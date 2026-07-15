############################################################
# BRITE C and BRITE D temporal trajectory clustering
#
# Input:
#   Results/BRITE_temporal_txt_only_significant_only/
#
# The input files are produced by:
#   Scripts/06_GAMS_and_Temporal_plots.R
#
# Workflow:
#   1. Analyze BRITE C and BRITE D separately.
#   2. Retain converged GAMs with FDR-adjusted p < 0.05.
#   3. Cluster z-scored GAM-fitted trajectories by temporal shape.
#   4. Use 1 - Pearson correlation and Ward.D2 clustering.
#   5. Evaluate k with an elbow plot and select the elbow automatically.
#   6. Allow the selected k to be overridden manually below.
#   7. Produce both an all-feature and reduced, readable heatmap.
############################################################

rm(list = ls())

# ---- User settings ----------------------------------------------------------

run_label <- "txt_only_significant_only"

# Leave as NA to use the automatically detected elbow.
# Enter an integer, for example 3, to override one result after reviewing
# the elbow plot.
manual_k <- c(
  brite_C = NA_integer_,
  brite_D = NA_integer_
)

# Number of strongest significant trajectories in the reduced heatmap.
top_n_heatmap <- c(
  brite_C = 60,
  brite_D = 75
)

# Maximum k evaluated in the elbow diagnostic.
max_k_to_test <- 12

# Show labels on the all-feature heatmap only when it remains readable.
max_rows_with_labels <- 100

# ---- Setup ------------------------------------------------------------------

required_packages <- c(
  "tidyverse",
  "dendextend",
  "pheatmap",
  "RColorBrewer",
  "magick"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(dendextend)
  library(pheatmap)
  library(RColorBrewer)
  library(magick)
})

theme_set(theme_bw(base_size = 12))

input_dir <- file.path(
  "Results",
  paste0("BRITE_temporal_", run_label)
)

results_dir <- file.path("Results", "BRITE_clustering_elbow")
figures_dir <- file.path("Figures", "BRITE_clustering_elbow")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Helper functions -------------------------------------------------------

shorten_brite_label <- function(x, level_name) {
  x <- as.character(x)

  if (level_name == "brite_C") {
    x <- str_remove(x, "^\\d{5}\\s+")
    x <- str_remove(x, "\\s+\\[(PATH|BR):[^]]+\\]$")
  } else {
    ko_code <- str_extract(x, "^K\\d{5}")
    gene_text <- x %>%
      str_remove("^K\\d{5}\\s+") %>%
      str_extract("^[^;]+") %>%
      str_squish()

    x <- if_else(
      is.na(gene_text) | gene_text == "",
      ko_code,
      paste(ko_code, gene_text)
    )
  }

  make.unique(x)
}

row_zscore <- function(x) {
  row_mean <- rowMeans(x, na.rm = TRUE)
  row_sd <- apply(x, 1, sd, na.rm = TRUE)
  sweep(sweep(x, 1, row_mean, "-"), 1, row_sd, "/")
}

# Elbow-style within-cluster dissimilarity based on the same correlation
# distance matrix used for clustering.
total_within_dissimilarity <- function(cluster_ids, distance_matrix) {
  distance_full <- as.matrix(distance_matrix)

  sum(vapply(unique(cluster_ids), function(cluster_id) {
    members <- names(cluster_ids)[cluster_ids == cluster_id]

    if (length(members) <= 1) {
      return(0)
    }

    cluster_distances <- distance_full[
      members,
      members,
      drop = FALSE
    ]

    sum(cluster_distances) / (2 * length(members))
  }, numeric(1)))
}

# Select the point farthest from the straight line joining the first and last
# normalized elbow points. The plot and selected value are both exported so
# the automatic choice can be reviewed and overridden above.
detect_elbow_k <- function(elbow_table) {
  if (nrow(elbow_table) <= 2) {
    return(max(elbow_table$k))
  }

  x <- (elbow_table$k - min(elbow_table$k)) /
    (max(elbow_table$k) - min(elbow_table$k))

  y_range <- diff(range(elbow_table$total_within_dissimilarity))

  if (y_range == 0) {
    return(2L)
  }

  y <- (elbow_table$total_within_dissimilarity -
          min(elbow_table$total_within_dissimilarity)) / y_range

  distance_from_line <- abs(y - (1 - x)) / sqrt(2)
  interior <- seq_len(nrow(elbow_table))[-c(1, nrow(elbow_table))]

  elbow_table$k[interior[which.max(distance_from_line[interior])]]
}

# cutree() cluster numbers are arbitrary. Relabel clusters according to the
# timing of their mean peak so Cluster_1 is the earliest-peaking group.
standardize_cluster_labels <- function(cluster_ids, zscore_matrix) {
  assignment_table <- tibble(
    brite_class = names(cluster_ids),
    raw_cluster = as.integer(cluster_ids)
  )

  cluster_order <- assignment_table %>%
    mutate(
      peak_index = apply(
        zscore_matrix[brite_class, , drop = FALSE],
        1,
        which.max
      )
    ) %>%
    group_by(raw_cluster) %>%
    summarise(
      mean_peak_index = mean(peak_index),
      .groups = "drop"
    ) %>%
    arrange(mean_peak_index, raw_cluster) %>%
    mutate(Cluster_ID = paste0("Cluster_", row_number()))

  standardized <- assignment_table %>%
    left_join(cluster_order, by = "raw_cluster")

  setNames(standardized$Cluster_ID, standardized$brite_class)
}

save_dendrogram <- function(
    hc,
    assignments,
    cluster_colors,
    level_name,
    output_stem,
    show_labels) {

  dend <- as.dendrogram(hc) %>%
    set("branches_lwd", 0.9)

  if (show_labels) {
    label_lookup <- setNames(
      assignments$display_label,
      assignments$brite_class
    )
    labels(dend) <- label_lookup[labels(dend)]
    dend <- dend %>% set("labels_cex", 0.55)
  } else {
    dend <- dend %>% set("labels", rep("", length(labels(dend))))
  }

  draw_dendrogram <- function() {
    plot(
      dend,
      main = paste0(
        str_replace(level_name, "_", " "),
        " temporal trajectory clustering"
      ),
      ylab = "Ward.D2 height",
      leaflab = if (show_labels) "perpendicular" else "none"
    )

    dend_order <- hc$labels[hc$order]

    left_to_right <- assignments %>%
      mutate(dend_position = match(brite_class, dend_order)) %>%
      group_by(Cluster_ID) %>%
      summarise(
        first_position = min(dend_position),
        .groups = "drop"
      ) %>%
      arrange(first_position) %>%
      pull(Cluster_ID)

    rect.hclust(
      hc,
      k = n_distinct(assignments$Cluster_ID),
      border = cluster_colors[left_to_right]
    )
  }

  png(
    paste0(output_stem, ".png"),
    width = 12,
    height = 7,
    units = "in",
    res = 300
  )
  draw_dendrogram()
  dev.off()

  pdf(paste0(output_stem, ".pdf"), width = 12, height = 7)
  draw_dendrogram()
  dev.off()
}

save_heatmap <- function(
    heatmap_matrix,
    assignments,
    cluster_colors,
    level_name,
    output_stem,
    title,
    show_labels,
    height) {

  ordered_features <- assignments$brite_class
  plot_matrix <- heatmap_matrix[ordered_features, , drop = FALSE]
  rownames(plot_matrix) <- assignments$display_label

  annotation_row <- assignments %>%
    select(display_label, Cluster_ID) %>%
    column_to_rownames("display_label")

  annotation_colors <- list(Cluster_ID = cluster_colors)

  draw_heatmap <- function(filename) {
    pheatmap(
      mat = plot_matrix,
      filename = filename,
      width = 11,
      height = height,
      scale = "none",
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      annotation_row = annotation_row,
      annotation_colors = annotation_colors,
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      breaks = seq(-3, 3, length.out = 101),
      show_rownames = show_labels,
      fontsize_row = 6,
      fontsize_col = 10,
      border_color = NA,
      main = title
    )
  }

  draw_heatmap(paste0(output_stem, ".png"))
  draw_heatmap(paste0(output_stem, ".pdf"))
}

# ---- Main analysis function -------------------------------------------------

analyze_brite_clustering <- function(level_name) {
  message("\n=== ", level_name, " clustering ===")

  level_results_dir <- file.path(results_dir, level_name)
  level_figures_dir <- file.path(figures_dir, level_name)

  dir.create(level_results_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(level_figures_dir, showWarnings = FALSE, recursive = TRUE)

  summary_file <- file.path(
    input_dir,
    paste0(level_name, "_GAM_summary_table.csv")
  )

  trajectory_file <- file.path(
    input_dir,
    paste0(level_name, "_GAM_fitted_trajectories.csv")
  )

  if (!file.exists(summary_file) || !file.exists(trajectory_file)) {
    stop(
      "Missing GAM inputs for ", level_name,
      ". Run Scripts/06_GAMS_and_Temporal_plots.R first."
    )
  }

  gam_summary <- read_csv(summary_file, show_col_types = FALSE)
  fitted_trajectories <- read_csv(trajectory_file, show_col_types = FALSE)

  significant_features <- gam_summary %>%
    filter(
      converged == "Yes",
      significant_time == "Yes",
      !is.na(p_adj_fdr),
      p_adj_fdr < 0.05
    ) %>%
    distinct(brite_class, .keep_all = TRUE)

  if (nrow(significant_features) < 3) {
    stop(
      level_name,
      " has fewer than three significant, converged trajectories."
    )
  }

  trajectory_long <- fitted_trajectories %>%
    semi_join(significant_features, by = "brite_class") %>%
    select(brite_class, Collection_date, time_numeric, fitted) %>%
    distinct() %>%
    mutate(Collection_date = as.Date(Collection_date)) %>%
    arrange(brite_class, Collection_date)

  duplicate_check <- trajectory_long %>%
    count(brite_class, Collection_date) %>%
    filter(n > 1)

  if (nrow(duplicate_check) > 0) {
    stop(
      "Duplicate fitted values were found for a BRITE class and date in ",
      level_name, "."
    )
  }

  date_lookup <- trajectory_long %>%
    distinct(Collection_date, time_numeric) %>%
    arrange(Collection_date) %>%
    mutate(
      time_label = format(Collection_date, "%b'%y"),
      matrix_col = paste0("day_", time_numeric)
    )

  trajectory_matrix <- trajectory_long %>%
    mutate(matrix_col = paste0("day_", time_numeric)) %>%
    select(brite_class, matrix_col, fitted) %>%
    pivot_wider(names_from = matrix_col, values_from = fitted) %>%
    column_to_rownames("brite_class") %>%
    as.matrix()

  trajectory_matrix <- trajectory_matrix[
    , date_lookup$matrix_col,
    drop = FALSE
  ]

  if (anyNA(trajectory_matrix)) {
    stop(
      "The fitted trajectory matrix contains missing values for ",
      level_name, "."
    )
  }

  trajectory_sd <- apply(trajectory_matrix, 1, sd)
  zero_variance <- names(trajectory_sd)[
    is.na(trajectory_sd) | trajectory_sd == 0
  ]

  if (length(zero_variance) > 0) {
    warning(
      "Removing zero-variance fitted trajectories from ", level_name,
      ": ", paste(zero_variance, collapse = ", ")
    )
    trajectory_matrix <- trajectory_matrix[
      !rownames(trajectory_matrix) %in% zero_variance,
      ,
      drop = FALSE
    ]
    significant_features <- significant_features %>%
      filter(!brite_class %in% zero_variance)
  }

  zscore_matrix <- row_zscore(trajectory_matrix)
  zscore_matrix[zscore_matrix < -3] <- -3
  zscore_matrix[zscore_matrix > 3] <- 3
  colnames(zscore_matrix) <- date_lookup$time_label

  write_csv(
    as.data.frame(zscore_matrix) %>%
      rownames_to_column("brite_class"),
    file.path(level_results_dir, "zscored_fitted_trajectories.csv")
  )

  correlation_matrix <- cor(
    t(zscore_matrix),
    method = "pearson",
    use = "pairwise.complete.obs"
  )

  if (anyNA(correlation_matrix)) {
    stop(
      "The correlation matrix contains missing values for ",
      level_name, "."
    )
  }

  dissimilarity_matrix <- 1 - correlation_matrix
  dissimilarity_matrix[dissimilarity_matrix < 0] <- 0
  diag(dissimilarity_matrix) <- 0
  distance_object <- as.dist(dissimilarity_matrix)

  hc <- hclust(distance_object, method = "ward.D2")

  n_features <- nrow(zscore_matrix)
  k_max <- min(max_k_to_test, n_features - 1)

  if (k_max < 2) {
    stop("Too few trajectories to evaluate clustering for ", level_name, ".")
  }

  elbow_table <- map_dfr(1:k_max, function(k) {
    cluster_ids <- cutree(hc, k = k)

    tibble(
      k = k,
      total_within_dissimilarity = total_within_dissimilarity(
        cluster_ids,
        distance_object
      )
    )
  })

  automatic_k <- detect_elbow_k(elbow_table)
  requested_k <- unname(manual_k[level_name])

  selected_k <- if (
    length(requested_k) == 1 && !is.na(requested_k)
  ) {
    as.integer(requested_k)
  } else {
    as.integer(automatic_k)
  }

  if (selected_k < 2 || selected_k > k_max) {
    stop(
      "Selected k for ", level_name,
      " must be between 2 and ", k_max, "."
    )
  }

  elbow_table <- elbow_table %>%
    mutate(
      automatic_elbow_k = automatic_k,
      selected_k = selected_k,
      selected = k == selected_k
    )

  write_csv(
    elbow_table,
    file.path(level_results_dir, "elbow_diagnostics.csv")
  )

  elbow_plot <- ggplot(
    elbow_table,
    aes(k, total_within_dissimilarity)
  ) +
    geom_line(linewidth = 0.8, color = "grey30") +
    geom_point(size = 2.5, color = "grey20") +
    geom_vline(
      xintercept = selected_k,
      linetype = "dashed",
      color = "#0072B2",
      linewidth = 0.7
    ) +
    scale_x_continuous(breaks = 1:k_max) +
    labs(
      title = paste0(str_replace(level_name, "_", " "), " elbow diagnostic"),
      subtitle = paste0(
        "Automatically detected elbow: k = ", automatic_k,
        "; selected k = ", selected_k
      ),
      x = "Number of clusters (k)",
      y = "Total within-cluster dissimilarity"
    ) +
    theme(panel.grid.minor = element_blank())

  ggsave(
    file.path(level_figures_dir, "elbow_diagnostic.png"),
    elbow_plot,
    width = 8,
    height = 5.5,
    dpi = 300
  )
  ggsave(
    file.path(level_figures_dir, "elbow_diagnostic.pdf"),
    elbow_plot,
    width = 8,
    height = 5.5
  )

  raw_cluster_ids <- cutree(hc, k = selected_k)
  standardized_ids <- standardize_cluster_labels(
    raw_cluster_ids,
    zscore_matrix
  )

  cluster_colors <- setNames(
    hcl.colors(selected_k, palette = "Dark 3"),
    paste0("Cluster_", seq_len(selected_k))
  )

  assignments <- tibble(
    brite_class = names(raw_cluster_ids),
    raw_cluster = as.integer(raw_cluster_ids),
    Cluster_ID = standardized_ids[names(raw_cluster_ids)]
  ) %>%
    left_join(
      significant_features %>%
        select(
          brite_class,
          p_value,
          p_adj_fdr,
          adj_r_squared,
          F_statistic,
          edf
        ),
      by = "brite_class"
    ) %>%
    mutate(
      display_label = shorten_brite_label(brite_class, level_name),
      dend_position = match(brite_class, hc$labels[hc$order])
    ) %>%
    arrange(Cluster_ID, dend_position)

  write_csv(
    assignments %>% select(-dend_position),
    file.path(level_results_dir, "cluster_assignments.csv")
  )

  cluster_summary <- assignments %>%
    count(Cluster_ID, name = "n_features") %>%
    mutate(
      brite_level = level_name,
      selected_k = selected_k,
      automatic_elbow_k = automatic_k
    ) %>%
    relocate(brite_level, selected_k, automatic_elbow_k)

  write_csv(
    cluster_summary,
    file.path(level_results_dir, "cluster_summary.csv")
  )

  dendrogram_stem <- file.path(
    level_figures_dir,
    paste0("dendrogram_k", selected_k)
  )

  save_dendrogram(
    hc = hc,
    assignments = assignments,
    cluster_colors = cluster_colors,
    level_name = level_name,
    output_stem = dendrogram_stem,
    show_labels = n_features <= max_rows_with_labels
  )

  trajectory_plot_data <- as.data.frame(zscore_matrix) %>%
    rownames_to_column("brite_class") %>%
    pivot_longer(
      -brite_class,
      names_to = "time_label",
      values_to = "zscore_fitted"
    ) %>%
    left_join(
      assignments %>% select(brite_class, Cluster_ID),
      by = "brite_class"
    ) %>%
    mutate(
      time_label = factor(time_label, levels = date_lookup$time_label),
      Cluster_ID = factor(
        Cluster_ID,
        levels = names(cluster_colors)
      )
    )

  trajectory_summary <- trajectory_plot_data %>%
    group_by(Cluster_ID, time_label) %>%
    summarise(
      mean_z = mean(zscore_fitted),
      se_z = if_else(
        n() > 1,
        sd(zscore_fitted) / sqrt(n()),
        0
      ),
      n_features = n(),
      .groups = "drop"
    )

  write_csv(
    trajectory_summary,
    file.path(level_results_dir, "cluster_average_trajectories.csv")
  )

  trajectory_plot <- ggplot(
    trajectory_summary,
    aes(
      time_label,
      mean_z,
      group = Cluster_ID,
      color = Cluster_ID,
      fill = Cluster_ID
    )
  ) +
    geom_ribbon(
      aes(ymin = mean_z - se_z, ymax = mean_z + se_z),
      alpha = 0.18,
      color = NA
    ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = cluster_colors, drop = FALSE) +
    scale_fill_manual(values = cluster_colors, drop = FALSE) +
    labs(
      title = paste0(
        str_replace(level_name, "_", " "),
        " cluster-average fitted trajectories"
      ),
      subtitle = paste0("Elbow-selected k = ", selected_k),
      x = NULL,
      y = "Mean z-scored fitted expression",
      color = "Cluster",
      fill = "Cluster"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    )

  trajectory_stem <- file.path(
    level_figures_dir,
    paste0("cluster_average_trajectories_k", selected_k)
  )

  ggsave(
    paste0(trajectory_stem, ".png"),
    trajectory_plot,
    width = 9,
    height = 5.5,
    dpi = 300
  )
  ggsave(
    paste0(trajectory_stem, ".pdf"),
    trajectory_plot,
    width = 9,
    height = 5.5
  )

  all_heatmap_height <- min(40, max(9, 0.13 * n_features))
  all_heatmap_stem <- file.path(
    level_figures_dir,
    paste0("heatmap_all_significant_k", selected_k)
  )

  save_heatmap(
    heatmap_matrix = zscore_matrix,
    assignments = assignments,
    cluster_colors = cluster_colors,
    level_name = level_name,
    output_stem = all_heatmap_stem,
    title = paste0(
      str_replace(level_name, "_", " "),
      ": all significant fitted trajectories"
    ),
    show_labels = n_features <= max_rows_with_labels,
    height = all_heatmap_height
  )

  n_top <- min(unname(top_n_heatmap[level_name]), n_features)

  top_assignments <- assignments %>%
    mutate(
      strength_rank = coalesce(adj_r_squared, -Inf),
      f_rank = coalesce(F_statistic, -Inf)
    ) %>%
    arrange(desc(strength_rank), desc(f_rank), p_adj_fdr) %>%
    slice_head(n = n_top) %>%
    arrange(Cluster_ID, dend_position)

  write_csv(
    top_assignments %>%
      select(-dend_position, -strength_rank, -f_rank),
    file.path(
      level_results_dir,
      paste0("top_", n_top, "_heatmap_features.csv")
    )
  )

  top_heatmap_height <- min(20, max(9, 0.22 * n_top))
  top_heatmap_stem <- file.path(
    level_figures_dir,
    paste0("heatmap_top_", n_top, "_k", selected_k)
  )

  save_heatmap(
    heatmap_matrix = zscore_matrix,
    assignments = top_assignments,
    cluster_colors = cluster_colors,
    level_name = level_name,
    output_stem = top_heatmap_stem,
    title = paste0(
      str_replace(level_name, "_", " "),
      ": strongest ", n_top, " significant temporal signals"
    ),
    show_labels = TRUE,
    height = top_heatmap_height
  )

  # The readable combined figure uses the reduced heatmap. The all-feature
  # heatmap remains available as a separate complete overview.
  dendrogram_image <- image_read(paste0(dendrogram_stem, ".png")) %>%
    image_resize("2400x")
  trajectory_image <- image_read(paste0(trajectory_stem, ".png")) %>%
    image_resize("2400x")
  top_heatmap_image <- image_read(paste0(top_heatmap_stem, ".png")) %>%
    image_resize("2400x")

  combined_image <- image_append(
    c(dendrogram_image, trajectory_image, top_heatmap_image),
    stack = TRUE
  )

  image_write(
    combined_image,
    file.path(
      level_figures_dir,
      paste0("combined_clustering_figure_k", selected_k, ".png")
    )
  )
  image_write(
    combined_image,
    file.path(
      level_figures_dir,
      paste0("combined_clustering_figure_k", selected_k, ".pdf")
    )
  )

  tibble(
    brite_level = level_name,
    n_significant_input = nrow(significant_features),
    n_clustered = n_features,
    automatic_elbow_k = automatic_k,
    selected_k = selected_k,
    manual_override_used = !is.na(requested_k),
    reduced_heatmap_n = n_top
  )
}

# ---- Run BRITE C and BRITE D separately ------------------------------------

run_summary <- bind_rows(
  analyze_brite_clustering("brite_C"),
  analyze_brite_clustering("brite_D")
)

write_csv(
  run_summary,
  file.path(results_dir, "BRITE_clustering_run_summary.csv")
)

cat("\nBRITE clustering complete.\n")
print(run_summary)
cat("\nResults:", results_dir, "\n")
cat("Figures:", figures_dir, "\n")
