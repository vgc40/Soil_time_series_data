############################################################
# MicroResp-relevant BRITE C temporal trajectory clustering
#
# Required first steps:
#   1. Run Scripts/09_BRITE_GAM_random_plot_all_pathways.R
#   2. Run Scripts/10_MicroResp_targeted_BRITE_C_GAM_analysis.R
#
# Scope:
#   - BRITE C pathways explicitly mapped in Data/Microresp_pathways.csv
#   - Converged models only
#   - All mapped pathways are clustered, regardless of temporal significance
#   - Targeted-FDR significance is retained as an annotation
#
# Clustering:
#   - Population GAM trajectories with Plot random effects excluded
#   - Row z-score to compare temporal shape rather than expression magnitude
#   - 1 - Pearson correlation dissimilarity
#   - Ward.D2 hierarchical clustering
#   - Elbow-selected number of clusters, with an optional manual override
############################################################

rm(list = ls())

# ---- User settings ----------------------------------------------------------

# Leave as NA to use the automatically detected elbow.
manual_k <- NA_integer_

max_k_to_test <- 10L

input_dir <- file.path("Results", "MicroResp_targeted_BRITE_C_gam")
results_dir <- file.path("Results", "MicroResp_BRITE_C_clustering_elbow")
figures_dir <- file.path("Figures", "MicroResp_BRITE_C_clustering_elbow")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Setup ------------------------------------------------------------------

required_packages <- c(
  "tidyverse",
  "dendextend",
  "pheatmap",
  "RColorBrewer"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org",
    dependencies = TRUE
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(dendextend)
  library(pheatmap)
  library(RColorBrewer)
})

theme_set(theme_bw(base_size = 12))

summary_file <- file.path(
  input_dir,
  "MicroResp_target_BRITE_C_GAM_summary.csv"
)

trajectory_file <- file.path(
  input_dir,
  "MicroResp_target_BRITE_C_fitted_trajectories.csv"
)

mapping_file <- file.path(
  input_dir,
  "MicroResp_to_BRITE_C_mapping_audit.csv"
)

required_inputs <- c(summary_file, trajectory_file, mapping_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs) > 0) {
  stop(
    "Missing required input(s): ",
    paste(missing_inputs, collapse = ", "),
    ". Run Scripts/09 and 10 first."
  )
}

# ---- Helper functions -------------------------------------------------------

row_zscore <- function(x) {
  row_mean <- rowMeans(x, na.rm = TRUE)
  row_sd <- apply(x, 1, sd, na.rm = TRUE)
  sweep(sweep(x, 1, row_mean, "-"), 1, row_sd, "/")
}

short_pathway_label <- function(x) {
  x %>%
    str_remove("^\\d{5}\\s+") %>%
    str_remove("\\s+\\[PATH:ko\\d{5}\\]$")
}

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

detect_elbow_k <- function(elbow_table) {
  if (nrow(elbow_table) <= 2) {
    return(max(elbow_table$k))
  }

  x <- (elbow_table$k - min(elbow_table$k)) /
    (max(elbow_table$k) - min(elbow_table$k))

  y_range <- diff(range(elbow_table$total_within_dissimilarity))

  if (!is.finite(y_range) || y_range == 0) {
    return(2L)
  }

  y <- (elbow_table$total_within_dissimilarity -
    min(elbow_table$total_within_dissimilarity)) / y_range

  distance_from_line <- abs(y - (1 - x)) / sqrt(2)
  interior <- seq_len(nrow(elbow_table))[-c(1, nrow(elbow_table))]

  as.integer(
    elbow_table$k[interior[which.max(distance_from_line[interior])]]
  )
}

standardize_cluster_labels <- function(cluster_ids, zscore_matrix) {
  assignment_table <- tibble(
    pathway_code = names(cluster_ids),
    raw_cluster = as.integer(cluster_ids)
  )

  cluster_order <- assignment_table %>%
    mutate(
      peak_index = apply(
        zscore_matrix[pathway_code, , drop = FALSE],
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

  setNames(standardized$Cluster_ID, standardized$pathway_code)
}

save_heatmap <- function(
    heatmap_matrix,
    assignments,
    cluster_colors,
    filename_stem,
    title) {

  if (nrow(assignments) == 0) {
    return(invisible(NULL))
  }

  plot_matrix <- heatmap_matrix[
    assignments$pathway_code,
    ,
    drop = FALSE
  ]

  rownames(plot_matrix) <- assignments$display_label

  annotation_row <- assignments %>%
    select(display_label, Cluster_ID, Temporal_significance) %>%
    column_to_rownames("display_label")

  significance_colors <- c(
    Significant = "#238b45",
    `Not significant` = "#bdbdbd"
  )

  annotation_colors <- list(
    Cluster_ID = cluster_colors,
    Temporal_significance = significance_colors
  )

  plot_height <- max(8, 0.27 * nrow(assignments))

  for (extension in c("png", "pdf")) {
    pheatmap(
      mat = plot_matrix,
      filename = paste0(filename_stem, ".", extension),
      width = 11,
      height = plot_height,
      scale = "none",
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      annotation_row = annotation_row,
      annotation_colors = annotation_colors,
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      breaks = seq(-3, 3, length.out = 101),
      show_rownames = TRUE,
      fontsize_row = 7,
      fontsize_col = 10,
      border_color = NA,
      main = title
    )
  }
}

# ---- Read targeted results --------------------------------------------------

target_summary <- read_csv(summary_file, show_col_types = FALSE) %>%
  filter(
    fit_status == "Converged",
    converged == "Yes",
    !is.na(pathway_code)
  ) %>%
  distinct(pathway_code, .keep_all = TRUE)

target_trajectories <- read_csv(
  trajectory_file,
  show_col_types = FALSE
) %>%
  mutate(Collection_date = as.Date(Collection_date)) %>%
  semi_join(target_summary, by = "pathway_code") %>%
  select(
    pathway_code,
    brite_class,
    Collection_date,
    time_numeric,
    fitted
  ) %>%
  distinct() %>%
  arrange(pathway_code, Collection_date)

mapping_audit <- read_csv(mapping_file, show_col_types = FALSE)

if (nrow(target_summary) < 3) {
  stop("Fewer than three converged MicroResp-linked pathways are available.")
}

duplicate_check <- target_trajectories %>%
  count(pathway_code, Collection_date) %>%
  filter(n > 1)

if (nrow(duplicate_check) > 0) {
  stop("Duplicate fitted values were found for a pathway and date.")
}

date_lookup <- target_trajectories %>%
  distinct(Collection_date, time_numeric) %>%
  arrange(Collection_date) %>%
  mutate(
    time_label = format(Collection_date, "%b'%y"),
    matrix_col = paste0("day_", time_numeric)
  )

trajectory_matrix <- target_trajectories %>%
  mutate(matrix_col = paste0("day_", time_numeric)) %>%
  select(pathway_code, matrix_col, fitted) %>%
  pivot_wider(names_from = matrix_col, values_from = fitted) %>%
  column_to_rownames("pathway_code") %>%
  as.matrix()

trajectory_matrix <- trajectory_matrix[
  , date_lookup$matrix_col,
  drop = FALSE
]

if (anyNA(trajectory_matrix)) {
  stop("The MicroResp-linked fitted trajectory matrix contains missing values.")
}

trajectory_sd <- apply(trajectory_matrix, 1, sd)
zero_variance <- names(trajectory_sd)[
  !is.finite(trajectory_sd) | trajectory_sd == 0
]

if (length(zero_variance) > 0) {
  warning(
    "Removing zero-variance fitted trajectories: ",
    paste(zero_variance, collapse = ", ")
  )

  trajectory_matrix <- trajectory_matrix[
    !rownames(trajectory_matrix) %in% zero_variance,
    ,
    drop = FALSE
  ]

  target_summary <- target_summary %>%
    filter(!pathway_code %in% zero_variance)
}

zscore_matrix <- row_zscore(trajectory_matrix)
zscore_matrix[zscore_matrix < -3] <- -3
zscore_matrix[zscore_matrix > 3] <- 3
colnames(zscore_matrix) <- date_lookup$time_label

write_csv(
  as.data.frame(zscore_matrix) %>%
    rownames_to_column("pathway_code"),
  file.path(results_dir, "zscored_MicroResp_BRITE_C_trajectories.csv")
)

# ---- Correlation clustering and elbow --------------------------------------

correlation_matrix <- cor(
  t(zscore_matrix),
  method = "pearson",
  use = "pairwise.complete.obs"
)

if (anyNA(correlation_matrix)) {
  stop("The pathway correlation matrix contains missing values.")
}

dissimilarity_matrix <- 1 - correlation_matrix
dissimilarity_matrix[dissimilarity_matrix < 0] <- 0
diag(dissimilarity_matrix) <- 0
distance_object <- as.dist(dissimilarity_matrix)

hc <- hclust(distance_object, method = "ward.D2")

n_pathways <- nrow(zscore_matrix)
k_max <- min(max_k_to_test, n_pathways - 1L)

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
selected_k <- if (!is.na(manual_k)) as.integer(manual_k) else automatic_k

if (selected_k < 2 || selected_k > k_max) {
  stop("Selected k must be between 2 and ", k_max, ".")
}

elbow_table <- elbow_table %>%
  mutate(
    automatic_elbow_k = automatic_k,
    selected_k = selected_k,
    selected = k == selected_k
  )

write_csv(
  elbow_table,
  file.path(results_dir, "MicroResp_BRITE_C_elbow_diagnostics.csv")
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
    title = "MicroResp-linked BRITE C elbow diagnostic",
    subtitle = paste0(
      "Automatically detected elbow: k = ", automatic_k,
      "; selected k = ", selected_k
    ),
    x = "Number of clusters (k)",
    y = "Total within-cluster dissimilarity"
  ) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(figures_dir, "MicroResp_BRITE_C_elbow_diagnostic.png"),
  elbow_plot,
  width = 8,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(figures_dir, "MicroResp_BRITE_C_elbow_diagnostic.pdf"),
  elbow_plot,
  width = 8,
  height = 5.5
)

# ---- Assignments and annotations -------------------------------------------

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
  pathway_code = names(raw_cluster_ids),
  raw_cluster = as.integer(raw_cluster_ids),
  Cluster_ID = standardized_ids[names(raw_cluster_ids)]
) %>%
  left_join(
    target_summary %>%
      select(
        pathway_code,
        brite_class,
        linked_substrates,
        n_linked_substrates,
        p_value,
        p_adj_fdr,
        p_adj_fdr_microresp_targets,
        significant_time,
        significant_time_microresp_targets,
        deviance_explained_full,
        time_variance_explained,
        edf
      ),
    by = "pathway_code"
  ) %>%
  mutate(
    Temporal_significance = if_else(
      significant_time_microresp_targets == "Yes",
      "Significant",
      "Not significant"
    ),
    display_label = paste0(
      str_remove(pathway_code, "^ko"),
      " ",
      short_pathway_label(brite_class)
    ),
    dend_position = match(pathway_code, hc$labels[hc$order])
  ) %>%
  arrange(Cluster_ID, dend_position)

write_csv(
  assignments %>% select(-dend_position),
  file.path(results_dir, "MicroResp_BRITE_C_cluster_assignments.csv")
)

substrate_cluster_map <- mapping_audit %>%
  select(
    Substrate,
    Substrate_category,
    Table_group,
    pathway_code,
    mapped_brite_class = brite_C
  ) %>%
  inner_join(
    assignments %>%
      select(
        pathway_code,
        Cluster_ID,
        Temporal_significance,
        p_adj_fdr_microresp_targets,
        time_variance_explained
      ),
    by = "pathway_code"
  ) %>%
  arrange(Substrate, Cluster_ID, pathway_code)

write_csv(
  substrate_cluster_map,
  file.path(results_dir, "MicroResp_substrate_to_BRITE_C_clusters.csv")
)

cluster_summary <- assignments %>%
  group_by(Cluster_ID) %>%
  summarise(
    n_pathways = n(),
    n_temporally_significant = sum(
      Temporal_significance == "Significant"
    ),
    mean_time_variance_explained = mean(
      time_variance_explained,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  mutate(
    selected_k = selected_k,
    automatic_elbow_k = automatic_k
  ) %>%
  relocate(selected_k, automatic_elbow_k)

write_csv(
  cluster_summary,
  file.path(results_dir, "MicroResp_BRITE_C_cluster_summary.csv")
)

# ---- Dendrogram -------------------------------------------------------------

dend <- as.dendrogram(hc) %>%
  set("branches_lwd", 0.9)

label_lookup <- setNames(assignments$display_label, assignments$pathway_code)
labels(dend) <- label_lookup[labels(dend)]
dend <- dend %>% set("labels_cex", 0.6)

draw_dendrogram <- function() {
  plot(
    dend,
    main = "MicroResp-linked BRITE C temporal clustering",
    ylab = "Ward.D2 height",
    leaflab = "perpendicular"
  )

  dend_order <- hc$labels[hc$order]

  left_to_right <- assignments %>%
    mutate(dend_position_plot = match(pathway_code, dend_order)) %>%
    group_by(Cluster_ID) %>%
    summarise(
      first_position = min(dend_position_plot),
      .groups = "drop"
    ) %>%
    arrange(first_position) %>%
    pull(Cluster_ID)

  rect.hclust(
    hc,
    k = selected_k,
    border = cluster_colors[left_to_right]
  )
}

png(
  file.path(figures_dir, "MicroResp_BRITE_C_dendrogram.png"),
  width = 12,
  height = 7,
  units = "in",
  res = 300
)
draw_dendrogram()
dev.off()

pdf(
  file.path(figures_dir, "MicroResp_BRITE_C_dendrogram.pdf"),
  width = 12,
  height = 7
)
draw_dendrogram()
dev.off()

# ---- Cluster-average trajectories ------------------------------------------

trajectory_plot_data <- as.data.frame(zscore_matrix) %>%
  rownames_to_column("pathway_code") %>%
  pivot_longer(
    -pathway_code,
    names_to = "time_label",
    values_to = "zscore_fitted"
  ) %>%
  left_join(
    assignments %>% select(pathway_code, Cluster_ID),
    by = "pathway_code"
  ) %>%
  mutate(
    time_label = factor(time_label, levels = date_lookup$time_label),
    Cluster_ID = factor(Cluster_ID, levels = names(cluster_colors))
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
    n_pathways = n(),
    .groups = "drop"
  )

write_csv(
  trajectory_summary,
  file.path(results_dir, "MicroResp_BRITE_C_cluster_average_trajectories.csv")
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
    title = "MicroResp-linked BRITE C cluster-average trajectories",
    subtitle = paste0(
      "All mapped pathways; elbow-selected k = ", selected_k
    ),
    x = NULL,
    y = "Mean z-scored fitted expression",
    color = "Cluster",
    fill = "Cluster"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(
    figures_dir,
    "MicroResp_BRITE_C_cluster_average_trajectories.png"
  ),
  trajectory_plot,
  width = 9,
  height = 5.5,
  dpi = 300
)

ggsave(
  file.path(
    figures_dir,
    "MicroResp_BRITE_C_cluster_average_trajectories.pdf"
  ),
  trajectory_plot,
  width = 9,
  height = 5.5
)

# ---- Heatmaps ---------------------------------------------------------------

save_heatmap(
  heatmap_matrix = zscore_matrix,
  assignments = assignments,
  cluster_colors = cluster_colors,
  filename_stem = file.path(
    figures_dir,
    "MicroResp_BRITE_C_heatmap_all_mapped"
  ),
  title = "All MicroResp-linked BRITE C fitted trajectories"
)

significant_assignments <- assignments %>%
  filter(Temporal_significance == "Significant")

save_heatmap(
  heatmap_matrix = zscore_matrix,
  assignments = significant_assignments,
  cluster_colors = cluster_colors,
  filename_stem = file.path(
    figures_dir,
    "MicroResp_BRITE_C_heatmap_target_FDR_significant"
  ),
  title = "Target-FDR significant MicroResp-linked BRITE C trajectories"
)

# ---- Run summary ------------------------------------------------------------

run_summary <- tibble(
  n_mapping_rows = nrow(mapping_audit),
  n_unique_mapped_converged_pathways = nrow(target_summary),
  n_clustered_pathways = nrow(zscore_matrix),
  n_target_FDR_significant = sum(
    assignments$Temporal_significance == "Significant"
  ),
  automatic_elbow_k = automatic_k,
  selected_k = selected_k,
  manual_override_used = !is.na(manual_k),
  clustering_distance = "1 - Pearson correlation",
  clustering_method = "Ward.D2",
  trajectory_scale = "Row z-score of population GAM fitted trajectories"
)

write_csv(
  run_summary,
  file.path(results_dir, "MicroResp_BRITE_C_clustering_run_summary.csv")
)

cat("\nMicroResp-linked BRITE C clustering complete.\n")
print(run_summary)
cat("\nResults: ", results_dir, "\n", sep = "")
cat("Figures: ", figures_dir, "\n", sep = "")

