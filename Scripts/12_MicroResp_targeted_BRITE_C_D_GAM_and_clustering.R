############################################################
# MicroResp-targeted BRITE C and BRITE D GAM review and
# temporal trajectory clustering
#
# Required first step:
#   Run Scripts/09_BRITE_GAM_random_plot_all_pathways.R
#
# This script uses the explicit MicroResp-to-BRITE C/D mapping in
# Data/MicroResp_all_BRITE_C_D_complete.csv. It does NOT refit the GAMs:
# it subsets the already-corrected random-Plot GAM results from Script 09.
#
# For each BRITE level it:
#   1. audits the MicroResp mapping;
#   2. exports GAM statistics (p, FDR, adjusted R2, full-model deviance,
#      time-only variance explained, EDF, and fit status);
#   3. calculates a separate FDR across the MicroResp-targeted features;
#   4. clusters every mapped feature with a converged GAM, regardless of
#      temporal significance;
#   5. produces GAM panels, dendrograms, cluster-average trajectories,
#      and heatmaps.
#
# BRITE D heatmaps are reduced for readability by selecting the strongest
# time-only trajectories within each cluster. That selection is only for
# display: all mapped, converged BRITE D features remain in clustering.
############################################################

rm(list = ls())

# ---- User settings ----------------------------------------------------------

mapping_file <- file.path("Data", "MicroResp_all_BRITE_C_D_complete.csv")
all_gam_dir <- file.path("Results", "BRITE_temporal_gam_random_plot_all_pathways")

# Leave as NA to use the elbow detected from correlation distance.
manual_k <- c(brite_C = NA_integer_, brite_D = NA_integer_)
max_k_to_test <- 10L

# BRITE D heatmap selection: at most this many features per cluster, ranked
# by time-only variance explained. The final number can be smaller than the
# cap when a cluster contains fewer features.
top_brite_D_per_cluster <- 15L

# Six GAM panels per PDF page, matching the earlier MicroResp workflow.
gam_panels_per_page <- 6L

results_dir <- file.path("Results", "MicroResp_targeted_BRITE_C_D")
figures_dir <- file.path("Figures", "MicroResp_targeted_BRITE_C_D")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Packages ---------------------------------------------------------------

required_packages <- c("tidyverse", "dendextend", "pheatmap", "RColorBrewer", "ggforce")
missing_packages <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org", dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(dendextend)
  library(pheatmap)
  library(RColorBrewer)
  library(ggforce)
})
theme_set(theme_bw(base_size = 11))

# ---- Helpers ----------------------------------------------------------------

standard_name <- function(x) str_to_lower(str_replace_all(x, "[^a-zA-Z0-9]", ""))

find_column <- function(data, options, required = TRUE) {
  names_standard <- standard_name(names(data))
  option_standard <- standard_name(options)
  match_index <- match(option_standard, names_standard, nomatch = 0)
  match_index <- match_index[match_index > 0]
  if (length(match_index) == 0) {
    if (required) stop("Could not find a required mapping column. Expected one of: ", paste(options, collapse = ", "))
    return(NULL)
  }
  names(data)[match_index[1]]
}

row_zscore <- function(x) {
  means <- rowMeans(x, na.rm = TRUE)
  sds <- apply(x, 1, sd, na.rm = TRUE)
  sweep(sweep(x, 1, means, "-"), 1, sds, "/")
}

total_within_dissimilarity <- function(cluster_ids, distance_matrix) {
  full <- as.matrix(distance_matrix)
  sum(vapply(unique(cluster_ids), function(cluster_id) {
    members <- names(cluster_ids)[cluster_ids == cluster_id]
    if (length(members) < 2) return(0)
    sum(full[members, members, drop = FALSE]) / (2 * length(members))
  }, numeric(1)))
}

detect_elbow_k <- function(tab) {
  if (nrow(tab) <= 2) return(max(tab$k))
  x <- (tab$k - min(tab$k)) / diff(range(tab$k))
  yrange <- diff(range(tab$total_within_dissimilarity))
  if (!is.finite(yrange) || yrange == 0) return(2L)
  y <- (tab$total_within_dissimilarity - min(tab$total_within_dissimilarity)) / yrange
  interior <- seq_len(nrow(tab))[-c(1, nrow(tab))]
  tab$k[interior[which.max(abs(y[interior] - (1 - x[interior])) / sqrt(2))]]
}

short_label <- function(x, level) {
  if (level == "brite_C") {
    x %>% str_remove("^\\d{5}\\s+") %>% str_remove("\\s+\\[PATH:ko\\d{5}\\]$")
  } else {
    code <- str_extract(x, "^K\\d{5}")
    text <- x %>% str_remove("^K\\d{5}\\s+") %>% str_extract("^[^;]+") %>% str_squish()
    if_else(is.na(text) | text == "", code, paste(code, text))
  }
}

make_paginated_gam_plots <- function(points, predictions, level, level_dir) {
  labels <- predictions %>% distinct(feature_code, plot_label)
  n_panels <- nrow(labels)
  n_pages <- ceiling(n_panels / gam_panels_per_page)
  if (n_panels == 0) return(tibble(brite_level = level, n_panels = 0L, n_pages = 0L))

  pdf_file <- file.path(level_dir, paste0(level, "_targeted_GAMs_all_features.pdf"))
  pdf(pdf_file, width = 14, height = 8.5, onefile = TRUE)
  for (page in seq_len(n_pages)) {
    plot_page <- ggplot(points, aes(Collection_date, Expression)) +
      geom_point(alpha = 0.35, size = 0.9, color = "grey35") +
      geom_ribbon(
        data = predictions,
        aes(x = Collection_date, ymin = lower, ymax = upper),
        inherit.aes = FALSE, fill = "#9ecae1", alpha = 0.35
      ) +
      geom_line(
        data = predictions,
        aes(x = Collection_date, y = fitted),
        inherit.aes = FALSE, color = "#08519c", linewidth = 0.8
      ) +
      ggforce::facet_wrap_paginate(~plot_label, scales = "free_y", ncol = 3, nrow = 2, page = page) +
      labs(title = paste0("MicroResp-linked ", str_replace(level, "_", " "), " GAMs"),
           subtitle = paste0("All mapped, converged features. Population time trend, page ", page, " of ", n_pages),
           x = "Collection date", y = "Expression") +
      theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(size = 7))
    print(plot_page)
    ggsave(file.path(level_dir, paste0(level, "_targeted_GAMs_page_", sprintf("%02d", page), ".png")), plot_page, width = 14, height = 8.5, dpi = 300)
  }
  dev.off()
  tibble(brite_level = level, n_panels = n_panels, panels_per_page = gam_panels_per_page, n_pages = n_pages)
}

save_heatmap <- function(mat, assignments, level, suffix, title, output_dir) {
  if (nrow(assignments) == 0) return(invisible(NULL))
  plot_mat <- mat[assignments$feature_code, , drop = FALSE]
  rownames(plot_mat) <- make.unique(assignments$display_label)
  annotation <- assignments %>%
    transmute(display_label = rownames(plot_mat), Cluster_ID, Temporal_significance) %>%
    column_to_rownames("display_label")
  colors <- setNames(hcl.colors(n_distinct(assignments$Cluster_ID), "Dark 3"), sort(unique(assignments$Cluster_ID)))
  annotation_colors <- list(Cluster_ID = colors, Temporal_significance = c(Significant = "#238b45", `Not significant` = "#bdbdbd"))
  for (extension in c("png", "pdf")) {
    pheatmap(plot_mat, filename = file.path(output_dir, paste0(level, "_", suffix, ".", extension)),
             width = 11, height = max(7, 0.27 * nrow(assignments)), scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
             annotation_row = annotation, annotation_colors = annotation_colors,
             color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100), breaks = seq(-3, 3, length.out = 101),
             fontsize_row = if (nrow(assignments) > 80) 5 else 7, border_color = NA, main = title)
  }
}

# ---- Validate and standardize the mapping -----------------------------------

if (!file.exists(mapping_file)) {
  candidates <- list.files("Data", pattern = "^MicroResp_all_BRITE_C_D_complete\\.(csv|CSV)$", full.names = TRUE)
  if (length(candidates) == 1) mapping_file <- candidates
  if (!file.exists(mapping_file)) stop("Mapping file not found: ", mapping_file)
}

raw_map <- read_csv(mapping_file, show_col_types = FALSE)
substrate_col <- find_column(raw_map, c("Substrate", "MicroResp_substrate"), required = FALSE)
category_col <- find_column(raw_map, c("Substrate_category", "MicroResp_substrate_category"), required = FALSE)
brite_c_col <- find_column(raw_map, c("brite_C", "BRITE_C"))
brite_d_col <- find_column(raw_map, c("brite_D", "BRITE_D"))

mapping <- raw_map %>%
  transmute(
    Substrate = if (is.null(substrate_col)) "Not provided" else as.character(.data[[substrate_col]]),
    Substrate_category = if (is.null(category_col)) "Not provided" else as.character(.data[[category_col]]),
    brite_C_mapping = as.character(.data[[brite_c_col]]),
    brite_D_mapping = as.character(.data[[brite_d_col]]),
    brite_C_code = str_to_lower(str_extract(brite_C_mapping, regex("ko\\d{5}", ignore_case = TRUE))),
    brite_D_code = str_to_upper(str_extract(brite_D_mapping, regex("K\\d{5}", ignore_case = TRUE)))
  ) %>%
  mutate(across(c(Substrate, Substrate_category), ~replace_na(.x, "Not provided")))

# ---- Analyze one BRITE level ------------------------------------------------

analyze_level <- function(level) {
  level_dir <- file.path(figures_dir, level)
  dir.create(level_dir, showWarnings = FALSE, recursive = TRUE)
  code_col <- if (level == "brite_C") "brite_C_code" else "brite_D_code"
  mapping_text_col <- if (level == "brite_C") "brite_C_mapping" else "brite_D_mapping"

  summary_file <- file.path(all_gam_dir, paste0(level, "_GAM_summary_table.csv"))
  trajectory_file <- file.path(all_gam_dir, paste0(level, "_GAM_fitted_trajectories.csv"))
  data_file <- file.path(all_gam_dir, paste0(level, "_GAM_long_model_data.csv"))
  required <- c(summary_file, trajectory_file, data_file)
  if (any(!file.exists(required))) stop("Missing Script 09 output for ", level, ": ", paste(required[!file.exists(required)], collapse = ", "))

  map_level <- mapping %>%
    transmute(Substrate, Substrate_category, mapped_brite_class = .data[[mapping_text_col]], feature_code = .data[[code_col]]) %>%
    filter(!is.na(feature_code)) %>% distinct()

  summary <- read_csv(summary_file, show_col_types = FALSE) %>%
    mutate(feature_code = if (level == "brite_C") str_to_lower(str_extract(brite_class, "ko\\d{5}")) else str_to_upper(str_extract(brite_class, "^K\\d{5}")))
  trajectories <- read_csv(trajectory_file, show_col_types = FALSE) %>%
    mutate(Collection_date = as.Date(Collection_date), feature_code = if (level == "brite_C") str_to_lower(str_extract(brite_class, "ko\\d{5}")) else str_to_upper(str_extract(brite_class, "^K\\d{5}")))
  model_data <- read_csv(data_file, show_col_types = FALSE) %>%
    mutate(Collection_date = as.Date(Collection_date), feature_code = if (level == "brite_C") str_to_lower(str_extract(brite_class, "ko\\d{5}")) else str_to_upper(str_extract(brite_class, "^K\\d{5}")))

  target_features <- map_level %>% group_by(feature_code) %>%
    summarise(mapped_brite_class = first(mapped_brite_class), linked_substrates = paste(sort(unique(Substrate)), collapse = "; "),
              substrate_categories = paste(sort(unique(Substrate_category)), collapse = "; "), n_linked_substrates = n_distinct(Substrate), .groups = "drop") %>%
    left_join(summary, by = "feature_code") %>%
    mutate(p_value_target_set = if_else(fit_status == "Converged" & converged == "Yes", p_value, NA_real_),
           p_adj_fdr_microresp_targets = p.adjust(p_value_target_set, method = "BH"),
           significant_time_microresp_targets = if_else(!is.na(p_adj_fdr_microresp_targets) & p_adj_fdr_microresp_targets < 0.05, "Yes", "No")) %>%
    select(-p_value_target_set) %>% arrange(p_adj_fdr_microresp_targets, p_value, feature_code)

  mapping_audit <- map_level %>% left_join(target_features %>% select(feature_code, brite_class, fit_status, converged, p_value, p_adj_fdr, p_adj_fdr_microresp_targets), by = "feature_code") %>%
    mutate(mapping_status = case_when(is.na(brite_class) ~ "No match in filtered BRITE data", fit_status == "Converged" & converged == "Yes" ~ "Matched; model converged", TRUE ~ paste0("Matched; model status: ", fit_status)))
  write_csv(mapping_audit, file.path(results_dir, paste0(level, "_MicroResp_mapping_audit.csv")))
  write_csv(target_features, file.path(results_dir, paste0(level, "_MicroResp_GAM_statistics.csv")))

  eligible <- target_features %>% filter(fit_status == "Converged", converged == "Yes", !is.na(brite_class))
  eligible_codes <- eligible$feature_code
  target_trajectories <- trajectories %>% filter(feature_code %in% eligible_codes)
  target_points <- model_data %>% filter(feature_code %in% eligible_codes)

  if (n_distinct(target_trajectories$feature_code) < 3) stop("Fewer than three mapped, converged ", level, " features are available for clustering.")
  duplicate_values <- target_trajectories %>% count(feature_code, Collection_date) %>% filter(n > 1)
  if (nrow(duplicate_values) > 0) stop("Duplicate fitted values detected for ", level)

  dates <- target_trajectories %>% distinct(Collection_date, time_numeric) %>% arrange(time_numeric)
  fitted_matrix <- target_trajectories %>% mutate(date_id = paste0("day_", time_numeric)) %>%
    select(feature_code, date_id, fitted) %>% pivot_wider(names_from = date_id, values_from = fitted) %>%
    column_to_rownames("feature_code") %>% as.matrix()
  fitted_matrix <- fitted_matrix[, paste0("day_", dates$time_numeric), drop = FALSE]
  zero_variance <- rownames(fitted_matrix)[!is.finite(apply(fitted_matrix, 1, sd)) | apply(fitted_matrix, 1, sd) == 0]
  if (length(zero_variance) > 0) {
    fitted_matrix <- fitted_matrix[!rownames(fitted_matrix) %in% zero_variance, , drop = FALSE]
    eligible <- eligible %>% filter(!feature_code %in% zero_variance)
    target_trajectories <- target_trajectories %>% filter(!feature_code %in% zero_variance)
    target_points <- target_points %>% filter(!feature_code %in% zero_variance)
  }
  zmat <- row_zscore(fitted_matrix); zmat[zmat < -3] <- -3; zmat[zmat > 3] <- 3
  colnames(zmat) <- format(dates$Collection_date, "%b'%y")
  write_csv(as.data.frame(zmat) %>% rownames_to_column("feature_code"), file.path(results_dir, paste0(level, "_MicroResp_zscored_fitted_trajectories.csv")))

  correlation <- cor(t(zmat), use = "pairwise.complete.obs")
  distance_matrix <- 1 - correlation
  distance_matrix[distance_matrix < 0] <- 0
  diag(distance_matrix) <- 0
  distance <- as.dist(distance_matrix)
  
  hc <- hclust(distance, method = "ward.D2")
  kmax <- min(max_k_to_test, nrow(zmat) - 1L)
  elbow <- map_dfr(1:kmax, function(k) {
    ids <- cutree(hc, k); tibble(k = k, total_within_dissimilarity = total_within_dissimilarity(ids, distance))
  })
  automatic_k <- detect_elbow_k(elbow)
  selected_k <- if (is.na(manual_k[[level]])) automatic_k else manual_k[[level]]
  if (selected_k < 2 || selected_k > kmax) stop("Selected k is outside the valid range for ", level)
  elbow <- elbow %>% mutate(automatic_elbow_k = automatic_k, selected_k = selected_k, selected = k == selected_k)
  write_csv(elbow, file.path(results_dir, paste0(level, "_MicroResp_elbow_diagnostics.csv")))
  elbow_plot <- ggplot(elbow, aes(k, total_within_dissimilarity)) + geom_line() + geom_point() + geom_vline(xintercept = selected_k, linetype = "dashed", color = "#0072B2") +
    scale_x_continuous(breaks = 1:kmax) + labs(title = paste0("MicroResp-linked ", str_replace(level, "_", " "), " elbow diagnostic"), x = "Number of clusters", y = "Within-cluster dissimilarity") + theme(panel.grid.minor = element_blank())
  ggsave(file.path(level_dir, paste0(level, "_elbow_diagnostic.png")), elbow_plot, width = 8, height = 5.5, dpi = 300)

  raw_ids <- cutree(hc, selected_k)
  raw_order <- tibble(feature_code = names(raw_ids), raw_cluster = as.integer(raw_ids)) %>% mutate(peak = apply(zmat[feature_code, , drop = FALSE], 1, which.max)) %>%
    group_by(raw_cluster) %>% summarise(mean_peak = mean(peak), .groups = "drop") %>% arrange(mean_peak, raw_cluster) %>% mutate(Cluster_ID = paste0("Cluster_", row_number()))
  assignments <- tibble(feature_code = names(raw_ids), raw_cluster = as.integer(raw_ids)) %>% left_join(raw_order %>% select(raw_cluster, Cluster_ID), by = "raw_cluster") %>%
    left_join(eligible, by = "feature_code") %>%
    mutate(Temporal_significance = if_else(significant_time_microresp_targets == "Yes", "Significant", "Not significant"),
           display_label = paste(feature_code, short_label(brite_class, level)), dend_position = match(feature_code, hc$labels[hc$order])) %>% arrange(Cluster_ID, dend_position)
  write_csv(assignments %>% select(-dend_position), file.path(results_dir, paste0(level, "_MicroResp_cluster_assignments.csv")))
  write_csv(map_level %>% inner_join(assignments %>% select(feature_code, Cluster_ID, Temporal_significance, p_adj_fdr_microresp_targets, time_variance_explained), by = "feature_code"), file.path(results_dir, paste0(level, "_MicroResp_substrate_to_cluster.csv")))

  cluster_colours <- setNames(hcl.colors(selected_k, "Dark 3"), paste0("Cluster_", seq_len(selected_k)))
  dend <- as.dendrogram(hc) %>% set("branches_lwd", 0.9)
  if (nrow(assignments) <= 80) { labels(dend) <- setNames(assignments$display_label, assignments$feature_code)[labels(dend)]; dend <- dend %>% set("labels_cex", 0.5) } else dend <- dend %>% set("labels", rep("", length(labels(dend))))
  draw_dend <- function() { plot(dend, main = paste0("MicroResp-linked ", str_replace(level, "_", " "), " clustering"), ylab = "Ward.D2 height", leaflab = "perpendicular"); rect.hclust(hc, k = selected_k, border = cluster_colours) }
  png(file.path(level_dir, paste0(level, "_dendrogram.png")), width = 12, height = 7, units = "in", res = 300); draw_dend(); dev.off()
  pdf(file.path(level_dir, paste0(level, "_dendrogram.pdf")), width = 12, height = 7); draw_dend(); dev.off()

  trend_data <- as.data.frame(zmat) %>% rownames_to_column("feature_code") %>% pivot_longer(-feature_code, names_to = "date_label", values_to = "z_fitted") %>% left_join(assignments %>% select(feature_code, Cluster_ID), by = "feature_code") %>%
    mutate(date_label = factor(date_label, levels = colnames(zmat)), Cluster_ID = factor(Cluster_ID, levels = names(cluster_colours)))
  trend_summary <- trend_data %>% group_by(Cluster_ID, date_label) %>% summarise(mean_z = mean(z_fitted), se_z = if_else(n() > 1, sd(z_fitted) / sqrt(n()), 0), n_features = n(), .groups = "drop")
  write_csv(trend_summary, file.path(results_dir, paste0(level, "_MicroResp_cluster_average_trajectories.csv")))
  trend_plot <- ggplot(trend_summary, aes(date_label, mean_z, group = Cluster_ID, color = Cluster_ID, fill = Cluster_ID)) + geom_ribbon(aes(ymin = mean_z - se_z, ymax = mean_z + se_z), alpha = .18, color = NA) + geom_line(linewidth = 1) + geom_point(size = 2) +
    scale_color_manual(values = cluster_colours) + scale_fill_manual(values = cluster_colours) + labs(title = paste0("MicroResp-linked ", str_replace(level, "_", " "), " cluster-average trajectories"), subtitle = paste0("All mapped, converged features, k = ", selected_k), x = NULL, y = "Mean z-scored fitted expression") + theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())
  ggsave(file.path(level_dir, paste0(level, "_cluster_average_trajectories.png")), trend_plot, width = 9, height = 5.5, dpi = 300)
  ggsave(file.path(level_dir, paste0(level, "_cluster_average_trajectories.pdf")), trend_plot, width = 9, height = 5.5)

  # BRITE C is small enough to show in full. BRITE D gets a readable,
  # cluster-balanced heatmap, with the selection audit written to Results.
  if (level == "brite_C") {
    heatmap_assignments <- assignments
    heatmap_description <- "All mapped, converged BRITE C pathways"
  } else {
    heatmap_assignments <- assignments %>% group_by(Cluster_ID) %>% arrange(desc(time_variance_explained), p_adj_fdr_microresp_targets, .by_group = TRUE) %>% slice_head(n = top_brite_D_per_cluster) %>% ungroup()
    heatmap_description <- paste0("Top ", top_brite_D_per_cluster, " BRITE D KOs per cluster, ranked by time-only variance explained")
    write_csv(heatmap_assignments, file.path(results_dir, "brite_D_MicroResp_heatmap_feature_selection.csv"))
  }
  save_heatmap(zmat, heatmap_assignments, level, "MicroResp_heatmap", heatmap_description, level_dir)

  # Keep BRITE C GAM panels complete. For BRITE D, show the same balanced,
  # top-trajectory selection used in the readable heatmap. The full BRITE D
  # GAM statistics and cluster assignments are still exported above.
  plotted_codes <- heatmap_assignments$feature_code
  labeled_predictions <- target_trajectories %>%
    filter(feature_code %in% plotted_codes) %>%
    left_join(assignments %>% select(feature_code, display_label), by = "feature_code")
  labeled_points <- target_points %>%
    filter(feature_code %in% plotted_codes) %>%
    left_join(assignments %>% select(feature_code, display_label), by = "feature_code")
  gam_page_summary <- make_paginated_gam_plots(
    labeled_points %>% rename(plot_label = display_label),
    labeled_predictions %>% rename(plot_label = display_label),
    level,
    level_dir
  )
  write_csv(gam_page_summary, file.path(results_dir, paste0(level, "_MicroResp_GAM_page_summary.csv")))

  summary_row <- tibble(brite_level = level, n_mapping_rows = nrow(map_level), n_unique_requested_features = n_distinct(map_level$feature_code), n_matched_features = sum(!is.na(target_features$brite_class)), n_converged_features = nrow(eligible), n_clustered_features = nrow(zmat), n_target_FDR_significant = sum(assignments$Temporal_significance == "Significant"), automatic_elbow_k = automatic_k, selected_k = selected_k, heatmap_display_rule = heatmap_description)
  write_csv(summary_row, file.path(results_dir, paste0(level, "_MicroResp_run_summary.csv")))
  summary_row
}

# ---- Run --------------------------------------------------------------------

run_summary <- bind_rows(analyze_level("brite_C"), analyze_level("brite_D"))
write_csv(run_summary, file.path(results_dir, "MicroResp_BRITE_C_D_run_summary.csv"))
write_csv(tibble(mapping_file = mapping_file, all_pathway_gam_directory = all_gam_dir,
                 model = "Expression ~ s(time_numeric, k = 6) + s(Plot, bs = 're') fit in Script 09",
                 clustering = "Population GAM trajectories, row z-score, 1 - Pearson correlation, Ward.D2",
                 inclusion = "All mapped features with converged GAMs, regardless of target FDR significance",
                 BRITE_D_heatmap = "Top time-only trajectories within each cluster; display selection only"),
          file.path(results_dir, "run_manifest.csv"))

cat("\nMicroResp-targeted BRITE C/D GAM review and clustering complete.\n")
print(run_summary)
