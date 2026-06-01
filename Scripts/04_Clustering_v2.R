############################################################
# Step 3: Trajectory Clustering
# Correlation distance + Ward.D2 comparison version
#
# Objective:
# Cluster substrates by temporal shape of fitted GAM trajectories,
# using substrates with significant temporal trends from Step 2.
#
# This version:
# 1. Uses GAM-fitted trajectories, not raw observations.
# 2. Z-scores each substrate trajectory to compare shape, not magnitude.
# 3. Uses correlation distance: 1 - Pearson r.
# 4. Uses Ward.D2 linkage for comparison with the original dendrogram.
# 5. Evaluates k = 2 through k = 8 using silhouette width.
# 6. Adds an elbow-style diagnostic using within-cluster dissimilarity.
# 7. Exports and plots k = 2, k = 3, and the silhouette-selected best k.
# 8. Uses fixed cluster colors:
#      Cluster_1 = grey
#      Cluster_2 = green
# 9. Standardizes cluster labels by temporal direction to prevent color flipping.
############################################################

# -----------------------------
# 0. Setup
# -----------------------------

rm(list = ls())

required_packages <- c(
  "tidyverse",
  "cluster",
  "dendextend",
  "pheatmap",
  "RColorBrewer",
  "patchwork",
  "magick"
)

installed <- rownames(installed.packages())

for (pkg in required_packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

library(tidyverse)
library(cluster)
library(dendextend)
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(magick)

dir.create("Results", showWarnings = FALSE, recursive = TRUE)
dir.create("Figures", showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1. Load Step 2 outputs
# -----------------------------

gam_summary <- read_csv("Results/Step2_GAM_summary_table.csv", show_col_types = FALSE)
fitted_traj <- read_csv("Results/Step2_GAM_fitted_trajectories.csv", show_col_types = FALSE)

# Expected columns:
# gam_summary:
#   Substrate, Substrate_category, p_adj_FDR, Significant, Converged
#
# fitted_traj:
#   Substrate, time_numeric, fitted_SIR, SE

# -----------------------------
# 2. Identify significant and stable/excluded substrates
# -----------------------------

sig_substrates <- gam_summary %>%
  filter(
    Converged == "Yes",
    p_adj_FDR < 0.05
  ) %>%
  arrange(Substrate_category, Substrate)

stable_or_excluded_substrates <- gam_summary %>%
  filter(
    is.na(p_adj_FDR) | p_adj_FDR >= 0.05 | Converged != "Yes"
  ) %>%
  arrange(Substrate_category, Substrate)

write_csv(
  stable_or_excluded_substrates,
  "Results/Step3_temporally_stable_or_excluded_substrates_correlation_wardD2.csv"
)

if (nrow(sig_substrates) < 3) {
  stop("Fewer than 3 significant substrates available for clustering. Step 3 cannot proceed reliably.")
}

# -----------------------------
# 3. Create fitted trajectory matrix
# -----------------------------

traj_sig <- fitted_traj %>%
  semi_join(sig_substrates, by = "Substrate") %>%
  left_join(
    sig_substrates %>% select(Substrate, Substrate_category),
    by = "Substrate"
  ) %>%
  arrange(Substrate, time_numeric)

time_grid <- sort(unique(traj_sig$time_numeric))

traj_matrix <- traj_sig %>%
  select(Substrate, time_numeric, fitted_SIR) %>%
  mutate(time_numeric = paste0("time_", time_numeric)) %>%
  pivot_wider(
    names_from = time_numeric,
    values_from = fitted_SIR
  ) %>%
  column_to_rownames("Substrate") %>%
  as.matrix()

if (any(is.na(traj_matrix))) {
  stop("Trajectory matrix contains missing fitted values. Check Step2_GAM_fitted_trajectories.csv.")
}

# -----------------------------
# 4. Row-wise z-score normalization
# -----------------------------
# This removes magnitude differences and keeps temporal shape.

zero_var_rows <- rownames(traj_matrix)[apply(traj_matrix, 1, sd, na.rm = TRUE) == 0]

if (length(zero_var_rows) > 0) {
  warning(
    "Removing substrates with zero variance in fitted trajectory: ",
    paste(zero_var_rows, collapse = ", ")
  )
  
  traj_matrix <- traj_matrix[!rownames(traj_matrix) %in% zero_var_rows, , drop = FALSE]
  
  sig_substrates <- sig_substrates %>%
    filter(!Substrate %in% zero_var_rows)
}

zscore_matrix <- t(scale(t(traj_matrix)))

zscore_df <- as.data.frame(zscore_matrix) %>%
  rownames_to_column("Substrate") %>%
  left_join(
    sig_substrates %>% select(Substrate, Substrate_category),
    by = "Substrate"
  ) %>%
  relocate(Substrate, Substrate_category)

write_csv(
  zscore_df,
  "Results/Step3_zscore_trajectory_matrix_correlation_wardD2.csv"
)

# -----------------------------
# 5. Correlation distance and Ward.D2 clustering
# -----------------------------
# Distance = 1 - Pearson correlation.
# Ward.D2 is used here for comparison with the original version.
# Note: Ward.D2 is formally designed for Euclidean distances.

cor_matrix <- cor(
  t(zscore_matrix),
  method = "pearson",
  use = "pairwise.complete.obs"
)

dist_matrix <- as.dist(1 - cor_matrix)

hc <- hclust(
  dist_matrix,
  method = "ward.D2"
)

write_csv(
  as.data.frame(cor_matrix) %>% rownames_to_column("Substrate"),
  "Results/Step3_correlation_matrix_correlation_wardD2.csv"
)

write_csv(
  as.data.frame(as.matrix(dist_matrix)) %>% rownames_to_column("Substrate"),
  "Results/Step3_correlation_distance_matrix_correlation_wardD2.csv"
)

# -----------------------------
# 6. Helper functions
# -----------------------------

# Elbow-style total within-cluster dissimilarity.
# Uses the same correlation distance matrix used for clustering.
calc_total_within_dissimilarity <- function(cluster_ids, dist_matrix) {
  
  dist_full <- as.matrix(dist_matrix)
  
  sum(sapply(unique(cluster_ids), function(cl) {
    
    members <- names(cluster_ids)[cluster_ids == cl]
    
    if (length(members) <= 1) {
      return(0)
    }
    
    d_sub <- dist_full[members, members, drop = FALSE]
    
    # Sum of pairwise distances within cluster, adjusted to avoid double counting.
    sum(d_sub) / (2 * length(members))
  }))
}

# Standardize cluster labels by temporal direction.
# Raw cutree() labels can be arbitrary.
#
# Cluster_1 = lower late-minus-early tendency
# Cluster_2 = higher late-minus-early tendency
#
# For k > 2, clusters are ordered from most early-high/decreasing
# to most late-high/increasing.
standardize_cluster_labels <- function(cluster_ids, zscore_matrix) {
  
  cluster_ids <- as.integer(cluster_ids)
  names(cluster_ids) <- rownames(zscore_matrix)
  
  early_cols <- seq_len(max(1, floor(ncol(zscore_matrix) / 3)))
  late_cols <- seq(
    from = max(1, ceiling(2 * ncol(zscore_matrix) / 3)),
    to = ncol(zscore_matrix)
  )
  
  cluster_trend <- tibble(
    Substrate = names(cluster_ids),
    raw_cluster = cluster_ids,
    early_mean = rowMeans(zscore_matrix[, early_cols, drop = FALSE], na.rm = TRUE),
    late_mean = rowMeans(zscore_matrix[, late_cols, drop = FALSE], na.rm = TRUE),
    late_minus_early = late_mean - early_mean
  ) %>%
    group_by(raw_cluster) %>%
    summarise(
      cluster_late_minus_early = mean(late_minus_early, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(cluster_late_minus_early) %>%
    mutate(
      Cluster_ID = paste0("Cluster_", row_number())
    )
  
  relabeled <- tibble(
    Substrate = names(cluster_ids),
    raw_cluster = cluster_ids
  ) %>%
    left_join(cluster_trend, by = "raw_cluster") %>%
    select(Substrate, Cluster_ID)
  
  standardized_ids <- relabeled$Cluster_ID
  names(standardized_ids) <- relabeled$Substrate
  
  standardized_ids
}

# -----------------------------
# 7. Evaluate candidate cluster numbers
# -----------------------------

n_substrates <- nrow(zscore_matrix)

k_min <- 2
k_max <- min(8, n_substrates - 1)

if (k_max < k_min) {
  stop("Not enough substrates to evaluate k = 2 through k = 8.")
}

# Silhouette diagnostics, valid for k >= 2
silhouette_results <- tibble()

for (k in k_min:k_max) {
  
  cluster_tmp <- cutree(hc, k = k)
  
  sil_tmp <- silhouette(cluster_tmp, dist_matrix)
  
  silhouette_results <- bind_rows(
    silhouette_results,
    tibble(
      k = k,
      avg_silhouette_width = mean(sil_tmp[, "sil_width"]),
      min_silhouette_width = min(sil_tmp[, "sil_width"]),
      n_negative_silhouette = sum(sil_tmp[, "sil_width"] < 0)
    )
  )
}

best_k_silhouette <- silhouette_results %>%
  arrange(desc(avg_silhouette_width)) %>%
  slice(1) %>%
  pull(k)

silhouette_results <- silhouette_results %>%
  mutate(
    selected_by_silhouette = k == best_k_silhouette
  )

# Elbow diagnostics, include k = 1 for visual comparison
elbow_results <- tibble()

for (k in 1:k_max) {
  
  cluster_tmp <- cutree(hc, k = k)
  
  elbow_results <- bind_rows(
    elbow_results,
    tibble(
      k = k,
      total_within_dissimilarity = calc_total_within_dissimilarity(
        cluster_ids = cluster_tmp,
        dist_matrix = dist_matrix
      )
    )
  )
}

clustering_diagnostics <- elbow_results %>%
  left_join(silhouette_results, by = "k") %>%
  mutate(
    selected_k_reference = k == 2
  )

write_csv(
  clustering_diagnostics,
  "Results/Step3_clustering_diagnostics_correlation_wardD2.csv"
)

message("Best k by average silhouette width, correlation + Ward.D2: ", best_k_silhouette)

# -----------------------------
# 8. Diagnostic plots: elbow + silhouette
# -----------------------------

p_elbow <- ggplot(
  elbow_results,
  aes(x = k, y = total_within_dissimilarity)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 2, linetype = "dashed", linewidth = 0.5) +
  scale_x_continuous(breaks = 1:k_max) +
  labs(
    x = "Number of clusters (k)",
    y = "Total within-cluster dissimilarity",
    title = "Elbow diagnostic",
    subtitle = "Dashed line marks k = 2"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
  )

p_silhouette <- ggplot(
  silhouette_results,
  aes(x = k, y = avg_silhouette_width)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 2, linetype = "dashed", linewidth = 0.5) +
  scale_x_continuous(breaks = k_min:k_max) +
  labs(
    x = "Number of clusters (k)",
    y = "Average silhouette width",
    title = "Silhouette diagnostic",
    subtitle = paste0("Best k by silhouette = ", best_k_silhouette)
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
  )

p_diagnostics <- p_elbow / p_silhouette

ggsave(
  filename = "Figures/FigureS_Step3_elbow_and_silhouette_diagnostics_correlation_wardD2.png",
  plot = p_diagnostics,
  width = 8,
  height = 9,
  dpi = 300
)

ggsave(
  filename = "Figures/FigureS_Step3_elbow_and_silhouette_diagnostics_correlation_wardD2.pdf",
  plot = p_diagnostics,
  width = 8,
  height = 9
)

# -----------------------------
# 9. Export k = 2 and k = 3 cluster solutions
# -----------------------------

# Always export k = 2 and k = 3 for comparison.
# Also export the silhouette-selected best k, in case it is not 2 or 3.

candidate_ks <- sort(unique(c(2, 3, best_k_silhouette)))

cluster_solution_list <- list()

for (k in candidate_ks) {
  
  cluster_ids_raw <- cutree(hc, k = k)
  
  cluster_ids_standardized <- standardize_cluster_labels(
    cluster_ids = cluster_ids_raw,
    zscore_matrix = zscore_matrix
  )
  
  sil <- silhouette(cluster_ids_raw, dist_matrix)
  
  cluster_assignments <- tibble(
    Substrate = names(cluster_ids_raw),
    Raw_cluster_ID = paste0("Raw_cluster_", cluster_ids_raw),
    Cluster_ID = cluster_ids_standardized[names(cluster_ids_raw)],
    Silhouette_width = sil[, "sil_width"]
  ) %>%
    left_join(
      sig_substrates %>% select(Substrate, Substrate_category),
      by = "Substrate"
    ) %>%
    arrange(Cluster_ID, Substrate_category, Substrate) %>%
    select(
      Substrate,
      Substrate_category,
      Raw_cluster_ID,
      Cluster_ID,
      Silhouette_width
    )
  
  cluster_solution_list[[paste0("k", k)]] <- cluster_assignments
  
  write_csv(
    cluster_assignments,
    paste0("Results/Step3_cluster_assignments_k", k, "_correlation_wardD2.csv")
  )
  
  cluster_category_summary <- cluster_assignments %>%
    count(Cluster_ID, Substrate_category, name = "n_substrates") %>%
    group_by(Cluster_ID) %>%
    mutate(cluster_total = sum(n_substrates)) %>%
    ungroup() %>%
    arrange(Cluster_ID, desc(n_substrates), Substrate_category)
  
  write_csv(
    cluster_category_summary,
    paste0("Results/Step3_cluster_category_summary_k", k, "_correlation_wardD2.csv")
  )
}

# Compare k = 2 and k = 3 assignments, if both are available.

if (all(c("k2", "k3") %in% names(cluster_solution_list))) {
  
  cluster_compare <- cluster_solution_list$k2 %>%
    select(
      Substrate,
      Substrate_category,
      Cluster_k2 = Cluster_ID,
      Silhouette_k2 = Silhouette_width
    ) %>%
    left_join(
      cluster_solution_list$k3 %>%
        select(
          Substrate,
          Cluster_k3 = Cluster_ID,
          Silhouette_k3 = Silhouette_width
        ),
      by = "Substrate"
    ) %>%
    mutate(
      delta_silhouette_k3_minus_k2 = Silhouette_k3 - Silhouette_k2
    ) %>%
    arrange(Silhouette_k2)
  
  write_csv(
    cluster_compare,
    "Results/Step3_k2_vs_k3_cluster_comparison_correlation_wardD2.csv"
  )
}

# -----------------------------
# 10. Prepare labels and palettes
# -----------------------------

month_labels <- c(
  "Jun'23", "Jul'23", "Aug'23", "Sep'23", "Oct'23", "Nov'23",
  "Dec'23", "Jan'24", "Feb'24", "Mar'24", "Apr'24", "May'24", "Jun'24"
)

if (length(time_grid) != length(month_labels)) {
  month_labels <- paste0("T", seq_along(time_grid))
}

time_lookup <- tibble(
  time_col = colnames(zscore_matrix),
  time_numeric = time_grid,
  Month = month_labels
)

category_levels <- sig_substrates %>%
  arrange(Substrate_category) %>%
  pull(Substrate_category) %>%
  unique()

category_palette <- setNames(
  colorRampPalette(brewer.pal(8, "Dark2"))(length(category_levels)),
  category_levels
)

# Fixed named cluster palette.
# For k = 2:
#   Cluster_1 = grey
#   Cluster_2 = green
#
# Additional colors are only used for k > 2.
base_cluster_colors <- c(
  "Cluster_1" = "grey65",
  "Cluster_2" = "#1B9E77",
  "Cluster_3" = "#D95F02",
  "Cluster_4" = "#7570B3",
  "Cluster_5" = "#E7298A",
  "Cluster_6" = "#66A61E",
  "Cluster_7" = "#E6AB02",
  "Cluster_8" = "#A6761D"
)

# -----------------------------
# 11. Plotting function for k-specific figures
# -----------------------------

make_cluster_figures <- function(k, cluster_assignments) {
  
  cluster_levels <- paste0("Cluster_", seq_len(k))
  
  cluster_assignments <- cluster_assignments %>%
    mutate(
      Cluster_ID = factor(Cluster_ID, levels = cluster_levels)
    )
  
  cluster_palette <- base_cluster_colors[cluster_levels]
  
  # -----------------------------
  # Figure A: Dendrogram
  # -----------------------------
  
  dend <- as.dendrogram(hc)
  
  leaf_substrates <- labels(dend)
  
  leaf_categories <- cluster_assignments$Substrate_category[
    match(leaf_substrates, cluster_assignments$Substrate)
  ]
  
  dend <- dend %>%
    set("labels_cex", 0.75) %>%
    set("branches_lwd", 1.1) %>%
    set("labels_col", category_palette[leaf_categories])
  
  png(
    filename = paste0("Figures/Figure2A_dendrogram_k", k, "_correlation_wardD2.png"),
    width = 10,
    height = 7,
    units = "in",
    res = 300
  )
  
  plot(
    dend,
    main = paste0("Substrate trajectory clustering, k = ", k, ", correlation + Ward.D2"),
    ylab = "Ward.D2 height"
  )
  
  # Match dendrogram rectangle colors to actual left-to-right cluster order.
  # This prevents Cluster_1 and Cluster_2 colors from appearing flipped.
  dend_order <- hc$labels[hc$order]
  
  cluster_order_left_to_right <- cluster_assignments %>%
    mutate(
      dend_position = match(Substrate, dend_order)
    ) %>%
    group_by(Cluster_ID) %>%
    summarise(
      first_dend_position = min(dend_position, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(first_dend_position) %>%
    pull(Cluster_ID) %>%
    as.character()
  
  rect.hclust(
    hc,
    k = k,
    border = cluster_palette[cluster_order_left_to_right]
  )
  
  legend(
    "topright",
    legend = names(category_palette),
    fill = category_palette,
    title = "Substrate category",
    cex = 0.75,
    bty = "n"
  )
  
  dev.off()
  
  # -----------------------------
  # Figure B: Cluster-averaged trajectories
  # -----------------------------
  
  zscore_long <- as.data.frame(zscore_matrix) %>%
    rownames_to_column("Substrate") %>%
    pivot_longer(
      cols = -Substrate,
      names_to = "time_col",
      values_to = "zscore_fitted_SIR"
    ) %>%
    left_join(cluster_assignments, by = "Substrate") %>%
    left_join(time_lookup, by = "time_col") %>%
    mutate(
      Month = factor(Month, levels = month_labels),
      Cluster_ID = factor(Cluster_ID, levels = cluster_levels)
    )
  
  cluster_traj_summary <- zscore_long %>%
    group_by(Cluster_ID, time_numeric, Month) %>%
    summarise(
      mean_z = mean(zscore_fitted_SIR, na.rm = TRUE),
      se_z = sd(zscore_fitted_SIR, na.rm = TRUE) / sqrt(n()),
      n_substrates = n(),
      .groups = "drop"
    )
  
  p_cluster_traj <- ggplot(
    cluster_traj_summary,
    aes(
      x = Month,
      y = mean_z,
      group = Cluster_ID,
      color = Cluster_ID,
      fill = Cluster_ID
    )
  ) +
    geom_ribbon(
      aes(ymin = mean_z - se_z, ymax = mean_z + se_z),
      alpha = 0.20,
      color = NA
    ) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2) +
    scale_color_manual(values = cluster_palette, drop = FALSE) +
    scale_fill_manual(values = cluster_palette, drop = FALSE) +
    labs(
      x = NULL,
      y = "Mean z-scored fitted SIR",
      color = "Cluster",
      fill = "Cluster"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  
  ggsave(
    filename = paste0("Figures/Figure2B_cluster_trajectories_k", k, "_correlation_wardD2.png"),
    plot = p_cluster_traj,
    width = 9,
    height = 5.5,
    dpi = 300
  )
  
  # -----------------------------
  # Figure C: Heatmap
  # -----------------------------
  
  dend_order <- hc$labels[hc$order]
  
  ordered_substrates <- cluster_assignments %>%
    mutate(dend_order = match(Substrate, dend_order)) %>%
    arrange(Cluster_ID, dend_order) %>%
    pull(Substrate)
  
  heatmap_matrix <- zscore_matrix[ordered_substrates, , drop = FALSE]
  colnames(heatmap_matrix) <- month_labels
  
  row_annotation <- cluster_assignments %>%
    select(Substrate, Cluster_ID, Substrate_category) %>%
    mutate(
      Cluster_ID = factor(Cluster_ID, levels = cluster_levels),
      Substrate_category = factor(Substrate_category, levels = names(category_palette))
    ) %>%
    column_to_rownames("Substrate")
  
  row_annotation <- row_annotation[ordered_substrates, , drop = FALSE]
  
  annotation_colors <- list(
    Cluster_ID = cluster_palette,
    Substrate_category = category_palette
  )
  
  pheatmap(
    mat = heatmap_matrix,
    filename = paste0("Figures/Figure2C_heatmap_k", k, "_correlation_wardD2.png"),
    width = 9,
    height = 8,
    scale = "none",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    annotation_row = row_annotation,
    annotation_colors = annotation_colors,
    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
    fontsize_row = 8,
    fontsize_col = 10,
    main = paste0("Z-scored fitted MicroResp trajectories, k = ", k, ", correlation + Ward.D2")
  )
  
  # -----------------------------
  # Combined stacked figure
  # -----------------------------
  
  fig_a <- image_read(paste0("Figures/Figure2A_dendrogram_k", k, "_correlation_wardD2.png"))
  fig_b <- image_read(paste0("Figures/Figure2B_cluster_trajectories_k", k, "_correlation_wardD2.png"))
  fig_c <- image_read(paste0("Figures/Figure2C_heatmap_k", k, "_correlation_wardD2.png"))
  
  fig_a <- image_resize(fig_a, "2400x")
  fig_b <- image_resize(fig_b, "2400x")
  fig_c <- image_resize(fig_c, "2400x")
  
  combined_png <- image_append(c(fig_a, fig_b, fig_c), stack = TRUE)
  
  image_write(
    combined_png,
    path = paste0("Figures/Figure2_combined_k", k, "_correlation_wardD2.png"),
    format = "png"
  )
  
  image_write(
    combined_png,
    path = paste0("Figures/Figure2_combined_k", k, "_correlation_wardD2.pdf"),
    format = "pdf"
  )
}

# Generate figures for candidate k values.

for (k in candidate_ks) {
  make_cluster_figures(
    k = k,
    cluster_assignments = cluster_solution_list[[paste0("k", k)]]
  )
}

# -----------------------------
# 12. Select final k from diagnostics
# -----------------------------
# Use k = 2 as the primary solution, while reporting silhouette and elbow diagnostics.
# Change this to best_k_silhouette only if you want automatic k selection.

final_k <- 2

final_cluster_assignments <- cluster_solution_list[[paste0("k", final_k)]]

write_csv(
  final_cluster_assignments,
  "Results/Step3_cluster_assignments_FINAL_correlation_wardD2.csv"
)

final_cluster_category_summary <- final_cluster_assignments %>%
  count(Cluster_ID, Substrate_category, name = "n_substrates") %>%
  group_by(Cluster_ID) %>%
  mutate(cluster_total = sum(n_substrates)) %>%
  ungroup() %>%
  arrange(Cluster_ID, desc(n_substrates), Substrate_category)

write_csv(
  final_cluster_category_summary,
  "Results/Step3_cluster_category_summary_FINAL_correlation_wardD2.csv"
)

message("Final k selected for primary output: ", final_k)

# -----------------------------
# 13. Console summary
# -----------------------------

cat("\nStep 3 complete: correlation distance + Ward.D2 version.\n")
cat("Significant substrates clustered:", nrow(sig_substrates), "\n")
cat("Temporally stable or excluded substrates:", nrow(stable_or_excluded_substrates), "\n")
cat("Best k by average silhouette width:", best_k_silhouette, "\n")
cat("Primary final k used for output:", final_k, "\n")
cat("Elbow and silhouette diagnostics written to:\n")
cat("Figures/FigureS_Step3_elbow_and_silhouette_diagnostics_correlation_wardD2.png\n")
cat("Figures/FigureS_Step3_elbow_and_silhouette_diagnostics_correlation_wardD2.pdf\n\n")

cat("Clustering diagnostics:\n")
print(clustering_diagnostics)

cat("\nFinal cluster sizes:\n")
print(final_cluster_assignments %>% count(Cluster_ID, name = "n_substrates"))

if ("k2" %in% names(cluster_solution_list)) {
  cat("\nCluster sizes, k = 2:\n")
  print(cluster_solution_list$k2 %>% count(Cluster_ID, name = "n_substrates"))
}

if ("k3" %in% names(cluster_solution_list)) {
  cat("\nCluster sizes, k = 3:\n")
  print(cluster_solution_list$k3 %>% count(Cluster_ID, name = "n_substrates"))
}

if (exists("cluster_compare")) {
  cat("\nSubstrates with weakest k = 2 silhouette values:\n")
  print(cluster_compare %>% arrange(Silhouette_k2) %>% slice(1:min(10, n())))
  
  cat("\nSubstrates with weakest k = 3 silhouette values:\n")
  print(cluster_compare %>% arrange(Silhouette_k3) %>% slice(1:min(10, n())))
}

cat("\nOutputs written to Results/ and Figures/.\n")