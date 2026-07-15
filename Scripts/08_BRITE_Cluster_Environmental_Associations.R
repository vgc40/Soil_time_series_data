############################################################
# BRITE cluster associations with environmental parameters
#
# This script uses sample-level expression to calculate a score for each
# BRITE temporal cluster in every metatranscriptome sample. It then relates
# those cluster scores to environmental and biogeochemical measurements.
#
# Inputs:
#   Data/Data_merged_cleaned.csv
#   Results/BRITE_temporal_txt_only_significant_only/
#     brite_C_GAM_long_model_data.csv
#     brite_D_GAM_long_model_data.csv
#   Results/BRITE_clustering_elbow/brite_C/cluster_assignments.csv
#   Results/BRITE_clustering_elbow/brite_D/cluster_assignments.csv
#
# Main outputs for BRITE C and BRITE D separately:
#   1. Sample-level cluster scores
#   2. Spearman correlation tables and heatmaps
#   3. Random-intercept mixed-model association tables
#   4. Scatterplot PDFs for the strongest associations
#
# Important interpretation:
# These are descriptive associations. Both cluster scores and environmental
# variables change through time, so an association does not establish that an
# environmental variable caused the transcript trajectory.
############################################################

rm(list = ls())

# ---- User settings ----------------------------------------------------------

run_label <- "txt_only_significant_only"

# Default variables chosen to connect directly with the manuscript storyboard.
# The script skips any listed variable that is not present in the metadata.
environmental_variables <- c(
  "GWC",
  "pH",
  "SoilTemp_mean",
  "SEOC_ug_gdw",
  "SETN_ug_gdw",
  "SEOC_SETN_ratio",
  "MBC_ug_gdw_bgc",
  "MBN_ug_gdw_bgc",
  "Roots_gdw",
  "Resp_rate_mg_C_m2_h"
)

environmental_labels <- c(
  "GWC" = "Gravimetric water content",
  "pH" = "pH",
  "SoilTemp_mean" = "Mean soil temperature",
  "SEOC_ug_gdw" = "Salt-extractable organic C",
  "SETN_ug_gdw" = "Salt-extractable total N",
  "SEOC_SETN_ratio" = "SEOC:SETN",
  "MBC_ug_gdw_bgc" = "Microbial biomass C",
  "MBN_ug_gdw_bgc" = "Microbial biomass N",
  "Roots_gdw" = "Root biomass",
  "Resp_rate_mg_C_m2_h" = "Field respiration"
)

# Maximum number of associations shown in the scatterplot PDF for each level.
max_scatter_associations <- 12

# Number of plots on each PDF page.
plots_per_page <- 6

# ---- Packages and directories ----------------------------------------------

required_packages <- c(
  "tidyverse",
  "lme4",
  "lmerTest",
  "patchwork"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(patchwork)
})

theme_set(theme_bw(base_size = 11))

gam_results_dir <- file.path(
  "Results",
  paste0("BRITE_temporal_", run_label)
)

cluster_results_dir <- file.path("Results", "BRITE_clustering_elbow")

results_dir <- file.path(
  "Results",
  "BRITE_cluster_environment_associations"
)

figures_dir <- file.path(
  "Figures",
  "BRITE_cluster_environment_associations"
)

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Metadata ---------------------------------------------------------------

metadata <- read.csv(
  "Data/Data_merged_cleaned.csv",
  check.names = FALSE
) %>%
  mutate(
    Collection_date = as.Date(Collection_date),
    Plot = factor(Plot),
    SEOC_SETN_ratio = if_else(
      !is.na(SETN_ug_gdw) & SETN_ug_gdw != 0,
      SEOC_ug_gdw / SETN_ug_gdw,
      NA_real_
    ),
    Month_label = format(Collection_date, "%b'%y")
  )

available_environmental_variables <- environmental_variables[
  environmental_variables %in% names(metadata)
]

missing_environmental_variables <- setdiff(
  environmental_variables,
  available_environmental_variables
)

if (length(missing_environmental_variables) > 0) {
  warning(
    "Skipping environmental variables that are absent from metadata: ",
    paste(missing_environmental_variables, collapse = ", ")
  )
}

if (length(available_environmental_variables) == 0) {
  stop("None of the requested environmental variables were found.")
}

# ---- Helper functions -------------------------------------------------------

safe_zscore <- function(x) {
  x <- as.numeric(x)
  x_sd <- sd(x, na.rm = TRUE)

  if (is.na(x_sd) || x_sd == 0) {
    return(rep(NA_real_, length(x)))
  }

  (x - mean(x, na.rm = TRUE)) / x_sd
}

significance_symbol <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

calculate_spearman <- function(data, cluster_id, environmental_variable) {
  d <- data %>%
    filter(Cluster_ID == cluster_id) %>%
    select(cluster_score, all_of(environmental_variable)) %>%
    drop_na()

  if (
    nrow(d) < 5 ||
    n_distinct(d$cluster_score) < 3 ||
    n_distinct(d[[environmental_variable]]) < 3
  ) {
    return(tibble(
      Cluster_ID = cluster_id,
      environmental_variable = environmental_variable,
      n_complete = nrow(d),
      rho = NA_real_,
      p_value = NA_real_
    ))
  }

  test <- suppressWarnings(
    cor.test(
      d$cluster_score,
      d[[environmental_variable]],
      method = "spearman",
      exact = FALSE
    )
  )

  tibble(
    Cluster_ID = cluster_id,
    environmental_variable = environmental_variable,
    n_complete = nrow(d),
    rho = unname(test$estimate),
    p_value = test$p.value
  )
}

fit_environment_mixed_model <- function(
    data,
    cluster_id,
    environmental_variable) {

  d <- data %>%
    filter(Cluster_ID == cluster_id) %>%
    transmute(
      Plot = droplevels(factor(Plot)),
      Collection_date,
      score_z = safe_zscore(cluster_score),
      env_z = safe_zscore(.data[[environmental_variable]])
    ) %>%
    drop_na(score_z, env_z, Plot)

  empty_result <- tibble(
    Cluster_ID = cluster_id,
    environmental_variable = environmental_variable,
    n_complete = nrow(d),
    n_plots = n_distinct(d$Plot),
    estimate = NA_real_,
    std_error = NA_real_,
    df = NA_real_,
    t_value = NA_real_,
    p_value = NA_real_,
    singular_fit = NA,
    model_error = NA_character_
  )

  if (
    nrow(d) < 10 ||
    n_distinct(d$Plot) < 3 ||
    n_distinct(d$score_z) < 3 ||
    n_distinct(d$env_z) < 3
  ) {
    empty_result$model_error <- "Insufficient complete variation"
    return(empty_result)
  }

  model_error <- NA_character_

  model <- tryCatch(
    lmer(
      score_z ~ env_z + (1 | Plot),
      data = d,
      REML = TRUE
    ),
    error = function(e) {
      model_error <<- conditionMessage(e)
      NULL
    }
  )

  if (is.null(model)) {
    empty_result$model_error <- model_error
    return(empty_result)
  }

  coefficient_table <- coef(summary(model))

  if (!"env_z" %in% rownames(coefficient_table)) {
    empty_result$model_error <- "Environmental coefficient unavailable"
    return(empty_result)
  }

  coefficient <- coefficient_table["env_z", , drop = FALSE]

  tibble(
    Cluster_ID = cluster_id,
    environmental_variable = environmental_variable,
    n_complete = nrow(d),
    n_plots = n_distinct(d$Plot),
    estimate = coefficient[1, "Estimate"],
    std_error = coefficient[1, "Std. Error"],
    df = coefficient[1, "df"],
    t_value = coefficient[1, "t value"],
    p_value = coefficient[1, "Pr(>|t|)"],
    singular_fit = isSingular(model, tol = 1e-4),
    model_error = model_error
  )
}

make_scatterplot <- function(
    data,
    correlation_row,
    mixed_model_row,
    level_name) {

  cluster_id <- correlation_row$Cluster_ID
  environmental_variable <- correlation_row$environmental_variable

  d <- data %>%
    filter(Cluster_ID == cluster_id) %>%
    filter(
      !is.na(cluster_score),
      !is.na(.data[[environmental_variable]])
    )

  model_p <- if (nrow(mixed_model_row) == 1) {
    mixed_model_row$p_adj_fdr
  } else {
    NA_real_
  }

  annotation <- paste0(
    "Spearman rho = ", round(correlation_row$rho, 2),
    "\nSpearman FDR = ", signif(correlation_row$p_adj_fdr, 2),
    "\nMixed-model FDR = ",
    if_else(is.na(model_p), "NA", as.character(signif(model_p, 2)))
  )

  ggplot(
    d,
    aes(
      x = .data[[environmental_variable]],
      y = cluster_score,
      color = Collection_date
    )
  ) +
    geom_point(size = 2.3, alpha = 0.8) +
    geom_smooth(
      aes(group = 1),
      method = "lm",
      se = TRUE,
      color = "grey25",
      fill = "grey75",
      linewidth = 0.7
    ) +
    scale_color_viridis_c(option = "C", end = 0.9) +
    annotate(
      "label",
      x = -Inf,
      y = Inf,
      label = annotation,
      hjust = -0.05,
      vjust = 1.05,
      size = 3.2,
      label.size = 0.2,
      fill = "white"
    ) +
    labs(
      title = paste0(
        str_replace(level_name, "_", " "),
        " ", cluster_id
      ),
      x = unname(environmental_labels[environmental_variable]),
      y = "Sample-level cluster score",
      color = "Collection date"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "bottom"
    )
}

# ---- Main function ----------------------------------------------------------

analyze_brite_environment <- function(level_name) {
  message("\n=== ", level_name, " environmental associations ===")

  level_results_dir <- file.path(results_dir, level_name)
  level_figures_dir <- file.path(figures_dir, level_name)

  dir.create(level_results_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(level_figures_dir, showWarnings = FALSE, recursive = TRUE)

  model_data_file <- file.path(
    gam_results_dir,
    paste0(level_name, "_GAM_long_model_data.csv")
  )

  assignment_file <- file.path(
    cluster_results_dir,
    level_name,
    "cluster_assignments.csv"
  )

  if (!file.exists(model_data_file)) {
    stop("Missing model data: ", model_data_file)
  }

  if (!file.exists(assignment_file)) {
    stop("Missing cluster assignments: ", assignment_file)
  }

  model_data <- read_csv(model_data_file, show_col_types = FALSE)
  assignments <- read_csv(assignment_file, show_col_types = FALSE)

  if (!level_name %in% names(model_data)) {
    stop(
      "Expected feature column '", level_name,
      "' was not found in ", model_data_file
    )
  }

  # Standardize expression within each BRITE feature so features with large
  # absolute expression do not dominate the cluster score.
  expression_z <- model_data %>%
    rename(brite_class = all_of(level_name)) %>%
    semi_join(assignments, by = "brite_class") %>%
    mutate(
      Collection_date = as.Date(Collection_date),
      Plot = factor(Plot),
      Expression = as.numeric(Expression)
    ) %>%
    group_by(brite_class) %>%
    mutate(feature_expression_z = safe_zscore(Expression)) %>%
    ungroup() %>%
    drop_na(feature_expression_z) %>%
    left_join(
      assignments %>% select(brite_class, Cluster_ID),
      by = "brite_class"
    )

  cluster_scores <- expression_z %>%
    group_by(
      Sample_cat,
      Collection_date,
      Plot,
      Cluster_ID
    ) %>%
    summarise(
      cluster_score = mean(feature_expression_z, na.rm = TRUE),
      n_features_contributing = n_distinct(brite_class),
      .groups = "drop"
    ) %>%
    left_join(
      metadata %>%
        select(
          Sample_cat,
          Collection_date,
          Plot,
          Irrigation_trmt_pct,
          Irrigation_trmt_cat,
          Month_label,
          all_of(available_environmental_variables)
        ),
      by = c("Sample_cat", "Collection_date", "Plot")
    ) %>%
    mutate(
      brite_level = level_name,
      Cluster_ID = factor(
        Cluster_ID,
        levels = sort(unique(Cluster_ID))
      )
    )

  write_csv(
    cluster_scores,
    file.path(level_results_dir, "sample_level_cluster_scores.csv")
  )

  cluster_ids <- levels(cluster_scores$Cluster_ID)

  association_grid <- expand_grid(
    Cluster_ID = cluster_ids,
    environmental_variable = available_environmental_variables
  )

  spearman_results <- pmap_dfr(
    association_grid,
    function(Cluster_ID, environmental_variable) {
      calculate_spearman(
        cluster_scores,
        Cluster_ID,
        environmental_variable
      )
    }
  ) %>%
    mutate(
      brite_level = level_name,
      p_adj_fdr = p.adjust(p_value, method = "BH"),
      significance = significance_symbol(p_adj_fdr),
      environmental_label = unname(
        environmental_labels[environmental_variable]
      )
    ) %>%
    relocate(brite_level)

  write_csv(
    spearman_results,
    file.path(level_results_dir, "spearman_environment_correlations.csv")
  )

  mixed_model_results <- pmap_dfr(
    association_grid,
    function(Cluster_ID, environmental_variable) {
      fit_environment_mixed_model(
        cluster_scores,
        Cluster_ID,
        environmental_variable
      )
    }
  ) %>%
    mutate(
      brite_level = level_name,
      p_adj_fdr = p.adjust(p_value, method = "BH"),
      significance = significance_symbol(p_adj_fdr),
      environmental_label = unname(
        environmental_labels[environmental_variable]
      )
    ) %>%
    relocate(brite_level)

  write_csv(
    mixed_model_results,
    file.path(level_results_dir, "mixed_model_environment_associations.csv")
  )

  # ---- Correlation heatmap --------------------------------------------------

  heatmap_data <- spearman_results %>%
    mutate(
      environmental_label = factor(
        environmental_label,
        levels = rev(unname(
          environmental_labels[available_environmental_variables]
        ))
      ),
      Cluster_ID = factor(Cluster_ID, levels = cluster_ids),
      tile_label = if_else(
        is.na(rho),
        "NA",
        paste0(sprintf("%.2f", rho), significance)
      )
    )

  correlation_heatmap <- ggplot(
    heatmap_data,
    aes(Cluster_ID, environmental_label, fill = rho)
  ) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = tile_label), size = 3.5) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      na.value = "grey90"
    ) +
    labs(
      title = paste0(
        str_replace(level_name, "_", " "),
        " cluster associations with environmental parameters"
      ),
      subtitle = paste0(
        "Spearman correlations using sample-level cluster scores; ",
        "* FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001"
      ),
      x = NULL,
      y = NULL,
      fill = "Spearman\nrho"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(level_figures_dir, "spearman_correlation_heatmap.png"),
    correlation_heatmap,
    width = 8.5,
    height = 6.5,
    dpi = 300
  )

  ggsave(
    file.path(level_figures_dir, "spearman_correlation_heatmap.pdf"),
    correlation_heatmap,
    width = 8.5,
    height = 6.5
  )

  # ---- Mixed-model coefficient heatmap -------------------------------------

  model_heatmap_data <- mixed_model_results %>%
    mutate(
      environmental_label = factor(
        environmental_label,
        levels = rev(unname(
          environmental_labels[available_environmental_variables]
        ))
      ),
      Cluster_ID = factor(Cluster_ID, levels = cluster_ids),
      tile_label = if_else(
        is.na(estimate),
        "NA",
        paste0(sprintf("%.2f", estimate), significance)
      )
    )

  coefficient_limit <- max(
    abs(model_heatmap_data$estimate),
    na.rm = TRUE
  )

  if (!is.finite(coefficient_limit) || coefficient_limit == 0) {
    coefficient_limit <- 1
  }

  model_heatmap <- ggplot(
    model_heatmap_data,
    aes(Cluster_ID, environmental_label, fill = estimate)
  ) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = tile_label), size = 3.5) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-coefficient_limit, coefficient_limit),
      na.value = "grey90"
    ) +
    labs(
      title = paste0(
        str_replace(level_name, "_", " "),
        " environmental mixed-model coefficients"
      ),
      subtitle = paste0(
        "Standardized coefficient from score ~ environment + (1 | Plot); ",
        "* FDR < 0.05"
      ),
      x = NULL,
      y = NULL,
      fill = "Standardized\ncoefficient"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )

  ggsave(
    file.path(level_figures_dir, "mixed_model_coefficient_heatmap.png"),
    model_heatmap,
    width = 8.5,
    height = 6.5,
    dpi = 300
  )

  ggsave(
    file.path(level_figures_dir, "mixed_model_coefficient_heatmap.pdf"),
    model_heatmap,
    width = 8.5,
    height = 6.5
  )

  # ---- Scatterplots ---------------------------------------------------------

  selected_associations <- spearman_results %>%
    filter(!is.na(rho), !is.na(p_adj_fdr)) %>%
    arrange(p_adj_fdr, desc(abs(rho)))

  significant_associations <- selected_associations %>%
    filter(p_adj_fdr < 0.05)

  if (nrow(significant_associations) > 0) {
    selected_associations <- significant_associations
  }

  selected_associations <- selected_associations %>%
    slice_head(n = max_scatter_associations)

  write_csv(
    selected_associations,
    file.path(level_results_dir, "scatterplot_associations_selected.csv")
  )

  scatterplots <- vector("list", nrow(selected_associations))

  if (nrow(selected_associations) > 0) {
    for (i in seq_len(nrow(selected_associations))) {
      correlation_row <- selected_associations[i, ]

      model_row <- mixed_model_results %>%
        filter(
          Cluster_ID == correlation_row$Cluster_ID,
          environmental_variable ==
            correlation_row$environmental_variable
        )

      scatterplots[[i]] <- make_scatterplot(
        cluster_scores,
        correlation_row,
        model_row,
        level_name
      )
    }

    pdf(
      file.path(level_figures_dir, "top_environment_association_scatterplots.pdf"),
      width = 11,
      height = 8.5,
      onefile = TRUE
    )

    page_starts <- seq(1, length(scatterplots), by = plots_per_page)

    for (page_start in page_starts) {
      page_end <- min(
        page_start + plots_per_page - 1,
        length(scatterplots)
      )

      print(
        wrap_plots(
          scatterplots[page_start:page_end],
          ncol = 3,
          nrow = 2
        )
      )
    }

    dev.off()

    first_page_plots <- scatterplots[
      seq_len(min(plots_per_page, length(scatterplots)))
    ]

    first_page <- wrap_plots(
      first_page_plots,
      ncol = 3,
      nrow = 2
    )

    ggsave(
      file.path(level_figures_dir, "top_environment_associations_first_page.png"),
      first_page,
      width = 11,
      height = 8.5,
      dpi = 300
    )
  }

  # ---- Monthly trajectories for context ------------------------------------

  monthly_cluster_scores <- cluster_scores %>%
    group_by(Cluster_ID, Collection_date, Month_label) %>%
    summarise(
      mean_cluster_score = mean(cluster_score, na.rm = TRUE),
      se_cluster_score = sd(cluster_score, na.rm = TRUE) / sqrt(n()),
      n_samples = n(),
      .groups = "drop"
    )

  write_csv(
    monthly_cluster_scores,
    file.path(level_results_dir, "monthly_cluster_score_summary.csv")
  )

  monthly_plot <- ggplot(
    monthly_cluster_scores,
    aes(
      Collection_date,
      mean_cluster_score,
      color = Cluster_ID,
      fill = Cluster_ID
    )
  ) +
    geom_ribbon(
      aes(
        ymin = mean_cluster_score - se_cluster_score,
        ymax = mean_cluster_score + se_cluster_score
      ),
      alpha = 0.18,
      color = NA
    ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(
      title = paste0(
        str_replace(level_name, "_", " "),
        " sample-level cluster scores through time"
      ),
      x = NULL,
      y = "Mean sample-level cluster score",
      color = "Cluster",
      fill = "Cluster"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )

  ggsave(
    file.path(level_figures_dir, "sample_level_cluster_scores_through_time.png"),
    monthly_plot,
    width = 9,
    height = 5.5,
    dpi = 300
  )

  ggsave(
    file.path(level_figures_dir, "sample_level_cluster_scores_through_time.pdf"),
    monthly_plot,
    width = 9,
    height = 5.5
  )

  tibble(
    brite_level = level_name,
    n_samples = n_distinct(cluster_scores$Sample_cat),
    n_dates = n_distinct(cluster_scores$Collection_date),
    n_clusters = n_distinct(cluster_scores$Cluster_ID),
    n_environmental_variables = length(
      available_environmental_variables
    ),
    n_spearman_fdr_significant = sum(
      spearman_results$p_adj_fdr < 0.05,
      na.rm = TRUE
    ),
    n_mixed_model_fdr_significant = sum(
      mixed_model_results$p_adj_fdr < 0.05,
      na.rm = TRUE
    )
  )
}

# ---- Run BRITE C and BRITE D separately ------------------------------------

run_summary <- bind_rows(
  analyze_brite_environment("brite_C"),
  analyze_brite_environment("brite_D")
)

write_csv(
  run_summary,
  file.path(results_dir, "BRITE_environment_association_run_summary.csv")
)

write_csv(
  tibble(
    environmental_variable = available_environmental_variables,
    environmental_label = unname(
      environmental_labels[available_environmental_variables]
    )
  ),
  file.path(results_dir, "environmental_variables_used.csv")
)

cat("\nBRITE environmental association analysis complete.\n")
print(run_summary)
cat("\nResults:", results_dir, "\n")
cat("Figures:", figures_dir, "\n")
