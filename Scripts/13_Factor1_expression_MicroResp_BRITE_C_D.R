############################################################
# Factor 1 gene expression, BRITE summaries, random-plot GAMs,
# MicroResp overlap review, and trajectory clustering
#
# This script:
#   1. filters the original expression matrix to Factor 1 genes;
#   2. summarizes Factor 1 genes across BRITE C and BRITE D;
#   3. fits random-plot temporal GAMs for Factor 1 BRITE C and D;
#   4. compares Factor 1 BRITE C and D features against the MicroResp mapping;
#   5. clusters the MicroResp-linked Factor 1 BRITE trajectories and
#      saves heatmaps, dendrograms, and cluster-average trajectories.
############################################################

rm(list = ls())

required_packages <- c(
  "mgcv",
  "tidyverse",
  "readr",
  "dendextend",
  "pheatmap",
  "RColorBrewer",
  "ggforce"
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
  library(mgcv)
  library(tidyverse)
  library(readr)
  library(dendextend)
  library(pheatmap)
  library(RColorBrewer)
  library(ggforce)
})

theme_set(theme_bw(base_size = 11))

# ---- User settings ----------------------------------------------------------

expr_file <- paste0(
  "C:/Users/gara009/OneDrive - PNNL/",
  "Soil Microbiome SFA - A3.1 Rhizosphere and Moisture Effects/",
  "Prosser_Monthly/Manuscripts/A3-1-1_Prosser_time_biogeochem/",
  "Transcriptomic_abundance/geTMM_normalized_expression_noZeroes.csv"
)

annotation_file <- paste0(
  "C:/Users/gara009/OneDrive - PNNL/",
  "Soil Microbiome SFA - A3.1 Rhizosphere and Moisture Effects/",
  "Prosser_Monthly/Manuscripts/A3-1-1_Prosser_time_biogeochem/",
  "Annotations/2.5kbSizeOut.emapper.annotations.kegg_annotated.csv"
)

factor1_file <- file.path("Data", "factor1_all_omic_genes.csv")
meta_file <- file.path("Data", "Data_merged_cleaned.csv")
microresp_mapping_file <- file.path("Data", "MicroResp_all_BRITE_C_D_complete.csv")
remove_text_file <- file.path("Data", "pathwaystoremove", "kegg_human_exclude_2026-07-13.txt")

results_dir <- file.path("Results", "Factor1_BRITE_MicroResp")
figures_dir <- file.path("Figures", "Factor1_BRITE_MicroResp")
gam_results_dir <- file.path(results_dir, "temporal_gam")
gam_figures_dir <- file.path(figures_dir, "temporal_gam")
clustering_results_dir <- file.path(results_dir, "clustering")
clustering_figures_dir <- file.path(figures_dir, "clustering")
diagnostics_dir <- file.path(results_dir, "diagnostics")

for (directory in c(
  results_dir,
  figures_dir,
  gam_results_dir,
  gam_figures_dir,
  clustering_results_dir,
  clustering_figures_dir,
  diagnostics_dir
)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

time_basis_k <- 6L
accepted_fit_status <- "Converged"
max_k_to_test <- 10L
manual_k <- c(brite_C = NA_integer_, brite_D = NA_integer_)
top_brite_D_per_cluster <- 15L
gam_panels_per_page <- 6L

# ---- Read source data -------------------------------------------------------

factor1_genes <- read_csv(factor1_file, show_col_types = FALSE) %>%
  transmute(
    Query = as.character(feature),
    factor = as.character(factor),
    loading_value = as.numeric(value),
    abs_loading_value = as.numeric(abs_value),
    view = as.character(view)
  ) %>%
  distinct(Query, .keep_all = TRUE)

expr <- read.csv(expr_file, check.names = FALSE)
names(expr)[1] <- "Query"

annotation <- read.csv(annotation_file, check.names = FALSE) %>%
  mutate(
    Query = as.character(Query),
    brite_C = as.character(brite_C),
    brite_D = as.character(brite_D)
  )

meta <- read.csv(meta_file, check.names = FALSE) %>%
  mutate(
    Collection_date = as.Date(Collection_date),
    Plot = factor(Plot),
    time_numeric = as.numeric(Collection_date - min(Collection_date, na.rm = TRUE))
  ) %>%
  select(
    Sample_cat,
    Collection_date,
    Plot,
    Irrigation_trmt_pct,
    Irrigation_trmt_cat,
    time_numeric
  )

factor1_expr <- expr %>%
  mutate(Query = as.character(Query)) %>%
  semi_join(factor1_genes, by = "Query")

sample_cols <- setdiff(names(factor1_expr), "Query")

factor1_annotated <- factor1_genes %>%
  left_join(annotation, by = "Query") %>%
  left_join(factor1_expr, by = "Query")

write_csv(factor1_annotated, file.path(results_dir, "factor1_gene_annotation_expression_table.csv"))

factor1_gene_expression <- factor1_annotated %>%
  select(Query, all_of(sample_cols)) %>%
  group_by(Query) %>%
  summarise(across(all_of(sample_cols), ~sum(.x, na.rm = TRUE)), .groups = "drop")

factor1_brite_c_expression <- factor1_annotated %>%
  filter(!is.na(brite_C), brite_C != "") %>%
  group_by(brite_C) %>%
  summarise(across(all_of(sample_cols), ~sum(.x, na.rm = TRUE)), .groups = "drop")

factor1_brite_d_expression <- factor1_annotated %>%
  filter(!is.na(brite_D), brite_D != "") %>%
  group_by(brite_D) %>%
  summarise(across(all_of(sample_cols), ~sum(.x, na.rm = TRUE)), .groups = "drop")

write_csv(factor1_gene_expression, file.path(results_dir, "factor1_gene_expression.csv"))
write_csv(factor1_brite_c_expression, file.path(results_dir, "factor1_brite_C_expression.csv"))
write_csv(factor1_brite_d_expression, file.path(results_dir, "factor1_brite_D_expression.csv"))

# ---- Gene/pathway summaries and plots --------------------------------------

factor1_gene_summary <- factor1_annotated %>%
  summarise(
    n_factor1_genes = n_distinct(Query),
    n_factor1_genes_with_expression = sum(Query %in% factor1_expr$Query),
    n_factor1_genes_with_annotation = sum(!is.na(brite_C) | !is.na(brite_D)),
    n_unique_brite_C = n_distinct(brite_C[!is.na(brite_C) & brite_C != ""]),
    n_unique_brite_D = n_distinct(brite_D[!is.na(brite_D) & brite_D != ""])
  )

write_csv(factor1_gene_summary, file.path(results_dir, "factor1_gene_summary.csv"))

factor1_brite_c_summary <- factor1_annotated %>%
  filter(!is.na(brite_C), brite_C != "") %>%
  group_by(brite_C) %>%
  summarise(
    n_factor1_genes = n_distinct(Query),
    mean_abs_loading = mean(abs_loading_value, na.rm = TRUE),
    max_abs_loading = max(abs_loading_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_factor1_genes), desc(mean_abs_loading), brite_C)

factor1_brite_d_summary <- factor1_annotated %>%
  filter(!is.na(brite_D), brite_D != "") %>%
  mutate(brite_D_code = str_extract(brite_D, "^K\\d{5}")) %>%
  group_by(brite_D, brite_D_code) %>%
  summarise(
    n_factor1_genes = n_distinct(Query),
    mean_abs_loading = mean(abs_loading_value, na.rm = TRUE),
    max_abs_loading = max(abs_loading_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_factor1_genes), desc(mean_abs_loading), brite_D)

write_csv(factor1_brite_c_summary, file.path(results_dir, "factor1_brite_C_summary.csv"))
write_csv(factor1_brite_d_summary, file.path(results_dir, "factor1_brite_D_summary.csv"))

brite_c_plot_data <- factor1_brite_c_summary %>%
  slice_max(order_by = n_factor1_genes, n = 20)

if (nrow(brite_c_plot_data) > 0) {
  brite_c_plot <- ggplot(
    brite_c_plot_data,
    aes(x = reorder(brite_C, n_factor1_genes), y = n_factor1_genes)
  ) +
    geom_col(fill = "#2b8cbe") +
    coord_flip() +
    labs(
      title = "Factor 1 genes across BRITE C pathways",
      subtitle = "Top BRITE C pathways by number of Factor 1 genes",
      x = NULL,
      y = "Number of Factor 1 genes"
    ) +
    theme(panel.grid.minor = element_blank())

  ggsave(
    file.path(figures_dir, "factor1_brite_C_gene_counts.png"),
    brite_c_plot,
    width = 10,
    height = 7,
    dpi = 300
  )
}

brite_d_plot_data <- factor1_brite_d_summary %>%
  slice_max(order_by = n_factor1_genes, n = 30)

if (nrow(brite_d_plot_data) > 0) {
  brite_d_plot <- ggplot(
    brite_d_plot_data,
    aes(x = reorder(brite_D, n_factor1_genes), y = n_factor1_genes)
  ) +
    geom_col(fill = "#d95f0e") +
    coord_flip() +
    labs(
      title = "Factor 1 genes across BRITE D features",
      subtitle = "Top BRITE D features by number of Factor 1 genes",
      x = NULL,
      y = "Number of Factor 1 genes"
    ) +
    theme(panel.grid.minor = element_blank(), axis.text.y = element_text(size = 7))

  ggsave(
    file.path(figures_dir, "factor1_brite_D_gene_counts.png"),
    brite_d_plot,
    width = 11,
    height = 9,
    dpi = 300
  )
}

# ---- Remove excluded human pathways ----------------------------------------

human_remove_table <- read.delim(remove_text_file, check.names = FALSE) %>%
  mutate(
    path_ko = str_extract(as.character(path_ko), "ko\\d{5}"),
    KO_entry = str_extract(as.character(KO_entry), "K\\d{5}")
  )

brite_c_remove_codes <- human_remove_table %>%
  filter(!is.na(path_ko)) %>%
  distinct(path_ko)

brite_d_remove_codes <- human_remove_table %>%
  filter(!is.na(KO_entry)) %>%
  distinct(KO_entry)

factor1_brite_d_model_input <- factor1_brite_d_expression %>%
  mutate(
    brite_D = as.character(brite_D),
    ko_code = str_extract(brite_D, "^K\\d{5}"),
    remove_from_text = ko_code %in% brite_d_remove_codes$KO_entry
  )

factor1_brite_c_model_input <- factor1_brite_c_expression %>%
  mutate(
    brite_C = as.character(brite_C),
    pathway_code = str_extract(brite_C, "ko\\d{5}"),
    remove_from_text = pathway_code %in% brite_c_remove_codes$path_ko
  )

factor1_exclusion_summary <- bind_rows(
  factor1_brite_d_model_input %>%
    summarise(
      brite_level = "brite_D",
      starting_rows = n(),
      removed_from_text = sum(remove_from_text, na.rm = TRUE),
      remaining_rows = sum(!remove_from_text, na.rm = TRUE)
    ),
  factor1_brite_c_model_input %>%
    summarise(
      brite_level = "brite_C",
      starting_rows = n(),
      removed_from_text = sum(remove_from_text, na.rm = TRUE),
      remaining_rows = sum(!remove_from_text, na.rm = TRUE)
    )
) %>%
  mutate(
    exclusion_source = basename(remove_text_file),
    model = paste0(
      "Expression ~ s(time_numeric, k = ",
      time_basis_k,
      ") + s(Plot, bs = 're')"
    )
  )

write_csv(factor1_exclusion_summary, file.path(results_dir, "factor1_exclusion_summary.csv"))

factor1_brite_d_model_input <- factor1_brite_d_model_input %>%
  filter(!remove_from_text) %>%
  select(-ko_code, -remove_from_text)

factor1_brite_c_model_input <- factor1_brite_c_model_input %>%
  filter(!remove_from_text) %>%
  select(-pathway_code, -remove_from_text)

# ---- Small helper functions used repeatedly --------------------------------

safe_number <- function(x) {
  if (length(x) == 1 && !is.null(x) && is.finite(as.numeric(x))) {
    return(as.numeric(x))
  }
  NA_real_
}

collapse_messages <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) {
    return(NA_character_)
  }
  paste(x, collapse = " | ")
}

make_summary_row <- function(
    level_name,
    brite_category,
    data,
    fit_status,
    model = NULL,
    gam_summary = NULL,
    time_row = NULL,
    plot_row = NULL,
    time_variance_explained = NA_real_,
    warning_msg = NA_character_,
    error_msg = NA_character_) {
  tibble(
    brite_level = level_name,
    brite_class = brite_category,
    n_obs = nrow(data),
    n_plots = n_distinct(data$Plot),
    n_timepoints = n_distinct(data$time_numeric),
    n_unique_expression = n_distinct(data$Expression),
    expression_variance = if (nrow(data) > 1) var(data$Expression, na.rm = TRUE) else NA_real_,
    edf = safe_number(time_row$edf),
    ref_df = safe_number(time_row$Ref.df),
    F_statistic = safe_number(time_row$F),
    p_value = safe_number(time_row[["p-value"]]),
    plot_random_effect_edf = safe_number(plot_row$edf),
    adj_r_squared = safe_number(gam_summary$r.sq),
    adj_r_squared_full = safe_number(gam_summary$r.sq),
    deviance_explained = safe_number(gam_summary$dev.expl),
    deviance_explained_full = safe_number(gam_summary$dev.expl),
    time_variance_explained = safe_number(time_variance_explained),
    model_deviance = if (!is.null(model)) safe_number(deviance(model)) else NA_real_,
    null_deviance = if (!is.null(model)) safe_number(model$null.deviance) else NA_real_,
    scale = safe_number(gam_summary$scale),
    AIC = if (!is.null(model)) safe_number(AIC(model)) else NA_real_,
    fit_status = fit_status,
    converged = if_else(fit_status == accepted_fit_status, "Yes", "No"),
    model_warning = warning_msg,
    model_error = error_msg
  )
}

row_zscore <- function(x) {
  row_means <- rowMeans(x, na.rm = TRUE)
  row_sds <- apply(x, 1, sd, na.rm = TRUE)
  sweep(sweep(x, 1, row_means, "-"), 1, row_sds, "/")
}

total_within_dissimilarity <- function(cluster_ids, distance_matrix) {
  distance_full <- as.matrix(distance_matrix)
  sum(vapply(unique(cluster_ids), function(cluster_id) {
    members <- names(cluster_ids)[cluster_ids == cluster_id]
    if (length(members) <= 1) {
      return(0)
    }
    sum(distance_full[members, members, drop = FALSE]) / (2 * length(members))
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
  interior <- seq_len(nrow(elbow_table))[-c(1, nrow(elbow_table))]
  elbow_table$k[interior[which.max(abs(y[interior] - (1 - x[interior])) / sqrt(2))]]
}

short_label <- function(x, level) {
  if (level == "brite_C") {
    x %>%
      str_remove("^\\d{5}\\s+") %>%
      str_remove("\\s+\\[PATH:ko\\d{5}\\]$")
  } else {
    code <- str_extract(x, "^K\\d{5}")
    text <- x %>%
      str_remove("^K\\d{5}\\s+") %>%
      str_extract("^[^;]+") %>%
      str_squish()
    if_else(is.na(text) | text == "", code, paste(code, text))
  }
}

save_heatmap <- function(mat, assignments, level, suffix, title, output_dir) {
  if (nrow(assignments) == 0) {
    return(invisible(NULL))
  }
  plot_mat <- mat[assignments$feature_code, , drop = FALSE]
  rownames(plot_mat) <- make.unique(assignments$display_label)
  annotation <- assignments %>%
    transmute(display_label = rownames(plot_mat), Cluster_ID, Temporal_significance) %>%
    column_to_rownames("display_label")
  cluster_colors <- setNames(
    hcl.colors(n_distinct(assignments$Cluster_ID), "Dark 3"),
    sort(unique(assignments$Cluster_ID))
  )
  annotation_colors <- list(
    Cluster_ID = cluster_colors,
    Temporal_significance = c(Significant = "#238b45", `Not significant` = "#bdbdbd")
  )
  for (extension in c("png", "pdf")) {
    pheatmap(
      plot_mat,
      filename = file.path(output_dir, paste0(level, "_", suffix, ".", extension)),
      width = 11,
      height = max(7, 0.27 * nrow(assignments)),
      scale = "none",
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      annotation_row = annotation,
      annotation_colors = annotation_colors,
      color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
      breaks = seq(-3, 3, length.out = 101),
      fontsize_row = if (nrow(assignments) > 80) 5 else 7,
      border_color = NA,
      main = title
    )
  }
}

# ---- Fit Factor 1 BRITE C and D GAMs ---------------------------------------

fit_brite_level <- function(brite_table, brite_col, level_name) {
  sample_cols <- setdiff(names(brite_table), brite_col)

  model_data <- brite_table %>%
    filter(!is.na(.data[[brite_col]]), .data[[brite_col]] != "") %>%
    pivot_longer(
      all_of(sample_cols),
      names_to = "Sample_mt",
      values_to = "Expression"
    ) %>%
    mutate(Sample_cat = str_remove(Sample_mt, "_MT$")) %>%
    left_join(meta, by = "Sample_cat") %>%
    mutate(
      !!brite_col := as.character(.data[[brite_col]]),
      Expression = as.numeric(Expression),
      Plot = droplevels(factor(Plot))
    ) %>%
    drop_na(Collection_date, Plot, time_numeric, Expression)

  model_data %>%
    summarise(
      n_rows = n(),
      n_brite_categories = n_distinct(.data[[brite_col]]),
      n_samples = n_distinct(Sample_cat),
      n_plots = n_distinct(Plot),
      n_timepoints = n_distinct(time_numeric),
      min_date = min(Collection_date),
      max_date = max(Collection_date)
    ) %>%
    write_csv(file.path(diagnostics_dir, paste0(level_name, "_input_summary.csv")))

  time_grid <- model_data %>%
    distinct(time_numeric, Collection_date) %>%
    arrange(time_numeric, Collection_date)

  brite_categories <- unique(model_data[[brite_col]])
  summary_rows <- vector("list", length(brite_categories))
  prediction_rows <- vector("list", length(brite_categories))
  names(summary_rows) <- brite_categories
  names(prediction_rows) <- brite_categories

  for (brite_category in brite_categories) {
    d <- model_data %>%
      filter(.data[[brite_col]] == brite_category) %>%
      arrange(time_numeric, Plot, Sample_cat) %>%
      mutate(Plot = droplevels(factor(Plot)))

    if (n_distinct(d$time_numeric) < 4 || n_distinct(d$Plot) < 2) {
      summary_rows[[brite_category]] <- make_summary_row(
        level_name,
        brite_category,
        d,
        fit_status = "Skipped",
        error_msg = "Insufficient distinct timepoints or plots"
      )
      next
    }

    expression_variance <- var(d$Expression, na.rm = TRUE)
    variance_tolerance <- .Machine$double.eps * max(1, mean(abs(d$Expression), na.rm = TRUE))^2

    if (
      n_distinct(d$Expression) < 3 ||
      !is.finite(expression_variance) ||
      expression_variance <= variance_tolerance
    ) {
      summary_rows[[brite_category]] <- make_summary_row(
        level_name,
        brite_category,
        d,
        fit_status = "Skipped",
        error_msg = "Zero or near-zero expression variance"
      )
      next
    }

    warnings_seen <- character()
    error_msg <- NA_character_
    k_used <- min(time_basis_k, n_distinct(d$time_numeric) - 1L)

    model <- tryCatch(
      withCallingHandlers(
        gam(
          Expression ~ s(time_numeric, k = k_used) + s(Plot, bs = "re"),
          data = d,
          method = "REML"
        ),
        warning = function(w) {
          warnings_seen <<- c(warnings_seen, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        error_msg <<- conditionMessage(e)
        NULL
      }
    )

    warning_msg <- collapse_messages(warnings_seen)

    if (is.null(model)) {
      summary_rows[[brite_category]] <- make_summary_row(
        level_name,
        brite_category,
        d,
        fit_status = "Failed",
        warning_msg = warning_msg,
        error_msg = error_msg
      )
      next
    }

    fit_status <- case_when(
      !isTRUE(model$converged) ~ "Not converged",
      !is.na(warning_msg) ~ "Warning",
      TRUE ~ "Converged"
    )

    gam_summary <- summary(model)
    smooth_table <- as.data.frame(gam_summary$s.table)
    smooth_table$smooth_term <- rownames(smooth_table)

    time_row <- smooth_table %>%
      filter(smooth_term == "s(time_numeric)") %>%
      slice_head(n = 1)

    plot_row <- smooth_table %>%
      filter(str_detect(smooth_term, "^s\\(Plot\\)")) %>%
      slice_head(n = 1)

    time_predictions_observed <- predict(
      model,
      newdata = d,
      type = "response",
      exclude = "s(Plot)"
    )

    total_sum_squares <- sum(
      (d$Expression - mean(d$Expression, na.rm = TRUE))^2,
      na.rm = TRUE
    )

    time_residual_sum_squares <- sum(
      (d$Expression - time_predictions_observed)^2,
      na.rm = TRUE
    )

    time_variance_explained <- if (is.finite(total_sum_squares) && total_sum_squares > 0) {
      1 - time_residual_sum_squares / total_sum_squares
    } else {
      NA_real_
    }

    summary_rows[[brite_category]] <- make_summary_row(
      level_name,
      brite_category,
      d,
      fit_status = fit_status,
      model = model,
      gam_summary = gam_summary,
      time_row = time_row,
      plot_row = plot_row,
      time_variance_explained = time_variance_explained,
      warning_msg = warning_msg,
      error_msg = error_msg
    )

    reference_plot <- levels(d$Plot)[1]
    newdata <- time_grid %>%
      mutate(Plot = factor(reference_plot, levels = levels(d$Plot)))

    prediction <- predict(
      model,
      newdata = newdata,
      type = "response",
      se.fit = TRUE,
      exclude = "s(Plot)"
    )

    prediction_rows[[brite_category]] <- newdata %>%
      mutate(
        brite_level = level_name,
        brite_class = brite_category,
        fitted = as.numeric(prediction$fit),
        se = as.numeric(prediction$se.fit),
        lower = fitted - 1.96 * se,
        upper = fitted + 1.96 * se,
        prediction_type = "Population temporal trajectory; Plot RE excluded",
        fit_status = fit_status
      )
  }

  gam_summary_table <- bind_rows(summary_rows) %>%
    mutate(
      model_method = "mgcv::gam REML",
      model_formula = paste0(
        "Expression ~ s(time_numeric, k = ",
        time_basis_k,
        ") + s(Plot, bs = 're')"
      ),
      p_value_for_fdr = if_else(fit_status == accepted_fit_status, p_value, NA_real_),
      p_adj_fdr = p.adjust(p_value_for_fdr, method = "BH"),
      significant_time = if_else(!is.na(p_adj_fdr) & p_adj_fdr < 0.05, "Yes", "No")
    ) %>%
    select(-p_value_for_fdr) %>%
    arrange(p_adj_fdr, p_value, brite_class)

  gam_predictions <- bind_rows(prediction_rows) %>%
    left_join(
      gam_summary_table %>%
        select(
          brite_level,
          brite_class,
          p_value,
          p_adj_fdr,
          significant_time,
          deviance_explained_full,
          time_variance_explained
        ),
      by = c("brite_level", "brite_class")
    )

  model_data_output <- model_data %>%
    rename(brite_class = all_of(brite_col))

  write_csv(gam_summary_table, file.path(gam_results_dir, paste0(level_name, "_GAM_summary_table.csv")))
  write_csv(gam_predictions, file.path(gam_results_dir, paste0(level_name, "_GAM_fitted_trajectories.csv")))
  write_csv(model_data_output, file.path(gam_results_dir, paste0(level_name, "_GAM_long_model_data.csv")))
  write_csv(
    filter(gam_summary_table, significant_time == "Yes"),
    file.path(gam_results_dir, paste0(level_name, "_GAM_significant_time_classes.csv"))
  )

  significant_categories <- gam_summary_table %>%
    filter(significant_time == "Yes")

  if (nrow(significant_categories) > 0) {
    plot_points <- model_data_output %>%
      semi_join(significant_categories, by = "brite_class")

    plot_predictions <- gam_predictions %>%
      filter(significant_time == "Yes")

    temporal_plot <- ggplot(plot_points, aes(Collection_date, Expression)) +
      geom_point(alpha = 0.4, size = 1, color = "grey35") +
      geom_ribbon(
        data = plot_predictions,
        aes(Collection_date, ymin = lower, ymax = upper),
        inherit.aes = FALSE,
        fill = "#9ecae1",
        alpha = 0.35
      ) +
      geom_line(
        data = plot_predictions,
        aes(Collection_date, fitted),
        inherit.aes = FALSE,
        color = "#08519c",
        linewidth = 0.8
      ) +
      facet_wrap(~ brite_class, scales = "free_y") +
      labs(
        title = paste0("Factor 1 ", level_name, " significant temporal GAMs"),
        subtitle = "Population temporal trajectory with Plot random effect excluded from fitted line",
        x = "Collection date",
        y = "Expression"
      ) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 7)
      )

    ggsave(
      file.path(gam_figures_dir, paste0(level_name, "_significant_temporal_GAMs.png")),
      temporal_plot,
      width = 14,
      height = 10,
      dpi = 300,
      limitsize = FALSE
    )
  }

  list(
    model_data = model_data_output,
    summary_table = gam_summary_table,
    predictions = gam_predictions
  )
}

gam_results <- list(
  brite_C = fit_brite_level(factor1_brite_c_model_input, "brite_C", "brite_C"),
  brite_D = fit_brite_level(factor1_brite_d_model_input, "brite_D", "brite_D")
)

factor1_gam_run_summary <- bind_rows(
  gam_results$brite_C$summary_table,
  gam_results$brite_D$summary_table
) %>%
  group_by(brite_level) %>%
  summarise(
    n_features = n(),
    n_converged = sum(fit_status == "Converged"),
    n_warning = sum(fit_status == "Warning"),
    n_failed = sum(fit_status == "Failed"),
    n_skipped = sum(fit_status == "Skipped"),
    n_significant_time = sum(significant_time == "Yes"),
    .groups = "drop"
  )

write_csv(factor1_gam_run_summary, file.path(gam_results_dir, "Factor1_GAM_run_summary.csv"))

# ---- MicroResp overlap and clustering inputs --------------------------------

microresp_mapping <- read_csv(microresp_mapping_file, show_col_types = FALSE) %>%
  transmute(
    Substrate = as.character(Substrate),
    Substrate_category = as.character(Substrate_category),
    brite_C_mapping = as.character(brite_C),
    brite_D_mapping = as.character(brite_D),
    brite_C_code = str_to_lower(str_extract(brite_C_mapping, regex("ko\\d{5}", ignore_case = TRUE))),
    brite_D_code = str_to_upper(str_extract(brite_D_mapping, regex("K\\d{5}", ignore_case = TRUE)))
  )

factor1_brite_c_codes <- factor1_brite_c_summary %>%
  transmute(
    brite_class = brite_C,
    feature_code = str_to_lower(str_extract(brite_C, "ko\\d{5}")),
    n_factor1_genes
  )

factor1_brite_d_codes <- factor1_brite_d_summary %>%
  transmute(
    brite_class = brite_D,
    feature_code = str_to_upper(str_extract(brite_D, "^K\\d{5}")),
    n_factor1_genes
  )

factor1_microresp_overlap <- bind_rows(
  factor1_brite_c_codes %>%
    left_join(
      microresp_mapping %>%
        filter(!is.na(brite_C_code)) %>%
        group_by(brite_C_code) %>%
        summarise(
          n_microresp_substrates = n_distinct(Substrate),
          linked_substrates = paste(sort(unique(Substrate)), collapse = "; "),
          .groups = "drop"
        ),
      by = c("feature_code" = "brite_C_code")
    ) %>%
    mutate(
      brite_level = "brite_C",
      matched_microresp = !is.na(n_microresp_substrates)
    ),
  factor1_brite_d_codes %>%
    left_join(
      microresp_mapping %>%
        filter(!is.na(brite_D_code)) %>%
        group_by(brite_D_code) %>%
        summarise(
          n_microresp_substrates = n_distinct(Substrate),
          linked_substrates = paste(sort(unique(Substrate)), collapse = "; "),
          .groups = "drop"
        ),
      by = c("feature_code" = "brite_D_code")
    ) %>%
    mutate(
      brite_level = "brite_D",
      matched_microresp = !is.na(n_microresp_substrates)
    )
) %>%
  relocate(brite_level, brite_class, feature_code)

write_csv(factor1_microresp_overlap, file.path(results_dir, "factor1_microresp_overlap_table.csv"))

factor1_microresp_overlap_summary <- factor1_microresp_overlap %>%
  group_by(brite_level) %>%
  summarise(
    n_factor1_features = n(),
    n_microresp_matched_features = sum(matched_microresp),
    pct_microresp_matched_features = 100 * mean(matched_microresp),
    .groups = "drop"
  )

write_csv(factor1_microresp_overlap_summary, file.path(results_dir, "factor1_microresp_overlap_summary.csv"))

# ---- Cluster Factor 1 MicroResp-linked BRITE C and D ------------------------

make_paginated_gam_plots <- function(points, predictions, level, level_dir) {
  labels <- predictions %>% distinct(feature_code, plot_label)
  n_panels <- nrow(labels)
  n_pages <- ceiling(n_panels / gam_panels_per_page)
  if (n_panels == 0) {
    return(tibble(brite_level = level, n_panels = 0L, n_pages = 0L))
  }

  pdf_file <- file.path(level_dir, paste0(level, "_Factor1_targeted_GAMs_all_features.pdf"))
  pdf(pdf_file, width = 14, height = 8.5, onefile = TRUE)
  for (page in seq_len(n_pages)) {
    plot_page <- ggplot(points, aes(Collection_date, Expression)) +
      geom_point(alpha = 0.35, size = 0.9, color = "grey35") +
      geom_ribbon(
        data = predictions,
        aes(x = Collection_date, ymin = lower, ymax = upper),
        inherit.aes = FALSE,
        fill = "#9ecae1",
        alpha = 0.35
      ) +
      geom_line(
        data = predictions,
        aes(x = Collection_date, y = fitted),
        inherit.aes = FALSE,
        color = "#08519c",
        linewidth = 0.8
      ) +
      ggforce::facet_wrap_paginate(
        ~plot_label,
        scales = "free_y",
        ncol = 3,
        nrow = 2,
        page = page
      ) +
      labs(
        title = paste0("Factor 1 MicroResp-linked ", str_replace(level, "_", " "), " GAMs"),
        subtitle = paste0("All mapped, converged features, page ", page, " of ", n_pages),
        x = "Collection date",
        y = "Expression"
      ) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 7)
      )
    print(plot_page)
    ggsave(
      file.path(level_dir, paste0(level, "_Factor1_targeted_GAMs_page_", sprintf("%02d", page), ".png")),
      plot_page,
      width = 14,
      height = 8.5,
      dpi = 300
    )
  }
  dev.off()

  tibble(
    brite_level = level,
    n_panels = n_panels,
    panels_per_page = gam_panels_per_page,
    n_pages = n_pages
  )
}

analyze_level <- function(level) {
  level_dir <- file.path(clustering_figures_dir, level)
  dir.create(level_dir, recursive = TRUE, showWarnings = FALSE)

  code_col <- if (level == "brite_C") "brite_C_code" else "brite_D_code"
  mapping_text_col <- if (level == "brite_C") "brite_C_mapping" else "brite_D_mapping"

  summary <- read_csv(
    file.path(gam_results_dir, paste0(level, "_GAM_summary_table.csv")),
    show_col_types = FALSE
  ) %>%
    mutate(
      feature_code = if (level == "brite_C") {
        str_to_lower(str_extract(brite_class, "ko\\d{5}"))
      } else {
        str_to_upper(str_extract(brite_class, "^K\\d{5}"))
      }
    )

  trajectories <- read_csv(
    file.path(gam_results_dir, paste0(level, "_GAM_fitted_trajectories.csv")),
    show_col_types = FALSE
  ) %>%
    mutate(
      Collection_date = as.Date(Collection_date),
      feature_code = if (level == "brite_C") {
        str_to_lower(str_extract(brite_class, "ko\\d{5}"))
      } else {
        str_to_upper(str_extract(brite_class, "^K\\d{5}"))
      }
    )

  model_data <- read_csv(
    file.path(gam_results_dir, paste0(level, "_GAM_long_model_data.csv")),
    show_col_types = FALSE
  ) %>%
    mutate(
      Collection_date = as.Date(Collection_date),
      feature_code = if (level == "brite_C") {
        str_to_lower(str_extract(brite_class, "ko\\d{5}"))
      } else {
        str_to_upper(str_extract(brite_class, "^K\\d{5}"))
      }
    )

  map_level <- microresp_mapping %>%
    transmute(
      Substrate,
      Substrate_category,
      mapped_brite_class = .data[[mapping_text_col]],
      feature_code = .data[[code_col]]
    ) %>%
    filter(!is.na(feature_code)) %>%
    distinct() %>%
    inner_join(
      factor1_microresp_overlap %>%
        filter(brite_level == level, matched_microresp) %>%
        select(feature_code),
      by = "feature_code"
    )

  target_features <- map_level %>%
    group_by(feature_code) %>%
    summarise(
      mapped_brite_class = first(mapped_brite_class),
      linked_substrates = paste(sort(unique(Substrate)), collapse = "; "),
      substrate_categories = paste(sort(unique(Substrate_category)), collapse = "; "),
      n_linked_substrates = n_distinct(Substrate),
      .groups = "drop"
    ) %>%
    left_join(summary, by = "feature_code") %>%
    left_join(
      if (level == "brite_C") factor1_brite_c_codes else factor1_brite_d_codes,
      by = c("feature_code", "brite_class")
    ) %>%
    mutate(
      p_value_target_set = if_else(fit_status == "Converged" & converged == "Yes", p_value, NA_real_),
      p_adj_fdr_microresp_targets = p.adjust(p_value_target_set, method = "BH"),
      significant_time_microresp_targets = if_else(
        !is.na(p_adj_fdr_microresp_targets) & p_adj_fdr_microresp_targets < 0.05,
        "Yes",
        "No"
      )
    ) %>%
    select(-p_value_target_set) %>%
    arrange(p_adj_fdr_microresp_targets, p_value, feature_code)

  write_csv(target_features, file.path(clustering_results_dir, paste0(level, "_Factor1_MicroResp_GAM_statistics.csv")))

  eligible <- target_features %>%
    filter(fit_status == "Converged", converged == "Yes", !is.na(brite_class))

  if (nrow(eligible) < 3) {
    write_csv(
      tibble(
        brite_level = level,
        note = "Fewer than three mapped, converged Factor 1 MicroResp-linked features. Clustering skipped."
      ),
      file.path(clustering_results_dir, paste0(level, "_Factor1_clustering_skipped.csv"))
    )
    return(tibble(
      brite_level = level,
      n_mapping_rows = nrow(map_level),
      n_unique_requested_features = n_distinct(map_level$feature_code),
      n_matched_features = sum(!is.na(target_features$brite_class)),
      n_converged_features = nrow(eligible),
      n_clustered_features = 0,
      n_target_FDR_significant = sum(target_features$significant_time_microresp_targets == "Yes", na.rm = TRUE),
      automatic_elbow_k = NA_integer_,
      selected_k = NA_integer_,
      heatmap_display_rule = "Clustering skipped"
    ))
  }

  eligible_codes <- eligible$feature_code
  target_trajectories <- trajectories %>% filter(feature_code %in% eligible_codes)
  target_points <- model_data %>% filter(feature_code %in% eligible_codes)

  dates <- target_trajectories %>%
    distinct(Collection_date, time_numeric) %>%
    arrange(time_numeric)

  fitted_matrix <- target_trajectories %>%
    mutate(date_id = paste0("day_", time_numeric)) %>%
    select(feature_code, date_id, fitted) %>%
    pivot_wider(names_from = date_id, values_from = fitted) %>%
    column_to_rownames("feature_code") %>%
    as.matrix()

  fitted_matrix <- fitted_matrix[, paste0("day_", dates$time_numeric), drop = FALSE]

  zero_variance <- rownames(fitted_matrix)[
    !is.finite(apply(fitted_matrix, 1, sd)) | apply(fitted_matrix, 1, sd) == 0
  ]

  if (length(zero_variance) > 0) {
    fitted_matrix <- fitted_matrix[!rownames(fitted_matrix) %in% zero_variance, , drop = FALSE]
    eligible <- eligible %>% filter(!feature_code %in% zero_variance)
    target_trajectories <- target_trajectories %>% filter(!feature_code %in% zero_variance)
    target_points <- target_points %>% filter(!feature_code %in% zero_variance)
  }

  zmat <- row_zscore(fitted_matrix)
  zmat[zmat < -3] <- -3
  zmat[zmat > 3] <- 3
  colnames(zmat) <- format(dates$Collection_date, "%b'%y")

  write_csv(
    as.data.frame(zmat) %>% rownames_to_column("feature_code"),
    file.path(clustering_results_dir, paste0(level, "_Factor1_MicroResp_zscored_fitted_trajectories.csv"))
  )

  correlation <- cor(t(zmat), use = "pairwise.complete.obs")
  distance_matrix <- 1 - correlation
  distance_matrix[distance_matrix < 0] <- 0
  diag(distance_matrix) <- 0
  distance <- as.dist(distance_matrix)

  hc <- hclust(distance, method = "ward.D2")
  kmax <- min(max_k_to_test, nrow(zmat) - 1L)
  elbow <- map_dfr(1:kmax, function(k) {
    ids <- cutree(hc, k)
    tibble(k = k, total_within_dissimilarity = total_within_dissimilarity(ids, distance))
  })

  automatic_k <- detect_elbow_k(elbow)
  selected_k <- if (is.na(manual_k[[level]])) automatic_k else manual_k[[level]]

  elbow <- elbow %>%
    mutate(
      automatic_elbow_k = automatic_k,
      selected_k = selected_k,
      selected = k == selected_k
    )

  write_csv(elbow, file.path(clustering_results_dir, paste0(level, "_Factor1_MicroResp_elbow_diagnostics.csv")))

  elbow_plot <- ggplot(elbow, aes(k, total_within_dissimilarity)) +
    geom_line() +
    geom_point() +
    geom_vline(xintercept = selected_k, linetype = "dashed", color = "#0072B2") +
    scale_x_continuous(breaks = 1:kmax) +
    labs(
      title = paste0("Factor 1 MicroResp-linked ", str_replace(level, "_", " "), " elbow diagnostic"),
      x = "Number of clusters",
      y = "Within-cluster dissimilarity"
    ) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(level_dir, paste0(level, "_Factor1_elbow_diagnostic.png")), elbow_plot, width = 8, height = 5.5, dpi = 300)

  raw_ids <- cutree(hc, selected_k)
  raw_order <- tibble(
    feature_code = names(raw_ids),
    raw_cluster = as.integer(raw_ids)
  ) %>%
    mutate(peak = apply(zmat[feature_code, , drop = FALSE], 1, which.max)) %>%
    group_by(raw_cluster) %>%
    summarise(mean_peak = mean(peak), .groups = "drop") %>%
    arrange(mean_peak, raw_cluster) %>%
    mutate(Cluster_ID = paste0("Cluster_", row_number()))

  assignments <- tibble(
    feature_code = names(raw_ids),
    raw_cluster = as.integer(raw_ids)
  ) %>%
    left_join(raw_order %>% select(raw_cluster, Cluster_ID), by = "raw_cluster") %>%
    left_join(eligible, by = "feature_code") %>%
    mutate(
      Temporal_significance = if_else(significant_time_microresp_targets == "Yes", "Significant", "Not significant"),
      display_label = paste(feature_code, short_label(brite_class, level)),
      dend_position = match(feature_code, hc$labels[hc$order])
    ) %>%
    arrange(Cluster_ID, dend_position)

  write_csv(assignments %>% select(-dend_position), file.path(clustering_results_dir, paste0(level, "_Factor1_MicroResp_cluster_assignments.csv")))

  cluster_colors <- setNames(hcl.colors(selected_k, "Dark 3"), paste0("Cluster_", seq_len(selected_k)))
  dend <- as.dendrogram(hc) %>% set("branches_lwd", 0.9)

  if (nrow(assignments) <= 80) {
    labels(dend) <- setNames(assignments$display_label, assignments$feature_code)[labels(dend)]
    dend <- dend %>% set("labels_cex", 0.5)
  } else {
    dend <- dend %>% set("labels", rep("", length(labels(dend))))
  }

  draw_dend <- function() {
    plot(
      dend,
      main = paste0("Factor 1 MicroResp-linked ", str_replace(level, "_", " "), " clustering"),
      ylab = "Ward.D2 height",
      leaflab = "perpendicular"
    )
    rect.hclust(hc, k = selected_k, border = cluster_colors)
  }

  png(file.path(level_dir, paste0(level, "_Factor1_dendrogram.png")), width = 12, height = 7, units = "in", res = 300)
  draw_dend()
  dev.off()

  pdf(file.path(level_dir, paste0(level, "_Factor1_dendrogram.pdf")), width = 12, height = 7)
  draw_dend()
  dev.off()

  trend_data <- as.data.frame(zmat) %>%
    rownames_to_column("feature_code") %>%
    pivot_longer(-feature_code, names_to = "date_label", values_to = "z_fitted") %>%
    left_join(assignments %>% select(feature_code, Cluster_ID), by = "feature_code") %>%
    mutate(
      date_label = factor(date_label, levels = colnames(zmat)),
      Cluster_ID = factor(Cluster_ID, levels = names(cluster_colors))
    )

  trend_summary <- trend_data %>%
    group_by(Cluster_ID, date_label) %>%
    summarise(
      mean_z = mean(z_fitted),
      se_z = if_else(n() > 1, sd(z_fitted) / sqrt(n()), 0),
      n_features = n(),
      .groups = "drop"
    )

  write_csv(trend_summary, file.path(clustering_results_dir, paste0(level, "_Factor1_MicroResp_cluster_average_trajectories.csv")))

  trend_plot <- ggplot(
    trend_summary,
    aes(date_label, mean_z, group = Cluster_ID, color = Cluster_ID, fill = Cluster_ID)
  ) +
    geom_ribbon(aes(ymin = mean_z - se_z, ymax = mean_z + se_z), alpha = 0.18, color = NA) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    scale_color_manual(values = cluster_colors) +
    scale_fill_manual(values = cluster_colors) +
    labs(
      title = paste0("Factor 1 MicroResp-linked ", str_replace(level, "_", " "), " cluster-average trajectories"),
      subtitle = paste0("All mapped, converged features, k = ", selected_k),
      x = NULL,
      y = "Mean z-scored fitted expression"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.minor = element_blank())

  ggsave(file.path(level_dir, paste0(level, "_Factor1_cluster_average_trajectories.png")), trend_plot, width = 9, height = 5.5, dpi = 300)
  ggsave(file.path(level_dir, paste0(level, "_Factor1_cluster_average_trajectories.pdf")), trend_plot, width = 9, height = 5.5)

  if (level == "brite_C") {
    heatmap_assignments <- assignments
    heatmap_description <- "All mapped, converged Factor 1 BRITE C pathways"
  } else {
    heatmap_assignments <- assignments %>%
      group_by(Cluster_ID) %>%
      arrange(desc(time_variance_explained), p_adj_fdr_microresp_targets, .by_group = TRUE) %>%
      slice_head(n = top_brite_D_per_cluster) %>%
      ungroup()
    heatmap_description <- paste0(
      "Top ",
      top_brite_D_per_cluster,
      " Factor 1 BRITE D KOs per cluster ranked by time-only variance explained"
    )
    write_csv(
      heatmap_assignments,
      file.path(clustering_results_dir, "brite_D_Factor1_MicroResp_heatmap_feature_selection.csv")
    )
  }

  save_heatmap(zmat, heatmap_assignments, level, "Factor1_MicroResp_heatmap", heatmap_description, level_dir)

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

  write_csv(gam_page_summary, file.path(clustering_results_dir, paste0(level, "_Factor1_MicroResp_GAM_page_summary.csv")))

  tibble(
    brite_level = level,
    n_mapping_rows = nrow(map_level),
    n_unique_requested_features = n_distinct(map_level$feature_code),
    n_matched_features = sum(!is.na(target_features$brite_class)),
    n_converged_features = nrow(eligible),
    n_clustered_features = nrow(zmat),
    n_target_FDR_significant = sum(assignments$Temporal_significance == "Significant"),
    automatic_elbow_k = automatic_k,
    selected_k = selected_k,
    heatmap_display_rule = heatmap_description
  )
}

clustering_summary <- bind_rows(analyze_level("brite_C"), analyze_level("brite_D"))
write_csv(clustering_summary, file.path(clustering_results_dir, "Factor1_MicroResp_clustering_run_summary.csv"))

write_csv(
  tibble(
    factor1_file = factor1_file,
    expression_file = expr_file,
    annotation_file = annotation_file,
    microresp_mapping_file = microresp_mapping_file,
    model = paste0("Expression ~ s(time_numeric, k = ", time_basis_k, ") + s(Plot, bs = 're')"),
    clustering = "Population GAM trajectories, row z-score, 1 - Pearson correlation, Ward.D2",
    run_date = as.character(Sys.Date())
  ),
  file.path(results_dir, "run_manifest.csv")
)

cat("\nFactor 1 expression, BRITE, GAM, and MicroResp clustering analysis complete.\n")
print(factor1_gene_summary)
print(factor1_gam_run_summary)
print(clustering_summary)
