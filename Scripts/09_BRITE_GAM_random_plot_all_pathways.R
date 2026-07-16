############################################################
# BRITE C and BRITE D temporal GAMs with Plot as a random effect
#
# This script refits all filtered BRITE C and BRITE D features using
# mgcv::gam() rather than the mgcv::gamm()/nlme wrapper. Plot remains a
# random effect through s(Plot, bs = "re"). The complete gam object provides
# valid full-model deviance explained and clearer convergence diagnostics.
#
# Outputs:
#   Results/BRITE_temporal_gam_random_plot_all_pathways/
#   Figures/BRITE_temporal_gam_random_plot_all_pathways/
############################################################

rm(list = ls())

required_packages <- c("mgcv", "tidyverse", "readr")
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
})

theme_set(theme_bw(base_size = 11))

# ---- User settings ----------------------------------------------------------

run_label <- "gam_random_plot_all_pathways"

# There are seven observed dates. k = 6 allows a nonlinear trajectory without
# using the maximum flexibility available from those seven dates.
time_basis_k <- 6L

# Only models with this status contribute p-values to the FDR correction.
accepted_fit_status <- "Converged"

results_dir <- file.path("Results", paste0("BRITE_temporal_", run_label))
figures_dir <- file.path("Figures", paste0("BRITE_temporal_", run_label))
diagnostics_dir <- file.path(results_dir, "diagnostics")

for (directory in c(results_dir, figures_dir, diagnostics_dir)) {
  dir.create(directory, showWarnings = FALSE, recursive = TRUE)
}

# ---- Input data -------------------------------------------------------------

meta <- read.csv("Data/Data_merged_cleaned.csv", check.names = FALSE) %>%
  mutate(
    Collection_date = as.Date(Collection_date),
    Plot = factor(Plot),
    time_numeric = as.numeric(
      Collection_date - min(Collection_date, na.rm = TRUE)
    )
  ) %>%
  select(
    Sample_cat,
    Collection_date,
    Plot,
    Irrigation_trmt_pct,
    Irrigation_trmt_cat,
    time_numeric
  )

remove_text_file <- file.path(
  "Data",
  "pathwaystoremove",
  "kegg_human_exclude_2026-07-13.txt"
)

human_remove_table <- read.delim(
  remove_text_file,
  check.names = FALSE
) %>%
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

brite_D <- read.csv("Data/brite_D_expression.csv", check.names = FALSE) %>%
  mutate(
    brite_D = as.character(brite_D),
    ko_code = str_extract(brite_D, "^K\\d{5}"),
    remove_from_text = ko_code %in% brite_d_remove_codes$KO_entry
  )

brite_C <- read.csv("Data/brite_C_expression.csv", check.names = FALSE) %>%
  mutate(
    brite_C = as.character(brite_C),
    pathway_code = str_extract(brite_C, "ko\\d{5}"),
    remove_from_text = pathway_code %in% brite_c_remove_codes$path_ko
  )

exclusion_summary <- bind_rows(
  brite_D %>%
    summarise(
      brite_level = "brite_D",
      starting_rows = n(),
      removed_from_text = sum(remove_from_text, na.rm = TRUE),
      remaining_rows = sum(!remove_from_text, na.rm = TRUE)
    ),
  brite_C %>%
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

write_csv(
  exclusion_summary,
  file.path(results_dir, "BRITE_exclusion_summary.csv")
)

write_csv(
  filter(brite_D, remove_from_text),
  file.path(results_dir, "brite_D_removed_classes.csv")
)

write_csv(
  filter(brite_C, remove_from_text),
  file.path(results_dir, "brite_C_removed_classes.csv")
)

brite_D <- brite_D %>%
  filter(!remove_from_text) %>%
  select(-ko_code, -remove_from_text)

brite_C <- brite_C %>%
  filter(!remove_from_text) %>%
  select(-pathway_code, -remove_from_text)

# ---- Helper functions -------------------------------------------------------

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
    expression_variance = if (nrow(data) > 1) {
      var(data$Expression, na.rm = TRUE)
    } else {
      NA_real_
    },
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
    model_deviance = if (!is.null(model)) {
      safe_number(deviance(model))
    } else {
      NA_real_
    },
    null_deviance = if (!is.null(model)) {
      safe_number(model$null.deviance)
    } else {
      NA_real_
    },
    scale = safe_number(gam_summary$scale),
    AIC = if (!is.null(model)) safe_number(AIC(model)) else NA_real_,
    fit_status = fit_status,
    converged = if_else(fit_status == accepted_fit_status, "Yes", "No"),
    model_warning = warning_msg,
    model_error = error_msg
  )
}

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
      irrigation_values = paste(
        sort(unique(Irrigation_trmt_pct)),
        collapse = ";"
      ),
      min_date = min(Collection_date),
      max_date = max(Collection_date)
    ) %>%
    write_csv(
      file.path(diagnostics_dir, paste0(level_name, "_input_summary.csv"))
    )

  time_grid <- model_data %>%
    distinct(time_numeric, Collection_date) %>%
    arrange(time_numeric, Collection_date)

  brite_categories <- unique(model_data[[brite_col]])
  summary_rows <- vector("list", length(brite_categories))
  prediction_rows <- vector("list", length(brite_categories))

  names(summary_rows) <- brite_categories
  names(prediction_rows) <- brite_categories

  cat(
    "\n=== FITTING ",
    level_name,
    " RANDOM-PLOT GAMS (",
    length(brite_categories),
    " categories) ===\n",
    sep = ""
  )

  for (i in seq_along(brite_categories)) {
    brite_category <- brite_categories[i]

    if (i %% 100 == 0) {
      cat("  ", level_name, ": ", i, " of ", length(brite_categories), "\n")
    }

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
    variance_tolerance <- .Machine$double.eps *
      max(1, mean(abs(d$Expression), na.rm = TRUE))^2

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
          Expression ~
            s(time_numeric, k = k_used) +
            s(Plot, bs = "re"),
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

    time_variance_explained <- if (
      is.finite(total_sum_squares) && total_sum_squares > 0
    ) {
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
      exclusion_rule = "txt file only",
      model_method = "mgcv::gam REML",
      model_formula = paste0(
        "Expression ~ s(time_numeric, k = ",
        time_basis_k,
        ") + s(Plot, bs = 're')"
      ),
      p_value_for_fdr = if_else(
        fit_status == accepted_fit_status,
        p_value,
        NA_real_
      ),
      p_adj_fdr = p.adjust(p_value_for_fdr, method = "BH"),
      significant_time = if_else(
        !is.na(p_adj_fdr) & p_adj_fdr < 0.05,
        "Yes",
        "No"
      )
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

  write_csv(
    gam_summary_table,
    file.path(results_dir, paste0(level_name, "_GAM_summary_table.csv"))
  )

  write_csv(
    gam_predictions,
    file.path(results_dir, paste0(level_name, "_GAM_fitted_trajectories.csv"))
  )

  write_csv(
    model_data_output,
    file.path(results_dir, paste0(level_name, "_GAM_long_model_data.csv"))
  )

  write_csv(
    filter(gam_summary_table, significant_time == "Yes"),
    file.path(
      results_dir,
      paste0(level_name, "_GAM_significant_time_classes.csv")
    )
  )

  fit_status_summary <- gam_summary_table %>%
    count(fit_status, name = "n_features") %>%
    mutate(brite_level = level_name) %>%
    relocate(brite_level)

  write_csv(
    fit_status_summary,
    file.path(diagnostics_dir, paste0(level_name, "_fit_status_summary.csv"))
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
        title = paste0(
          level_name,
          ": significant temporal GAMs with random Plot effect"
        ),
        subtitle = "Lines show the population temporal trajectory; Plot random effects are excluded",
        x = "Collection date",
        y = "Expression"
      ) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 7)
      )

    ggsave(
      file.path(
        figures_dir,
        paste0(level_name, "_significant_temporal_GAMs.png")
      ),
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
    predictions = gam_predictions,
    fit_status_summary = fit_status_summary
  )
}

# ---- Run BRITE C and BRITE D separately ------------------------------------

results <- list(
  brite_C = fit_brite_level(brite_C, "brite_C", "brite_C"),
  brite_D = fit_brite_level(brite_D, "brite_D", "brite_D")
)

combined_summary <- bind_rows(
  results$brite_C$summary_table,
  results$brite_D$summary_table
)

write_csv(
  combined_summary,
  file.path(results_dir, "BRITE_GAM_summary_table_all_levels.csv")
)

run_summary <- combined_summary %>%
  group_by(brite_level) %>%
  summarise(
    n_features = n(),
    n_converged = sum(fit_status == "Converged"),
    n_warning = sum(fit_status == "Warning"),
    n_failed = sum(fit_status == "Failed"),
    n_skipped = sum(fit_status == "Skipped"),
    n_significant_time = sum(significant_time == "Yes"),
    n_deviance_explained_available = sum(
      !is.na(deviance_explained_full)
    ),
    .groups = "drop"
  )

write_csv(
  run_summary,
  file.path(results_dir, "BRITE_GAM_run_summary.csv")
)

run_manifest <- tibble(
  run_label = run_label,
  exclusion_source = basename(remove_text_file),
  model_method = "mgcv::gam REML",
  model_formula = paste0(
    "Expression ~ s(time_numeric, k = ",
    time_basis_k,
    ") + s(Plot, bs = 're')"
  ),
  full_deviance_definition = "Time smooth plus Plot random effect",
  time_variance_definition = "Population time smooth only; Plot random effect excluded",
  fdr_scope = "Separately within all converged BRITE C and BRITE D features",
  results_dir = results_dir,
  figures_dir = figures_dir,
  run_date = as.character(Sys.Date())
)

write_csv(run_manifest, file.path(results_dir, "run_manifest.csv"))

cat("\nRandom-Plot GAM analysis complete.\n")
print(run_summary)
cat("\nResults: ", results_dir, "\n", sep = "")
cat("Figures: ", figures_dir, "\n", sep = "")

