rm(list = ls())

required_packages <- c("mgcv", "tidyverse", "broom", "readr")
missing <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(mgcv)
  library(tidyverse)
  library(broom)
  library(readr)
})

theme_set(theme_bw(base_size = 11))

run_label <- "txt_only_significant_only"
results_dir <- file.path("Results", paste0("BRITE_temporal_", run_label))
figures_dir <- file.path("Figures", paste0("BRITE_temporal_", run_label))
diagnostics_dir <- file.path(results_dir, "diagnostics")

for (d in c(results_dir, figures_dir, diagnostics_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

meta <- read.csv("Data/Data_merged_cleaned.csv", check.names = FALSE) %>%
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

remove_text_file <- file.path("Data", "pathwaystoremove", "kegg_human_exclude_2026-07-13.txt")

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
      exclusion_rule = run_label,
      source_file = basename(remove_text_file),
      starting_rows = n(),
      removed_from_text = sum(remove_from_text, na.rm = TRUE),
      remaining_rows = sum(!remove_from_text, na.rm = TRUE)
    ),
  brite_C %>%
    summarise(
      brite_level = "brite_C",
      exclusion_rule = run_label,
      source_file = basename(remove_text_file),
      starting_rows = n(),
      removed_from_text = sum(remove_from_text, na.rm = TRUE),
      remaining_rows = sum(!remove_from_text, na.rm = TRUE)
    )
)

write_csv(exclusion_summary, file.path(results_dir, "BRITE_exclusion_summary.csv"))
write_csv(filter(brite_D, remove_from_text), file.path(results_dir, "brite_D_removed_classes.csv"))
write_csv(filter(brite_C, remove_from_text), file.path(results_dir, "brite_C_removed_classes.csv"))

brite_D <- brite_D %>%
  filter(!remove_from_text) %>%
  select(-ko_code, -remove_from_text)

brite_C <- brite_C %>%
  filter(!remove_from_text) %>%
  select(-pathway_code, -remove_from_text)

write_csv(brite_D, file.path(results_dir, "brite_D_expression_filtered_for_GAMs.csv"))
write_csv(brite_C, file.path(results_dir, "brite_C_expression_filtered_for_GAMs.csv"))

safe_val <- function(x) if (length(x) == 1 && !is.null(x)) x else NA_real_

make_summary_row <- function(level_name, brite_category, data,
                             gam_summary = NULL, smooth_row = NULL,
                             model = NULL, converged = "No",
                             warning_msg = NA_character_, error_msg = NA_character_) {
  tibble(
    brite_level = level_name,
    brite_class = brite_category,
    n_obs = nrow(data),
    n_plots = n_distinct(data$Plot),
    n_timepoints = n_distinct(data$time_numeric),
    edf = safe_val(smooth_row$edf),
    ref_df = safe_val(smooth_row$Ref.df),
    F_statistic = safe_val(smooth_row$F),
    p_value = safe_val(smooth_row[["p-value"]]),
    adj_r_squared = safe_val(gam_summary$r.sq),
    deviance_explained = safe_val(gam_summary$dev.expl),
    scale = safe_val(gam_summary$scale),
    AIC = if (!is.null(model)) safe_val(AIC(model$lme)) else NA_real_,
    converged = converged,
    model_warning = warning_msg,
    model_error = error_msg
  )
}

analyze_brite_level <- function(brite_table, brite_col, level_name) {
  sample_cols <- setdiff(names(brite_table), brite_col)

  model_data <- brite_table %>%
    filter(!is.na(.data[[brite_col]]), .data[[brite_col]] != "") %>%
    pivot_longer(all_of(sample_cols), names_to = "Sample_mt", values_to = "Expression") %>%
    mutate(Sample_cat = str_remove(Sample_mt, "_MT$")) %>%
    left_join(meta, by = "Sample_cat") %>%
    mutate(
      !!brite_col := as.character(.data[[brite_col]]),
      Expression = as.numeric(Expression)
    ) %>%
    drop_na(Collection_date, Plot, time_numeric, Expression)

  model_data %>%
    summarise(
      exclusion_rule = run_label,
      n_rows = n(),
      n_brite_categories = n_distinct(.data[[brite_col]]),
      n_samples = n_distinct(Sample_cat),
      min_date = min(Collection_date),
      max_date = max(Collection_date)
    ) %>%
    write_csv(file.path(diagnostics_dir, paste0(level_name, "_input_summary.csv")))

  time_grid <- model_data %>%
    distinct(time_numeric, Collection_date) %>%
    arrange(time_numeric, Collection_date)

  brite_categories <- unique(model_data[[brite_col]])
  reference_plot <- levels(model_data$Plot)[1]

  summary_rows <- list()
  prediction_list <- list()

  cat("\n=== FITTING", level_name, "GAMMS (", length(brite_categories), "categories ) ===\n")

  for (i in seq_along(brite_categories)) {
    brite_category <- brite_categories[i]

    if (i %% 100 == 0) {
      cat("  ", level_name, ":", i, "of", length(brite_categories), "\n")
    }

    d <- model_data %>%
      filter(.data[[brite_col]] == brite_category) %>%
      arrange(time_numeric, Sample_cat)

    if (n_distinct(d$time_numeric) < 4 || n_distinct(d$Plot) < 2) {
      summary_rows[[brite_category]] <- make_summary_row(
        level_name, brite_category, d,
        error_msg = "Insufficient distinct timepoints or plots for GAMM"
      )
      next
    }

    warning_msg <- NA_character_
    error_msg <- NA_character_

    model <- tryCatch(
      withCallingHandlers(
        gamm(
          Expression ~ s(time_numeric, k = min(10, n_distinct(d$time_numeric))),
          random = list(Plot = ~1),
          data = d,
          method = "REML"
        ),
        warning = function(w) {
          warning_msg <<- conditionMessage(w)
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        error_msg <<- conditionMessage(e)
        NULL
      }
    )

    if (is.null(model)) {
      summary_rows[[brite_category]] <- make_summary_row(
        level_name, brite_category, d,
        warning_msg = warning_msg,
        error_msg = error_msg
      )
      next
    }

    gam_summary <- summary(model$gam)
    smooth_tbl <- as.data.frame(gam_summary$s.table)
    smooth_row <- smooth_tbl[1, , drop = FALSE]

    summary_rows[[brite_category]] <- make_summary_row(
      level_name, brite_category, d,
      gam_summary = gam_summary,
      smooth_row = smooth_row,
      model = model,
      converged = "Yes",
      warning_msg = warning_msg,
      error_msg = error_msg
    )

    newdata <- time_grid %>%
      mutate(Plot = factor(reference_plot, levels = levels(model_data$Plot)))

    pred <- predict(model$gam, newdata = newdata, type = "response", se.fit = TRUE)

    prediction_list[[brite_category]] <- newdata %>%
      mutate(
        brite_level = level_name,
        brite_class = brite_category,
        fitted = as.numeric(pred$fit),
        se = as.numeric(pred$se.fit),
        lower = fitted - 1.96 * se,
        upper = fitted + 1.96 * se
      )
  }

  gam_summary_table <- bind_rows(summary_rows) %>%
    mutate(
      exclusion_rule = run_label,
      p_adj_fdr = p.adjust(p_value, method = "BH"),
      significant_time = if_else(!is.na(p_adj_fdr) & p_adj_fdr < 0.05, "Yes", "No")
    ) %>%
    arrange(p_adj_fdr, p_value, brite_class)

  gam_predictions <- bind_rows(prediction_list) %>%
    left_join(
      gam_summary_table %>%
        select(brite_level, brite_class, p_value, p_adj_fdr, significant_time),
      by = c("brite_level", "brite_class")
    ) %>%
    mutate(exclusion_rule = run_label)

  significant_categories <- gam_summary_table %>%
    filter(significant_time == "Yes")

  write_csv(gam_summary_table, file.path(results_dir, paste0(level_name, "_GAM_summary_table.csv")))
  write_csv(gam_predictions, file.path(results_dir, paste0(level_name, "_GAM_fitted_trajectories.csv")))
  write_csv(model_data, file.path(results_dir, paste0(level_name, "_GAM_long_model_data.csv")))
  write_csv(significant_categories, file.path(results_dir, paste0(level_name, "_GAM_significant_time_classes.csv")))

  if (nrow(significant_categories) > 0) {
    brite_order_sig <- significant_categories$brite_class

    plot_points_sig <- model_data %>%
      rename(brite_class = all_of(brite_col)) %>%
      semi_join(significant_categories, by = "brite_class") %>%
      mutate(brite_class = factor(brite_class, levels = brite_order_sig))

    plot_pred_sig <- gam_predictions %>%
      filter(significant_time == "Yes") %>%
      mutate(brite_class = factor(brite_class, levels = brite_order_sig))

    temporal_plot_sig <- ggplot(plot_points_sig, aes(Collection_date, Expression)) +
      geom_point(alpha = 0.45, size = 1.1, color = "grey35") +
      geom_ribbon(
        data = plot_pred_sig,
        aes(Collection_date, ymin = lower, ymax = upper),
        inherit.aes = FALSE,
        fill = "#9ecae1",
        alpha = 0.35
      ) +
      geom_line(
        data = plot_pred_sig,
        aes(Collection_date, fitted),
        inherit.aes = FALSE,
        color = "#08519c",
        linewidth = 0.8
      ) +
      facet_wrap(~ brite_class, scales = "free_y") +
      labs(
        title = paste0(level_name, ": significant temporal GAMs after txt-only exclusion"),
        subtitle = "Exclusion rule: txt_only_significant_only",
        x = "Collection date",
        y = "Expression"
      ) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 8)
      )

    ggsave(
      file.path(figures_dir, paste0(level_name, "_significant_temporal_GAMs.png")),
      temporal_plot_sig,
      width = 14,
      height = 10,
      dpi = 300
    )
    ggsave(
      file.path(figures_dir, paste0(level_name, "_significant_temporal_GAMs.pdf")),
      temporal_plot_sig,
      width = 14,
      height = 10
    )
  }

  cat(level_name, "classes modeled:", nrow(gam_summary_table), "\n")
  cat(level_name, "significant temporal classes (FDR < 0.05):", nrow(significant_categories), "\n")

  list(
    model_data = model_data,
    summary_table = gam_summary_table,
    predictions = gam_predictions
  )
}

results <- list(
  brite_D = analyze_brite_level(brite_D, "brite_D", "brite_D"),
  brite_C = analyze_brite_level(brite_C, "brite_C", "brite_C")
)

combined_summary <- bind_rows(
  results$brite_D$summary_table,
  results$brite_C$summary_table
)

write_csv(combined_summary, file.path(results_dir, "BRITE_GAM_summary_table_all_levels.csv"))

run_manifest <- tibble(
  run_label = run_label,
  exclusion_source = basename(remove_text_file),
  exclusion_logic = "txt file only",
  plot_logic = "significant classes only, FDR < 0.05",
  results_dir = results_dir,
  figures_dir = figures_dir,
  run_date = as.character(Sys.Date())
)

write_csv(run_manifest, file.path(results_dir, "run_manifest.csv"))

cat("\nFiles saved in:\n")
cat(results_dir, "\n")
cat(figures_dir, "\n")
