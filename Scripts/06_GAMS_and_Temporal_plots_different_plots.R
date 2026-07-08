rm(list = ls())
## ---- Setup -----------------------------------------------------------------
required_packages <- c("mgcv", "tidyverse", "broom", "readr")
missing <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(mgcv); library(tidyverse); library(broom); library(readr)
})
theme_set(theme_bw(base_size = 11))
for (d in c("Figures", "Results", "Results/Step2_diagnostics")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

## ---- Metadata --------------------------------------------------------------
meta <- read.csv("Data/Data_merged_cleaned.csv", check.names = FALSE) %>%
  mutate(
    Collection_date = as.Date(Collection_date),
    Plot            = factor(Plot),
    time_numeric    = as.numeric(Collection_date - min(Collection_date, na.rm = TRUE))
  ) %>%
  select(Sample_cat, Collection_date, Plot,
         Irrigation_trmt_pct, Irrigation_trmt_cat, time_numeric)

## ---- Helpers ---------------------------------------------------------------
## Return x only if it is a single, non-NULL value; otherwise NA.
## This guards against BOTH NULL and length-0 values (e.g. gamm's dev.expl),
## which would otherwise collapse the whole tibble to 0 rows.
safe_val <- function(x) if (length(x) == 1 && !is.null(x)) x else NA_real_

## Build one summary row.
## gam_summary/smooth_row/model can be NULL when the model was not fitted.
make_summary_row <- function(level_name, brite_category, data,
                             gam_summary = NULL, smooth_row = NULL,
                             model = NULL, converged = "No",
                             warning_msg = NA_character_, error_msg = NA_character_) {
  tibble(
    brite_level        = level_name,
    brite_class        = brite_category,
    n_obs              = nrow(data),
    n_plots            = n_distinct(data$Plot),
    n_timepoints       = n_distinct(data$time_numeric),
    edf                = safe_val(smooth_row$edf),
    ref_df             = safe_val(smooth_row$Ref.df),
    F_statistic        = safe_val(smooth_row$F),
    p_value            = safe_val(smooth_row[["p-value"]]),
    adj_r_squared      = safe_val(gam_summary$r.sq),
    deviance_explained = safe_val(gam_summary$dev.expl),
    scale              = safe_val(gam_summary$scale),
    AIC                = if (!is.null(model)) safe_val(AIC(model$lme)) else NA_real_,
    converged          = converged,
    model_warning      = warning_msg,
    model_error        = error_msg
  )
}

## ---- Main analysis ---------------------------------------------------------
analyze_brite_level <- function(brite_table, brite_col, level_name) {
  
  sample_cols <- setdiff(names(brite_table), brite_col)
  
  model_data <- brite_table %>%
    filter(!is.na(.data[[brite_col]]), .data[[brite_col]] != "") %>%
    pivot_longer(all_of(sample_cols), names_to = "Sample_mt", values_to = "Expression") %>%
    mutate(Sample_cat = str_remove(Sample_mt, "_MT$")) %>%
    left_join(meta, by = "Sample_cat") %>%
    mutate(
      !!brite_col := as.character(.data[[brite_col]]),
      Expression  = as.numeric(Expression)
    ) %>%
    drop_na(Collection_date, Plot, time_numeric, Expression)
  
  ## Diagnostic summary of the input
  model_data %>%
    summarise(
      n_rows             = n(),
      n_brite_categories = n_distinct(.data[[brite_col]]),
      n_samples          = n_distinct(Sample_cat),
      min_date           = min(Collection_date),
      max_date           = max(Collection_date)
    ) %>%
    write_csv(paste0("Results/Step2_diagnostics/", level_name, "_input_summary_2.csv"))
  
  time_grid <- model_data %>%
    distinct(time_numeric, Collection_date) %>%
    arrange(time_numeric, Collection_date)
  
  brite_categories <- unique(model_data[[brite_col]])
  reference_plot   <- levels(model_data$Plot)[1]
  
  summary_rows    <- list()
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
    
    ## Skip categories with too little data
    if (n_distinct(d$time_numeric) < 4 || n_distinct(d$Plot) < 2) {
      summary_rows[[brite_category]] <- make_summary_row(
        level_name, brite_category, d,
        error_msg = "Insufficient distinct timepoints or plots for GAMM"
      )
      next
    }
    
    ## Fit GAMM, capturing warnings/errors
    warning_msg <- NA_character_
    error_msg   <- NA_character_
    model <- tryCatch(
      withCallingHandlers(
        gamm(
          Expression ~ s(time_numeric, k = min(10, n_distinct(d$time_numeric))),
          random = list(Plot = ~1), data = d, method = "REML"
        ),
        warning = function(w) { warning_msg <<- conditionMessage(w); invokeRestart("muffleWarning") }
      ),
      error = function(e) { error_msg <<- conditionMessage(e); NULL }
    )
    
    if (is.null(model)) {
      summary_rows[[brite_category]] <- make_summary_row(
        level_name, brite_category, d,
        warning_msg = warning_msg, error_msg = error_msg
      )
      next
    }
    
    ## Successful fit: extract stats
    gam_summary <- summary(model$gam)
    smooth_tbl  <- as.data.frame(gam_summary$s.table)
    smooth_row  <- smooth_tbl[1, , drop = FALSE]   # only one smooth term
    
    summary_rows[[brite_category]] <- make_summary_row(
      level_name, brite_category, d,
      gam_summary = gam_summary, smooth_row = smooth_row,
      model = model, converged = "Yes",
      warning_msg = warning_msg, error_msg = error_msg
    )
    
    ## Predictions on the shared time grid
    newdata <- time_grid %>%
      mutate(Plot = factor(reference_plot, levels = levels(model_data$Plot)))
    pred <- predict(model$gam, newdata = newdata, type = "response", se.fit = TRUE)
    
    prediction_list[[brite_category]] <- newdata %>%
      mutate(
        brite_level = level_name,
        brite_class = brite_category,
        fitted      = as.numeric(pred$fit),
        se          = as.numeric(pred$se.fit),
        lower       = fitted - 1.96 * se,
        upper       = fitted + 1.96 * se
      )
  }
  
  ## ---- Assemble results ----------------------------------------------------
  gam_summary_table <- bind_rows(summary_rows) %>%
    mutate(
      p_adj_fdr        = p.adjust(p_value, method = "BH"),
      significant_time = if_else(!is.na(p_adj_fdr) & p_adj_fdr < 0.05, "Yes", "No")
    ) %>%
    arrange(p_adj_fdr, p_value, brite_class)
  
  gam_predictions <- bind_rows(prediction_list) %>%
    left_join(
      gam_summary_table %>%
        select(brite_level, brite_class, p_value, p_adj_fdr, significant_time),
      by = c("brite_level", "brite_class")
    )
  
  significant_categories <- filter(gam_summary_table, significant_time == "Yes")
  
  ## ---- Write outputs -------------------------------------------------------
  write_csv(gam_summary_table,      paste0("Results/", level_name, "_GAM_summary_table_2.csv"))
  write_csv(gam_predictions,        paste0("Results/", level_name, "_GAM_fitted_trajectories_2.csv"))
  write_csv(model_data,             paste0("Results/", level_name, "_GAM_long_model_data_2.csv"))
  write_csv(significant_categories, paste0("Results/", level_name, "_GAM_significant_time_classes_2.csv"))
  
  ## ---- Plot significant classes --------------------------------------------
  if (nrow(significant_categories) > 0) {
    brite_order <- significant_categories$brite_class  # already ordered by p_adj_fdr
    
    plot_points <- model_data %>%
      rename(brite_class = all_of(brite_col)) %>%
      semi_join(significant_categories, by = "brite_class") %>%
      mutate(brite_class = factor(brite_class, levels = brite_order))
    
    plot_pred <- gam_predictions %>%
      filter(significant_time == "Yes") %>%
      mutate(brite_class = factor(brite_class, levels = brite_order))
    
    ## Split significant classes into pages of 6.
    ## Each PDF page has a fixed 2-row by 3-column facet layout.
    max_plots_per_page <- 6
    brite_pages <- split(
      brite_order,
      ceiling(seq_along(brite_order) / max_plots_per_page)
    )
    n_pages <- length(brite_pages)
    pdf_file <- paste0("Figures/", level_name, "_significant_temporal_GAMs_v2.pdf")
    
    grDevices::pdf(pdf_file, width = 14, height = 8.5, onefile = TRUE)
    for (page_i in seq_along(brite_pages)) {
      page_classes <- brite_pages[[page_i]]
      
      page_points <- plot_points %>%
        filter(brite_class %in% page_classes) %>%
        mutate(brite_class = factor(brite_class, levels = page_classes))
      
      page_pred <- plot_pred %>%
        filter(brite_class %in% page_classes) %>%
        mutate(brite_class = factor(brite_class, levels = page_classes))
      
      temporal_plot <- ggplot(page_points, aes(Collection_date, Expression)) +
        geom_point(alpha = 0.45, size = 1.1, color = "grey35") +
        geom_ribbon(data = page_pred,
                    aes(Collection_date, ymin = lower, ymax = upper),
                    inherit.aes = FALSE, fill = "#9ecae1", alpha = 0.35) +
        geom_line(data = page_pred, aes(Collection_date, fitted),
                  inherit.aes = FALSE, color = "#08519c", linewidth = 0.8) +
        facet_wrap(~ brite_class, scales = "free_y", nrow = 2, ncol = 3) +
        labs(
          title = paste0(
            level_name,
            ": BRITE classes with significant temporal GAM smooths, page ",
            page_i, " of ", n_pages
          ),
          x = "Collection date", y = "Expression"
        ) +
        theme(
          panel.grid.minor = element_blank(),
          axis.text.x      = element_text(angle = 45, hjust = 1),
          strip.text       = element_text(size = 8)
        )
      
      print(temporal_plot)
      
      ## Optional page-wise PNGs using the same 2 x 3 layout.
      ggsave(
        paste0(
          "Figures/", level_name,
          "_significant_temporal_GAMs_page_",
          sprintf("%02d", page_i), ".png"
        ),
        temporal_plot, width = 14, height = 8.5, dpi = 300
      )
    }
    grDevices::dev.off()
  }
  
  cat(level_name, "classes modeled:", nrow(gam_summary_table), "\n")
  cat(level_name, "significant temporal classes (FDR < 0.05):",
      nrow(significant_categories), "\n")
  
  list(
    model_data    = model_data,
    summary_table = gam_summary_table,
    predictions   = gam_predictions
  )
}

## ---- Run both BRITE levels -------------------------------------------------
brite_D <- read.csv("Data/brite_D_expression.csv", check.names = FALSE)
brite_C <- read.csv("Data/brite_C_expression.csv", check.names = FALSE)

results <- list(
  brite_D = analyze_brite_level(brite_D, "brite_D", "brite_D"),
  brite_C = analyze_brite_level(brite_C, "brite_C", "brite_C")
)

## ---- Combined summary ------------------------------------------------------
combined_summary <- bind_rows(
  results$brite_D$summary_table,
  results$brite_C$summary_table
)
write_csv(combined_summary, "Results/BRITE_GAM_summary_table_all_levels_2.csv")

cat("\nFiles saved:\n",
    "Results/brite_D_GAM_summary_table.csv\n",
    "Results/brite_C_GAM_summary_table.csv\n",
    "Results/BRITE_GAM_summary_table_all_levels.csv\n",
    "Results/brite_D_GAM_significant_time_classes.csv\n",
    "Results/brite_C_GAM_significant_time_classes.csv\n",
    "Figures/brite_D_significant_temporal_GAMs.pdf and page-wise PNGs\n",
    "Figures/brite_C_significant_temporal_GAMs.pdf and page-wise PNGs\n",
    sep = "")