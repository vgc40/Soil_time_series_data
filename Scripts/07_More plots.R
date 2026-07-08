## =============================================================================
## Standalone plotting script for BRITE GAM results
## Reads pre-generated CSVs from Results/ and (re)builds figures only.
## No model fitting — safe to run/tweak repeatedly.
## Each panel is annotated with the GAM fit statistics.
## =============================================================================
rm(list = ls())

## ---- Setup -----------------------------------------------------------------
required_packages <- c("tidyverse", "readr", "ggforce")
missing <- required_packages[!required_packages %in% rownames(installed.packages())]
if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(tidyverse); library(readr); library(ggforce)
})

theme_set(theme_bw(base_size = 11))
dir.create("Figures", showWarnings = FALSE, recursive = TRUE)

## ---- Plotting settings -----------------------------------------------------
facets_per_page <- 4          # panels per page (2 x 2 grid)
facet_ncol      <- 2
facet_nrow      <- 2
page_width      <- 10
page_height     <- 8
save_pngs       <- TRUE        # set FALSE to only produce the multi-page PDF

## ---- Function: plot one BRITE level ----------------------------------------
plot_brite_level <- function(level_name) {
  
  long_file    <- paste0("Results/", level_name, "_GAM_long_model_data.csv")
  pred_file    <- paste0("Results/", level_name, "_GAM_fitted_trajectories.csv")
  summary_file <- paste0("Results/", level_name, "_GAM_summary_table.csv")
  
  ## Guard: make sure the inputs exist
  if (!file.exists(long_file) || !file.exists(pred_file) || !file.exists(summary_file)) {
    cat("Skipping", level_name, "- missing input CSVs.\n")
    return(invisible(NULL))
  }
  
  model_data  <- read_csv(long_file, show_col_types = FALSE)
  predictions <- read_csv(pred_file, show_col_types = FALSE)
  summary_tbl <- read_csv(summary_file, show_col_types = FALSE)
  
  ## The brite class column in model_data is named after the level (e.g. brite_D)
  brite_col <- level_name
  
  ## Keep only significant categories, ordered by FDR then p-value
  sig <- summary_tbl %>%
    filter(significant_time == "Yes") %>%
    arrange(p_adj_fdr, p_value, brite_class)
  
  if (nrow(sig) == 0) {
    cat(level_name, ": no significant categories to plot.\n")
    return(invisible(NULL))
  }
  
  brite_order <- sig$brite_class
  
  plot_points <- model_data %>%
    rename(brite_class = all_of(brite_col)) %>%
    filter(brite_class %in% brite_order) %>%
    mutate(
      brite_class     = factor(brite_class, levels = brite_order),
      Collection_date = as.Date(Collection_date)
    )
  
  plot_pred <- predictions %>%
    filter(significant_time == "Yes", brite_class %in% brite_order) %>%
    mutate(
      brite_class     = factor(brite_class, levels = brite_order),
      Collection_date = as.Date(Collection_date)
    )
  
  ## ---- Build a stats label for each significant category -------------------
  stat_labels <- sig %>%
    transmute(
      brite_class = factor(brite_class, levels = brite_order),
      label = sprintf(
        "edf = %.2f\nF = %.2f\np(FDR) = %.3g\nR-sq(adj) = %.3f\nDev. expl = %s",
        edf, F_statistic, p_adj_fdr, adj_r_squared,
        ifelse(is.na(deviance_explained), "NA",
               sprintf("%.1f%%", 100 * deviance_explained))
      )
    )
  
  n_pages <- ceiling(nlevels(plot_points$brite_class) / facets_per_page)
  
  base_plot <- ggplot(plot_points, aes(Collection_date, Expression)) +
    geom_point(alpha = 0.45, size = 1.1, color = "grey35") +
    geom_ribbon(data = plot_pred, aes(Collection_date, ymin = lower, ymax = upper),
                inherit.aes = FALSE, fill = "#9ecae1", alpha = 0.35) +
    geom_line(data = plot_pred, aes(Collection_date, fitted),
              inherit.aes = FALSE, color = "#08519c", linewidth = 0.8) +
    ## --- GAM stats annotation, top-left of each panel -----------------------
  geom_text(
    data = stat_labels,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.1,
    size = 2.6, color = "grey20", lineheight = 0.95
  ) +
    labs(
      title = paste0(level_name, ": significant temporal GAM smooths (with fit stats)"),
      x = "Collection date", y = "Expression"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      strip.text       = element_text(size = 8)
    )
  
  ## ---- Multi-page PDF (one page per group of 4) ----------------------------
  pdf_path <- paste0("Figures/", level_name, "_significant_temporal_GAMs.pdf")
  pdf(pdf_path, width = page_width, height = page_height)
  for (pg in seq_len(n_pages)) {
    print(
      base_plot +
        facet_wrap_paginate(~ brite_class, ncol = facet_ncol, nrow = facet_nrow,
                            scales = "free_y", page = pg)
    )
  }
  dev.off()
  
  ## ---- Optional PNGs (one file per page) -----------------------------------
  if (save_pngs) {
    for (pg in seq_len(n_pages)) {
      ggsave(
        paste0("Figures/", level_name, "_significant_temporal_GAMs_page", pg, ".png"),
        base_plot +
          facet_wrap_paginate(~ brite_class, ncol = facet_ncol, nrow = facet_nrow,
                              scales = "free_y", page = pg),
        width = page_width, height = page_height, dpi = 300
      )
    }
  }
  
  cat(level_name, ": plotted", nrow(sig), "significant categories across",
      n_pages, "page(s).\n")
  invisible(TRUE)
}

## ---- Run for both levels ---------------------------------------------------
plot_brite_level("brite_D")
plot_brite_level("brite_C")

cat("\nDone. Figures written to the Figures/ folder.\n")