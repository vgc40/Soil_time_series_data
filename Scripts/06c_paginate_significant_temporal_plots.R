rm(list = ls())

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

theme_set(theme_bw(base_size = 11))

run_label <- "txt_only_significant_only"
results_dir <- file.path("Results", paste0("BRITE_temporal_", run_label))
figures_dir <- file.path("Figures", paste0("BRITE_temporal_", run_label))

dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

plot_significant_pages <- function(level_name, brite_col) {
  model_data <- read_csv(file.path(results_dir, paste0(level_name, "_GAM_long_model_data.csv")), show_col_types = FALSE) %>%
    mutate(Collection_date = as.Date(Collection_date))

  gam_predictions <- read_csv(file.path(results_dir, paste0(level_name, "_GAM_fitted_trajectories.csv")), show_col_types = FALSE) %>%
    mutate(Collection_date = as.Date(Collection_date))

  significant_categories <- read_csv(file.path(results_dir, paste0(level_name, "_GAM_significant_time_classes.csv")), show_col_types = FALSE)

  if (nrow(significant_categories) == 0) {
    write_csv(
      tibble(level_name = level_name, note = "No significant classes to plot"),
      file.path(results_dir, paste0(level_name, "_significant_plot_page_index.csv"))
    )
    return(invisible(NULL))
  }

  brite_order <- significant_categories$brite_class
  page_size <- 6

  page_index <- tibble(
    brite_class = brite_order,
    page = ceiling(seq_along(brite_order) / page_size)
  )

  write_csv(
    page_index,
    file.path(results_dir, paste0(level_name, "_significant_plot_page_index.csv"))
  )

  plot_points <- model_data %>%
    rename(brite_class = all_of(brite_col)) %>%
    semi_join(significant_categories, by = "brite_class") %>%
    mutate(brite_class = factor(brite_class, levels = brite_order))

  plot_pred <- gam_predictions %>%
    filter(significant_time == "Yes") %>%
    mutate(brite_class = factor(brite_class, levels = brite_order))

  pdf(
    file.path(figures_dir, paste0(level_name, "_significant_temporal_GAMs.pdf")),
    width = 11,
    height = 8.5,
    onefile = TRUE
  )

  for (page_number in seq_len(max(page_index$page))) {
    page_classes <- page_index %>%
      filter(page == page_number) %>%
      pull(brite_class)

    page_points <- plot_points %>%
      filter(brite_class %in% page_classes) %>%
      mutate(brite_class = factor(brite_class, levels = page_classes))

    page_pred <- plot_pred %>%
      filter(brite_class %in% page_classes) %>%
      mutate(brite_class = factor(brite_class, levels = page_classes))

    page_plot <- ggplot(page_points, aes(Collection_date, Expression)) +
      geom_point(alpha = 0.45, size = 1, color = "grey35") +
      geom_ribbon(
        data = page_pred,
        aes(Collection_date, ymin = lower, ymax = upper),
        inherit.aes = FALSE,
        fill = "#9ecae1",
        alpha = 0.35
      ) +
      geom_line(
        data = page_pred,
        aes(Collection_date, fitted),
        inherit.aes = FALSE,
        color = "#08519c",
        linewidth = 0.8
      ) +
      facet_wrap(~ brite_class, scales = "free_y", ncol = 2) +
      labs(
        title = paste0(level_name, ": significant temporal GAMs after txt-only exclusion"),
        subtitle = paste0("Page ", page_number, " of ", max(page_index$page), " | 6 panels per page"),
        x = "Collection date",
        y = "Expression"
      ) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 8)
      )

    print(page_plot)
  }

  dev.off()
}

plot_significant_pages("brite_D", "brite_D")
plot_significant_pages("brite_C", "brite_C")
