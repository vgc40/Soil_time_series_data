args <- commandArgs(trailingOnly = TRUE)
level_name <- if (length(args) > 0) args[1] else "brite_C"
write_png <- if (length(args) > 1) args[2] else "yes"
write_png <- tolower(write_png) %in% c("yes", "y", "true", "1")

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

theme_set(theme_bw(base_size = 11))

dir.create("Figures", showWarnings = FALSE, recursive = TRUE)
dir.create("Results", showWarnings = FALSE, recursive = TRUE)

model_data <- read_csv(paste0("Results/", level_name, "_GAM_long_model_data.csv"), show_col_types = FALSE) %>%
  mutate(Collection_date = as.Date(Collection_date))

gam_predictions <- read_csv(paste0("Results/", level_name, "_GAM_fitted_trajectories.csv"), show_col_types = FALSE) %>%
  mutate(Collection_date = as.Date(Collection_date))

gam_summary_table <- read_csv(paste0("Results/", level_name, "_GAM_summary_table.csv"), show_col_types = FALSE)

brite_order_all <- gam_summary_table$brite_class
page_size <- 24
page_index <- tibble(
  brite_class = brite_order_all,
  page = ceiling(seq_along(brite_order_all) / page_size)
)

total_pages <- max(page_index$page)
write_csv(page_index, paste0("Results/", level_name, "_GAM_plot_page_index.csv"))

if (level_name == "brite_D") {
  plot_points_all <- model_data %>% rename(brite_class = brite_D)
} else {
  plot_points_all <- model_data %>% rename(brite_class = brite_C)
}

plot_points_all <- plot_points_all %>%
  left_join(
    gam_summary_table %>% select(brite_class, p_value, p_adj_fdr, significant_time),
    by = "brite_class"
  )

plot_pred_all <- gam_predictions

plot_points_all$brite_class <- factor(plot_points_all$brite_class, levels = brite_order_all)
plot_pred_all$brite_class <- factor(plot_pred_all$brite_class, levels = brite_order_all)

pdf(
  paste0("Figures/", level_name, "_all_filtered_temporal_GAMs.pdf"),
  width = 14,
  height = 10,
  onefile = TRUE
)

for (page_number in seq_len(total_pages)) {
  page_classes <- page_index %>%
    filter(page == page_number) %>%
    pull(brite_class)

  page_points <- plot_points_all %>%
    filter(brite_class %in% page_classes) %>%
    mutate(brite_class = factor(brite_class, levels = page_classes))

  page_pred <- plot_pred_all %>%
    filter(brite_class %in% page_classes) %>%
    mutate(brite_class = factor(brite_class, levels = page_classes))

  page_plot <- ggplot(page_points, aes(Collection_date, Expression)) +
    geom_point(alpha = 0.4, size = 0.9, color = "grey35") +
    geom_ribbon(
      data = page_pred,
      aes(Collection_date, ymin = lower, ymax = upper),
      inherit.aes = FALSE,
      fill = "#9ecae1",
      alpha = 0.3
    ) +
    geom_line(
      data = page_pred,
      aes(Collection_date, fitted),
      inherit.aes = FALSE,
      color = "#08519c",
      linewidth = 0.7
    ) +
    facet_wrap(~ brite_class, scales = "free_y", ncol = 4) +
    labs(
      title = paste0(level_name, " filtered temporal GAMs: page ", page_number, " of ", total_pages),
      subtitle = "All remaining BRITE classes after exclusions",
      x = "Collection date",
      y = "Expression"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 7)
    )

  print(page_plot)

  if (write_png) {
    ggsave(
      paste0("Figures/", level_name, "_all_filtered_temporal_GAMs_page_", str_pad(page_number, 3, pad = "0"), ".png"),
      page_plot,
      width = 14,
      height = 10,
      dpi = 300
    )
  }

  cat(level_name, "page", page_number, "of", total_pages, "written\n")
}

dev.off()
