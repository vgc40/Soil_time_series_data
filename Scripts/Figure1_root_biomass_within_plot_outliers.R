# =============================================================================
# Figure 1 supplement: root biomass outlier screening WITHIN each plot
# Standalone diagnostic. Does NOT remove observations or refit Figure 1 models.
# =============================================================================
rm(list = ls())

required_packages <- c("tidyverse", "lubridate")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages)) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

input_file <- file.path("Data", "Data_merged_cleaned.csv")
result_dir <- file.path("Results", "manuscript_final", "Figure1")
figure_dir <- file.path("Figures", "manuscript_final", "Figure1")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) stop("Input file not found: ", input_file)

data <- readr::read_csv(input_file, show_col_types = FALSE) %>%
  filter(as.character(Irrigation_trmt_pct) == "100") %>%
  transmute(
    observation_id = paste(
      Sample_cat, paste0("plot", Plot), as.Date(Collection_date), sep = "_"
    ),
    Sample_cat,
    Sampling_Date = as.Date(Collection_date),
    Plot_ID = factor(Plot),
    Block = factor(Block),
    root_biomass = as.numeric(Roots_gdw),
    total_soil_respiration = as.numeric(Resp_rate_mg_C_m2_h),
    Season = if_else(
      month(Sampling_Date) %in% 4:10,
      "Plant-active, Apr-Oct", "Plant-dormant, Nov-Mar"
    )
  ) %>%
  filter(!is.na(root_biomass), !is.na(Sampling_Date), !is.na(Plot_ID)) %>%
  arrange(Plot_ID, Sampling_Date)

if (anyDuplicated(data$observation_id)) {
  warning("Some observation IDs repeat. Inspect duplicates before interpreting flags.")
}

# Modified z-score uses plot-specific median and unscaled MAD.
# Threshold > 3.5 matches your original root biomass rule.
# Also compute an exploratory plot-specific 1.5 * IQR flag (not automatic exclusion).
# With roughly 13 dates per plot, neither statistic establishes measurement error.
plot_stats <- data %>%
  group_by(Plot_ID) %>%
  summarise(
    n_plot = n(),
    plot_median = median(root_biomass),
    plot_mad = mad(root_biomass, constant = 1),
    q1 = quantile(root_biomass, 0.25, names = FALSE),
    q3 = quantile(root_biomass, 0.75, names = FALSE),
    plot_iqr = IQR(root_biomass),
    .groups = "drop"
  ) %>%
  mutate(
    iqr_lower = q1 - 1.5 * plot_iqr,
    iqr_upper = q3 + 1.5 * plot_iqr
  )

# Global modified z-score calculated on the SAME root-nonmissing observations
# to show how flags differ when the reference distribution changes.
global_median <- median(data$root_biomass)
global_mad <- mad(data$root_biomass, constant = 1)

audit <- data %>%
  left_join(plot_stats, by = "Plot_ID") %>%
  group_by(Plot_ID) %>%
  arrange(Sampling_Date, .by_group = TRUE) %>%
  mutate(
    within_plot_modified_z = if_else(
      plot_mad > 0,
      0.6745 * (root_biomass - plot_median) / plot_mad,
      NA_real_
    ),
    # If MAD is zero, leave z-score undefined, rather than dividing by zero.
    within_plot_modified_z_flag = !is.na(within_plot_modified_z) &
      abs(within_plot_modified_z) > 3.5,
    within_plot_iqr_flag = root_biomass < iqr_lower |
      root_biomass > iqr_upper,
    previous_root = lag(root_biomass),
    next_root = lead(root_biomass),
    difference_previous = root_biomass - previous_root,
    difference_next = next_root - root_biomass,
    # Visual diagnostic only: not a statistical outlier test.
    differs_from_both_neighbors = !is.na(previous_root) &
      !is.na(next_root) &
      sign(root_biomass - previous_root) == sign(root_biomass - next_root),
    any_within_plot_flag = within_plot_modified_z_flag | within_plot_iqr_flag,
    within_plot_criteria = case_when(
      within_plot_modified_z_flag & within_plot_iqr_flag ~ "modified_z; IQR",
      within_plot_modified_z_flag ~ "modified_z",
      within_plot_iqr_flag ~ "IQR",
      TRUE ~ "None"
    )
  ) %>%
  ungroup() %>%
  mutate(
    global_modified_z = if (global_mad > 0) {
      0.6745 * (root_biomass - global_median) / global_mad
    } else {
      NA_real_
    },
    global_modified_z_flag = !is.na(global_modified_z) &
      abs(global_modified_z) > 3.5,
    flag_comparison = case_when(
      global_modified_z_flag & within_plot_modified_z_flag ~ "Both modified-z screens",
      global_modified_z_flag ~ "Global modified-z only",
      within_plot_modified_z_flag ~ "Within-plot modified-z only",
      TRUE ~ "Neither modified-z screen"
    )
  )

# OPTIONAL: compare with ALL original global-screening criteria if that CSV exists.
original_file <- file.path(result_dir, "Figure1_root_biomass_outlier_review.csv")
if (file.exists(original_file)) {
  original <- readr::read_csv(original_file, show_col_types = FALSE)
  if (all(c("observation_id", "any_statistical_flag", "triggered_criteria") %in% names(original))) {
    original <- original %>%
      select(observation_id,
             original_any_statistical_flag = any_statistical_flag,
             original_triggered_criteria = triggered_criteria) %>%
      distinct(observation_id, .keep_all = TRUE)
    audit <- audit %>% left_join(original, by = "observation_id")
  }
}

summary_by_plot <- audit %>%
  group_by(Plot_ID) %>%
  summarise(
    n = n(),
    median_root = first(plot_median),
    mad_root = first(plot_mad),
    iqr_root = first(plot_iqr),
    n_within_plot_modified_z = sum(within_plot_modified_z_flag),
    n_within_plot_iqr = sum(within_plot_iqr_flag),
    n_either_within_plot = sum(any_within_plot_flag),
    n_global_modified_z = sum(global_modified_z_flag),
    .groups = "drop"
  )

flagged_observations <- audit %>%
  filter(any_within_plot_flag | global_modified_z_flag) %>%
  arrange(Plot_ID, Sampling_Date)

readr::write_csv(audit, file.path(result_dir, "Figure1_root_outliers_within_plot_all.csv"))
readr::write_csv(flagged_observations, file.path(result_dir, "Figure1_root_outliers_within_plot_flagged.csv"))
readr::write_csv(summary_by_plot, file.path(result_dir, "Figure1_root_outliers_within_plot_summary.csv"))

season_colors <- c(
  "Plant-active, Apr-Oct" = "#2E8B57",
  "Plant-dormant, Nov-Mar" = "#5B6F8C"
)
season_periods <- tibble(
  Season = c("Plant-active, Apr-Oct", "Plant-dormant, Nov-Mar", "Plant-active, Apr-Oct"),
  xmin = as.Date(c("2023-06-01", "2023-11-01", "2024-04-01")),
  xmax = as.Date(c("2023-10-31", "2024-03-31", "2024-06-30"))
)

p <- ggplot(audit, aes(Sampling_Date, root_biomass)) +
  geom_rect(
    data = season_periods,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Season),
    inherit.aes = FALSE, alpha = 0.17, color = NA
  ) +
  geom_line(aes(group = Plot_ID), color = "grey45", linewidth = 0.6) +
  geom_point(aes(color = any_within_plot_flag, shape = any_within_plot_flag), size = 2.5) +
  facet_wrap(~Plot_ID, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = season_colors) +
  scale_color_manual(values = c(`FALSE` = "grey30", `TRUE` = "#D55E00"),
                     labels = c(`FALSE` = "Not flagged", `TRUE` = "Within-plot flag")) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     labels = c(`FALSE` = "Not flagged", `TRUE` = "Within-plot flag")) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y") +
  labs(title = "Root biomass: within-plot outlier screening",
       subtitle = "Plot-specific modified z-score (>3.5) or 1.5 x IQR; flags do not imply exclusion",
       x = NULL, y = "Root biomass (g dry weight)", color = "Screening", shape = "Screening") +
  guides(fill = "none") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Fixed y scale is useful for comparing absolute biomass across plots.
p_same_y <- p + facet_wrap(~Plot_ID, ncol = 3, scales = "fixed")

ggsave(file.path(figure_dir, "Figure1_root_outliers_within_plot_free_y.png"),
       p, width = 10, height = 7, dpi = 320)
ggsave(file.path(figure_dir, "Figure1_root_outliers_within_plot_free_y.pdf"),
       p, width = 10, height = 7)
ggsave(file.path(figure_dir, "Figure1_root_outliers_within_plot_same_y.png"),
       p_same_y, width = 10, height = 7, dpi = 320)

cat("\nWithin-plot root biomass screening summary:\n")
print(summary_by_plot, n = Inf)
cat("\nFlagged observations (global or within-plot modified z, or within-plot IQR):\n")
print(flagged_observations %>%
        select(observation_id, Plot_ID, Sampling_Date, root_biomass,
               global_modified_z_flag, within_plot_modified_z_flag,
               within_plot_iqr_flag, flag_comparison), n = Inf)
cat("\nNo observations excluded. Output saved to: ", result_dir, "\n", sep = "")
