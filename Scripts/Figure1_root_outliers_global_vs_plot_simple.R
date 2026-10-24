# =============================================================================
# ROOT BIOMASS OUTLIERS: ONE METHOD, TWO COMPARISONS
# Whole 100% irrigation dataset vs. within each plot
# Does not exclude any observations or change Figure 1 models.
# =============================================================================

library(tidyverse)

input_file <- file.path("Data", "Data_merged_cleaned.csv")
output_dir <- file.path("Results", "manuscript_final", "Figure1", "simple_outliers")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# One pre-specified screening criterion: absolute modified z-score > 3.5.
# MAD is the unscaled median absolute deviation (constant = 1).
cutoff <- 3.5

if (!file.exists(input_file)) stop("File not found: ", input_file)

# Match the original Figure 1 analytic sample: 100% irrigation with both
# root biomass and respiration present. To include every measured root biomass
# value instead, remove the !is.na(total_soil_respiration) condition below.
dat <- read_csv(input_file, show_col_types = FALSE) %>%
  filter(as.character(Irrigation_trmt_pct) == "100") %>%
  transmute(
    Sample_cat = as.character(Sample_cat),
    Plot_ID = as.character(Plot),
    Sampling_Date = as.Date(Collection_date),
    root_biomass = as.numeric(Roots_gdw),
    total_soil_respiration = as.numeric(Resp_rate_mg_C_m2_h)
  ) %>%
  filter(!is.na(Plot_ID), !is.na(Sampling_Date),
         !is.na(root_biomass), !is.na(total_soil_respiration)) %>%
  arrange(Plot_ID, Sampling_Date) %>%
  mutate(observation_id = paste(Sample_cat, paste0("plot", Plot_ID),
                                Sampling_Date, sep = "_"))

if (anyDuplicated(dat$observation_id)) {
  warning("Duplicated observation IDs detected; check sampling records.")
}

# Whole-dataset median and MAD (all plots and dates pooled).
whole_median <- median(dat$root_biomass)
whole_mad <- mad(dat$root_biomass, constant = 1)

# Each plot's median and MAD (all dates within that plot).
plot_reference <- dat %>%
  group_by(Plot_ID) %>%
  summarise(
    n_plot = n(),
    plot_median = median(root_biomass),
    plot_mad = mad(root_biomass, constant = 1),
    .groups = "drop"
  )

if (whole_mad == 0 || any(plot_reference$plot_mad == 0)) {
  warning("MAD is zero for the whole dataset or some plots. ",
          "Those modified z-scores will be NA, not flagged.")
}

results <- dat %>%
  left_join(plot_reference, by = "Plot_ID") %>%
  mutate(
    whole_median = whole_median,
    whole_mad = whole_mad,
    z_whole = if_else(whole_mad > 0,
                      0.6745 * (root_biomass - whole_median) / whole_mad,
                      NA_real_),
    z_plot = if_else(plot_mad > 0,
                     0.6745 * (root_biomass - plot_median) / plot_mad,
                     NA_real_),
    flagged_whole = !is.na(z_whole) & abs(z_whole) > cutoff,
    flagged_plot = !is.na(z_plot) & abs(z_plot) > cutoff,
    classification = case_when(
      flagged_whole & flagged_plot ~ "Both",
      flagged_whole ~ "Whole dataset only",
      flagged_plot ~ "Within plot only",
      TRUE ~ "Neither"
    )
  ) %>%
  select(observation_id, Sample_cat, Plot_ID, Sampling_Date, root_biomass,
         whole_median, whole_mad, plot_median, plot_mad,
         z_whole, z_plot, flagged_whole, flagged_plot, classification)

summary_counts <- results %>%
  count(classification, name = "n_observations") %>%
  complete(classification = c("Both", "Whole dataset only",
                              "Within plot only", "Neither"),
           fill = list(n_observations = 0))

plot_summary <- results %>%
  group_by(Plot_ID) %>%
  summarise(n = n(),
            n_flagged_whole = sum(flagged_whole),
            n_flagged_within_plot = sum(flagged_plot),
            n_flagged_either = sum(flagged_whole | flagged_plot),
            .groups = "drop")

write_csv(results, file.path(output_dir, "root_outliers_all.csv"))
write_csv(filter(results, flagged_whole | flagged_plot),
          file.path(output_dir, "root_outliers_flagged.csv"))
write_csv(summary_counts, file.path(output_dir, "root_outliers_summary.csv"))
write_csv(plot_summary, file.path(output_dir, "root_outliers_plot_summary.csv"))

# Visualization: each plot over time, colored by which test flagged it.
p <- ggplot(results, aes(Sampling_Date, root_biomass, group = Plot_ID)) +
  geom_line(color = "grey65") +
  geom_point(aes(color = classification), size = 2.6) +
  facet_wrap(~ Plot_ID, ncol = 3) +
  scale_color_manual(
    values = c("Neither" = "grey45", "Both" = "#C0392B",
               "Whole dataset only" = "#0072B2",
               "Within plot only" = "#E69F00"),
    breaks = c("Both", "Whole dataset only", "Within plot only", "Neither")
  ) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") +
  labs(title = "Root biomass outlier screening",
       subtitle = "Modified z-score threshold: absolute value > 3.5",
       x = NULL, y = "Root biomass (g dry weight)", color = "Flag") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(output_dir, "root_outliers_global_vs_plot.png"), p,
       width = 11, height = 7, dpi = 300)
ggsave(file.path(output_dir, "root_outliers_global_vs_plot.pdf"), p,
       width = 11, height = 7)

cat("\nObservations assessed:", nrow(results), "\n")
cat("Whole-dataset median:", whole_median, "MAD:", whole_mad, "\n\n")
print(summary_counts)
cat("\nFlagged samples:\n")
print(results %>% filter(flagged_whole | flagged_plot) %>%
        select(Sample_cat, Plot_ID, Sampling_Date, root_biomass,
               z_whole, z_plot, classification), n = Inf)
cat("\nSaved results to: ", output_dir, "\n", sep = "")
