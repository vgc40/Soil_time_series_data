# =============================================================================
# Figure 1: root biomass outlier review, additional temporal and same-date checks
# Standalone script. All observations retained in primary analyses.
# R packages: tidyverse, lubridate, mgcv, lme4, lmerTest, emmeans
# =============================================================================
rm(list = ls())

required <- c("tidyverse", "lubridate", "mgcv", "lme4", "lmerTest", "emmeans")
missing <- setdiff(required, rownames(installed.packages()))
if (length(missing)) stop("Install packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(mgcv)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

input_file <- file.path("Data", "Data_merged_cleaned.csv")
result_dir <- file.path("Results", "manuscript_final", "Figure1")
figure_dir <- file.path("Figures", "manuscript_final", "Figure1")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_file)) stop("Missing source file: ", input_file)

all_data <- readr::read_csv(input_file, show_col_types = FALSE) %>%
  filter(as.character(Irrigation_trmt_pct) == "100") %>%
  transmute(
    observation_id = paste(Sample_cat, paste0("plot", Plot), as.Date(Collection_date), sep = "_"),
    Sample_cat,
    Plot_ID = factor(Plot),
    Block = factor(Block),
    Sampling_Date = as.Date(Collection_date),
    root_biomass = as.numeric(Roots_gdw),
    total_soil_respiration = as.numeric(Resp_rate_mg_C_m2_h),
    Season = factor(if_else(month(Sampling_Date) %in% 4:10,
                            "Plant-active, Apr-Oct", "Plant-dormant, Nov-Mar"),
                    levels = c("Plant-active, Apr-Oct", "Plant-dormant, Nov-Mar"))
  ) %>%
  filter(!is.na(root_biomass), !is.na(Plot_ID), !is.na(Sampling_Date)) %>%
  arrange(Plot_ID, Sampling_Date) %>%
  mutate(time_days = as.numeric(Sampling_Date - min(Sampling_Date)),
         Sampling_Date_factor = factor(Sampling_Date))
if (anyDuplicated(all_data$observation_id)) stop("Duplicate observation IDs, review source data.")
if (any(all_data$root_biomass < 0)) warning("Negative root biomass found, review these values.")

# 1. Recalculate original within-plot and pooled screens.
plot_stats <- all_data %>% group_by(Plot_ID) %>%
  summarise(n_plot = n(), plot_median = median(root_biomass),
            plot_mad = mad(root_biomass, constant = 1),
            q1 = quantile(root_biomass, 0.25, names = FALSE),
            q3 = quantile(root_biomass, 0.75, names = FALSE),
            .groups = "drop") %>%
  mutate(plot_iqr = q3 - q1,
         iqr_lower = q1 - 1.5 * plot_iqr,
         iqr_upper = q3 + 1.5 * plot_iqr)

pooled_median <- median(all_data$root_biomass)
pooled_mad <- mad(all_data$root_biomass, constant = 1)

audit <- all_data %>% left_join(plot_stats, by = "Plot_ID") %>%
  group_by(Plot_ID) %>% arrange(Sampling_Date, .by_group = TRUE) %>%
  mutate(
    plot_modified_z = if_else(plot_mad > 0,
                              0.6745 * (root_biomass - plot_median) / plot_mad, NA_real_),
    plot_modified_z_flag = !is.na(plot_modified_z) & abs(plot_modified_z) > 3.5,
    plot_iqr_flag = root_biomass < iqr_lower | root_biomass > iqr_upper,
    previous_root = lag(root_biomass),
    next_root = lead(root_biomass),
    previous_date = lag(Sampling_Date),
    next_date = lead(Sampling_Date),
    within_plot_any_flag = plot_modified_z_flag | plot_iqr_flag
  ) %>% ungroup() %>%
  mutate(pooled_modified_z = if (pooled_mad > 0) {
           0.6745 * (root_biomass - pooled_median) / pooled_mad
         } else NA_real_,
         pooled_modified_z_flag = !is.na(pooled_modified_z) &
           abs(pooled_modified_z) > 3.5)

# 2. Plot-specific local temporal check.
# Predict each interior month by LINEAR INTERPOLATION between preceding and
# following months of the same plot, accounting for actual elapsed days.
# The difference is a descriptive local departure, NOT a formal hypothesis test.
# We deliberately do not call a point an outlier based on interpolation alone.
audit <- audit %>% group_by(Plot_ID) %>%
  arrange(Sampling_Date, .by_group = TRUE) %>%
  mutate(
    local_expected = if_else(
      !is.na(previous_root) & !is.na(next_root) &
        as.numeric(next_date - previous_date) > 0,
      previous_root + (next_root - previous_root) *
        as.numeric(Sampling_Date - previous_date) /
        as.numeric(next_date - previous_date), NA_real_
    ),
    local_deviation = root_biomass - local_expected,
    local_abs_deviation = abs(local_deviation)
  ) %>% ungroup()

# 3. Common seasonal trajectory with plot effects, fitted leaving each
# observation out in turn. This is a POOL-OF-PLOTS model that accounts for
# annual time rather than incorrectly fitting a flexible GAM to ~13 points/plot.
# The standardized deviations are predictive diagnostics, NOT calibrated tests;
# threshold 3 is a REVIEW flag only. Check sensitivity to the smooth below.
fit_seasonal <- function(training, k = 5) {
  mgcv::gam(root_biomass ~ s(time_days, k = k, bs = "cr") + Plot_ID,
            data = training, method = "REML")
}

k_time <- min(5L, n_distinct(audit$time_days) - 1L)
if (k_time < 3L) stop("Insufficient dates for seasonal smoother")
loo <- lapply(seq_len(nrow(audit)), function(i) {
  train <- audit[-i, , drop = FALSE]
  test <- audit[i, , drop = FALSE]
  result <- tryCatch({
    model <- fit_seasonal(train, k = k_time)
    estimate <- as.numeric(predict(model, newdata = test))
    residual_scale <- sqrt(mean(residuals(model)^2))
    tibble(observation_id = test$observation_id,
           temporal_loo_expected = estimate,
           temporal_loo_deviation = test$root_biomass - estimate,
           temporal_loo_train_rmse = residual_scale,
           temporal_loo_standardized = if (residual_scale > 0) {
             (test$root_biomass - estimate) / residual_scale
           } else NA_real_,
           temporal_loo_review_flag = residual_scale > 0 &&
             abs((test$root_biomass - estimate) / residual_scale) > 3,
           temporal_fit_error = NA_character_)
  }, error = function(e) tibble(
    observation_id = test$observation_id,
    temporal_loo_expected = NA_real_, temporal_loo_deviation = NA_real_,
    temporal_loo_train_rmse = NA_real_, temporal_loo_standardized = NA_real_,
    temporal_loo_review_flag = FALSE,
    temporal_fit_error = conditionMessage(e)))
  result
}) %>% bind_rows()
audit <- audit %>% left_join(loo, by = "observation_id")

# 4. Compare each root measurement to the OTHER plots on that date.
# Peer n is usually only 5, so median/MAD scores are fragile and not used
# as formal outlier tests. This screen is descriptive and context-dependent.
peer <- lapply(seq_len(nrow(audit)), function(i) {
  row <- audit[i, ]
  others <- audit %>% filter(Sampling_Date == row$Sampling_Date,
                             observation_id != row$observation_id)
  n_peer <- nrow(others)
  peer_median <- if (n_peer) median(others$root_biomass) else NA_real_
  peer_mad <- if (n_peer > 1) mad(others$root_biomass, constant = 1) else NA_real_
  tibble(observation_id = row$observation_id,
         peer_n = n_peer,
         peer_median = peer_median,
         peer_mad = peer_mad,
         peer_difference = row$root_biomass - peer_median,
         peer_modified_z = if (is.finite(peer_mad) && peer_mad > 0) {
           0.6745 * (row$root_biomass - peer_median) / peer_mad
         } else NA_real_)
}) %>% bind_rows()
audit <- audit %>% left_join(peer, by = "observation_id") %>%
  mutate(
    peer_descriptive_flag = peer_n >= 4 & !is.na(peer_modified_z) &
      abs(peer_modified_z) > 3.5,
    review_any_flag = within_plot_any_flag | pooled_modified_z_flag |
      temporal_loo_review_flag | peer_descriptive_flag,
    review_reason = paste0(
      if_else(plot_modified_z_flag, "plot_modified_z; ", ""),
      if_else(plot_iqr_flag, "plot_IQR; ", ""),
      if_else(pooled_modified_z_flag, "pooled_modified_z; ", ""),
      if_else(temporal_loo_review_flag, "temporal_LOO_diagnostic; ", ""),
      if_else(peer_descriptive_flag, "same_date_peer_diagnostic; ", "")
    )
  )

# 5. Sensitivity of root-biomass x season relationship to exclusion scenarios.
# Refit original crossed random-intercept model; if plot variance collapses
# to zero, use date-only random intercept as in original Figure 1 code.
# Keep failed/singular fits visible rather than silently switching to lm.
model_data <- audit %>% filter(is.finite(total_soil_respiration)) %>%
  mutate(Sampling_Date_factor = factor(Sampling_Date))

fit_sensitivity <- function(d, scenario) {
  if (nrow(d) < 20 || n_distinct(d$Sampling_Date) < 3) {
    return(tibble(scenario = scenario, n = nrow(d), model = NA_character_,
                  season = NA_character_, slope = NA_real_, lower_95 = NA_real_,
                  upper_95 = NA_real_, slope_p = NA_real_,
                  interaction_estimate = NA_real_, interaction_p = NA_real_,
                  singular = NA, error = "Insufficient observations"))
  }
  tryCatch({
    d <- droplevels(d)
    ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
    m <- lmer(total_soil_respiration ~ root_biomass * Season +
                (1 | Plot_ID) + (1 | Sampling_Date_factor),
              data = d, REML = TRUE, control = ctrl)
    model_name <- "plot_and_date_random_intercepts"
    v <- as.data.frame(VarCorr(m))
    plot_v <- v$vcov[v$grp == "Plot_ID"]
    date_v <- v$vcov[v$grp == "Sampling_Date_factor"]
    if (isSingular(m, tol = 1e-5) && length(plot_v) == 1 &&
        length(date_v) == 1 && plot_v < 1e-6 && date_v >= 1e-6) {
      m <- lmer(total_soil_respiration ~ root_biomass * Season +
                  (1 | Sampling_Date_factor), data = d, REML = TRUE, control = ctrl)
      model_name <- "date_random_intercept_only"
    }
    trends <- as.data.frame(summary(emmeans::emtrends(
      m, specs = ~ Season, var = "root_biomass"), infer = c(TRUE, TRUE)))
    cc <- coef(summary(m))
    int_row <- grep("root_biomass:Season", rownames(cc))
    tibble(scenario = scenario, n = nrow(d), model = model_name,
           season = as.character(trends$Season),
           slope = trends$root_biomass.trend,
           lower_95 = trends$lower.CL,
           upper_95 = trends$upper.CL,
           slope_p = trends$p.value,
           interaction_estimate = if (length(int_row) == 1) cc[int_row, "Estimate"] else NA_real_,
           interaction_p = if (length(int_row) == 1) cc[int_row, "Pr(>|t|)"] else NA_real_,
           singular = isSingular(m, tol = 1e-5), error = NA_character_)
  }, error = function(e) tibble(
    scenario = scenario, n = nrow(d), model = NA_character_,
    season = NA_character_, slope = NA_real_, lower_95 = NA_real_,
    upper_95 = NA_real_, slope_p = NA_real_,
    interaction_estimate = NA_real_, interaction_p = NA_real_,
    singular = NA, error = conditionMessage(e)))
}

scenario_flags <- list(
  all_observations = rep(FALSE, nrow(model_data)),
  exclude_plot_modified_z = model_data$plot_modified_z_flag,
  exclude_plot_IQR = model_data$plot_iqr_flag,
  exclude_either_plot_screen = model_data$within_plot_any_flag,
  exclude_both_plot_screens = model_data$plot_modified_z_flag & model_data$plot_iqr_flag,
  exclude_temporal_LOO_review = model_data$temporal_loo_review_flag,
  exclude_all_diagnostic_flags = model_data$review_any_flag
)
sensitivity <- purrr::imap_dfr(scenario_flags, function(flag, name) {
  fit_sensitivity(model_data[!flag, , drop = FALSE], name) %>%
    mutate(n_excluded = sum(flag), .after = n)
})

# Also assess influence of each individual flagged point (one at a time).
ids_to_review <- model_data$observation_id[model_data$review_any_flag]
leave_one_review <- purrr::map_dfr(ids_to_review, function(id) {
  fit_sensitivity(model_data %>% filter(observation_id != id),
                  paste0("leave_one_", id)) %>%
    mutate(excluded_observation_id = id)
})

# 6. Output. The review list includes the original and new diagnostics.
summary_by_plot <- audit %>% group_by(Plot_ID) %>%
  summarise(n = n(), plot_median = first(plot_median), plot_mad = first(plot_mad),
            n_plot_modified_z = sum(plot_modified_z_flag),
            n_plot_IQR = sum(plot_iqr_flag),
            n_temporal_LOO_review = sum(temporal_loo_review_flag),
            n_peer_descriptive = sum(peer_descriptive_flag),
            n_review_any = sum(review_any_flag), .groups = "drop")

write_csv(audit, file.path(result_dir, "Figure1_root_robust_all_observations.csv"))
write_csv(audit %>% filter(review_any_flag),
          file.path(result_dir, "Figure1_root_robust_review_candidates.csv"))
write_csv(summary_by_plot, file.path(result_dir, "Figure1_root_robust_plot_summary.csv"))
write_csv(sensitivity, file.path(result_dir, "Figure1_root_robust_sensitivity_scenarios.csv"))
write_csv(leave_one_review, file.path(result_dir, "Figure1_root_robust_leave_one_out_models.csv"))

p1 <- ggplot(audit, aes(Sampling_Date, root_biomass)) +
  geom_line(aes(group = Plot_ID), color = "grey50") +
  geom_point(aes(color = review_any_flag, shape = review_any_flag), size = 2.5) +
  facet_wrap(~Plot_ID, ncol = 3, scales = "free_y") +
  scale_color_manual(values = c(`FALSE` = "grey30", `TRUE` = "#D55E00")) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
  labs(title = "Root biomass observations for further review",
       subtitle = "Flag from any screening/diagnostic method, not exclusion decision",
       x = NULL, y = "Root biomass (g dry weight)",
       color = "Review flag", shape = "Review flag") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(figure_dir, "Figure1_root_robust_within_plot_review.png"),
       p1, width = 11, height = 7, dpi = 300)
ggsave(file.path(figure_dir, "Figure1_root_robust_within_plot_review.pdf"),
       p1, width = 11, height = 7)

p2 <- ggplot(audit, aes(temporal_loo_expected, root_biomass)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(color = temporal_loo_review_flag), size = 2) +
  facet_wrap(~Plot_ID, ncol = 3) +
  scale_color_manual(values = c(`FALSE` = "grey35", `TRUE` = "#D55E00")) +
  labs(title = "Observed root biomass vs leave-one-out seasonal expectation",
       subtitle = "Expected values from pooled seasonal GAM adjusted for plot",
       x = "Seasonal model prediction leaving sample out (g)",
       y = "Observed root biomass (g)", color = "Temporal review flag") +
  theme_bw(base_size = 11)
ggsave(file.path(figure_dir, "Figure1_root_robust_temporal_LOO.png"),
       p2, width = 11, height = 7, dpi = 300)

cat("\nPlot-level review summary:\n")
print(summary_by_plot, n = Inf)
cat("\nSensitivity scenarios:\n")
print(sensitivity, n = Inf)
cat("\nAll records RETAINED in full-data model. No automatic deletion.\n")
cat("\nNOTE: temporal LOO and same-date peer flags are exploratory diagnostics,\n",
    "NOT validated outlier tests or reasons to exclude observations by themselves.\n", sep = "")
