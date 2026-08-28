# =============================================================================
# Figure 1: field seasonal context and root-respiration relationship
# =============================================================================
#
# Uses only the 100% irrigation treatment. Statistical outlier flags prompt
# review and sensitivity analysis; observations are removed only when an
# independent measurement or data-processing error is verified.

rm(list = ls())

required_packages <- c(
  "tidyverse", "lubridate", "lme4", "lmerTest", "emmeans",
  "performance", "MASS", "patchwork"
)
missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]
if (length(missing_packages) > 0) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(performance)
  library(patchwork)
})

set.seed(42)

season_colors <- c(
  "Plant-active, Apr-Oct" = "#2E8B57",
  "Plant-dormant, Nov-Mar" = "#5B6F8C"
)
figure_dir <- file.path("Figures", "manuscript_final", "Figure1")
result_dir <- file.path("Results", "manuscript_final", "Figure1")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

theme_manuscript <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.key = element_blank()
  )

field_data <- read_csv(
  file.path("Data", "Data_merged_cleaned.csv"),
  show_col_types = FALSE
) %>%
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
    Season = factor(
      if_else(
        month(Sampling_Date) %in% 4:10,
        "Plant-active, Apr-Oct",
        "Plant-dormant, Nov-Mar"
      ),
      levels = names(season_colors)
    )
  ) %>%
  arrange(Sampling_Date, Plot_ID) %>%
  mutate(
    Sampling_Date_factor = factor(Sampling_Date),
    sequential_month = dense_rank(Sampling_Date)
  )

duplicate_sample_ids <- field_data %>%
  count(Sample_cat, name = "record_count") %>%
  filter(record_count > 1) %>%
  pull(Sample_cat)

analysis_data <- field_data %>%
  filter(!is.na(root_biomass), !is.na(total_soil_respiration))

monthly_respiration <- analysis_data %>%
  group_by(Sampling_Date, Season) %>%
  summarise(
    n = n(),
    mean_respiration = mean(total_soil_respiration),
    se_respiration = sd(total_soil_respiration) / sqrt(n),
    ci95_lower = mean_respiration - qt(0.975, n - 1) * se_respiration,
    ci95_upper = mean_respiration + qt(0.975, n - 1) * se_respiration,
    .groups = "drop"
  )

season_periods <- tribble(
  ~Season, ~xmin, ~xmax,
  "Plant-active, Apr-Oct", as.Date("2023-06-01"), as.Date("2023-10-31"),
  "Plant-dormant, Nov-Mar", as.Date("2023-11-01"), as.Date("2024-03-31"),
  "Plant-active, Apr-Oct", as.Date("2024-04-01"), as.Date("2024-06-30")
) %>%
  mutate(Season = factor(Season, levels = names(season_colors)))

panel_a <- ggplot(
  analysis_data,
  aes(Sampling_Date, total_soil_respiration)
) +
  geom_rect(
    data = season_periods,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Season),
    inherit.aes = FALSE,
    alpha = 0.20,
    color = NA
  ) +
  geom_point(
    aes(shape = Plot_ID),
    color = "grey35",
    alpha = 0.65,
    size = 2,
    position = position_jitter(width = 1.8, height = 0, seed = 42)
  ) +
  geom_ribbon(
    data = monthly_respiration,
    aes(x = Sampling_Date, ymin = ci95_lower, ymax = ci95_upper),
    inherit.aes = FALSE,
    fill = "grey55",
    alpha = 0.25
  ) +
  geom_line(
    data = monthly_respiration,
    aes(x = Sampling_Date, y = mean_respiration),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.8
  ) +
  geom_point(
    data = monthly_respiration,
    aes(x = Sampling_Date, y = mean_respiration),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    size = 3
  ) +
  scale_fill_manual(values = season_colors, drop = FALSE) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b\n%Y") +
  labs(
    title = "A",
    subtitle = "Points are plots; line and ribbon are monthly mean and 95% CI",
    x = NULL,
    y = expression(
      "Total soil respiration (mg CO"[2] * "-C " * m^{-2} * " " * h^{-1} * ")"
    ),
    shape = "Plot"
  ) +
  guides(fill = "none", shape = guide_legend(nrow = 2)) +
  theme_manuscript +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

# Exact harvest and individual irrigation dates are not marked because no
# repository event-date table was found that verifies those dates.

# ---- Combined interaction model --------------------------------------------

crossed_model <- lmer(
  total_soil_respiration ~ root_biomass * Season +
    (1 | Plot_ID) + (1 | Sampling_Date_factor),
  data = analysis_data,
  REML = TRUE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 100000)
  )
)

crossed_variances <- as.data.frame(VarCorr(crossed_model)) %>%
  as_tibble()
plot_variance <- crossed_variances %>%
  filter(grp == "Plot_ID") %>%
  pull(vcov)
date_variance <- crossed_variances %>%
  filter(grp == "Sampling_Date_factor") %>%
  pull(vcov)

primary_model <- crossed_model
primary_model_formula <- paste(
  "total_soil_respiration ~ root_biomass * Season +",
  "(1 | Plot_ID) + (1 | Sampling_Date_factor)"
)
model_simplification <- "None"

if (
  isSingular(crossed_model, tol = 1e-5) &&
    plot_variance < 1e-6 &&
    date_variance >= 1e-6
) {
  primary_model <- lmer(
    total_soil_respiration ~ root_biomass * Season +
      (1 | Sampling_Date_factor),
    data = analysis_data,
    REML = TRUE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 100000)
    )
  )
  primary_model_formula <- paste(
    "total_soil_respiration ~ root_biomass * Season +",
    "(1 | Sampling_Date_factor)"
  )
  model_simplification <- paste(
    "The requested crossed model was singular because Plot_ID variance was",
    "estimated as zero. The zero-variance Plot_ID term was removed and the",
    "sampling-date random intercept was retained."
  )
}

if (isSingular(primary_model, tol = 1e-5)) {
  stop(
    "The simplified mixed model remains singular; review it rather than ",
    "silently falling back to a simple linear regression."
  )
}

season_slopes <- emtrends(
  primary_model,
  specs = ~ Season,
  var = "root_biomass"
) %>%
  summary(infer = c(TRUE, TRUE)) %>%
  as_tibble() %>%
  rename(slope = root_biomass.trend, p_value = p.value)

fixed_effects <- as.data.frame(summary(primary_model)$coefficients) %>%
  rownames_to_column("term") %>%
  as_tibble() %>%
  rename(
    estimate = Estimate,
    std_error = `Std. Error`,
    statistic = `t value`,
    p_value = `Pr(>|t|)`
  )

fixed_effect_ci <- as.data.frame(
  confint(primary_model, parm = "beta_", method = "Wald")
) %>%
  rownames_to_column("term") %>%
  as_tibble() %>%
  rename(conf_low = `2.5 %`, conf_high = `97.5 %`)

fixed_effects <- fixed_effects %>%
  left_join(fixed_effect_ci, by = "term")

interaction_result <- fixed_effects %>%
  filter(str_detect(term, "root_biomass:Season"))
model_r2 <- performance::r2_nakagawa(primary_model) %>%
  as.data.frame() %>%
  as_tibble()

model_structure <- tibble(
  model = c("Requested crossed model", "Primary model"),
  formula = c(
    paste(
      "total_soil_respiration ~ root_biomass * Season +",
      "(1 | Plot_ID) + (1 | Sampling_Date_factor)"
    ),
    primary_model_formula
  ),
  singular = c(
    isSingular(crossed_model, tol = 1e-5),
    isSingular(primary_model, tol = 1e-5)
  ),
  model_change = c("Initial requested structure", model_simplification),
  n_observations = nrow(analysis_data),
  n_plots = n_distinct(analysis_data$Plot_ID),
  n_dates = n_distinct(analysis_data$Sampling_Date)
)

write_csv(
  crossed_variances,
  file.path(result_dir, "Figure1_initial_random_effect_variances.csv")
)
write_csv(
  model_structure,
  file.path(result_dir, "Figure1_model_structure_summary.csv")
)
write_csv(
  fixed_effects,
  file.path(result_dir, "Figure1_primary_fixed_effects.csv")
)
write_csv(
  season_slopes,
  file.path(result_dir, "Figure1_seasonal_slopes.csv")
)
write_csv(
  model_r2,
  file.path(result_dir, "Figure1_marginal_conditional_R2.csv")
)

# ---- Root-biomass outlier and influence review -----------------------------

root_median <- median(analysis_data$root_biomass)
root_mad <- mad(analysis_data$root_biomass, constant = 1)
screening_lm <- lm(
  total_soil_respiration ~ root_biomass * Season,
  data = analysis_data
)

analysis_data <- analysis_data %>%
  mutate(
    modified_z_root = 0.6745 * (root_biomass - root_median) / root_mad,
    standardized_mixed_residual = residuals(primary_model) / sigma(primary_model),
    studentized_screening_residual = rstudent(screening_lm),
    cooks_distance = cooks.distance(screening_lm),
    leverage = hatvalues(screening_lm),
    bonferroni_outlier_p = pmin(
      1,
      2 * pt(
        abs(studentized_screening_residual),
        df = df.residual(screening_lm) - 1,
        lower.tail = FALSE
      ) * nrow(analysis_data)
    )
  ) %>%
  group_by(Plot_ID) %>%
  arrange(Sampling_Date, .by_group = TRUE) %>%
  mutate(
    root_change_from_previous = root_biomass - lag(root_biomass),
    root_change_to_next = lead(root_biomass) - root_biomass
  ) %>%
  ungroup() %>%
  mutate(
    modified_z_flag = abs(modified_z_root) > 3.5,
    mixed_residual_flag = abs(standardized_mixed_residual) > 3,
    studentized_residual_flag = abs(studentized_screening_residual) > 3,
    cooks_distance_threshold = 4 / nrow(analysis_data),
    cooks_distance_flag = cooks_distance > cooks_distance_threshold,
    leverage_threshold = 2 * length(coef(screening_lm)) / nrow(analysis_data),
    leverage_flag = leverage > leverage_threshold,
    bonferroni_outlier_flag = bonferroni_outlier_p < 0.05,
    impossible_or_negative_flag = root_biomass < 0,
    duplicate_record_flag = Sample_cat %in% duplicate_sample_ids,
    any_statistical_flag = modified_z_flag |
      mixed_residual_flag |
      studentized_residual_flag |
      cooks_distance_flag |
      leverage_flag |
      bonferroni_outlier_flag
  ) %>%
  rowwise() %>%
  mutate(
    triggered_criteria = paste(
      c(
        "modified_z_root_gt_3.5",
        "absolute_mixed_residual_gt_3",
        "absolute_studentized_residual_gt_3",
        "cooks_distance_gt_4_over_n",
        "leverage_gt_2p_over_n",
        "bonferroni_outlier_test"
      )[c_across(c(
        modified_z_flag,
        mixed_residual_flag,
        studentized_residual_flag,
        cooks_distance_flag,
        leverage_flag,
        bonferroni_outlier_flag
      ))],
      collapse = "; "
    ),
    triggered_criteria = if_else(
      triggered_criteria == "",
      "None",
      triggered_criteria
    )
  ) %>%
  ungroup() %>%
  mutate(
    unit_consistency_review = paste(
      "All values use the repository column Roots_gdw;",
      "no conflicting unit field was found."
    ),
    sample_match_review = if_else(
      is.na(Sample_cat) | is.na(Plot_ID) | is.na(Sampling_Date),
      "Potential sample or plot mismatch",
      "Sample, plot, block, and date metadata are present"
    ),
    laboratory_or_field_note_review = paste(
      "No repository note independently verifies measurement failure",
      "for this observation."
    ),
    verified_measurement_or_processing_problem = FALSE,
    final_decision = "Retain",
    written_justification = if_else(
      any_statistical_flag,
      paste(
        "Retained: statistically unusual, but no independent measurement,",
        "unit, transcription, or sample-identification error was verified."
      ),
      "Retained: no screening flag or independent error evidence."
    )
  )

outlier_review <- analysis_data %>%
  select(
    observation_id, Sample_cat, Plot = Plot_ID, Block, Sampling_Date, Season,
    root_biomass, total_soil_respiration, modified_z_root,
    standardized_mixed_residual, studentized_screening_residual,
    cooks_distance, cooks_distance_threshold, leverage, leverage_threshold,
    bonferroni_outlier_p, modified_z_flag, mixed_residual_flag,
    studentized_residual_flag, cooks_distance_flag, leverage_flag,
    bonferroni_outlier_flag, impossible_or_negative_flag,
    duplicate_record_flag, root_change_from_previous, root_change_to_next,
    triggered_criteria, unit_consistency_review, sample_match_review,
    laboratory_or_field_note_review,
    verified_measurement_or_processing_problem, final_decision,
    written_justification
  )

write_csv(
  outlier_review,
  file.path(result_dir, "Figure1_root_biomass_outlier_review.csv")
)

root_trajectory_plot <- ggplot(
  analysis_data,
  aes(Sampling_Date, root_biomass, group = Plot_ID)
) +
  geom_rect(
    data = season_periods,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Season),
    inherit.aes = FALSE,
    alpha = 0.20,
    color = NA
  ) +
  geom_line(color = "grey45", linewidth = 0.6) +
  geom_point(
    aes(color = any_statistical_flag, shape = any_statistical_flag),
    size = 2.5
  ) +
  facet_wrap(~ Plot_ID, ncol = 3) +
  scale_fill_manual(values = season_colors) +
  scale_color_manual(
    values = c(`FALSE` = "grey30", `TRUE` = "#D55E00"),
    labels = c(`FALSE` = "Not flagged", `TRUE` = "Candidate flag")
  ) +
  scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 17),
    labels = c(`FALSE` = "Not flagged", `TRUE` = "Candidate flag")
  ) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y") +
  labs(
    title = "Root biomass trajectories used for outlier review",
    subtitle = "Flags prompt investigation and do not imply exclusion",
    x = NULL,
    y = "Root biomass (g dry weight)",
    color = "Screening result",
    shape = "Screening result"
  ) +
  guides(fill = "none") +
  theme_manuscript +
  theme(legend.position = "bottom")

ggsave(
  file.path(figure_dir, "Figure1_root_biomass_trajectory_diagnostics.png"),
  root_trajectory_plot,
  width = 10,
  height = 7,
  dpi = 320
)
ggsave(
  file.path(figure_dir, "Figure1_root_biomass_trajectory_diagnostics.pdf"),
  root_trajectory_plot,
  width = 10,
  height = 7
)

# ---- Leave-one-flag sensitivity models -------------------------------------

flagged_ids <- analysis_data %>%
  filter(any_statistical_flag) %>%
  pull(observation_id)

sensitivity_results <- tibble()

for (excluded_id in flagged_ids) {
  sensitivity_data <- analysis_data %>%
    filter(observation_id != excluded_id)

  sensitivity_model <- lmer(
    total_soil_respiration ~ root_biomass * Season +
      (1 | Sampling_Date_factor),
    data = sensitivity_data,
    REML = TRUE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 100000)
    )
  )

  sensitivity_slopes <- emtrends(
    sensitivity_model,
    specs = ~ Season,
    var = "root_biomass"
  ) %>%
    summary(infer = c(TRUE, TRUE)) %>%
    as_tibble()

  sensitivity_coefficients <- as.data.frame(
    summary(sensitivity_model)$coefficients
  ) %>%
    rownames_to_column("term") %>%
    as_tibble()

  sensitivity_ci <- as.data.frame(
    confint(sensitivity_model, parm = "beta_", method = "Wald")
  ) %>%
    rownames_to_column("term") %>%
    as_tibble()

  sensitivity_interaction <- sensitivity_coefficients %>%
    filter(str_detect(term, "root_biomass:Season")) %>%
    left_join(sensitivity_ci, by = "term") %>%
    transmute(
      estimate = Estimate,
      conf_low = `2.5 %`,
      conf_high = `97.5 %`,
      p_value = `Pr(>|t|)`
    )

  sensitivity_results <- bind_rows(
    sensitivity_results,
    sensitivity_slopes %>%
      transmute(
        model = paste0("Exclude_", excluded_id),
        excluded_observation_id = excluded_id,
        Season,
        slope = root_biomass.trend,
        slope_se = SE,
        slope_df = df,
        slope_lower_95 = lower.CL,
        slope_upper_95 = upper.CL,
        slope_p_value = p.value,
        interaction_estimate = sensitivity_interaction$estimate,
        interaction_lower_95 = sensitivity_interaction$conf_low,
        interaction_upper_95 = sensitivity_interaction$conf_high,
        interaction_p_value = sensitivity_interaction$p_value,
        n_observations = nrow(sensitivity_data),
        n_plots = n_distinct(sensitivity_data$Plot_ID),
        singular = isSingular(sensitivity_model, tol = 1e-5)
      )
  )
}

full_data_comparison <- season_slopes %>%
  transmute(
    model = "Full_valid_data",
    excluded_observation_id = NA_character_,
    Season,
    slope,
    slope_se = SE,
    slope_df = df,
    slope_lower_95 = lower.CL,
    slope_upper_95 = upper.CL,
    slope_p_value = p_value,
    interaction_estimate = interaction_result$estimate,
    interaction_lower_95 = interaction_result$conf_low,
    interaction_upper_95 = interaction_result$conf_high,
    interaction_p_value = interaction_result$p_value,
    n_observations = nrow(analysis_data),
    n_plots = n_distinct(analysis_data$Plot_ID),
    singular = isSingular(primary_model, tol = 1e-5)
  )

# No observation is independently verified as erroneous. Consequently the
# justified-error sensitivity model is identical to the complete-data model.
verified_error_sensitivity <- full_data_comparison %>%
  mutate(
    model = "Verified_errors_excluded",
    excluded_observation_id = "None"
  )

model_comparison <- bind_rows(
  full_data_comparison,
  verified_error_sensitivity,
  sensitivity_results
)
write_csv(
  model_comparison,
  file.path(result_dir, "Figure1_full_and_sensitivity_model_comparison.csv")
)

# Robust marginal sensitivity; it does not replace the repeated-measures model.
robust_model <- MASS::rlm(
  total_soil_respiration ~ root_biomass * Season,
  data = analysis_data,
  maxit = 100
)
robust_coefficients <- summary(robust_model)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  as_tibble() %>%
  rename(
    estimate = Value,
    standard_error = `Std. Error`,
    statistic = `t value`
  ) %>%
  mutate(
    lower_95 = estimate - 1.96 * standard_error,
    upper_95 = estimate + 1.96 * standard_error,
    analysis_note = paste(
      "Marginal robust-regression sensitivity only;",
      "primary inference uses the mixed model."
    )
  )
write_csv(
  robust_coefficients,
  file.path(result_dir, "Figure1_robust_regression_sensitivity.csv")
)

# ---- Figure 1B predictions from the combined model -------------------------

root_prediction_grid <- crossing(
  root_biomass = seq(
    min(analysis_data$root_biomass),
    max(analysis_data$root_biomass),
    length.out = 200
  ),
  Season = factor(names(season_colors), levels = names(season_colors))
)
fixed_design_matrix <- model.matrix(
  ~ root_biomass * Season,
  data = root_prediction_grid
)
fixed_vcov <- as.matrix(vcov(primary_model))
fixed_beta <- fixef(primary_model)

root_prediction_grid <- root_prediction_grid %>%
  mutate(
    fitted = as.numeric(fixed_design_matrix %*% fixed_beta),
    prediction_se = sqrt(
      diag(fixed_design_matrix %*% fixed_vcov %*% t(fixed_design_matrix))
    ),
    lower_95 = fitted - 1.96 * prediction_se,
    upper_95 = fitted + 1.96 * prediction_se
  )


# ------------------------------------------------------------
# PANEL B: Root biomass vs total soil respiration by season
# ------------------------------------------------------------
# ------------------------------------------------------------
# PANEL B: Root biomass vs total soil respiration by season
# Raw data only, with model statistics
# ------------------------------------------------------------

# Create labels for season-specific slope statistics
# Plant-active label goes bottom-right
# Plant-dormant label stays top-left
slope_labels <- season_slopes %>%
  mutate(
    label = sprintf(
      "Slope = %.2f\n95%% CI: %.2f to %.2f\np = %.3f",
      slope, lower.CL, upper.CL, p_value
    ),
    
    x = case_when(
      Season == "Plant-active, Apr-Oct" ~ Inf,
      TRUE ~ -Inf
    ),
    
    y = case_when(
      Season == "Plant-active, Apr-Oct" ~ -Inf,
      TRUE ~ Inf
    ),
    
    hjust = case_when(
      Season == "Plant-active, Apr-Oct" ~ 1.05,
      TRUE ~ -0.05
    ),
    
    vjust = case_when(
      Season == "Plant-active, Apr-Oct" ~ -0.25,
      TRUE ~ 1.1
    )
  )

# Make Panel B
panel_b <- ggplot(
  analysis_data,
  aes(
    x = root_biomass,
    y = total_soil_respiration
   # color = as.factor(analysis_data$Sampling_Date),
    
  )
) +
  
  # Observed data only
  geom_point(
    aes(shape = Plot_ID),
    color = "grey35",
    size = 2.4,
    alpha = 0.75
  ) +
  
  # Season-specific model statistics
  geom_text(
    data = slope_labels,
    aes(
      x = x,
      y = y,
      label = label,
      hjust = hjust,
      vjust = vjust
    ),
    inherit.aes = FALSE,
    size = 3.2
  ) +
  
  # One season per panel
  # Allow independent x-axis ranges
  facet_wrap(
    ~ Season,
    ncol = 1,
    scales = "free_x"
  ) +
  
  labs(
    title = "B",
    
    subtitle = paste0(
      "Root biomass × season interaction: p = ",
      format.pval(
        interaction_result$p_value,
        digits = 2
      )
    ),
    
    x = "Root biomass (g dry weight)",
    
    y = expression(
      "Total soil respiration (mg CO"[2] *
        "-C " * m^{-2} * " " * h^{-1} * ")"
    )
  ) +
  
  guides(
    shape = "none"
  ) +
  
  theme_manuscript +
  
  theme(
    strip.text = element_text(face = "bold")
  )
# Display plot
panel_b

figure_1 <- panel_a + panel_b +
  plot_layout(widths = c(1.18, 1), guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  file.path(figure_dir, "Figure1_main.png"),
  figure_1,
  width = 13,
  height = 6.2,
  dpi = 320
)
ggsave(
  file.path(figure_dir, "Figure1_main.pdf"),
  figure_1,
  width = 13,
  height = 6.2
)
ggsave(
  file.path(figure_dir, "Figure1A_field_respiration_time_series.png"),
  panel_a,
  width = 7.3,
  height = 6,
  dpi = 320
)
ggsave(
  file.path(figure_dir, "Figure1A_field_respiration_time_series.pdf"),
  panel_a,
  width = 7.3,
  height = 6
)
ggsave(
  file.path(figure_dir, "Figure1B_root_biomass_vs_respiration.png"),
  panel_b,
  width = 6.5,
  height = 6,
  dpi = 320
)
ggsave(
  file.path(figure_dir, "Figure1B_root_biomass_vs_respiration.pdf"),
  panel_b,
  width = 6.5,
  height = 6
)

# ---- Model diagnostics and final summary -----------------------------------

diagnostic_data <- analysis_data %>%
  transmute(
    observation_id,
    fitted = fitted(primary_model),
    residual = residuals(primary_model),
    standardized_residual = residuals(primary_model) / sigma(primary_model)
  )
shapiro_result <- shapiro.test(diagnostic_data$residual)

diagnostic_summary <- tibble(
  diagnostic = c(
    "Primary model singular",
    "Maximum absolute standardized residual",
    "Shapiro-Wilk residual W",
    "Shapiro-Wilk residual p-value",
    "Correlation of absolute residuals with fitted values",
    "Initial Plot_ID random-intercept variance",
    "Initial sampling-date random-intercept variance",
    "Candidate outliers",
    "Independently verified errors excluded"
  ),
  value = c(
    isSingular(primary_model, tol = 1e-5),
    max(abs(diagnostic_data$standardized_residual)),
    unname(shapiro_result$statistic),
    shapiro_result$p.value,
    cor(abs(diagnostic_data$residual), diagnostic_data$fitted),
    plot_variance,
    date_variance,
    length(flagged_ids),
    0
  )
)

write_csv(
  diagnostic_data,
  file.path(result_dir, "Figure1_primary_model_residuals.csv")
)
write_csv(
  diagnostic_summary,
  file.path(result_dir, "Figure1_model_diagnostics_summary.csv")
)
write_csv(
  monthly_respiration,
  file.path(result_dir, "Figure1A_monthly_mean_95CI.csv")
)

residual_plot <- ggplot(
  diagnostic_data,
  aes(fitted, standardized_residual)
) +
  geom_hline(
    yintercept = c(-3, 0, 3),
    linetype = c("dashed", "solid", "dashed"),
    color = c("#D55E00", "grey50", "#D55E00")
  ) +
  geom_point(color = "grey30", alpha = 0.75) +
  labs(
    title = "Residuals versus fitted values",
    x = "Fitted total soil respiration",
    y = "Standardized conditional residual"
  ) +
  theme_manuscript

qq_plot <- ggplot(
  diagnostic_data,
  aes(sample = standardized_residual)
) +
  stat_qq(color = "grey30") +
  stat_qq_line(color = "#0072B2") +
  labs(
    title = "Normal Q-Q plot",
    x = "Theoretical quantile",
    y = "Standardized conditional residual"
  ) +
  theme_manuscript

diagnostic_plot <- residual_plot + qq_plot
ggsave(
  file.path(figure_dir, "Figure1_model_diagnostics.png"),
  diagnostic_plot,
  width = 10,
  height = 4.8,
  dpi = 320
)
ggsave(
  file.path(figure_dir, "Figure1_model_diagnostics.pdf"),
  diagnostic_plot,
  width = 10,
  height = 4.8
)

analysis_summary <- tibble(
  item = c(
    "Observations in primary model",
    "Plots in primary model",
    "Sampling dates in primary model",
    "Plant-active observations",
    "Plant-dormant observations",
    "Candidate outlier observations",
    "Observations excluded from primary model",
    "Primary model formula",
    "Model simplification",
    "Root biomass x Season estimate",
    "Root biomass x Season p-value",
    "Marginal R2",
    "Conditional R2",
    "Harvest or irrigation event markers"
  ),
  value = c(
    nrow(analysis_data),
    n_distinct(analysis_data$Plot_ID),
    n_distinct(analysis_data$Sampling_Date),
    sum(analysis_data$Season == "Plant-active, Apr-Oct"),
    sum(analysis_data$Season == "Plant-dormant, Nov-Mar"),
    length(flagged_ids),
    0,
    primary_model_formula,
    model_simplification,
    interaction_result$estimate,
    interaction_result$p_value,
    model_r2$R2_marginal,
    model_r2$R2_conditional,
    paste(
      "Not shown: dates could not be verified from a repository",
      "harvest or irrigation event-date table."
    )
  )
)
write_csv(
  analysis_summary,
  file.path(result_dir, "Figure1_analysis_summary.csv")
)

cat("\nFigure 1 analysis complete.\n")
print(analysis_summary)
cat("\nSeasonal slopes:\n")
print(season_slopes)
cat("\nCandidate outliers retained for primary analysis:\n")
print(
  outlier_review %>%
    filter(
      modified_z_flag |
        mixed_residual_flag |
        studentized_residual_flag |
        cooks_distance_flag |
        leverage_flag |
        bonferroni_outlier_flag
    ) %>%
    select(
      observation_id,
      root_biomass,
      total_soil_respiration,
      triggered_criteria,
      final_decision
    )
)
