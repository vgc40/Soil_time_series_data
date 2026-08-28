# =============================================================================
# Figure 1 SI sensitivity analysis
# Season-specific models + log1p sensitivity
# =============================================================================

rm(list = ls())

library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)
library(emmeans)
library(patchwork)

result_dir <- file.path("Results", "manuscript_final", "Figure1_SI")
figure_dir <- file.path("Figures", "manuscript_final", "Figure1_SI")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

theme_manuscript <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------

analysis_data <- read_csv(
  file.path("Data", "Data_merged_cleaned.csv"),
  show_col_types = FALSE
) %>%
  filter(as.character(Irrigation_trmt_pct) == "100") %>%
  transmute(
    Sampling_Date = as.Date(Collection_date),
    Plot_ID = factor(Plot),
    root_biomass = as.numeric(Roots_gdw),
    total_soil_respiration = as.numeric(Resp_rate_mg_C_m2_h),
    Season = factor(
      if_else(
        month(Sampling_Date) %in% 4:10,
        "Plant-active, Apr-Oct",
        "Plant-dormant, Nov-Mar"
      ),
      levels = c(
        "Plant-active, Apr-Oct",
        "Plant-dormant, Nov-Mar"
      )
    ),
    Sampling_Date_factor = factor(Sampling_Date)
  ) %>%
  filter(
    !is.na(root_biomass),
    !is.na(total_soil_respiration)
  )

active_data <- analysis_data %>%
  filter(Season == "Plant-active, Apr-Oct")

dormant_data <- analysis_data %>%
  filter(Season == "Plant-dormant, Nov-Mar")

# -----------------------------------------------------------------------------
# 1. Separate seasonal models, raw respiration
# -----------------------------------------------------------------------------

active_raw <- lmer(
  total_soil_respiration ~ root_biomass +
    (1 | Sampling_Date_factor),
  data = active_data,
  REML = TRUE
)

dormant_raw <- lmer(
  total_soil_respiration ~ root_biomass +
    (1 | Sampling_Date_factor),
  data = dormant_data,
  REML = TRUE
)

# -----------------------------------------------------------------------------
# 2. Separate seasonal models, log1p respiration
# -----------------------------------------------------------------------------

active_log <- lmer(
  log1p(total_soil_respiration) ~ root_biomass +
    (1 | Sampling_Date_factor),
  data = active_data,
  REML = TRUE
)

dormant_log <- lmer(
  log1p(total_soil_respiration) ~ root_biomass +
    (1 | Sampling_Date_factor),
  data = dormant_data,
  REML = TRUE
)

# -----------------------------------------------------------------------------
# 3. Extract slopes and confidence intervals
# -----------------------------------------------------------------------------

active_raw_coef <- summary(active_raw)$coefficients["root_biomass", ]
active_log_coef <- summary(active_log)$coefficients["root_biomass", ]
dormant_raw_coef <- summary(dormant_raw)$coefficients["root_biomass", ]
dormant_log_coef <- summary(dormant_log)$coefficients["root_biomass", ]

active_raw_ci <- confint(
  active_raw,
  parm = "root_biomass",
  method = "Wald"
)

active_log_ci <- confint(
  active_log,
  parm = "root_biomass",
  method = "Wald"
)

dormant_raw_ci <- confint(
  dormant_raw,
  parm = "root_biomass",
  method = "Wald"
)

dormant_log_ci <- confint(
  dormant_log,
  parm = "root_biomass",
  method = "Wald"
)

season_model_results <- tibble(
  Season = c(
    "Plant-active",
    "Plant-active",
    "Plant-dormant",
    "Plant-dormant"
  ),
  response = c(
    "raw",
    "log1p",
    "raw",
    "log1p"
  ),
  n = c(
    nrow(active_data),
    nrow(active_data),
    nrow(dormant_data),
    nrow(dormant_data)
  ),
  slope = c(
    active_raw_coef["Estimate"],
    active_log_coef["Estimate"],
    dormant_raw_coef["Estimate"],
    dormant_log_coef["Estimate"]
  ),
  SE = c(
    active_raw_coef["Std. Error"],
    active_log_coef["Std. Error"],
    dormant_raw_coef["Std. Error"],
    dormant_log_coef["Std. Error"]
  ),
  p_value = c(
    active_raw_coef["Pr(>|t|)"],
    active_log_coef["Pr(>|t|)"],
    dormant_raw_coef["Pr(>|t|)"],
    dormant_log_coef["Pr(>|t|)"]
  ),
  lower_95 = c(
    active_raw_ci[1],
    active_log_ci[1],
    dormant_raw_ci[1],
    dormant_log_ci[1]
  ),
  upper_95 = c(
    active_raw_ci[2],
    active_log_ci[2],
    dormant_raw_ci[2],
    dormant_log_ci[2]
  )
)

write_csv(
  season_model_results,
  file.path(
    result_dir,
    "FigureS1_season_specific_model_results.csv"
  )
)

print(season_model_results)

# -----------------------------------------------------------------------------
# 4. Combined interaction models
# -----------------------------------------------------------------------------

combined_raw <- lmer(
  total_soil_respiration ~ root_biomass * Season +
    (1 | Sampling_Date_factor),
  data = analysis_data,
  REML = TRUE
)

combined_log <- lmer(
  log1p(total_soil_respiration) ~ root_biomass * Season +
    (1 | Sampling_Date_factor),
  data = analysis_data,
  REML = TRUE
)

raw_slopes <- emtrends(
  combined_raw,
  ~ Season,
  var = "root_biomass"
) %>%
  summary(infer = c(TRUE, TRUE)) %>%
  as_tibble() %>%
  mutate(model = "raw")

log_slopes <- emtrends(
  combined_log,
  ~ Season,
  var = "root_biomass"
) %>%
  summary(infer = c(TRUE, TRUE)) %>%
  as_tibble() %>%
  mutate(model = "log1p")

combined_slopes <- bind_rows(
  raw_slopes,
  log_slopes
)

write_csv(
  combined_slopes,
  file.path(
    result_dir,
    "FigureS1_combined_model_slopes.csv"
  )
)

# Interaction terms
raw_interaction <- summary(combined_raw)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  filter(str_detect(term, "root_biomass:Season"))

log_interaction <- summary(combined_log)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  filter(str_detect(term, "root_biomass:Season"))

interaction_results <- bind_rows(
  raw_interaction %>% mutate(model = "raw"),
  log_interaction %>% mutate(model = "log1p")
)

write_csv(
  interaction_results,
  file.path(
    result_dir,
    "FigureS1_interaction_results.csv"
  )
)

# -----------------------------------------------------------------------------
# 5. Residual diagnostics
# -----------------------------------------------------------------------------

active_raw_diag <- tibble(
  model = "Active raw",
  fitted = fitted(active_raw),
  residual = residuals(active_raw),
  standardized_residual = residuals(active_raw) / sigma(active_raw)
)

active_log_diag <- tibble(
  model = "Active log1p",
  fitted = fitted(active_log),
  residual = residuals(active_log),
  standardized_residual = residuals(active_log) / sigma(active_log)
)

dormant_raw_diag <- tibble(
  model = "Dormant raw",
  fitted = fitted(dormant_raw),
  residual = residuals(dormant_raw),
  standardized_residual = residuals(dormant_raw) / sigma(dormant_raw)
)

dormant_log_diag <- tibble(
  model = "Dormant log1p",
  fitted = fitted(dormant_log),
  residual = residuals(dormant_log),
  standardized_residual = residuals(dormant_log) / sigma(dormant_log)
)

diagnostic_data <- bind_rows(
  active_raw_diag,
  active_log_diag,
  dormant_raw_diag,
  dormant_log_diag
)

# Residual vs fitted
p_resid <- ggplot(
  diagnostic_data,
  aes(fitted, standardized_residual)
) +
  geom_hline(yintercept = 0) +
  geom_hline(
    yintercept = c(-3, 3),
    linetype = "dashed"
  ) +
  geom_point(alpha = 0.7) +
  facet_wrap(~ model, scales = "free") +
  labs(
    x = "Fitted value",
    y = "Standardized residual",
    title = "Residuals versus fitted values"
  ) +
  theme_manuscript

# QQ plots
p_qq <- ggplot(
  diagnostic_data,
  aes(sample = standardized_residual)
) +
  stat_qq() +
  stat_qq_line() +
  facet_wrap(~ model, scales = "free") +
  labs(
    x = "Theoretical quantile",
    y = "Standardized residual",
    title = "Normal Q-Q plots"
  ) +
  theme_manuscript

diagnostic_plot <- p_resid / p_qq

ggsave(
  file.path(
    figure_dir,
    "FigureS1_season_model_diagnostics.pdf"
  ),
  diagnostic_plot,
  width = 10,
  height = 9
)

ggsave(
  file.path(
    figure_dir,
    "FigureS1_season_model_diagnostics.png"
  ),
  diagnostic_plot,
  width = 10,
  height = 9,
  dpi = 320
)

# -----------------------------------------------------------------------------
# 6. Simple numerical diagnostic summary
# -----------------------------------------------------------------------------

diagnostic_summary <- diagnostic_data %>%
  group_by(model) %>%
  summarise(
    max_abs_standardized_residual =
      max(abs(standardized_residual)),
    residual_fitted_correlation =
      cor(abs(residual), fitted),
    shapiro_p =
      shapiro.test(residual)$p.value,
    .groups = "drop"
  )

write_csv(
  diagnostic_summary,
  file.path(
    result_dir,
    "FigureS1_diagnostic_summary.csv"
  )
)

print(diagnostic_summary)

cat("\nSI sensitivity analysis complete.\n")