# =============================================================================
# FIGURE 1: Influence assessment for the three extreme root-biomass observations
# Uses exactly the primary Figure 1 mixed model and compares all data to
# removal of SM1621, SM1078, SM1080. Also performs single-observation deletion
# influence diagnostics for EVERY observation using the same mixed model.
# =============================================================================

packages <- c("tidyverse", "lme4", "lmerTest", "emmeans")
missing <- setdiff(packages, rownames(installed.packages()))
if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

output_dir <- file.path("Results", "manuscript_final", "Figure1", "three_extremes_influence")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Prefer the original raw analysis source. If it is unavailable, use the
# previously exported robust all-observations CSV (which contains 77 model rows).
raw_path <- file.path("Data", "Data_merged_cleaned.csv")
export_path <- file.path("Results", "manuscript_final", "Figure1",
                         "Figure1_root_robust_all_observations.csv")
if (file.exists(raw_path)) {
  dat <- readr::read_csv(raw_path, show_col_types = FALSE) %>%
    filter(as.character(Irrigation_trmt_pct) == "100") %>%
    transmute(
      Sample_cat = as.character(Sample_cat),
      Plot_ID = factor(Plot),
      Sampling_Date = as.Date(Collection_date),
      root_biomass = as.numeric(Roots_gdw),
      total_soil_respiration = as.numeric(Resp_rate_mg_C_m2_h)
    ) %>%
    filter(!is.na(root_biomass), !is.na(total_soil_respiration),
           !is.na(Sampling_Date), !is.na(Plot_ID))
} else if (file.exists(export_path)) {
  message("Original merged file not found; using exported 77-row analytic dataset.")
  dat <- readr::read_csv(export_path, show_col_types = FALSE) %>%
    transmute(
      Sample_cat = as.character(Sample_cat),
      Plot_ID = factor(Plot_ID),
      Sampling_Date = as.Date(Sampling_Date),
      root_biomass = as.numeric(root_biomass),
      total_soil_respiration = as.numeric(total_soil_respiration)
    ) %>%
    filter(!is.na(root_biomass), !is.na(total_soil_respiration))
} else {
  stop("Cannot find Data/Data_merged_cleaned.csv or ", export_path)
}

dat <- dat %>%
  mutate(
    Sampling_Date_factor = factor(Sampling_Date),
    Season = factor(if_else(lubridate::month(Sampling_Date) %in% 4:10,
                            "Plant-active, Apr-Oct", "Plant-dormant, Nov-Mar"),
                    levels = c("Plant-active, Apr-Oct", "Plant-dormant, Nov-Mar"))
  ) %>%
  arrange(Sampling_Date, Plot_ID, Sample_cat) %>%
  mutate(observation_id = paste(Sample_cat, paste0("plot", Plot_ID),
                                Sampling_Date, sep = "_"))

if (anyDuplicated(dat$observation_id)) stop("Nonunique observation_id in data.")
excluded_samples <- c("SM1621", "SM1078", "SM1080")
if (any(table(factor(dat$Sample_cat, levels = excluded_samples)) != 1)) {
  stop("Each target sample must occur exactly once. Check IDs and source file.")
}

# First fit the same crossed model as in the original script. If its plot
# variance is effectively zero, use date random intercept only, matching the
# original primary-model simplification. Hold this structure fixed in ALL
# refits, so comparisons change only because observations are removed.
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000))
full_crossed <- lmer(total_soil_respiration ~ root_biomass * Season +
                       (1 | Plot_ID) + (1 | Sampling_Date_factor),
                     data = dat, REML = TRUE, control = ctrl)
vv <- as.data.frame(VarCorr(full_crossed))
plot_v <- vv$vcov[vv$grp == "Plot_ID"]
date_v <- vv$vcov[vv$grp == "Sampling_Date_factor"]
if (isSingular(full_crossed, tol = 1e-5) && plot_v < 1e-6 && date_v >= 1e-6) {
  formula_used <- total_soil_respiration ~ root_biomass * Season +
    (1 | Sampling_Date_factor)
} else {
  formula_used <- total_soil_respiration ~ root_biomass * Season +
    (1 | Plot_ID) + (1 | Sampling_Date_factor)
}
message("Fixed model for comparisons: ", paste(deparse(formula_used), collapse = " "))
fit_model <- function(d) lmer(formula_used, data = d, REML = TRUE, control = ctrl)
full_fit <- fit_model(dat)
if (isSingular(full_fit, tol = 1e-5)) {
  warning("Primary model is singular. Inspect model structure before interpreting.")
}

get_stats <- function(fit, scenario, n) {
  ss <- summary(emmeans::emtrends(fit, specs = ~ Season,
                                  var = "root_biomass"), infer = c(TRUE, TRUE)) %>%
    as_tibble()
  coef_tab <- as.data.frame(summary(fit)$coefficients) %>%
    rownames_to_column("term") %>% as_tibble()
  interaction <- coef_tab %>% filter(grepl("^root_biomass:Season", term))
  tibble(scenario = scenario, n = n,
         singular = isSingular(fit, tol = 1e-5),
         season = as.character(ss$Season),
         slope = ss$root_biomass.trend,
         lower_95 = ss$lower.CL, upper_95 = ss$upper.CL,
         slope_p = ss$p.value,
         interaction_estimate = interaction$Estimate,
         interaction_p = interaction$`Pr(>|t|)`)
}

reduced_dat <- dat %>% filter(!Sample_cat %in% excluded_samples)
reduced_fit <- fit_model(reduced_dat)
comparison <- bind_rows(
  get_stats(full_fit, "All observations", nrow(dat)),
  get_stats(reduced_fit, "Exclude SM1621, SM1078, SM1080", nrow(reduced_dat))
)
write_csv(comparison, file.path(output_dir, "three_point_model_comparison.csv"))

# Single-deletion model-based influence. A Cook-like fixed-effects distance:
# (beta_minus_i - beta_full)' V_full^-1 (beta_minus_i - beta_full) / p
# where V_full is covariance of full-model fixed effects, p = number of fixed
# effects. This is NOT ordinary Cook's distance from lm, and does not have
# a validated p-value or universal 4/n threshold for mixed models.
beta_full <- fixef(full_fit)
V_full <- as.matrix(vcov(full_fit))
p <- length(beta_full)
refit_one <- function(i) {
  tryCatch({
    fit <- fit_model(dat[-i, , drop = FALSE])
    delta <- fixef(fit)[names(beta_full)] - beta_full
    cook_like <- as.numeric(t(delta) %*% solve(V_full, delta)) / p
    slopes <- get_stats(fit, "Exclude one", nrow(dat) - 1)
    act <- slopes %>% filter(season == "Plant-active, Apr-Oct")
    tibble(observation_id = dat$observation_id[i],
           Sample_cat = dat$Sample_cat[i],
           Plot_ID = as.character(dat$Plot_ID[i]),
           Sampling_Date = dat$Sampling_Date[i],
           root_biomass = dat$root_biomass[i],
           total_soil_respiration = dat$total_soil_respiration[i],
           is_target = dat$Sample_cat[i] %in% excluded_samples,
           cook_like_fixed_effect_distance = cook_like,
           active_slope_if_excluded = act$slope,
           change_in_active_slope = act$slope -
             (get_stats(full_fit, "Full", nrow(dat)) %>%
                filter(season == "Plant-active, Apr-Oct"))$slope,
           active_slope_p_if_excluded = act$slope_p,
           interaction_p_if_excluded = act$interaction_p,
           refit_singular = isSingular(fit, tol = 1e-5),
           error = NA_character_)
  }, error = function(e) {
    tibble(observation_id = dat$observation_id[i],
           Sample_cat = dat$Sample_cat[i],
           Plot_ID = as.character(dat$Plot_ID[i]),
           Sampling_Date = dat$Sampling_Date[i],
           root_biomass = dat$root_biomass[i],
           total_soil_respiration = dat$total_soil_respiration[i],
           is_target = dat$Sample_cat[i] %in% excluded_samples,
           cook_like_fixed_effect_distance = NA_real_,
           active_slope_if_excluded = NA_real_,
           change_in_active_slope = NA_real_,
           active_slope_p_if_excluded = NA_real_,
           interaction_p_if_excluded = NA_real_,
           refit_singular = NA,
           error = conditionMessage(e))
  })
}
message("Computing leave-one-observation-out influence for ", nrow(dat), " observations...")
influence <- bind_rows(lapply(seq_len(nrow(dat)), refit_one)) %>%
  arrange(desc(cook_like_fixed_effect_distance)) %>%
  mutate(influence_rank = row_number())
write_csv(influence, file.path(output_dir, "mixed_model_leave_one_out_influence.csv"))

p_influence <- ggplot(influence,
                      aes(x = reorder(observation_id, cook_like_fixed_effect_distance),
                          y = cook_like_fixed_effect_distance,
                          fill = is_target)) +
  geom_col() + coord_flip() +
  scale_fill_manual(values = c(`FALSE` = "grey65", `TRUE` = "#C8372D"),
                    labels = c(`FALSE` = "Other observations", `TRUE` = "Three extreme root values")) +
  labs(x = NULL, y = "Fixed-effect influence (Cook-like distance)", fill = NULL,
       title = "Mixed-model influence: remove one observation at a time") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 5))
ggsave(file.path(output_dir, "mixed_model_influence.pdf"), p_influence,
       width = 8, height = 12)

cat("\n===== MODEL COMPARISON =====\n")
print(comparison, n = Inf)
cat("\n===== THREE EXTREME OBSERVATIONS: INFLUENCE RANKS =====\n")
print(influence %>% filter(is_target) %>%
        select(Sample_cat, root_biomass, influence_rank,
               cook_like_fixed_effect_distance,
               change_in_active_slope, active_slope_p_if_excluded), n = Inf)
cat("\nOutputs saved to: ", output_dir, "\n", sep = "")
