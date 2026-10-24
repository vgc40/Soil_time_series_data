# =============================================================================
# Figure 2: seasonal substrate-induced respiration (Q1; H1a, H1b)
# =============================================================================
#
# Q1: How does substrate-induced respiration (SIR) differ between plant-active
#     and plant-dormant periods, and do seasonal responses differ among
#     substrates?
# H1a: SIR profiles differ between periods.
# H1b: Relative to each substrate's mean across sampled dates, readily
#      available substrates respond more during the plant-active period and
#      complex residue-associated substrates respond more during the
#      plant-dormant period.
#
# Data: MicroResp v3, 100% irrigation, 5-h, carbon-normalized SIR
# (ngCO2-C mgC-1 h-1), 31 substrates, 77 samples on 13 dates.
#
# Main figure
#   A. NMDS of the 31-substrate SIR profile; date-level PERMANOVA (H1a).
#   B. Per-substrate seasonal difference (dormant - active, log2 fold change).
#   C. Seasonal difference by availability class; class x season test (H1b).
#   D. Class trajectories relative to each substrate's annual mean.
# Supplement
#   S2. GAM trajectories of all 31 substrates, grouped by substrate group.
#   S2b. Sensitivity of the H1b contrast to analysis choices.
#   Table: environmental vector fitting on date-mean ordination scores.
#
# Response scale: log2(SIR) centred on the substrate's mean of date means,
# so 0 = the substrate's annual mean and +1 = twice that mean. Seasonal
# differences on this scale are log2 fold changes.

rm(list = ls())
source(file.path("Scripts", "helpers.R"))

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(mgcv)
  library(patchwork)
  library(ggrepel)
})

set.seed(20231023)
emm_options(lmer.df = "satterthwaite")
dirs <- output_dirs("Figure2")

classes <- read_csv(input_files$substrate_classes, show_col_types = FALSE)
availability_levels <- c("Readily available", "Complex", "Other")

# -----------------------------------------------------------------------------
# Data preparation
# -----------------------------------------------------------------------------

prepare_sir <- function(sir, availability_col = "availability_class",
                        standardize = c("centred_log2", "z")) {
  standardize <- match.arg(standardize)
  sir %>%
    left_join(classes, by = "Substrate") %>%
    mutate(
      availability = factor(.data[[availability_col]], levels = availability_levels),
      log2_sir = log2(SIR)
    ) %>%
    group_by(Substrate) %>%
    mutate(
      annual_mean = mean(tapply(log2_sir, Sampling_Date, mean)),
      response = if (standardize == "centred_log2") {
        log2_sir - annual_mean
      } else {
        (log2_sir - mean(log2_sir)) / sd(log2_sir)
      }
    ) %>%
    ungroup() %>%
    mutate(Date_f = factor(Sampling_Date))
}

sir_raw <- load_microresp_100()
stopifnot(n_distinct(sir_raw$Sample_cat) == 77, n_distinct(sir_raw$Substrate) == 31)
sir <- prepare_sir(sir_raw)

sir_wide <- sir_raw %>%
  select(Sample_cat, Sampling_Date, Plot_ID, Season, Substrate, SIR) %>%
  pivot_wider(names_from = Substrate, values_from = SIR) %>%
  arrange(Sampling_Date, Plot_ID)
sir_matrix <- as.matrix(select(sir_wide, -c(Sample_cat, Sampling_Date, Plot_ID, Season)))
stopifnot(!anyNA(sir_matrix))

# -----------------------------------------------------------------------------
# A. Whole-profile ordination and date-level season test (H1a)
# -----------------------------------------------------------------------------

nmds <- metaMDS(sir_matrix, distance = "bray", k = 2, trymax = 200,
                autotransform = FALSE, trace = 0)

nmds_scores <- as_tibble(scores(nmds, display = "sites")) %>%
  bind_cols(select(sir_wide, Sample_cat, Sampling_Date, Plot_ID, Season))

date_centroids <- nmds_scores %>%
  group_by(Sampling_Date, Season) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = "drop") %>%
  arrange(Sampling_Date) %>%
  mutate(label = format(Sampling_Date, "%b %y"))

h1a_test <- date_level_season_test(sir_matrix, sir_wide$Sampling_Date)
write_csv(h1a_test, file.path(dirs$results, "Figure2A_season_PERMANOVA_dispersion.csv"))
write_csv(nmds_scores, file.path(dirs$results, "Figure2A_NMDS_sample_scores.csv"))

primary_h1a <- filter(h1a_test, level == "Date means (primary)")

panel_a <- ggplot(nmds_scores, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour = Season), alpha = 0.35, size = 1.4) +
  geom_path(data = date_centroids, colour = "grey45", linewidth = 0.4,
            arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
  geom_point(data = date_centroids, aes(fill = Season), shape = 21,
             size = 3.2, colour = "black", stroke = 0.4) +
  geom_text_repel(data = date_centroids, aes(label = label), size = 2.5,
                  min.segment.length = 0, seed = 1, max.overlaps = 20) +
  scale_colour_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  labs(
    title = "SIR profile",
    subtitle = sprintf(
      "Date-level PERMANOVA R² = %.2f, p = %s; dispersion p = %s; stress = %.2f",
      primary_h1a$R2, format_p(primary_h1a$p_permanova),
      format_p(primary_h1a$p_dispersion), nmds$stress
    ),
    colour = NULL, fill = NULL
  ) +
  theme_manuscript +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 7.5))

# Environmental vector fitting (supplement): date-mean NMDS scores against
# date-mean field conditions, permuting whole dates.
field <- load_field_100()
env_dates <- field %>%
  group_by(Sampling_Date) %>%
  summarise(
    `Soil temperature` = mean(SoilTemp_mean),
    `Gravimetric water content` = mean(GWC),
    `Root biomass` = mean(Roots_gdw),
    `Field soil respiration` = mean(Resp_rate_mg_C_m2_h),
    `Microbial biomass C` = mean(MBC_ug_gdw_bgc),
    .groups = "drop"
  )
envfit_input <- date_centroids %>% left_join(env_dates, by = "Sampling_Date")
envfit_result <- envfit(
  select(envfit_input, NMDS1, NMDS2),
  select(envfit_input, -c(Sampling_Date, Season, NMDS1, NMDS2, label)),
  permutations = 9999
)
envfit_table <- tibble(
  variable = names(envfit_result$vectors$r),
  r2 = envfit_result$vectors$r,
  p = envfit_result$vectors$pvals,
  NMDS1_direction = envfit_result$vectors$arrows[, 1],
  NMDS2_direction = envfit_result$vectors$arrows[, 2],
  n_dates = nrow(envfit_input),
  note = "Exploratory; date-mean scores and conditions; 9999 free permutations of dates."
)
write_csv(envfit_table, file.path(dirs$results, "Figure2_envfit_date_level.csv"))

# -----------------------------------------------------------------------------
# B. Per-substrate seasonal differences
# -----------------------------------------------------------------------------

fit_season_effect <- function(d) {
  m <- lmer(response ~ Season + (1 | Date_f) + (1 | Plot_ID), data = d)
  co <- summary(m)$coefficients["SeasonPlant-dormant", ]
  tibble(
    estimate = co[["Estimate"]], se = co[["Std. Error"]], df = co[["df"]],
    p = co[["Pr(>|t|)"]],
    lower = estimate - qt(0.975, df) * se,
    upper = estimate + qt(0.975, df) * se,
    singular = isSingular(m)
  )
}

substrate_effects <- sir %>%
  group_by(Substrate, display_name, availability, Table_group, compound_class) %>%
  group_modify(~ fit_season_effect(.x)) %>%
  ungroup() %>%
  mutate(
    p_fdr = p.adjust(p, method = "BH"),
    fold_change = 2^estimate
  )

# GAM trajectories (supplement) identify substrates without detectable
# temporal structure ("steady" substrates, KH comment 47).
fit_gam <- function(d) {
  d <- mutate(d, time_days = as.numeric(Sampling_Date - min(sir$Sampling_Date)))
  m <- gam(log2_sir ~ s(time_days, k = 7) + s(Plot_ID, bs = "re"),
           data = d, method = "REML")
  s_tab <- summary(m)$s.table
  grid <- tibble(
    Sampling_Date = seq(min(d$Sampling_Date), max(d$Sampling_Date), by = "3 days"),
    Plot_ID = d$Plot_ID[1]
  ) %>%
    mutate(time_days = as.numeric(Sampling_Date - min(sir$Sampling_Date)))
  pred <- predict(m, grid, exclude = "s(Plot_ID)", se.fit = TRUE)
  list(
    stats = tibble(
      edf = s_tab["s(time_days)", "edf"],
      F = s_tab["s(time_days)", "F"],
      p = s_tab["s(time_days)", "p-value"],
      deviance_explained = summary(m)$dev.expl
    ),
    fit = grid %>% mutate(fit = pred$fit, se = pred$se.fit)
  )
}

gam_fits <- sir %>%
  group_by(Substrate) %>%
  group_map(~ c(list(Substrate = .y$Substrate), fit_gam(.x)))

gam_stats <- map_dfr(gam_fits, ~ mutate(.x$stats, Substrate = .x$Substrate)) %>%
  mutate(p_fdr = p.adjust(p, method = "BH"))
gam_curves <- map_dfr(gam_fits, ~ mutate(.x$fit, Substrate = .x$Substrate))

substrate_effects <- substrate_effects %>%
  left_join(select(gam_stats, Substrate, gam_deviance_explained = deviance_explained,
                   gam_p_fdr = p_fdr), by = "Substrate") %>%
  mutate(steady = gam_p_fdr > 0.05)

write_csv(substrate_effects, file.path(dirs$results, "Figure2B_substrate_season_effects.csv"))
write_csv(gam_stats %>% left_join(classes, by = "Substrate"),
          file.path(dirs$results, "FigureS2_GAM_temporal_statistics.csv"))

plot_b_data <- substrate_effects %>%
  arrange(availability, estimate) %>%
  mutate(
    label = if_else(steady, paste0(display_name, " (steady)"), display_name),
    label = factor(label, levels = label),
    signif = if_else(p_fdr < 0.05, "FDR < 0.05", "FDR ≥ 0.05")
  )

panel_b <- ggplot(plot_b_data, aes(estimate, label, colour = availability)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, linewidth = 0.5, orientation = "y") +
  geom_point(aes(shape = signif), size = 2) +
  scale_colour_manual(values = class_colors[availability_levels], drop = FALSE) +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "FDR ≥ 0.05" = 1)) +
  facet_grid(availability ~ ., scales = "free_y", space = "free_y") +
  labs(
    title = "Seasonal difference by substrate",
    x = "Dormant − active (log2 fold change, 95% CI)", y = NULL,
    colour = NULL, shape = NULL
  ) +
  theme_manuscript +
  theme(legend.position = "bottom", strip.text.y = element_text(size = 7),
        axis.text.y = element_text(size = 7)) +
  guides(colour = "none")

# -----------------------------------------------------------------------------
# C. Class x season test (H1b)
# -----------------------------------------------------------------------------
# Mixed model on all 31 substrates. Random effects: date (season is a date-
# level treatment), date x class (the error term for class-specific seasonal
# change), sample (substrates measured on the same soil), and substrate x
# season (substrates within a class differ in their seasonal response).

fit_class_model <- function(d) {
  m <- lmer(
    response ~ Season * availability +
      (1 | Date_f) + (1 | Date_f:availability) + (1 | Sample_cat) +
      (1 | Substrate:Season),
    data = d,
    control = lmerControl(optimizer = "bobyqa")
  )
  season_contrasts <- emmeans(m, ~ Season | availability) %>%
    contrast(method = "revpairwise")
  season_by_class <- season_contrasts %>%
    summary(infer = TRUE) %>%
    as_tibble() %>%
    transmute(availability, contrast = "Dormant - active", estimate, SE, df,
              lower = lower.CL, upper = upper.CL, p = p.value)
  # Difference of the seasonal contrasts, complex minus readily available
  # (availability levels: Readily available, Complex, Other).
  interaction <- contrast(
    season_contrasts, by = NULL,
    method = list("(Dormant - active): Complex - Readily available" =
                    c(-1, 1, 0))
  ) %>%
    summary(infer = TRUE) %>%
    as_tibble() %>%
    transmute(contrast = as.character(contrast), estimate, SE, df,
              lower = lower.CL, upper = upper.CL, p = p.value)
  list(model = m, season_by_class = season_by_class, interaction = interaction,
       singular = isSingular(m))
}

# Distribution-free check: compare substrate-level seasonal effects between
# the two H1b groups by exact enumeration of group labels.
class_permutation_test <- function(effects) {
  e <- filter(effects, availability %in% c("Readily available", "Complex"))
  obs <- mean(e$estimate[e$availability == "Complex"]) -
    mean(e$estimate[e$availability == "Readily available"])
  n_complex <- sum(e$availability == "Complex")
  perm <- combn(nrow(e), n_complex, function(idx) {
    mean(e$estimate[idx]) - mean(e$estimate[-idx])
  })
  tibble(
    observed_difference = obs,
    p_two_sided = mean(abs(perm) >= abs(obs) - 1e-12),
    p_one_sided_H1b = mean(perm >= obs - 1e-12),
    n_label_sets = length(perm),
    n_complex = n_complex,
    n_readily_available = nrow(e) - n_complex
  )
}

primary_class <- fit_class_model(sir)
class_perm <- class_permutation_test(substrate_effects)

# Shared seasonal shift across all 31 substrates: each sample's mean
# response over substrates, with date as the replication unit for season.
overall_shift <- sir %>%
  group_by(Sample_cat, Date_f, Plot_ID, Season) %>%
  summarise(response = mean(response), .groups = "drop") %>%
  fit_season_effect() %>%
  mutate(contrast = "All 31 substrates, dormant - active")
write_csv(overall_shift, file.path(dirs$results, "Figure2C_overall_seasonal_shift.csv"))

write_csv(primary_class$season_by_class, file.path(dirs$results, "Figure2C_season_effect_by_class.csv"))
write_csv(bind_rows(primary_class$interaction %>% mutate(test = "Mixed model"),
                    class_perm %>% transmute(test = "Exact label permutation (substrate effects)",
                                             contrast = "(Dormant - active): Complex - Readily available",
                                             estimate = observed_difference,
                                             p = p_two_sided, p_one_sided_H1b)),
          file.path(dirs$results, "Figure2C_H1b_class_by_season_tests.csv"))

h1b <- primary_class$interaction

panel_c <- ggplot(primary_class$season_by_class,
                  aes(availability, estimate, colour = availability)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_jitter(data = substrate_effects, aes(availability, estimate),
              width = 0.12, alpha = 0.4, size = 1.2, inherit.aes = FALSE,
              colour = "grey40") +
  geom_pointrange(aes(ymin = lower, ymax = upper), size = 0.5, linewidth = 0.8) +
  scale_colour_manual(values = class_colors[availability_levels]) +
  labs(
    title = "Seasonal difference by class",
    subtitle = sprintf(
      "All substrates: %.2f (95%% CI %.2f to %.2f), p = %s\nComplex − readily available = %.2f (%.2f to %.2f), p = %s; exact permutation p = %s",
      overall_shift$estimate, overall_shift$lower, overall_shift$upper,
      format_p(overall_shift$p),
      h1b$estimate, h1b$lower, h1b$upper, format_p(h1b$p),
      format_p(class_perm$p_two_sided)
    ),
    x = NULL, y = "Dormant − active (log2 fold change)"
  ) +
  theme_manuscript +
  theme(legend.position = "none", plot.subtitle = element_text(size = 7.5))

# -----------------------------------------------------------------------------
# D. Class trajectories
# -----------------------------------------------------------------------------

class_trajectories <- sir %>%
  group_by(Sample_cat, Sampling_Date, Plot_ID, Season, availability) %>%
  summarise(response = mean(response), .groups = "drop") %>%
  group_by(Sampling_Date, Season, availability) %>%
  summarise(
    mean = mean(response), n = n(), se = sd(response) / sqrt(n),
    lower = mean - qt(0.975, n - 1) * se, upper = mean + qt(0.975, n - 1) * se,
    .groups = "drop"
  )
write_csv(class_trajectories, file.path(dirs$results, "Figure2D_class_trajectories.csv"))

panel_d <- ggplot(class_trajectories, aes(Sampling_Date, mean, colour = availability)) +
  season_band_layers(sir$Sampling_Date) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = availability),
              alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  scale_fill_manual(values = class_colors[availability_levels], guide = "none") +
  scale_colour_manual(values = class_colors[availability_levels]) +
  date_axis +
  labs(
    title = "Class trajectories",
    x = NULL, y = "log2 SIR relative to annual mean", colour = NULL
  ) +
  theme_manuscript +
  theme(legend.position = "bottom")

figure2 <- (panel_a + panel_d) / (panel_b + panel_c) +
  plot_layout(heights = c(1, 1.5), widths = c(1, 1)) +
  plot_annotation(tag_levels = "A")
save_figure(figure2, dirs$figures, "Figure2_main", width = 11, height = 11)

# -----------------------------------------------------------------------------
# Supplement S2: all 31 GAM trajectories grouped by substrate group
# -----------------------------------------------------------------------------

facet_order <- gam_stats %>%
  left_join(classes, by = "Substrate") %>%
  arrange(Table_group, desc(deviance_explained)) %>%
  mutate(
    facet = sprintf("%s\n%s · %.0f%%, FDR p %s", Table_group, display_name,
                    100 * deviance_explained, format_p(p_fdr))
  )

s2_points <- sir %>%
  left_join(select(facet_order, Substrate, facet), by = "Substrate") %>%
  mutate(facet = factor(facet, levels = facet_order$facet))
s2_curves <- gam_curves %>%
  left_join(select(facet_order, Substrate, facet), by = "Substrate") %>%
  mutate(facet = factor(facet, levels = facet_order$facet))

figure_s2 <- ggplot(s2_points, aes(Sampling_Date, log2_sir)) +
  season_band_layers(sir$Sampling_Date) +
  geom_point(size = 0.6, alpha = 0.5) +
  geom_ribbon(data = s2_curves, aes(Sampling_Date, ymin = fit - 1.96 * se,
                                    ymax = fit + 1.96 * se),
              inherit.aes = FALSE, fill = "grey30", alpha = 0.2) +
  geom_line(data = s2_curves, aes(Sampling_Date, fit), inherit.aes = FALSE,
            linewidth = 0.6) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b\n%y") +
  facet_wrap(~ facet, scales = "free_y", ncol = 6) +
  labs(x = NULL, y = expression(log[2]~SIR~(ngCO[2]*"-C"~mgC^{-1}~h^{-1}))) +
  theme_manuscript +
  theme(strip.text = element_text(size = 5.5), axis.text = element_text(size = 6))
save_figure(figure_s2, dirs$figures, "FigureS2_substrate_GAM_trajectories",
            width = 12, height = 12)

# -----------------------------------------------------------------------------
# Supplement S2b: sensitivity of the H1b contrast
# -----------------------------------------------------------------------------

run_sensitivity <- function(label, d) {
  fit <- fit_class_model(d)
  fit$interaction %>% mutate(analysis = label, singular = fit$singular,
                             n_obs = nrow(d))
}

sir_24h <- prepare_sir(load_microresp_100(timepoint = 24))
sir_gdw <- prepare_sir(load_microresp_100(normalization = "ngCO2-C_gdw_h"))

sensitivity <- bind_rows(
  run_sensitivity("Primary (5 h, per mg C, centred log2)", sir),
  run_sensitivity("Amino acids counted as readily available",
                  prepare_sir(sir_raw, "availability_class_sensitivity")),
  run_sensitivity("Within-substrate z-score", prepare_sir(sir_raw, standardize = "z")),
  run_sensitivity("Transition dates (Oct 2023, Apr 2024) excluded",
                  prepare_sir(filter(sir_raw, !Sampling_Date %in% transition_dates))),
  run_sensitivity("24-h incubation", sir_24h),
  run_sensitivity("Per g dry soil", sir_gdw)
)

# The H1a test under the same alternatives.
h1a_sensitivity <- list(
  "Transition dates excluded" = filter(sir_raw, !Sampling_Date %in% transition_dates),
  "24-h incubation" = load_microresp_100(timepoint = 24),
  "Per g dry soil" = load_microresp_100(normalization = "ngCO2-C_gdw_h")
) %>%
  imap_dfr(function(d, label) {
    w <- d %>% select(Sample_cat, Sampling_Date, Plot_ID, Substrate, SIR) %>%
      pivot_wider(names_from = Substrate, values_from = SIR) %>%
      arrange(Sampling_Date, Plot_ID)
    date_level_season_test(as.matrix(select(w, -c(Sample_cat, Sampling_Date, Plot_ID))),
                           w$Sampling_Date) %>%
      mutate(analysis = label)
  })

write_csv(sensitivity, file.path(dirs$results, "FigureS2b_H1b_sensitivity.csv"))
write_csv(h1a_sensitivity, file.path(dirs$results, "FigureS2b_H1a_sensitivity.csv"))

figure_s2b <- sensitivity %>%
  mutate(analysis = fct_rev(fct_inorder(analysis))) %>%
  ggplot(aes(estimate, analysis)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, orientation = "y") +
  geom_point(size = 2) +
  labs(
    x = "Complex − readily available difference in seasonal effect\n(dormant − active; positive supports H1b)",
    y = NULL
  ) +
  theme_manuscript
save_figure(figure_s2b, dirs$figures, "FigureS2b_H1b_sensitivity", width = 7, height = 3.5)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

summary_table <- tribble(
  ~item, ~value,
  "Samples", as.character(n_distinct(sir$Sample_cat)),
  "Sampling dates (active / dormant)", sprintf("%d (%d / %d)", primary_h1a$n_units,
                                               primary_h1a$n_active_dates, primary_h1a$n_dormant_dates),
  "Substrates", as.character(n_distinct(sir$Substrate)),
  "NMDS stress", sprintf("%.3f", nmds$stress),
  "H1a date-level PERMANOVA R2", sprintf("%.3f", primary_h1a$R2),
  "H1a date-level PERMANOVA p (exact, 1287 label sets)", format_p(primary_h1a$p_permanova),
  "H1a date-level dispersion p", format_p(primary_h1a$p_dispersion),
  "Substrates higher in dormancy (FDR < 0.05)",
  as.character(sum(substrate_effects$estimate > 0 & substrate_effects$p_fdr < 0.05)),
  "Substrates higher in active period (FDR < 0.05)",
  as.character(sum(substrate_effects$estimate < 0 & substrate_effects$p_fdr < 0.05)),
  "Substrates with dormant > active point estimate",
  sprintf("%d of %d", sum(substrate_effects$estimate > 0), nrow(substrate_effects)),
  "Overall shift, all substrates (dormant - active, log2)",
  sprintf("%.3f [%.3f, %.3f], p = %s", overall_shift$estimate, overall_shift$lower,
          overall_shift$upper, format_p(overall_shift$p)),
  "Steady substrates (GAM FDR > 0.05)",
  paste(substrate_effects$display_name[substrate_effects$steady], collapse = "; "),
  "H1b interaction (complex - readily available), log2",
  sprintf("%.3f [%.3f, %.3f], p = %s", h1b$estimate, h1b$lower, h1b$upper, format_p(h1b$p)),
  "H1b exact permutation p (two-sided)", format_p(class_perm$p_two_sided),
  "Class model singular", as.character(primary_class$singular)
)
write_csv(summary_table, file.path(dirs$results, "Figure2_analysis_summary.csv"))
print(summary_table, n = Inf)
print(primary_class$season_by_class)
print(sensitivity)
