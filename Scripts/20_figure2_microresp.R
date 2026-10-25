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
# A season difference in an SIR profile can be of two kinds, tested
# separately:
#   level - all substrates respond more (or less) in one period;
#   mix   - the relative pattern across substrates changes (some substrates
#           go up relative to others).
#
# Main figure
#   A. Level: mean SIR across all 31 substrates through time (H1a, level).
#   B. Whole profile: NMDS; date-level tests of the whole profile and of the
#      mix alone (H1a).
#   C. Seasonal difference for each substrate, grouped by chemistry.
#   D. H1b: complex vs readily available overall, and polymer vs its
#      building blocks within each chemical class.
# Supplement
#   S2.  GAM trajectories of all 31 substrates, grouped by chemistry.
#   S2b. Sensitivity of the H1b contrast to analysis choices.
#   S2c. Trajectories of complex and readily available substrates.
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
availability_levels <- c("Readily available", "Complex")
group_levels <- c("Carbohydrate", "Organic acid", "Amino acid", "Amino sugar",
                  "Aromatic", "Lipid-like", "Aliphatic amine")

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
      substrate_group = factor(substrate_group, levels = group_levels),
      log2_sir = log2(SIR)
    ) %>%
    filter(!is.na(availability)) %>%
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
stopifnot(!anyNA(sir$substrate_group))

to_wide <- function(d) {
  w <- d %>%
    select(Sample_cat, Sampling_Date, Plot_ID, Substrate, SIR) %>%
    pivot_wider(names_from = Substrate, values_from = SIR) %>%
    arrange(Sampling_Date, Plot_ID)
  list(meta = select(w, Sample_cat, Sampling_Date, Plot_ID),
       matrix = as.matrix(select(w, -c(Sample_cat, Sampling_Date, Plot_ID))))
}
sir_wide <- to_wide(sir_raw)
stopifnot(!anyNA(sir_wide$matrix))

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

bands_layers <- function() season_band_layers(sir$Sampling_Date)

# -----------------------------------------------------------------------------
# A. Level: mean response across all substrates (H1a, level)
# -----------------------------------------------------------------------------
# Each sample's mean over the 31 substrates; the season effect is tested with
# date as the replication unit (mixed model with date and plot intercepts).

sample_level <- sir %>%
  group_by(Sample_cat, Sampling_Date, Date_f, Plot_ID, Season) %>%
  summarise(response = mean(response), .groups = "drop")

overall_shift <- fit_season_effect(sample_level) %>%
  mutate(contrast = "All 31 substrates, dormant - active",
         percent_change = 100 * (2^estimate - 1))
write_csv(overall_shift, file.path(dirs$results, "Figure2A_overall_level_shift.csv"))

# How much of the variation in a sample's overall SIR level is shared by all
# plots on a date (why dates, not plots, are the replicates for season).
variance_model <- lmer(response ~ 1 + (1 | Date_f) + (1 | Plot_ID), data = sample_level)
variance_shares <- as.data.frame(VarCorr(variance_model)) %>%
  transmute(source = recode(grp, Date_f = "Sampling date", Plot_ID = "Plot",
                            Residual = "Plot x date (residual)"),
            variance = vcov, share = vcov / sum(vcov))
write_csv(variance_shares, file.path(dirs$results, "Figure2A_variance_shares_date_vs_plot.csv"))

level_dates <- sample_level %>%
  group_by(Sampling_Date, Season) %>%
  summarise(mean = mean(response), n = n(), se = sd(response) / sqrt(n),
            lower = mean - qt(0.975, n - 1) * se, upper = mean + qt(0.975, n - 1) * se,
            .groups = "drop")
write_csv(level_dates, file.path(dirs$results, "Figure2A_overall_level_by_date.csv"))
write_csv(sample_level, file.path(dirs$results, "Figure2A_overall_level_by_sample.csv"))

panel_a <- ggplot(level_dates, aes(Sampling_Date, mean)) +
  bands_layers() +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_point(data = sample_level, aes(y = response), colour = "grey40",
             alpha = 0.45, size = 1) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey20", alpha = 0.15) +
  geom_line(linewidth = 0.7) +
  geom_point(aes(colour = Season), size = 2) +
  scale_colour_manual(values = season_colors) +
  date_axis +
  labs(
    title = "Level: SIR averaged over all 31 substrates",
    subtitle = sprintf("Dormant vs active: %+.0f%% (log2 %.2f, 95%% CI %.2f to %.2f), %s",
                       overall_shift$percent_change, overall_shift$estimate,
                       overall_shift$lower, overall_shift$upper, p_label(overall_shift$p)),
    x = NULL, y = "log2 SIR relative to each\nsubstrate's annual mean", colour = NULL
  ) +
  theme_manuscript +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 7.5))

# -----------------------------------------------------------------------------
# B. Whole profile and mix (H1a)
# -----------------------------------------------------------------------------
# Whole profile: Bray-Curtis on the 31 SIR rates (sensitive to level and mix).
# Mix only: each sample's rates divided by its own total, so a sample that is
# uniformly higher has the same profile; differences are in relative pattern.

nmds <- metaMDS(sir_wide$matrix, distance = "bray", k = 2, trymax = 200,
                autotransform = FALSE, trace = 0)
nmds_scores <- as_tibble(scores(nmds, display = "sites")) %>%
  bind_cols(sir_wide$meta) %>%
  mutate(Season = assign_season(Sampling_Date))
date_centroids <- nmds_scores %>%
  group_by(Sampling_Date, Season) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = "drop") %>%
  mutate(label = format(Sampling_Date, "%b %y"))

profile_test <- date_level_season_test(sir_wide$matrix, sir_wide$meta$Sampling_Date)
mix_test <- date_level_season_test(sir_wide$matrix / rowSums(sir_wide$matrix),
                                   sir_wide$meta$Sampling_Date)
h1a_tests <- bind_rows(
  profile_test %>% mutate(analysis = "Whole profile (level + mix)"),
  mix_test %>% mutate(analysis = "Mix only (profiles scaled to each sample's total)")
)
write_csv(h1a_tests, file.path(dirs$results, "Figure2B_season_PERMANOVA_profile_and_mix.csv"))
write_csv(nmds_scores, file.path(dirs$results, "Figure2B_NMDS_sample_scores.csv"))

primary_profile <- filter(profile_test, level == "Date means (primary)")
primary_mix <- filter(mix_test, level == "Date means (primary)")

panel_b <- ggplot(nmds_scores, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour = Season), alpha = 0.45, size = 1.5) +
  geom_point(data = date_centroids, aes(fill = Season), shape = 21,
             size = 3.4, colour = "black", stroke = 0.4) +
  geom_text_repel(data = date_centroids, aes(label = label), size = 2.5,
                  min.segment.length = 0.2, seed = 1, max.overlaps = 20) +
  scale_colour_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  labs(
    title = "Whole profile of 31 substrates (NMDS)",
    subtitle = sprintf(
      "Small points: samples; large: date means. Season tests with dates as replicates:\nwhole profile R² = %.2f, %s; mix only R² = %.2f, %s",
      primary_profile$R2, p_label(primary_profile$p_permanova),
      primary_mix$R2, p_label(primary_mix$p_permanova)
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
  n_dates = nrow(envfit_input),
  note = "Exploratory; date-mean scores and conditions; 9999 free permutations of dates."
)
write_csv(envfit_table, file.path(dirs$results, "Figure2_envfit_date_level.csv"))

# -----------------------------------------------------------------------------
# C. Seasonal difference for each substrate
# -----------------------------------------------------------------------------

substrate_effects <- sir %>%
  group_by(Substrate, display_name, availability, substrate_group, compound_class,
           polymer_monomer_class) %>%
  group_modify(~ fit_season_effect(.x)) %>%
  ungroup() %>%
  mutate(p_fdr = p.adjust(p, method = "BH"), fold_change = 2^estimate)

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
      edf = s_tab["s(time_days)", "edf"], F = s_tab["s(time_days)", "F"],
      p = s_tab["s(time_days)", "p-value"], deviance_explained = summary(m)$dev.expl
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
write_csv(substrate_effects, file.path(dirs$results, "Figure2C_substrate_season_effects.csv"))
write_csv(gam_stats %>% left_join(classes, by = "Substrate"),
          file.path(dirs$results, "FigureS2_GAM_temporal_statistics.csv"))

plot_c_data <- substrate_effects %>%
  arrange(substrate_group, availability, estimate) %>%
  mutate(
    label = if_else(steady, paste0(display_name, " (steady)"), display_name),
    label = factor(label, levels = label),
    signif = if_else(p_fdr < 0.05, "Significant after FDR", "Not significant after FDR")
  )

panel_c <- ggplot(plot_c_data, aes(estimate, label, colour = availability)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, linewidth = 0.5,
                orientation = "y") +
  geom_point(aes(shape = signif), size = 2) +
  scale_colour_manual(values = class_colors[availability_levels]) +
  scale_shape_manual(values = c("Significant after FDR" = 16, "Not significant after FDR" = 1)) +
  facet_grid(substrate_group ~ ., scales = "free_y", space = "free_y",
             labeller = label_wrap_gen(width = 9)) +
  labs(
    title = "Seasonal difference for each substrate",
    subtitle = "Right of 0 = higher in dormancy; +1 = twice as high",
    x = "Dormant − active (log2 fold change, 95% CI)", y = NULL,
    colour = NULL, shape = NULL
  ) +
  theme_manuscript +
  theme(legend.position = "bottom", legend.box = "vertical",
        strip.text.y = element_text(size = 6.5, angle = 0),
        axis.text.y = element_text(size = 7), plot.subtitle = element_text(size = 7.5))

# -----------------------------------------------------------------------------
# D. H1b: complex vs readily available
# -----------------------------------------------------------------------------
# Mixed model. Random effects: date (season is a date-level treatment), date
# x availability (the error term for availability-specific seasonal change),
# sample (substrates measured on the same soil), and substrate x season
# (substrates within a group differ in their seasonal response).

fit_availability_model <- function(d, extra_fixed = NULL) {
  fixed <- paste(c("Season * availability", extra_fixed), collapse = " + ")
  m <- lmer(
    as.formula(paste("response ~", fixed,
                     "+ (1 | Date_f) + (1 | Date_f:availability) + (1 | Sample_cat) + (1 | Substrate:Season)")),
    data = d, control = lmerControl(optimizer = "bobyqa")
  )
  season_contrasts <- emmeans(m, ~ Season | availability) %>% contrast(method = "revpairwise")
  season_by_class <- season_contrasts %>%
    summary(infer = TRUE) %>%
    as_tibble() %>%
    transmute(availability, estimate, SE, df, lower = lower.CL, upper = upper.CL, p = p.value)
  interaction <- contrast(season_contrasts, by = NULL,
                          method = list("Complex - readily available" = c(-1, 1))) %>%
    summary(infer = TRUE) %>%
    as_tibble() %>%
    transmute(estimate, SE, df, lower = lower.CL, upper = upper.CL, p = p.value)
  list(season_by_class = season_by_class, interaction = interaction,
       singular = isSingular(m), n_obs = nrow(d))
}

primary_h1b <- fit_availability_model(sir)

# Within-class versions: the polymer (complex) against its own building
# blocks (readily available) in the same chemical class, so chemistry is held
# constant and only availability differs.
pm_data <- filter(sir, !is.na(polymer_monomer_class))
within_class <- bind_rows(
  map_dfr(sort(unique(pm_data$polymer_monomer_class)), function(cls) {
    d <- filter(pm_data, polymer_monomer_class == cls)
    fit <- fit_availability_model(d)
    fit$interaction %>%
      mutate(comparison = cls, singular = fit$singular,
             complex_substrates = paste(sort(unique(d$display_name[d$availability == "Complex"])), collapse = ", "),
             readily_substrates = paste(sort(unique(d$display_name[d$availability == "Readily available"])), collapse = ", "))
  }),
  {
    fit <- fit_availability_model(pm_data, extra_fixed = "Season * polymer_monomer_class")
    fit$interaction %>% mutate(comparison = "Pooled within-class", singular = fit$singular,
                               complex_substrates = "", readily_substrates = "")
  }
)

# Distribution-free check: random reassignment of the 8 "complex" labels
# among the 31 substrate-level seasonal effects.
class_permutation_test <- function(effects, n_perm = 99999) {
  is_complex <- effects$availability == "Complex"
  obs <- mean(effects$estimate[is_complex]) - mean(effects$estimate[!is_complex])
  perm <- replicate(n_perm, {
    idx <- sample(nrow(effects), sum(is_complex))
    mean(effects$estimate[idx]) - mean(effects$estimate[-idx])
  })
  tibble(observed_difference = obs,
         p_two_sided = (sum(abs(perm) >= abs(obs) - 1e-12) + 1) / (n_perm + 1),
         n_permutations = n_perm)
}
class_perm <- class_permutation_test(substrate_effects)

h1b <- primary_h1b$interaction
h1b_table <- bind_rows(
  h1b %>% mutate(comparison = "All substrates: complex vs readily available"),
  within_class
) %>%
  mutate(test = "Mixed model")
write_csv(primary_h1b$season_by_class, file.path(dirs$results, "Figure2D_season_effect_by_availability.csv"))
write_csv(h1b_table, file.path(dirs$results, "Figure2D_H1b_tests.csv"))
write_csv(class_perm, file.path(dirs$results, "Figure2D_H1b_permutation_check.csv"))

d_rows <- h1b_table %>%
  mutate(
    label = case_when(
      comparison == "All substrates: complex vs readily available" ~ "All substrates\n(8 complex vs 23 readily available)",
      comparison == "Pooled within-class" ~ "Pooled within-class",
      TRUE ~ sprintf("%s: %s\nvs %d small molecules", comparison, complex_substrates,
                     str_count(readily_substrates, ",") + 1)
    ),
    label = fct_rev(fct_inorder(label)),
    kind = if_else(str_detect(comparison, "All substrates|Pooled"), "Overall", "Within one class")
  )

panel_d <- ggplot(d_rows, aes(estimate, label, colour = kind)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, linewidth = 0.7,
                orientation = "y") +
  geom_point(size = 2.4) +
  geom_text(aes(x = upper, label = p_label(p)), hjust = -0.15, size = 2.5,
            show.legend = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  scale_colour_manual(values = c("Overall" = "grey10", "Within one class" = class_colors[["Complex"]])) +
  labs(
    title = "H1b: did complex substrates gain more?",
    subtitle = "Extra dormancy increase of complex over readily available substrates\n(right of 0 supports H1b)",
    x = "Difference in dormant − active change (log2)", y = NULL, colour = NULL
  ) +
  theme_manuscript +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 7),
        plot.subtitle = element_text(size = 7.5))

figure2 <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_layout(heights = c(1, 1.7)) +
  plot_annotation(tag_levels = "A")
save_figure(figure2, dirs$figures, "Figure2_main", width = 11.5, height = 12)

# -----------------------------------------------------------------------------
# Supplement S2: GAM trajectories grouped by chemistry
# -----------------------------------------------------------------------------

facet_order <- gam_stats %>%
  left_join(classes, by = "Substrate") %>%
  mutate(substrate_group = factor(substrate_group, levels = group_levels)) %>%
  arrange(substrate_group, desc(deviance_explained)) %>%
  mutate(facet = sprintf("%s: %s\n%.0f%% explained, FDR p %s", substrate_group, display_name,
                         100 * deviance_explained, format_p(p_fdr)))

s2_points <- sir %>%
  left_join(select(facet_order, Substrate, facet), by = "Substrate") %>%
  mutate(facet = factor(facet, levels = facet_order$facet))
s2_curves <- gam_curves %>%
  left_join(select(facet_order, Substrate, facet), by = "Substrate") %>%
  mutate(facet = factor(facet, levels = facet_order$facet))

figure_s2 <- ggplot(s2_points, aes(Sampling_Date, log2_sir)) +
  bands_layers() +
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
# Supplement S2b: sensitivity of the H1b contrast; S2c: group trajectories
# -----------------------------------------------------------------------------

run_sensitivity <- function(label, d) {
  fit <- fit_availability_model(d)
  fit$interaction %>% mutate(analysis = label, singular = fit$singular, n_obs = fit$n_obs)
}

sensitivity <- bind_rows(
  run_sensitivity("Primary (5 h, per mg C, centred log2)", sir),
  run_sensitivity("Root-exudate rule (amino sugars, phenolics, amines excluded)",
                  prepare_sir(sir_raw, "availability_class_sensitivity")),
  run_sensitivity("Within-substrate z-score", prepare_sir(sir_raw, standardize = "z")),
  run_sensitivity("Transition dates (Oct 2023, Apr 2024) excluded",
                  prepare_sir(filter(sir_raw, !Sampling_Date %in% transition_dates))),
  run_sensitivity("24-h incubation", prepare_sir(load_microresp_100(timepoint = 24))),
  run_sensitivity("Per g dry soil", prepare_sir(load_microresp_100(normalization = "ngCO2-C_gdw_h")))
)

h1a_sensitivity <- list(
  "Transition dates excluded" = filter(sir_raw, !Sampling_Date %in% transition_dates),
  "24-h incubation" = load_microresp_100(timepoint = 24),
  "Per g dry soil" = load_microresp_100(normalization = "ngCO2-C_gdw_h")
) %>%
  imap_dfr(function(d, label) {
    w <- to_wide(d)
    date_level_season_test(w$matrix, w$meta$Sampling_Date) %>% mutate(analysis = label)
  })

write_csv(sensitivity, file.path(dirs$results, "FigureS2b_H1b_sensitivity.csv"))
write_csv(h1a_sensitivity, file.path(dirs$results, "FigureS2b_H1a_sensitivity.csv"))

figure_s2b <- sensitivity %>%
  mutate(analysis = fct_rev(fct_inorder(analysis))) %>%
  ggplot(aes(estimate, analysis)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, orientation = "y") +
  geom_point(size = 2) +
  geom_text(aes(x = upper, label = p_label(p)), hjust = -0.15, size = 2.6) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  labs(x = "Extra dormancy increase of complex over readily available substrates (log2)\npositive supports H1b",
       y = NULL) +
  theme_manuscript
save_figure(figure_s2b, dirs$figures, "FigureS2b_H1b_sensitivity", width = 8, height = 3.5)

group_trajectories <- sir %>%
  group_by(Sample_cat, Sampling_Date, Plot_ID, Season, availability) %>%
  summarise(response = mean(response), .groups = "drop") %>%
  group_by(Sampling_Date, Season, availability) %>%
  summarise(mean = mean(response), n = n(), se = sd(response) / sqrt(n),
            lower = mean - qt(0.975, n - 1) * se, upper = mean + qt(0.975, n - 1) * se,
            .groups = "drop")
write_csv(group_trajectories, file.path(dirs$results, "FigureS2c_availability_trajectories.csv"))

figure_s2c <- ggplot(group_trajectories, aes(Sampling_Date, mean, colour = availability)) +
  bands_layers() +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = availability), alpha = 0.12, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = class_colors[availability_levels]) +
  scale_fill_manual(values = class_colors[availability_levels], guide = "none") +
  date_axis +
  labs(x = NULL, y = "log2 SIR relative to each substrate's annual mean", colour = NULL) +
  theme_manuscript +
  theme(legend.position = "bottom")
save_figure(figure_s2c, dirs$figures, "FigureS2c_availability_trajectories", width = 7, height = 4)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

summary_table <- tribble(
  ~item, ~value,
  "Samples", as.character(n_distinct(sir$Sample_cat)),
  "Sampling dates (active / dormant)", sprintf("%d (%d / %d)", primary_profile$n_units,
                                               primary_profile$n_active_dates, primary_profile$n_dormant_dates),
  "Substrates (complex / readily available)", sprintf("%d (%d / %d)", n_distinct(sir$Substrate),
                                                      sum(substrate_effects$availability == "Complex"),
                                                      sum(substrate_effects$availability == "Readily available")),
  "Level: all substrates, dormant - active (log2)",
  sprintf("%.3f [%.3f, %.3f], %s", overall_shift$estimate, overall_shift$lower,
          overall_shift$upper, p_label(overall_shift$p)),
  "Level: percent change", sprintf("%+.0f%%", overall_shift$percent_change),
  "Share of overall-level variation due to sampling date",
  sprintf("%.0f%%", 100 * variance_shares$share[variance_shares$source == "Sampling date"]),
  "Share due to plot", sprintf("%.0f%%", 100 * variance_shares$share[variance_shares$source == "Plot"]),
  "Whole profile PERMANOVA R2 / p (date level)", sprintf("%.3f / %s", primary_profile$R2, p_label(primary_profile$p_permanova)),
  "Mix only PERMANOVA R2 / p (date level)", sprintf("%.3f / %s", primary_mix$R2, p_label(primary_mix$p_permanova)),
  "Whole profile dispersion p (date level)", p_label(primary_profile$p_dispersion),
  "NMDS stress", sprintf("%.3f", nmds$stress),
  "Substrates with dormant > active point estimate",
  sprintf("%d of %d", sum(substrate_effects$estimate > 0), nrow(substrate_effects)),
  "Substrates significant after FDR", as.character(sum(substrate_effects$p_fdr < 0.05)),
  "Steady substrates (GAM FDR > 0.05)",
  paste(substrate_effects$display_name[substrate_effects$steady], collapse = "; "),
  "H1b: complex - readily available (log2)",
  sprintf("%.3f [%.3f, %.3f], %s", h1b$estimate, h1b$lower, h1b$upper, p_label(h1b$p)),
  "H1b permutation check p", p_label(class_perm$p_two_sided),
  "H1b model singular", as.character(primary_h1b$singular)
)
write_csv(summary_table, file.path(dirs$results, "Figure2_analysis_summary.csv"))
print(summary_table, n = Inf)
print(primary_h1b$season_by_class)
print(h1b_table)
print(sensitivity)
