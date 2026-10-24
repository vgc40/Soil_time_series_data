# =============================================================================
# Figure S1: root biomass and field soil respiration (supporting analysis)
# =============================================================================
#
# Not a storyboard question (decision D3); kept as context for interpreting
# field respiration, which integrates root and microbial respiration.
#
# Model: respiration ~ root biomass x Season + (1 | Date) [+ (1 | Plot) when
# its variance is estimable], fitted to all 78 plot-dates (decision D11).
# Sensitivity: exclude the three observations flagged as root-biomass
# outliers by both a global and a within-plot modified z-score screen
# (|z| > 3.5; median/MAD): SM1621, SM1078, SM1080.

rm(list = ls())
source(file.path("Scripts", "helpers.R"))

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
})

emm_options(lmer.df = "satterthwaite")
dirs <- output_dirs("FigureS1")

roots <- load_field_100() %>%
  transmute(Sample_cat, Sampling_Date, Date_f = factor(Sampling_Date), Plot_ID,
            Season, root_biomass = Roots_gdw,
            respiration = Resp_rate_mg_C_m2_h)
stopifnot(nrow(roots) == 78, !anyNA(roots$root_biomass), !anyNA(roots$respiration))

# -----------------------------------------------------------------------------
# Outlier screens (modified z-score, Iglewicz & Hoaglin threshold 3.5)
# -----------------------------------------------------------------------------

modified_z <- function(x) 0.6745 * (x - median(x)) / mad(x, constant = 1)

screened <- roots %>%
  mutate(z_global = modified_z(root_biomass)) %>%
  group_by(Plot_ID) %>%
  mutate(z_within_plot = modified_z(root_biomass)) %>%
  ungroup() %>%
  mutate(
    flag_global = abs(z_global) > 3.5,
    flag_within_plot = abs(z_within_plot) > 3.5,
    flag_both = flag_global & flag_within_plot
  )
write_csv(screened, file.path(dirs$results, "FigureS1_root_outlier_screens.csv"))

both_flagged <- screened$Sample_cat[screened$flag_both]
expected_flags <- c("SM1078", "SM1080", "SM1621")
if (!setequal(both_flagged, expected_flags)) {
  warning("Both-screen flags differ from the 3 agreed in D11: ",
          paste(both_flagged, collapse = ", "))
}

# -----------------------------------------------------------------------------
# Season-specific slopes
# -----------------------------------------------------------------------------

fit_slopes <- function(d, label) {
  m <- lmer(respiration ~ root_biomass * Season + (1 | Date_f) + (1 | Plot_ID), data = d)
  random <- "(1 | Date) + (1 | Plot)"
  if (isSingular(m)) {
    m <- lmer(respiration ~ root_biomass * Season + (1 | Date_f), data = d)
    random <- "(1 | Date); plot variance estimated as zero"
  }
  slopes <- emtrends(m, ~ Season, var = "root_biomass") %>%
    summary(infer = TRUE) %>%
    as_tibble() %>%
    transmute(Season, slope = root_biomass.trend, lower = lower.CL,
              upper = upper.CL, p = p.value)
  interaction <- summary(m)$coefficients["root_biomass:SeasonPlant-dormant", ]
  slopes %>% mutate(
    analysis = label, n = nrow(d), random_effects = random,
    interaction_estimate = interaction[["Estimate"]],
    interaction_p = interaction[["Pr(>|t|)"]]
  )
}

slope_results <- bind_rows(
  fit_slopes(screened, "All observations (primary)"),
  fit_slopes(filter(screened, !flag_both),
             paste("Excluding", paste(sort(both_flagged), collapse = ", ")))
)
write_csv(slope_results, file.path(dirs$results, "FigureS1_root_respiration_slopes.csv"))

# -----------------------------------------------------------------------------
# Figure
# -----------------------------------------------------------------------------

line_data <- slope_results %>%
  distinct(analysis) %>%
  pull(analysis) %>%
  map_dfr(function(a) {
    d <- if (str_detect(a, "primary")) screened else filter(screened, !flag_both)
    m <- lm(respiration ~ root_biomass * Season, data = d)
    grid <- d %>%
      group_by(Season) %>%
      reframe(root_biomass = seq(min(root_biomass), max(root_biomass),
                                 length.out = 50))
    grid %>% mutate(fit = predict(m, newdata = grid), analysis = a)
  })

slope_labels <- slope_results %>%
  mutate(label = sprintf("%s: %.2f (%.2f to %.2f), %s", Season, slope, lower,
                         upper, p_label(p))) %>%
  group_by(analysis) %>%
  summarise(label = paste(label, collapse = "\n"), .groups = "drop")

figure_s1 <- ggplot(screened, aes(root_biomass, respiration, colour = Season)) +
  geom_point(aes(shape = flag_both), size = 2, alpha = 0.8) +
  geom_line(data = line_data, aes(y = fit), linewidth = 0.8) +
  geom_text(data = slope_labels, aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = 1.02, vjust = 1.2, size = 2.6) +
  facet_wrap(~ analysis) +
  scale_colour_manual(values = season_colors) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                     labels = c("Retained", "Flagged by both screens")) +
  labs(x = "Root biomass (mg g-1 dry soil, 0-15 cm)",
       y = "Field soil respiration (mg CO2-C m-2 h-1)",
       colour = NULL, shape = NULL,
       caption = "Lines: fixed-effect slopes for display; reported slopes from the mixed model.") +
  theme_manuscript +
  theme(legend.position = "bottom")
save_figure(figure_s1, dirs$figures, "FigureS1_root_respiration", width = 10, height = 5)

print(slope_results)
