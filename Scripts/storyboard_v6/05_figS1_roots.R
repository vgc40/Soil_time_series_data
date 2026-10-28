# =============================================================================
# Figure S1. Root biomass vs field soil CO2 by season
# =============================================================================
#
# Reuses the v4 Figure S1 model (Scripts/11_figureS1_root_respiration.R):
# CO2 ~ roots x season + (1 | date) [+ (1 | plot) when estimable], fitted to
# all 78 plot-dates and again without the points flagged by both root screens
# (Step 0.3), which are drawn as crosses.

source(file.path("Scripts", "storyboard_v6", "00_setup.R"))
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
})
emm_options(lmer.df = "satterthwaite")

roots <- read_csv(file.path(dir_v6_data, "field_100.csv"), show_col_types = FALSE) |>
  transmute(Sample_cat, Collection_date, date_f = factor(Collection_date),
            Plot = factor(Plot), season = factor(season, season_levels),
            root_biomass = Roots_gdw, respiration = Resp_rate_mg_C_m2_h, flag_both)
stopifnot(nrow(roots) == 78, !anyNA(roots$root_biomass), !anyNA(roots$respiration))
flagged <- sort(roots$Sample_cat[roots$flag_both])

fit_slopes <- function(d, label) {
  m <- lmer(respiration ~ root_biomass * season + (1 | date_f) + (1 | Plot), data = d)
  random <- "(1 | date) + (1 | plot)"
  if (isSingular(m)) {
    m <- lmer(respiration ~ root_biomass * season + (1 | date_f), data = d)
    random <- "(1 | date); plot variance estimated as zero"
  }
  emtrends(m, ~ season, var = "root_biomass") |>
    summary(infer = TRUE) |>
    as_tibble() |>
    transmute(season, slope = root_biomass.trend, lower = lower.CL, upper = upper.CL,
              p = p.value, analysis = label, n = nrow(d), random_effects = random)
}

slopes <- bind_rows(
  fit_slopes(roots, "All observations"),
  fit_slopes(filter(roots, !flag_both), paste("Without", paste(flagged, collapse = ", ")))
)

line_data <- distinct(slopes, analysis) |>
  pull() |>
  map_dfr(\(a) {
    d <- if (a == "All observations") roots else filter(roots, !flag_both)
    m <- lm(respiration ~ root_biomass * season, data = d)
    grid <- d |>
      group_by(season) |>
      reframe(root_biomass = seq(min(root_biomass), max(root_biomass), length.out = 50))
    grid |> mutate(fit = predict(m, newdata = grid), analysis = a)
  })
labels <- slopes |>
  mutate(label = sprintf("%s: %.2f (%.2f to %.2f), %s", season, slope, lower, upper,
                         fmt_p(p))) |>
  group_by(analysis) |>
  summarise(label = paste(label, collapse = "\n"), .groups = "drop")
panel_points <- bind_rows(
  roots |> mutate(analysis = "All observations"),
  roots |> filter(!flag_both) |> mutate(analysis = unique(slopes$analysis)[2])
)

s1 <- ggplot(panel_points, aes(root_biomass, respiration, colour = season)) +
  geom_point(aes(shape = flag_both), size = 1.6, alpha = 0.8) +
  geom_point(data = roots |> filter(flag_both) |> mutate(analysis = unique(slopes$analysis)[2]),
             shape = 4, size = 1.6, colour = "grey60") +
  geom_line(data = line_data, aes(y = fit), linewidth = 0.6) +
  geom_text(data = labels, aes(x = Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = 1.03, vjust = 1.2, size = 2.1) +
  facet_wrap(~ analysis) +
  scale_colour_manual(values = season_colors, name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                     labels = c("Kept", "Flagged by both root screens"), name = NULL) +
  labs(x = "Root biomass (mg g⁻¹ dry soil, 0-15 cm)",
       y = "Field soil CO2 (mg CO2-C m⁻² h⁻¹)",
       title = "Root biomass and field soil CO2 by season",
       caption = "Slopes (mg CO2-C m⁻² h⁻¹ per mg g⁻¹) from the mixed model with date as a random effect; lines are fixed-effect fits for display.\nRight panel: flagged points shown in grey, not used.") +
  theme_v6 +
  theme(legend.position = "bottom")
save_v6(s1, "FigS1", 180, 100)

write_parts(result_rows("FigS1", "",
                        paste0("root biomass slope | ", slopes$season, " | ", slopes$analysis),
                        slopes$slope, slopes$lower, slopes$upper, slopes$p,
                        n_dates = 13, n_samples = slopes$n,
                        notes = paste0("lmer CO2 ~ roots x season + ", slopes$random_effects,
                                       "; emtrends Satterthwaite")),
            "FigS1")
print(slopes)
