# =============================================================================
# Figure 1: study design and seasonal field context
# =============================================================================
#
# Supporting measurements used to interpret Q1-Q3 (100% irrigation, six
# plots, 13 monthly dates). Root biomass is hand-picked from pooled triplicate
# 0-15 cm cores, not rhizosphere material.
#
#   A. Sampling timeline: dates with MicroResp/field data and with
#      metatranscriptomes, season shading, and the April-October drip
#      irrigation window. Harvest and fall-saturation markers are pending
#      dates from S. Bell (decision D10).
#   B-G. Soil temperature, gravimetric water content, field soil
#      respiration, root biomass, microbial biomass C, and summed C-, N-, and
#      P-acquiring extracellular enzyme activity.
#
# Season differences are tested with date as the replication unit:
#   log(y) ~ Season + (1 | Date) + (1 | Plot) for plot-level variables;
#   a two-sample t-test on date means for soil temperature, which has one
#   value per date.

rm(list = ls())
source(file.path("Scripts", "helpers.R"))

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(patchwork)
})

dirs <- output_dirs("Figure1")

field <- load_field_100()
stopifnot(nrow(field) == 78)

availability <- read_csv(
  file.path("Results", "manuscript_v4", "sample_availability",
            "sample_availability_by_date.csv"),
  show_col_types = FALSE
)

# -----------------------------------------------------------------------------
# A. Sampling timeline
# -----------------------------------------------------------------------------

timeline <- bind_rows(
  availability %>% filter(microresp > 0) %>%
    transmute(Sampling_Date, stream = "Field, EEA, MicroResp", n = microresp),
  availability %>% filter(metatranscriptomes > 0) %>%
    transmute(Sampling_Date, stream = "Metatranscriptomes", n = metatranscriptomes)
) %>%
  mutate(stream = factor(stream, levels = c("Metatranscriptomes",
                                            "Field, EEA, MicroResp")))

irrigation_windows <- tibble(
  start = as.Date(c("2023-06-01", "2024-04-01")),
  end = as.Date(c("2023-10-31", "2024-06-05"))
)

panel_a <- ggplot(timeline, aes(Sampling_Date, stream)) +
  season_band_layers(field$Sampling_Date) +
  geom_segment(data = irrigation_windows,
               aes(x = start, xend = end, y = 0.45, yend = 0.45),
               inherit.aes = FALSE, linewidth = 2, colour = "#4F81BD") +
  annotate("text", x = as.Date("2023-06-05"), y = 0.45, label = "Drip irrigation (Apr-Oct)",
           hjust = 0, vjust = -0.8, size = 2.6, colour = "#4F81BD") +
  geom_point(aes(size = n), shape = 21, fill = "grey25", colour = "white") +
  geom_text(aes(label = n), colour = "white", size = 2.2) +
  annotate("text", x = as.Date("2023-08-15"), y = 2.55, label = "Plant-active",
           size = 3, colour = season_colors[["Plant-active"]], fontface = "bold") +
  annotate("text", x = as.Date("2024-01-15"), y = 2.55, label = "Plant-dormant",
           size = 3, colour = season_colors[["Plant-dormant"]], fontface = "bold") +
  scale_size_continuous(range = c(4, 5.5), guide = "none") +
  scale_y_discrete(expand = expansion(add = c(0.9, 0.7))) +
  date_axis +
  labs(title = "Sampling design (100% irrigation; numbers = plots sampled)",
       x = NULL, y = NULL) +
  theme_manuscript +
  theme(panel.grid.major.y = element_blank())

# -----------------------------------------------------------------------------
# B-G. Field context variables
# -----------------------------------------------------------------------------

variables <- tribble(
  ~column, ~label, ~unit, ~log_scale,
  "SoilTemp_mean", "Soil temperature", "°C", FALSE,
  "GWC", "Gravimetric water content", "%", TRUE,
  "Resp_rate_mg_C_m2_h", "Field soil respiration", "mg CO2-C m-2 h-1", TRUE,
  "Roots_gdw", "Root biomass (0-15 cm)", "mg g-1 dry soil", TRUE,
  "MBC_ug_gdw_bgc", "Microbial biomass C", "ug C g-1 dry soil", TRUE,
  "Sum_Carbon", "C-acquiring enzymes (BG + BX + CB)", "nmol g-1 h-1", TRUE,
  "Sum_Nitrogen", "N-acquiring enzymes (NAG + Leu + Ala)", "nmol g-1 h-1", TRUE,
  "Sum_Phosphorus", "P-acquiring enzyme (AP)", "nmol g-1 h-1", TRUE
)

field_long <- field %>%
  select(Sample_cat, Sampling_Date, Plot_ID, Season, all_of(variables$column)) %>%
  pivot_longer(all_of(variables$column), names_to = "column", values_to = "value") %>%
  left_join(variables, by = "column") %>%
  filter(!is.na(value))

date_summary <- field_long %>%
  group_by(column, label, unit, Sampling_Date, Season) %>%
  summarise(
    mean = mean(value), n = n(), se = sd(value) / sqrt(n),
    lower = if_else(n > 1, mean - qt(0.975, pmax(n - 1, 1)) * se, mean),
    upper = if_else(n > 1, mean + qt(0.975, pmax(n - 1, 1)) * se, mean),
    .groups = "drop"
  )
write_csv(date_summary, file.path(dirs$results, "Figure1_date_means_95CI.csv"))

season_test <- function(d, var, use_log) {
  if (var == "SoilTemp_mean") {
    dm <- d %>% group_by(Sampling_Date, Season) %>% summarise(value = mean(value), .groups = "drop")
    tt <- t.test(value ~ Season, data = dm)
    return(tibble(
      model = "Welch t-test on date means", scale = "original",
      estimate = tt$estimate[[2]] - tt$estimate[[1]], lower = -tt$conf.int[2],
      upper = -tt$conf.int[1], p = tt$p.value,
      ratio = NA_real_, n_dates = nrow(dm), singular = NA
    ))
  }
  d <- mutate(d, y = if (use_log) log(value) else value, Date_f = factor(Sampling_Date))
  m <- lmer(y ~ Season + (1 | Date_f) + (1 | Plot_ID), data = d)
  co <- summary(m)$coefficients["SeasonPlant-dormant", ]
  est <- co[["Estimate"]]
  half <- qt(0.975, co[["df"]]) * co[["Std. Error"]]
  tibble(
    model = "lmer: y ~ Season + (1|Date) + (1|Plot)",
    scale = if (use_log) "log" else "original",
    estimate = est, lower = est - half, upper = est + half, p = co[["Pr(>|t|)"]],
    ratio = if (use_log) exp(est) else NA_real_,
    n_dates = n_distinct(d$Sampling_Date), singular = isSingular(m)
  )
}

season_tests <- variables %>%
  mutate(test = pmap(list(column, log_scale), function(col, ls) {
    season_test(filter(field_long, column == col), col, ls)
  })) %>%
  unnest(test) %>%
  mutate(
    effect_label = if_else(
      is.na(ratio),
      sprintf("Dormant - active = %.1f %s, %s", estimate, unit, p_label(p)),
      sprintf("Dormant/active = %.2f×, %s", ratio, p_label(p))
    )
  )
write_csv(season_tests, file.path(dirs$results, "Figure1_season_tests.csv"))

# Log-scaled variables are summarised as geometric means with 95% CIs, and
# plotted as log10 values on a linear axis with back-transformed labels so
# the season shading can span the full panel.
context_panel <- function(col) {
  info <- filter(variables, column == col)
  test <- filter(season_tests, column == col)
  tf <- if (info$log_scale) log10 else identity
  pts <- filter(field_long, column == col) %>% mutate(y = tf(value))
  dm <- pts %>%
    group_by(Sampling_Date, Season) %>%
    summarise(mean = mean(y), n = n(), se = sd(y) / sqrt(n), .groups = "drop") %>%
    mutate(lower = if_else(n > 1, mean - qt(0.975, pmax(n - 1, 1)) * se, mean),
           upper = if_else(n > 1, mean + qt(0.975, pmax(n - 1, 1)) * se, mean))
  p <- ggplot(dm, aes(Sampling_Date, mean)) +
    season_band_layers(field$Sampling_Date) +
    geom_point(data = pts, aes(y = y), colour = "grey45", alpha = 0.45, size = 0.9) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey20", alpha = 0.15) +
    geom_line(linewidth = 0.6) +
    geom_point(aes(colour = Season), size = 1.6) +
    scale_colour_manual(values = season_colors, guide = "none") +
    scale_x_date(date_breaks = "2 months", date_labels = "%b\n%y") +
    labs(title = info$label, subtitle = test$effect_label, x = NULL,
         y = if (info$log_scale) paste(info$unit, "(log scale)") else info$unit) +
    theme_manuscript +
    theme(plot.title = element_text(size = 9), plot.subtitle = element_text(size = 7))
  if (info$log_scale) {
    p <- p + scale_y_continuous(labels = function(b) format(signif(10^b, 2),
                                                            scientific = FALSE,
                                                            drop0trailing = TRUE))
  }
  p
}

context_panels <- map(variables$column, context_panel)

figure1 <- panel_a / wrap_plots(context_panels, ncol = 4) +
  plot_layout(heights = c(0.45, 2)) +
  plot_annotation(tag_levels = "A")
save_figure(figure1, dirs$figures, "Figure1_main", width = 12, height = 8)

print(select(season_tests, label, model, estimate, lower, upper, ratio, p, singular))
