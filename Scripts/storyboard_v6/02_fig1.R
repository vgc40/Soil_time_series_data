# =============================================================================
# Figure 1. Study context (brief: Figure 1)
# =============================================================================
#
# A: sampling calendar. B-F: soil temperature, gravimetric water content,
# microbial biomass C (Step 0 rule applied), field soil CO2, root biomass.
# Season tests reuse the v4 Figure 1 method: a Welch t-test on date means for
# soil temperature (one value per date), and
# log(y) ~ season + (1 | date) + (1 | plot) for plot-level variables.

source(file.path("Scripts", "storyboard_v6", "00_setup.R"))
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
})

field <- read_csv(file.path(dir_v6_data, "field_100.csv"), show_col_types = FALSE) |>
  mutate(season = factor(season, season_levels))
mr <- read_csv(file.path(dir_v6_data, "microresp_long.csv"), show_col_types = FALSE)
transcript_samples <- read_csv(input_files$brite_d, n_max = 0,
                               show_col_types = FALSE) |>
  names() |> str_subset("_MT$") |> str_remove("_MT$")

# -----------------------------------------------------------------------------
# A. Sampling calendar
# -----------------------------------------------------------------------------

calendar <- bind_rows(
  field |> count(Collection_date, name = "n") |> mutate(stream = "Field and enzymes"),
  mr |> distinct(Sample_cat, Collection_date) |> count(Collection_date, name = "n") |>
    mutate(stream = "MicroResp"),
  field |> filter(Sample_cat %in% transcript_samples) |>
    count(Collection_date, name = "n") |> mutate(stream = "Transcripts")
) |>
  mutate(stream = factor(stream, c("Transcripts", "MicroResp", "Field and enzymes")))

season_labels <- tibble(
  x = as.Date(c("2023-08-10", "2024-01-01", "2024-05-05")),
  label = c("Plant-active", "Plant-dormant", "Plant-active"),
  season = c("Active", "Dormant", "Active")
)

panel_a <- ggplot(calendar, aes(Collection_date, stream)) +
  timeline_layers(marker_labels = FALSE) +
  annotate("text", x = plant_markers$date, y = 3.75, label = plant_markers$label,
           size = 2.2, colour = "grey30") +
  geom_text(data = season_labels, aes(x, y = 0.45, label = label, colour = season),
            inherit.aes = FALSE, size = 2.4, fontface = "bold") +
  geom_point(size = 4.2, colour = "grey25") +
  geom_text(aes(label = n), colour = "white", size = 2.1) +
  scale_colour_manual(values = season_colors, guide = "none") +
  scale_y_discrete(expand = expansion(add = c(0.8, 0.9))) +
  timeline_x +
  labs(title = "Sampling calendar (100% irrigation; numbers = plots)", x = NULL, y = NULL) +
  theme_v6

# -----------------------------------------------------------------------------
# B-F. Field context
# -----------------------------------------------------------------------------

variables <- tribble(
  ~column, ~label, ~unit, ~log_scale,
  "SoilTemp_mean", "Soil temperature", "°C", FALSE,
  "GWC", "Gravimetric water content", "%", TRUE,
  "MBC_used", "Microbial biomass C", "µg C g⁻¹ dry soil", TRUE,
  "Resp_rate_mg_C_m2_h", "Field soil CO2 (context only)", "mg CO2-C m⁻² h⁻¹", TRUE,
  "Roots_gdw", "Root biomass (0-15 cm)", "mg g⁻¹ dry soil", TRUE
)

field_long <- field |>
  select(Sample_cat, Collection_date, Plot, season, all_of(variables$column)) |>
  pivot_longer(all_of(variables$column), names_to = "column", values_to = "value") |>
  filter(!is.na(value))
stopifnot(all(field_long$value[field_long$column != "SoilTemp_mean"] > 0))

season_effect <- function(col, use_log) {
  d <- filter(field_long, column == col)
  if (col == "SoilTemp_mean") {
    dm <- d |> group_by(Collection_date, season) |>
      summarise(value = mean(value), .groups = "drop")
    st <- season_test(dm$value, dm$season)
    return(tibble(model = "Welch t-test on date means", estimate = st$diff,
                  ci_low = st$ci_low, ci_high = st$ci_high, p = st$p_welch,
                  ratio = FALSE, n_dates = nrow(dm), n_samples = nrow(d),
                  singular = NA))
  }
  d <- mutate(d, y = log(value), date_f = factor(Collection_date))
  m <- lmer(y ~ season + (1 | date_f) + (1 | Plot), data = d)
  co <- summary(m)$coefficients["seasonDormant", ]
  half <- qt(0.975, co[["df"]]) * co[["Std. Error"]]
  tibble(model = "lmer: log(y) ~ season + (1|date) + (1|plot)",
         estimate = exp(co[["Estimate"]]), ci_low = exp(co[["Estimate"]] - half),
         ci_high = exp(co[["Estimate"]] + half), p = co[["Pr(>|t|)"]],
         ratio = TRUE, n_dates = n_distinct(d$Collection_date), n_samples = nrow(d),
         singular = isSingular(m))
}

effects <- variables |>
  mutate(test = map2(column, log_scale, season_effect)) |>
  unnest(test) |>
  mutate(subtitle = if_else(
    ratio,
    sprintf("Dormant/active = %.2f×\n(95%% CI %.2f–%.2f), %s", estimate, ci_low,
            ci_high, fmt_p(p)),
    sprintf("Dormant − active = %.1f %s\n(95%% CI %.1f to %.1f), %s", estimate,
            unit, ci_low, ci_high, fmt_p(p))
  ))

context_panel <- function(col) {
  info <- filter(variables, column == col)
  eff <- filter(effects, column == col)
  tf <- if (info$log_scale) log10 else identity
  pts <- filter(field_long, column == col) |> mutate(y = tf(value))
  dm <- date_mean_ci(pts, "y")
  p <- ggplot(dm, aes(Collection_date, mean)) +
    timeline_layers(label_size = 1.8) +
    geom_point(data = pts, aes(y = y), colour = "grey55", alpha = 0.6, size = 0.7) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 6, linewidth = 0.35) +
    geom_line(linewidth = 0.4) +
    geom_point(aes(colour = season), size = 1.5) +
    scale_colour_manual(values = season_colors, guide = "none") +
    timeline_x_small +
    labs(title = info$label, subtitle = eff$subtitle, x = NULL,
         y = if (info$log_scale) paste(info$unit, "(log scale)") else info$unit) +
    theme_v6
  if (info$log_scale) {
    p <- p + scale_y_continuous(labels = \(b) format(signif(10^b, 2), scientific = FALSE,
                                                     drop0trailing = TRUE, trim = TRUE))
  }
  list(plot = p, dm = dm |> mutate(column = col, log_scale = info$log_scale))
}

panels <- map(variables$column, context_panel)

fig1 <- panel_a / wrap_plots(map(panels, "plot"), ncol = 3) +
  plot_layout(heights = c(0.55, 2)) +
  plot_annotation(tag_levels = "A")
save_v6(fig1, "Fig1", 180, 165)

# -----------------------------------------------------------------------------
# Results rows
# -----------------------------------------------------------------------------

panel_letters <- setNames(LETTERS[2:6], variables$column)

rows <- bind_rows(
  result_rows("Fig1", "A", paste0("plots sampled | ", calendar$stream, " | ",
                                  calendar$Collection_date),
              calendar$n, n_dates = 1, n_samples = calendar$n),
  map_dfr(panels, \(x) {
    dm <- x$dm
    back <- if (dm$log_scale[1]) \(v) 10^v else identity
    info <- filter(variables, column == dm$column[1])
    result_rows("Fig1", panel_letters[[dm$column[1]]],
                paste0(info$label, " date mean | ", dm$Collection_date),
                back(dm$mean), back(dm$ci_low), back(dm$ci_high), n_dates = 1,
                n_samples = dm$n,
                notes = if_else(dm$log_scale, "geometric mean, t-based 95% CI on log10 scale",
                                "mean, t-based 95% CI"))
  }),
  result_rows("Fig1", panel_letters[effects$column],
              if_else(effects$ratio, paste0(effects$label, " dormant/active ratio"),
                      paste0(effects$label, " dormant - active (°C)")),
              effects$estimate, effects$ci_low, effects$ci_high, effects$p,
              n_dates = effects$n_dates, n_samples = effects$n_samples,
              notes = paste0(effects$model,
                             if_else(effects$singular %in% TRUE, "; singular fit", ""),
                             if_else(effects$column == "MBC_used",
                                     "; SM1625 and SM1539 set to missing (Step 0)", "")))
)
write_parts(rows, "Fig1")
print(effects |> select(label, estimate, ci_low, ci_high, p, n_samples, singular))
