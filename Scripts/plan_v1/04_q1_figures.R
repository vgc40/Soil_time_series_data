# =============================================================================
# 04. Q1 figures (plan §3.4)
# =============================================================================
#
# Main figures 1-3 and supplementary figures, each saved as PDF and 300 dpi
# PNG with a caption file. Reads the tables written by scripts 01-03.

source(file.path("Scripts", "plan_v1", "00_setup.R"))
suppressPackageStartupMessages({
  library(patchwork)
  library(mgcv)
})

read_res <- function(name) {
  read_csv(file.path(dir_results, name), show_col_types = FALSE)
}

date_scores <- read_res("date_scores_100pct_5h.csv") |>
  mutate(season = factor(season, levels = c("Active", "Dormant")))
sample_scores <- read_res("sample_scores_100pct_5h.csv") |>
  mutate(season = factor(season, levels = c("Active", "Dormant")))
season_tests <- read_res("cp2_t1_t3_t6_season_tests.csv")
roots <- read_res("root_flags_active_100pct.csv")
compound_effects <- read_res("cp2_t5_compound_effects.csv")

fmt_p <- function(p) if_else(p < 0.001, "p < 0.001", paste0("p = ", signif(p, 2)))

# Mean and Welch-style 95% CI for each season's date means
season_summary <- function(df, value) {
  df |>
    group_by(season) |>
    summarise(n = n(), mean = mean(.data[[value]]), sd = sd(.data[[value]]),
              .groups = "drop") |>
    mutate(half = qt(0.975, n - 1) * sd / sqrt(n),
           ci_low = mean - half, ci_high = mean + half)
}

# Dates whose mean includes a sample with a replaced or lowest MBC (D1)
mbc_flag_dates <- sample_scores |>
  filter(normalization == "per_MBC", MBC_flag != "none") |>
  distinct(Collection_date) |>
  pull()

# -----------------------------------------------------------------------------
# Figure 1. Season comparison of the simple score
# -----------------------------------------------------------------------------
# Plan: add the biomass-normalized result if it differs from the per-C result
# in direction or significance (p < 0.05, Welch).

t1_C <- season_tests |> filter(test == "1_H1a_simple_score")
t1_MBC <- season_tests |> filter(test == "6_H1a_simple_score")
add_mbc <- sign(t1_C$diff) != sign(t1_MBC$diff) |
  (t1_C$p_welch < 0.05) != (t1_MBC$p_welch < 0.05)

fig1_norms <- if (add_mbc) c("per_C", "per_MBC") else "per_C"
norm_labels <- c(per_C = "Per mg added C", per_MBC = "Per mg added C per µg MBC")

fig1_points <- date_scores |>
  filter(normalization %in% fig1_norms) |>
  mutate(mbc_flag = normalization == "per_MBC" & Collection_date %in% mbc_flag_dates,
         panel = norm_labels[normalization])
fig1_summary <- fig1_points |>
  group_by(normalization, panel) |>
  group_modify(\(d, k) season_summary(d, "simple_score")) |>
  ungroup()
fig1_labels <- season_tests |>
  filter(test %in% c("1_H1a_simple_score", "6_H1a_simple_score"),
         normalization %in% fig1_norms) |>
  mutate(panel = norm_labels[normalization],
         label = sprintf("Dormant − active = %.2f\n(95%% CI %.2f to %.2f)\nWelch %s; Wilcoxon %s",
                         diff, ci_low, ci_high, fmt_p(p_welch), fmt_p(p_wilcox)))

fig1 <- ggplot(fig1_points, aes(season, simple_score)) +
  geom_point(aes(colour = season, shape = mbc_flag),
             position = position_jitter(width = 0.08, seed = 2026),
             size = 1.8, alpha = 0.8) +
  geom_errorbar(data = fig1_summary, aes(y = mean, ymin = ci_low, ymax = ci_high),
                width = 0.15, linewidth = 0.5, position = position_nudge(x = 0.28)) +
  geom_point(data = fig1_summary, aes(y = mean), size = 2.4,
             position = position_nudge(x = 0.28)) +
  geom_text(data = fig1_labels, aes(x = 0.5, y = Inf, label = label),
            hjust = 0, vjust = 1.2, size = 2.5, lineheight = 0.9) +
  scale_colour_manual(values = season_colors, guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1),
                     labels = c(`FALSE` = "Date mean",
                                `TRUE` = "Date includes a replaced or lowest MBC sample"),
                     name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.45))) +
  facet_wrap(~panel, scales = "free_y") +
  labs(x = NULL, y = "Simple-compound score (mean z)") +
  theme_plan +
  theme(legend.position = "bottom")

fig1_width <- if (add_mbc) 180 else 85
save_plan_figure(
  fig1, "Fig1_season_simple_score", fig1_width, 85,
  paste0(
    "Figure 1. Simple-compound score by season (H1a, test 1). ",
    "Points are date means of the sample-level simple score (mean z-score of ",
    "15 simple compounds; water-subtracted, 5 h, 100% irrigation). ",
    "n = 13 dates (8 active, 5 dormant), 77 samples. ",
    "Black points and bars: season mean of date means and its 95% t-interval. ",
    "Annotation: dormant minus active difference with Welch 95% CI, Welch and ",
    "Wilcoxon rank-sum p-values. ",
    if (add_mbc) paste0(
      "Right panel: the same score after dividing each value by microbial ",
      "biomass C, shown because its result differs from the per-C result in ",
      "significance (plan §3.4). Open circles: dates whose mean includes a ",
      "sample with MBC replaced (0 set to 0.07 µg C/g) or the lowest measured ",
      "MBC (0.14 µg C/g); see deviation log D1."
    ) else ""
  )
)

# -----------------------------------------------------------------------------
# Figure 2. Simple score vs root biomass (active season)
# -----------------------------------------------------------------------------

root_tests <- read_res("cp2_t2_root_biomass_score.csv")
fig2_data <- sample_scores |>
  filter(normalization == "per_C", season == "Active") |>
  inner_join(roots |> select(Sample_cat, Roots_gdw, root_flag), by = "Sample_cat")

fig2_label <- sprintf(
  "All cores: rho = %.2f, %s (n = %d)\nFlagged removed: rho = %.2f, %s (n = %d)",
  root_tests$rho[1], fmt_p(root_tests$p[1]), root_tests$n[1],
  root_tests$rho[2], fmt_p(root_tests$p[2]), root_tests$n[2]
)

fig2 <- ggplot(fig2_data, aes(Roots_gdw, simple_score)) +
  geom_smooth(aes(linetype = "All cores"), method = "lm", formula = y ~ x,
              se = FALSE, colour = "grey40", linewidth = 0.5) +
  geom_smooth(data = filter(fig2_data, !root_flag),
              aes(linetype = "Flagged cores removed"), method = "lm",
              formula = y ~ x, se = FALSE, colour = group_colors[["Simple"]],
              linewidth = 0.5) +
  geom_point(aes(shape = root_flag), colour = group_colors[["Simple"]], size = 1.6) +
  annotate("text", x = Inf, y = Inf, label = fig2_label, hjust = 1.05,
           vjust = 1.3, size = 2.4, lineheight = 0.9) +
  scale_x_log10() +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4),
                     labels = c(`FALSE` = "Core", `TRUE` = "Root-dense (|M| > 3.5)"),
                     name = NULL) +
  scale_linetype_manual(values = c("All cores" = "dashed",
                                   "Flagged cores removed" = "solid"),
                        name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  labs(x = "Root biomass (Roots_gdw, log scale)",
       y = "Simple-compound score (mean z)") +
  theme_plan +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.margin = margin(0, 0, 0, 0))

save_plan_figure(
  fig2, "Fig2_simple_score_root_biomass", 85, 95,
  paste0(
    "Figure 2. Simple-compound score against root biomass in the plant-active ",
    "season (H1a, test 2). Points are samples (100% irrigation, 5 h, ",
    "water-subtracted per-C values); n = 48 active-season samples on 8 dates. ",
    "Crosses: root-dense cores flagged by modified z-score |M| > 3.5 (n = 6; ",
    "SM1080, SM1078, SM1621, SM1126, SM1054, SM1150). Lines are least-squares ",
    "fits for display, with all cores (dashed) and without flagged cores ",
    "(solid). Statistics are Spearman rank correlations. x axis is log10. ",
    "No error bars."
  )
)

# -----------------------------------------------------------------------------
# Figure 3. Forest plot of compound-level season effects
# -----------------------------------------------------------------------------

summary_rows <- season_tests |>
  filter(normalization == "per_C") |>
  transmute(
    Substrate = case_when(
      test == "1_H1a_simple_score" ~ "Simple-group mean",
      test == "polymer_score (for forest plot)" ~ "Polymer-group mean",
      test == "3_H1b_polymer_minus_simple" ~ "Polymer − simple (paired)"
    ),
    group = "Summary", diff, ci_low, ci_high, p_welch
  )

forest <- bind_rows(
  compound_effects |> select(Substrate, group, diff, ci_low, ci_high, p_welch),
  summary_rows
) |>
  mutate(
    block = factor(group, levels = c("Simple", "Polymer", "Reference", "Summary"),
                   labels = c("Simple", "Polymer", "Ref.", "Summary")),
    label = str_replace_all(Substrate, c("_" = " ", "^D\\." = "D-", "^L\\." = "L-",
                                         "^N\\.Acetyl\\.D\\." = "N-acetyl-D-"))
  ) |>
  group_by(block) |>
  arrange(if_else(block == "Summary", 0, diff), .by_group = TRUE) |>
  ungroup() |>
  mutate(label = factor(label, levels = unique(label)))

summary_order <- c("Polymer − simple (paired)", "Polymer-group mean", "Simple-group mean")
forest <- forest |>
  mutate(label = fct_relevel(label, summary_order, after = 0))

fig3 <- ggplot(forest, aes(diff, label, colour = group)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, linewidth = 0.5,
                orientation = "y") +
  geom_point(size = 1.8) +
  facet_grid(block ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = group_colors, guide = "none") +
  labs(x = "Dormant − active (z-score units, Welch 95% CI)", y = NULL) +
  theme_plan +
  theme(strip.text.y = element_text(angle = 0))

save_plan_figure(
  fig3, "Fig3_forest_compound_effects", 85, 120,
  paste0(
    "Figure 3. Season effect for each compound (test 5) and for the group ",
    "scores (tests 1 and 3). Each compound's water-subtracted per-C value (5 h, ",
    "100% irrigation) was z-scored across 77 samples and averaged to 13 date ",
    "means (8 active, 5 dormant). Points: dormant minus active difference of ",
    "date means; bars: Welch 95% t-interval. Blue: simple compounds (n = 15); ",
    "orange: polymers (n = 5); dark grey: glucose (reference). Summary rows: ",
    "simple-group mean score (test 1), polymer-group mean score, and the paired ",
    "polymer − simple difference computed per date (test 3). Compounds are ",
    "sorted by effect within group. No multiple-comparison correction (see ",
    "cp2_t5_compound_effects.csv for FDR-adjusted p-values)."
  )
)

# -----------------------------------------------------------------------------
# Figure S1. Enzymes: season comparison and enzyme-substrate correlations
# -----------------------------------------------------------------------------

enzyme_dates <- read_res("cp2_t7_enzyme_date_means.csv") |>
  mutate(season = factor(season, levels = c("Active", "Dormant")))
enzyme_season <- read_res("cp2_t7_enzyme_season.csv")
enzyme_substrate <- read_res("cp2_t7_enzyme_substrate.csv")
compound_dates <- read_res("cp2_t5_compound_date_means.csv")

enz_summary <- enzyme_dates |>
  group_by(enzyme) |>
  group_modify(\(d, k) season_summary(d, "activity")) |>
  ungroup()
enz_labels <- enzyme_season |>
  mutate(label = sprintf("ratio %.2f\nWelch %s", ratio, fmt_p(p_welch)))

s1a <- ggplot(enzyme_dates, aes(season, activity)) +
  geom_point(aes(colour = season), position = position_jitter(width = 0.08, seed = 2026),
             size = 1.5) +
  geom_errorbar(data = enz_summary, aes(y = mean, ymin = ci_low, ymax = ci_high),
                width = 0.15, position = position_nudge(x = 0.28)) +
  geom_point(data = enz_summary, aes(y = mean), size = 2,
             position = position_nudge(x = 0.28)) +
  geom_text(data = enz_labels, aes(x = 0.5, y = Inf, label = label), hjust = 0,
            vjust = 1.2, size = 2.3, lineheight = 0.9) +
  facet_wrap(~enzyme, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = season_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(x = NULL, y = "Activity (nmol/g dry soil/h)", tag = "A") +
  theme_plan

enz_pairs <- enzyme_substrate |>
  pmap(\(enzyme, Substrate, ...) {
    inner_join(enzyme_dates[enzyme_dates$enzyme == enzyme, ],
               compound_dates[compound_dates$Substrate == Substrate,
                              c("Collection_date", "z")],
               by = "Collection_date") |>
      mutate(pair = paste0(enzyme, " vs ", Substrate))
  }) |>
  list_rbind()
pair_labels <- enzyme_substrate |>
  mutate(pair = paste0(enzyme, " vs ", Substrate),
         label = sprintf("rho = %.2f\n%s", rho, fmt_p(p)))

s1b <- ggplot(enz_pairs, aes(activity, z)) +
  geom_point(aes(colour = season), size = 1.5) +
  geom_text(data = pair_labels, aes(x = -Inf, y = Inf, label = label), hjust = -0.1,
            vjust = 1.2, size = 2.3, lineheight = 0.9) +
  facet_wrap(~pair, scales = "free", nrow = 1) +
  scale_colour_manual(values = season_colors, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(x = "Enzyme activity, date mean (nmol/g dry soil/h)",
       y = "Substrate z, date mean", tag = "B") +
  theme_plan +
  theme(legend.position = "bottom")

save_plan_figure(
  s1a / s1b, "FigS1_enzymes", 180, 120,
  paste0(
    "Figure S1. Potential enzyme activities and their matching substrates ",
    "(test 7). (A) Date means of CB, BG, BX and NAG activity by season (100% ",
    "irrigation; n = 13 dates, 8 active and 5 dormant, 78 samples); black ",
    "points and bars are the season mean of date means and its 95% t-interval; ",
    "annotations give the dormant/active ratio and Welch p-value. (B) Date ",
    "means of enzyme activity against the date-mean z-score of the matching ",
    "substrate (5 h, water-subtracted per C); rho is the Spearman correlation ",
    "across 13 dates. No error bars in B."
  )
)

# -----------------------------------------------------------------------------
# Figure S2. Later-to-early rate ratio by season and group
# -----------------------------------------------------------------------------

ratio_dates <- read_res("cp2_t8_rate_ratio_dates.csv") |>
  mutate(season = factor(season, levels = c("Active", "Dormant")),
         group = factor(group, levels = c("Simple", "Polymer")))
ratio_tests <- read_res("cp2_t8_rate_ratio_season.csv") |>
  mutate(group = factor(group, levels = c("Simple", "Polymer")),
         label = sprintf("Dormant − active = %.2f\nWelch %s",
                         diff, fmt_p(p_welch)))
ratio_summary <- ratio_dates |>
  group_by(group) |>
  group_modify(\(d, k) season_summary(d, "ratio_24_5")) |>
  ungroup()

s2 <- ggplot(ratio_dates, aes(season, ratio_24_5)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_point(aes(colour = group), position = position_jitter(width = 0.08, seed = 2026),
             size = 1.6) +
  geom_errorbar(data = ratio_summary, aes(y = mean, ymin = ci_low, ymax = ci_high),
                width = 0.15, position = position_nudge(x = 0.28)) +
  geom_point(data = ratio_summary, aes(y = mean), size = 2,
             position = position_nudge(x = 0.28)) +
  geom_text(data = ratio_tests, aes(x = 0.5, y = Inf, label = label), hjust = 0,
            vjust = 1.2, size = 2.4, lineheight = 0.9) +
  facet_wrap(~group) +
  scale_colour_manual(values = group_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  labs(x = NULL, y = "Rate at 24 h / rate at 5 h") +
  theme_plan

save_plan_figure(
  s2, "FigS2_rate_ratio_24h_5h", 85, 75,
  paste0(
    "Figure S2. Ratio of the 24 h to the 5 h water-subtracted per-C rate by ",
    "season and group (test 8). For each date and group, the mean rate over ",
    "all plots and group compounds at 24 h was divided by the same mean at 5 h ",
    "(deviation log D6). Whether the 24 h value is the rate over 5-24 h or ",
    "0-24 h is not stated in the data package (open item O1). Points: dates ",
    "(n = 13; 8 active, 5 dormant; 100% irrigation). Black points and bars: ",
    "season mean and 95% t-interval. Dashed line: ratio of 1."
  )
)

# -----------------------------------------------------------------------------
# Figure S3. Time series of the simple and polymer scores with GAM smooths
# -----------------------------------------------------------------------------
# Description only: s(day, k = 6) with plot as a random effect.

ts_samples <- sample_scores |>
  filter(normalization == "per_C") |>
  pivot_longer(c(simple_score, polymer_score), names_to = "score", values_to = "value") |>
  mutate(group = if_else(score == "simple_score", "Simple", "Polymer"),
         group = factor(group, levels = c("Simple", "Polymer")),
         day = as.numeric(Collection_date - min(study_dates)),
         Plot = factor(Plot))

grid <- tibble(Collection_date = seq(min(study_dates), max(study_dates), by = "3 days")) |>
  mutate(day = as.numeric(Collection_date - min(study_dates)))

gam_fits <- ts_samples |>
  group_by(group) |>
  group_modify(\(d, k) {
    fit <- gam(value ~ s(day, k = 6) + s(Plot, bs = "re"), data = d, method = "REML")
    nd <- grid |> mutate(Plot = d$Plot[1])
    pr <- predict(fit, nd, exclude = "s(Plot)", se.fit = TRUE)
    nd |> mutate(fit = pr$fit, low = pr$fit - 1.96 * pr$se.fit,
                 high = pr$fit + 1.96 * pr$se.fit)
  }) |>
  ungroup()

gam_summary <- ts_samples |>
  group_by(group) |>
  group_modify(\(d, k) {
    fit <- gam(value ~ s(day, k = 6) + s(Plot, bs = "re"), data = d, method = "REML")
    st <- summary(fit)$s.table
    tibble(edf_day = st["s(day)", "edf"], p_day = st["s(day)", "p-value"],
           dev_expl = summary(fit)$dev.expl, n = nrow(d))
  }) |>
  ungroup()
write_csv(gam_summary, file.path(dir_results, "cp2_gam_time_series.csv"))

ts_dates <- date_scores |>
  filter(normalization == "per_C") |>
  pivot_longer(c(simple_score, polymer_score), names_to = "score", values_to = "value") |>
  mutate(group = factor(if_else(score == "simple_score", "Simple", "Polymer"),
                        levels = c("Simple", "Polymer")))

s3 <- ggplot() +
  dormant_rect() +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(data = ts_samples, aes(Collection_date, value, colour = group),
             size = 0.9, alpha = 0.35) +
  geom_ribbon(data = gam_fits, aes(Collection_date, ymin = low, ymax = high, fill = group),
              alpha = 0.2) +
  geom_line(data = gam_fits, aes(Collection_date, fit, colour = group), linewidth = 0.7) +
  geom_point(data = ts_dates, aes(Collection_date, value, colour = group), size = 1.8) +
  scale_colour_manual(values = group_colors, name = NULL) +
  scale_fill_manual(values = group_colors, guide = "none") +
  date_axis +
  labs(x = NULL, y = "Score (mean z)") +
  theme_plan +
  theme(legend.position = "bottom")

save_plan_figure(
  s3, "FigS3_score_time_series", 180, 85,
  paste0(
    "Figure S3. Simple and polymer scores over the 13 sampling dates ",
    "(description only). Small points: samples (n = 77; 100% irrigation, 5 h, ",
    "water-subtracted per C); large points: date means. Lines and bands: GAM ",
    "smooth s(day, k = 6) with plot as a random effect, fitted to samples and ",
    "shown at the population level with approximate 95% confidence bands. ",
    "Grey shading: plant-dormant period (5 dates)."
  )
)

# -----------------------------------------------------------------------------
# Figure S4 (diagnostic, not in plan): biomass-normalized sample scores
# -----------------------------------------------------------------------------
# Shows the samples with replaced or lowest MBC (D1) so their weight in the
# per-MBC scores is visible.

mbc_samples <- sample_scores |>
  filter(normalization == "per_MBC") |>
  pivot_longer(c(simple_score, polymer_score), names_to = "score", values_to = "value") |>
  mutate(group = factor(if_else(score == "simple_score", "Simple", "Polymer"),
                        levels = c("Simple", "Polymer")),
         MBC_flag = factor(MBC_flag, levels = c("none", "lowest_measured", "replaced_zero")))

s4 <- ggplot(mbc_samples, aes(Collection_date, value)) +
  dormant_rect() +
  geom_point(aes(shape = MBC_flag, colour = group), size = 1.6) +
  geom_text(data = filter(mbc_samples, MBC_flag != "none"),
            aes(label = Sample_cat), size = 2.2, nudge_y = 0.5) +
  facet_wrap(~group, nrow = 1) +
  scale_colour_manual(values = group_colors, guide = "none") +
  scale_shape_manual(values = c(none = 16, lowest_measured = 2, replaced_zero = 1),
                     labels = c(none = "MBC measured",
                                lowest_measured = "Lowest measured MBC (0.14)",
                                replaced_zero = "MBC 0 replaced with 0.07"),
                     name = NULL) +
  date_axis +
  labs(x = NULL, y = "Biomass-normalized score (mean z)") +
  theme_plan +
  theme(legend.position = "bottom")

save_plan_figure(
  s4, "FigS4_mbc_normalized_scores_flags", 180, 80,
  paste0(
    "Figure S4 (diagnostic; not in the analysis plan). Biomass-normalized ",
    "simple and polymer scores for each sample (n = 77; 100% irrigation, 5 h). ",
    "Values were divided by microbial biomass C before z-scoring (plan §2.7). ",
    "Triangle: SM1539, the lowest measured MBC (0.14 µg C/g). Open circle: ",
    "SM1625, MBC of 0 replaced with 0.07 µg C/g (deviation log D1). ",
    "Grey shading: plant-dormant period. No error bars."
  )
)

cat("Figure 1 includes per-MBC panel:", add_mbc, "\n")
print(gam_summary)
