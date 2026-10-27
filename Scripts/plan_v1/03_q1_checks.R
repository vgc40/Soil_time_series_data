# =============================================================================
# 03. Q1 checks, exploratory analyses, and 25% replication (plan §3.3-3.6)
# =============================================================================
#
# Mixed models, sliding-window permutation, robustness variants for tests 1
# and 3, exploratory timing and field correlations, and tests 1 and 3 on the
# 25% irrigation plots.

source(file.path("Scripts", "plan_v1", "00_setup.R"))
suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
})

sir <- read_csv(file.path(dir_results, "sir_prepared_long.csv"),
                show_col_types = FALSE) |>
  mutate(season = factor(season, levels = c("Active", "Dormant")),
         # Uncorrected per-C value, recomputed with the mgC_per_well column (D5)
         SIR_C = SIR_gdw * gdw_well / mgC_per_well)
sir100_5 <- sir |> filter(Irrigation_trmt_pct == 100, Timepoint == 5)

date_scores <- read_csv(file.path(dir_results, "date_scores_100pct_5h.csv"),
                        show_col_types = FALSE)
dates_C <- date_scores |> filter(normalization == "per_C") |> arrange(Collection_date)
water_summary <- read_csv(file.path(dir_results, "cp1_water_subtraction_summary.csv"),
                          show_col_types = FALSE)

# -----------------------------------------------------------------------------
# §3.3 Mixed models
# -----------------------------------------------------------------------------

z_long <- make_z_long(sir100_5, "induced_C",
                      c(simple_compounds, polymer_compounds)) |>
  left_join(compound_groups |> rename(cgroup = group), by = "Substrate") |>
  mutate(group = factor(cgroup, levels = c("Simple", "Polymer")),
         date = factor(Collection_date), Plot = factor(Plot))

fit_h1a <- lmer(z ~ season + (1 | Substrate) + (1 | Plot) + (1 | date),
                data = filter(z_long, group == "Simple"))
fit_h1b <- lmer(z ~ season * group + (1 | Substrate) + (1 | Plot) + (1 | date),
                data = z_long)

tidy_fixed <- function(fit, model) {
  co <- summary(fit)$coefficients
  ci <- confint(fit, parm = "beta_", method = "Wald")
  tibble(model = model, term = rownames(co), estimate = co[, "Estimate"],
         se = co[, "Std. Error"], df = co[, "df"], t = co[, "t value"],
         p = co[, "Pr(>|t|)"], ci_low = ci[, 1], ci_high = ci[, 2],
         singular = isSingular(fit), n_obs = nobs(fit))
}
mixed <- bind_rows(tidy_fixed(fit_h1a, "H1a: simple z ~ season"),
                   tidy_fixed(fit_h1b, "H1b: z ~ season * group"))
write_csv(mixed, file.path(dir_results, "cp2_mixed_models.csv"))

mixed_var <- bind_rows(
  as.data.frame(VarCorr(fit_h1a)) |> mutate(model = "H1a"),
  as.data.frame(VarCorr(fit_h1b)) |> mutate(model = "H1b")
) |> select(model, grp, vcov, sdcor)
write_csv(mixed_var, file.path(dir_results, "cp2_mixed_models_variance.csv"))

# -----------------------------------------------------------------------------
# §3.3 Sliding-window permutation (no wrap-around)
# -----------------------------------------------------------------------------

window_shift <- function(values, dates, width = 5) {
  n <- length(values)
  starts <- seq_len(n - width + 1)
  obs_in <- dates %in% dormant_dates
  stat <- function(inside) mean(values[inside]) - mean(values[!inside])
  d <- map_dbl(starts, \(s) stat(seq_len(n) %in% s:(s + width - 1)))
  observed <- stat(obs_in)
  tibble(n_positions = length(starts),
         observed_start = dates[which(obs_in)[1]],
         observed_diff = observed,
         rank_abs = sum(abs(d) >= abs(observed)),
         p = mean(abs(d) >= abs(observed)),
         min_possible_p = 1 / length(starts),
         window_diffs = paste(round(d, 3), collapse = "; "))
}
window_results <- bind_rows(
  window_shift(dates_C$simple_score, dates_C$Collection_date) |>
    mutate(score = "simple_score"),
  window_shift(dates_C$diff_score, dates_C$Collection_date) |>
    mutate(score = "polymer_minus_simple")
) |> relocate(score)
write_csv(window_results, file.path(dir_results, "cp2_window_permutation.csv"))

# -----------------------------------------------------------------------------
# §3.3 Robustness variants for tests 1 and 3
# -----------------------------------------------------------------------------
# Each variant rebuilds the date means of the simple score and the
# polymer - simple difference, then reruns the season test.

variant_tests <- function(df, value, variant, simple = simple_compounds,
                          polymer = polymer_compounds, summed = FALSE) {
  d <- df |> filter(Substrate %in% c(simple, polymer))
  if (summed) {
    s <- d |>
      mutate(group = if_else(Substrate %in% simple, "simple", "polymer")) |>
      group_by(Sample_cat, Collection_date, season, group) |>
      summarise(v = sum(.data[[value]]), .groups = "drop") |>
      pivot_wider(names_from = group, values_from = v, names_glue = "{group}_score")
  } else {
    s <- d |>
      group_by(Substrate) |>
      mutate(z = (.data[[value]] - mean(.data[[value]])) / sd(.data[[value]])) |>
      ungroup() |>
      mutate(group = if_else(Substrate %in% simple, "simple", "polymer")) |>
      group_by(Sample_cat, Collection_date, season, group) |>
      summarise(v = mean(z), .groups = "drop") |>
      pivot_wider(names_from = group, values_from = v, names_glue = "{group}_score")
  }
  dm <- s |>
    group_by(Collection_date, season) |>
    summarise(simple_score = mean(simple_score), polymer_score = mean(polymer_score),
              .groups = "drop") |>
    mutate(diff_score = polymer_score - simple_score)
  bind_rows(
    season_test(dm$simple_score, dm$season) |> mutate(test = "1_simple_score"),
    season_test(dm$diff_score, dm$season) |> mutate(test = "3_polymer_minus_simple")
  ) |> mutate(variant = variant, n_simple = length(simple),
              n_polymer = length(polymer))
}

sir100_24 <- sir |> filter(Irrigation_trmt_pct == 100, Timepoint == 24)
flagged <- water_summary |>
  filter(flag_below_75pct, Substrate %in% c(simple_compounds, polymer_compounds)) |>
  pull(Substrate)

robustness <- bind_rows(
  variant_tests(sir100_5, "induced_C", "main (per C, water-subtracted, 5 h, z)"),
  variant_tests(sir100_5, "SIR_C", "a: uncorrected (no water subtraction)"),
  variant_tests(sir100_5, "induced_gdw", "b: per gram dry soil"),
  variant_tests(sir100_24, "induced_C", "c: 24 h readings"),
  variant_tests(sir100_5, "induced_C", "d: raw summed induced_C", summed = TRUE),
  map(c(simple_compounds, polymer_compounds), \(cmp) variant_tests(
    sir100_5, "induced_C", paste0("e: leave out ", cmp),
    simple = setdiff(simple_compounds, cmp), polymer = setdiff(polymer_compounds, cmp)
  )) |> list_rbind(),
  variant_tests(sir100_5, "induced_C",
                paste0("f: drop flagged (", paste(flagged, collapse = ", "), ")"),
                simple = setdiff(simple_compounds, flagged),
                polymer = setdiff(polymer_compounds, flagged)),
  variant_tests(sir100_5, "induced_C", "g: drop D.Cellobiose",
                simple = setdiff(simple_compounds, "D.Cellobiose"))
) |>
  relocate(test, variant) |>
  group_by(test) |>
  mutate(p_welch_fdr = p.adjust(p_welch, "BH"),
         p_wilcox_fdr = p.adjust(p_wilcox, "BH")) |>
  ungroup() |>
  arrange(test)
write_csv(robustness, file.path(dir_results, "cp2_robustness.csv"))

# -----------------------------------------------------------------------------
# §3.5 Exploratory: timing of polymer peaks
# -----------------------------------------------------------------------------

compound_dates <- read_csv(file.path(dir_results, "cp2_t5_compound_date_means.csv"),
                           show_col_types = FALSE)
peaks <- compound_dates |>
  filter(Substrate %in% c("Cellulose", "Hemicellulose", "Chitin", "Chitosan")) |>
  group_by(Substrate) |>
  slice_max(z, n = 1) |>
  ungroup() |>
  select(Substrate, peak_date = Collection_date, season, peak_z = z)
write_csv(peaks, file.path(dir_results, "cp2_explore_polymer_peaks.csv"))

# §3.5 Exploratory: date-level correlations with field variables
bgc <- read_csv(input_files$bgc, na = c("", "NA", "na"), show_col_types = FALSE) |>
  mutate(Collection_date = as.Date(Collection_date))
field_dates <- bgc |>
  filter(Irrigation_trmt_pct == 100, Collection_date %in% study_dates) |>
  group_by(Collection_date) |>
  summarise(SoilTemp_mean = mean(SoilTemp_mean), MBC_ug_gdw = mean(MBC_ug_gdw),
            SEOC_ug_gdw = mean(SEOC_ug_gdw), .groups = "drop")

field_cor <- expand_grid(score = c("simple_score", "polymer_score"),
                         variable = c("SoilTemp_mean", "MBC_ug_gdw", "SEOC_ug_gdw")) |>
  pmap(\(score, variable) {
    j <- inner_join(dates_C, field_dates, by = "Collection_date")
    cor_both(j[[score]], j[[variable]]) |> mutate(score = score, variable = variable)
  }) |>
  list_rbind() |>
  relocate(score, variable)
write_csv(field_cor, file.path(dir_results, "cp2_explore_field_correlations.csv"))

# -----------------------------------------------------------------------------
# §3.6 Replication in the 25% plots (tests 1 and 3, same pipeline)
# -----------------------------------------------------------------------------

sir25_5 <- sir |> filter(Irrigation_trmt_pct == 25, Timepoint == 5)
scores25 <- make_scores(sir25_5, "induced_C")
dates25 <- date_means(scores25)
write_csv(dates25, file.path(dir_results, "cp2_25pct_date_scores.csv"))

replication <- bind_rows(
  season_test(dates25$simple_score, dates25$season) |>
    mutate(test = "1_H1a_simple_score"),
  season_test(dates25$diff_score, dates25$season) |>
    mutate(test = "3_H1b_polymer_minus_simple")
) |> relocate(test) |> mutate(treatment = "25%", n_samples = nrow(scores25))
write_csv(replication, file.path(dir_results, "cp2_25pct_replication.csv"))

# -----------------------------------------------------------------------------
# Console summary
# -----------------------------------------------------------------------------

options(width = 200)
print(mixed)
print(mixed_var)
print(window_results |> select(-window_diffs))
print(robustness |> select(test, variant, diff, ci_low, ci_high, p_welch, p_wilcox,
                           p_welch_fdr), n = Inf)
print(peaks)
print(field_cor)
print(replication)
