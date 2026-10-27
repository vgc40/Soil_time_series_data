# =============================================================================
# 06. Q2 analyses (plan §5, Checkpoint 4)
# =============================================================================
#
# Expression values from the approved gene sets (D8), potential values
# (per-C SIR as main, D11; per-MBC SIR as robustness; enzyme activities), the bridge to Q1,
# H2a, H2b, H2c (lag and pathway step), plot-level lag models, and the two
# robustness variants. Expression is geTMM summed per KO (O2).

source(file.path("Scripts", "plan_v1", "00_setup.R"))
suppressPackageStartupMessages({
  library(lme4)
})

set.seed(2026)

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------

gene_sets_all <- read_csv(file.path(dir_results, "KO_gene_sets_RECOMMENDED.csv"),
                          show_col_types = FALSE)
write_csv(gene_sets_all, file.path(dir_results, "KO_gene_sets_APPROVED.csv"))
gene_sets <- gene_sets_all |> filter(status == "use in Q2")
set_info <- gene_sets |> distinct(set, set_type, group, pathway_step)

bgc <- read_csv(input_files$bgc, na = c("", "NA", "na"), show_col_types = FALSE) |>
  mutate(Collection_date = as.Date(Collection_date))
samples <- bgc |>
  filter(Irrigation_trmt_pct == 100, Collection_date %in% study_dates) |>
  select(Sample_cat, Collection_date, Plot) |>
  mutate(season = assign_season(Collection_date))

sir <- read_csv(file.path(dir_results, "sir_prepared_long.csv"), show_col_types = FALSE) |>
  mutate(season = factor(season, levels = c("Active", "Dormant")))
sir100_5 <- sir |> filter(Irrigation_trmt_pct == 100, Timepoint == 5)

eea <- read_csv(input_files$eea, na = c("", "NA", "na", "#N/A"), show_col_types = FALSE) |>
  select(Sample_cat, BG_nmol_gdw_h:Ala_nmol_gdw_h) |>
  inner_join(samples, by = "Sample_cat")

# -----------------------------------------------------------------------------
# §5.1 Expression values: zero replacement, set sums, centred log ratio
# -----------------------------------------------------------------------------
# Per sample: zeros -> half the smallest nonzero value in that sample. Set sum
# over the replaced values; CLR = log(set sum) - mean log expression over all
# KOs in that sample. Rows without a KO identifier are not KOs and are left
# out of the per-sample mean.

expr <- read_csv(input_files$brite_d, show_col_types = FALSE)
names(expr)[1] <- "brite_D"
expr <- expr |>
  mutate(KO = str_extract(brite_D, "^K\\d{5}")) |>
  filter(!is.na(KO))
sample_cols <- str_subset(names(expr), "_MT$")

mat <- as.matrix(expr[, sample_cols])
rownames(mat) <- expr$KO
half_min <- apply(mat, 2, \(x) min(x[x > 0]) / 2)
mat_rep <- sweep(mat, 2, half_min, \(x, h) ifelse(x == 0, h, x))
mean_log <- colMeans(log(mat_rep))

zero_replacement <- tibble(sample = sample_cols, half_min = half_min,
                           n_zero = colSums(mat == 0), n_KOs = nrow(mat),
                           mean_log = mean_log)
write_csv(zero_replacement, file.path(dir_results, "cp4_zero_replacement.csv"))

clr_for_sets <- function(sets) {
  sets |>
    distinct(set, KO) |>
    group_by(set) |>
    group_modify(\(d, k) {
      s <- colSums(mat_rep[d$KO, , drop = FALSE])
      tibble(sample = names(s), set_sum = s, clr = log(s) - mean_log[names(s)])
    }) |>
    ungroup() |>
    mutate(Sample_cat = str_remove(sample, "_MT$")) |>
    inner_join(samples, by = "Sample_cat")
}

expr_clr <- clr_for_sets(gene_sets)
write_csv(expr_clr, file.path(dir_results, "cp4_expression_clr_samples.csv"))

transcript_dates <- sort(unique(expr_clr$Collection_date))
transcript_plots <- expr_clr |> distinct(Collection_date, Plot)

# -----------------------------------------------------------------------------
# §5.2 Potential values (sample level)
# -----------------------------------------------------------------------------
# Compounds: water-subtracted per-C SIR (main, D11) and per-MBC SIR (§2.7,
# robustness a). Enzymes: measured activity.

compound_sets <- set_info |> filter(set_type == "compound") |> pull(set)
enzyme_sets <- set_info |> filter(set_type == "enzyme") |> pull(set)

potential_samples <- function(sir_value) {
  bind_rows(
    sir100_5 |>
      filter(Substrate %in% compound_sets) |>
      transmute(set = Substrate, Sample_cat, Collection_date, Plot, season,
                potential = .data[[sir_value]], MBC_flag),
    eea |>
      pivot_longer(BG_nmol_gdw_h:Ala_nmol_gdw_h, names_to = "set",
                   values_to = "potential") |>
      mutate(set = str_remove(set, "_nmol_gdw_h"), MBC_flag = "none") |>
      filter(set %in% enzyme_sets) |>
      select(set, Sample_cat, Collection_date, Plot, season, potential, MBC_flag)
  )
}

# Main potential is per mg C (D11); the plan's per-MBC version is robustness (a)
pot_main <- potential_samples("induced_C")
pot_mbc <- potential_samples("induced_C_per_MBC")
write_csv(pot_main, file.path(dir_results, "cp4_potential_samples.csv"))
write_csv(pot_mbc, file.path(dir_results, "cp4_potential_samples_per_MBC.csv"))

# -----------------------------------------------------------------------------
# §5.3.1 Bridge to Q1: tests 1 and 3 on the 7 transcript dates
# -----------------------------------------------------------------------------
# Uses the Q1 date-mean scores (z-scores from all 77 samples), restricted to
# the transcript dates.

date_scores <- read_csv(file.path(dir_results, "date_scores_100pct_5h.csv"),
                        show_col_types = FALSE) |>
  filter(normalization == "per_C")
q1_tests <- read_csv(file.path(dir_results, "cp2_t1_t3_t6_season_tests.csv"),
                     show_col_types = FALSE) |>
  filter(normalization == "per_C", str_detect(test, "^[13]_"))

bridge_dates <- date_scores |> filter(Collection_date %in% transcript_dates)
bridge <- bind_rows(
  season_test(bridge_dates$simple_score, bridge_dates$season) |>
    mutate(test = "1_H1a_simple_score"),
  season_test(bridge_dates$diff_score, bridge_dates$season) |>
    mutate(test = "3_H1b_polymer_minus_simple")
) |>
  mutate(dates = "7 transcript dates") |>
  left_join(q1_tests |> select(test, diff_13 = diff, p_welch_13 = p_welch),
            by = "test") |>
  mutate(same_direction = sign(diff) == sign(diff_13)) |>
  relocate(test, dates)
write_csv(bridge, file.path(dir_results, "cp4_bridge_q1.csv"))

# -----------------------------------------------------------------------------
# §5.3.2 H2a: simple and polymer expression scores (descriptive)
# -----------------------------------------------------------------------------

zscore <- function(x) (x - mean(x)) / sd(x)

welch_ci <- function(value, season) {
  d <- value[season == "Dormant"]
  a <- value[season == "Active"]
  w <- t.test(d, a)
  tibble(n_dormant = length(d), n_active = length(a), mean_dormant = mean(d),
         mean_active = mean(a), diff = mean(d) - mean(a),
         ci_low = w$conf.int[1], ci_high = w$conf.int[2])
}

h2a_samples <- expr_clr |>
  inner_join(set_info, by = "set") |>
  filter(group %in% c("Simple", "Polymer")) |>
  group_by(set) |>
  mutate(z = zscore(clr)) |>
  ungroup() |>
  group_by(Sample_cat, Collection_date, Plot, season, group) |>
  summarise(score = mean(z), n_sets = n(), .groups = "drop")

h2a_dates <- h2a_samples |>
  group_by(Collection_date, season, group) |>
  summarise(score = mean(score), n_plots = n(), .groups = "drop")
write_csv(h2a_dates, file.path(dir_results, "cp4_h2a_expression_scores_dates.csv"))

h2a <- h2a_dates |>
  pivot_wider(names_from = group, values_from = score) |>
  mutate(Polymer_minus_Simple = Polymer - Simple) |>
  pivot_longer(c(Simple, Polymer, Polymer_minus_Simple), names_to = "score",
               values_to = "value") |>
  group_by(score) |>
  group_modify(\(d, k) welch_ci(d$value, d$season)) |>
  ungroup() |>
  mutate(n_sets = case_when(score == "Simple" ~ sum(set_info$group == "Simple"),
                            score == "Polymer" ~ sum(set_info$group == "Polymer"),
                            .default = NA_integer_))
write_csv(h2a, file.path(dir_results, "cp4_h2a_expression_season.csv"))

# -----------------------------------------------------------------------------
# §5.3.3 H2b: seasonal concordance (compounds, 7 transcript dates)
# -----------------------------------------------------------------------------
# Season effect in z-units: z-score across the samples on the transcript
# dates, average to date means, dormant - active.

season_effect_z <- function(df, value) {
  df |>
    group_by(set) |>
    mutate(z = zscore(.data[[value]])) |>
    group_by(set, Collection_date, season) |>
    summarise(z = mean(z), n_plots = n(), .groups = "drop") |>
    group_by(set) |>
    group_modify(\(d, k) welch_ci(d$z, d$season)) |>
    ungroup()
}

h2b_effects <- function(pot, expr_values) {
  e <- season_effect_z(expr_values |> filter(set %in% compound_sets), "clr")
  p <- season_effect_z(pot |> filter(set %in% compound_sets,
                                     Collection_date %in% transcript_dates),
                       "potential")
  inner_join(e |> select(set, expr_effect = diff, expr_ci_low = ci_low,
                         expr_ci_high = ci_high),
             p |> select(set, pot_effect = diff, pot_ci_low = ci_low,
                         pot_ci_high = ci_high),
             by = "set") |>
    left_join(set_info, by = "set") |>
    mutate(concordant = sign(expr_effect) == sign(pot_effect))
}

h2b_summary <- function(eff, variant) {
  b <- binom.test(sum(eff$concordant), nrow(eff), p = 0.5)
  s <- suppressWarnings(cor.test(eff$pot_effect, eff$expr_effect,
                                 method = "spearman", exact = TRUE))
  tibble(variant = variant, n_compounds = nrow(eff),
         n_concordant = sum(eff$concordant), share_concordant = mean(eff$concordant),
         binom_p = b$p.value, spearman_rho = unname(s$estimate),
         spearman_p = s$p.value)
}

h2b <- h2b_effects(pot_main, expr_clr)
write_csv(h2b, file.path(dir_results, "cp4_h2b_compound_effects.csv"))
h2b_main <- h2b_summary(h2b, "main (potential per C; D11)")

# -----------------------------------------------------------------------------
# §5.3.4 H2c: lag (date means, plot-matched)
# -----------------------------------------------------------------------------

lag_pairs <- tibble(
  t = as.Date(c("2023-06-07", "2023-07-05", "2023-09-06", "2023-10-04",
                "2024-01-02", "2024-03-06", "2024-05-08")),
  t_next = as.Date(c("2023-07-05", "2023-08-02", "2023-10-04", "2023-11-01",
                     "2024-02-12", "2024-04-03", "2024-06-05"))
)
stopifnot(all(lag_pairs$t %in% transcript_dates))

# Plot-level table: expression at t, potential at t and at t_next, same plot
plot_pairs <- function(pot, expr_values) {
  expr_values |>
    select(set, Plot, t = Collection_date, clr) |>
    inner_join(lag_pairs, by = "t") |>
    left_join(pot |> select(set, Plot, t = Collection_date, pot_same = potential,
                            flag_same = MBC_flag),
              by = c("set", "Plot", "t")) |>
    left_join(pot |> select(set, Plot, t_next = Collection_date,
                            pot_next = potential, flag_next = MBC_flag),
              by = c("set", "Plot", "t_next"))
}

lag_date_means <- function(pp) {
  bind_rows(
    pp |> filter(!is.na(pot_same)) |>
      group_by(set, t, t_next) |>
      summarise(expr = mean(clr), pot = mean(pot_same), n_plots = n(),
                .groups = "drop") |> mutate(comparison = "same"),
    pp |> filter(!is.na(pot_next)) |>
      group_by(set, t, t_next) |>
      summarise(expr = mean(clr), pot = mean(pot_next), n_plots = n(),
                .groups = "drop") |> mutate(comparison = "lag")
  )
}

lag_stats <- function(dm, n_boot = 2000) {
  set.seed(2026)  # same draws for every variant, so results do not depend on run order
  dm |>
    group_by(set) |>
    group_modify(\(d, k) {
      same <- d |> filter(comparison == "same") |> arrange(t)
      lag <- d |> filter(comparison == "lag") |> arrange(t)
      j <- inner_join(same |> select(t, e_same = expr, p_same = pot),
                      lag |> select(t, e_lag = expr, p_lag = pot), by = "t")
      r_same <- cor(j$e_same, j$p_same)
      r_lag <- cor(j$e_lag, j$p_lag)
      boot <- replicate(n_boot, {
        i <- sample(nrow(j), replace = TRUE)
        suppressWarnings(cor(j$e_lag[i], j$p_lag[i]) - cor(j$e_same[i], j$p_same[i]))
      })
      boot <- boot[is.finite(boot)]
      tibble(n_pairs = nrow(j), r_same = r_same, r_lag = r_lag,
             lag_advantage = r_lag - r_same,
             ci_low = unname(quantile(boot, 0.025)),
             ci_high = unname(quantile(boot, 0.975)),
             n_boot_valid = length(boot))
    }) |>
    ungroup() |>
    left_join(set_info, by = "set")
}

pp_main <- plot_pairs(pot_main, expr_clr)
write_csv(pp_main, file.path(dir_results, "cp4_lag_plot_pairs.csv"))
dm_main <- lag_date_means(pp_main)
write_csv(dm_main, file.path(dir_results, "cp4_lag_date_means.csv"))
lag_main <- lag_stats(dm_main)
write_csv(lag_main, file.path(dir_results, "cp4_h2c_lag_advantage.csv"))

# -----------------------------------------------------------------------------
# §5.3.5 H2c: extracellular vs intracellular lag advantage (permutation)
# -----------------------------------------------------------------------------

step_permutation <- function(lag, variant, n_perm = 10000) {
  set.seed(2026)
  x <- lag$lag_advantage
  lab <- lag$pathway_step
  stat <- function(l) mean(x[l == "extracellular"]) - mean(x[l == "intracellular"])
  obs <- stat(lab)
  perm <- replicate(n_perm, stat(sample(lab)))
  tibble(variant = variant,
         n_extracellular = sum(lab == "extracellular"),
         n_intracellular = sum(lab == "intracellular"),
         mean_extracellular = mean(x[lab == "extracellular"]),
         mean_intracellular = mean(x[lab == "intracellular"]),
         diff = obs, p_perm = mean(abs(perm) >= abs(obs)), n_perm = n_perm)
}

step_main <- step_permutation(lag_main, "main (potential per C; D11)")

# -----------------------------------------------------------------------------
# §5.4 Plot-level lag model
# -----------------------------------------------------------------------------
# For each gene set: potential(t+1) ~ expression(t) + (1|date) and
# potential(t) ~ expression(t) + (1|date), plot-matched pairs. Response and
# predictor z-scored within each model's data so the fits are on one scale.
# ML fits for AIC. Marginal R2 = var(fixed) / (var(fixed) + var(date) + var(resid)).

fit_lag_model <- function(d) {
  d <- d |> mutate(y = zscore(y), x = zscore(x), date = factor(t))
  fit <- lmer(y ~ x + (1 | date), data = d, REML = FALSE)
  vc <- as.data.frame(VarCorr(fit))
  var_f <- var(as.vector(model.matrix(fit) %*% fixef(fit)))
  var_date <- vc$vcov[vc$grp == "date"]
  var_res <- vc$vcov[vc$grp == "Residual"]
  tibble(n = nobs(fit), slope = unname(fixef(fit)["x"]),
         marginal_R2 = var_f / (var_f + var_date + var_res),
         AIC = AIC(fit), singular = isSingular(fit))
}

# Both models use the same rows (plot-dates with potential at t and t+1), so
# their AICs are comparable (D9).
lag_models <- pp_main |>
  filter(!is.na(pot_same), !is.na(pot_next)) |>
  group_by(set) |>
  group_modify(\(d, k) bind_rows(
    fit_lag_model(d |> filter(!is.na(pot_next)) |> transmute(t, x = clr, y = pot_next)) |>
      mutate(model = "lag: potential(t+1) ~ expression(t)"),
    fit_lag_model(d |> filter(!is.na(pot_same)) |> transmute(t, x = clr, y = pot_same)) |>
      mutate(model = "same: potential(t) ~ expression(t)")
  )) |>
  ungroup() |>
  left_join(set_info, by = "set")
write_csv(lag_models, file.path(dir_results, "cp4_lag_models_plot_level.csv"))

lag_model_compare <- lag_models |>
  mutate(model = if_else(str_starts(model, "lag"), "lag", "same")) |>
  select(set, pathway_step, model, n, marginal_R2, AIC, singular) |>
  pivot_wider(names_from = model, values_from = c(n, marginal_R2, AIC, singular)) |>
  mutate(R2_lag_minus_same = marginal_R2_lag - marginal_R2_same,
         AIC_lag_minus_same = AIC_lag - AIC_same)
write_csv(lag_model_compare, file.path(dir_results, "cp4_lag_models_compare.csv"))

# -----------------------------------------------------------------------------
# §5.4 Robustness: (a) potential without biomass normalization;
# (b) high-confidence KOs only
# -----------------------------------------------------------------------------

# (a) Potential per MBC, the plan's original main version (D11)
h2b_mbc <- h2b_effects(pot_mbc, expr_clr)
lag_mbc <- lag_stats(lag_date_means(plot_pairs(pot_mbc, expr_clr)))

# (b) High-confidence KOs only; sets that become near-zero (sum 0 in 9 or more
# samples, D8 rule 4) are left out.
high_sets <- gene_sets |> filter(confidence == "high")
high_zero <- high_sets |>
  distinct(set, KO) |>
  group_by(set) |>
  summarise(n_zero = sum(colSums(mat[KO, , drop = FALSE]) == 0), .groups = "drop")
high_keep <- high_zero |> filter(n_zero < ceiling(length(sample_cols) / 4)) |> pull(set)
expr_clr_high <- clr_for_sets(high_sets |> filter(set %in% high_keep))

h2b_high <- h2b_effects(pot_main, expr_clr_high)
lag_high <- lag_stats(lag_date_means(plot_pairs(pot_main, expr_clr_high)))

robust_h2b <- bind_rows(
  h2b_main,
  h2b_summary(h2b_mbc, "a: potential per MBC (plan original)"),
  h2b_summary(h2b_high, "b: high-confidence KOs only")
)
robust_step <- bind_rows(
  step_main,
  step_permutation(lag_mbc, "a: potential per MBC (plan original)"),
  step_permutation(lag_high, "b: high-confidence KOs only")
)
write_csv(robust_h2b, file.path(dir_results, "cp4_h2b_summary_and_robustness.csv"))
write_csv(robust_step, file.path(dir_results, "cp4_h2c_step_permutation_and_robustness.csv"))
write_csv(bind_rows(lag_main |> mutate(variant = "main"),
                    lag_mbc |> mutate(variant = "a: potential per MBC"),
                    lag_high |> mutate(variant = "b: high-confidence KOs")),
          file.path(dir_results, "cp4_h2c_lag_advantage_all_variants.csv"))
write_csv(bind_rows(h2b |> mutate(variant = "main"),
                    h2b_mbc |> mutate(variant = "a: potential per MBC"),
                    h2b_high |> mutate(variant = "b: high-confidence KOs")),
          file.path(dir_results, "cp4_h2b_compound_effects_all_variants.csv"))
write_csv(high_zero |> mutate(kept = set %in% high_keep),
          file.path(dir_results, "cp4_high_confidence_sets_kept.csv"))

# -----------------------------------------------------------------------------
# Console summary
# -----------------------------------------------------------------------------

options(width = 200)
cat("Transcript dates:", format(transcript_dates), "\n")
print(summary(zero_replacement$half_min)); print(summary(zero_replacement$n_zero))
print(bridge)
print(h2a)
print(h2b |> select(set, pathway_step, expr_effect, pot_effect, concordant), n = Inf)
print(robust_h2b)
print(lag_main |> select(set, pathway_step, n_pairs, r_same, r_lag, lag_advantage,
                         ci_low, ci_high) |> arrange(pathway_step, desc(lag_advantage)),
      n = Inf)
print(robust_step)
print(lag_model_compare |> select(set, pathway_step, n_lag, n_same, marginal_R2_lag,
                                  marginal_R2_same, AIC_lag_minus_same), n = Inf)
print(high_zero |> filter(!set %in% high_keep))
