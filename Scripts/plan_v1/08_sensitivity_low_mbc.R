# =============================================================================
# 08. Sensitivity: biomass-normalized analyses without the low-MBC samples
# =============================================================================
#
# Requested by VGC after Checkpoint 4 (deviation log D10). Reruns the
# biomass-normalized analyses without SM1625 (MBC 0 replaced with 0.07) and
# SM1539 (MBC 0.14):
#   Q1 test 6 (tests 1, 3, 4 on per-MBC scores, 13 dates)
#   Q2 H2b (seasonal concordance) and H2c (lag advantage, pathway step)
# Per-MBC z-scores are recomputed without the two samples. Run after 06.

source(file.path("Scripts", "plan_v1", "06_q2_analyses.R"))
set.seed(2026)

low_mbc <- c("SM1625", "SM1539")

# -----------------------------------------------------------------------------
# Q1 test 6 without the two samples
# -----------------------------------------------------------------------------

sir100_5_sens <- sir100_5 |> filter(!Sample_cat %in% low_mbc)
scores_MBC_sens <- make_scores(sir100_5_sens, "induced_C_per_MBC")
dates_MBC_sens <- date_means(scores_MBC_sens)
write_csv(dates_MBC_sens, file.path(dir_results, "sens_q1_date_scores_per_MBC.csv"))

q1_sens <- bind_rows(
  season_test(dates_MBC_sens$simple_score, dates_MBC_sens$season) |>
    mutate(test = "6_H1a_simple_score"),
  season_test(dates_MBC_sens$diff_score, dates_MBC_sens$season) |>
    mutate(test = "6_H1b_polymer_minus_simple"),
  season_test(dates_MBC_sens$polymer_score, dates_MBC_sens$season) |>
    mutate(test = "6_polymer_score")
) |>
  mutate(variant = "without SM1625, SM1539", n_samples = nrow(scores_MBC_sens)) |>
  relocate(test, variant)

q1_orig <- read_csv(file.path(dir_results, "cp2_t1_t3_t6_season_tests.csv"),
                    show_col_types = FALSE) |>
  filter(normalization == "per_MBC") |>
  mutate(variant = "all samples (Checkpoint 2)", n_samples = 77L) |>
  select(-normalization)
write_csv(bind_rows(q1_orig, q1_sens) |> relocate(test, variant) |> arrange(test),
          file.path(dir_results, "sens_q1_test6_season.csv"))

coupling_sens <- bind_rows(
  cor_both(dates_MBC_sens$polymer_score, dates_MBC_sens$simple_score) |>
    mutate(dates = "all 13"),
  with(filter(dates_MBC_sens, season == "Dormant"), cor_both(polymer_score, simple_score)) |>
    mutate(dates = "dormant 5")
) |> mutate(variant = "without SM1625, SM1539") |> relocate(variant, dates)
write_csv(coupling_sens, file.path(dir_results, "sens_q1_test6_coupling.csv"))

# -----------------------------------------------------------------------------
# Q2 H2b and H2c without the two samples
# -----------------------------------------------------------------------------
# Enzyme potential does not use MBC, so enzyme sets are unchanged; they stay
# in the pathway-step comparison.

pot_sens <- pot_main |> filter(!Sample_cat %in% low_mbc)

h2b_sens <- h2b_effects(pot_sens, expr_clr)
write_csv(h2b_sens, file.path(dir_results, "sens_q2_h2b_compound_effects.csv"))
h2b_sens_summary <- bind_rows(
  h2b_main,
  h2b_summary(h2b_sens, "without SM1625, SM1539 (potential per MBC)")
)
write_csv(h2b_sens_summary, file.path(dir_results, "sens_q2_h2b_summary.csv"))

lag_sens <- lag_stats(lag_date_means(plot_pairs(pot_sens, expr_clr)))
write_csv(lag_sens, file.path(dir_results, "sens_q2_h2c_lag_advantage.csv"))
step_sens <- bind_rows(
  step_main,
  step_permutation(lag_sens, "without SM1625, SM1539 (potential per MBC)")
)
write_csv(step_sens, file.path(dir_results, "sens_q2_h2c_step_permutation.csv"))

# -----------------------------------------------------------------------------
# Console summary
# -----------------------------------------------------------------------------

options(width = 200)
cat("\n===== SENSITIVITY =====\n")
print(bind_rows(q1_orig, q1_sens) |> select(test, variant, n_samples, diff, ci_low,
                                            ci_high, p_welch, p_wilcox) |> arrange(test))
print(dates_MBC_sens |> select(Collection_date, season, n_plots, simple_score,
                               polymer_score))
print(coupling_sens)
print(h2b_sens |> select(set, pathway_step, expr_effect, pot_effect, concordant), n = Inf)
print(h2b_sens_summary)
print(lag_sens |> select(set, pathway_step, r_same, r_lag, lag_advantage, ci_low, ci_high) |>
        arrange(pathway_step, desc(lag_advantage)), n = Inf)
print(step_sens)
