# =============================================================================
# 02. Q1 hypothesis tests and interpretation steps (plan §3.1, §3.2)
# =============================================================================
#
# Tests 1-5 (hypothesis tests) and 6-8 (biomass normalization, enzyme
# support, later-to-early rate ratio). 100% plots, 5 h, water-subtracted
# per-C values unless stated. Date means are the unit of replication.

source(file.path("Scripts", "plan_v1", "00_setup.R"))

sir <- read_csv(file.path(dir_results, "sir_prepared_long.csv"),
                show_col_types = FALSE) |>
  mutate(season = factor(season, levels = c("Active", "Dormant")))
sir100_5 <- sir |> filter(Irrigation_trmt_pct == 100, Timepoint == 5)

sample_scores <- read_csv(file.path(dir_results, "sample_scores_100pct_5h.csv"),
                          show_col_types = FALSE)
date_scores <- read_csv(file.path(dir_results, "date_scores_100pct_5h.csv"),
                        show_col_types = FALSE)
roots <- read_csv(file.path(dir_results, "root_flags_active_100pct.csv"),
                  show_col_types = FALSE)

bgc <- read_csv(input_files$bgc, na = c("", "NA", "na"), show_col_types = FALSE) |>
  mutate(Collection_date = as.Date(Collection_date))
eea <- read_csv(input_files$eea, na = c("", "NA", "na", "#N/A"),
                show_col_types = FALSE) |>
  select(Sample_cat, BG_nmol_gdw_h:Ala_nmol_gdw_h) |>
  inner_join(bgc |> select(Sample_cat, Collection_date, Plot, Irrigation_trmt_pct),
             by = "Sample_cat") |>
  filter(Irrigation_trmt_pct == 100, Collection_date %in% study_dates) |>
  mutate(season = assign_season(Collection_date))

dates_C <- date_scores |> filter(normalization == "per_C")
dates_MBC <- date_scores |> filter(normalization == "per_MBC")

# -----------------------------------------------------------------------------
# Tests 1 and 3 (and 6: the same on biomass-normalized scores)
# -----------------------------------------------------------------------------

season_tests <- bind_rows(
  season_test(dates_C$simple_score, dates_C$season) |>
    mutate(test = "1_H1a_simple_score", normalization = "per_C"),
  season_test(dates_C$diff_score, dates_C$season) |>
    mutate(test = "3_H1b_polymer_minus_simple", normalization = "per_C"),
  season_test(dates_C$polymer_score, dates_C$season) |>
    mutate(test = "polymer_score (for forest plot)", normalization = "per_C"),
  season_test(dates_MBC$simple_score, dates_MBC$season) |>
    mutate(test = "6_H1a_simple_score", normalization = "per_MBC"),
  season_test(dates_MBC$diff_score, dates_MBC$season) |>
    mutate(test = "6_H1b_polymer_minus_simple", normalization = "per_MBC"),
  season_test(dates_MBC$polymer_score, dates_MBC$season) |>
    mutate(test = "6_polymer_score", normalization = "per_MBC")
) |>
  relocate(test, normalization)
write_csv(season_tests, file.path(dir_results, "cp2_t1_t3_t6_season_tests.csv"))

# -----------------------------------------------------------------------------
# Test 2. Simple score vs root biomass (active season, samples)
# -----------------------------------------------------------------------------

active_samples <- sample_scores |>
  filter(normalization == "per_C", season == "Active") |>
  inner_join(roots |> select(Sample_cat, Roots_gdw, root_flag), by = "Sample_cat")

spearman <- function(x, y) {
  s <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  tibble(n = sum(complete.cases(x, y)), rho = unname(s$estimate),
         p = s$p.value)
}

root_tests <- bind_rows(
  spearman(active_samples$simple_score, active_samples$Roots_gdw) |>
    mutate(set = "all active cores"),
  with(filter(active_samples, !root_flag),
       spearman(simple_score, Roots_gdw)) |>
    mutate(set = "flagged root cores removed")
) |>
  relocate(set)
write_csv(root_tests, file.path(dir_results, "cp2_t2_root_biomass_score.csv"))

root_compounds <- sir100_5 |>
  filter(Substrate %in% simple_compounds, season == "Active") |>
  inner_join(roots |> select(Sample_cat, Roots_gdw, root_flag), by = "Sample_cat") |>
  group_by(Substrate) |>
  group_modify(\(d, k) bind_rows(
    spearman(d$induced_C, d$Roots_gdw) |> mutate(set = "all active cores"),
    with(filter(d, !root_flag), spearman(induced_C, Roots_gdw)) |>
      mutate(set = "flagged root cores removed")
  )) |>
  ungroup() |>
  group_by(set) |>
  mutate(p_fdr = p.adjust(p, "BH")) |>
  ungroup()
write_csv(root_compounds, file.path(dir_results, "cp2_t2_root_biomass_compounds.csv"))

root_compound_counts <- root_compounds |>
  group_by(set) |>
  summarise(n_compounds = n(), n_positive = sum(rho > 0),
            n_positive_p05 = sum(rho > 0 & p < 0.05), .groups = "drop")
write_csv(root_compound_counts,
          file.path(dir_results, "cp2_t2_root_biomass_compound_counts.csv"))

# -----------------------------------------------------------------------------
# Test 4 (and 6). Coupling of polymer and simple date means
# -----------------------------------------------------------------------------

coupling <- bind_rows(
  cor_both(dates_C$polymer_score, dates_C$simple_score) |>
    mutate(normalization = "per_C", dates = "all 13"),
  with(filter(dates_C, season == "Dormant"), cor_both(polymer_score, simple_score)) |>
    mutate(normalization = "per_C", dates = "dormant 5"),
  cor_both(dates_MBC$polymer_score, dates_MBC$simple_score) |>
    mutate(normalization = "per_MBC", dates = "all 13"),
  with(filter(dates_MBC, season == "Dormant"), cor_both(polymer_score, simple_score)) |>
    mutate(normalization = "per_MBC", dates = "dormant 5")
) |>
  relocate(normalization, dates)
write_csv(coupling, file.path(dir_results, "cp2_t4_t6_coupling.csv"))

# -----------------------------------------------------------------------------
# Test 5. Compound-level season effects
# -----------------------------------------------------------------------------

forest_compounds <- c(simple_compounds, polymer_compounds, "D.Glucose")

compound_dates <- make_z_long(sir100_5, "induced_C", forest_compounds) |>
  group_by(Substrate, Collection_date, season) |>
  summarise(z = mean(z), induced_C = mean(induced_C), n_plots = n(),
            .groups = "drop") |>
  left_join(compound_groups, by = "Substrate")
write_csv(compound_dates, file.path(dir_results, "cp2_t5_compound_date_means.csv"))

compound_effects <- compound_dates |>
  group_by(Substrate, group) |>
  group_modify(\(d, k) season_test(d$z, d$season) |>
                 mutate(ratio_induced_C = mean(d$induced_C[d$season == "Dormant"]) /
                          mean(d$induced_C[d$season == "Active"]))) |>
  ungroup() |>
  mutate(p_welch_fdr = p.adjust(p_welch, "BH"),
         p_wilcox_fdr = p.adjust(p_wilcox, "BH")) |>
  arrange(factor(group, c("Simple", "Polymer", "Reference")), desc(diff))
write_csv(compound_effects, file.path(dir_results, "cp2_t5_compound_effects.csv"))

# -----------------------------------------------------------------------------
# Test 7. Enzyme support for H1b
# -----------------------------------------------------------------------------

enzyme_dates <- eea |>
  group_by(Collection_date, season) |>
  summarise(n_plots = n(),
            across(c(CB_nmol_gdw_h, BG_nmol_gdw_h, BX_nmol_gdw_h, NAG_nmol_gdw_h),
                   mean),
            .groups = "drop") |>
  pivot_longer(ends_with("_nmol_gdw_h"), names_to = "enzyme",
               values_to = "activity") |>
  mutate(enzyme = str_remove(enzyme, "_nmol_gdw_h"))
write_csv(enzyme_dates, file.path(dir_results, "cp2_t7_enzyme_date_means.csv"))

enzyme_season <- enzyme_dates |>
  group_by(enzyme) |>
  group_modify(\(d, k) season_test(d$activity, d$season) |>
                 mutate(ratio = mean_dormant / mean_active)) |>
  ungroup() |>
  mutate(p_welch_fdr = p.adjust(p_welch, "BH"),
         p_wilcox_fdr = p.adjust(p_wilcox, "BH"))
write_csv(enzyme_season, file.path(dir_results, "cp2_t7_enzyme_season.csv"))

enzyme_pairs <- tribble(
  ~enzyme, ~Substrate,
  "CB", "Cellulose", "BG", "Cellulose", "BX", "Hemicellulose", "NAG", "Chitin"
)

enzyme_substrate <- enzyme_pairs |>
  pmap(\(enzyme, Substrate) {
    e <- enzyme_dates[enzyme_dates$enzyme == enzyme, ]
    s <- compound_dates[compound_dates$Substrate == Substrate, ]
    j <- inner_join(e, s, by = "Collection_date")
    spearman(j$activity, j$z) |> mutate(enzyme = enzyme, Substrate = Substrate)
  }) |>
  list_rbind() |>
  relocate(enzyme, Substrate) |>
  mutate(p_fdr = p.adjust(p, "BH"))
write_csv(enzyme_substrate, file.path(dir_results, "cp2_t7_enzyme_substrate.csv"))

# -----------------------------------------------------------------------------
# Test 8. Later (24 h) to early (5 h) rate ratio (see D6, O1)
# -----------------------------------------------------------------------------
# Per date and group: mean induced_C at 24 h over all plots and group
# compounds, divided by the same mean at 5 h.

rate_ratio_dates <- sir |>
  filter(Irrigation_trmt_pct == 100, group %in% c("Simple", "Polymer")) |>
  group_by(Collection_date, season, group, Timepoint) |>
  summarise(mean_induced_C = mean(induced_C), n = n(), .groups = "drop") |>
  pivot_wider(names_from = Timepoint, values_from = c(mean_induced_C, n),
              names_prefix = "h") |>
  mutate(ratio_24_5 = mean_induced_C_h24 / mean_induced_C_h5)
write_csv(rate_ratio_dates, file.path(dir_results, "cp2_t8_rate_ratio_dates.csv"))

rate_ratio_tests <- rate_ratio_dates |>
  group_by(group) |>
  group_modify(\(d, k) season_test(d$ratio_24_5, d$season)) |>
  ungroup()
write_csv(rate_ratio_tests, file.path(dir_results, "cp2_t8_rate_ratio_season.csv"))

# -----------------------------------------------------------------------------
# Console summary
# -----------------------------------------------------------------------------

options(width = 200)
print(season_tests)
print(root_tests)
print(root_compound_counts)
print(coupling)
print(compound_effects |> select(Substrate, group, diff, ci_low, ci_high, p_welch,
                                 p_wilcox, ratio_induced_C), n = Inf)
print(enzyme_season)
print(enzyme_substrate)
print(rate_ratio_tests)
