# =============================================================================
# 01. Shared data preparation (plan §2, Checkpoint 1)
# =============================================================================
#
# Reads the DataHub March tables, subtracts the water control, converts to per
# mg added C, assigns season, builds the simple and polymer scores (per C and
# per MBC), and flags low above-water compounds and root-dense cores.
# Writes all prepared tables and the Checkpoint 1 summaries to Results/plan_v1.

source(file.path("Scripts", "plan_v1", "00_setup.R"))

# -----------------------------------------------------------------------------
# Read inputs
# -----------------------------------------------------------------------------

# MicroResp is not UTF-8 (plan §1). Water rows have SIR = "na" in the per-C
# normalization because no C is added.
microresp <- read_csv(input_files$microresp,
                      locale = locale(encoding = "latin1"),
                      na = c("", "NA", "na"), show_col_types = FALSE) |>
  mutate(Collection_date = mdy(Collection_date))

bgc <- read_csv(input_files$bgc, na = c("", "NA", "na"),
                show_col_types = FALSE) |>
  mutate(Collection_date = as.Date(Collection_date))

# EEA dates are M/D/YY and one row (SM1502) has "#N/A". Dates are taken from
# the biogeochemistry table through Sample_cat, so the EEA date is not used.
eea <- read_csv(input_files$eea, na = c("", "NA", "na", "#N/A"),
                show_col_types = FALSE) |>
  select(Sample_cat, BG_nmol_gdw_h:Ala_nmol_gdw_h, MBC_ug_gdw, MBN_ug_gdw)

transcript_samples <- read_csv(input_files$brite_d, n_max = 0,
                               show_col_types = FALSE) |>
  names() |>
  str_subset("_MT$") |>
  str_remove("_MT$")

# Sample table: one row per plot and date on the 13 study dates
samples <- bgc |>
  filter(Collection_date %in% study_dates) |>
  select(Sample_cat, Collection_date, Plot, Irrigation_trmt_pct) |>
  mutate(season = assign_season(Collection_date))

# -----------------------------------------------------------------------------
# §2.2 Per-C conversion check
# -----------------------------------------------------------------------------
# Recompute the file's per-C values from the per-gram values:
# SIR_C = SIR_gdw * gdw_well / mgC_per_well

sir_wide <- microresp |>
  select(Sample_cat, Collection_date, Plot, Irrigation_trmt_pct, Timepoint,
         Substrate, gdw_well, mgC_per_well, Normalization, SIR) |>
  pivot_wider(names_from = Normalization, values_from = SIR) |>
  rename(SIR_gdw = `ngCO2-C_gdw_h`, SIR_C_file = `ngCO2-C_mgC_h`)

percarbon_check <- sir_wide |>
  filter(Substrate != "Water") |>
  mutate(SIR_C_recomputed = SIR_gdw * gdw_well / mgC_per_well,
         abs_diff = abs(SIR_C_recomputed - SIR_C_file),
         rel_diff = abs_diff / abs(SIR_C_file))

percarbon_summary <- percarbon_check |>
  summarise(n_values = n(),
            max_abs_diff = max(abs_diff, na.rm = TRUE),
            median_abs_diff = median(abs_diff, na.rm = TRUE),
            max_rel_diff = max(rel_diff[SIR_C_file != 0], na.rm = TRUE),
            share_within_0.01 = mean(abs_diff <= 0.01, na.rm = TRUE),
            share_within_1pct = mean(rel_diff <= 0.01, na.rm = TRUE),
            correlation = cor(SIR_C_recomputed, SIR_C_file, use = "complete.obs"))
write_csv(percarbon_summary, file.path(dir_results, "cp1_percarbon_check.csv"))

# Per substrate: the mg C per well implied by the file's per-C values
# (open item O3; interim decision D5 uses the mgC_per_well column)
percarbon_by_substrate <- percarbon_check |>
  mutate(implied_mgC = SIR_gdw * gdw_well / SIR_C_file) |>
  group_by(Substrate) |>
  summarise(mgC_per_well = first(mgC_per_well),
            implied_mgC_median = median(implied_mgC, na.rm = TRUE),
            implied_mgC_min = min(implied_mgC, na.rm = TRUE),
            implied_mgC_max = max(implied_mgC, na.rm = TRUE),
            median_rel_diff = median(rel_diff, na.rm = TRUE),
            max_rel_diff = max(rel_diff, na.rm = TRUE)) |>
  mutate(match = max_rel_diff <= 0.05) |>
  arrange(match, Substrate)
write_csv(percarbon_by_substrate,
          file.path(dir_results, "cp1_percarbon_check_by_substrate.csv"))
write_csv(percarbon_check |> slice_max(rel_diff, n = 20),
          file.path(dir_results, "cp1_percarbon_check_largest.csv"))

# -----------------------------------------------------------------------------
# §2.2 Water subtraction (both timepoints, both treatments; negatives kept)
# -----------------------------------------------------------------------------

water <- sir_wide |>
  filter(Substrate == "Water") |>
  select(Sample_cat, Timepoint, water_gdw = SIR_gdw)

sir <- sir_wide |>
  filter(Substrate != "Water") |>
  left_join(water, by = c("Sample_cat", "Timepoint")) |>
  mutate(induced_gdw = SIR_gdw - water_gdw,
         induced_C = induced_gdw * gdw_well / mgC_per_well,
         season = assign_season(Collection_date)) |>
  left_join(compound_groups, by = "Substrate") |>
  mutate(group = replace_na(group, "Excluded"))

# -----------------------------------------------------------------------------
# §2.7 MBC, with zeros replaced by half the lowest nonzero value (D1)
# -----------------------------------------------------------------------------

mbc_half_min <- min(eea$MBC_ug_gdw[eea$MBC_ug_gdw > 0], na.rm = TRUE) / 2

mbc <- eea |>
  select(Sample_cat, MBC_ug_gdw) |>
  mutate(
    MBC_flag = case_when(
      MBC_ug_gdw == 0 ~ "replaced_zero",
      MBC_ug_gdw == mbc_half_min * 2 ~ "lowest_measured",
      .default = "none"
    ),
    MBC_used = if_else(MBC_ug_gdw == 0, mbc_half_min, MBC_ug_gdw)
  )

mbc_replacements <- mbc |>
  filter(MBC_flag != "none") |>
  left_join(samples, by = "Sample_cat") |>
  arrange(Irrigation_trmt_pct, Collection_date)
write_csv(mbc_replacements, file.path(dir_results, "cp1_mbc_replacements.csv"))

sir <- sir |>
  left_join(mbc, by = "Sample_cat") |>
  mutate(induced_C_per_MBC = induced_C / MBC_used)

write_csv(sir, file.path(dir_results, "sir_prepared_long.csv"))

# -----------------------------------------------------------------------------
# §2.1 + §2.5 Above-water summary and flag (100%, 5 h)
# -----------------------------------------------------------------------------

sir100_5 <- sir |> filter(Irrigation_trmt_pct == 100, Timepoint == 5)

water_summary <- sir100_5 |>
  group_by(Substrate, group) |>
  summarise(n_samples = n(),
            median_induced_gdw = median(induced_gdw, na.rm = TRUE),
            share_above_zero = mean(induced_gdw > 0, na.rm = TRUE),
            .groups = "drop") |>
  mutate(flag_below_75pct = share_above_zero < 0.75) |>
  arrange(factor(group, c("Simple", "Polymer", "Reference", "Q2", "Excluded")),
          Substrate)
write_csv(water_summary, file.path(dir_results, "cp1_water_subtraction_summary.csv"))

flagged_compounds <- water_summary |>
  filter(flag_below_75pct, group %in% c("Simple", "Polymer", "Reference"))

# -----------------------------------------------------------------------------
# §2.6 + §2.7 Scores (100%, 5 h)
# -----------------------------------------------------------------------------

scores_C <- make_scores(sir100_5, "induced_C") |>
  mutate(normalization = "per_C")
scores_MBC <- make_scores(sir100_5, "induced_C_per_MBC") |>
  mutate(normalization = "per_MBC")

sample_scores <- bind_rows(scores_C, scores_MBC) |>
  left_join(mbc |> select(Sample_cat, MBC_flag), by = "Sample_cat")
write_csv(sample_scores, file.path(dir_results, "sample_scores_100pct_5h.csv"))

date_scores <- sample_scores |>
  group_by(normalization) |>
  group_modify(\(d, k) date_means(d)) |>
  ungroup()
write_csv(date_scores, file.path(dir_results, "date_scores_100pct_5h.csv"))

# -----------------------------------------------------------------------------
# §2.8 Root-dense cores (100%, active season)
# -----------------------------------------------------------------------------
# M = 0.6745 * (x - median) / MAD, with MAD the raw median absolute deviation
# (constant = 1, not R's default 1.4826).

roots <- bgc |>
  filter(Irrigation_trmt_pct == 100, Collection_date %in% study_dates) |>
  mutate(season = assign_season(Collection_date)) |>
  filter(season == "Active") |>
  select(Sample_cat, Collection_date, Plot, Roots_gdw) |>
  mutate(M = 0.6745 * (Roots_gdw - median(Roots_gdw, na.rm = TRUE)) /
           mad(Roots_gdw, constant = 1, na.rm = TRUE),
         root_flag = abs(M) > 3.5) |>
  arrange(desc(abs(M)))
write_csv(roots, file.path(dir_results, "root_flags_active_100pct.csv"))

# -----------------------------------------------------------------------------
# §2.9 Sample availability and unmatched samples
# -----------------------------------------------------------------------------

mr_samples <- unique(microresp$Sample_cat)
field_vars <- c("Roots_gdw", "SoilTemp_mean", "GWC", "SEOC_ug_gdw")

availability <- samples |>
  left_join(bgc |> select(Sample_cat, all_of(field_vars)), by = "Sample_cat") |>
  left_join(eea |> select(Sample_cat, BG_nmol_gdw_h, MBC_ug_gdw),
            by = "Sample_cat") |>
  mutate(microresp = Sample_cat %in% mr_samples,
         eea = !is.na(BG_nmol_gdw_h),
         mbc = !is.na(MBC_ug_gdw),
         field = if_all(all_of(field_vars), \(x) !is.na(x)),
         missing_field = pmap_chr(pick(all_of(field_vars)), \(...) {
           paste(field_vars[is.na(c(...))], collapse = ", ")
         }),
         transcripts = Sample_cat %in% transcript_samples)

sample_counts <- availability |>
  group_by(Irrigation_trmt_pct, Collection_date, season) |>
  summarise(plots = n(), across(c(microresp, eea, mbc, field, transcripts), sum),
            .groups = "drop")
write_csv(sample_counts, file.path(dir_results, "cp1_sample_counts.csv"))

unmatched <- bind_rows(
  availability |>
    filter(!(microresp & eea & mbc & field)) |>
    mutate(issue = paste0(
      if_else(!microresp, "no MicroResp; ", ""),
      if_else(!eea, "no EEA; ", ""),
      if_else(eea & !mbc, "MBC missing; ", ""),
      if_else(!field, paste0("missing field: ", missing_field), "")
    )) |>
    select(Sample_cat, Collection_date, Plot, Irrigation_trmt_pct, issue),
  tibble(Sample_cat = setdiff(mr_samples, samples$Sample_cat),
         issue = "MicroResp sample not in field table on study dates"),
  tibble(Sample_cat = setdiff(transcript_samples, samples$Sample_cat),
         issue = "transcript sample not in field table on study dates")
) |>
  arrange(Irrigation_trmt_pct, Collection_date)
write_csv(unmatched, file.path(dir_results, "cp1_unmatched_samples.csv"))

# Transcripts: plots without a sample on transcript dates (100%)
transcript_dates <- availability |>
  filter(transcripts) |>
  distinct(Collection_date) |>
  pull()
transcript_gaps <- availability |>
  filter(Irrigation_trmt_pct == 100, Collection_date %in% transcript_dates,
         !transcripts) |>
  select(Sample_cat, Collection_date, Plot)
write_csv(transcript_gaps, file.path(dir_results, "cp1_transcript_gaps.csv"))

# -----------------------------------------------------------------------------
# Console summary
# -----------------------------------------------------------------------------

print(percarbon_summary)
print(percarbon_by_substrate, n = Inf)
print(water_summary, n = Inf)
print(flagged_compounds)
print(mbc_replacements)
print(filter(roots, root_flag))
print(sample_counts, n = Inf)
print(unmatched, n = Inf)
print(transcript_gaps)
