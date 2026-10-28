# =============================================================================
# 01. Step 0 checks and shared MicroResp tables (brief Step 0, rules 4-5)
# =============================================================================
#
# 100% plots, 13 dates. Writes:
#   data/microresp_long.csv   one row per sample x substrate x timepoint, with
#                             the primary, gross, fold-over-water, and per-MBC
#                             responses
#   data/z_5h.csv             z-scores and relative responses (5 h)
#   data/step0_*.csv          assay QC, MBC rule, root screens, negatives
#   data/field_100.csv        field context table for Figure 1 and S1

source(file.path("Scripts", "storyboard_v6", "00_setup.R"))

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------

microresp <- read_csv(input_files$microresp, locale = locale(encoding = "latin1"),
                      na = c("", "NA", "na"), show_col_types = FALSE) |>
  mutate(Collection_date = mdy(Collection_date))

bgc <- read_csv(input_files$bgc, na = c("", "NA", "na"), show_col_types = FALSE) |>
  mutate(Collection_date = as.Date(Collection_date))

# Field soil respiration rate exists only in the v4 cleaned table (computed by
# 01_Data_merging_and_cleaning.R); other field variables match DataHub exactly.
field_clean <- read_csv(file.path("Data", "Data_merged_cleaned.csv"),
                        show_col_types = FALSE) |>
  filter(Irrigation_trmt_pct == 100) |>
  select(Sample_cat, Resp_rate_mg_C_m2_h)

samples <- bgc |>
  filter(Irrigation_trmt_pct == 100, Collection_date %in% study_dates) |>
  select(Sample_cat, Collection_date, Plot, SoilTemp_mean, GWC, Roots_gdw,
         MBC_ug_gdw) |>
  left_join(field_clean, by = "Sample_cat") |>
  mutate(season = assign_season(Collection_date))
stopifnot(nrow(samples) == 78)

stopifnot(setequal(setdiff(unique(microresp$Substrate), "Water"),
                   substrate_classes$substrate))

# -----------------------------------------------------------------------------
# MicroResp responses (rule 4)
# -----------------------------------------------------------------------------
# primary = (substrate CO2 - water CO2, same sample) / C added
# gross   = substrate CO2 / C added
# fold    = substrate CO2 / water CO2
# Per-gram rates are converted to per mg C with gdw_well / mgC_per_well, as in
# plan_v1 (D5: mgC_per_well column; open item O3).

mr_gdw <- microresp |>
  filter(Irrigation_trmt_pct == 100, Collection_date %in% study_dates,
         Normalization == "ngCO2-C_gdw_h") |>
  select(Sample_cat, Collection_date, Plot, Timepoint, Substrate, gdw_well,
         mgC_per_well, SIR_gdw = SIR)

water <- mr_gdw |>
  filter(Substrate == "Water") |>
  select(Sample_cat, Timepoint, water_gdw = SIR_gdw)

mr_long <- mr_gdw |>
  filter(Substrate != "Water") |>
  left_join(water, by = c("Sample_cat", "Timepoint")) |>
  left_join(samples |> select(Sample_cat, season, MBC_ug_gdw), by = "Sample_cat") |>
  mutate(primary = (SIR_gdw - water_gdw) * gdw_well / mgC_per_well,
         gross = SIR_gdw * gdw_well / mgC_per_well,
         fold = SIR_gdw / water_gdw,
         mgC_per_gdw = mgC_per_well / gdw_well)

# Agreement with the plan_v1 prepared table (same formula)
plan_v1_sir <- read_csv(file.path(dir_results, "sir_prepared_long.csv"),
                        show_col_types = FALSE) |>
  filter(Irrigation_trmt_pct == 100, Collection_date %in% study_dates) |>
  select(Sample_cat, Timepoint, Substrate, induced_C)
check_v1 <- mr_long |>
  inner_join(plan_v1_sir, by = c("Sample_cat", "Timepoint", "Substrate"))
stopifnot(nrow(check_v1) == nrow(mr_long),
          isTRUE(all.equal(check_v1$primary, check_v1$induced_C)))

# -----------------------------------------------------------------------------
# Step 0.2 Microbial biomass: at or near zero -> missing
# -----------------------------------------------------------------------------
# Rule: MBC = 0 (below the process blank, reported as zero in DataHub) or
# below 1% of the 100% median MBC is set to missing.

mbc_median <- median(samples$MBC_ug_gdw, na.rm = TRUE)
mbc_threshold <- 0.01 * mbc_median
mbc_check <- samples |>
  select(Sample_cat, Collection_date, Plot, MBC_ug_gdw) |>
  arrange(MBC_ug_gdw) |>
  mutate(rank = row_number(),
         times_below_median = mbc_median / MBC_ug_gdw,
         set_missing = MBC_ug_gdw < mbc_threshold)
write_csv(mbc_check |> slice_head(n = 6), file.path(dir_v6_data, "step0_mbc_lowest.csv"))
mbc_missing <- mbc_check$Sample_cat[mbc_check$set_missing]

samples <- samples |>
  mutate(MBC_used = if_else(Sample_cat %in% mbc_missing, NA_real_, MBC_ug_gdw))
mr_long <- mr_long |>
  mutate(MBC_used = if_else(Sample_cat %in% mbc_missing, NA_real_, MBC_ug_gdw),
         primary_per_MBC = primary / MBC_used)

# -----------------------------------------------------------------------------
# Step 0.3 Root outliers: whole-dataset and within-plot modified z screens
# -----------------------------------------------------------------------------

modified_z <- function(x) 0.6745 * (x - median(x)) / mad(x, constant = 1)

root_screen <- samples |>
  select(Sample_cat, Collection_date, Plot, season, Roots_gdw) |>
  mutate(z_whole = modified_z(Roots_gdw)) |>
  group_by(Plot) |>
  mutate(z_within_plot = modified_z(Roots_gdw)) |>
  ungroup() |>
  mutate(flag_whole = abs(z_whole) > 3.5,
         flag_within_plot = abs(z_within_plot) > 3.5,
         flag_both = flag_whole & flag_within_plot)
write_csv(root_screen |> filter(flag_whole | flag_within_plot) |> arrange(desc(z_whole)),
          file.path(dir_v6_data, "step0_root_flags.csv"))
root_flags <- root_screen$Sample_cat[root_screen$flag_both]
if (!setequal(root_flags, c("SM1621", "SM1078", "SM1080"))) {
  warning("Root flags differ from the expected SM1621, SM1078, SM1080: ",
          paste(root_flags, collapse = ", "))
}

samples <- samples |> left_join(root_screen |> select(Sample_cat, flag_both),
                                by = "Sample_cat")
write_csv(samples, file.path(dir_v6_data, "field_100.csv"))

# -----------------------------------------------------------------------------
# Step 0.1 Assay QC by date
# -----------------------------------------------------------------------------
# The MicroResp file has no run date or incubation temperature. The DataHub
# README states incubation at room temperature and a CO2 calculation at 20 C,
# and MicroResp subsamples were run within a week of collection. QC is
# therefore reported by collection date.

water_5h <- mr_gdw |>
  filter(Substrate == "Water", Timepoint == 5) |>
  left_join(samples |> select(Sample_cat, season, GWC), by = "Sample_cat")

assay_qc <- mr_long |>
  filter(Timepoint == 5) |>
  group_by(Collection_date, season) |>
  summarise(n_samples = n_distinct(Sample_cat),
            gdw_well_mean = mean(gdw_well[!duplicated(Sample_cat)]),
            gdw_well_min = min(gdw_well), gdw_well_max = max(gdw_well),
            glucose_mgC_per_gdw_mean = mean(mgC_per_gdw[Substrate == "D.Glucose"]),
            glucose_mgC_per_gdw_min = min(mgC_per_gdw[Substrate == "D.Glucose"]),
            glucose_mgC_per_gdw_max = max(mgC_per_gdw[Substrate == "D.Glucose"]),
            .groups = "drop") |>
  left_join(water_5h |>
              group_by(Collection_date) |>
              summarise(water_CO2_mean = mean(SIR_gdw), water_CO2_sd = sd(SIR_gdw),
                        GWC_mean = mean(GWC), .groups = "drop"),
            by = "Collection_date") |>
  mutate(incubation_temperature = "not recorded (README: room temperature; CO2 calc assumes 20 C)")
write_csv(assay_qc, file.path(dir_v6_data, "step0_assay_qc.csv"))
write_csv(water_5h, file.path(dir_v6_data, "water_control_5h.csv"))

# Same amount of C per well on every date? (mgC_per_well is one value per substrate)
mgc_constant <- microresp |>
  filter(Substrate != "Water") |>
  group_by(Substrate) |>
  summarise(n_values = n_distinct(mgC_per_well), .groups = "drop")
stopifnot(all(mgc_constant$n_values == 1))

# Season comparison of soil per well (date means; per g dry soil depends on it)
gdw_dates <- mr_gdw |>
  distinct(Sample_cat, Collection_date, gdw_well) |>
  mutate(season = assign_season(Collection_date)) |>
  group_by(Collection_date, season) |>
  summarise(gdw = mean(gdw_well), .groups = "drop")
gdw_season <- season_test(gdw_dates$gdw, gdw_dates$season)
write_csv(gdw_season, file.path(dir_v6_data, "step0_gdw_well_season.csv"))

# -----------------------------------------------------------------------------
# Rule 4: negative primary values by substrate (kept as they are)
# -----------------------------------------------------------------------------

negatives <- mr_long |>
  group_by(Timepoint, Substrate) |>
  summarise(n_values = sum(!is.na(primary)), n_negative = sum(primary < 0, na.rm = TRUE),
            .groups = "drop") |>
  mutate(share_negative = n_negative / n_values,
         Substrate = factor(Substrate, substrate_levels)) |>
  arrange(Timepoint, Substrate)
write_csv(negatives, file.path(dir_v6_data, "step0_negative_values.csv"))

# Fold over water: water CO2 at or below zero makes the ratio undefined
fold_problems <- water |> filter(water_gdw <= 0)
write_csv(fold_problems, file.path(dir_v6_data, "step0_water_nonpositive.csv"))

write_csv(mr_long, file.path(dir_v6_data, "microresp_long.csv"))

# -----------------------------------------------------------------------------
# Rules 4-5: z-scores and relative response (5 h and 24 h)
# -----------------------------------------------------------------------------
# z-score each substrate across all samples of that timepoint; relative =
# z minus the mean z of all substrates in the same sample.

make_z <- function(df, value) {
  df |>
    filter(!is.na(.data[[value]])) |>
    group_by(Substrate) |>
    mutate(z = zscore(.data[[value]])) |>
    group_by(Sample_cat) |>
    mutate(level = mean(z), relative = z - level) |>
    ungroup() |>
    left_join(substrate_classes |> select(Substrate = substrate, chem_class),
              by = "Substrate") |>
    select(Sample_cat, Collection_date, Plot, season, Substrate, chem_class,
           value = all_of(value), z, level, relative)
}

z_all <- bind_rows(
  make_z(filter(mr_long, Timepoint == 5), "primary") |> mutate(version = "primary_5h"),
  make_z(filter(mr_long, Timepoint == 24), "primary") |> mutate(version = "primary_24h"),
  make_z(filter(mr_long, Timepoint == 5), "gross") |> mutate(version = "gross_5h"),
  make_z(filter(mr_long, Timepoint == 5, is.finite(fold)), "fold") |>
    mutate(version = "fold_5h"),
  make_z(filter(mr_long, Timepoint == 5), "primary_per_MBC") |>
    mutate(version = "per_MBC_5h")
)
write_csv(z_all, file.path(dir_v6_data, "microresp_z.csv"))

# -----------------------------------------------------------------------------
# Console summary for the log
# -----------------------------------------------------------------------------

print(assay_qc, width = Inf)
print(gdw_season)
print(mbc_check |> slice_head(n = 5))
print(root_screen |> filter(flag_whole | flag_within_plot))
print(negatives |> filter(Timepoint == 5), n = Inf)
print(fold_problems)
print(z_all |> group_by(version) |> summarise(n_samples = n_distinct(Sample_cat)))
