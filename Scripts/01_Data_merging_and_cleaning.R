# ============================================================================
# 01 Data merging and cleaning
# Soils SFA A3.1 Prosser Time Series
# ============================================================================
#
# Purpose:
#   Build a reproducible cleaned biogeochemistry + EEA dataset.
#
# Cleaning decisions carried forward from the existing exploration workflow:
#   1. Merge soil biogeochemistry and EEA data by Sample_cat.
#   2. Parse Collection_date robustly and create Month, Year, Month_name,
#      irrigation-based Season, and Irrigation_active.
#   3. Remove the July 2023 outlier samples identified in the prior
#      exploration workflow: SM1082, SM1085, SM1087, SM1088.
#   4. Treat zeros in selected C, N, P, and microbial biomass variables as
#      below detection limit values. Replace those zeros with LOD/2, where
#      LOD/2 is half of the minimum non-zero value for that variable.
#   5. If negative field respiration values exist, replace them with half of
#      the minimum non-negative respiration value.
#
# Decision rationale:
#   The source README notes that negative analytical values in carbon,
#   nitrogen, and phosphorus analyses were converted to zeros and can be
#   interpreted as values below detection limit. LOD/2 preserves samples while
#   avoiding true zeros in downstream transformations and multivariate models.
#   The fixed outlier list is used so this script reproduces the existing
#   documented cleaning decision instead of changing results whenever the data
#   distribution changes.

rm(list = ls())

library(tidyverse)
library(lubridate)

# ----------------------------------------------------------------------------
# 1. File paths and output folders
# ----------------------------------------------------------------------------

bgc_source_candidates <- c(
  "DataHub/Soil_Biogeochemistry/SoilSFA_A3.1_Soil_Biogeochemistry.csv",
  "DataHub March/Soil_Biogeochemistry/SoilMicrobiomeSFA_prosser-time-series_Soil_Biogeochemistry.csv"
)

eea_source_candidates <- c(
  "DataHub/Soil_Phenotypes/SoilSFA_A3.1_EEA.csv",
  "DataHub March/Soil_Phenotypes/SoilMicrobiomeSFA_prosser-time-series_EEA.csv"
)

bgc_source <- bgc_source_candidates[file.exists(bgc_source_candidates)][1]
eea_source <- eea_source_candidates[file.exists(eea_source_candidates)][1]

if (is.na(bgc_source)) {
  stop(
    "Could not find the soil biogeochemistry source file. Checked: ",
    paste(bgc_source_candidates, collapse = "; ")
  )
}

if (is.na(eea_source)) {
  stop(
    "Could not find the EEA source file. Checked: ",
    paste(eea_source_candidates, collapse = "; ")
  )
}

dir.create("Data", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Data_merging_and_cleaning", showWarnings = FALSE, recursive = TRUE)

cat("\n=== SOURCE FILES ===\n")
cat("Biogeochemistry:", bgc_source, "\n")
cat("EEA:", eea_source, "\n")

# ----------------------------------------------------------------------------
# 2. Load source data
# ----------------------------------------------------------------------------

bgc_data <- read_csv(bgc_source, show_col_types = FALSE, guess_max = 10000)
eea_data <- read_csv(eea_source, show_col_types = FALSE, guess_max = 10000)

if ("Irrigaiton_trmt_cat" %in% names(eea_data) &&
    !"Irrigation_trmt_cat" %in% names(eea_data)) {
  eea_data <- eea_data %>%
    rename(Irrigation_trmt_cat = Irrigaiton_trmt_cat)
}

duplicate_bgc_samples <- bgc_data %>%
  count(Sample_cat, name = "n") %>%
  filter(n > 1)

duplicate_eea_samples <- eea_data %>%
  count(Sample_cat, name = "n") %>%
  filter(n > 1)

if (nrow(duplicate_bgc_samples) > 0) {
  write_csv(
    duplicate_bgc_samples,
    "Results/Data_merging_and_cleaning/duplicate_bgc_sample_ids.csv"
  )
  stop("Duplicate Sample_cat values found in the biogeochemistry source file.")
}

if (nrow(duplicate_eea_samples) > 0) {
  write_csv(
    duplicate_eea_samples,
    "Results/Data_merging_and_cleaning/duplicate_eea_sample_ids.csv"
  )
  stop("Duplicate Sample_cat values found in the EEA source file.")
}

cat("\n=== SOURCE DATA CHECK ===\n")
cat("BGC rows:", nrow(bgc_data), "\n")
cat("EEA rows:", nrow(eea_data), "\n")

# ----------------------------------------------------------------------------
# 3. Parse dates and create common metadata fields
# ----------------------------------------------------------------------------

bgc_data <- bgc_data %>%
  mutate(
    Collection_date = parse_date_time(
      as.character(Collection_date),
      orders = c("ymd", "mdy", "dmy", "Ymd", "mdY", "dmY")
    ) %>%
      as.Date(),
    Month = month(Collection_date),
    Year = year(Collection_date),
    Month_name = month(Collection_date, label = TRUE, abbr = TRUE),
    Season = case_when(
      Month %in% c(4, 5, 6, 7, 8, 9, 10) ~ "Growing",
      Month %in% c(11, 12, 1, 2, 3) ~ "Non-Growing",
      TRUE ~ NA_character_
    ),
    Season = factor(Season, levels = c("Growing", "Non-Growing")),
    Irrigation_active = if_else(Month %in% c(4, 5, 6, 7, 8, 9, 10), 1, 0),
    Resp_rate_mg_C_m2_h = CO2_C_60_m2 / 1000,
    Irrigation_trmt_pct = factor(Irrigation_trmt_pct, levels = c("25", "100"))
  )

date_parse_failures <- bgc_data %>%
  filter(is.na(Collection_date)) %>%
  select(Sample_cat)

if (nrow(date_parse_failures) > 0) {
  write_csv(
    date_parse_failures,
    "Results/Data_merging_and_cleaning/date_parse_failures.csv"
  )
  stop("Some Collection_date values could not be parsed. See date_parse_failures.csv.")
}

# ----------------------------------------------------------------------------
# 4. Merge BGC and EEA measurements
# ----------------------------------------------------------------------------

eea_measure_cols <- c(
  "Sample_cat",
  "BG_nmol_gdw_h",
  "BX_nmol_gdw_h",
  "CB_nmol_gdw_h",
  "AP_nmol_gdw_h",
  "NAG_nmol_gdw_h",
  "Leu_nmol_gdw_h",
  "Ala_nmol_gdw_h",
  "Sum_Carbon",
  "Sum_Nitrogen",
  "Sum_Phosphorus",
  "EEA_C_N",
  "EEA_C_P",
  "EEA_N_P",
  "MBC_ug_gdw",
  "MBN_ug_gdw"
)

if (!"Irrigation_trmt_cat" %in% names(bgc_data) &&
    "Irrigation_trmt_cat" %in% names(eea_data)) {
  eea_measure_cols <- c(eea_measure_cols, "Irrigation_trmt_cat")
}

eea_for_merge <- eea_data %>%
  select(any_of(eea_measure_cols))

samples_missing_eea <- bgc_data %>%
  anti_join(eea_for_merge, by = "Sample_cat") %>%
  select(Sample_cat, Collection_date, Plot, Irrigation_trmt_pct)

eea_without_bgc <- eea_for_merge %>%
  anti_join(bgc_data, by = "Sample_cat") %>%
  select(Sample_cat)

write_csv(samples_missing_eea, "Results/Data_merging_and_cleaning/bgc_samples_without_eea.csv")
write_csv(eea_without_bgc, "Results/Data_merging_and_cleaning/eea_samples_without_bgc.csv")

rf_data <- bgc_data %>%
  left_join(
    eea_for_merge,
    by = "Sample_cat",
    suffix = c("_bgc", "_eea")
  )

if (!"MBC_ug_gdw_bgc" %in% names(rf_data) && "MBC_ug_gdw" %in% names(rf_data)) {
  rf_data <- rf_data %>%
    rename(MBC_ug_gdw_bgc = MBC_ug_gdw)
}

if (!"MBN_ug_gdw_bgc" %in% names(rf_data) && "MBN_ug_gdw" %in% names(rf_data)) {
  rf_data <- rf_data %>%
    rename(MBN_ug_gdw_bgc = MBN_ug_gdw)
}

cat("\n=== MERGE CHECK ===\n")
cat("Merged rows:", nrow(rf_data), "\n")
cat("BGC samples without EEA:", nrow(samples_missing_eea), "\n")
cat("EEA samples without BGC:", nrow(eea_without_bgc), "\n")
cat("Samples with respiration data:", sum(!is.na(rf_data$Resp_rate_mg_C_m2_h)), "\n")

# ----------------------------------------------------------------------------
# 5. Define analysis variables
# ----------------------------------------------------------------------------

soil_chemistry <- c("GWC", "pH", "SEOC_ug_gdw", "SETN_ug_gdw")
microbial_biomass <- c("MBC_ug_gdw_bgc", "MBN_ug_gdw_bgc")
extractable_p <- c("Olson_IP_ug_gdw", "Oslon_OP_ug_gdw",
                   "MBP_IP_ug_gdw", "MBP_OP_ug_gdw")
enzymes_c <- c("BG_nmol_gdw_h", "BX_nmol_gdw_h", "CB_nmol_gdw_h")
enzymes_n <- c("NAG_nmol_gdw_h", "Leu_nmol_gdw_h", "Ala_nmol_gdw_h")
enzymes_p <- c("AP_nmol_gdw_h")
enzyme_sums_ratios <- c("Sum_Carbon", "Sum_Nitrogen", "Sum_Phosphorus",
                        "EEA_C_N", "EEA_C_P", "EEA_N_P")
environmental <- c("SoilTemp_mean", "SoilTemp_min", "SoilTemp_max",
                   "Roots_gdw")

all_predictors <- c(
  soil_chemistry,
  microbial_biomass,
  extractable_p,
  enzymes_c,
  enzymes_n,
  enzymes_p,
  enzyme_sums_ratios,
  environmental
)

response_var <- "Resp_rate_mg_C_m2_h"

required_cols <- c(
  "Sample_cat",
  "Collection_date",
  "Month",
  "Month_name",
  "Season",
  "Year",
  "Plot",
  "Irrigation_trmt_pct",
  "Irrigation_trmt_cat",
  "Block",
  "Rep",
  "Irrigation_active",
  all_predictors,
  response_var
)

missing_required_cols <- setdiff(required_cols, names(rf_data))

if (length(missing_required_cols) > 0) {
  stop(
    "The merged dataset is missing required columns: ",
    paste(missing_required_cols, collapse = ", ")
  )
}

rf_analysis <- rf_data %>%
  select(all_of(required_cols)) %>%
  filter(!is.na(.data[[response_var]]))

cat("\n=== ANALYSIS DATASET BEFORE CLEANING ===\n")
cat("Samples with respiration:", nrow(rf_analysis), "\n")
cat("Date range:",
    format(min(rf_analysis$Collection_date), "%Y-%m-%d"), "to",
    format(max(rf_analysis$Collection_date), "%Y-%m-%d"), "\n")

# ----------------------------------------------------------------------------
# 6. Document missing data before cleaning
# ----------------------------------------------------------------------------

missing_summary <- rf_analysis %>%
  select(all_of(all_predictors), all_of(response_var)) %>%
  summarise(across(
    everything(),
    list(
      n_missing = ~sum(is.na(.)),
      pct_missing = ~round(sum(is.na(.)) / n() * 100, 1)
    )
  )) %>%
  pivot_longer(
    everything(),
    names_to = c("Variable", ".value"),
    names_pattern = "(.+)_(n_missing|pct_missing)"
  ) %>%
  arrange(desc(pct_missing), Variable)

write_csv(missing_summary, "Results/Data_merging_and_cleaning/missing_data_summary_before_cleaning.csv")

# ----------------------------------------------------------------------------
# 7. Remove documented outlier samples
# ----------------------------------------------------------------------------

outlier_samples <- c("SM1082", "SM1085", "SM1087", "SM1088")

outliers_removed <- rf_analysis %>%
  filter(Sample_cat %in% outlier_samples) %>%
  select(Sample_cat, Collection_date, Month_name, Plot, Irrigation_trmt_pct,
         SEOC_ug_gdw, SETN_ug_gdw, MBP_IP_ug_gdw, MBP_OP_ug_gdw,
         Resp_rate_mg_C_m2_h)

rf_no_outliers <- rf_analysis %>%
  filter(!Sample_cat %in% outlier_samples)

write_csv(outliers_removed, "Results/Data_merging_and_cleaning/outlier_samples_removed.csv")

cat("\n=== OUTLIER REMOVAL ===\n")
cat("Outlier samples requested for removal:", paste(outlier_samples, collapse = ", "), "\n")
cat("Outlier samples found and removed:", nrow(outliers_removed), "\n")
cat("Samples after outlier removal:", nrow(rf_no_outliers), "\n")

# ----------------------------------------------------------------------------
# 8. Analyze and replace zeros with LOD/2
# ----------------------------------------------------------------------------

vars_with_zeros_to_treat <- c(
  "MBP_OP_ug_gdw",
  "Oslon_OP_ug_gdw",
  "Olson_IP_ug_gdw",
  "MBP_IP_ug_gdw",
  "SEOC_ug_gdw",
  "MBC_ug_gdw_bgc",
  "SETN_ug_gdw",
  "MBN_ug_gdw_bgc"
)

vars_with_zeros_to_treat <- intersect(vars_with_zeros_to_treat, names(rf_no_outliers))

zero_summary <- rf_no_outliers %>%
  select(all_of(vars_with_zeros_to_treat)) %>%
  summarise(across(
    everything(),
    list(
      n_zeros = ~sum(. == 0, na.rm = TRUE),
      pct_zeros = ~round(sum(. == 0, na.rm = TRUE) / n() * 100, 1),
      min_nonzero = ~min(.[. > 0], na.rm = TRUE)
    )
  )) %>%
  pivot_longer(
    everything(),
    names_to = c("Variable", ".value"),
    names_pattern = "(.+)_(n_zeros|pct_zeros|min_nonzero)"
  ) %>%
  mutate(
    min_nonzero = if_else(is.infinite(min_nonzero), NA_real_, min_nonzero),
    LOD_half = min_nonzero / 2
  ) %>%
  arrange(desc(n_zeros), Variable)

zero_replacement_values <- zero_summary %>%
  filter(n_zeros > 0, !is.na(LOD_half))

zero_values_detail <- rf_no_outliers %>%
  select(Sample_cat, Collection_date, Month_name, Season,
         Irrigation_trmt_pct, all_of(vars_with_zeros_to_treat)) %>%
  pivot_longer(
    cols = all_of(vars_with_zeros_to_treat),
    names_to = "Variable",
    values_to = "Value"
  ) %>%
  filter(Value == 0) %>%
  arrange(Sample_cat, Variable)

write_csv(zero_summary, "Results/Data_merging_and_cleaning/zero_values_summary_before_replacement.csv")
write_csv(zero_values_detail, "Results/Data_merging_and_cleaning/zero_values_detail_before_replacement.csv")
write_csv(zero_replacement_values, "Results/Data_merging_and_cleaning/zero_replacement_values_lod_half.csv")

rf_cleaned <- rf_no_outliers

cat("\n=== ZERO REPLACEMENT ===\n")
cat("Variables checked for below-detection zeros:",
    paste(vars_with_zeros_to_treat, collapse = ", "), "\n")
cat("Variables with zeros replaced:", nrow(zero_replacement_values), "\n")

for (i in seq_len(nrow(zero_replacement_values))) {
  var <- zero_replacement_values$Variable[i]
  replacement_value <- zero_replacement_values$LOD_half[i]
  n_zeros_before <- sum(rf_cleaned[[var]] == 0, na.rm = TRUE)

  rf_cleaned[[var]][rf_cleaned[[var]] == 0] <- replacement_value

  cat(sprintf(
    "  %s: %d zeros replaced with %.6f\n",
    var,
    n_zeros_before,
    replacement_value
  ))
}

remaining_zeros_treated_vars <- rf_cleaned %>%
  select(all_of(vars_with_zeros_to_treat)) %>%
  summarise(across(everything(), ~sum(. == 0, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "n_zeros_remaining")

write_csv(
  remaining_zeros_treated_vars,
  "Results/Data_merging_and_cleaning/remaining_zeros_after_replacement.csv"
)

# ----------------------------------------------------------------------------
# 9. Replace negative respiration values if present
# ----------------------------------------------------------------------------

negative_respiration_samples <- rf_cleaned %>%
  filter(Resp_rate_mg_C_m2_h < 0) %>%
  select(Sample_cat, Collection_date, Month_name, Season,
         Irrigation_trmt_pct, Resp_rate_mg_C_m2_h)

n_negative_respiration <- nrow(negative_respiration_samples)
respiration_replacement_value <- NA_real_

cat("\n=== NEGATIVE RESPIRATION CHECK ===\n")
cat("Negative respiration values found:", n_negative_respiration, "\n")

if (n_negative_respiration > 0) {
  min_non_negative_respiration <- min(
    rf_cleaned$Resp_rate_mg_C_m2_h[rf_cleaned$Resp_rate_mg_C_m2_h >= 0],
    na.rm = TRUE
  )

  respiration_replacement_value <- min_non_negative_respiration / 2

  rf_cleaned <- rf_cleaned %>%
    mutate(
      Resp_rate_mg_C_m2_h = if_else(
        Resp_rate_mg_C_m2_h < 0,
        respiration_replacement_value,
        Resp_rate_mg_C_m2_h
      )
    )

  cat("Replacement value:", respiration_replacement_value, "\n")
}

write_csv(
  negative_respiration_samples,
  "Results/Data_merging_and_cleaning/negative_respiration_values_before_replacement.csv"
)

# ----------------------------------------------------------------------------
# 10. Create complete-case version for analyses that require no missing values
# ----------------------------------------------------------------------------

rf_complete_cleaned <- rf_cleaned %>%
  drop_na(all_of(all_predictors), all_of(response_var))

cat("\n=== CLEANED DATASET SUMMARY ===\n")
cat("Samples before cleaning:", nrow(rf_analysis), "\n")
cat("Samples after outlier removal:", nrow(rf_no_outliers), "\n")
cat("Samples in cleaned dataset:", nrow(rf_cleaned), "\n")
cat("Complete cases after cleaning:", nrow(rf_complete_cleaned), "\n")

final_stratification <- rf_cleaned %>%
  group_by(Season, Irrigation_trmt_pct) %>%
  summarise(
    N = n(),
    N_months = n_distinct(Month_name),
    Mean_respiration = round(mean(Resp_rate_mg_C_m2_h, na.rm = TRUE), 3),
    Min_respiration = round(min(Resp_rate_mg_C_m2_h, na.rm = TRUE), 3),
    Max_respiration = round(max(Resp_rate_mg_C_m2_h, na.rm = TRUE), 3),
    .groups = "drop"
  )

write_csv(final_stratification, "Results/Data_merging_and_cleaning/final_stratification_summary.csv")

# ----------------------------------------------------------------------------
# 11. Save cleaned data and decision logs
# ----------------------------------------------------------------------------

cleaning_decisions <- tibble(
  Decision = c(
    "Merge BGC and EEA",
    "Parse Collection_date",
    "Remove documented outliers",
    "Replace zeros with LOD/2",
    "Replace negative respiration",
    "Keep complete-case output as a separate file"
  ),
  Rule = c(
    "Left join EEA measurements to BGC records by Sample_cat.",
    "Parse common ymd/mdy/dmy formats, then derive Month, Year, Month_name, Season, and Irrigation_active.",
    paste("Drop fixed sample list:", paste(outlier_samples, collapse = ", ")),
    "For selected below-detection variables, replace zero with half of that variable's minimum non-zero value.",
    "If field respiration is negative, replace it with half of the minimum non-negative respiration value.",
    "Do not drop incomplete rows from the main cleaned output; save a second complete-case dataset."
  ),
  Rationale = c(
    "Sample_cat is the shared sample identifier used by the existing analysis scripts.",
    "Date formats can differ across source exports; robust parsing avoids silent date failures.",
    "These July 2023 samples were identified in the prior exploration workflow as extreme outliers.",
    "The README documents converted negative analytical values as zeros/below detection; LOD/2 preserves sample size.",
    "Negative field respiration is not meaningful for downstream positive-scale analyses.",
    "Some analyses can use partial data, while RF/RDA-style workflows may require complete predictors."
  ),
  Affected = c(
    paste0("BGC rows: ", nrow(bgc_data), "; EEA rows: ", nrow(eea_data),
           "; merged rows: ", nrow(rf_data)),
    paste0("Parsed dates for ", nrow(bgc_data), " BGC records."),
    paste0(nrow(outliers_removed), " samples removed."),
    paste0(sum(zero_replacement_values$n_zeros), " zero values replaced across ",
           nrow(zero_replacement_values), " variables."),
    paste0(n_negative_respiration, " values replaced."),
    paste0(nrow(rf_complete_cleaned), " complete cases saved separately.")
  )
)

cleaned_data_summary <- tibble(
  Metric = c(
    "BGC source",
    "EEA source",
    "BGC rows",
    "EEA rows",
    "Merged rows",
    "Rows with respiration before cleaning",
    "Outliers removed",
    "Rows in cleaned data",
    "Complete cases after cleaning",
    "Zero values replaced",
    "Negative respiration values replaced",
    "Date range start",
    "Date range end"
  ),
  Value = c(
    bgc_source,
    eea_source,
    as.character(nrow(bgc_data)),
    as.character(nrow(eea_data)),
    as.character(nrow(rf_data)),
    as.character(nrow(rf_analysis)),
    as.character(nrow(outliers_removed)),
    as.character(nrow(rf_cleaned)),
    as.character(nrow(rf_complete_cleaned)),
    as.character(sum(zero_replacement_values$n_zeros)),
    as.character(n_negative_respiration),
    format(min(rf_cleaned$Collection_date), "%Y-%m-%d"),
    format(max(rf_cleaned$Collection_date), "%Y-%m-%d")
  )
)

write_csv(rf_cleaned, "Data/Data_merged_cleaned.csv")
write_csv(rf_complete_cleaned, "Data/Data_merged_cleaned_complete_cases.csv")
write_csv(cleaning_decisions, "Results/Data_merging_and_cleaning/cleaning_decisions.csv")
write_csv(cleaned_data_summary, "Results/Data_merging_and_cleaning/cleaned_data_summary.csv")

