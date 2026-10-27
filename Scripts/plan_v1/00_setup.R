# =============================================================================
# Shared definitions for analysis plan v1 (analysis_plan_v1.md)
# =============================================================================
#
# Sourced by every plan_v1 script so that input paths, season dates, compound
# groups, and score construction are stated once. Run scripts from the
# repository root.

suppressPackageStartupMessages({
  library(tidyverse)
})

# -----------------------------------------------------------------------------
# Input files (plan §1; source-of-truth rule: DataHub March first, see D2)
# -----------------------------------------------------------------------------

input_files <- list(
  microresp = file.path(
    "DataHub March", "Soil_Phenotypes",
    "SoilMicrobiomeSFA_prosser-time-series_MicroResp_v3.csv"
  ),
  eea = file.path(
    "DataHub March", "Soil_Phenotypes",
    "SoilMicrobiomeSFA_prosser-time-series_EEA.csv"
  ),
  bgc = file.path(
    "DataHub March", "Soil_Biogeochemistry",
    "SoilMicrobiomeSFA_prosser-time-series_Soil_Biogeochemistry.csv"
  ),
  brite_d = file.path("Data", "brite_D_expression.csv"),
  chitin_kos = file.path(
    "Data", "merged_chitin_chitosan_methyl_oleate_gene_categories.csv"
  )
)

dir_results <- file.path("Results", "plan_v1")
dir_figures <- file.path("Figures", "plan_v1")
dir.create(dir_results, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_figures, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Season (plan §2.3)
# -----------------------------------------------------------------------------

dormant_dates <- as.Date(c(
  "2023-11-01", "2023-12-04", "2024-01-02", "2024-02-12", "2024-03-06"
))

study_dates <- as.Date(c(
  "2023-06-07", "2023-07-05", "2023-08-02", "2023-09-06", "2023-10-04",
  "2023-11-01", "2023-12-04", "2024-01-02", "2024-02-12", "2024-03-06",
  "2024-04-03", "2024-05-08", "2024-06-05"
))

assign_season <- function(date) {
  factor(if_else(date %in% dormant_dates, "Dormant", "Active"),
         levels = c("Active", "Dormant"))
}

# -----------------------------------------------------------------------------
# Compound groups (plan §2.4; fixed)
# -----------------------------------------------------------------------------

compound_groups <- tribble(
  ~Substrate, ~group,
  "D.Fructose", "Simple", "D.Xylose", "Simple", "D.Trehalose", "Simple",
  "D.Cellobiose", "Simple", "Citric_acid", "Simple", "Malic_acid", "Simple",
  "Oxalic_acid", "Simple", "Succinic_acid", "Simple", "Lactic_acid", "Simple",
  "Glutamic_acid", "Simple", "L.Arginine", "Simple", "L.Glutamine", "Simple",
  "Glycine", "Simple", "L.Alanine", "Simple", "L.Leucine", "Simple",
  "Cellulose", "Polymer", "Hemicellulose", "Polymer", "Lignin", "Polymer",
  "Chitin", "Polymer", "Chitosan", "Polymer",
  "D.Glucose", "Reference",
  "N.Acetyl.D.glucosamine", "Q2", "Glucosamine_hydrochloride", "Q2"
)

simple_compounds <- compound_groups$Substrate[compound_groups$group == "Simple"]
polymer_compounds <- compound_groups$Substrate[compound_groups$group == "Polymer"]

# -----------------------------------------------------------------------------
# Standardized scores (plan §2.6)
# -----------------------------------------------------------------------------
# z-score each compound across the samples supplied, then average z within
# the simple and polymer groups for each sample. `value` names the column to
# standardize, so the same function serves per-C, per-MBC, and the
# robustness variants.

make_scores <- function(df, value) {
  z <- df |>
    filter(Substrate %in% c(simple_compounds, polymer_compounds)) |>
    group_by(Substrate) |>
    mutate(z = (.data[[value]] - mean(.data[[value]], na.rm = TRUE)) /
             sd(.data[[value]], na.rm = TRUE)) |>
    ungroup() |>
    select(-any_of("group")) |>
    left_join(compound_groups, by = "Substrate")
  z |>
    group_by(Sample_cat, Collection_date, Plot, season, group) |>
    summarise(score = mean(z, na.rm = TRUE), n_compounds = sum(!is.na(z)),
              .groups = "drop") |>
    pivot_wider(names_from = group, values_from = c(score, n_compounds)) |>
    rename(simple_score = score_Simple, polymer_score = score_Polymer)
}

date_means <- function(scores) {
  scores |>
    group_by(Collection_date, season) |>
    summarise(n_plots = n(),
              across(c(simple_score, polymer_score), \(x) mean(x, na.rm = TRUE)),
              .groups = "drop") |>
    mutate(diff_score = polymer_score - simple_score)
}
