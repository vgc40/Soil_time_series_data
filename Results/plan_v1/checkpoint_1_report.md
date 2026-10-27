# Checkpoint 1 report: shared data preparation (plan §2)

Script: `Scripts/plan_v1/01_data_prep.R` (shared definitions in
`00_setup.R`). Inputs: DataHub March tables only (see deviation log D2).
Decisions and open items: `deviation_log.md` (D1–D5, O1–O3).

## 1. Units confirmed (READMEs and data dictionary)

| Variable | Source | Unit |
|---|---|---|
| SIR, per gram | MicroResp v3, `Normalization = ngCO2-C_gdw_h` | ng CO2-C per g dry soil per h |
| SIR, per C | MicroResp v3, `Normalization = ngCO2-C_mgC_h` | ng CO2-C per mg added substrate C per h |
| `gdw_well` | MicroResp v3 | g dry soil per well |
| `mgC_per_well` | MicroResp v3 | mg substrate C added per well (stored to 2 decimals) |
| `Timepoint` | MicroResp v3 | end point of the measurement period, h (5 or 24); see O1 |
| BG, BX, CB, AP, NAG, Leu, Ala | EEA | nmol per g dry soil per h |
| MBC, MBN | EEA (identical to biogeochemistry file for all 156 shared samples) | µg C (N) per g dry soil |
| `Roots_gdw` | Biogeochemistry | README: g per g dry soil; values 0.58–53 suggest mg per g (D2) |
| `SoilTemp_mean` | Biogeochemistry | °C, 72 h mean at 18 cm before sampling |
| `SEOC_ug_gdw` | Biogeochemistry | µg C per g dry soil |
| `CO2_C_*_m2` | Biogeochemistry | not given in README or dictionary |

The data dictionary has no entry for root biomass or the field respiration
columns.

## 2. Sample counts per date

100% irrigation plots. Field = root biomass, soil temperature, GWC, and SEOC
all present. MBC is from the EEA file.

| Date | Season | MicroResp | EEA | MBC | Field | Transcripts |
|---|---|---|---|---|---|---|
| 2023-06-07 | Active | 6 | 6 | 6 | 6 | 5 |
| 2023-07-05 | Active | 6 | 6 | 6 | 6 | 5 |
| 2023-08-02 | Active | 6 | 6 | 6 | 6 | 0 |
| 2023-09-06 | Active | 6 | 6 | 6 | 6 | 5 |
| 2023-10-04 | Active | 6 | 6 | 6 | 6 | 5 |
| 2023-11-01 | Dormant | 6 | 6 | 6 | 6 | 0 |
| 2023-12-04 | Dormant | 6 | 6 | 6 | 6 | 0 |
| 2024-01-02 | Dormant | 6 | 6 | 6 | 6 | 5 |
| 2024-02-12 | Dormant | 5 | 6 | 6 | 6 | 0 |
| 2024-03-06 | Dormant | 6 | 6 | 6 | 6 | 5 |
| 2024-04-03 | Active | 6 | 6 | 6 | 6 | 0 |
| 2024-05-08 | Active | 6 | 6 | 6 | 6 | 5 |
| 2024-06-05 | Active | 6 | 6 | 6 | 6 | 0 |
| **Total** | 8 active, 5 dormant | **77** | **78** | **78** | **78** | **35** |

25% plots: 6 samples on every date for MicroResp, EEA, MBC, and field
variables (78 each); no transcripts. Full table: `cp1_sample_counts.csv`.

## 3. Per-C conversion check (§2.2)

Recomputed `SIR_gdw × gdw_well / mgC_per_well` against the file's per-C
values for 9,610 values (31 substrates, 155 samples, 2 timepoints; Water has
no per-C value).

- **28 of 31 substrates match** within rounding. Largest relative difference
  4.1% (oxalic acid), explained by `mgC_per_well` being stored to two
  decimals (oxalic acid: 0.10 in the column, 0.104 implied).
- **3 substrates do not match** (open item O3):

| Substrate | Group | `mgC_per_well` | Implied mg C, median (range) | Median relative difference |
|---|---|---|---|---|
| Succinic acid | Simple | 0.16 | 2.246 (2.243–2.249) | 13.04 |
| Chitin | Polymer | 1.30 | 1.430 (0.434–2.205) | 0.28 |
| Methyl oleate | Excluded | 2.84 | 2.841 (2.836–17.037) | 0.00 (up to 5.0) |

- Overall: correlation 0.911; 66.9% of values within 1%.
- Per-gram values are consistent: the implied `gdw_well` matches the column
  for every sample.
- **Interim decision (D5):** `induced_C` uses the `mgC_per_well` column for
  all substrates while Sheryl Bell confirms the correct amounts.

Files: `cp1_percarbon_check.csv`, `cp1_percarbon_check_by_substrate.csv`,
`cp1_percarbon_check_largest.csv`.

## 4. Water subtraction summary (100%, 5 h; n = 77 samples each)

`induced_gdw` = substrate minus Water, same sample and timepoint, in ng
CO2-C per g dry soil per h. Negative values kept.

| Substrate | Group | Median `induced_gdw` | Share > 0 | Flag (< 75%) |
|---|---|---|---|---|
| Citric_acid | Simple | 302.3 | 1.000 | |
| D.Cellobiose | Simple | 1.6 | 0.519 | **flag** |
| D.Fructose | Simple | 77.4 | 1.000 | |
| D.Trehalose | Simple | 50.2 | 1.000 | |
| D.Xylose | Simple | 76.0 | 0.974 | |
| Glutamic_acid | Simple | 135.5 | 1.000 | |
| Glycine | Simple | 34.4 | 0.922 | |
| L.Alanine | Simple | 30.1 | 0.922 | |
| L.Arginine | Simple | 21.5 | 0.792 | |
| L.Glutamine | Simple | 74.2 | 1.000 | |
| L.Leucine | Simple | 15.9 | 0.857 | |
| Lactic_acid | Simple | 272.4 | 0.987 | |
| Malic_acid | Simple | 246.8 | 1.000 | |
| Oxalic_acid | Simple | 135.6 | 1.000 | |
| Succinic_acid | Simple | 150.1 | 1.000 | |
| Cellulose | Polymer | 7.3 | 0.701 | **flag** |
| Chitin | Polymer | 9.0 | 0.714 | **flag** |
| Chitosan | Polymer | 58.9 | 0.857 | |
| Hemicellulose | Polymer | 87.2 | 1.000 | |
| Lignin | Polymer | 75.3 | 1.000 | |
| D.Glucose | Reference | 85.0 | 0.987 | |
| Glucosamine_hydrochloride | Q2 | 25.3 | 0.857 | |
| N.Acetyl.D.glucosamine | Q2 | 41.2 | 0.896 | |
| BSA | Excluded | 23.8 | 0.805 | |
| Diethylamine | Excluded | 91.0 | 1.000 | |
| GABA | Excluded | 18.1 | 0.779 | |
| Methyl_oleate | Excluded | 18.3 | 0.844 | |
| Phenylethylamine | Excluded | −13.3 | 0.156 | below 75% |
| Salicylic_acid | Excluded | 32.5 | 0.909 | |
| Tween40 | Excluded | 25.6 | 0.870 | |
| Vanillic_acid | Excluded | 58.3 | 1.000 | |

File: `cp1_water_subtraction_summary.csv`.

## 5. Flags

### Compounds below 75% above water (§2.5)
Within the analysis groups: **D.Cellobiose** (Simple, 51.9%), **Cellulose**
(Polymer, 70.1%), **Chitin** (Polymer, 71.4%). All remain in their groups for
the main analyses and are dropped only in robustness variant (f).

### Root-dense cores (§2.8)
100% plots, active season, n = 48 cores. Median `Roots_gdw` 5.21, MAD
(unscaled) 2.33. Six cores have |M| > 3.5:

| Sample | Date | Plot | `Roots_gdw` | M |
|---|---|---|---|---|
| SM1080 | 2023-07-05 | 40 | 53.02 | 13.87 |
| SM1078 | 2023-07-05 | 25 | 39.21 | 9.85 |
| SM1621 | 2024-05-08 | 1 | 34.13 | 8.39 |
| SM1126 | 2023-09-06 | 25 | 20.30 | 4.38 |
| SM1054 | 2023-06-07 | 25 | 18.67 | 3.91 |
| SM1150 | 2023-10-04 | 25 | 17.67 | 3.62 |

The plan expected four (SM1078, SM1080, SM1621, SM1126). The rule as written
also flags SM1054 and SM1150. All six are dropped in the H1a root-biomass
sensitivity analysis. File: `root_flags_active_100pct.csv`.

### MBC replacements (§2.7, D1)
Zeros replaced with 0.07 µg C per g (half the lowest nonzero value, 0.14).
Column `MBC_flag` marks them in every output table.

| Sample | Date | Plot | Treatment | MBC | Flag |
|---|---|---|---|---|---|
| SM1625 | 2024-05-08 | 46 | 100% | 0 → 0.07 | replaced_zero |
| SM1539 | 2024-02-12 | 57 | 100% | 0.14 (kept) | lowest_measured |
| SM1062 | 2023-06-07 | 33 | 25% | 0 → 0.07 | replaced_zero |
| SM1083 | 2023-07-05 | 8 | 25% | 0 → 0.07 | replaced_zero |
| SM1632 | 2024-05-08 | 53 | 25% | 0 → 0.07 | replaced_zero |

File: `cp1_mbc_replacements.csv`.

## 6. Unmatched and missing samples (§2.9)

- **SM1536** (100%, 2024-02-12, plot 27): no MicroResp measurement. Present
  in EEA and field tables. February 2024 date means use 5 plots.
- **Transcripts:** plot 57 has no transcript sample on any of the 7
  transcript dates (SM1058, SM1082, SM1130, SM1154, SM1507, SM1551, SM1626).
  Transcript date means use 5 plots.
- No MicroResp or transcript sample lacks a field-table match. SM1082 and
  SM1508, named in the plan, are both in the DataHub biogeochemistry file
  (D4).
- No missing values in EEA, MBC, root biomass, soil temperature, GWC, or
  SEOC on the 13 study dates. The EEA file records SM1502's date as `#N/A`;
  dates are taken from the biogeochemistry table by `Sample_cat`, so this has
  no effect.

Files: `cp1_unmatched_samples.csv`, `cp1_transcript_gaps.csv`.

## 7. Prepared tables for later checkpoints

| File | Contents |
|---|---|
| `sir_prepared_long.csv` | Every sample × substrate × timepoint, both treatments: per-gram SIR, Water, `induced_gdw`, `induced_C`, MBC used, `MBC_flag`, `induced_C_per_MBC`, season, group |
| `sample_scores_100pct_5h.csv` | Simple and polymer scores per sample, per C and per MBC |
| `date_scores_100pct_5h.csv` | Date means of both scores and polymer − simple, with n plots |

## 8. Open items carried forward

- **O1:** whether the 24 h reading is the rate over 5–24 h or 0–24 h. Needed
  for Checkpoint 2, test 8.
- **O2:** units of `brite_D_expression.csv`. Needed for Checkpoint 4.
- **O3:** correct mg C per well for succinic acid, chitin, and methyl oleate.
  Awaiting Sheryl Bell's answer.
