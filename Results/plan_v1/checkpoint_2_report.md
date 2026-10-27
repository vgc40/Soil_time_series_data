# Checkpoint 2 report: Q1 analyses (plan §3)

Scripts: `Scripts/plan_v1/02_q1_tests.R` (tests 1–8),
`03_q1_checks.R` (checks, exploratory, 25% replication), `04_q1_figures.R`.
All tests use the 100% plots, 5 h, water-subtracted per-C values unless
stated. Season comparisons use date means: n = 13 dates (8 active, 5
dormant), 77 samples. Differences are dormant − active; CIs are Welch 95%
t-intervals. Decisions affecting this checkpoint: D1, D5, D6, D7 in
`deviation_log.md`. Results are reported without interpretation.

## 3.1 Hypothesis tests

### Test 1. H1a, season effect on the simple score

| Normalization | Mean dormant | Mean active | Difference | 95% CI | Welch t (df) | p Welch | W | p Wilcoxon |
|---|---|---|---|---|---|---|---|---|
| Per C | 0.325 | −0.196 | **0.521** | 0.218 to 0.825 | 3.79 (10.7) | 0.0031 | 39 | 0.0031 |

### Test 2. H1a, simple score vs root biomass (active season)

Spearman correlation, sample level.

| Set | n samples | rho | p |
|---|---|---|---|
| All active cores | 48 | −0.007 | 0.960 |
| 6 flagged root cores removed | 42 | 0.087 | 0.584 |

Per compound (`induced_C` vs `Roots_gdw`, 15 simple compounds):

| Set | Positive rho | Positive with p < 0.05 |
|---|---|---|
| All active cores | 8 of 15 | 1 |
| Flagged cores removed | 9 of 15 | 1 |

Per-compound values: `cp2_t2_root_biomass_compounds.csv` (with FDR).

### Test 3. H1b, polymer − simple difference

| Normalization | Mean dormant | Mean active | Difference | 95% CI | Welch t (df) | p Welch | W | p Wilcoxon |
|---|---|---|---|---|---|---|---|---|
| Per C | 0.071 | −0.038 | **0.108** | −0.337 to 0.553 | 0.58 (6.6) | 0.580 | 26 | 0.435 |

Polymer score alone (for the forest plot): difference 0.630 (95% CI 0.039 to
1.221), p Welch 0.041, p Wilcoxon 0.045.

### Test 4. Coupling of polymer and simple date means

| Dates | n | Pearson r | p | Spearman rho | p |
|---|---|---|---|---|---|
| All 13 | 13 | 0.788 | 0.0014 | 0.742 | 0.0053 |
| Dormant only | 5 | 0.737 | 0.156 | 0.700 | 0.233 |

### Test 5. Compound-level season effects

z-score date means; dormant − active with Welch 95% CI. Ratio: dormant/active
mean of date-mean `induced_C` (real units). p-values unadjusted; FDR-adjusted
values in `cp2_t5_compound_effects.csv`.

| Compound | Group | Difference (z) | 95% CI | p Welch | p Wilcoxon | Ratio induced_C |
|---|---|---|---|---|---|---|
| Glycine | Simple | 1.15 | 0.59 to 1.71 | 0.0010 | 0.0016 | 2.02 |
| D.Fructose | Simple | 0.89 | 0.09 to 1.68 | 0.034 | 0.065 | 1.54 |
| D.Xylose | Simple | 0.87 | 0.09 to 1.64 | 0.032 | 0.045 | 1.61 |
| L.Leucine | Simple | 0.81 | −0.23 to 1.84 | 0.108 | 0.171 | 2.22 |
| L.Glutamine | Simple | 0.77 | −0.11 to 1.64 | 0.080 | 0.127 | 1.53 |
| Oxalic_acid | Simple | 0.70 | 0.02 to 1.38 | 0.046 | 0.093 | 1.55 |
| D.Cellobiose | Simple | 0.61 | −0.53 to 1.75 | 0.221 | 0.284 | −45.6 |
| D.Trehalose | Simple | 0.43 | −0.41 to 1.26 | 0.276 | 0.354 | 1.32 |
| L.Alanine | Simple | 0.38 | −0.86 to 1.61 | 0.469 | 0.622 | 1.38 |
| Lactic_acid | Simple | 0.34 | −0.49 to 1.17 | 0.387 | 0.622 | 1.25 |
| Succinic_acid | Simple | 0.30 | −0.66 to 1.26 | 0.490 | 0.435 | 1.25 |
| Citric_acid | Simple | 0.24 | −0.52 to 1.00 | 0.504 | 0.354 | 1.15 |
| Malic_acid | Simple | 0.14 | −0.42 to 0.69 | 0.603 | 0.833 | 1.09 |
| L.Arginine | Simple | 0.13 | −0.34 to 0.61 | 0.543 | 0.622 | 1.17 |
| Glutamic_acid | Simple | 0.09 | −0.67 to 0.84 | 0.801 | 0.524 | 1.06 |
| Chitin | Polymer | 1.13 | 0.24 to 2.03 | 0.022 | 0.0031 | 6.96 |
| Lignin | Polymer | 0.90 | 0.17 to 1.62 | 0.022 | 0.019 | 1.55 |
| Cellulose | Polymer | 0.81 | −0.22 to 1.84 | 0.100 | 0.065 | 4.57 |
| Chitosan | Polymer | 0.37 | −1.19 to 1.93 | 0.578 | 0.833 | 1.44 |
| Hemicellulose | Polymer | −0.06 | −0.68 to 0.55 | 0.827 | 0.833 | 0.96 |
| D.Glucose | Reference | 0.84 | 0.15 to 1.52 | 0.021 | 0.045 | 1.54 |

19 of 20 group compounds have a positive difference; hemicellulose is the
exception.
The cellobiose ratio is negative because its active-season mean `induced_C`
is below zero.

## 3.2 Interpretation steps

### Test 6. Biomass-normalized scores (tests 1, 3, 4 repeated)

| Test | Difference | 95% CI | p Welch | p Wilcoxon |
|---|---|---|---|---|
| 1. Simple score | 0.048 | −0.484 to 0.580 | 0.841 | 0.435 |
| 3. Polymer − simple | 0.166 | −0.096 to 0.428 | 0.190 | 0.222 |
| Polymer score | 0.214 | −0.487 to 0.916 | 0.451 | 0.222 |

| Coupling (test 4) | n | Pearson r | p | Spearman rho | p |
|---|---|---|---|---|---|
| All 13 dates | 13 | 0.795 | 0.0012 | 0.643 | 0.021 |
| Dormant only | 5 | 1.000 | < 0.001 | 0.900 | 0.083 |

**Note on these values (factual).** The two samples with replaced or lowest
MBC (SM1625, 0.07; SM1539, 0.14; D1) have biomass-normalized simple scores
of 6.66 and 4.49 and polymer scores of 2.09 and 6.27. The other 75 samples
have simple scores between −0.154 and −0.019 and polymer scores between
−0.116 and 0.003 (the highest is SM1420, MBC 7.78). The date means for
2024-02-12 (0.78) and 2024-05-08 (0.98) carry the two flagged samples; the
other 11 date means lie between −0.152 and −0.129. The dormant-only Pearson
r = 1.000 reflects the 2024-02-12 date. See Figure S4.

### Test 7. Enzyme support for H1b

Season comparison of date means (n = 13 dates, 78 samples):

| Enzyme | Mean dormant | Mean active | Ratio | Difference | 95% CI | p Welch | p Wilcoxon |
|---|---|---|---|---|---|---|---|
| BG | 319.5 | 329.6 | 0.97 | −10.1 | −124.7 to 104.4 | 0.848 | 1.000 |
| BX | 36.8 | 34.5 | 1.07 | 2.4 | −8.2 to 12.9 | 0.630 | 0.524 |
| CB | 15.6 | 14.8 | 1.06 | 0.8 | −6.8 to 8.4 | 0.808 | 0.833 |
| NAG | 41.1 | 35.7 | 1.15 | 5.4 | −8.0 to 18.7 | 0.371 | 0.524 |

Units nmol per g dry soil per h. FDR-adjusted p: all ≥ 0.85.

Spearman correlation across 13 date means, enzyme vs matching substrate z:

| Enzyme | Substrate | rho | p | p FDR |
|---|---|---|---|---|
| CB | Cellulose | 0.236 | 0.437 | 0.583 |
| BG | Cellulose | 0.258 | 0.394 | 0.583 |
| BX | Hemicellulose | 0.060 | 0.845 | 0.845 |
| NAG | Chitin | 0.319 | 0.289 | 0.583 |

### Test 8. Later (24 h) to early (5 h) rate ratio

**What the 24 h value represents:** the rate over 5–24 h. The indicator
plate was replaced at 5 h and the new plate read at 24 h (earlier manuscript
methods; see deviation log, resolution of O1). Whether the 5–24 h CO2 was
divided by 19 h or 24 h is not stated; that constant does not affect the
season comparison. Ratio computed per date and group as mean `induced_C` at
24 h over mean at 5 h (D6).

| Group | Mean ratio dormant | Mean ratio active | Difference | 95% CI | p Welch | p Wilcoxon |
|---|---|---|---|---|---|---|
| Simple | 0.906 | 1.080 | −0.174 | −0.754 to 0.406 | 0.522 | 0.943 |
| Polymer | 1.100 | 1.229 | −0.129 | −0.664 to 0.407 | 0.598 | 0.724 |

## 3.3 Checks

### Mixed models

Satterthwaite df (lmerTest). n = 1,155 observations (H1a: 15 compounds × 77
samples) and 1,540 (H1b: 20 × 77). Both fits are singular: the compound
variance is zero because each compound's z-scores are centred.

| Model | Term | Estimate | SE | df | p | 95% CI (Wald) |
|---|---|---|---|---|---|---|
| H1a: simple z ~ season | season (dormant) | 0.519 | 0.208 | 3.2 | 0.083 | 0.110 to 0.927 |
| H1b: z ~ season × group | season (dormant) | 0.522 | 0.213 | 3.0 | 0.091 | 0.105 to 0.939 |
| H1b | group (polymer) | −0.038 | 0.068 | 1520 | 0.579 | −0.171 to 0.095 |
| H1b | season × group (interaction) | 0.100 | 0.111 | 1520 | 0.366 | −0.117 to 0.316 |

Variance components (SD): date 0.35, plot 0.21 (H1a) / 0.20 (H1b),
residual 0.91.

**Agreement with date-level tests:**
- H1a: the direction and size agree (0.519 vs 0.521). Significance does not:
  p = 0.083 in the mixed model vs p = 0.0031 for the date-level Welch test.
- H1b: both not significant (interaction p = 0.366; date-level p = 0.580).

### Sliding-window permutation (no wrap-around)

9 window positions; minimum possible p = 1/9 = 0.111. The true dormant
window starts at 2023-11-01.

| Score | Observed window − outside | Positions with \|diff\| ≥ observed | p |
|---|---|---|---|
| Simple score | 0.521 | 1 of 9 | 0.111 |
| Polymer − simple | 0.108 | 4 of 9 | 0.444 |

All window differences: `cp2_window_permutation.csv`.

### Robustness variants (tests 1 and 3)

FDR (Benjamini–Hochberg) within each test across all 27 variants. Full table
including each leave-one-out: `cp2_robustness.csv`.

| Variant | Test 1 diff (95% CI) | p Welch (FDR) | Test 3 diff (95% CI) | p Welch (FDR) |
|---|---|---|---|---|
| Main | 0.521 (0.218 to 0.825) | 0.0031 (0.0054) | 0.108 (−0.337 to 0.553) | 0.580 (0.809) |
| a. Uncorrected (no water subtraction) | 0.534 (0.019 to 1.048) | 0.043 (0.045) | 0.089 (−0.338 to 0.516) | 0.616 (0.809) |
| b. Per gram dry soil | 0.470 (0.153 to 0.786) | 0.0075 (0.0085) | 0.100 (−0.345 to 0.545) | 0.608 (0.809) |
| c. 24 h readings | 0.410 (−0.338 to 1.158) | 0.227 (0.227) | 0.075 (−0.395 to 0.546) | 0.731 (0.822) |
| d. Raw summed induced_C (ng CO2-C/mg C/h) | 786 (62 to 1510) | 0.036 (0.039) | −763 (−1484 to −42) | 0.040 (0.809) |
| e. Leave one out, range over 20 | 0.476 to 0.552 | max 0.0074 | −0.018 to 0.281 | min 0.221 |
| f. Drop flagged (cellobiose, cellulose, chitin) | 0.515 (0.228 to 0.802) | 0.0023 (0.0054) | −0.114 (−0.786 to 0.558) | 0.686 (0.809) |
| g. Drop cellobiose | 0.515 (0.228 to 0.802) | 0.0023 (0.0054) | 0.115 (−0.361 to 0.591) | 0.584 (0.809) |

Leave-one-out extremes for test 3: dropping hemicellulose gives 0.281
(p = 0.221); dropping chitin gives −0.018 (p = 0.940). Variant (d) depends on
the `mgC_per_well` values under open item O3.

## 3.5 Exploratory (descriptive)

Date of the highest date-mean z-score:

| Compound | Peak date | Season | Peak z |
|---|---|---|---|
| Hemicellulose | 2023-09-06 | Active | 1.15 |
| Chitin | 2023-12-04 | Dormant | 1.94 |
| Cellulose | 2024-01-02 | Dormant | 1.46 |
| Chitosan | 2024-03-06 | Dormant | 1.71 |

Date-level correlations (n = 13 dates, 100% plots):

| Score | Variable | Pearson r | p | Spearman rho | p |
|---|---|---|---|---|---|
| Simple | Soil temperature | −0.711 | 0.0064 | −0.703 | 0.0096 |
| Simple | MBC | 0.578 | 0.038 | 0.302 | 0.315 |
| Simple | SEOC | −0.608 | 0.028 | −0.154 | 0.617 |
| Polymer | Soil temperature | −0.624 | 0.023 | −0.533 | 0.064 |
| Polymer | MBC | 0.370 | 0.213 | 0.181 | 0.554 |
| Polymer | SEOC | −0.378 | 0.203 | −0.066 | 0.835 |

GAM (Figure S3, description only): s(day) edf 3.15, p = 0.0023, deviance
explained 32% (simple); edf 3.52, p < 0.001, 31% (polymer).

## 3.6 Replication in the 25% plots

Same pipeline, z-scores computed within the 25% samples (n = 78 samples, 13
dates, 6 plots per date).

| Test | Mean dormant | Mean active | Difference | 95% CI | p Welch | p Wilcoxon |
|---|---|---|---|---|---|---|
| 1. Simple score | 0.179 | −0.112 | 0.290 | −0.327 to 0.907 | 0.291 | 0.093 |
| 3. Polymer − simple | 0.069 | −0.043 | 0.113 | −0.400 to 0.625 | 0.616 | 0.354 |

## Figures (`Figures/plan_v1/`, each with PDF, 300 dpi PNG, and caption .txt)

| File | Content |
|---|---|
| `Fig1_season_simple_score` | Test 1. Includes the per-MBC panel because its result differs in significance (plan rule); dates containing D1-flagged samples are open circles |
| `Fig2_simple_score_root_biomass` | Test 2, flagged cores marked, trend lines with and without them |
| `Fig3_forest_compound_effects` | Test 5 with summary rows for tests 1 and 3 |
| `FigS1_enzymes` | Test 7 |
| `FigS2_rate_ratio_24h_5h` | Test 8 |
| `FigS3_score_time_series` | Time series with GAM smooths, dormant period shaded |
| `FigS4_mbc_normalized_scores_flags` | Extra diagnostic (D7): per-MBC sample scores with flagged samples labelled |

The robustness table is `Results/plan_v1/cp2_robustness.csv`.

## Items for review

1. **Biomass-normalized results (test 6, Figure 1 right panel)** are driven
   by two samples (see note under test 6).
2. **Test 8 computation (D6)** was chosen during implementation and needs
   your approval.
3. **Open items:** O1 (24 h meaning), O2 (transcript units), O3 (mg C for
   succinic acid, chitin, methyl oleate).
