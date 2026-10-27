# Checkpoint 4 report: Q2 analyses (plan §5)

Scripts: `Scripts/plan_v1/06_q2_analyses.R` (tests),
`07_q2_figures.R` (figures). Gene sets: approved D8 sets
(`KO_gene_sets_APPROVED.csv`; 19 compound sets and 7 enzyme sets used).
Implementation choices not fixed by the plan: deviation log D9. Results are
reported without interpretation.

Design: 7 transcript dates (2023-06-07, 07-05, 09-06, 10-04; 2024-01-02,
03-06, 05-08), 2 dormant and 5 active, 5 plots per date (35 samples; plot 57
has none). 100% irrigation.

## 5.1 Expression values

- **Units:** geTMM (gene-length-corrected TMM counts per million) summed per
  KO (deviation log, O2).
- **Zero replacement:** per sample, zeros set to half the smallest nonzero
  value. Half-minimum ranged 0.081 to 0.301 (median 0.198); zeros per sample
  4,220 to 6,478 of 9,022 KOs.
- **Centred log ratio:** log(set sum) − mean log over all 9,022 KOs in that
  sample. Files: `cp4_expression_clr_samples.csv`, `cp4_zero_replacement.csv`.

## 5.2 Potential values

- Compounds: water-subtracted 5 h per-C SIR divided by MBC (§2.7).
- Enzymes: measured activity (nmol per g dry soil per h).

**Note (factual).** On the transcript dates, SM1625 (2024-05-08, plot 46,
MBC replaced with 0.07; D1) has biomass-normalized values about 1,000 times
the next-highest sample (D-fructose: 881.7 vs 0.94). This sample has
transcripts and is in every per-MBC analysis below.

## 5.3 Tests

### 5.3.1 Bridge to Q1 (per-C scores, 7 transcript dates)

| Test | 7 dates: difference (95% CI) | p Welch | p Wilcoxon | 13 dates: difference | p Welch | Same direction |
|---|---|---|---|---|---|---|
| 1. Simple score | 0.676 (0.146 to 1.206) | 0.023 | 0.095 | 0.521 | 0.0031 | yes |
| 3. Polymer − simple | 0.307 (−1.077 to 1.692) | 0.371 | 0.381 | 0.108 | 0.580 | yes |

With 2 dormant and 5 active dates, the smallest possible Wilcoxon p is 0.095.

### 5.3.2 H2a, expression scores (descriptive; no p-values)

z-score of each set's CLR across 35 samples, averaged within group, date
means.

| Score | Sets | Mean dormant | Mean active | Dormant − active | 95% CI (Welch) |
|---|---|---|---|---|---|
| Simple | 14 | −0.515 | 0.206 | −0.721 | −1.090 to −0.352 |
| Polymer | 3 | −0.058 | 0.023 | −0.082 | −3.040 to 2.876 |
| Polymer − simple | — | 0.457 | −0.183 | 0.640 | −2.340 to 3.620 |

### 5.3.3 H2b, seasonal concordance (19 compounds, 7 transcript dates)

| Variant | n | Concordant | Share | Sign test p | Spearman rho | p |
|---|---|---|---|---|---|---|
| Main (potential per MBC) | 19 | 13 | 0.684 | 0.167 | 0.335 | 0.161 |
| a. Potential per C | 19 | 11 | 0.579 | 0.648 | 0.156 | 0.522 |
| b. High-confidence KOs only | 17 | 13 | 0.765 | 0.049 | 0.218 | 0.399 |

**Note (factual).** In the main variant, 17 of 19 compounds have a potential
season effect between −0.219 and −0.212; the exceptions are cellulose
(0.221), cellobiose (0.225) and chitin (0.039). The spread of the per-C
effects is −0.59 to 1.59. Per-compound effects for every variant:
`cp4_h2b_compound_effects_all_variants.csv`.

Expression season effects (z, dormant − active) range from −1.70
(L-arginine) to 0.86 (D-cellobiose). 12 of 19 are negative.

### 5.3.4 H2c, lag advantage (26 gene sets, 7 date pairs, plot-matched)

Lag advantage = r_lag − r_same (Pearson, date means); 95% percentile
bootstrap CI, 2,000 resamples of date pairs.

| Set | Step | r_same | r_lag | Lag advantage | 95% CI |
|---|---|---|---|---|---|
| Ala (enzyme) | extracellular | −0.542 | 0.585 | 1.127 | 0.318 to 1.793 |
| AP (enzyme) | extracellular | −0.332 | 0.792 | 1.124 | 0.287 to 1.622 |
| Cellulose | extracellular | −0.219 | 0.666 | 0.884 | −0.175 to 1.690 |
| Leu (enzyme) | extracellular | −0.576 | 0.119 | 0.695 | −0.189 to 1.571 |
| BG (enzyme) | extracellular | −0.422 | −0.198 | 0.224 | −0.905 to 1.389 |
| Chitin | extracellular | −0.136 | −0.044 | 0.091 | −0.604 to 1.120 |
| Hemicellulose | extracellular | −0.175 | −0.169 | 0.005 | −0.653 to 0.855 |
| CB (enzyme) | extracellular | −0.222 | −0.261 | −0.039 | −0.464 to 0.598 |
| BX (enzyme) | extracellular | 0.329 | −0.060 | −0.389 | −1.370 to 1.010 |
| NAG (enzyme) | extracellular | 0.440 | −0.019 | −0.459 | −1.490 to 0.997 |
| Lactic acid | intracellular | −0.140 | 0.807 | 0.946 | −0.418 to 1.710 |
| D-Fructose | intracellular | −0.167 | 0.249 | 0.416 | −0.573 to 1.260 |
| D-Xylose | intracellular | −0.303 | 0.081 | 0.384 | −0.560 to 1.310 |
| D-Cellobiose | intracellular | 0.308 | 0.540 | 0.232 | −0.732 to 0.752 |
| L-Glutamine | intracellular | −0.251 | −0.108 | 0.144 | −0.771 to 0.611 |
| Malic acid | intracellular | 0.330 | 0.437 | 0.107 | −1.490 to 1.020 |
| Oxalic acid | intracellular | 0.297 | 0.370 | 0.073 | −1.210 to 0.775 |
| Glutamic acid | intracellular | 0.235 | 0.190 | −0.045 | −1.350 to 0.674 |
| L-Leucine | intracellular | 0.108 | 0.053 | −0.054 | −0.483 to 0.963 |
| D-Trehalose | intracellular | 0.417 | 0.327 | −0.089 | −0.607 to 1.150 |
| Citric acid | intracellular | 0.220 | 0.051 | −0.169 | −1.300 to 0.607 |
| L-Alanine | intracellular | 0.126 | −0.052 | −0.178 | −0.766 to 1.180 |
| Glucosamine | intracellular | 0.017 | −0.318 | −0.334 | −1.620 to 0.106 |
| Glycine | intracellular | 0.255 | −0.365 | −0.620 | −1.300 to 0.609 |
| L-Arginine | intracellular | 0.176 | −0.516 | −0.692 | −1.290 to 0.971 |
| NAG (compound) | intracellular | 0.544 | −0.417 | −0.960 | −1.580 to −0.127 |

Positive lag advantage: 7 of 10 extracellular, 7 of 16 intracellular. CI
excludes 0: Ala and AP (above), NAG compound (below).

### 5.3.5 H2c, pathway step (permutation, 10,000 shuffles)

| Variant | n ext | n int | Mean ext | Mean int | Difference | p |
|---|---|---|---|---|---|---|
| Main | 10 | 16 | 0.326 | −0.053 | 0.379 | 0.077 |
| a. Potential per C | 10 | 16 | 0.340 | 0.033 | 0.307 | 0.093 |
| b. High-confidence KOs only | 10 | 14 | 0.278 | −0.083 | 0.361 | 0.131 |

## 5.4 Checks

### Plot-level lag models

Per set: potential(t+1) ~ expression(t) + (1|date) vs potential(t) ~
expression(t) + (1|date); ML fits on the same rows (34 plot-dates for SIR
sets, 35 for enzyme sets); variables z-scored. 17 of 52 fits are singular.
Full table: `cp4_lag_models_compare.csv`.

| Step | Sets | Lag model lower AIC | Lag AIC lower by > 2 | Same-date AIC lower by > 2 | Lag marginal R² higher |
|---|---|---|---|---|---|
| Extracellular | 10 | 2 | 1 | 6 | 3 |
| Intracellular | 16 | 14 | 7 | 1 | 11 |

Marginal R² ranges: lag 0.0000 to 0.227 (highest: lactic acid); same-date
0.0000 to 0.130 (highest: chitin).

### Robustness

- (a) Potential without biomass normalization and (b) high-confidence KOs
  only: rows in the H2b and pathway-step tables above. Per-set lag
  advantages for all variants: `cp4_h2c_lag_advantage_all_variants.csv`.
- In (b), glucosamine and malic acid are left out (high-confidence sets zero
  in 27 and 35 samples; D9). Lactic acid's lag advantage CI excludes 0 in
  (b) (1.103, 0.132 to 1.823).

## Figures (`Figures/plan_v1/`; PDF, 300 dpi PNG, caption .txt)

| File | Content |
|---|---|
| `Fig4_H2b_quadrant` | H2b quadrant scatter; panel A per MBC (plan), panel B per C (robustness a; D9) |
| `Fig5_H2c_lag_advantage` | Lag advantage per gene set with 95% CI, grouped by pathway step |
| `Fig6_example_chitin_NAG_pathway` | Chitin set expression with chitin potential and NAG activity; NAG set expression with NAG potential |
| `FigS5_H2a_expression_scores` | Expression group scores by season |
| `FigS6_bridge_Q1` | Tests 1 and 3 on 13 vs 7 dates |
| `FigS7_lag_models_plot_level` | Marginal R² of lag vs same-date models |

Approved gene-set table: `Results/plan_v1/KO_gene_sets_APPROVED.csv`.

## Items for review

1. **Per-MBC potential.** SM1625 dominates every compound's biomass-
   normalized z-scores on the transcript dates (5.2, 5.3.3 notes). The main
   H2b result and the compound-set r_same values carry this sample. Not run:
   a version without SM1625, or with per-C potential as the main variant.
2. **D9 choices** made during implementation.
3. **O3** (mg C per well for succinic acid, chitin, methyl oleate) is still
   open. In Q2 it affects only the per-C variant; z-scores and correlations
   do not depend on a constant per-substrate amount unless chitin's amount
   varied by well.
