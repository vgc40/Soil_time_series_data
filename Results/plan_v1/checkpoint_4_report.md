# Checkpoint 4 report: Q2 analyses (plan §5)

Scripts: `Scripts/plan_v1/06_q2_analyses.R` (tests),
`07_q2_figures.R` (figures), `08_sensitivity_low_mbc.R` (D10). Gene sets:
approved D8 sets (`KO_gene_sets_APPROVED.csv`; 19 compound sets and 7 enzyme
sets used). Implementation choices not fixed by the plan: deviation log D9.

**Revised 2026-10-26 (D11):** the main potential for compound gene sets is
now per mg C (water-subtracted 5 h SIR, without dividing by MBC). The plan's
per-MBC version is robustness variant (a). This change was made after seeing
results and is recorded in the deviation log. Results are reported without
interpretation.

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

- Compounds (main, D11): water-subtracted 5 h SIR per mg C
  (`cp4_potential_samples.csv`).
- Compounds (robustness a): the same divided by MBC (§2.7)
  (`cp4_potential_samples_per_MBC.csv`). On the transcript dates, SM1625
  (2024-05-08, plot 46, MBC replaced with 0.07; D1) has per-MBC values about
  1,000 times the next-highest sample (D-fructose: 881.7 vs 0.94).
- Enzymes: measured activity (nmol per g dry soil per h), in both versions.

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
| **Main (potential per C; D11)** | 19 | 11 | 0.579 | 0.648 | 0.156 | 0.522 |
| a. Potential per MBC (plan original) | 19 | 13 | 0.684 | 0.167 | 0.335 | 0.161 |
| b. High-confidence KOs only (per C) | 17 | 8 | 0.471 | 1.000 | −0.007 | 0.981 |

Per-compound effects, main variant (dormant − active, z):

| Compound | Step | Expression | Potential (95% CI) | Concordant |
|---|---|---|---|---|
| Cellulose | extracellular | 0.40 | 1.05 (−4.61 to 6.71) | yes |
| Chitin | extracellular | −1.05 | 1.44 (0.65 to 2.24) | no |
| Hemicellulose | extracellular | 0.41 | 0.17 (−0.70 to 1.03) | yes |
| Citric acid | intracellular | −1.25 | −0.53 (−1.48 to 0.42) | yes |
| D-Cellobiose | intracellular | 0.86 | 1.04 (−5.23 to 7.31) | yes |
| D-Fructose | intracellular | −0.99 | 1.48 (0.19 to 2.77) | no |
| D-Trehalose | intracellular | 0.16 | 0.93 (−0.81 to 2.67) | yes |
| D-Xylose | intracellular | 0.37 | 1.30 (−1.45 to 4.05) | yes |
| Glucosamine | intracellular | 0.57 | 1.12 (−3.80 to 6.05) | yes |
| Glutamic acid | intracellular | −1.33 | −0.59 (−1.41 to 0.23) | yes |
| Glycine | intracellular | −1.26 | 1.51 (0.98 to 2.04) | no |
| L-Alanine | intracellular | −1.61 | 1.34 (−4.08 to 6.77) | no |
| L-Arginine | intracellular | −1.70 | 0.46 (−0.46 to 1.37) | no |
| L-Glutamine | intracellular | 0.58 | 1.51 (0.74 to 2.27) | yes |
| L-Leucine | intracellular | −1.07 | 1.59 (−1.94 to 5.11) | no |
| Lactic acid | intracellular | −0.40 | −0.21 (−1.33 to 0.90) | yes |
| Malic acid | intracellular | −1.00 | −0.22 (−2.41 to 1.98) | yes |
| NAG (compound) | intracellular | −0.01 | 0.10 (−2.89 to 3.10) | no |
| Oxalic acid | intracellular | −1.47 | 0.06 (−0.89 to 1.01) | no |

15 of 19 potential effects are positive; 12 of 19 expression effects are
negative. CIs use 2 dormant dates. In variant (a), 16 of 19 potential
effects lie between −0.219 and −0.212, because SM1625 dominates the per-MBC
z-scores.

### 5.3.4 H2c, lag advantage (26 gene sets, 7 date pairs, plot-matched; main variant)

Lag advantage = r_lag − r_same (Pearson, date means); 95% percentile
bootstrap CI, 2,000 resamples of date pairs. Enzyme sets are the same in all
variants.

| Set | Step | r_same | r_lag | Lag advantage | 95% CI |
|---|---|---|---|---|---|
| Ala (enzyme) | extracellular | −0.542 | 0.585 | 1.127 | 0.318 to 1.793 |
| AP (enzyme) | extracellular | −0.332 | 0.792 | 1.124 | 0.287 to 1.622 |
| Leu (enzyme) | extracellular | −0.576 | 0.119 | 0.695 | −0.189 to 1.574 |
| Chitin | extracellular | −0.702 | −0.117 | 0.585 | −0.955 to 1.336 |
| Cellulose | extracellular | 0.351 | 0.623 | 0.273 | −0.203 to 1.690 |
| Hemicellulose | extracellular | −0.139 | 0.120 | 0.259 | −0.739 to 1.488 |
| BG (enzyme) | extracellular | −0.422 | −0.198 | 0.224 | −0.905 to 1.395 |
| CB (enzyme) | extracellular | −0.222 | −0.261 | −0.039 | −0.464 to 0.598 |
| BX (enzyme) | extracellular | 0.329 | −0.060 | −0.389 | −1.369 to 1.005 |
| NAG (enzyme) | extracellular | 0.440 | −0.018 | −0.459 | −1.489 to 0.997 |
| L-Leucine | intracellular | −0.788 | −0.106 | 0.681 | −0.009 to 1.074 |
| D-Trehalose | intracellular | −0.087 | 0.373 | 0.461 | −0.399 to 1.107 |
| L-Alanine | intracellular | −0.872 | −0.544 | 0.328 | −0.584 to 0.877 |
| NAG (compound) | intracellular | 0.020 | 0.320 | 0.300 | −0.316 to 1.092 |
| Malic acid | intracellular | −0.050 | 0.196 | 0.246 | −1.110 to 0.918 |
| D-Fructose | intracellular | −0.690 | −0.490 | 0.200 | −0.456 to 0.927 |
| L-Arginine | intracellular | −0.778 | −0.642 | 0.135 | −1.009 to 0.981 |
| Glycine | intracellular | −0.712 | −0.585 | 0.126 | −1.275 to 0.844 |
| Oxalic acid | intracellular | 0.049 | 0.094 | 0.045 | −0.701 to 0.897 |
| D-Cellobiose | intracellular | 0.561 | 0.531 | −0.030 | −0.182 to 0.650 |
| Lactic acid | intracellular | 0.011 | −0.116 | −0.127 | −1.123 to 1.025 |
| L-Glutamine | intracellular | 0.313 | 0.185 | −0.128 | −0.808 to 0.786 |
| Citric acid | intracellular | 0.273 | −0.062 | −0.334 | −1.096 to 0.536 |
| D-Xylose | intracellular | 0.366 | −0.068 | −0.434 | −1.117 to 0.910 |
| Glutamic acid | intracellular | 0.416 | −0.025 | −0.441 | −1.465 to 0.476 |
| Glucosamine | intracellular | 0.483 | −0.025 | −0.508 | −1.515 to 0.175 |

Positive lag advantage: 7 of 10 extracellular, 9 of 16 intracellular. CI
excludes 0: Ala and AP (both above 0). Per set for all variants:
`cp4_h2c_lag_advantage_all_variants.csv`.

### 5.3.5 H2c, pathway step (permutation, 10,000 shuffles)

| Variant | n ext | n int | Mean ext | Mean int | Difference | p |
|---|---|---|---|---|---|---|
| **Main (potential per C; D11)** | 10 | 16 | 0.340 | 0.033 | 0.307 | 0.090 |
| a. Potential per MBC (plan original) | 10 | 16 | 0.326 | −0.053 | 0.379 | 0.079 |
| b. High-confidence KOs only (per C) | 10 | 14 | 0.302 | 0.074 | 0.229 | 0.173 |

## 5.4 Checks

### Plot-level lag models (main variant)

Per set: potential(t+1) ~ expression(t) + (1|date) vs potential(t) ~
expression(t) + (1|date); ML fits on the same rows (34 plot-dates for SIR
sets, 35 for enzyme sets); variables z-scored. 5 of 52 fits are singular
(date variance estimated at zero). Full table: `cp4_lag_models_compare.csv`.

| Step | Sets | Lag model lower AIC | Lag AIC lower by > 2 | Same-date AIC lower by > 2 | Lag marginal R² higher |
|---|---|---|---|---|---|
| Extracellular | 10 | 0 | 0 | 7 | 4 |
| Intracellular | 16 | 5 | 3 | 9 | 5 |

Marginal R² ranges: lag 0.0000 to 0.074 (highest: CB); same-date 0.0001 to
0.114 (highest: BX).

### Robustness

- (a) Potential per MBC and (b) high-confidence KOs only: rows in the H2b
  and pathway-step tables above.
- In (b), glucosamine and malic acid are left out (high-confidence sets zero
  in 27 and 35 samples; D9).

## Figures (`Figures/plan_v1/`; PDF, 300 dpi PNG, caption .txt)

| File | Content |
|---|---|
| `Fig4_H2b_quadrant` | H2b quadrant scatter; panel A per C (main, D11), panel B per MBC (robustness a) |
| `Fig5_H2c_lag_advantage` | Lag advantage per gene set with 95% CI, grouped by pathway step (main) |
| `Fig6_example_chitin_NAG_pathway` | Chitin set expression with chitin potential (per C and per MBC) and NAG activity; NAG set expression with NAG potential |
| `FigS5_H2a_expression_scores` | Expression group scores by season |
| `FigS6_bridge_Q1` | Tests 1 and 3 on 13 vs 7 dates |
| `FigS7_lag_models_plot_level` | Marginal R² of lag vs same-date models (main) |

Approved gene-set table: `Results/plan_v1/KO_gene_sets_APPROVED.csv`.

## Items for review

1. **D11** changes the plan's main Q2 potential after results were seen.
2. **D9 choices** made during implementation.
3. **O3** (mg C per well for succinic acid, chitin, methyl oleate) is still
   open. With per-C potential as the main version, z-scores and
   correlations still do not depend on a constant per-substrate amount; they
   would change only if chitin's amount varied by well.

---

## Addendum: per-MBC analyses without SM1625 and SM1539 (D10)

Since D11 these compare with robustness variant (a), the per-MBC version.
Potential per MBC with the two samples removed and z-scores recomputed.
Expression and enzyme potential are unchanged. Script
`08_sensitivity_low_mbc.R`.

### H2b

| Variant | n | Concordant | Share | Sign test p | Spearman rho | p |
|---|---|---|---|---|---|---|
| Per MBC, all samples | 19 | 13 | 0.684 | 0.167 | 0.335 | 0.161 |
| Per MBC, without the two | 19 | 12 | 0.632 | 0.359 | 0.286 | 0.234 |

Potential season effects without the two range from −0.739 (glutamic acid)
to 1.27 (chitin); 7 of 19 are negative. Per compound:
`sens_q2_h2b_compound_effects.csv`.

### H2c

| Variant | n ext | n int | Mean ext | Mean int | Difference | Permutation p |
|---|---|---|---|---|---|---|
| Per MBC, all samples | 10 | 16 | 0.326 | −0.053 | 0.379 | 0.079 |
| Per MBC, without the two | 10 | 16 | 0.296 | 0.107 | 0.189 | 0.421 |

- Only r_same changes for compound sets: SM1625 is on a transcript date
  (2024-05-08); SM1539 is plot 57, which has no transcripts, so it never
  entered the plot-matched lag pairs.
- Positive lag advantage: 7 of 10 extracellular, 10 of 16 intracellular.
- CI excludes 0: AP (1.16, 0.10 to 1.68), Ala (1.12, 0.27 to 1.80), lactic
  acid (0.99, 0.24 to 1.51) above; NAG compound (−1.15, −1.64 to −0.20)
  below.
- Per set: `sens_q2_h2c_lag_advantage.csv`; summary
  `sens_q2_h2c_step_permutation.csv`.

**Note (factual).** SM1420 (2023-12-04, MBC 7.78), which dominates the
per-MBC Q1 scores once the two samples are removed (Checkpoint 2 addendum),
does not enter H2b or H2c: 2023-12-04 is neither a transcript date nor a
next date in the lag pairs.
