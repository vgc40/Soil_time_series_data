# Analysis instructions: Prosser time series, Q1 and Q2

These instructions implement a locked analysis plan. Follow them in order. **Do not change any analysis decision after seeing results.** If something cannot be done as written (missing column, unexpected format, ambiguous rule), stop, record it in `Results/plan_v1/deviation_log.md`, and ask before proceeding.

---

## 0. Ground rules

- **Language:** R. Use tidyverse, ggplot2, lme4 + lmerTest, and base stats. Install packages if missing.
- **Branch:** create `analysis-plan-v1` from `manuscript-v4-analysis`. Work only on this branch.
- **Existing scripts:** ignore all existing scripts in `Scripts/` for analysis logic. The only exception: you may read them to copy the **ggplot theme and styling** (see section 7).
- **Do not modify any input data file.**
- **Outputs:**
  - Scripts: `Scripts/plan_v1/` (numbered in run order, e.g., `01_data_prep.R`)
  - Tables: `Results/plan_v1/`
  - Figures: `Figures/plan_v1/`
- **Checkpoints:** there are four (sections 2, 3, 4, 5). At each one: write a short report (`Results/plan_v1/checkpoint_N_report.md`), commit with a clear message, **push**, then stop and wait for approval before continuing.
- **Reporting:** report results factually (estimates, CIs, p-values, n). Do not write interpretations of what results mean biologically.
- **Save the plan:** before any analysis, commit this file as `analysis_plan_v1.md` at the repository root, so the plan is timestamped.

---

## 1. Inputs

**Source-of-truth rule:** use `DataHub March/` for every variable it contains. Use `Data/` only for variables not in DataHub March.

| Variable | Source |
|---|---|
| MicroResp SIR | `DataHub March/Soil_Phenotypes/SoilMicrobiomeSFA_prosser-time-series_MicroResp_v3.csv` |
| Enzyme activities (BG, BX, CB, AP, NAG, Leu, Ala), MBC, MBN | `DataHub March/Soil_Phenotypes/SoilMicrobiomeSFA_prosser-time-series_EEA.csv` |
| Biogeochemistry (GWC, pH, extractable C and N, etc.) | `DataHub March/Soil_Biogeochemistry/SoilMicrobiomeSFA_prosser-time-series_Soil_Biogeochemistry.csv` |
| Root biomass, field respiration, soil temperature | `Data/Data_merged_cleaned.csv` (only if not in DataHub March) |
| KO-level expression | `Data/brite_D_expression.csv` (columns `SM####_MT`) |
| Existing curated chitin/chitosan KOs | `Data/merged_chitin_chitosan_methyl_oleate_gene_categories.csv` |

Notes:
- The MicroResp file is not UTF-8. Read it with `locale(encoding = "latin1")`.
- Dates in DataHub files may be `M/D/YY`. Parse explicitly and convert to `YYYY-MM-DD`.
- Join tables on `Sample_cat`. Strip the `_MT` suffix from transcript column names to match.
- Read all DataHub READMEs and the data dictionary (`DataHub March/SoilSFA_A3-1_prosser_dictionary.xlsx`) before starting, and confirm the units for every variable you use.

---

## 2. Shared data preparation → Checkpoint 1

### 2.1 Filter
- 100% irrigation plots only (`Irrigation_trmt_pct == 100`), for everything except the 25% replication check (section 3.6).

### 2.2 Water subtraction and per-C conversion
For each sample and substrate, using the **per gram dry soil** values (`Normalization == "ngCO2-C_gdw_h"`):

```
induced_gdw = SIR_gdw(substrate) − SIR_gdw(Water, same Sample_cat, same Timepoint)
induced_C   = induced_gdw × gdw_well / mgC_per_well
```

- **Verify the conversion:** recompute the file's existing per-C values from the uncorrected per-gram values with the same formula (`SIR_gdw × gdw_well / mgC_per_well`). Report the agreement (it should match within rounding). If it does not, stop and report.
- **Keep negative values.** Do not truncate them to zero.
- Each value in the file is already a mean of technical well replicates. Do not average again.
- Do this for both `Timepoint` 5 and 24.

### 2.3 Season
- **Dormant:** 2023-11-01, 2023-12-04, 2024-01-02, 2024-02-12, 2024-03-06 (5 dates).
- **Active:** all other dates (8 dates, including 2024-04-03 and 2024-06-05).

### 2.4 Compound groups (fixed, do not change)
- **Simple (15):** D.Fructose, D.Xylose, D.Trehalose, D.Cellobiose, Citric_acid, Malic_acid, Oxalic_acid, Succinic_acid, Lactic_acid, Glutamic_acid, L.Arginine, L.Glutamine, Glycine, L.Alanine, L.Leucine
- **Polymer (5):** Cellulose, Hemicellulose, Lignin, Chitin, Chitosan
- **Reference:** D.Glucose
- **Held for Q2:** N.Acetyl.D.glucosamine, Glucosamine_hydrochloride
- **Excluded from group analyses:** all other substrates

### 2.5 Above-water flag
- For each compound, compute the share of samples (100%, 5 h) with `induced_gdw > 0`.
- **Flag** compounds below 75%. Keep all compounds in their groups for the main analyses. Flagged compounds are dropped only in a robustness check.

### 2.6 Standardized scores
- For each compound, z-score `induced_C` across all 100% samples at 5 h.
- **Simple score** and **polymer score** for each sample: mean z-score across the group's compounds.
- **Date means:** average sample scores to one value per date.

### 2.7 Biomass-normalized values
- `induced_C_per_MBC = induced_C / MBC` (MBC in µg C per g dry soil, from the EEA file, matched by `Sample_cat`).
- Recompute the z-scores and group scores from these values, exactly as in 2.6.

### 2.8 Root-dense core flag
- On root biomass, 100% plots, active-season samples only.
- Modified z-score: `M = 0.6745 × (x − median) / MAD`, where MAD is the median absolute deviation.
- **Flag** cores with `|M| > 3.5`. Expected flags are SM1078, SM1080, SM1621, and SM1126. Report what you actually get.
- Flagged cores are dropped only in the H1a root-biomass sensitivity analysis.

### 2.9 Missing and unmatched samples
- Use whatever samples exist. **Do not impute.**
- Date means use the plots available on that date (note February 2024 has 5 plots).
- Samples without a field-table match (e.g., SM1082, SM1508) are dropped only from analyses that need field variables. They stay in MicroResp-only analyses.

### Checkpoint 1 report
- Sample counts per date (MicroResp, EEA, field variables, transcripts).
- Per-C conversion check result.
- Water-subtraction summary per compound: median `induced_gdw` and share above zero.
- Flagged compounds (2.5) and flagged root cores (2.8).
- A list of all unmatched or missing samples.

**Commit, push, stop.**

---

## 3. Q1 analyses → Checkpoint 2

**Q1.** Does plant dormancy shift microbial substrate-use potential away from directly assimilable compounds and toward polymers that require extracellular breakdown?

- **H1a:** simple-compound responses lower in dormancy; within the active period, they increase with root biomass.
- **H1b:** polymer responses increase relative to simple responses in dormancy.
- **Alternative (priming/coupling):** the groups move together, with no relative shift.

All tests use the 100% plots, 5 h, water-subtracted per-C values, unless stated otherwise.

### 3.1 Hypothesis tests
1. **H1a, season effect.** Compare the 13 date means of the simple score: dormant (n = 5) vs. active (n = 8). Welch t-test and Wilcoxon rank-sum. Report the difference (dormant − active), its 95% CI, and both p-values.
2. **H1a, root biomass.** Active-season samples only. Spearman correlation between the sample-level simple score and root biomass, with and without the flagged root cores. Also compute the Spearman correlation for each of the 15 simple compounds individually (using `induced_C`), and report how many are positive.
3. **H1b, relative shift.** For each date, compute `polymer score − simple score`. Compare dormant vs. active date means of this difference, with the same tests as 3.1.1.
4. **Priming/coupling.** Pearson and Spearman correlations between the polymer and simple date means, across all 13 dates and within the dormant dates only.
5. **Compound-level effects.** For each simple and polymer compound, plus glucose: date means of the z-scores, then dormant − active difference with a 95% CI (Welch t-interval). Also report the dormant/active ratio of mean `induced_C` in real units.

### 3.2 Interpretation steps
6. **Biomass normalization.** Repeat tests 1, 3, and 4 using the biomass-normalized scores (2.7).
7. **Enzyme support for H1b.** For CB, BG, BX, and NAG: compare dormant vs. active date means (same tests as 3.1.1). Then compute Spearman correlations across the 13 date means between each enzyme and its matching substrate: CB and BG with cellulose, BX with hemicellulose, NAG with chitin.
8. **Existing capacity vs. growth.** First confirm from the README what the 24 h value represents (the rate over 5–24 h, or over 0–24 h), and report it. Then, for each group and season, compute the ratio of the later rate to the 5 h rate, and compare it between seasons.

### 3.3 Checks
- **Mixed models:** `z ~ season + (1|compound) + (1|Plot) + (1|date)` for the simple group (H1a), and `z ~ season * group + (1|compound) + (1|Plot) + (1|date)` for H1b (the interaction term is the test). Report whether conclusions agree with the date-level tests.
- **Circular-shift permutation:** order the 13 dates chronologically. Slide a window of 5 consecutive dates through every possible start position, without wrapping around the series ends. For each position, compute the difference in window mean minus outside mean, for the simple score and for the polymer − simple difference. p-value = the share of positions with an absolute difference at least as large as the observed (true dormant window). Report the number of positions, and note the minimum possible p-value.
- **Robustness variants** for tests 1 and 3: (a) uncorrected values (no water subtraction); (b) per gram soil instead of per mg C; (c) 24 h readings; (d) raw summed `induced_C` instead of z-scores; (e) leave out one compound at a time; (f) drop the flagged compounds from 2.5; (g) drop cellobiose.

### 3.4 Figures
**Main:**
1. Season comparison of the simple score: date means by season, with each season's mean and 95% CI. If the biomass-normalized result differs in direction or significance, add it as a second set of points.
2. Simple score vs. root biomass, active-season samples, with the flagged cores marked. Show trend lines with and without them.
3. Forest plot of the compound-level season effects (3.1.5), colored by group and sorted within group, with three summary rows at the bottom: simple-group mean, polymer-group mean, and the paired polymer − simple difference (3.1.3).

**Supplement:**
- Enzyme season comparisons and enzyme–substrate correlations (3.2.7).
- The later-to-early rate ratio by season and group (3.2.8).
- Time series of the simple and polymer scores over all 13 dates, with the dormant period shaded and a GAM smooth (mgcv, `s(day, k = 6)`, with plot as a random effect). For description only.
- Robustness table, with all variants and FDR-adjusted p-values.

### 3.5 Exploratory (report descriptively)
- Timing: date of peak z-score for cellulose and hemicellulose vs. chitin and chitosan.
- Date-level correlations of the simple and polymer scores with soil temperature, MBC, and extractable organic C.

### 3.6 Replication in the 25% plots
- Repeat tests 1 and 3 on the 25% irrigation plots, using the same pipeline. Report them separately.

### Checkpoint 2 report
- Results tables for 3.1–3.6, and all figures.

**Commit, push, stop.**

---

## 4. KO gene sets → Checkpoint 3 (requires approval)

Draft gene sets from the KEGG names in `brite_D_expression.csv`.

**Inclusion rule:**
- **Polymers:** the enzymes that break down that polymer (e.g., chitinases for chitin; endoglucanases and cellobiohydrolases for cellulose; xylanases for hemicellulose). Start from the existing curated chitin and chitosan list.
- **Simple compounds, NAG, and glucosamine:** a specific transporter plus the first committed catabolic step for that compound.
- **Exclude** shared central metabolism (glycolysis, TCA cycle, pentose phosphate pathway, general amino acid metabolism not specific to the compound).
- **Enzyme gene sets:** the KOs encoding each assayed enzyme (CB, BG, BX, NAG, Leu, Ala, AP).
- Label each gene set **extracellular** (polymer breakdown and all assayed enzymes) or **intracellular** (uptake and catabolism).
- If a compound has no specific KOs, mark it "no specific gene set." Do not substitute general pathways.

**Deliverable:** `Results/plan_v1/KO_gene_sets_DRAFT.csv`, with columns: compound or enzyme, KO, KEGG name, pathway step, reason for inclusion, and confidence (high/medium/low). Also report the summed expression of each gene set, so near-zero sets are visible.

**Commit, push, stop. Do not run any Q2 analysis until the gene sets are approved.**

---

## 5. Q2 analyses → Checkpoint 4

**Q2.** How do seasonal changes in field gene expression relate to changes in substrate-use potential?

- **H2a (descriptive):** in dormancy, expression of simple-compound genes decreases and polymer-degrading genes increase relative to them.
- **H2b, seasonal concordance:** expression and potential shift in the same direction between seasons.
- **H2c, monthly timing:** potential agrees better with the previous month's expression than with same-date expression, and more so for extracellular than intracellular steps.

### 5.1 Expression values
- First confirm the units of `brite_D_expression.csv` (e.g., TPM or raw counts), and report them.
- For each sample: replace zeros with half the smallest nonzero value in that sample. Compute each gene set's summed expression. Then take its centered log ratio: log(gene-set sum) minus that sample's mean log expression across all KOs.
- Average to date means (7 dates × 5 plots).

### 5.2 Potential values
- The biomass-normalized per-C SIR from 2.7, for each compound with an approved gene set.
- Measured enzyme activities, from the EEA file, for the enzyme gene sets.

### 5.3 Tests
1. **Bridge to Q1.** Repeat Q1 tests 1 and 3 on the 7 transcript dates only. Report whether the direction matches the 13-date result.
2. **H2a.** Build simple and polymer expression scores (z-score each gene set across samples, then average within group). Compute dormant − active differences with 95% CIs. Report as descriptive (2 dormant dates); no p-values.
3. **H2b.** For each compound with an approved gene set, calculate its season effect (dormant − active, z-units) in expression and in potential, both on the 7 transcript dates. Report the share of compounds in the concordant quadrants, a two-sided binomial sign test against 50%, and the Spearman correlation between the effects.
4. **H2c, lag.** Pair each transcript date with the next sampling date: Jun→Jul, Jul→Aug, Sep→Oct, Oct→Nov (2023), Jan→Feb, Mar→Apr, May→Jun (2024), matching plots. For each gene set, using date means:
   - r_same = Pearson correlation between expression and potential on the same date (7 pairs)
   - r_lag = Pearson correlation between expression at t and potential at the paired next date (7 pairs)
   - lag advantage = r_lag − r_same, with a 95% bootstrap CI (resample date pairs, 2,000 iterations)
   Do this for the SIR gene sets and the enzyme gene sets.
5. **H2c, pathway step.** Compare the mean lag advantage between extracellular and intracellular gene sets. Test with a permutation test (shuffle the pathway-step labels, 10,000 iterations).

### 5.4 Checks
- **Plot-level lag model:** `potential(t+1) ~ expression(t) + (1|date)` vs. `potential(t) ~ expression(t) + (1|date)`, using plot-matched pairs. Compare the fits (marginal R², AIC).
- **Robustness:** potential without biomass normalization; gene sets restricted to high-confidence KOs only.

### 5.5 Figures
**Main:**
1. **Quadrant scatter (H2b):** one point per compound, potential season effect (x) against expression season effect (y), colored by pathway step, with lines at zero.
2. **Lag advantage (H2c):** one point per gene set (compounds and enzymes), with 95% CIs, grouped as extracellular vs. intracellular.
3. **Example pathway:** over the 7 transcript dates, chitinase expression with chitin potential and NAG enzyme activity (all 13 dates), plus NAG uptake and catabolism expression with NAG potential.

**Supplement:**
- Expression group scores by season (H2a).
- The bridge check (5.3.1).
- The approved KO gene-set table.
- Plot-level lag model results (5.4).

### Checkpoint 4 report
- Results tables for 5.3–5.4, and all figures.

**Commit, push, stop.**

---

## 6. Statistical conventions

- Significance level: 0.05.
- 95% CIs for season effects: Welch t-intervals on date means.
- The unit of replication for season comparisons is the **date**, not the sample.
- No multiple-comparison correction in the forest plot (it is descriptive). Report FDR-adjusted (Benjamini–Hochberg) p-values in the supplementary tables.
- Always report n (dates and samples) alongside each result.
- Set a random seed (`set.seed(2026)`) for all bootstrap and permutation tests.

---

## 7. Figure style

- ggplot2, base theme `theme_bw()`. Match the visual style (fonts, sizes, grid lines, facet labels, legend placement) of the existing figures in `Figures/manuscript_final/`. You may read the existing scripts **only** to copy theme and styling code.
- Colors: Okabe–Ito palette. Simple compounds blue (`#0072B2`), polymers orange (`#E69F00`), dormant-period shading light gray. Glucose and other references in dark gray.
- Text: sans-serif, 8–10 pt.
- Export each figure as PDF and as PNG at 300 dpi, at journal widths: 85 mm (single column) or 180 mm (double column).
- Every figure gets a matching caption file (`.txt`) stating what is plotted, n, and what the error bars represent.
