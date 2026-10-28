# Storyboard v6 figure log

Branch `storyboard-v6` (from `analysis-plan-v1`). Scripts: `Scripts/storyboard_v6/`
(run `00_run_all.R` from the repository root). All numbers shown in figures are
in `results_table.csv`; intermediate tables are in `data/`. This log reports
checks, exclusions, choices, and oddities only.

## Step 0 checks

### 1. Assay QC (Fig S8; `data/step0_assay_qc.csv`)
- **Incubation temperature: not available.** The MicroResp file has no
  temperature column. The DataHub README says plates were "incubated at room
  temperature" and the CO2 calculation uses 20 °C.
- **Run dates: not available.** The file has collection dates only. The README
  says MicroResp subsamples were run within one week of collection. QC is
  therefore reported by collection date, not by run date.
- **Substrate added:** `mgC_per_well` is a single value per substrate on every
  date (checked). Per gram of dry soil it varies only with soil per well
  (`gdw_well`, 0.23–0.44 g on the 100% plots). Date means of `gdw_well`:
  dormant − active = 0.013 g (95% CI −0.009 to 0.035; Welch p = 0.21; 13
  dates). Glucose as an example: 1.56–1.96 mg C g⁻¹ dry soil (date means). One
  2023-06-07 sample had 0.23 g per well (2.70 mg C g⁻¹ for glucose).
- **Water-control CO2 (5 h)**, date means 46.0–96.8 ng CO2-C g⁻¹ h⁻¹ (Fig S8A;
  season comparison in Fig S9A).

### 2. Microbial biomass (`data/step0_mbc_lowest.csv`)
- **Rule:** MBC = 0 (below the process blank, reported as zero in DataHub) or
  below 1% of the 100% median (median 165.9 µg C g⁻¹, so below 1.66) is set to
  **missing**, not zero.
- **Set to missing:** SM1625 (2024-05-08, plot 46, MBC 0) and SM1539
  (2024-02-12, plot 57, MBC 0.14).
- This replaces plan_v1 D1 (half-minimum replacement) for v6. Next lowest,
  kept: SM1420 (2023-12-04, plot 1, 7.78; 21× below the median).
- Affects Fig 1D (76 samples) and Fig S4 (75 samples).

### 3. Root outliers (`data/step0_root_flags.csv`)
- Modified z-score (0.6745 × (x − median) / MAD, |z| > 3.5) over all 78
  100% plot-dates, and again within each plot.
- **Flagged by both screens: SM1078, SM1080 (both 2023-07-05) and SM1621
  (2024-05-08)**, as expected.
- Flagged by one screen only, not used: SM1126 (whole-dataset only); SM1081,
  SM1398, SM1422 (within-plot only).
- Used only in Fig S1 (crosses; slopes with and without them). Root biomass is
  shown unfiltered in Fig 1F.

### Rule 4: negative primary values (5 h, 77 samples each; `data/step0_negative_values.csv`)
Kept as they are. Counts out of 77:
phenylethylamine 65, cellobiose 37, cellulose 23, chitin 22, GABA 17,
arginine 16, BSA 15, methyl oleate 12, chitosan 11, glucosamine 11, leucine 11,
Tween 40 10, N-acetylglucosamine 8, salicylic acid 7, alanine 6, glycine 6,
xylose 2, glucose 1, lactic acid 1. All other substrates 0.

## Exclusions
| What | IDs | Rule |
|---|---|---|
| MicroResp sample missing | SM1536 (2024-02-12, plot 27) | not measured; that date has 5 MicroResp plots |
| MBC set to missing | SM1625, SM1539 | Step 0.2 |
| Root points left out of the second S1 fit | SM1078, SM1080, SM1621 | Step 0.3 |
| Fold-over-water | none | no sample had water CO2 ≤ 0 |
| Gene sets not detected (Fig 3B grey, left out of class means) | Lignin, Chitosan | D8: set sum 0 in 9 or more of 35 samples |
| Classes left out of Fig 3 and Fig 4A (no gene sets) | Aromatic / phenolic acids; Amines; Protein; Fatty acid esters / surfactants | no v5 gene sets |
| Substrates in Fig 3 classes with no usable gene set | Glucose, Succinic acid, GABA (no set); Lignin, Chitosan (not detected) | v5 sets |

## Choices and deviations from the brief
1. **Season test:** plan_v1 `season_test()` (dormant − active of date means,
   Welch 95% CI, Welch p). The Wilcoxon p from the same function is in the
   `notes` column. `p` in the results table is the Welch p.
2. **Season shading:** bands change at the midpoints between the last date of
   one season and the first of the next (18 Oct 2023, 20 Mar 2024), so each
   sampling date sits inside its own band. Season membership is still
   Apr–Oct active / Nov–Mar dormant.
3. **Plant calendar dates are approximate:** harvests drawn at 15 Jul 2023 and
   15 Oct 2023, emergence at 1 Apr 2024. Exact dates were not supplied.
   On heatmaps (Fig 2C, S3) the markers are dashed lines between tiles; Fig 3B
   (7 date columns, no time axis) has season strips only.
4. **`substrate_classes.csv`:** the file row order is the heatmap row order.
   I have not seen Kirsten's grid. Order used: plant structural polymers →
   sugars (cellobiose, xylose first, to line up with cellulose, hemicellulose)
   → aromatic acids (vanillic first, for lignin) → chitinous polymers → amino
   sugars (NAG, glucosamine lined up with chitin, chitosan) → organic acids →
   amino acids → amines → protein → fatty acid esters. Both sugars and
   aromatic acids hold plant-polymer products, so only one of them can sit
   directly below the plant polymers. **The `eco_association` labels were
   assigned by me and need checking** (all rows `provisional = TRUE`).
   Amines have no label.
5. **Fig 1:** method reused from v4 Figure 1 (Welch on date means for
   temperature; `log(y) ~ season + (1|date) + (1|plot)` for the others,
   reported as dormant/active ratios). Field soil CO2 (`Resp_rate_mg_C_m2_h`)
   comes from `Data/Data_merged_cleaned.csv` because DataHub has only chamber
   concentrations. Temperature, GWC, roots, and MBC values in that table
   match DataHub exactly. The GWC model is singular (date variance only).
   Panel A shows field/enzyme, MicroResp, and transcript plots per date.
6. **Fig 2C/D:** class trajectories are the mean relative response of the
   substrates in the class, per sample, then averaged by date. Holm is over the
   4 tested classes only.
7. **Fig 2E:** per pair and date, the mean of (polymer z − product z) per
   sample (equal to the difference of relative responses). The per-pair rows
   show a Welch CI only. Pooled test: per date, the mean over the 5 pairs, then
   the season test.
8. **Fig 3:** expression is the plan_v1/v5 centred log ratio (log share of
   expression against the sample's geometric mean over all KOs; values match
   `cp4_expression_clr_samples.csv`). Each gene set is z-scored across the 35
   samples, class = mean of its sets per sample, then date means.
9. **Fig 4:** every series is z-scored at the sample level (MicroResp per
   substrate across 77 samples, level included, per mg C; enzymes across 78;
   gene sets across 35), then averaged by date. Fig 4A uses all MicroResp and
   enzyme plots on the 7 transcript dates, not only the 5 plots with
   transcripts. **Enzyme-to-class matching is my choice:** CB, BX → plant
   structural polymers; BG → simple sugars; NAG → chitinous polymers; Leu, Ala
   → amino acids; none for amino sugars and organic acids. In Fig 4C the
   "cellulase/cellobiohydrolase genes" are the v5 Cellulose set (endoglucanase
   + CBH); "β-glucosidase genes" are the v5 BG enzyme set.
10. **Fig S7:** Spearman ρ across the 7 transcript dates between class
    expression at t and the MicroResp class z (level included) at t or at the
    next sampling date. Leave-one-date-out ranges drop each transcript date.
    No p-values.
11. **Fig S9 gross version** uses `mgC_per_well` from the column, as for the
    primary version (open item O3, below), not the file's own per-C values.
12. **Fig S1:** the plot random effect was estimated as zero in both fits, so
    the models use date only (the same fallback as v4).
13. **Fig S6:** NMDS shown for display, as in v4 (allowed for S6).

## Things that look odd (no interpretation)
- **O3 is still open.** For succinic acid, chitin, and methyl oleate, the file's
  per-C values imply different mg C than the `mgC_per_well` column (plan_v1
  deviation log). Per-substrate z-scores do not depend on this constant, but
  any value in real per-C units does (only the unscaled responses in
  `data/microresp_long.csv`).
- Phenylethylamine is below the water control in 65 of 77 samples.
- Fig 3C: the plant structural polymer CI is −7.4 to 8.2 (two dormant dates;
  Welch df is small).
- Glucosamine's gene set is retained as in v5 but is zero in many samples
  (27 of 35 under the plan_v1 high-confidence variant).
- SM1420 (2023-12-04, plot 1) has field soil CO2 of 0.15 mg m⁻² h⁻¹, the
  lowest value; the other five plots that date are 8.5–14.4. The same sample
  has the third-lowest MBC (7.78).
- 2023-07-05 has the lowest MBC date mean; it includes the two root-flagged
  cores SM1078 and SM1080 (MBC 26 and 25).

## Results table
`results_table.csv` has 2,001 rows. `ci_low`/`ci_high` hold a 95% CI except
where `notes` says otherwise (Fig S7: leave-one-date-out min and max; Fig S8B:
min and max across samples). Date-mean rows have `n_dates = 1` and
`n_samples` = plots on that date.

## Round 2

Changes from `Figure_brief_v6_fixes_for_ClaudeCode.md`. Figures 1, 2, 4A, S1,
S2, S3, S5, S6, S7, S8, and S9 have no method or layout changes. Fig 2C has new
annotation text, and Fig 2 and Fig 3 were regenerated from the updated class
file.

1. **Fig S4 rebuilt at the date level.** For each substrate and date: the
   date-mean response per mg C added (all MicroResp plots on that date), divided
   by the date-mean MBC (plots with MBC after the Step 0 rule; SM1625 and SM1539
   missing). Each substrate is z-scored across the 13 dates, and the z-scores
   are averaged across the 31 substrates per date. Panel A has no plot points or
   CIs because the values exist only at the date level. Panel B uses the Fig 2B
   season test. Panel C repeats the test without SM1420 for:
   - the primary version, which reruns the Fig 2B pipeline with SM1420 removed
     before z-scoring (76 samples). The brief calls this "per gram soil";
     the primary response is per mg C added, so that is what is shown;
   - the date-level per-biomass version, which drops SM1420 from both date
     means.
   Results: per-biomass 0.36 (95% CI −0.22 to 0.94), p = 0.179; per-biomass
   without SM1420 0.28 (−0.36 to 0.92), p = 0.315; primary without SM1420 0.46
   (0.10 to 0.83), p = 0.018.
   - The primary-without-SM1420 test is also a row under Fig 2, panel B, marked
     as a sensitivity. It is not drawn in Fig 2.
   - The earlier sample-level per-MBC z-scores (`per_MBC_5h` in
     `data/microresp_z.csv`) are still written but no longer used in any figure.
2. **Phenylethylamine sensitivity:** the Fig 2B overall-level test with
   phenylethylamine removed (mean z of the other 30 substrates; the
   per-substrate z-scores are unchanged): 0.48 (0.13 to 0.84), p = 0.014. This
   is one row under Fig 2, panel B. There is no figure.
3. **Fig 3C:** the x-axis is limited to −3 to 3, and CIs that go past the limit
   end in arrows. The figure caption gives the full CIs:
   - plant structural polymers 0.40 (−7.4 to 8.2);
   - chitinous structural polymers −1.05 (−3.2 to 1.1).

   The chitinous CI also crosses −3, so it gets an arrow too.
4. **eco_association updated** in `Data/substrate_classes.csv` (still
   provisional):
   - Trehalose: "microbial / plant storage–stress sugar"
   - Lactic acid: "microbial metabolite; root-exudate-associated"
   - Salicylic acid: "plant signaling; root-exudate-associated"
   - Phenylethylamine and Diethylamine: "microbial / other"

   Fig 3B has no annotation column, so its content is unchanged.
5. **Plant calendar markers: NOT updated.** The round 2 brief has blank date
   fields for both harvests and emergence. The approximate dates (15 Jul 2023,
   15 Oct 2023, 1 Apr 2024) and the "approximate" note above remain until the
   dates are supplied.
6. **Fig 4B/C:**
   - x-axis ticks every second sampling date (no overlapping labels);
   - lines between dates more than 45 days apart are dashed. This affects only
     the transcript series (Oct→Jan and Mar→May). MicroResp and enzyme series
     have at most 41 days between dates (Jan→Feb).

The results table now has 2,005 rows.
