# Deviation log: analysis plan v1

Each entry records a point where the plan could not be followed exactly as
written, the decision taken, who approved it, and the date. Entries are added
in the order they arise and are never edited after the fact.

---

## D1. MBC values of zero (plan §2.7)

- **Issue:** `induced_C_per_MBC = induced_C / MBC` is undefined when MBC = 0.
  Four samples have MBC = 0 in the EEA file. The DataHub README states that
  values below the process blank "were converted to zeros" and gives no
  numeric detection limit.
  - 100%: SM1625 (2024-05-08, plot 46)
  - 25%: SM1062 (2023-06-07), SM1083 (2023-07-05), SM1632 (2024-05-08)
- **Decision:** replace each zero with half the lowest nonzero MBC across
  both irrigation treatments. The lowest nonzero value is SM1539
  (2024-02-12, 100%) at 0.14 µg C/g dry soil, so zeros become 0.07.
  SM1539 itself is a measured value and is kept at 0.14.
- **Flag:** every sample carries an `MBC_flag` column:
  `replaced_zero` (MBC was 0, set to 0.07), `lowest_measured` (SM1539),
  or `none`. Figures that use biomass-normalized values mark flagged samples
  with a distinct symbol.
- **Note:** both flagged 100% samples have MBC about 1,000 times below the
  100% median (166 µg C/g), so their biomass-normalized values will be
  extreme. No sensitivity analysis was added; that was not requested.
- **Approved by:** VGC, 2026-10-26.

## D2. Source for root biomass, soil temperature, and field respiration (plan §1)

- **Issue:** the input table lists `Data/Data_merged_cleaned.csv` for these
  variables "only if not in DataHub March". All three are in the DataHub
  biogeochemistry file (`Roots_gdw`, `Root_mass`, `SoilTemp_mean`,
  `CO2_C_*_m2`).
- **Decision:** following the source-of-truth rule, take them from DataHub
  March. Root biomass is `Roots_gdw` (root mass per unit soil), not
  `Root_mass` (g per sample).
- **Unit note:** the README gives `Roots_gdw` as g per g dry soil, but values
  range 0.58 to 53, which fits mg per g better. The variable is used only in
  rank-based statistics (Spearman, modified z on the median), so the unit
  does not change any result. Flagged for confirmation.
- **Consequence:** values are raw DataHub values, without the outlier
  cleaning applied by the earlier pipeline. For example, SM1082 SEOC (set to
  NA in `Data_merged_cleaned.csv`) is used as reported.
- **Approved by:** VGC, 2026-10-26.

## D3. Figure style source (plan §7)

- **Issue:** `Figures/manuscript_final/` does not exist.
- **Decision:** match `Figures/manuscript_v4/` and the theme in
  `Scripts/helpers.R` (`theme_manuscript`).
- **Approved by:** VGC, 2026-10-26.

## D4. Samples listed as unmatched in the plan (plan §2.9)

- **Issue:** the plan names SM1082 and SM1508 as samples without a field
  match. Both are present in the DataHub biogeochemistry file (SM1082: 100%,
  2023-07-05, plot 57; SM1508: 25%, 2024-01-02, plot 8). Their absence was
  specific to the earlier cleaned table.
- **Decision:** keep both in all analyses. The actual list of unmatched
  samples is reported in the Checkpoint 1 report.
- **Approved by:** VGC, 2026-10-26.

## D5. Interim mg C values while O3 is open (plan §2.2)

- **Decision:** compute `induced_C` with the `mgC_per_well` column for all
  substrates, as the plan's formula states, and continue while Sheryl Bell
  confirms the correct mg C for succinic acid, chitin, and methyl oleate.
  If she reports different values, `01_data_prep.R` is rerun and every
  downstream result is regenerated.
- **Scope of the effect:** `mgC_per_well` is one constant per substrate, so
  per-compound z-scores (§2.6, §2.7) and dormant/active ratios of means
  (§3.1.5) do not depend on it. Only values in real per-C units are affected
  (robustness variant (d), summed raw `induced_C`). If chitin amounts truly
  varied by well, chitin z-scores would change.
- **Approved by:** VGC, 2026-10-26.

## D6. How the later-to-early rate ratio is computed (plan §3.2.8)

- **Issue:** the plan asks for "the ratio of the later rate to the 5 h
  rate" for each group and season, without saying at which level it is
  formed. Water-subtracted rates can be zero or negative for single
  compound-sample pairs (cellobiose median 1.6), so ratios formed per sample
  and compound are unstable.
- **Decision:** for each date and group, divide the mean `induced_C` at 24 h
  (over all plots and group compounds) by the same mean at 5 h. The 13 date
  ratios are then compared between seasons with the tests of §3.1.1.
- **Note:** O1 (whether 24 h is a 5–24 h or 0–24 h rate) is still open; the
  ratio uses the 24 h values as reported.
- **Approved by:** taken by Claude while implementing; for VGC review at
  Checkpoint 2.

## D7. Extra diagnostic figure (not in plan)

- **Addition:** Figure S4 shows the biomass-normalized sample scores with
  the D1-flagged samples marked, as VGC asked for a way to see them in
  figures. It adds no test.
- **Approved by:** requested by VGC, 2026-10-26.

## Resolutions of open items (2026-10-26)

Source folder: `OneDrive - PNNL/Soil Microbiome SFA - A3.1 Rhizosphere and
Moisture Effects/Prosser_Monthly/Manuscripts/A3-1-1_Prosser_time_biogeochem/`
(pointed to by VGC).

- **O1 resolved: the 24 h value is the rate over 5–24 h.** The MicroResp
  methods in `Archive/SFA3-1-1_StoichiometryManuscript.docx` state that the
  indicator plate was removed after 5 h and replaced with a new one, which
  was read again at 24 h "to calculate the rates of respiration early and
  late in the incubation". The text does not say whether the 5–24 h CO2 was
  divided by 19 h or 24 h. Either way the factor is the same constant for
  every sample, so the season comparison in test 8 (D6) is unaffected; only
  the absolute ratio would shift.
- **O2 resolved (strong evidence, generating code not found): expression
  values are geTMM.** `Transcriptomic_abundance/geTMM_normalization.R`
  converts featureCounts (counts < 5 set to 0) to reads per kilobase, then
  applies edgeR TMM and `cpm()`: gene-length-corrected TMM counts per million
  (geTMM). Per-sample totals of `brite_D_expression.csv` are 1.13 to 1.79
  times the gene-level geTMM totals (median 1.57; Pearson r = 0.98 across 35
  samples), consistent with gene-level geTMM summed per KO, with genes that
  carry more than one KO counted once per KO. The script that rolled genes up
  to KOs is not in either location. The centred log ratio in plan §5.1 is
  scale-free, so the per-sample total does not affect it.
- **Support for D1.** `Archive/SFA3-1-1_Storyboard_VGC.docx` records the
  same rule used earlier in the project: zeros below detection limit
  replaced with half the minimum detected value for each variable.
- **O3 not resolved.** The methods state that per-C values are per-gram
  rates divided by mg C added per g dry soil, which matches the formula used
  here, but give no per-substrate amounts. Still waiting on Sheryl Bell.

## D8. Recommended gene sets (plan §4; proposed, awaiting approval)

VGC asked Claude to make its recommended choices for the Checkpoint 3
decisions. Rules, applied in order, to `KO_gene_sets_DRAFT.csv`:

1. Drop every low-confidence KO (15 KOs).
2. Leucine: drop the shared branched-chain amino acid transporters (LivKHMGF,
   BrnQ; 6 KOs). Set becomes leucine dehydrogenase, Ivd, MccA, MccB.
3. Succinic acid: "no specific gene set" (only shared transporters; first
   catabolic step is in the TCA cycle).
4. Near-zero rule: a set whose sum is 0 in 9 or more of 35 samples (more
   than a quarter) is marked "insufficient expression" and is not used in
   Q2 tests. This marks lignin (0 in 33 samples after rule 1) and chitosan
   (0 in 20).

Kept as they are: KOs shared between sets (cellulose/CB, Leu/Ala,
glycine/alanine); the eukaryotic hexosaminidase in the NAG enzyme set; the
"extracellular" label for all enzyme sets, as the plan states.

Output: `KO_gene_sets_RECOMMENDED.csv`, `..._expression.csv`,
`..._dropped.csv`.

- **Approved by:** proposed by Claude, 2026-10-26; awaiting VGC approval at
  Checkpoint 3.

## D9. Q2 implementation choices (plan §5; for VGC review at Checkpoint 4)

Choices taken by Claude where plan §5 does not specify:

- **Rows without a KO identifier** in `brite_D_expression.csv` are not KOs
  and are left out of the per-sample mean log expression (§5.1).
- **Potential on transcript dates (H2b)** uses every 100% plot with
  MicroResp on those dates (6 per date), not only the 5 plots with
  transcripts. Plots are matched only in the lag analysis, where the plan
  asks for it.
- **Bridge (§5.3.1)** subsets the Q1 date-mean scores (z-scores from all 77
  samples) to the 7 transcript dates; z-scores are not recomputed.
- **H2c lag and same-date correlations** use, for each date pair, the plots
  with transcripts at t; potential at t and at t+1 is averaged over those
  plots (4 plots for 2024-01-02 to 2024-02-12, SIR sets).
- **Bootstrap CI** for the lag advantage: percentile interval; resamples in
  which a correlation is undefined (zero variance) are dropped (count in
  `n_boot_valid`).
- **Plot-level lag models:** both models fitted by maximum likelihood on the
  same rows (plot-dates with potential at t and t+1), with response and
  predictor z-scored, so AIC and marginal R² are comparable. Marginal R² is
  var(fixed) / (var(fixed) + var(date) + var(residual)).
- **High-confidence robustness (§5.4b):** sets that become near-zero under
  D8 rule 4 are left out (glucosamine: 0 in 27 samples; malic acid: 0 in
  35).
- **Figure 4** shows the per-C robustness result as a second panel next to
  the plan's per-MBC result, because the per-MBC effects are nearly
  identical across compounds (see Checkpoint 4 report).

- **Approved by:** taken by Claude while implementing; for VGC review at
  Checkpoint 4.

## Open items (not yet decided)

- **O3. Per-C conversion check fails for three substrates (plan §2.2; STOP).**
  Recomputing `SIR_gdw × gdw_well / mgC_per_well` reproduces the file's
  per-C values for 28 of 31 substrates, with small constant offsets (up to
  4%) explained by `mgC_per_well` being stored to two decimals (for example
  oxalic acid: column 0.10, implied 0.104). Three substrates do not match.
  The mg C implied by the file's own per-C values is:
  - **Succinic acid** (Simple group): 2.25 on every date, against 0.16 in the
    column. 2.25 is the salicylic acid value, so the file's per-C succinic
    acid values are about 14 times lower than the column implies.
  - **Chitin** (Polymer group): 1.73 on 2023-06-07 to 2023-10-04, then
    sample-specific values from 0.43 to 2.20 on 2023-11-01 to 2024-06-05,
    against 1.3 in the column.
  - **Methyl oleate** (excluded from groups): 16.8 on 2023-07-05 to
    2023-10-04, 2.84 on all other dates, against 2.84 in the column.
  Per-gram dry soil values are not affected (implied `gdw_well` matches the
  column for every sample). Decision needed on which mg C values are correct
  before `induced_C` is final.

- **O1. Meaning of the 24 h reading (plan §3.2.8).** The README defines
  `Timepoint` only as the "end point of the measurement period in hours
  (5 or 24)". It does not say whether the 24 h rate covers 5–24 h or
  0–24 h. Needed before Checkpoint 2, test 8.
- **O2. Units of `Data/brite_D_expression.csv` (plan §5.1).** Not documented
  in the repository. Sample totals are 1.2 to 2.4 million and values are not
  integers, so the file is neither TPM nor raw counts. Two rows have no KO
  name. Needed before Checkpoint 4.
