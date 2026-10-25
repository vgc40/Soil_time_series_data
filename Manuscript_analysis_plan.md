# Manuscript analysis plan (draft v1, 2026-10-23)

> **Status (2026-10-23):** implemented on branch `manuscript-v4-analysis`. Scripts are listed in `Scripts/README.md`, results are summarised in `Docs/Storyboard_v4.docx`, and open items are in `Docs/Analysis_notes_v4.docx`. Where they differ, the decision log in section 6 overrides earlier sections.

Goal: bring the analysis pipeline in line with the current storyboard (Q1–Q3), fix known data issues, and produce one main figure per question plus a context figure. Draft for review; nothing below has been implemented yet.

---

## 0. Where the current pipeline stands vs. the storyboard

| Storyboard element | Current state | Gap |
|---|---|---|
| Data input | Scripts 00, 02, 03 read the old `..._MicroResp.csv` | Old file has the January 2024 ID offset + duplicate (SM1503 copy of SM1502). `_v3` fixes it. |
| Q1 / H1a (SIR profile differs by season) | Fig 2A NMDS + PERMANOVA (R² = 0.10, p = 0.002) | Permutes samples within plot, so the 6 plots per date count as independent replicates of season. Dispersion also differs (betadisper p = 0.007), and that isn't reported alongside the result. |
| Q1 / H1b (readily available vs complex substrates) | Fig 2B–D show "top 3 substrates by GAM deviance" (methyl oleate, chitosan, L-glutamine) | Not tied to H1b. No a priori readily-available/complex grouping and no test of a class × season difference. |
| Q2 / H2 (transcript composition differs by season) | Fig 3A NMDS + PERMANOVA (R² = 0.23, p = 0.001) | Same pseudo-replication issue. With 5 active + 2 dormant dates, a date-level permutation has only 21 distinct arrangements (smallest possible p ≈ 0.048). |
| Q2 descriptive (uniquely mapped categories per class) | Fig 3B heatmap of 31 "linked" pathways | Includes non-specific pathways mapped to many substrates (TCA cycle → 10 groups, glycolysis → 4). The storyboard asks for categories mapped *uniquely* to 4 compound classes. |
| Q3 (class SIR vs class transcripts, plus one-month lag) | Fig 3C top-3 z-score overlay; Fig S3 global Procrustes | No class-level comparison and no lag analysis. |
| Supporting data (EEA, MBC, respiration, roots, environment) | Fig 1 = root biomass × respiration mixed model + outlier review (9 candidates in the committed script, 3 in the untracked Oct-16 scripts) | No storyboard question uses root–respiration. Fig 1 should become seasonal context. |

---

## 1. Data fixes (do first; everything downstream depends on them)

1. **Switch to `SoilMicrobiomeSFA_prosser-time-series_MicroResp_v3.csv`** in `00_sample_availability.R`, `02_figure2_microresp.R` and `03_figure3_brite_expression.R`. Read `SIR` with `na = c("na", "NA")` (310 Water rows are `na`). Expected 100% MicroResp n changes from 78 to 77: Jan 2024 has 6 (the duplicate is gone) and Feb 2024 has 5 (plot 27 / SM1536 was never measured).
2. **Re-run the sample-availability audit** and confirm: 77 MicroResp samples, 35 transcript samples (plots 1, 25, 27, 40, 46 on all 7 dates; plot 57 never sequenced), and lag pairs n = 5, 5, 5, 5, 4, 5, 5.
3. **July 2023 sample removal** in `01_Data_merging_and_cleaning.R` (SM1082, SM1085, SM1087, SM1088). SM1082 is the 100% plot 57 sample, which is why field n = 77. **Decided (D6):** keep SM1082 and set only its SEOC to NA, so field n becomes 78.
4. **EEA metadata:** SM1502 has `Collection_date = #N/A` in the EEA file. The merge is by `Sample_cat`, so it is harmless, but report it to the data owner along with the README errors (collection-date line; "six MUB / four MUC" vs the 5 + 2 listed).
5. **Transcript data provenance:** no script in the repo creates `brite_C_expression.csv` / `brite_D_expression.csv`. Units and normalization (counts? TPM? and which human/host pathways were removed via `Data/pathwaystoremove/`) need to be documented for Methods, and the decision affects Q2/Q3 (D5).
6. **Housekeeping:** commit or archive the untracked Oct-16 Figure 1 outlier scripts; add a `00_run_all.R` with the run order; write `sessionInfo()` to Results.

---

## 2. A priori classification tables (needed before any figure)

### 2a. Compound classes for Q2/Q3 (4 classes)
Derived from `Data/Microresp_pathways.csv`. A BRITE C category counts as "unique" when it maps to substrates in only one class. **Decided (D2):** BSA and GABA belong to the amino-acid class, and the rule is applied with no exceptions, giving **20 categories**.

| Class | MicroResp substrates | Uniquely mapped BRITE C |
|---|---|---|
| Lipid-like (2 substrates) | Methyl oleate, Tween40 | Fatty acid degradation; Glycerolipid; Glycerophospholipid; Cutin, suberine & wax (low: nonzero in 30/35, median 3.7). **4** |
| Aromatic (4) | Salicylic acid, Vanillic acid, Lignin, Phenylethylamine | Benzoate deg.; Aminobenzoate deg.; Naphthalene deg.; Phenylalanine; Tyrosine; Phenylpropanoid biosynthesis. **6** |
| Amino sugar (4) | Chitin, Chitosan, Glucosamine-HCl, N-acetyl-D-glucosamine | Amino sugar & nucleotide sugar metabolism. **1** |
| Amino acid (8) | Glutamic acid, L-Arg, L-Gln, Gly, L-Ala, L-Leu, BSA, GABA | Ala/Asp/Glu; Gly/Ser/Thr; Arg & Pro; Arg biosynthesis; Val/Leu/Ile deg.; Cys & Met; One carbon pool by folate; β-Alanine; Butanoate. **9** |

Substrates outside the 4 classes (sugars, organic acids, cellulosic, diethylamine) are not used in Q2/Q3. The script will recompute the unique set from the mapping file and stop if it is not these 20.

### 2b. Substrate availability classes for H1b
**Decided (D1), strict sets:**
- **Readily available (10):** D-glucose, D-fructose, D-trehalose, D-xylose, D-cellobiose, citric, lactic, malic, oxalic and succinic acids.
- **Complex, plant/microbial residue-associated (7):** cellulose, hemicellulose, lignin, chitin, chitosan, BSA, methyl oleate.
- **Not in the H1b contrast (14):** the 6 amino acids, diethylamine, phenylethylamine, GABA, glucosamine, NAG, salicylic acid, vanillic acid and Tween40. These are shown in the figures labelled "other".
- **Sensitivity:** amino acids moved into "readily available".

Both tables get saved as `Data/compound_class_mapping.csv`. Every script reads from that file, so the mapping is stated once and is easy to cite.

---

## 3. Statistical design principles (apply to all figures)

- **Season is a date-level treatment.** The replication unit for season is sampling date (13 for SIR; 7 for transcripts: 5 active, 2 dormant). Plots are repeated measures.
- **Primary season tests:**
  - Multivariate: PERMANOVA on date-mean profiles (SIR: 8 active vs 5 dormant dates, 1,287 possible permutations; transcripts: 21 possible, reported as exact).
  - Sensitivity: sample-level PERMANOVA with season labels permuted across whole dates and plot as strata.
  - Univariate: mixed models `y ~ Season + (1 | Date) + (1 | Plot)`.
- **Always report dispersion (betadisper)** next to each PERMANOVA. Where dispersion differs, describe the result as "location and/or dispersion".
- **H1b is a relative test.** Responses are z-scored (or log-ratioed) within substrate relative to that substrate's mean across dates, as the storyboard specifies. A shared seasonal shift in basal respiration then cancels out, and the class × season interaction is the actual test.
- **Q2 categories and all of Q3 are descriptive:** effect sizes, trajectories and rank correlations, with no per-category p-values. The storyboard already says this.
- **SIR settings:** 5 h, carbon-normalized (per mg C added) for the main figures. Sensitivity in SI: per g dry soil, and 24 h.

---

## 4. Figures

### Figure 1: Study design and seasonal field context (supporting measurements)
- **A.** Sampling timeline: 13 dates with plant-active/dormant shading and markers for the 7 transcript dates. Shows what was measured and when.
- **B.** Soil temperature and gravimetric water content.
- **C.** Field soil respiration and root biomass.
- **D.** Microbial biomass C/N and EEA summed C-, N- and P-acquiring activity (BG, BX, CB / NAG, Leu, Ala / AP).
- **Stats:** a season mixed model per variable, given in an SI table.
- **Root–respiration model and outlier sensitivity:** move to SI or drop (D3).

### Figure 2: Q1, seasonal substrate-induced respiration (H1a, H1b)
- **A. H1a:** NMDS/PCoA of 31-substrate SIR profiles, coloured by season, with date centroids joined in time order. Date-level PERMANOVA + betadisper.
- **B. H1b, per substrate:** dormant − active season effect on within-substrate z-scores, with 95% CI, from mixed models. Ordered and coloured by availability class; FDR across 31 substrates.
- **C. H1b, group test:** one mixed model `z ~ Season * AvailabilityClass + (1 | Substrate) + (1 | Date) + (1 | Plot)`. The interaction estimate is the H1b result. Show the class means by season.
- **D.** Monthly trajectories (date mean ± CI of z) per availability class.
- **SI:**
  - All 31 substrate trajectories (reuse the current GAM panels) and the heatmap/clusters.
  - 24 h and per-g-soil sensitivity.
  - Envfit with environment, as exploratory.
- **Current top-3 panels:** drop from the main figure (D4).

### Figure 3: Q2, seasonal transcript expression (H2 + descriptive class patterns)
- **A. H2:** NMDS of BRITE C composition (35 samples, 7 dates), coloured by season. Date-level PERMANOVA (exact, 21 permutations) and a sample-level sensitivity; betadisper.
- **B.** Heatmap of the 17 (or 18) uniquely mapped categories: date-mean z-scores within category, rows grouped by compound class, columns ordered by date with season bar.
- **C.** Class-level trajectories: mean z of each class's unique categories over the 7 dates. Descriptive.
- **SI:** full BRITE C heatmap; low-abundance caveat for Cutin/suberine/wax.

### Figure 4: Q3, correspondence between SIR and transcripts by class
For each of the 4 classes:
- **A. Concurrent:** the 7 transcript dates, with paired trajectories of class SIR (mean z of the class's substrates) and class transcripts (mean z of its unique categories). Plus a scatter with Spearman ρ (n = 7 dates; descriptive).
- **B. One-month lag:** transcripts at date t vs SIR at the next sampling date. Pairs are Jun→Jul, Jul→Aug, Sep→Oct, Oct→Nov, Jan→Feb, Mar→Apr, May→Jun, with n plots = 5, 5, 5, 5, 4, 5, 5. Compare lagged with concurrent ρ per class. Note that the gaps are not equal (Jan 2 → Feb 12 is about 6 weeks).
- **SI:** the current global Procrustes/PROTEST (Fig S3), kept as a whole-profile check, plus plot-level versions of A and B.

---

## 5. Order of work
1. Data fixes (section 1) → re-run availability audit.
2. Write `Data/compound_class_mapping.csv` after the D1/D2 decisions.
3. Figure 2 (Q1) → Figure 3 (Q2) → Figure 4 (Q3) → Figure 1 (context).
4. Captions, SI tables, results summary per figure; check every reported n against the audit.
5. Commit in logical steps (data fixes / mapping / each figure).

---

## 6. Decisions (resolved 2026-10-23; this log overrides earlier sections where they differ)

| # | Decision |
|---|---|
| D1 | *Revised 2026-10-25:* one rule. Complex = polymers and esters that must be broken apart outside the cell (cellulose, hemicellulose, lignin, chitin, chitosan, BSA, methyl oleate, Tween 40); readily available = the other 23 small molecules; no "other". Sensitivity: only root-exudate types (sugars, organic acids, amino acids) count as readily available. Polymer vs monomer comparisons added within chemical classes. (Section 2b shows the original strict sets.) |
| D2 | Amino-acid class includes BSA and GABA. Strict uniqueness rule, giving 20 categories (section 2a). |
| D3 | Root × respiration becomes a short SI analysis. Main Fig 1 shows root biomass only as seasonal context. |
| D4 | Top-3 panels dropped from the main figures. All 31 GAMs go in Fig S2, and the top 3 are named in the text. |
| D5 | Transcripts converted to within-sample proportions of total expression over the analysed BRITE C set. Values as provided kept as an SI sensitivity. **Still needed:** who produced `brite_C_expression.csv`, and its upstream unit (values are non-integer; per-sample totals range 0.88–4.08 M). |
| D6 | SM1082 (Jul 2023, plot 57): set only its SEOC (2,410 µg/g, ~60× July median) to NA and keep the sample. The 100% field data become 78/78 plot-dates. Removal of the three 25% samples is unchanged (not used here). |
| D7 | Apr–Oct / Nov–Mar stays primary. SI sensitivity drops the transition dates (Oct 2023, Apr 2024); for transcripts that leaves 4 active + 2 dormant dates. Citation to come from Sheryl. |
| D8 | Class aggregation. Transcripts: sum of the class's unique categories as a share of expression, then z-scored across dates. SIR: mean of within-substrate z. Sensitivity: transcripts as mean within-category z. |
| D9 | Fig 3A uses all 325 categories. SI sensitivity uses metabolism only (KEGG map numbers 00xxx/01xxx). |
| D10 | Fig 1 timeline shows season shading plus the Apr–Oct irrigation window. Harvest/saturation markers and the rhizotron data are on hold until Sheryl sends dates/files. |
| D11 | All root observations kept in the primary analysis. Sensitivity excludes only the 3 points flagged by both screens (SM1621, SM1078, SM1080). The 9-point version is dropped. |
| D12 | *2026-10-25:* H1a reported as two results, level (mean over substrates) and mix (profiles scaled to each sample's total). NMDS drawn without connecting lines. Q3 adds a leave-one-date-out check. Storyboard rewritten for non-specialist readers. |

**To ask Sheryl:** season-definition reference (D7); harvest, fall-saturation and irrigation start/stop dates, plus rhizotron data (D10); whether field substrate concentrations can be compared with the assay dose (KH 96).

## 6a. Original decision list (for reference)
- **D1:** Readily available vs complex substrate lists for H1b (section 2b), including BSA and Tween40.
- **D2:** Whether BSA and GABA belong in the amino-acid class, and whether One carbon pool by folate is included (section 2a).
- **D3:** Root–respiration analysis: SI, or drop?
- **D4:** Drop the "top 3 by GAM deviance" panels from the main figures?
- **D5:** Transcript normalization for Q2/Q3: as provided, relative abundance (proportion of total BRITE C), or another? And can you confirm how `brite_C_expression.csv` was produced?
- **D6:** Keep the July 2023 removal of SM1082/1085/1087/1088? On what basis?
- **D7:** Season boundaries: storyboard uses Apr (emergence) to Oct (harvest), with the citation still marked "[basis for season definitions]". OK to keep the current month-based split? (Still waiting on Sheryl, per v3 comments 31/33.)
- **D8:** How to aggregate a class (KH, comment 90). SIR is a bulk sum of activity, while metaT is per function. Options: sum of the class's categories as relative abundance, or the mean of within-category z-scores. Plan default: use the sum as primary and the mean z as a sensitivity.
- **D9:** Fig 3A input (KH, comments 74/77). *Mostly settled by your v3 replies (VG36/VG42):* keep all 325 BRITE C categories in Fig 3A. Open question: should a metabolism-only version go in the SI as a sensitivity?
- **D10:** Do event dates exist (irrigation on/off, harvest, fall saturation) for the Fig 1 timeline? Sheryl's comments 49/60 suggest she has them. Also, should the rhizotron data (comment 54) go into Fig 1?
- **D11:** Root outliers. VGC (comment 40) proposed removing any point flagged by both the global and within-plot screens (3 points). KH (comment 39) asked whether they are hot spots/moments rather than outliers. Plan default: keep all points and show the both-flagged exclusion only as a sensitivity.

---

## 7. Changes prompted by the v3 storyboard comments (`VGC_Storyboard_v3_comments.json`)

| Comment | Plan change |
|---|---|
| KH 64: "why are we highlighting plots… date and substrate are the story" | Fig 2A/3A colour by season, show date centroids joined in time order, and drop plot highlighting. Fits the date-level design. |
| AZ 65/94: potential activity is higher in dormancy for all substrates; this is the headline | Add a panel to Fig 2 with the absolute season effect (dormant/active ratio of raw SIR per substrate), alongside the relative H1b test. |
| KH 47/71, VGC 70: contrast seasonal vs year-round-steady substrates | Fig S2 stays ordered by deviance explained. Fig 2B also marks substrates with no detectable seasonality (e.g. hemicellulose, malic acid). |
| AZ 89: are some months more aligned than others? | Fig 4 / SI: per-date alignment (Procrustes residual per date, or per-date agreement of class z-scores). |
| KH 90/91: sum vs average for a comparable measure | Decision D8. |
| KH 74/77/78: what goes into the BRITE C NMDS; remove basic cell functions? | Decision D9, and Methods text on how metaT was rolled up to BRITE C (ties to D5 provenance). |
| KH 46: is the mapping based on annotations? | Add `rationale` and `source` columns to `compound_class_mapping.csv` (KEGG pathway membership of substrate catabolism). |
| KH 79: "just 3 – include the threshold" | Moot if top-3 is dropped (D4). Otherwise state the selection rule. |
| KH 59: wet–dry / Birch effect | SIR sensitivity: add field GWC at sampling as a covariate (rewetting to 60% WHC may cause a pulse in dry months). Discussion point. |
| KH 96: benchmark field substrate concentrations vs assay dose | Not possible from repo data. Ask Sheryl; flag for Discussion. |
| SB 51/53: roots are root-picking from pooled 0–15 cm cores, not rhizosphere | Fig 1 labels and Methods say "root biomass (0–15 cm, hand-picked)". |
| Title/framing/intro comments (KH 0, 2, 5; AZ 1, 7, 20, 26) | Writing, not analysis. Not in this plan. |

---

## 8. What to carry over from the full v3 storyboard text

**Keep / reuse:**
- **Fig S2, all 31 GAM trajectories.** The methods carry over unchanged (non-cyclic smooth, plot random effect, FDR). Layout: VG32 asks for grouping by compound class, KH33 liked ordering by % explained. Do both: group panels by class and order by deviance explained within each class.
- **Envfit result:** GWC and soil temperature were associated with the SIR ordination. Keep it as exploratory SI for Fig 2, and re-run it with the date-level permutation. It also gives an angle on the wet–dry/Birch question (KH 59).
- **The main v3 observation for Q3:** methyl oleate, chitosan and L-glutamine SIR rose in winter/early spring while the mapped transcripts hit their low around January. Fig 4 must therefore show the *sign* of correspondence (opposite trajectories), not only its strength. The concurrent vs one-month-lag comparison tests whether the mismatch is a timing offset.
- **Fig 1C, field respiration:** keep the v3 presentation (plot points, monthly mean ± 95% CI, October 2023 peak, sharp drop into November). Per SB15, annotate October as pre-saturation/harvest and November as post-saturation.
- **Procrustes/PROTEST with cyclic-shift permutations:** keep as SI (Fig S3), plus per-date residuals for AZ 48.
- **Methods text already written:**
  - why expression is rolled up to taxon-independent BRITE C (AZ 8)
  - the description of MicroResp as potential capacity after a resource spike
  - the z-score caveat that timing and direction are comparable but magnitudes are not

**Drop or demote (no longer in the current storyboard):**
- **Old H1a (root biomass → field respiration, seasonal slopes):** demote to SI context (D3). The sensitivity numbers in v3 (9 removed: active slope 0.49 → 2.55, p = 0.0035) predate the 3-point both-screens rule (2.20, p = 0.005). If this stays in SI, report only one rule (D11).
- **"Top 3 substrates by deviance explained" as the bridge between Fig 2 and Fig 3:** replaced by the a priori compound classes. The top 3 can still be named in the text, since methyl oleate, chitosan and L-glutamine sit in the lipid-like, amino-sugar and amino-acid classes.
- **Old H2 (positive correspondence expected):** the current storyboard treats Q3 as descriptive, so no directional prediction.

**Housekeeping from v3:**
- The ChatGPT-use note needs a matching AI-use statement in the manuscript (check SBB policy).
- The code-availability statement should point to the GitHub repo (`vgc40/Soil_time_series_data`).
