# Manuscript v4 scripts

Run everything from the repository root with `Rscript Scripts/00_run_all.R`
(R 4.6.1; about 5 minutes). The raw DataHub files must be present in
`DataHub March/` (see `Data/README.md` for origins and checksums).

| Order | Script | Purpose | Main outputs |
|---|---|---|---|
| 0 | `00_run_all.R` | Runs the steps below in fresh R sessions | `Results/manuscript_v4/run_logs/`, `sessionInfo.txt` |
| 1 | `01_Data_merging_and_cleaning.R` | Merge biogeochemistry and EEA; documented cleaning | `Data/Data_merged_cleaned.csv`, `Results/Data_merging_and_cleaning/` |
| 2 | `02_sample_availability.R` | Sample audit, MicroResp label check, one-month lag pairs | `Results/manuscript_v4/sample_availability/` |
| 3 | `03_compound_class_mapping.R` | A priori substrate classes and uniquely mapped BRITE C categories | `Data/compound_class_mapping.csv`, `Data/compound_class_unique_BRITE_C.csv` |
| 4 | `10_figure1_context.R` | Figure 1: design and seasonal field context | `Figures|Results/manuscript_v4/Figure1/` |
| 5 | `11_figureS1_root_respiration.R` | Figure S1: root biomass vs field respiration | `.../FigureS1/` |
| 6 | `20_figure2_microresp.R` | Figure 2 and S2: seasonal SIR (Q1, H1a, H1b) | `.../Figure2/` |
| 7 | `30_figure3_transcripts.R` | Figure 3 and S3: seasonal BRITE C expression (Q2, H2) | `.../Figure3/` |
| 8 | `40_figure4_correspondence.R` | Figure 4 and S4: SIR-transcript correspondence (Q3) | `.../Figure4/` |

`helpers.R` holds what every script shares: input paths, the season
definition, data loaders, the exact date-level season permutation test, and
the plotting theme.

`docs/` writes the storyboard and analysis notes from the results:
`Rscript Scripts/docs/write_storyboard_html.R` and then
`powershell -ExecutionPolicy Bypass -File Scripts/docs/html_to_docx.ps1`
(needs Microsoft Word). Output goes to `Docs/`.

Superseded scripts (the v3 pipeline and the October 2026 outlier exploration)
are in `Archive/2026-10-23_pre_v4/`, and in git history up to commit `80e0df9`.
