# =============================================================================
# Run the full manuscript v4 pipeline in order
# =============================================================================
#
# Run from the repository root:  Rscript Scripts/00_run_all.R
# Requires the untracked `DataHub March/` folder (see Data/README.md).
# Each step runs in a fresh R process so scripts cannot depend on objects
# left behind by an earlier step.

steps <- c(
  "01_Data_merging_and_cleaning.R",   # cleaned field + EEA table
  "02_sample_availability.R",         # sample audit, label check, lag pairs
  "03_compound_class_mapping.R",      # a priori classes (D1, D2)
  "10_figure1_context.R",             # Figure 1
  "11_figureS1_root_respiration.R",   # Figure S1
  "20_figure2_microresp.R",           # Figure 2, S2 (Q1)
  "30_figure3_transcripts.R",         # Figure 3, S3 (Q2)
  "40_figure4_correspondence.R"       # Figure 4, S4 (Q3)
)

rscript <- file.path(R.home("bin"), "Rscript")
log_dir <- file.path("Results", "manuscript_v4", "run_logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

for (step in steps) {
  message("Running ", step)
  log_file <- file.path(log_dir, sub("\\.R$", ".log", step))
  status <- system2(rscript, file.path("Scripts", step), stdout = log_file,
                    stderr = log_file)
  if (status != 0) stop(step, " failed; see ", log_file)
}

writeLines(capture.output(sessionInfo()),
           file.path("Results", "manuscript_v4", "sessionInfo.txt"))
message("All steps completed.")
