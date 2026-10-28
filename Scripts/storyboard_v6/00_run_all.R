# =============================================================================
# Run every storyboard v6 script and assemble Results/storyboard_v6
# =============================================================================
#
# Run from the repository root. Each script runs in a fresh R session; the
# per-figure result parts are then combined into results_table.csv, and the
# substrate class file used is copied next to the figures.

rscript <- file.path(R.home("bin"), "Rscript")
scripts <- c("01_step0_data.R", "02_fig1.R", "03_fig2_microresp.R",
             "04_fig3_fig4_transcripts.R", "05_figS1_roots.R",
             "06_figS6_permanova.R")
for (s in scripts) {
  message("Running ", s)
  status <- system2(rscript, file.path("Scripts", "storyboard_v6", s))
  if (status != 0) stop(s, " failed")
}

suppressPackageStartupMessages(library(tidyverse))
dir_v6 <- file.path("Results", "storyboard_v6")
figure_order <- c("Fig1", "Fig2", "Fig3", "Fig4", paste0("FigS", 1:9))
parts <- list.files(file.path(dir_v6, "results_parts"), pattern = "\\.csv$",
                    full.names = TRUE)
results <- map_dfr(parts, \(f) read_csv(f, col_types = cols(.default = col_character()))) |>
  mutate(figure = factor(figure, figure_order)) |>
  arrange(figure)
stopifnot(!anyNA(results$figure))
write_csv(results, file.path(dir_v6, "results_table.csv"), na = "")
file.copy(file.path("Data", "substrate_classes.csv"),
          file.path(dir_v6, "substrate_classes.csv"), overwrite = TRUE)
message("results_table.csv: ", nrow(results), " rows")
