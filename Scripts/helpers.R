# =============================================================================
# Shared definitions for the manuscript v4 pipeline
# =============================================================================
#
# Sourced by every analysis script so that the season definition, input file
# paths, data loaders, permutation tests, and plotting theme are stated once.
# Run scripts from the repository root (see Scripts/00_run_all.R).

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(vegan)
})

# -----------------------------------------------------------------------------
# Input files (see Data/README.md for provenance and checksums)
# -----------------------------------------------------------------------------

input_files <- list(
  microresp = file.path(
    "DataHub March", "Soil_Phenotypes",
    "SoilMicrobiomeSFA_prosser-time-series_MicroResp_v3.csv"
  ),
  field = file.path("Data", "Data_merged_cleaned.csv"),
  brite_c = file.path("Data", "brite_C_expression.csv"),
  human_exclusion = file.path(
    "Data", "pathwaystoremove", "kegg_human_exclude_2026-07-13.txt"
  ),
  microresp_pathways = file.path("Data", "Microresp_pathways.csv"),
  substrate_classes = file.path("Data", "compound_class_mapping.csv"),
  unique_brite_c = file.path("Data", "compound_class_unique_BRITE_C.csv")
)

# -----------------------------------------------------------------------------
# Season definition (decision D7)
# -----------------------------------------------------------------------------
# Plant-active: April (perennial emergence) through October (harvest).
# Plant-dormant: November through March.
# Transition dates excluded in the boundary sensitivity analysis.

season_levels <- c("Plant-active", "Plant-dormant")
season_colors <- c("Plant-active" = "#2E8B57", "Plant-dormant" = "#5B6F8C")
season_fill <- c("Plant-active" = "#E3F1E8", "Plant-dormant" = "#E4E8EF")
transition_dates <- as.Date(c("2023-10-04", "2024-04-03"))

assign_season <- function(date) {
  factor(
    if_else(month(date) %in% 4:10, "Plant-active", "Plant-dormant"),
    levels = season_levels
  )
}

# Shaded season bands spanning the sampling period, for time-series panels.
season_bands <- function(dates, pad_days = 15) {
  all_days <- seq(min(dates) - pad_days, max(dates) + pad_days, by = "day")
  tibble(day = all_days, Season = assign_season(all_days)) %>%
    mutate(run = cumsum(Season != lag(Season, default = first(Season)))) %>%
    group_by(run, Season) %>%
    summarise(start = min(day), end = max(day) + 1, .groups = "drop")
}

# Season shading as fixed-colour layers, leaving the fill scale free for data.
season_band_layers <- function(dates, alpha = 0.8) {
  bands <- season_bands(dates)
  map(season_levels, function(s) {
    geom_rect(
      data = filter(bands, Season == s),
      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE, fill = season_fill[[s]], alpha = alpha
    )
  })
}

# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

class_colors <- c(
  "Readily available" = "#D08C2E",
  "Complex" = "#6B4C9A",
  "Other" = "#9A9A9A",
  "Lipid-like" = "#C0504D",
  "Aromatic" = "#4F81BD",
  "Amino sugar" = "#9BBB59",
  "Amino acid" = "#8064A2"
)

measure_colors <- c("SIR" = "#B8562B", "Transcripts" = "#2B6CB8")

theme_manuscript <- theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    plot.title = element_text(face = "bold", size = 11),
    plot.tag = element_text(face = "bold", size = 12),
    legend.key.size = unit(0.4, "cm")
  )

date_axis <- scale_x_date(date_breaks = "1 month", date_labels = "%b\n%y")

save_figure <- function(plot, dir, name, width, height) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(dir, paste0(name, ".png")), plot, width = width,
         height = height, dpi = 300, bg = "white")
  ggsave(file.path(dir, paste0(name, ".pdf")), plot, width = width,
         height = height, bg = "white")
}

output_dirs <- function(figure) {
  dirs <- list(
    figures = file.path("Figures", "manuscript_v4", figure),
    results = file.path("Results", "manuscript_v4", figure)
  )
  walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  dirs
}

# -----------------------------------------------------------------------------
# Data loaders
# -----------------------------------------------------------------------------

# 100% irrigation MicroResp SIR, one row per sample x substrate.
load_microresp_100 <- function(timepoint = 5, normalization = "ngCO2-C_mgC_h") {
  read_csv(
    input_files$microresp,
    col_types = cols(.default = col_character()),
    na = c("", "NA", "na")
  ) %>%
    filter(
      Irrigation_trmt_pct == "100",
      as.numeric(Timepoint) == timepoint,
      Normalization == normalization,
      Substrate != "Water"
    ) %>%
    transmute(
      Sample_cat,
      Sampling_Date = mdy(Collection_date),
      Plot_ID = factor(Plot),
      Block = factor(Block),
      Substrate,
      SIR = as.numeric(SIR),
      Season = assign_season(Sampling_Date)
    )
}

# 100% irrigation field data (biogeochemistry + EEA), one row per sample.
load_field_100 <- function() {
  read_csv(input_files$field, show_col_types = FALSE) %>%
    filter(as.character(Irrigation_trmt_pct) == "100") %>%
    mutate(
      Sampling_Date = as.Date(Collection_date),
      Plot_ID = factor(Plot),
      Block = factor(Block),
      Season = assign_season(Sampling_Date)
    )
}

# BRITE C transcript expression for the 100% treatment, long format.
# Rows without a KEGG pathway code, curated human/non-soil pathways, and
# categories with zero expression in every sample are removed (325 remain). `proportion` is each category's share of a
# sample's total expression over the retained categories (decision D5).
load_brite_c <- function() {
  wide <- read_csv(input_files$brite_c, show_col_types = FALSE) %>%
    rename(brite_C = 1) %>%
    filter(!is.na(brite_C), brite_C != "")

  human_codes <- read.delim(input_files$human_exclusion, check.names = FALSE) %>%
    transmute(pathway_code = str_extract(as.character(path_ko), "ko\\d{5}")) %>%
    filter(!is.na(pathway_code)) %>%
    distinct() %>%
    pull(pathway_code)

  metadata <- load_field_100() %>%
    select(Sample_cat, Sampling_Date, Plot_ID, Block, Season)

  wide %>%
    mutate(pathway_code = str_extract(brite_C, "ko\\d{5}")) %>%
    filter(!is.na(pathway_code), !pathway_code %in% human_codes) %>%
    pivot_longer(-c(brite_C, pathway_code), names_to = "Sample_mt",
                 values_to = "expression") %>%
    mutate(Sample_cat = str_remove(Sample_mt, "_MT$")) %>%
    left_join(metadata, by = "Sample_cat") %>%
    group_by(brite_C) %>%
    filter(sum(expression) > 0) %>%
    group_by(Sample_cat) %>%
    mutate(proportion = expression / sum(expression)) %>%
    ungroup()
}

# -----------------------------------------------------------------------------
# Season tests with sampling date as the replication unit
# -----------------------------------------------------------------------------
# Season is assigned to whole sampling dates, so permutations reassign season
# labels to dates, never to individual plots. With few dates the full set of
# label assignments is enumerated and the p-value is exact.

season_label_sets <- function(date_season) {
  n_dormant <- sum(date_season == "Plant-dormant")
  combn(length(date_season), n_dormant, simplify = FALSE) %>%
    map(function(idx) {
      labels <- rep("Plant-active", length(date_season))
      labels[idx] <- "Plant-dormant"
      factor(labels, levels = season_levels)
    })
}

pseudo_f <- function(dist_matrix, groups) {
  fit <- suppressMessages(adonis2(dist_matrix ~ groups, permutations = 0))
  c(F = fit$F[1], R2 = fit$R2[1])
}

dispersion_f <- function(dist_matrix, groups) {
  bd <- betadisper(dist_matrix, groups)
  anova(bd)$`F value`[1]
}

# Exact date-level PERMANOVA and dispersion test.
# `profile` is a samples x variables matrix; `sample_dates` gives each row's
# sampling date. Returns tests on date-mean profiles (primary) and on samples
# with labels permuted by whole date (sensitivity).
date_level_season_test <- function(profile, sample_dates, method = "bray") {
  dates <- sort(unique(sample_dates))
  date_season <- assign_season(dates)
  label_sets <- season_label_sets(date_season)

  date_means <- rowsum(profile, as.character(sample_dates)) /
    as.vector(table(as.character(sample_dates)))
  date_means <- date_means[as.character(dates), , drop = FALSE]
  d_date <- vegdist(date_means, method = method)
  d_sample <- vegdist(profile, method = method)
  sample_index <- match(sample_dates, dates)

  run <- function(d, expand) {
    obs_groups <- if (expand) date_season[sample_index] else date_season
    obs <- pseudo_f(d, obs_groups)
    obs_disp <- dispersion_f(d, obs_groups)
    perm <- map(label_sets, function(lbl) {
      g <- if (expand) lbl[sample_index] else lbl
      c(pseudo_f(d, g)["F"], disp = dispersion_f(d, g))
    }) %>% do.call(rbind, .)
    tibble(
      F = obs[["F"]], R2 = obs[["R2"]],
      p_permanova = mean(perm[, "F"] >= obs[["F"]] - 1e-12),
      F_dispersion = obs_disp,
      p_dispersion = mean(perm[, "disp"] >= obs_disp - 1e-12),
      n_label_sets = length(label_sets),
      min_possible_p = 1 / length(label_sets)
    )
  }

  bind_rows(
    run(d_date, expand = FALSE) %>%
      mutate(level = "Date means (primary)", n_units = length(dates)),
    run(d_sample, expand = TRUE) %>%
      mutate(level = "Samples, labels permuted by date", n_units = nrow(profile))
  ) %>%
    mutate(
      n_active_dates = sum(date_season == "Plant-active"),
      n_dormant_dates = sum(date_season == "Plant-dormant")
    ) %>%
    relocate(level, n_units)
}

format_p <- function(p) {
  if_else(p < 0.001, "< 0.001", formatC(p, format = "f", digits = 3))
}
