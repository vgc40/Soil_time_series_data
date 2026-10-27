# =============================================================================
# Shared definitions for analysis plan v1 (analysis_plan_v1.md)
# =============================================================================
#
# Sourced by every plan_v1 script so that input paths, season dates, compound
# groups, and score construction are stated once. Run scripts from the
# repository root.

suppressPackageStartupMessages({
  library(tidyverse)
})

# -----------------------------------------------------------------------------
# Input files (plan §1; source-of-truth rule: DataHub March first, see D2)
# -----------------------------------------------------------------------------

input_files <- list(
  microresp = file.path(
    "DataHub March", "Soil_Phenotypes",
    "SoilMicrobiomeSFA_prosser-time-series_MicroResp_v3.csv"
  ),
  eea = file.path(
    "DataHub March", "Soil_Phenotypes",
    "SoilMicrobiomeSFA_prosser-time-series_EEA.csv"
  ),
  bgc = file.path(
    "DataHub March", "Soil_Biogeochemistry",
    "SoilMicrobiomeSFA_prosser-time-series_Soil_Biogeochemistry.csv"
  ),
  brite_d = file.path("Data", "brite_D_expression.csv"),
  chitin_kos = file.path(
    "Data", "merged_chitin_chitosan_methyl_oleate_gene_categories.csv"
  )
)

dir_results <- file.path("Results", "plan_v1")
dir_figures <- file.path("Figures", "plan_v1")
dir.create(dir_results, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_figures, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Season (plan §2.3)
# -----------------------------------------------------------------------------

dormant_dates <- as.Date(c(
  "2023-11-01", "2023-12-04", "2024-01-02", "2024-02-12", "2024-03-06"
))

study_dates <- as.Date(c(
  "2023-06-07", "2023-07-05", "2023-08-02", "2023-09-06", "2023-10-04",
  "2023-11-01", "2023-12-04", "2024-01-02", "2024-02-12", "2024-03-06",
  "2024-04-03", "2024-05-08", "2024-06-05"
))

assign_season <- function(date) {
  factor(if_else(date %in% dormant_dates, "Dormant", "Active"),
         levels = c("Active", "Dormant"))
}

# -----------------------------------------------------------------------------
# Compound groups (plan §2.4; fixed)
# -----------------------------------------------------------------------------

compound_groups <- tribble(
  ~Substrate, ~group,
  "D.Fructose", "Simple", "D.Xylose", "Simple", "D.Trehalose", "Simple",
  "D.Cellobiose", "Simple", "Citric_acid", "Simple", "Malic_acid", "Simple",
  "Oxalic_acid", "Simple", "Succinic_acid", "Simple", "Lactic_acid", "Simple",
  "Glutamic_acid", "Simple", "L.Arginine", "Simple", "L.Glutamine", "Simple",
  "Glycine", "Simple", "L.Alanine", "Simple", "L.Leucine", "Simple",
  "Cellulose", "Polymer", "Hemicellulose", "Polymer", "Lignin", "Polymer",
  "Chitin", "Polymer", "Chitosan", "Polymer",
  "D.Glucose", "Reference",
  "N.Acetyl.D.glucosamine", "Q2", "Glucosamine_hydrochloride", "Q2"
)

simple_compounds <- compound_groups$Substrate[compound_groups$group == "Simple"]
polymer_compounds <- compound_groups$Substrate[compound_groups$group == "Polymer"]

# -----------------------------------------------------------------------------
# Standardized scores (plan §2.6)
# -----------------------------------------------------------------------------
# z-score each compound across the samples supplied, then average z within
# the simple and polymer groups for each sample. `value` names the column to
# standardize, so the same function serves per-C, per-MBC, and the
# robustness variants.

make_scores <- function(df, value) {
  z <- df |>
    filter(Substrate %in% c(simple_compounds, polymer_compounds)) |>
    group_by(Substrate) |>
    mutate(z = (.data[[value]] - mean(.data[[value]], na.rm = TRUE)) /
             sd(.data[[value]], na.rm = TRUE)) |>
    ungroup() |>
    select(-any_of("group")) |>
    left_join(compound_groups, by = "Substrate")
  z |>
    group_by(Sample_cat, Collection_date, Plot, season, group) |>
    summarise(score = mean(z, na.rm = TRUE), n_compounds = sum(!is.na(z)),
              .groups = "drop") |>
    pivot_wider(names_from = group, values_from = c(score, n_compounds)) |>
    rename(simple_score = score_Simple, polymer_score = score_Polymer)
}

# Long table of per-compound z-scores, for compound-level effects and the
# mixed models
make_z_long <- function(df, value, compounds) {
  df |>
    filter(Substrate %in% compounds) |>
    group_by(Substrate) |>
    mutate(z = (.data[[value]] - mean(.data[[value]], na.rm = TRUE)) /
             sd(.data[[value]], na.rm = TRUE)) |>
    ungroup()
}

date_means <- function(scores) {
  scores |>
    group_by(Collection_date, season) |>
    summarise(n_plots = n(),
              across(c(simple_score, polymer_score), \(x) mean(x, na.rm = TRUE)),
              .groups = "drop") |>
    mutate(diff_score = polymer_score - simple_score)
}

# -----------------------------------------------------------------------------
# Season comparison on date means (plan §3.1.1, §6)
# -----------------------------------------------------------------------------
# Dormant minus active, Welch 95% CI, Welch and Wilcoxon rank-sum p-values.
# The date is the unit of replication.

season_test <- function(value, season) {
  d <- value[season == "Dormant"]
  a <- value[season == "Active"]
  welch <- t.test(d, a, var.equal = FALSE)
  wil <- wilcox.test(d, a, exact = TRUE)
  tibble(n_dormant = length(d), n_active = length(a),
         mean_dormant = mean(d), mean_active = mean(a),
         diff = mean(d) - mean(a),
         ci_low = welch$conf.int[1], ci_high = welch$conf.int[2],
         t = unname(welch$statistic), df = unname(welch$parameter),
         p_welch = welch$p.value, W = unname(wil$statistic),
         p_wilcox = wil$p.value)
}

cor_both <- function(x, y) {
  p <- cor.test(x, y, method = "pearson")
  s <- suppressWarnings(cor.test(x, y, method = "spearman", exact = TRUE))
  tibble(n = sum(complete.cases(x, y)),
         pearson_r = unname(p$estimate), pearson_p = p$p.value,
         spearman_rho = unname(s$estimate), spearman_p = s$p.value)
}

# -----------------------------------------------------------------------------
# Figure style (plan §7; matches Figures/manuscript_v4, see D3)
# -----------------------------------------------------------------------------

group_colors <- c(Simple = "#0072B2", Polymer = "#E69F00",
                  Reference = "grey30", Summary = "black")
season_colors <- c(Active = "#009E73", Dormant = "#56B4E9")
dormant_shade <- "grey90"

theme_plan <- theme_bw(base_size = 9, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 10),
    plot.tag = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.key.size = unit(0.4, "cm")
  )

date_axis <- scale_x_date(date_breaks = "1 month", date_labels = "%b\n%y")

# Saves PDF and 300 dpi PNG at a journal width (85 or 180 mm) plus a caption
save_plan_figure <- function(plot, name, width_mm, height_mm, caption) {
  stopifnot(width_mm %in% c(85, 180))
  for (ext in c("png", "pdf")) {
    ggsave(file.path(dir_figures, paste0(name, ".", ext)), plot,
           width = width_mm, height = height_mm, units = "mm", dpi = 300,
           bg = "white")
  }
  writeLines(caption, file.path(dir_figures, paste0(name, ".txt")))
}

dormant_rect <- function() {
  annotate("rect", xmin = as.Date("2023-10-19"), xmax = as.Date("2024-03-20"),
           ymin = -Inf, ymax = Inf, fill = dormant_shade)
}
