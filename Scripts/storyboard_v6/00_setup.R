# =============================================================================
# Shared definitions for storyboard v6 (Figure brief v6)
# =============================================================================
#
# Sourced by every storyboard_v6 script. Reuses the plan_v1 inputs, season
# dates, and season test (Scripts/plan_v1/00_setup.R), and adds the v6
# substrate classes, timeline styling, and the results-table recorder.
# Run scripts from the repository root.

source(file.path("Scripts", "plan_v1", "00_setup.R"))
suppressPackageStartupMessages({
  library(patchwork)
})

dir_v6 <- file.path("Results", "storyboard_v6")
dir_v6_data <- file.path(dir_v6, "data")
dir_v6_parts <- file.path(dir_v6, "results_parts")
walk(c(dir_v6, dir_v6_data, dir_v6_parts), dir.create,
     showWarnings = FALSE, recursive = TRUE)

season_levels <- c("Active", "Dormant")

# -----------------------------------------------------------------------------
# Substrate classes (brief rule 8): every grouping comes from this file only.
# File row order is the heatmap row order; class order is first appearance.
# -----------------------------------------------------------------------------

substrate_classes <- read_csv(file.path("Data", "substrate_classes.csv"),
                              show_col_types = FALSE,
                              col_types = cols(.default = col_character())) |>
  mutate(tested_class = as.logical(tested_class),
         provisional = as.logical(provisional),
         across(c(eco_association, polymer_pair), \(x) replace_na(x, "")))
stopifnot(nrow(substrate_classes) == 31,
          !anyDuplicated(substrate_classes$substrate))

class_levels <- unique(substrate_classes$chem_class)
substrate_levels <- substrate_classes$substrate
tested_classes <- substrate_classes |>
  filter(tested_class) |>
  distinct(chem_class) |>
  pull()

# Plain names for figure labels, derived from the MicroResp substrate names
substrate_label <- function(x) {
  lab <- x |>
    str_replace("^D\\.|^L\\.", "") |>
    str_replace("N\\.Acetyl\\.D\\.glucosamine", "N-acetylglucosamine") |>
    str_replace("_hydrochloride", "") |>
    str_replace_all("_", " ") |>
    str_replace("Tween40", "Tween 40")
  if_else(lab %in% c("BSA", "GABA"), lab, str_to_sentence(lab))
}

# Polymer -> breakdown product pairs, read from the polymer_pair column
polymer_pairs <- substrate_classes |>
  filter(polymer_pair != "") |>
  mutate(role = if_else(str_detect(chem_class, "polymers"), "polymer", "product")) |>
  select(polymer_pair, role, substrate) |>
  pivot_wider(names_from = role, values_from = substrate)
stopifnot(nrow(polymer_pairs) == 5, !anyNA(polymer_pairs))

class_colors <- setNames(
  c("#8C510A", "#D08C2E", "#B2182B", "#6B4C9A", "#9BBB59",
    "#2166AC", "#1B9E77", "#7F7F7F", "#E7298A", "#4D4D4D"),
  class_levels
)

# -----------------------------------------------------------------------------
# Timeline styling (brief rules 2-3)
# -----------------------------------------------------------------------------

season_colors <- c(Active = "#2E8B57", Dormant = "#5B6F8C")
season_fill <- c(Active = "#E3F1E8", Dormant = "#E4E8EF")

# Bands run to the midpoints between the last date of one season and the
# first of the next, so every sampling date sits inside its own band.
season_bands <- tibble(
  season = c("Active", "Dormant", "Active"),
  start = as.Date(c("2023-05-24", "2023-10-18", "2024-03-20")),
  end = as.Date(c("2023-10-18", "2024-03-20", "2024-06-19"))
)

# Approximate plant calendar: exact dates were not supplied in the brief
plant_markers <- tibble(
  date = as.Date(c("2023-07-15", "2023-10-15", "2024-04-01")),
  label = c("Harvest", "Harvest", "Emergence")
)

timeline_layers <- function(marker_labels = TRUE, label_size = 2.2) {
  # One fixed-fill layer per season, so faceted panels and data fill scales work
  bands <- map(season_levels, \(s) {
    geom_rect(data = filter(season_bands, season == s),
              aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
              fill = season_fill[[s]], inherit.aes = FALSE)
  })
  layers <- c(bands, list(
    geom_vline(xintercept = plant_markers$date, linewidth = 0.3,
               colour = "grey35", linetype = "22")
  ))
  if (marker_labels) {
    layers <- c(layers, list(
      annotate("text", x = plant_markers$date, y = Inf, label = plant_markers$label,
               angle = 90, hjust = 1.05, vjust = -0.4, size = label_size,
               colour = "grey30")
    ))
  }
  layers
}

timeline_x <- scale_x_date(limits = range(season_bands$start, season_bands$end),
                           breaks = study_dates, date_labels = "%b\n%y",
                           expand = expansion(0))

# For small panels: every other sampling date labelled
timeline_x_small <- scale_x_date(limits = range(season_bands$start, season_bands$end),
                                 breaks = study_dates[c(1, 3, 5, 7, 9, 11, 13)],
                                 date_labels = "%b\n%y", expand = expansion(0))

theme_v6 <- theme_bw(base_size = 8, base_family = "sans") +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    strip.text = element_text(size = 7),
    plot.title = element_text(face = "bold", size = 8.5),
    plot.subtitle = element_text(size = 7),
    plot.tag = element_text(face = "bold", size = 10),
    axis.text = element_text(size = 6.5),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7.5),
    legend.key.size = unit(0.35, "cm")
  )

save_v6 <- function(plot, name, width_mm, height_mm) {
  # cairo_pdf keeps Unicode labels (superscripts, minus signs) in the PDF
  ggsave(file.path(dir_v6, paste0(name, ".png")), plot, width = width_mm,
         height = height_mm, units = "mm", dpi = 300, bg = "white")
  ggsave(file.path(dir_v6, paste0(name, ".pdf")), plot, width = width_mm,
         height = height_mm, units = "mm", bg = "white", device = cairo_pdf)
}

# -----------------------------------------------------------------------------
# Summaries
# -----------------------------------------------------------------------------

zscore <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

# Mean and t-based 95% CI of plot values within each date
date_mean_ci <- function(df, value, ...) {
  df |>
    filter(!is.na(.data[[value]])) |>
    group_by(Collection_date, season, ...) |>
    summarise(n = n(), mean = mean(.data[[value]]), sd = sd(.data[[value]]),
              .groups = "drop") |>
    mutate(half = if_else(n > 1, qt(0.975, pmax(n - 1, 1)) * sd / sqrt(n), NA_real_),
           ci_low = mean - half, ci_high = mean + half) |>
    select(-half)
}

# Dormant - active difference with Welch 95% CI only (descriptive)
welch_ci <- function(value, season) {
  d <- value[season == "Dormant"]
  a <- value[season == "Active"]
  w <- t.test(d, a)
  tibble(n_dormant = length(d), n_active = length(a), diff = mean(d) - mean(a),
         ci_low = w$conf.int[1], ci_high = w$conf.int[2])
}

fmt_p <- function(p) if_else(p < 0.001, "p < 0.001", sprintf("p = %.3f", p))

# -----------------------------------------------------------------------------
# Results table recorder (brief deliverable 2)
# -----------------------------------------------------------------------------

results_columns <- c("figure", "panel", "quantity", "estimate", "ci_low",
                     "ci_high", "p", "p_holm", "n_dates", "n_samples", "notes")

result_rows <- function(figure, panel, quantity, estimate, ci_low = NA,
                        ci_high = NA, p = NA, p_holm = NA, n_dates = NA,
                        n_samples = NA, notes = "") {
  tibble(figure, panel, quantity, estimate = as.numeric(estimate),
         ci_low = as.numeric(ci_low), ci_high = as.numeric(ci_high),
         p = as.numeric(p), p_holm = as.numeric(p_holm),
         n_dates = as.integer(n_dates), n_samples = as.integer(n_samples),
         notes = as.character(notes))
}

write_parts <- function(rows, figure) {
  stopifnot(identical(names(rows), results_columns))
  write_csv(rows, file.path(dir_v6_parts, paste0(figure, ".csv")), na = "")
}

# Rows for date means from date_mean_ci()
date_rows <- function(dm, figure, panel, quantity_prefix, notes = "") {
  result_rows(figure, panel,
              paste0(quantity_prefix, " | ", dm$Collection_date),
              dm$mean, dm$ci_low, dm$ci_high, n_dates = 1, n_samples = dm$n,
              notes = notes)
}

# Rows for a season_test() result
season_rows <- function(st, figure, panel, quantity, n_samples = NA,
                        p_holm = NA, notes = "") {
  result_rows(figure, panel, quantity, st$diff, st$ci_low, st$ci_high,
              st$p_welch, p_holm, st$n_dormant + st$n_active, n_samples,
              paste0("dormant - active of date means; Welch CI and p; Wilcoxon p = ",
                     signif(st$p_wilcox, 3),
                     sprintf("; n dormant = %d, n active = %d", st$n_dormant,
                             st$n_active),
                     if_else(notes == "", "", paste0("; ", notes))))
}
