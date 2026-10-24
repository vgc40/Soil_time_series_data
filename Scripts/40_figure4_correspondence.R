# =============================================================================
# Figure 4: correspondence between SIR and transcript expression by class (Q3)
# =============================================================================
#
# Q3: How do seasonal changes in SIR within each compound class relate to
#     transcript abundance of the BRITE C categories mapped uniquely to that
#     class, and does this differ among classes? Also explored: transcripts
#     versus SIR measured at the next sampling date (one-month lag).
# Q3 is descriptive: rank correlations are reported without p-values because
# they rest on seven sampling dates.
#
# Class values (decision D8)
#   SIR: mean over the class's substrates of within-substrate z-scores of
#        log2 SIR (5-h, carbon-normalized).
#   Transcripts: summed share of total expression of the class's unique
#        BRITE C categories, z-scored. Sensitivity: mean of within-category
#        z-scores.
# z-scores are computed over the samples entering each comparison, so both
# measurements are on a common dimensionless scale; magnitudes are not
# comparable, only timing and direction.
#
# Main figure
#   A. Concurrent date-mean trajectories of class SIR and class transcripts.
#   B. Spearman correlations, concurrent vs one-month lag, date and plot level.
#   C. Per-date divergence (SIR z - transcript z) by class.
# Supplement
#   S4a. Concurrent and lagged scatter plots.
#   S4b. Whole-profile Procrustes comparison with per-date residuals.

rm(list = ls())
source(file.path("Scripts", "helpers.R"))

suppressPackageStartupMessages({
  library(patchwork)
})

set.seed(20231023)
dirs <- output_dirs("Figure4")

class_levels <- c("Lipid-like", "Aromatic", "Amino sugar", "Amino acid")

substrate_classes <- read_csv(input_files$substrate_classes, show_col_types = FALSE) %>%
  filter(!is.na(compound_class))
unique_brite <- read_csv(input_files$unique_brite_c, show_col_types = FALSE)
lag_pairs <- read_csv(
  file.path("Results", "manuscript_v4", "sample_availability",
            "lag_pairs_transcript_to_next_microresp.csv"),
  show_col_types = FALSE
) %>%
  filter(!is.na(next_Sample_cat)) %>%
  select(Sample_cat, next_Sample_cat)

zscore <- function(x) (x - mean(x)) / sd(x)

# -----------------------------------------------------------------------------
# Class values per sample
# -----------------------------------------------------------------------------

sir <- load_microresp_100() %>%
  inner_join(select(substrate_classes, Substrate, compound_class), by = "Substrate")

brite <- load_brite_c()

transcript_class <- function(samples) {
  d <- brite %>%
    filter(Sample_cat %in% samples) %>%
    inner_join(select(unique_brite, brite_C, compound_class), by = "brite_C")
  summed <- d %>%
    group_by(Sample_cat, compound_class) %>%
    summarise(value = sum(proportion), .groups = "drop") %>%
    group_by(compound_class) %>%
    mutate(transcript_z = zscore(value)) %>%
    ungroup() %>%
    select(-value)
  mean_z <- d %>%
    group_by(brite_C) %>%
    mutate(z = zscore(proportion)) %>%
    group_by(Sample_cat, compound_class) %>%
    summarise(transcript_z_meanz = mean(z), .groups = "drop")
  left_join(summed, mean_z, by = c("Sample_cat", "compound_class"))
}

sir_class <- function(samples) {
  sir %>%
    filter(Sample_cat %in% samples) %>%
    group_by(Substrate) %>%
    mutate(z = zscore(log2(SIR))) %>%
    group_by(Sample_cat, Sampling_Date, Plot_ID, compound_class) %>%
    summarise(sir_z = mean(z), .groups = "drop")
}

transcript_samples <- unique(brite$Sample_cat)
stopifnot(all(transcript_samples %in% sir$Sample_cat))

concurrent <- sir_class(transcript_samples) %>%
  inner_join(transcript_class(transcript_samples),
             by = c("Sample_cat", "compound_class")) %>%
  mutate(comparison = "Concurrent")

lagged <- lag_pairs %>%
  inner_join(transcript_class(lag_pairs$Sample_cat), by = "Sample_cat",
             relationship = "many-to-many") %>%
  inner_join(
    sir_class(lag_pairs$next_Sample_cat) %>%
      rename(next_Sample_cat = Sample_cat, next_date = Sampling_Date),
    by = c("next_Sample_cat", "compound_class")
  ) %>%
  left_join(distinct(brite, Sample_cat, Sampling_Date), by = "Sample_cat") %>%
  mutate(comparison = "One-month lag")

stopifnot(nrow(concurrent) == 35 * 4, nrow(lagged) == 34 * 4)

paired <- bind_rows(concurrent, lagged) %>%
  mutate(compound_class = factor(compound_class, levels = class_levels),
         comparison = factor(comparison, levels = c("Concurrent", "One-month lag")))
write_csv(paired, file.path(dirs$results, "Figure4_class_values_by_sample.csv"))

# Date means indexed by the transcript sampling date.
paired_dates <- paired %>%
  group_by(comparison, compound_class, Sampling_Date) %>%
  summarise(across(c(sir_z, transcript_z, transcript_z_meanz), mean),
            n_plots = n(), .groups = "drop") %>%
  mutate(Season = assign_season(Sampling_Date),
         divergence = sir_z - transcript_z)
write_csv(paired_dates, file.path(dirs$results, "Figure4_class_values_by_date.csv"))

# -----------------------------------------------------------------------------
# Spearman correlations (descriptive)
# -----------------------------------------------------------------------------

spearman <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))

correlations <- bind_rows(
  paired_dates %>%
    group_by(comparison, compound_class) %>%
    summarise(
      `Date level, summed share` = spearman(sir_z, transcript_z),
      `Date level, mean z (sensitivity)` = spearman(sir_z, transcript_z_meanz),
      n = n(), .groups = "drop"
    ),
  paired %>%
    group_by(comparison, compound_class) %>%
    summarise(
      `Plot level, summed share` = spearman(sir_z, transcript_z),
      `Plot level, mean z (sensitivity)` = spearman(sir_z, transcript_z_meanz),
      n = n(), .groups = "drop"
    )
) %>%
  pivot_longer(starts_with(c("Date", "Plot")), names_to = "level",
               values_to = "rho") %>%
  filter(!is.na(rho)) %>%
  mutate(
    scale = if_else(str_detect(level, "Date"), "Date means (n = 7)",
                    "Plot-dates (n = 34-35)"),
    transcript_metric = if_else(str_detect(level, "sensitivity"),
                                "Mean of category z (sensitivity)",
                                "Summed class share (primary)")
  )
write_csv(correlations, file.path(dirs$results, "Figure4B_spearman_correlations.csv"))

# -----------------------------------------------------------------------------
# A. Concurrent trajectories
# -----------------------------------------------------------------------------

concurrent_dates <- paired_dates %>%
  filter(comparison == "Concurrent") %>%
  select(compound_class, Sampling_Date, Season, SIR = sir_z,
         Transcripts = transcript_z) %>%
  pivot_longer(c(SIR, Transcripts), names_to = "measure", values_to = "z")

rho_labels <- correlations %>%
  filter(comparison == "Concurrent", scale == "Date means (n = 7)",
         transcript_metric == "Summed class share (primary)") %>%
  mutate(label = sprintf("rho = %.2f", rho))

panel_a <- ggplot(concurrent_dates, aes(Sampling_Date, z, colour = measure)) +
  season_band_layers(concurrent_dates$Sampling_Date) +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.6) +
  geom_text(data = rho_labels, aes(x = min(concurrent_dates$Sampling_Date), y = Inf,
                                   label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1.4, size = 3) +
  facet_wrap(~ compound_class, nrow = 1) +
  scale_colour_manual(values = measure_colors) +
  scale_x_date(breaks = sort(unique(concurrent_dates$Sampling_Date)),
               date_labels = "%b %y") +
  labs(title = "Concurrent class trajectories at the seven transcript dates",
       x = NULL, y = "Date-mean z-score", colour = NULL) +
  theme_manuscript +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))

# -----------------------------------------------------------------------------
# B. Correlation summary
# -----------------------------------------------------------------------------

panel_b <- correlations %>%
  mutate(compound_class = factor(compound_class, levels = rev(class_levels))) %>%
  ggplot(aes(rho, compound_class, colour = comparison, shape = transcript_metric)) +
  geom_vline(xintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_point(size = 2.4, position = position_dodge(width = 0.6)) +
  facet_wrap(~ scale, nrow = 1) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_colour_manual(values = c("Concurrent" = "grey15", "One-month lag" = "#D08C2E")) +
  scale_shape_manual(values = c("Summed class share (primary)" = 16,
                                "Mean of category z (sensitivity)" = 1)) +
  labs(title = "SIR–transcript rank correlation by class",
       x = "Spearman rho (descriptive)", y = NULL, colour = NULL, shape = NULL) +
  theme_manuscript +
  theme(legend.position = "bottom", legend.box = "vertical")

# -----------------------------------------------------------------------------
# C. Per-date divergence
# -----------------------------------------------------------------------------

divergence <- paired_dates %>%
  filter(comparison == "Concurrent") %>%
  mutate(date_label = factor(format(Sampling_Date, "%b %y"),
                             levels = format(sort(unique(Sampling_Date)), "%b %y")),
         compound_class = factor(compound_class, levels = rev(class_levels)))
date_seasons <- divergence %>% distinct(date_label, Season) %>% arrange(date_label)

panel_c <- ggplot(divergence, aes(date_label, compound_class, fill = divergence)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f", divergence)), size = 2.8) +
  scale_fill_gradient2(low = "#2B6CB8", mid = "white", high = "#B8562B",
                       midpoint = 0, name = "SIR z −\ntranscript z") +
  labs(title = "Where SIR and transcripts diverge",
       subtitle = "Positive: SIR relatively higher than expression; dormant dates in bold",
       x = NULL, y = NULL) +
  theme_manuscript +
  theme(panel.grid = element_blank(), plot.subtitle = element_text(size = 7.5),
        axis.text.x = element_text(
          colour = season_colors[as.character(date_seasons$Season)],
          face = if_else(date_seasons$Season == "Plant-dormant", "bold", "plain")
        ))

figure4 <- panel_a / (panel_b | panel_c) +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = "A")
save_figure(figure4, dirs$figures, "Figure4_main", width = 11, height = 8)

# -----------------------------------------------------------------------------
# Supplement S4a: scatter plots
# -----------------------------------------------------------------------------

figure_s4a <- ggplot(paired_dates, aes(transcript_z, sir_z)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_vline(xintercept = 0, colour = "grey80") +
  geom_point(aes(colour = Season), size = 2.2) +
  geom_text(aes(label = format(Sampling_Date, "%b")), size = 2.3, vjust = -0.9) +
  facet_grid(comparison ~ compound_class) +
  scale_colour_manual(values = season_colors) +
  labs(x = "Transcript class z (date mean; transcript date)",
       y = "SIR class z (date mean; same or next date)", colour = NULL) +
  theme_manuscript +
  theme(legend.position = "bottom")
save_figure(figure_s4a, dirs$figures, "FigureS4a_class_scatter", width = 10, height = 5.5)

# -----------------------------------------------------------------------------
# Supplement S4b: whole-profile Procrustes with per-date residuals
# -----------------------------------------------------------------------------

sir_matrix <- load_microresp_100() %>%
  filter(Sample_cat %in% transcript_samples) %>%
  select(Sample_cat, Sampling_Date, Plot_ID, Substrate, SIR) %>%
  pivot_wider(names_from = Substrate, values_from = SIR) %>%
  arrange(Plot_ID, Sampling_Date)
transcript_matrix <- brite %>%
  select(Sample_cat, brite_C, proportion) %>%
  pivot_wider(names_from = brite_C, values_from = proportion)
transcript_matrix <- transcript_matrix[match(sir_matrix$Sample_cat,
                                             transcript_matrix$Sample_cat), ]

nmds_sir <- metaMDS(as.matrix(select(sir_matrix, -c(Sample_cat, Sampling_Date, Plot_ID))),
                    distance = "bray", k = 2, trymax = 200, autotransform = FALSE, trace = 0)
nmds_tx <- metaMDS(as.matrix(select(transcript_matrix, -Sample_cat)),
                   distance = "bray", k = 2, trymax = 200, autotransform = FALSE, trace = 0)

# Cyclic shifts of the date order within plot keep each plot's series intact.
procrustes_perm <- how(
  within = Within(type = "series", mirror = FALSE),
  blocks = sir_matrix$Plot_ID,
  nperm = 999
)
protest_fit <- protest(scores(nmds_sir, display = "sites"),
                       scores(nmds_tx, display = "sites"),
                       permutations = procrustes_perm)

procrustes_residuals <- sir_matrix %>%
  select(Sample_cat, Sampling_Date, Plot_ID) %>%
  mutate(residual = residuals(protest_fit), Season = assign_season(Sampling_Date))
residuals_by_date <- procrustes_residuals %>%
  group_by(Sampling_Date, Season) %>%
  summarise(mean_residual = mean(residual), se = sd(residual) / sqrt(n()),
            .groups = "drop")

procrustes_summary <- tibble(
  procrustes_correlation = protest_fit$t0,
  protest_p = protest_fit$signif,
  sum_of_squares = protest_fit$ss,
  n_samples = nrow(sir_matrix),
  permutations = "Cyclic shifts within plot (series, no mirroring); exploratory",
  n_distinct_permutations = protest_fit$permutations
)
write_csv(procrustes_summary, file.path(dirs$results, "FigureS4b_PROTEST.csv"))
write_csv(procrustes_residuals, file.path(dirs$results, "FigureS4b_procrustes_residuals.csv"))

figure_s4b <- ggplot(residuals_by_date, aes(Sampling_Date, mean_residual, colour = Season)) +
  geom_point(data = procrustes_residuals, aes(y = residual), alpha = 0.35, size = 1.2) +
  geom_pointrange(aes(ymin = mean_residual - se, ymax = mean_residual + se)) +
  scale_colour_manual(values = season_colors) +
  scale_x_date(breaks = sort(unique(residuals_by_date$Sampling_Date)),
               date_labels = "%b %y") +
  labs(
    title = sprintf("Procrustes r = %.2f, PROTEST p = %s (exploratory)",
                    protest_fit$t0, format_p(protest_fit$signif)),
    x = NULL, y = "Procrustes residual (higher = poorer agreement)", colour = NULL
  ) +
  theme_manuscript +
  theme(legend.position = "bottom")
save_figure(figure_s4b, dirs$figures, "FigureS4b_procrustes_by_date", width = 6.5, height = 4)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

summary_table <- correlations %>%
  filter(transcript_metric == "Summed class share (primary)") %>%
  select(comparison, compound_class, scale, rho, n) %>%
  arrange(compound_class, scale, comparison)
write_csv(summary_table, file.path(dirs$results, "Figure4_analysis_summary.csv"))
print(summary_table, n = Inf)
print(procrustes_summary)
print(divergence %>% select(compound_class, date_label, sir_z, transcript_z, divergence), n = Inf)
