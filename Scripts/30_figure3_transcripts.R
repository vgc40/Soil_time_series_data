# =============================================================================
# Figure 3: seasonal field transcript expression (Q2; H2)
# =============================================================================
#
# Q2: How does functional transcript expression differ between plant-active
#     and plant-dormant periods, and do functions associated with different
#     compound classes show distinct seasonal patterns?
# H2: Functional transcript composition differs between periods.
# The class-level part of Q2 is descriptive: only two transcript dates fall
# in the plant-dormant period, so no per-category tests are made.
#
# Data: BRITE C expression for 35 samples (5 plots x 7 dates), 325 categories
# after removing rows without a KEGG pathway code and curated human/non-soil
# pathways. Expression is converted to within-sample proportions (D5).
#
# Main figure
#   A. NMDS of BRITE C composition; exact date-level PERMANOVA (H2).
#   B. Date-mean z-scores of the 20 BRITE C categories uniquely mapped to one
#      compound class.
#   C. Share of total expression for each compound class over time.
# Supplement
#   S3. H2 sensitivity: metabolism-only categories (D9), values as provided
#       (D5), and transition date excluded (D7).

rm(list = ls())
source(file.path("Scripts", "helpers.R"))

suppressPackageStartupMessages({
  library(patchwork)
  library(ggrepel)
})

set.seed(20231023)
dirs <- output_dirs("Figure3")

class_levels <- c("Lipid-like", "Aromatic", "Amino sugar", "Amino acid")
unique_brite <- read_csv(input_files$unique_brite_c, show_col_types = FALSE) %>%
  mutate(compound_class = factor(compound_class, levels = class_levels))

brite <- load_brite_c()
stopifnot(
  n_distinct(brite$Sample_cat) == 35,
  n_distinct(brite$brite_C) == 325,
  !anyNA(brite$Sampling_Date),
  all(unique_brite$brite_C %in% brite$brite_C)
)

to_matrix <- function(d, value) {
  w <- d %>%
    select(Sample_cat, Sampling_Date, Plot_ID, brite_C, value = all_of(value)) %>%
    pivot_wider(names_from = brite_C, values_from = value) %>%
    arrange(Sampling_Date, Plot_ID)
  list(meta = select(w, Sample_cat, Sampling_Date, Plot_ID),
       matrix = as.matrix(select(w, -c(Sample_cat, Sampling_Date, Plot_ID))))
}

# -----------------------------------------------------------------------------
# A. Composition and date-level season test (H2)
# -----------------------------------------------------------------------------

primary <- to_matrix(brite, "proportion")
nmds <- metaMDS(primary$matrix, distance = "bray", k = 2, trymax = 200,
                autotransform = FALSE, trace = 0)

nmds_scores <- as_tibble(scores(nmds, display = "sites")) %>%
  bind_cols(primary$meta) %>%
  mutate(Season = assign_season(Sampling_Date))

date_centroids <- nmds_scores %>%
  group_by(Sampling_Date, Season) %>%
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = "drop") %>%
  arrange(Sampling_Date) %>%
  mutate(label = format(Sampling_Date, "%b %y"))

h2_test <- date_level_season_test(primary$matrix, primary$meta$Sampling_Date)
write_csv(h2_test, file.path(dirs$results, "Figure3A_season_PERMANOVA_dispersion.csv"))
write_csv(nmds_scores, file.path(dirs$results, "Figure3A_NMDS_sample_scores.csv"))
primary_h2 <- filter(h2_test, level == "Date means (primary)")
sample_h2 <- filter(h2_test, level != "Date means (primary)")

panel_a <- ggplot(nmds_scores, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour = Season), alpha = 0.4, size = 1.6) +
  geom_path(data = date_centroids, colour = "grey45", linewidth = 0.4,
            arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
  geom_point(data = date_centroids, aes(fill = Season), shape = 21,
             size = 3.4, colour = "black", stroke = 0.4) +
  geom_text_repel(data = date_centroids, aes(label = label), size = 2.6,
                  min.segment.length = 0, seed = 1) +
  scale_colour_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  labs(
    title = "BRITE C transcript composition (325 categories)",
    subtitle = sprintf(
      "Date-level PERMANOVA R² = %.2f, exact p = %s (smallest possible %.3f);\nsample-level p = %s; dispersion p = %s; stress = %.2f",
      primary_h2$R2, format_p(primary_h2$p_permanova), primary_h2$min_possible_p,
      format_p(sample_h2$p_permanova), format_p(primary_h2$p_dispersion), nmds$stress
    ),
    colour = NULL, fill = NULL
  ) +
  theme_manuscript +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 7.5))

# -----------------------------------------------------------------------------
# B. Uniquely mapped categories: date-mean z-scores
# -----------------------------------------------------------------------------

category_dates <- brite %>%
  inner_join(unique_brite, by = "brite_C") %>%
  group_by(compound_class, brite_C_name, brite_C, Sampling_Date) %>%
  summarise(mean_proportion = mean(proportion), .groups = "drop") %>%
  group_by(brite_C) %>%
  mutate(z = (mean_proportion - mean(mean_proportion)) / sd(mean_proportion)) %>%
  ungroup() %>%
  mutate(Season = assign_season(Sampling_Date))

category_summary <- brite %>%
  inner_join(unique_brite, by = "brite_C") %>%
  group_by(compound_class, brite_C_name, brite_C) %>%
  summarise(
    mean_share_pct = 100 * mean(proportion),
    n_samples_nonzero = sum(expression > 0),
    .groups = "drop"
  ) %>%
  left_join(
    category_dates %>%
      group_by(brite_C, Season) %>%
      summarise(m = mean(mean_proportion), .groups = "drop") %>%
      pivot_wider(names_from = Season, values_from = m) %>%
      transmute(brite_C,
                dormant_to_active_ratio = `Plant-dormant` / `Plant-active`),
    by = "brite_C"
  ) %>%
  mutate(note = "Descriptive: 5 active and 2 dormant dates; no per-category test.")

write_csv(category_dates, file.path(dirs$results, "Figure3B_unique_category_date_means.csv"))
write_csv(category_summary, file.path(dirs$results, "Figure3B_unique_category_summary.csv"))

row_order <- category_summary %>%
  arrange(compound_class, desc(mean_share_pct)) %>%
  mutate(row_label = sprintf("%s (%s%%)", brite_C_name, format(signif(mean_share_pct, 2), scientific = FALSE)))

heat_data <- category_dates %>%
  left_join(select(row_order, brite_C, row_label), by = "brite_C") %>%
  mutate(
    row_label = factor(row_label, levels = rev(row_order$row_label)),
    date_label = factor(format(Sampling_Date, "%b %y"),
                        levels = format(sort(unique(Sampling_Date)), "%b %y"))
  )

# Date labels coloured by season (dormant dates in bold).
date_seasons <- heat_data %>% distinct(date_label, Season) %>% arrange(date_label)

panel_b <- ggplot(heat_data, aes(date_label, row_label, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-2.3, 2.3), oob = scales::squish,
                       name = "z (within\ncategory)") +
  scale_x_discrete(position = "top") +
  facet_grid(compound_class ~ ., scales = "free_y", space = "free_y",
             labeller = label_wrap_gen(width = 6)) +
  labs(title = "Categories mapped uniquely to one compound class",
       subtitle = "Row label: mean share of total expression; dormant dates in bold",
       x = NULL, y = NULL) +
  theme_manuscript +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x.top = element_text(
      colour = season_colors[as.character(date_seasons$Season)],
      face = if_else(date_seasons$Season == "Plant-dormant", "bold", "plain")
    ),
    strip.text.y = element_text(size = 7, angle = 0),
    panel.grid = element_blank(), plot.subtitle = element_text(size = 7.5)
  )

# -----------------------------------------------------------------------------
# C. Class share of total expression (decision D8: sum of unique categories)
# -----------------------------------------------------------------------------

class_share <- brite %>%
  inner_join(unique_brite, by = "brite_C") %>%
  group_by(Sample_cat, Sampling_Date, Plot_ID, Season, compound_class) %>%
  summarise(share_pct = 100 * sum(proportion), .groups = "drop")

class_share_dates <- class_share %>%
  group_by(Sampling_Date, Season, compound_class) %>%
  summarise(
    mean = mean(share_pct), n = n(), se = sd(share_pct) / sqrt(n),
    lower = mean - qt(0.975, n - 1) * se, upper = mean + qt(0.975, n - 1) * se,
    .groups = "drop"
  )
write_csv(class_share, file.path(dirs$results, "Figure3C_class_share_by_sample.csv"))
write_csv(class_share_dates, file.path(dirs$results, "Figure3C_class_share_by_date.csv"))

panel_c <- ggplot(class_share_dates, aes(Sampling_Date, mean, colour = compound_class)) +
  season_band_layers(class_share$Sampling_Date) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = compound_class),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.5) +
  geom_point(data = class_share, aes(Sampling_Date, share_pct), size = 0.6,
             alpha = 0.35) +
  scale_colour_manual(values = class_colors[class_levels]) +
  scale_fill_manual(values = class_colors[class_levels]) +
  facet_wrap(~ compound_class, scales = "free_y", ncol = 2) +
  scale_x_date(breaks = sort(unique(class_share$Sampling_Date)),
               date_labels = "%b %y") +
  labs(title = "Class share of total transcript expression",
       x = NULL, y = "% of total BRITE C expression (mean, 95% CI)") +
  theme_manuscript +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))

figure3 <- (free(panel_a) | panel_c) / panel_b +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "A")
save_figure(figure3, dirs$figures, "Figure3_main", width = 11, height = 10.5)

# -----------------------------------------------------------------------------
# Supplement: H2 sensitivity
# -----------------------------------------------------------------------------

metabolism_codes <- brite %>%
  distinct(brite_C, pathway_code) %>%
  filter(as.integer(str_remove(pathway_code, "ko")) < 2000) %>%
  pull(brite_C)

metabolism_only <- brite %>%
  filter(brite_C %in% metabolism_codes) %>%
  group_by(Sample_cat) %>%
  mutate(proportion = expression / sum(expression)) %>%
  ungroup() %>%
  to_matrix("proportion")

as_provided <- to_matrix(brite, "expression")

no_transition <- brite %>%
  filter(!Sampling_Date %in% transition_dates) %>%
  to_matrix("proportion")

h2_sensitivity <- bind_rows(
  h2_test %>% mutate(analysis = "Primary: 325 categories, proportions"),
  date_level_season_test(metabolism_only$matrix, metabolism_only$meta$Sampling_Date) %>%
    mutate(analysis = sprintf("Metabolism only (%d categories, KEGG map < 02000)",
                              length(metabolism_codes))),
  date_level_season_test(as_provided$matrix, as_provided$meta$Sampling_Date) %>%
    mutate(analysis = "Values as provided (no within-sample normalization)"),
  date_level_season_test(no_transition$matrix, no_transition$meta$Sampling_Date) %>%
    mutate(analysis = "Transition date (Oct 2023) excluded")
)
write_csv(h2_sensitivity, file.path(dirs$results, "FigureS3_H2_sensitivity.csv"))

total_expression <- brite %>%
  group_by(Sample_cat, Sampling_Date, Season) %>%
  summarise(total = sum(expression), .groups = "drop")
write_csv(total_expression, file.path(dirs$results, "FigureS3_total_expression_by_sample.csv"))

figure_s3 <- h2_sensitivity %>%
  mutate(analysis = fct_rev(fct_inorder(analysis))) %>%
  ggplot(aes(R2, analysis, colour = level)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  geom_text(aes(label = paste("p =", format_p(p_permanova))), size = 2.6,
            position = position_dodge(width = 0.5), hjust = -0.25,
            show.legend = FALSE) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.35))) +
  scale_colour_manual(values = c("grey20", "grey60")) +
  labs(x = "PERMANOVA R² for season", y = NULL, colour = NULL) +
  theme_manuscript +
  theme(legend.position = "bottom")
save_figure(figure_s3, dirs$figures, "FigureS3_H2_sensitivity", width = 8, height = 3.5)

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

summary_table <- tribble(
  ~item, ~value,
  "Samples (plots x dates)", sprintf("%d (%d x %d)", nrow(primary$matrix),
                                     n_distinct(primary$meta$Plot_ID),
                                     n_distinct(primary$meta$Sampling_Date)),
  "Dates active / dormant", sprintf("%d / %d", primary_h2$n_active_dates,
                                    primary_h2$n_dormant_dates),
  "BRITE C categories", as.character(ncol(primary$matrix)),
  "Per-sample total expression range (as provided)",
  sprintf("%.0f to %.0f", min(total_expression$total), max(total_expression$total)),
  "NMDS stress", sprintf("%.3f", nmds$stress),
  "H2 date-level PERMANOVA R2", sprintf("%.3f", primary_h2$R2),
  "H2 date-level p (exact, 21 label sets)", format_p(primary_h2$p_permanova),
  "H2 sample-level p (labels permuted by date)", format_p(sample_h2$p_permanova),
  "H2 date-level dispersion p", format_p(primary_h2$p_dispersion),
  "Unique categories", as.character(nrow(unique_brite)),
  "Unique categories with dormant/active ratio < 1",
  as.character(sum(category_summary$dormant_to_active_ratio < 1))
)
write_csv(summary_table, file.path(dirs$results, "Figure3_analysis_summary.csv"))
print(summary_table, n = Inf)
print(h2_sensitivity %>% select(analysis, level, R2, p_permanova, p_dispersion))
print(category_summary %>% select(compound_class, brite_C_name, mean_share_pct,
                                  n_samples_nonzero, dormant_to_active_ratio), n = Inf)
