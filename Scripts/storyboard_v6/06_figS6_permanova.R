# =============================================================================
# Figure S6. Overall transcript composition (reuses v4 Figure 3A)
# =============================================================================
#
# BRITE C expression (325 categories, within-sample proportions) for the 35
# samples; NMDS for display and the exact date-level PERMANOVA on date means
# from Scripts/helpers.R. Run in its own session because helpers.R and the
# plan_v1 setup define functions with the same names.

source(file.path("Scripts", "helpers.R"))
suppressPackageStartupMessages(library(ggrepel))
set.seed(20231023)

dir_v6 <- file.path("Results", "storyboard_v6")
brite <- load_brite_c()
stopifnot(n_distinct(brite$Sample_cat) == 35, n_distinct(brite$brite_C) == 325)

w <- brite |>
  select(Sample_cat, Sampling_Date, Plot_ID, brite_C, proportion) |>
  pivot_wider(names_from = brite_C, values_from = proportion) |>
  arrange(Sampling_Date, Plot_ID)
meta <- select(w, Sample_cat, Sampling_Date, Plot_ID)
profile <- as.matrix(select(w, -c(Sample_cat, Sampling_Date, Plot_ID)))

test <- date_level_season_test(profile, meta$Sampling_Date)
primary <- filter(test, level == "Date means (primary)")

nmds <- metaMDS(profile, distance = "bray", k = 2, trymax = 200,
                autotransform = FALSE, trace = 0)
scores_df <- as_tibble(scores(nmds, display = "sites")) |>
  bind_cols(meta) |>
  mutate(Season = assign_season(Sampling_Date))
centroids <- scores_df |>
  group_by(Sampling_Date, Season) |>
  summarise(NMDS1 = mean(NMDS1), NMDS2 = mean(NMDS2), .groups = "drop") |>
  mutate(label = format(Sampling_Date, "%b %y"))

v6_colors <- c("Plant-active" = "#2E8B57", "Plant-dormant" = "#5B6F8C")
s6 <- ggplot(scores_df, aes(NMDS1, NMDS2)) +
  geom_point(aes(colour = Season), alpha = 0.45, size = 1.4) +
  geom_point(data = centroids, aes(fill = Season), shape = 21, size = 3,
             colour = "black", stroke = 0.4) +
  geom_text_repel(data = centroids, aes(label = label), size = 2.3, seed = 1) +
  scale_colour_manual(values = v6_colors, name = NULL) +
  scale_fill_manual(values = v6_colors, name = NULL) +
  labs(title = "Overall transcript composition (325 BRITE C categories)",
       subtitle = sprintf(
         "PERMANOVA on 7 date means (exact, 21 label sets): R² = %.2f, F = %.2f, p = %.4f\n(smallest possible p = 1/21 = %.4f); dispersion p = %.3f. Small points: samples; large: date means.",
         primary$R2, primary$F, primary$p_permanova, primary$min_possible_p,
         primary$p_dispersion),
       caption = sprintf("NMDS (Bray-Curtis, k = 2) for display only; stress = %.3f.", nmds$stress)) +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 8.5),
        plot.subtitle = element_text(size = 6.5))
ggsave(file.path(dir_v6, "FigS6.png"), s6, width = 120, height = 115, units = "mm",
       dpi = 300, bg = "white")
ggsave(file.path(dir_v6, "FigS6.pdf"), s6, width = 120, height = 115, units = "mm",
       bg = "white", device = cairo_pdf)

rows <- tibble(
  figure = "FigS6", panel = "",
  quantity = c("PERMANOVA R2, season, date means", "PERMANOVA dispersion F, date means",
               "PERMANOVA R2, samples with labels permuted by date"),
  estimate = c(primary$R2, primary$F_dispersion, test$R2[test$level != "Date means (primary)"]),
  ci_low = NA_real_, ci_high = NA_real_,
  p = c(primary$p_permanova, primary$p_dispersion,
        test$p_permanova[test$level != "Date means (primary)"]),
  p_holm = NA_real_, n_dates = 7L, n_samples = 35L,
  notes = c(sprintf("exact permutation over %d label sets; F = %.3f; p is the smallest possible",
                    primary$n_label_sets, primary$F),
            "betadisper on date means", "sensitivity (not shown in the figure)")
)
write_csv(rows, file.path(dir_v6, "results_parts", "FigS6.csv"), na = "")
print(test)
