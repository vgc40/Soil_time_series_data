# =============================================================================
# 07. Q2 figures (plan §5.5)
# =============================================================================
#
# Main figures 4-6 and supplementary figures S5-S7, each saved as PDF and
# 300 dpi PNG with a caption file. Reads the tables written by script 06.

source(file.path("Scripts", "plan_v1", "00_setup.R"))
suppressPackageStartupMessages(library(patchwork))

read_res <- function(name) read_csv(file.path(dir_results, name), show_col_types = FALSE)

step_colors <- c(extracellular = "#E69F00", intracellular = "#0072B2")
nice_label <- function(x) {
  str_replace_all(x, c("_hydrochloride" = "", "_acid" = " acid", "_" = " ",
                       "^D\\." = "D-", "^L\\." = "L-",
                       "^N\\.Acetyl\\.D\\.glucosamine" = "NAG (compound)"))
}

# -----------------------------------------------------------------------------
# Figure 4. Quadrant scatter (H2b)
# -----------------------------------------------------------------------------

h2b_all <- read_res("cp4_h2b_compound_effects_all_variants.csv") |>
  filter(variant %in% c("main", "a: potential per MBC")) |>
  mutate(panel = if_else(variant == "main",
                         "A. Potential per mg C (main; D11)",
                         "B. Potential per µg MBC (plan original; robustness a)"),
         label = nice_label(set))
h2b_sum <- read_res("cp4_h2b_summary_and_robustness.csv") |>
  filter(!str_starts(variant, "b")) |>
  mutate(panel = if_else(str_starts(variant, "main"),
                         "A. Potential per mg C (main; D11)",
                         "B. Potential per µg MBC (plan original; robustness a)"),
         label = sprintf("%d of %d concordant, sign test p = %.2f\nSpearman rho = %.2f, p = %.2f",
                         n_concordant, n_compounds, binom_p, spearman_rho, spearman_p))

fig4 <- ggplot(h2b_all, aes(pot_effect, expr_effect)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_point(aes(colour = pathway_step), size = 1.8) +
  ggrepel::geom_text_repel(aes(label = label, colour = pathway_step), size = 2,
                           max.overlaps = Inf, min.segment.length = 0.2,
                           segment.size = 0.2, seed = 2026, show.legend = FALSE) +
  geom_text(data = h2b_sum, aes(x = -Inf, y = Inf, label = label), hjust = -0.03,
            vjust = 1.2, size = 2.3, lineheight = 0.9) +
  facet_wrap(~panel, scales = "free_x") +
  scale_colour_manual(values = step_colors, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.25))) +
  labs(x = "Potential season effect (dormant − active, z)",
       y = "Expression season effect (dormant − active, z)") +
  theme_plan +
  theme(legend.position = "bottom")

save_plan_figure(
  fig4, "Fig4_H2b_quadrant", 180, 95,
  paste0(
    "Figure 4. Seasonal concordance of expression and substrate-use potential ",
    "(H2b). One point per compound with an approved gene set (n = 19; 3 ",
    "extracellular polymer sets, 16 intracellular). Effects are dormant minus ",
    "active differences of date means of z-scores on the 7 transcript dates ",
    "(2 dormant, 5 active; 100% irrigation). Expression: centred log ratio of ",
    "the gene-set sum (35 samples, 5 plots per date). Potential: water-",
    "subtracted 5 h SIR (MicroResp samples on the same dates, 6 plots per ",
    "date). (A) Potential per mg C, the main version since deviation D11. ",
    "(B) Potential per µg MBC, the plan's original version (robustness a); the ",
    "effects are nearly identical across compounds because one sample (SM1625, ",
    "2024-05-08, MBC replaced with 0.07 µg C/g) dominates every compound's ",
    "z-scores. Concordant compounds fall in ",
    "the upper-right and lower-left quadrants. Sign test: two-sided binomial ",
    "test against 50%. No error bars."
  )
)

# -----------------------------------------------------------------------------
# Figure 5. Lag advantage (H2c)
# -----------------------------------------------------------------------------

lag <- read_res("cp4_h2c_lag_advantage.csv") |>
  mutate(label = nice_label(set),
         label = if_else(set_type == "enzyme", paste0(label, " (enzyme)"), label)) |>
  arrange(pathway_step, lag_advantage) |>
  mutate(label = factor(label, levels = unique(label)))
step_res <- read_res("cp4_h2c_step_permutation_and_robustness.csv") |>
  filter(str_starts(variant, "main"))
step_means <- lag |>
  group_by(pathway_step) |>
  summarise(mean = mean(lag_advantage), .groups = "drop")

fig5 <- ggplot(lag, aes(lag_advantage, label, colour = pathway_step)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(data = step_means, aes(xintercept = mean, colour = pathway_step),
             linewidth = 0.4, alpha = 0.6) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, orientation = "y") +
  geom_point(aes(shape = set_type), size = 1.8) +
  facet_grid(pathway_step ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = step_colors, guide = "none") +
  scale_shape_manual(values = c(compound = 16, enzyme = 17),
                     labels = c(compound = "SIR gene set", enzyme = "Enzyme gene set"),
                     name = NULL) +
  labs(x = "Lag advantage (r lag − r same, 95% bootstrap CI)", y = NULL,
       subtitle = sprintf("Extracellular − intracellular mean = %.2f\npermutation p = %.3f",
                          step_res$diff, step_res$p_perm)) +
  theme_plan +
  theme(legend.position = "bottom", strip.text.y = element_text(angle = 0),
        plot.subtitle = element_text(size = 8))

save_plan_figure(
  fig5, "Fig5_H2c_lag_advantage", 85, 140,
  paste0(
    "Figure 5. Lag advantage for each gene set (H2c). For each set, r_same is ",
    "the Pearson correlation between date-mean expression (centred log ratio) ",
    "and date-mean potential on the same date, and r_lag between expression at ",
    "a transcript date and potential at the next sampling date, both over 7 ",
    "date pairs with plots matched (5 plots per pair; 4 for 2024-01-02 to ",
    "2024-02-12 for SIR sets). Points: r_lag − r_same; bars: 95% percentile ",
    "bootstrap CI from 2,000 resamples of date pairs. Potential: water-",
    "subtracted 5 h SIR per mg C for compound sets (D11; n = 16 intracellular, 3 ",
    "extracellular) and measured activity for enzyme sets (n = 7, ",
    "extracellular). Vertical coloured lines: group means. Subtitle: ",
    "permutation test of the difference in means (10,000 label shuffles)."
  )
)

# -----------------------------------------------------------------------------
# Figure 6. Example pathway: chitin and NAG
# -----------------------------------------------------------------------------

expr_dates <- read_res("cp4_expression_clr_samples.csv") |>
  group_by(set, Collection_date) |>
  summarise(value = mean(clr), .groups = "drop")
pot_dates <- read_res("cp4_potential_samples_per_MBC.csv") |>
  group_by(set, Collection_date) |>
  summarise(value = mean(potential), .groups = "drop")
pot_perC_dates <- read_res("sir_prepared_long.csv") |>
  filter(Irrigation_trmt_pct == 100, Timepoint == 5,
         Substrate %in% c("Chitin", "N.Acetyl.D.glucosamine")) |>
  group_by(set = Substrate, Collection_date) |>
  summarise(value = mean(induced_C), .groups = "drop")

z_dates <- function(d) d |> group_by(series) |> mutate(z = (value - mean(value)) / sd(value)) |> ungroup()

panel_a <- bind_rows(
  expr_dates |> filter(set == "Chitin") |> mutate(series = "Chitinase expression (CLR)"),
  pot_dates |> filter(set == "Chitin") |> mutate(series = "Chitin potential (per µg MBC)"),
  pot_perC_dates |> filter(set == "Chitin") |> mutate(series = "Chitin potential (per mg C)"),
  pot_dates |> filter(set == "NAG") |> mutate(series = "NAG enzyme activity")
) |> z_dates()
panel_b <- bind_rows(
  expr_dates |> filter(set == "N.Acetyl.D.glucosamine") |>
    mutate(series = "NAG uptake + catabolism expression (CLR)"),
  pot_dates |> filter(set == "N.Acetyl.D.glucosamine") |>
    mutate(series = "NAG potential (per µg MBC)"),
  pot_perC_dates |> filter(set == "N.Acetyl.D.glucosamine") |>
    mutate(series = "NAG potential (per mg C)")
) |> z_dates()

series_colors <- c(
  "Chitinase expression (CLR)" = "#D55E00",
  "Chitin potential (per µg MBC)" = "#E69F00",
  "Chitin potential (per mg C)" = "#F0E442",
  "NAG enzyme activity" = "#009E73",
  "NAG uptake + catabolism expression (CLR)" = "#D55E00",
  "NAG potential (per µg MBC)" = "#0072B2",
  "NAG potential (per mg C)" = "#56B4E9"
)

pathway_panel <- function(d, tag) {
  ggplot(d, aes(Collection_date, z, colour = series)) +
    dormant_rect() +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 1.4) +
    scale_colour_manual(values = series_colors, name = NULL) +
    date_axis +
    labs(x = NULL, y = "z (date means)", tag = tag) +
    guides(colour = guide_legend(ncol = 2)) +
    theme_plan +
    theme(legend.position = "bottom")
}

fig6 <- pathway_panel(panel_a, "A") / pathway_panel(panel_b, "B")
save_plan_figure(
  fig6, "Fig6_example_chitin_NAG_pathway", 180, 150,
  paste0(
    "Figure 6. Example pathway: chitin and N-acetylglucosamine (NAG). Each ",
    "series is z-scored across its own date means for display. (A) Chitin ",
    "gene-set expression (chitinases and chitin-binding protein; centred log ",
    "ratio; 7 transcript dates, 5 plots each), chitin potential (water-",
    "subtracted 5 h SIR per mg C, the main version, and per µg MBC for comparison; 13 dates), and ",
    "NAG enzyme activity (13 dates). (B) NAG uptake and catabolism expression ",
    "(7 dates) with NAG potential per mg C and per µg MBC (13 dates). 100% ",
    "irrigation; 5-6 plots per date. The per-MBC series peak on 2024-02-12 and ",
    "2024-05-08 because of samples SM1539 and SM1625 (lowest and replaced MBC; ",
    "deviation log D1). Grey shading: plant-dormant period. No error bars."
  )
)

# -----------------------------------------------------------------------------
# Figure S5. H2a expression group scores by season
# -----------------------------------------------------------------------------

h2a_dates <- read_res("cp4_h2a_expression_scores_dates.csv") |>
  mutate(season = factor(season, levels = c("Active", "Dormant")),
         group = factor(group, levels = c("Simple", "Polymer")))
h2a <- read_res("cp4_h2a_expression_season.csv") |>
  filter(score %in% c("Simple", "Polymer")) |>
  mutate(group = factor(score, levels = c("Simple", "Polymer")),
         label = sprintf("Dormant − active = %.2f\n(95%% CI %.2f to %.2f)",
                         diff, ci_low, ci_high))

s5 <- ggplot(h2a_dates, aes(season, score)) +
  geom_point(aes(colour = group), size = 1.8,
             position = position_jitter(width = 0.08, seed = 2026)) +
  geom_text(data = h2a, aes(x = 0.5, y = Inf, label = label), hjust = 0, vjust = 1.2,
            size = 2.4, lineheight = 0.9) +
  facet_wrap(~group) +
  scale_colour_manual(values = group_colors, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.35))) +
  labs(x = NULL, y = "Expression score (mean z of set CLR)") +
  theme_plan

save_plan_figure(
  s5, "FigS5_H2a_expression_scores", 85, 75,
  paste0(
    "Figure S5. Expression scores by season (H2a, descriptive). Each gene ",
    "set's centred log ratio was z-scored across 35 samples and averaged within ",
    "group (simple: 14 sets; polymer: 3 sets: cellulose, hemicellulose, ",
    "chitin). Points: date means (7 dates: 5 active, 2 dormant; 5 plots each). ",
    "Annotations: dormant minus active difference with Welch 95% CI. No ",
    "p-values, as specified in the plan."
  )
)

# -----------------------------------------------------------------------------
# Figure S6. Bridge to Q1
# -----------------------------------------------------------------------------

bridge <- read_res("cp4_bridge_q1.csv") |>
  select(test, `7 transcript dates` = diff, ci_low, ci_high, `13 dates` = diff_13)
q1 <- read_res("cp2_t1_t3_t6_season_tests.csv") |>
  filter(normalization == "per_C", str_detect(test, "^[13]_")) |>
  select(test, diff, ci_low, ci_high) |>
  mutate(dates = "13 dates")
bridge_long <- bind_rows(
  read_res("cp4_bridge_q1.csv") |> select(test, diff, ci_low, ci_high) |>
    mutate(dates = "7 transcript dates"),
  q1
) |>
  mutate(test = recode(test, "1_H1a_simple_score" = "Test 1: simple score",
                       "3_H1b_polymer_minus_simple" = "Test 3: polymer − simple"))

s6 <- ggplot(bridge_long, aes(diff, test, colour = dates)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0, orientation = "y",
                position = position_dodge(width = 0.5)) +
  geom_point(size = 1.8, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("13 dates" = "black", "7 transcript dates" = "#CC79A7"),
                      name = NULL) +
  labs(x = "Dormant − active (z, Welch 95% CI)", y = NULL) +
  theme_plan +
  theme(legend.position = "bottom")

save_plan_figure(
  s6, "FigS6_bridge_Q1", 85, 60,
  paste0(
    "Figure S6. Bridge between Q1 and Q2 (test 5.3.1). Q1 tests 1 and 3 on ",
    "all 13 dates (8 active, 5 dormant; 77 samples) and on the 7 transcript ",
    "dates only (5 active, 2 dormant; 41 samples), using the same date-mean ",
    "scores (100% irrigation, 5 h, water-subtracted per C). Points: dormant ",
    "minus active difference; bars: Welch 95% CI."
  )
)

# -----------------------------------------------------------------------------
# Figure S7. Plot-level lag models
# -----------------------------------------------------------------------------

lm_cmp <- read_res("cp4_lag_models_compare.csv") |>
  mutate(label = nice_label(set))

s7 <- ggplot(lm_cmp, aes(marginal_R2_same, marginal_R2_lag, colour = pathway_step)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 1.8) +
  ggrepel::geom_text_repel(aes(label = label), size = 2, max.overlaps = Inf,
                           segment.size = 0.2, seed = 2026, show.legend = FALSE) +
  scale_colour_manual(values = step_colors, name = NULL) +
  labs(x = "Marginal R², potential(t) ~ expression(t)",
       y = "Marginal R², potential(t+1) ~ expression(t)") +
  theme_plan +
  theme(legend.position = "bottom")

save_plan_figure(
  s7, "FigS7_lag_models_plot_level", 85, 90,
  paste0(
    "Figure S7. Plot-level lag models (plan §5.4). For each gene set (n = 26), ",
    "linear mixed models with a random intercept for date, fitted by maximum ",
    "likelihood on plot-matched pairs (34 plot-dates for SIR sets, 35 for ",
    "enzyme sets; both models on the same rows). Response and predictor ",
    "z-scored. Points above the dashed 1:1 line have a higher marginal R² for ",
    "next-date potential than for same-date potential. AIC values are in ",
    "cp4_lag_models_compare.csv. No error bars."
  )
)
