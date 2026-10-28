# =============================================================================
# Figure 2 (Q1) and the MicroResp supplements S2, S3, S4, S5, S8, S9
# =============================================================================
#
# Reads data/microresp_z.csv (01_step0_data.R). Dates are the replicates:
# plots are averaged within each date before any season test (rule 1), and
# season tests use season_test() from plan_v1 (Welch CI and p; Wilcoxon p in
# the notes).

source(file.path("Scripts", "storyboard_v6", "00_setup.R"))
suppressPackageStartupMessages({
  library(ggh4x)
})

z_all <- read_csv(file.path(dir_v6_data, "microresp_z.csv"), show_col_types = FALSE) |>
  mutate(season = factor(season, season_levels),
         chem_class = factor(chem_class, class_levels),
         Substrate = factor(Substrate, substrate_levels))
get_z <- function(v) filter(z_all, version == v)

eco <- substrate_classes |>
  transmute(Substrate = factor(substrate, substrate_levels),
            chem_class = factor(chem_class, class_levels), eco_association)

# -----------------------------------------------------------------------------
# Shared builders
# -----------------------------------------------------------------------------

# Overall level: per sample, the mean z over all substrates; then date means
sample_level <- function(z) distinct(z, Sample_cat, Collection_date, Plot, season, level)

level_dates <- function(z) date_mean_ci(sample_level(z), "level")

timeline_panel <- function(pts, dm, value, title, ylab, colour = "black") {
  ggplot(dm, aes(Collection_date, mean)) +
    timeline_layers() +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey50") +
    geom_point(data = pts, aes(y = .data[[value]]), colour = "grey60", alpha = 0.6,
               size = 0.7) +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 6, linewidth = 0.35,
                  colour = colour) +
    geom_line(linewidth = 0.45, colour = colour) +
    geom_point(aes(fill = season), shape = 21, size = 1.8, colour = colour) +
    scale_fill_manual(values = season_colors, guide = "none") +
    timeline_x +
    labs(title = title, x = NULL, y = ylab) +
    theme_v6
}

# Season comparison of date means: dates as points, season mean with t CI,
# and the dormant - active difference in the subtitle
season_panel <- function(dm, title, ylab, st) {
  summ <- dm |>
    group_by(season) |>
    summarise(n = n(), m = mean(mean), s = sd(mean), .groups = "drop") |>
    mutate(half = qt(0.975, n - 1) * s / sqrt(n))
  ggplot(dm, aes(season, mean)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey50") +
    geom_point(aes(colour = season), position = position_jitter(width = 0.08, seed = 1),
               size = 1.6, alpha = 0.85) +
    geom_errorbar(data = summ, aes(y = m, ymin = m - half, ymax = m + half),
                  width = 0.12, linewidth = 0.45, position = position_nudge(x = 0.28)) +
    geom_point(data = summ, aes(y = m), size = 2.2, position = position_nudge(x = 0.28)) +
    scale_colour_manual(values = season_colors, guide = "none") +
    scale_x_discrete(labels = c(Active = sprintf("Active\n(%d dates)", st$n_active),
                                Dormant = sprintf("Dormant\n(%d dates)", st$n_dormant))) +
    labs(title = title,
         subtitle = sprintf("Dormant − active = %.2f\n95%% CI %.2f to %.2f\n%s",
                            st$diff, st$ci_low, st$ci_high, fmt_p(st$p_welch)),
         x = NULL, y = ylab) +
    theme_v6
}

# Heatmap of date means, rows grouped by chem_class with eco_association
# annotation; columns split into season runs with coloured strips.
season_run <- function(date) {
  factor(case_when(date < as.Date("2023-10-18") ~ "Plant-active",
                   date < as.Date("2024-03-20") ~ "Plant-dormant",
                   .default = "Plant-active "),
         levels = c("Plant-active", "Plant-dormant", "Plant-active "))
}

heatmap_panel <- function(cells, fill_label, title) {
  lim <- max(abs(cells$mean), na.rm = TRUE)
  markers <- plant_markers |> mutate(run = season_run(date))
  strips <- strip_themed(
    background_x = elem_list_rect(fill = season_fill[c("Active", "Dormant", "Active")]),
    background_y = elem_list_rect(fill = "grey95"),
    text_y = elem_list_text(angle = 0, hjust = 1, size = 5.5)
  )
  heat <- ggplot(cells, aes(Collection_date, fct_rev(Substrate), fill = mean)) +
    geom_tile(width = 27, colour = "white", linewidth = 0.2) +
    geom_vline(data = markers, aes(xintercept = date), linewidth = 0.3,
               colour = "grey30", linetype = "22") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-lim, lim), name = fill_label) +
    scale_x_date(breaks = study_dates, date_labels = "%b\n%y",
                 expand = expansion(add = 14)) +
    scale_y_discrete(labels = \(x) substrate_label(x)) +
    facet_grid2(chem_class ~ run, scales = "free", space = "free", switch = "y",
                strip = strips, labeller = labeller(chem_class = label_wrap_gen(30))) +
    labs(title = title, x = NULL, y = NULL) +
    theme_v6 +
    theme(strip.placement = "outside", panel.spacing.y = unit(0.6, "mm"),
          panel.spacing.x = unit(0.6, "mm"), legend.position = "bottom",
          legend.key.width = unit(1, "cm"), axis.text.y = element_text(size = 6))
  eco_col <- ggplot(eco |> filter(Substrate %in% cells$Substrate),
                    aes(1, fct_rev(Substrate), label = eco_association)) +
    geom_text(hjust = 0, size = 1.8, colour = "grey25") +
    scale_x_continuous(limits = c(1, 2.6), expand = expansion(0)) +
    facet_grid(chem_class ~ "Ecological\nassociation", scales = "free_y",
               space = "free_y") +
    labs(x = NULL, y = NULL) +
    theme_void(base_size = 7) +
    theme(strip.text.y = element_blank(), strip.text.x = element_text(size = 6),
          panel.spacing.y = unit(0.6, "mm"), plot.margin = margin(t = 13, b = 22))
  heat + eco_col + plot_layout(widths = c(5, 1.9))
}

heat_cells <- function(z, value) {
  date_mean_ci(z, value, Substrate, chem_class) |>
    mutate(run = season_run(Collection_date))
}

cell_rows <- function(cells, figure, panel, what, notes = "") {
  result_rows(figure, panel,
              paste0(what, " | ", cells$Substrate, " | ", cells$Collection_date),
              cells$mean, cells$ci_low, cells$ci_high, n_dates = 1, n_samples = cells$n,
              notes = notes)
}

# =============================================================================
# Figure 2
# =============================================================================

z5 <- get_z("primary_5h")
n_samples_5h <- n_distinct(z5$Sample_cat)

# A-B. Overall level
lvl_pts <- sample_level(z5)
lvl_dm <- level_dates(z5)
lvl_test <- season_test(lvl_dm$mean, lvl_dm$season)

panel_a <- timeline_panel(lvl_pts, lvl_dm, "level",
                          "Overall level: mean z-score across all 31 substrates",
                          "Mean z-score")
panel_b <- season_panel(lvl_dm, "Season comparison",
                        "Date mean of overall level (z)", lvl_test)

# C. Relative response heatmap
rel_cells <- heat_cells(z5, "relative")
panel_c <- heatmap_panel(rel_cells, "Relative response\n(z − sample mean z)",
                         "Relative response by substrate and date (date means)")

# D. Class trajectories of relative response
class_samples <- z5 |>
  group_by(Sample_cat, Collection_date, Plot, season, chem_class) |>
  summarise(relative = mean(relative), n_substrates = n(), .groups = "drop")
class_dm <- date_mean_ci(class_samples, "relative", chem_class)

class_tests <- class_dm |>
  filter(chem_class %in% tested_classes) |>
  group_by(chem_class) |>
  group_modify(\(d, k) season_test(d$mean, d$season)) |>
  ungroup() |>
  mutate(p_holm = p.adjust(p_welch, "holm"))

class_labels <- tibble(chem_class = factor(class_levels, class_levels)) |>
  left_join(class_tests, by = "chem_class") |>
  mutate(label = if_else(
    is.na(diff), "descriptive",
    sprintf("D − A = %.2f (%.2f to %.2f)\n%s; Holm %s", diff, ci_low, ci_high,
            fmt_p(p_welch), str_remove(fmt_p(p_holm), "^p "))
  ))

panel_d <- ggplot(class_dm, aes(Collection_date, mean, colour = chem_class)) +
  timeline_layers(marker_labels = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey50") +
  geom_point(data = class_samples, aes(y = relative), colour = "grey65", size = 0.4,
             alpha = 0.6) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  geom_text(data = class_labels, aes(x = as.Date("2023-05-28"), y = Inf, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1.15, size = 1.8, lineheight = 0.9) +
  scale_colour_manual(values = class_colors, guide = "none") +
  timeline_x_small +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.45))) +
  facet_wrap(~ chem_class, ncol = 5, labeller = label_wrap_gen(24)) +
  labs(title = "Class trajectories of relative response (mean over substrates in class; date means)",
       x = NULL, y = "Relative response") +
  theme_v6 +
  theme(axis.text.x = element_text(size = 5.5))

# E. Polymer vs breakdown product
pair_samples <- polymer_pairs |>
  mutate(data = map2(polymer, product, \(pol, pro) {
    z5 |>
      filter(Substrate %in% c(pol, pro)) |>
      mutate(role = if_else(Substrate == pol, "polymer", "product")) |>
      select(Sample_cat, Collection_date, season, role, z) |>
      pivot_wider(names_from = role, values_from = z) |>
      mutate(diff = polymer - product)
  })) |>
  select(polymer_pair, data) |>
  unnest(data)
pair_dm <- pair_samples |>
  group_by(polymer_pair, Collection_date, season) |>
  summarise(diff = mean(diff), n = n(), .groups = "drop")
pair_ci <- pair_dm |>
  group_by(polymer_pair) |>
  group_modify(\(d, k) welch_ci(d$diff, d$season)) |>
  ungroup()
pooled_dm <- pair_dm |>
  group_by(Collection_date, season) |>
  summarise(diff = mean(diff), .groups = "drop")
pooled_test <- season_test(pooled_dm$diff, pooled_dm$season)

pair_plot <- bind_rows(
  pair_ci |> mutate(row = polymer_pair, kind = "Pair (descriptive)"),
  pooled_test |> transmute(n_dormant, n_active, diff, ci_low, ci_high,
                           row = "Pooled (5 pairs)", kind = "Pooled test")
) |>
  mutate(row = factor(row, rev(c(polymer_pairs$polymer_pair, "Pooled (5 pairs)"))))

panel_e <- ggplot(pair_plot, aes(diff, row, colour = kind)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey50") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), width = 0.2, orientation = "y", linewidth = 0.45) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = c("Pair (descriptive)" = "grey30", "Pooled test" = "black"),
                      guide = "none") +
  labs(title = "Polymer minus breakdown product",
       subtitle = sprintf("Season shift (D − A) of polymer z − product z, 95%% CI\nPooled: %.2f (%.2f to %.2f), %s",
                          pooled_test$diff, pooled_test$ci_low, pooled_test$ci_high,
                          fmt_p(pooled_test$p_welch)),
       x = "Dormant − active (z units)", y = NULL) +
  theme_v6

fig2 <- (panel_a + panel_b + plot_layout(widths = c(3, 1))) /
  wrap_elements(full = panel_c) /
  panel_d /
  (panel_e + plot_spacer() + plot_layout(widths = c(1.2, 1))) +
  plot_layout(heights = c(0.8, 3.2, 1.4, 0.55)) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "D", "E")))
save_v6(fig2, "Fig2", 180, 320)

# Results rows
fig2_rows <- bind_rows(
  date_rows(lvl_dm, "Fig2", "A", "overall level (mean z, 31 substrates) date mean"),
  season_rows(lvl_test, "Fig2", "B", "overall level dormant - active",
              n_samples = n_samples_5h),
  cell_rows(rel_cells, "Fig2", "C", "relative response date mean"),
  result_rows("Fig2", "D", paste0("class relative response date mean | ",
                                  class_dm$chem_class, " | ", class_dm$Collection_date),
              class_dm$mean, class_dm$ci_low, class_dm$ci_high, n_dates = 1,
              n_samples = class_dm$n),
  pmap_dfr(class_tests, \(...) {
    st <- tibble(...)
    season_rows(st, "Fig2", "D", paste0("class relative response dormant - active | ",
                                        st$chem_class),
                n_samples = n_samples_5h, p_holm = st$p_holm,
                notes = "Holm over the 4 tested classes")
  }),
  result_rows("Fig2", "D", paste0("class relative response | ",
                                  setdiff(class_levels, tested_classes)),
              NA, notes = "descriptive (fewer than 3 substrates); no test"),
  result_rows("Fig2", "E", paste0("polymer z - product z, dormant - active | ",
                                  pair_ci$polymer_pair),
              pair_ci$diff, pair_ci$ci_low, pair_ci$ci_high,
              n_dates = pair_ci$n_dormant + pair_ci$n_active, n_samples = n_samples_5h,
              notes = "descriptive; Welch CI on date means"),
  season_rows(pooled_test, "Fig2", "E",
              "polymer z - product z, dormant - active | pooled over 5 pairs",
              n_samples = n_samples_5h,
              notes = "per date: mean over the 5 pairs of (polymer z - product z)")
)
write_parts(fig2_rows, "Fig2")

# =============================================================================
# Figure S2. Each substrate over the 13 dates (z, level included)
# =============================================================================

sub_dm <- heat_cells(z5, "z")
s2 <- ggplot(sub_dm, aes(Collection_date, mean, colour = chem_class)) +
  timeline_layers(marker_labels = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey50") +
  geom_point(data = z5, aes(y = z), colour = "grey65", size = 0.35, alpha = 0.6) +
  geom_line(linewidth = 0.45) +
  geom_point(size = 0.8) +
  scale_colour_manual(values = class_colors, name = "Chemical class") +
  timeline_x_small +
  facet_wrap(~ Substrate, ncol = 6, labeller = as_labeller(substrate_label)) +
  labs(title = "Each substrate through the year (z-score, level included; date means, plots in grey)",
       subtitle = "Panels are ordered and coloured by chemical class (Data/substrate_classes.csv)",
       x = NULL, y = "z-score") +
  theme_v6 +
  theme(legend.position = "bottom", axis.text.x = element_text(size = 5)) +
  guides(colour = guide_legend(ncol = 3, title.position = "top"))
save_v6(s2, "FigS2", 180, 230)
write_parts(cell_rows(sub_dm, "FigS2", "", "z date mean"), "FigS2")

# =============================================================================
# Figure S3. Heatmap of raw z-scores (level included)
# =============================================================================

s3 <- heatmap_panel(sub_dm, "z-score\n(level included)",
                    "z-score by substrate and date (date means; level included)")
save_v6(s3, "FigS3", 180, 220)
write_parts(cell_rows(sub_dm, "FigS3", "", "z date mean"), "FigS3")

# =============================================================================
# Figure S4. Per-biomass version of Figure 2A-B
# =============================================================================

# Round 2: date-level version. For each substrate and date, the date-mean
# response per mg C (all MicroResp plots) is divided by the date-mean MBC
# (plots with MBC after the Step 0 rule). Each substrate is z-scored across
# the 13 dates; the overall level is the mean z across substrates per date.
# A sample-level per-MBC z is dominated by SM1420 (MBC 7.78).

mr5 <- read_csv(file.path(dir_v6_data, "microresp_long.csv"), show_col_types = FALSE) |>
  filter(Timepoint == 5) |>
  mutate(season = factor(season, season_levels))

mbc_date_level <- function(drop = character()) {
  d <- filter(mr5, !Sample_cat %in% drop)
  mbc <- d |>
    distinct(Sample_cat, Collection_date, MBC_used) |>
    group_by(Collection_date) |>
    summarise(mbc = mean(MBC_used, na.rm = TRUE), n_mbc = sum(!is.na(MBC_used)),
              .groups = "drop")
  d |>
    group_by(Substrate, Collection_date, season) |>
    summarise(resp = mean(primary), n = n(), .groups = "drop") |>
    left_join(mbc, by = "Collection_date") |>
    mutate(value = resp / mbc) |>
    group_by(Substrate) |>
    mutate(z = zscore(value)) |>
    group_by(Collection_date, season) |>
    summarise(mean = mean(z), n = first(n), n_mbc = first(n_mbc), .groups = "drop") |>
    mutate(ci_low = NA_real_, ci_high = NA_real_)
}

# Fig 2B pipeline (sample-level z per substrate) without some samples
primary_level <- function(drop = character()) {
  mr5 |>
    filter(!Sample_cat %in% drop) |>
    group_by(Substrate) |>
    mutate(z = zscore(primary)) |>
    group_by(Sample_cat, Collection_date, season) |>
    summarise(level = mean(z), .groups = "drop") |>
    date_mean_ci("level")
}

mbc_dm <- mbc_date_level()
mbc_test <- season_test(mbc_dm$mean, mbc_dm$season)
mbc_dm_no1420 <- mbc_date_level("SM1420")
mbc_test_no1420 <- season_test(mbc_dm_no1420$mean, mbc_dm_no1420$season)
prim_dm_no1420 <- primary_level("SM1420")
prim_test_no1420 <- season_test(prim_dm_no1420$mean, prim_dm_no1420$season)

no_points <- tibble(Collection_date = as.Date(character()), level = numeric())
s4 <- (timeline_panel(no_points, mbc_dm, "level",
                      "Overall level per unit biomass (date level): mean z across 31 substrates",
                      "Mean z-score (across 13 dates)") +
         season_panel(mbc_dm, "Season comparison", "Overall level (z)", mbc_test) +
         plot_layout(widths = c(3, 1))) /
  (season_panel(prim_dm_no1420, "Without SM1420:\nprimary (per mg C)",
                "Date mean of overall level (z)", prim_test_no1420) +
     season_panel(mbc_dm_no1420, "Without SM1420:\nper biomass (date level)",
                  "Overall level (z)", mbc_test_no1420) +
     plot_spacer() + plot_layout(widths = c(1, 1, 1))) +
  plot_annotation(tag_levels = list(c("A", "B", "C", "")),
                  caption = "Per biomass: date-mean response per mg C added / date-mean MBC (SM1625 and SM1539 MBC set to missing, Step 0),\nz-scored per substrate across the 13 dates. No plot-level points or CIs in A: values exist only at the date level.")
save_v6(s4, "FigS4", 180, 150)
write_parts(bind_rows(
  result_rows("FigS4", "A", paste0("overall level per MBC (date level) | ", mbc_dm$Collection_date),
              mbc_dm$mean, n_dates = 1, n_samples = mbc_dm$n,
              notes = paste0("mean z across 31 substrates; date-mean MBC from ",
                             mbc_dm$n_mbc, " plots")),
  season_rows(mbc_test, "FigS4", "B", "overall level per MBC (date level) dormant - active",
              n_samples = n_distinct(mr5$Sample_cat),
              notes = "SM1625 and SM1539 MBC missing (Step 0)"),
  season_rows(prim_test_no1420, "FigS4", "C",
              "overall level (primary, per mg C) dormant - active, without SM1420",
              n_samples = n_distinct(mr5$Sample_cat) - 1,
              notes = "Fig 2B pipeline with SM1420 removed before z-scoring"),
  season_rows(mbc_test_no1420, "FigS4", "C",
              "overall level per MBC (date level) dormant - active, without SM1420",
              n_samples = n_distinct(mr5$Sample_cat) - 1,
              notes = "SM1625 and SM1539 MBC missing (Step 0)")
), "FigS4")

# -----------------------------------------------------------------------------
# Round 2 sensitivity rows for Fig 2B (no new figure)
# -----------------------------------------------------------------------------

no_pea_dm <- z5 |>
  filter(Substrate != "Phenylethylamine") |>
  group_by(Sample_cat, Collection_date, Plot, season) |>
  summarise(level = mean(z), n_substrates = n(), .groups = "drop") |>
  date_mean_ci("level")
stopifnot(all(count(filter(z5, Substrate != "Phenylethylamine"), Sample_cat)$n == 30))
no_pea_test <- season_test(no_pea_dm$mean, no_pea_dm$season)

write_parts(bind_rows(
  season_rows(prim_test_no1420, "Fig2", "B",
              "overall level dormant - active, without SM1420 (sensitivity)",
              n_samples = n_samples_5h - 1,
              notes = "not drawn in Fig 2; same value as Fig S4C primary"),
  season_rows(no_pea_test, "Fig2", "B",
              "overall level dormant - active, without phenylethylamine (30 substrates; sensitivity)",
              n_samples = n_samples_5h, notes = "not drawn in any figure")
), "Fig2_sensitivity")

# =============================================================================
# Figure S5. Figure 2B with the 5-24 h reading
# =============================================================================

z24 <- get_z("primary_24h")
dm24 <- level_dates(z24)
test24 <- season_test(dm24$mean, dm24$season)
s5 <- season_panel(dm24, "Overall level, 5–24 h reading", "Date mean of overall level (z)",
                   test24) +
  labs(caption = "Same as Figure 2B, using the 5–24 h reading (water subtracted, per mg C added).")
save_v6(s5, "FigS5", 85, 85)
write_parts(bind_rows(
  date_rows(dm24, "FigS5", "", "overall level 5-24 h date mean"),
  season_rows(test24, "FigS5", "", "overall level 5-24 h dormant - active",
              n_samples = n_distinct(z24$Sample_cat))
), "FigS5")

# =============================================================================
# Figures S8 (assay QC) and S9 (water correction)
# =============================================================================

water <- read_csv(file.path(dir_v6_data, "water_control_5h.csv"), show_col_types = FALSE) |>
  mutate(season = factor(season, season_levels))
water_dm <- date_mean_ci(water, "SIR_gdw")
water_test <- season_test(water_dm$mean, water_dm$season)
water_title <- "Water-only (basal) CO2, 5 h"
water_ylab <- "ng CO2-C g⁻¹ dry soil h⁻¹"

qc <- read_csv(file.path(dir_v6_data, "step0_assay_qc.csv"), show_col_types = FALSE) |>
  mutate(season = factor(season, season_levels))
s8_b <- ggplot(qc, aes(Collection_date, glucose_mgC_per_gdw_mean)) +
  timeline_layers() +
  geom_linerange(aes(ymin = glucose_mgC_per_gdw_min, ymax = glucose_mgC_per_gdw_max),
                 linewidth = 0.4) +
  geom_point(aes(fill = season), shape = 21, size = 1.8) +
  scale_fill_manual(values = season_colors, guide = "none") +
  timeline_x +
  labs(title = "Glucose C added per g dry soil (mean, range)",
       subtitle = "C per well is fixed for each substrate; per gram varies only with soil per well",
       x = NULL, y = "mg C g⁻¹ dry soil") +
  theme_v6
s8 <- (timeline_panel(water, water_dm, "SIR_gdw", "Water-control CO2 by date (5 h)",
                      water_ylab) /
         s8_b) +
  plot_annotation(tag_levels = "A",
                  caption = "Incubation temperature is not recorded in the MicroResp file (README: room temperature; CO2 calculation assumes 20 °C). Run dates are not in the file; QC is by collection date.")
save_v6(s8, "FigS8", 180, 130)
write_parts(bind_rows(
  date_rows(water_dm, "FigS8", "A", "water-control CO2 5 h (ng CO2-C/gdw/h) date mean"),
  result_rows("FigS8", "B", paste0("glucose mg C per g dry soil | ", qc$Collection_date),
              qc$glucose_mgC_per_gdw_mean, qc$glucose_mgC_per_gdw_min,
              qc$glucose_mgC_per_gdw_max, n_dates = 1, n_samples = qc$n_samples,
              notes = "mean; ci columns hold the min and max across samples")
), "FigS8")

versions <- tribble(
  ~version, ~label,
  "primary_5h", "Primary: (substrate − water) / C added",
  "gross_5h", "Gross: substrate / C added",
  "fold_5h", "Fold over water: substrate / water"
)
version_tests <- versions |>
  mutate(dm = map(version, \(v) level_dates(get_z(v))),
         test = map(dm, \(d) season_test(d$mean, d$season)),
         n_samples = map_int(version, \(v) n_distinct(get_z(v)$Sample_cat)))
s9_b <- pmap(list(version_tests$dm, version_tests$label, version_tests$test),
             \(d, l, t) season_panel(d, str_wrap(l, 24), "Date mean of overall level (z)", t))
s9 <- (timeline_panel(water, water_dm, "SIR_gdw", water_title, water_ylab) +
         season_panel(water_dm, "Basal CO2 by season", water_ylab, water_test) +
         plot_layout(widths = c(3, 1))) /
  wrap_plots(s9_b, nrow = 1) +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = list(c("A", "", "B", "", "")))
save_v6(s9, "FigS9", 180, 160)

write_parts(bind_rows(
  date_rows(water_dm, "FigS9", "A", "water-only CO2 5 h (ng CO2-C/gdw/h) date mean"),
  season_rows(water_test, "FigS9", "A", "water-only CO2 5 h dormant - active (ng CO2-C/gdw/h)",
              n_samples = nrow(water)),
  pmap_dfr(version_tests, \(version, label, dm, test, n_samples) {
    season_rows(test, "FigS9", "B", paste0("overall level dormant - active | ", label),
                n_samples = n_samples,
                notes = if_else(version == "primary_5h", "same as Fig2 B", ""))
  })
), "FigS9")

print(lvl_test)
print(class_tests |> select(chem_class, diff, ci_low, ci_high, p_welch, p_holm))
print(pair_ci)
print(pooled_test)
print(mbc_test); print(test24); print(water_test)
print(version_tests |> select(version, test) |> unnest(test) |> select(version, diff, ci_low, ci_high, p_welch))
print(bind_rows(mbc = mbc_test, mbc_no1420 = mbc_test_no1420,
                primary_no1420 = prim_test_no1420, no_phenylethylamine = no_pea_test,
                .id = "test") |> select(test, diff, ci_low, ci_high, p_welch, p_wilcox))
