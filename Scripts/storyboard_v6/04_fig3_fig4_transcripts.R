# =============================================================================
# Figure 3 (Q2), Figure 4 (Q3), and Figure S7 (exploratory one-month offset)
# =============================================================================
#
# Gene sets: Results/plan_v1/KO_gene_sets_APPROVED.csv (v5 sets, D8).
# Expression: share of total expression on a log scale, as in v5 (centred log
# ratio: log(set sum) - mean log expression over all KOs in the sample, after
# replacing zeros by half the smallest nonzero value in that sample).
# Compound gene sets are named after their MicroResp substrate, so they take
# their chem_class from Data/substrate_classes.csv.
# Transcripts are descriptive: differences with 95% CI, no p-values (rule 7).

source(file.path("Scripts", "storyboard_v6", "00_setup.R"))
suppressPackageStartupMessages({
  library(ggh4x)
})

field <- read_csv(file.path(dir_v6_data, "field_100.csv"), show_col_types = FALSE)
samples <- field |>
  select(Sample_cat, Collection_date, Plot) |>
  mutate(season = factor(assign_season(Collection_date), season_levels))

# -----------------------------------------------------------------------------
# Expression (CLR) for every approved set, including the near-zero sets
# -----------------------------------------------------------------------------

gene_sets <- read_csv(file.path(dir_results, "KO_gene_sets_APPROVED.csv"),
                      show_col_types = FALSE)
set_info <- gene_sets |>
  distinct(set, set_type, status) |>
  left_join(substrate_classes |> select(set = substrate, chem_class), by = "set")
stopifnot(all(!is.na(set_info$chem_class[set_info$set_type == "compound"])))

expr <- read_csv(input_files$brite_d, show_col_types = FALSE)
names(expr)[1] <- "brite_D"
expr <- expr |>
  mutate(KO = str_extract(brite_D, "^K\\d{5}")) |>
  filter(!is.na(KO))
sample_cols <- str_subset(names(expr), "_MT$")
mat <- as.matrix(expr[, sample_cols])
rownames(mat) <- expr$KO
half_min <- apply(mat, 2, \(x) min(x[x > 0]) / 2)
mat_rep <- sweep(mat, 2, half_min, \(x, h) ifelse(x == 0, h, x))
mean_log <- colMeans(log(mat_rep))

expr_sets <- gene_sets |>
  distinct(set, KO) |>
  group_by(set) |>
  group_modify(\(d, k) {
    s_raw <- colSums(mat[d$KO, , drop = FALSE])
    s <- colSums(mat_rep[d$KO, , drop = FALSE])
    tibble(sample = names(s), clr = log(s) - mean_log[names(s)],
           zero_sum = s_raw == 0)
  }) |>
  ungroup() |>
  mutate(Sample_cat = str_remove(sample, "_MT$")) |>
  inner_join(samples, by = "Sample_cat") |>
  left_join(set_info, by = "set") |>
  group_by(set) |>
  mutate(z = zscore(clr)) |>
  ungroup()
stopifnot(n_distinct(expr_sets$Sample_cat) == 35)

# Check against plan_v1 values for the sets used there
v1_clr <- read_csv(file.path(dir_results, "cp4_expression_clr_samples.csv"),
                   show_col_types = FALSE)
chk <- inner_join(expr_sets, v1_clr, by = c("set", "Sample_cat"))
stopifnot(max(abs(chk$clr.x - chk$clr.y)) < 1e-10)

transcript_dates <- sort(unique(expr_sets$Collection_date))
compound_sets <- expr_sets |> filter(set_type == "compound")
usable_sets <- compound_sets |> filter(status == "use in Q2")
fig3_classes <- class_levels[class_levels %in% usable_sets$chem_class]
classes_without_sets <- setdiff(class_levels, fig3_classes)

set_dm <- date_mean_ci(compound_sets, "z", set, chem_class, status) |>
  mutate(chem_class = factor(chem_class, class_levels),
         set = factor(set, substrate_levels))

class_samples <- usable_sets |>
  group_by(Sample_cat, Collection_date, Plot, season, chem_class) |>
  summarise(z = mean(z), n_sets = n(), .groups = "drop") |>
  mutate(chem_class = factor(chem_class, fig3_classes))
class_dm <- date_mean_ci(class_samples, "z", chem_class)
class_ci <- class_dm |>
  group_by(chem_class) |>
  group_modify(\(d, k) welch_ci(d$mean, d$season)) |>
  ungroup()
class_n_sets <- usable_sets |> distinct(chem_class, set) |> count(chem_class, name = "n_sets")

# =============================================================================
# Figure 3
# =============================================================================

panel_a <- ggplot(class_dm, aes(Collection_date, mean, colour = chem_class)) +
  timeline_layers(marker_labels = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey50") +
  geom_point(data = class_samples, aes(y = z), colour = "grey65", size = 0.5,
             alpha = 0.7) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 6, linewidth = 0.3) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = class_colors, guide = "none") +
  timeline_x_small +
  facet_wrap(~ chem_class, nrow = 2) +
  labs(title = "Class expression through the year (mean of gene-set z-scores; date means, plots in grey)",
       x = NULL, y = "Expression (z)") +
  theme_v6 +
  theme(axis.text.x = element_text(size = 5.5))

heat_data <- set_dm |>
  filter(chem_class %in% fig3_classes) |>
  mutate(chem_class = factor(chem_class, fig3_classes),
         detected = status == "use in Q2",
         fill_value = if_else(detected, mean, NA_real_),
         run = factor(if_else(season == "Dormant", "Plant-dormant", "Plant-active"),
                      c("Plant-active", "Plant-dormant")))
heat_lim <- max(abs(heat_data$fill_value), na.rm = TRUE)
panel_b <- ggplot(heat_data, aes(factor(format(Collection_date, "%b\n%y"),
                                        format(transcript_dates, "%b\n%y")),
                                 fct_rev(set))) +
  geom_tile(aes(fill = fill_value), colour = "white", linewidth = 0.3) +
  geom_text(data = filter(heat_data, !detected), label = "n.d.", size = 1.8,
            colour = "grey30") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-heat_lim, heat_lim), na.value = "grey80",
                       name = "Expression\n(z, date mean)") +
  scale_y_discrete(labels = \(x) substrate_label(x)) +
  facet_grid2(chem_class ~ run, scales = "free", space = "free", switch = "y",
              labeller = labeller(chem_class = label_wrap_gen(22)),
              strip = strip_themed(
                background_x = elem_list_rect(fill = season_fill),
                text_y = elem_list_text(angle = 0, hjust = 1, size = 5.5))) +
  labs(title = "Gene sets by date (rows ordered as Figure 2C)",
       subtitle = "Grey n.d. = not detected (set sum 0 in 9 or more of 35 samples; D8)",
       x = NULL, y = NULL) +
  theme_v6 +
  theme(strip.placement = "outside", panel.spacing = unit(0.6, "mm"))

# Axis limited to -3..3; CIs past the limit end in an arrow (round 2)
ci_limit <- 3
class_ci_plot <- class_ci |>
  mutate(lo = pmax(ci_low, -ci_limit), hi = pmin(ci_high, ci_limit),
         cut_lo = ci_low < -ci_limit, cut_hi = ci_high > ci_limit)
clipped <- filter(class_ci_plot, cut_lo | cut_hi)
clip_note <- if (nrow(clipped) > 0) {
  paste0("Arrows: CI extends past ±", ci_limit, ".\n",
         paste(sprintf("%s: %.2f (95%% CI %.1f to %.1f)", clipped$chem_class, clipped$diff,
                       clipped$ci_low, clipped$ci_high), collapse = "\n"))
} else ""
arrow_end <- arrow(length = unit(1.2, "mm"), type = "closed")

panel_c <- ggplot(class_ci_plot, aes(diff, fct_rev(chem_class), colour = chem_class)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey50") +
  geom_segment(aes(x = lo, xend = hi, yend = fct_rev(chem_class)), linewidth = 0.45) +
  geom_segment(data = filter(class_ci_plot, !cut_lo),
               aes(x = lo, xend = lo, y = as.numeric(fct_rev(chem_class)) - 0.1,
                   yend = as.numeric(fct_rev(chem_class)) + 0.1), linewidth = 0.45) +
  geom_segment(data = filter(class_ci_plot, !cut_hi),
               aes(x = hi, xend = hi, y = as.numeric(fct_rev(chem_class)) - 0.1,
                   yend = as.numeric(fct_rev(chem_class)) + 0.1), linewidth = 0.45) +
  geom_segment(data = filter(class_ci_plot, cut_lo),
               aes(x = 0.9 * lo, xend = lo, yend = fct_rev(chem_class)),
               linewidth = 0.45, arrow = arrow_end) +
  geom_segment(data = filter(class_ci_plot, cut_hi),
               aes(x = 0.9 * hi, xend = hi, yend = fct_rev(chem_class)),
               linewidth = 0.45, arrow = arrow_end) +
  geom_point(size = 1.8) +
  scale_colour_manual(values = class_colors, guide = "none") +
  scale_y_discrete(labels = label_wrap_gen(22)) +
  scale_x_continuous(limits = c(-ci_limit, ci_limit)) +
  labs(title = "Dormant − active",
       subtitle = "Descriptive: date means,\nWelch 95% CI; 2 dormant,\n5 active dates; no p",
       x = "D − A expression (z)", y = NULL) +
  theme_v6

fig3 <- wrap_elements(full = panel_a) /
  (panel_b + panel_c + plot_layout(widths = c(1.7, 1))) +
  plot_layout(heights = c(1.2, 2)) +
  plot_annotation(tag_levels = "A", caption = paste("C:", clip_note))
save_v6(fig3, "Fig3", 180, 230)

write_parts(bind_rows(
  result_rows("Fig3", "A", paste0("class expression z date mean | ", class_dm$chem_class,
                                  " | ", class_dm$Collection_date),
              class_dm$mean, class_dm$ci_low, class_dm$ci_high, n_dates = 1,
              n_samples = class_dm$n),
  result_rows("Fig3", "B", paste0("gene set expression z date mean | ", heat_data$set,
                                  " | ", heat_data$Collection_date),
              heat_data$mean, heat_data$ci_low, heat_data$ci_high, n_dates = 1,
              n_samples = heat_data$n,
              notes = if_else(heat_data$detected, "",
                              "shown as not detected (insufficient expression, D8)")),
  result_rows("Fig3", "C", paste0("class expression z dormant - active | ",
                                  class_ci$chem_class),
              class_ci$diff, class_ci$ci_low, class_ci$ci_high,
              n_dates = class_ci$n_dormant + class_ci$n_active, n_samples = 35,
              notes = paste0("descriptive; Welch CI on date means; ",
                             class_n_sets$n_sets[match(class_ci$chem_class,
                                                       class_n_sets$chem_class)],
                             " gene sets"))
), "Fig3")

# =============================================================================
# Figure 4
# =============================================================================
# All series are z-scores formed at the sample level (MicroResp per substrate
# across 77 samples, level included; enzymes across the 78 samples;
# transcripts per gene set across 35 samples), then averaged by date.

z5 <- read_csv(file.path(dir_v6_data, "microresp_z.csv"), show_col_types = FALSE) |>
  filter(version == "primary_5h") |>
  mutate(season = factor(season, season_levels))

eea <- read_csv(input_files$eea, na = c("", "NA", "na", "#N/A"), show_col_types = FALSE) |>
  select(Sample_cat, BG_nmol_gdw_h:Ala_nmol_gdw_h) |>
  inner_join(samples, by = "Sample_cat") |>
  pivot_longer(BG_nmol_gdw_h:Ala_nmol_gdw_h, names_to = "enzyme", values_to = "activity") |>
  mutate(enzyme = str_remove(enzyme, "_nmol_gdw_h")) |>
  filter(!is.na(activity)) |>
  group_by(enzyme) |>
  mutate(z = zscore(activity)) |>
  ungroup()

# Enzyme assigned to a class where it acts on that class's compounds
enzyme_class <- tribble(
  ~enzyme, ~chem_class,
  "CB", "Plant structural polymers",
  "BX", "Plant structural polymers",
  "BG", "Simple sugars and oligosaccharides",
  "NAG", "Chitinous structural polymers",
  "Leu", "Amino acids and related compounds",
  "Ala", "Amino acids and related compounds"
)

measure_colors <- c(MicroResp = "#B8562B", Enzyme = "#6B4C9A", Transcripts = "#2B6CB8")

mr_class <- z5 |>
  filter(chem_class %in% fig3_classes) |>
  group_by(Sample_cat, Collection_date, season, chem_class) |>
  summarise(z = mean(z), .groups = "drop") |>
  date_mean_ci("z", chem_class) |>
  mutate(measure = "MicroResp")
enz_class <- eea |>
  inner_join(enzyme_class, by = "enzyme") |>
  group_by(Sample_cat, Collection_date, season, chem_class) |>
  summarise(z = mean(z), .groups = "drop") |>
  date_mean_ci("z", chem_class) |>
  mutate(measure = "Enzyme")
tx_class <- class_dm |> mutate(chem_class = as.character(chem_class), measure = "Transcripts")

fig4a_data <- bind_rows(mr_class, enz_class, tx_class) |>
  filter(Collection_date %in% transcript_dates) |>
  mutate(chem_class = factor(chem_class, fig3_classes),
         measure = factor(measure, names(measure_colors)))
enz_labels <- enzyme_class |>
  group_by(chem_class) |>
  summarise(enz = paste(enzyme, collapse = ", ")) |>
  right_join(tibble(chem_class = fig3_classes), by = "chem_class") |>
  mutate(strip = paste0(chem_class, "\n(enzyme: ", replace_na(enz, "none"), ")"))
strip_lookup <- setNames(enz_labels$strip, enz_labels$chem_class)

dodge <- position_dodge(width = 10)
fig4a <- ggplot(fig4a_data, aes(Collection_date, mean, colour = measure)) +
  timeline_layers(marker_labels = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey50") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, linewidth = 0.3,
                position = dodge) +
  geom_line(linewidth = 0.5, position = dodge) +
  geom_point(size = 1, position = dodge) +
  scale_colour_manual(values = measure_colors, name = NULL) +
  timeline_x_small +
  facet_wrap(~ chem_class, nrow = 2,
             labeller = as_labeller(\(x) strip_lookup[x])) +
  labs(title = "Potential, enzymes, and expression on the 7 transcript dates (z-scores, date means ± 95% CI)",
       x = NULL, y = "z-score") +
  theme_v6 +
  theme(legend.position = "bottom", axis.text.x = element_text(size = 5.5))

# Pathway case studies: one row per step, in pathway order
pathway_series <- function(steps) {
  map_dfr(seq_len(nrow(steps)), \(i) {
    s <- steps[i, ]
    d <- switch(s$measure,
                MicroResp = filter(z5, Substrate == s$id),
                Enzyme = filter(eea, enzyme == s$id),
                Transcripts = filter(expr_sets, set == s$id))
    date_mean_ci(d, "z") |> mutate(step = s$step, measure = s$measure)
  }) |>
    mutate(step = factor(step, steps$step),
           measure = factor(measure, names(measure_colors)))
}

# Lines between consecutive dates; dashed where dates are more than 45 days
# apart (transcripts: Oct -> Jan and Mar -> May), as there is no data between
pathway_panel <- function(d, title) {
  segs <- d |>
    arrange(step, Collection_date) |>
    group_by(step, measure) |>
    mutate(x_end = lead(Collection_date), y_end = lead(mean)) |>
    ungroup() |>
    filter(!is.na(x_end)) |>
    mutate(gap = if_else(as.numeric(x_end - Collection_date) > 45, "gap", "adjacent"))
  ggplot(d, aes(Collection_date, mean, colour = measure)) +
    timeline_layers(marker_labels = FALSE) +
    geom_vline(xintercept = transcript_dates, linewidth = 0.2, colour = "grey80") +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey50") +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0, linewidth = 0.3) +
    geom_segment(data = segs, aes(xend = x_end, yend = y_end, linetype = gap),
                 linewidth = 0.5) +
    geom_point(size = 1) +
    scale_colour_manual(values = measure_colors, name = NULL, drop = FALSE) +
    scale_linetype_manual(values = c(adjacent = "solid", gap = "22"), guide = "none") +
    timeline_x_small +
    facet_grid(step ~ ., labeller = label_wrap_gen(22)) +
    labs(title = title, x = NULL, y = "z-score (date mean ± 95% CI)") +
    theme_v6 +
    theme(strip.text.y = element_text(angle = 0, size = 6), legend.position = "none",
          axis.text.x = element_text(size = 5.5))
}

chitin_steps <- tribble(
  ~step, ~measure, ~id,
  "Chitinase genes", "Transcripts", "Chitin",
  "Chitin response", "MicroResp", "Chitin",
  "NAG enzyme activity", "Enzyme", "NAG",
  "N-acetylglucosamine response", "MicroResp", "N.Acetyl.D.glucosamine",
  "NAG uptake genes", "Transcripts", "N.Acetyl.D.glucosamine"
)
cellulose_steps <- tribble(
  ~step, ~measure, ~id,
  "Cellulase / cellobiohydrolase genes", "Transcripts", "Cellulose",
  "Cellulose response", "MicroResp", "Cellulose",
  "CB enzyme activity", "Enzyme", "CB",
  "Cellobiose response", "MicroResp", "D.Cellobiose",
  "BG enzyme activity", "Enzyme", "BG",
  "β-glucosidase genes", "Transcripts", "BG",
  "Glucose response", "MicroResp", "D.Glucose"
)
chitin_data <- pathway_series(chitin_steps)
cellulose_data <- pathway_series(cellulose_steps)

fig4 <- fig4a /
  (pathway_panel(chitin_data, "Chitin pathway") +
     pathway_panel(cellulose_data, "Cellulose pathway")) +
  plot_layout(heights = c(1, 1.5)) +
  plot_annotation(tag_levels = "A",
                  caption = "MicroResp: water-subtracted response per mg C added (per-biomass version in Figure S4).\nTranscripts on 7 dates (thin grey lines); MicroResp and enzymes on 13.\nIn B and C, dashed lines join dates more than 45 days apart (no data between).")
save_v6(fig4, "Fig4", 180, 240)

series_rows <- function(d, panel, prefix) {
  result_rows("Fig4", panel,
              paste0(prefix, d$measure, " | ", d$group, " | ", d$Collection_date),
              d$mean, d$ci_low, d$ci_high, n_dates = 1, n_samples = d$n,
              notes = "z-score date mean")
}
write_parts(bind_rows(
  series_rows(fig4a_data |> mutate(group = chem_class), "A", "class z | "),
  series_rows(chitin_data |> mutate(group = step), "B", "chitin pathway | "),
  series_rows(cellulose_data |> mutate(group = step), "C", "cellulose pathway | ")
), "Fig4")

# =============================================================================
# Figure S7. Exploratory one-month offset (no significance claims)
# =============================================================================
# For each class: Spearman correlation across the transcript dates between
# class expression at date t and MicroResp class z (level included) at t
# ("same date") or at the next sampling date ("next date"). Leave-one-date-out
# ranges drop each transcript date in turn.

mr_all <- z5 |>
  filter(chem_class %in% fig3_classes) |>
  group_by(Sample_cat, Collection_date, chem_class) |>
  summarise(z = mean(z), .groups = "drop") |>
  group_by(Collection_date, chem_class) |>
  summarise(mr = mean(z), .groups = "drop")
next_date <- setNames(c(study_dates[-1], NA), as.character(study_dates))

offset_pairs <- tx_class |>
  select(chem_class, Collection_date, tx = mean) |>
  mutate(next_dt = as.Date(next_date[as.character(Collection_date)])) |>
  left_join(mr_all |> rename(mr_same = mr), by = c("chem_class", "Collection_date")) |>
  left_join(mr_all |> rename(next_dt = Collection_date, mr_next = mr),
            by = c("chem_class", "next_dt"))
stopifnot(!anyNA(offset_pairs))

rho <- function(x, y) cor(x, y, method = "spearman")
offset <- offset_pairs |>
  pivot_longer(c(mr_same, mr_next), names_to = "offset", values_to = "mr") |>
  mutate(offset = recode(offset, mr_same = "Same date", mr_next = "Next date")) |>
  group_by(chem_class, offset) |>
  summarise(rho = rho(tx, mr),
            lodo = list(map_dbl(seq_along(tx), \(i) rho(tx[-i], mr[-i]))),
            n_dates = n(), .groups = "drop") |>
  mutate(lodo_min = map_dbl(lodo, min), lodo_max = map_dbl(lodo, max),
         chem_class = factor(chem_class, fig3_classes),
         offset = factor(offset, c("Same date", "Next date")))

s7 <- ggplot(offset, aes(rho, fct_rev(chem_class), colour = offset)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey50") +
  geom_errorbar(aes(xmin = lodo_min, xmax = lodo_max), width = 0.2, orientation = "y",
                linewidth = 0.4, position = position_dodge(width = 0.5)) +
  geom_point(size = 1.8, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = c("Same date" = "grey40", "Next date" = "#B8562B"),
                      name = "MicroResp date") +
  scale_y_discrete(labels = label_wrap_gen(24)) +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(title = "EXPLORATORY: one-month offset",
       caption = "Transcripts at date t vs MicroResp class z (level included) at t or at the next sampling date.",
       subtitle = "Spearman ρ across the 7 transcript dates; bars = leave-one-date-out range.\nNo significance claims.",
       x = "Spearman ρ", y = NULL) +
  theme_v6 +
  theme(legend.position = "bottom")
save_v6(s7, "FigS7", 120, 100)
write_parts(result_rows("FigS7", "",
                        paste0("Spearman rho transcripts vs MicroResp | ", offset$chem_class,
                               " | ", offset$offset),
                        offset$rho, offset$lodo_min, offset$lodo_max,
                        n_dates = offset$n_dates,
                        notes = "exploratory; ci columns hold the leave-one-date-out min and max; no p-value"),
            "FigS7")

writeLines(c(
  paste("Classes without gene sets (left out of Figure 3):",
        paste(classes_without_sets, collapse = "; ")),
  paste("Substrates in Figure 3 classes without a usable gene set:",
        paste(setdiff(substrate_classes$substrate[substrate_classes$chem_class %in% fig3_classes],
                      usable_sets$set), collapse = "; "))
), file.path(dir_v6_data, "fig3_coverage.txt"))

print(class_ci)
print(offset |> select(-lodo))
cat(readLines(file.path(dir_v6_data, "fig3_coverage.txt")), sep = "\n")
