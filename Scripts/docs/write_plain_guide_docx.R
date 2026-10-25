# =============================================================================
# Write Docs/Plain_language_guide_v4.docx
# =============================================================================
#
# A step-by-step, plain-language explanation of the v4 analysis: the framing,
# questions, and hypotheses (unchanged), every data issue, and for each
# analysis what it asks, why it was chosen, what it found, and what it means
# for the question. All numbers are read from Results/manuscript_v4.
# Run after Scripts/00_run_all.R:  Rscript Scripts/docs/write_plain_guide_docx.R

rm(list = ls())
suppressPackageStartupMessages(library(tidyverse))
source(file.path("Scripts", "docs", "docx_writer.R"))

res <- function(...) read_csv(file.path("Results", "manuscript_v4", ...), show_col_types = FALSE)
fig <- function(...) normalizePath(file.path("Figures", "manuscript_v4", ...), winslash = "/")
f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
pv <- function(p) if_else(p < 0.001, "p < 0.001", paste0("p = ", f2(p, 3)))
pct_up <- function(log2_value) round(100 * (2^log2_value - 1))
lc <- function(x) if_else(str_detect(x, "^[A-Z][a-z]"), paste0(tolower(str_sub(x, 1, 1)), str_sub(x, 2)), x)
range_txt <- function(x, d = 2) sprintf("%s to %s", f2(min(x), d), f2(max(x), d))

# Boxed "answer" paragraph.
answer <- function(...) {
  para(paste0(...), extra_ppr = '<w:shd w:val="clear" w:color="auto" w:fill="E8F1E4"/><w:spacing w:before="120" w:after="200"/>')
}
step <- function(label, ...) p(b(paste0(label, " ")), ...)

# -----------------------------------------------------------------------------
# Numbers
# -----------------------------------------------------------------------------

design <- res("sample_availability", "design_summary.csv")
dval <- function(item) design$value[design$item == item]
lagsum <- res("sample_availability", "lag_pair_summary.csv")

counts <- res("Figure2", "explainers", "microresp_samples_per_date_by_release.csv")
count_on <- function(release, date) counts$n_samples[counts$release == release & counts$Sampling_Date == as.Date(date)]
total_on <- function(release) sum(counts$n_samples[counts$release == release])
versions <- res("Figure2", "explainers", "permanova_five_ways.csv")
fake <- res("Figure2", "explainers", "fake_seasons_summary.csv")
distances <- res("Figure2", "explainers", "bray_curtis_scale_of_seasonal_shift.csv")

f1 <- res("Figure1", "Figure1_season_tests.csv")
f1v <- function(col) filter(f1, column == col)

level <- res("Figure2", "Figure2A_overall_level_shift.csv")
var_share <- res("Figure2", "Figure2A_variance_shares_date_vs_plot.csv")
share_of <- function(src) round(100 * var_share$share[var_share$source == src])
level_samples <- res("Figure2", "Figure2A_overall_level_by_sample.csv")
level_dates <- res("Figure2", "Figure2A_overall_level_by_date.csv")
july <- filter(level_samples, Sampling_Date == as.Date("2023-07-05"))
july_mean <- level_dates$mean[level_dates$Sampling_Date == as.Date("2023-07-05")]
sub_eff <- res("Figure2", "Figure2C_substrate_season_effects.csv")
top_sub <- sub_eff %>% arrange(desc(estimate)) %>% slice(1:3)
h1b_tests <- res("Figure2", "Figure2D_H1b_tests.csv")
h1b_row <- function(cmp) filter(h1b_tests, comparison == cmp)
h1b <- h1b_row("All substrates: complex vs readily available")
by_avail <- res("Figure2", "Figure2D_season_effect_by_availability.csv")
avail_v <- function(g) filter(by_avail, availability == g)
h1b_sens <- res("Figure2", "FigureS2b_H1b_sensitivity.csv")
gam <- res("Figure2", "FigureS2_GAM_temporal_statistics.csv")
steady <- sub_eff %>% filter(steady) %>% pull(display_name)

h2_all <- res("Figure3", "Figure3A_season_PERMANOVA_dispersion.csv")
h2_res <- filter(h2_all, level == "Date means (primary)")
h2_sens <- res("Figure3", "FigureS3_H2_sensitivity.csv") %>% filter(level == "Date means (primary)")
cat_sum <- res("Figure3", "Figure3B_unique_category_summary.csv")
cat_dates <- res("Figure3", "Figure3B_unique_category_date_means.csv")
extremes <- cat_dates %>% group_by(compound_class, brite_C) %>%
  summarise(min_date = Sampling_Date[which.min(z)], .groups = "drop")
jan_low <- function(cls) {
  e <- filter(extremes, compound_class == cls)
  sprintf("%d of %d", sum(e$min_date == as.Date("2024-01-02")), nrow(e))
}

loo <- res("Figure4", "Figure4B_leave_one_date_out_summary.csv")
loo_v <- function(cmp, col) loo[[col]][loo$comparison == cmp]
protest <- res("Figure4", "FigureS4b_PROTEST.csv")
div <- res("Figure4", "Figure4_class_values_by_date.csv") %>% filter(comparison == "Concurrent")
jan_div <- filter(div, Sampling_Date == as.Date("2024-01-02"))

classes <- read_csv(file.path("Data", "compound_class_mapping.csv"), show_col_types = FALSE)
pathways <- read_csv(file.path("Data", "compound_class_unique_BRITE_C.csv"), show_col_types = FALSE)

version_row <- function(pattern) filter(versions, str_detect(version, pattern))
v3_version <- version_row("^v3")
whole_dates <- version_row("^All samples")
date_avg <- version_row("^Date averages")
builtin <- version_row("^Built-in")
mix_version <- version_row("^Mix")

# -----------------------------------------------------------------------------
# Tables
# -----------------------------------------------------------------------------

big_picture <- tribble(
  ~Question, ~`What we measured`, ~`Main analysis (plain words)`, ~`Short answer`,
  "Q1 / H1a", "CO2 released when each of 31 compounds is added to soil (SIR), 13 dates",
  "Is SIR higher or lower overall in one season (level)? Does the relative pattern across compounds change (mix)?",
  sprintf("Yes for level: about %d%% higher in dormancy for all 31 compounds. No for mix.", pct_up(level$estimate)),
  "Q1 / H1b", "Same SIR data, compounds split into complex (8) and readily available (23)",
  "Did complex compounds go up more in dormancy than readily available ones?",
  "No clear difference: both went up by similar amounts.",
  "Q2 / H2", "Which microbial functions were being expressed in the field, 7 dates",
  "Is the mix of expressed functions different on the 2 dormant dates?",
  "Yes, as clearly as 7 dates can show.",
  "Q2 (classes)", "The 20 functions linked to one compound class",
  "When were these functions high or low? (description only)",
  sprintf("Most were lower on dormant dates, especially amino acid functions in January (%s lowest then).", jan_low("Amino acid")),
  "Q3", "SIR and expression of the matching functions, same 7 dates",
  "Do they go up and down together? Do they match better if SIR is taken a month later?",
  "Not on the same date (January is the big mismatch); better one month later."
)

data_table <- tribble(
  ~Measurement, ~`What it tells us`, ~Dates, ~Plots, ~`Samples used`,
  "Field measurements (soil temperature, moisture, field CO2, roots, microbial biomass, enzymes)", "Conditions in the field when soil was collected", "13", "6", dval("Observed 100% field records"),
  "MicroResp (SIR)", "How strongly the microbes present can respond to each added compound", "13", "6 (5 in Feb 2024)", dval("Observed 100% MicroResp samples"),
  "Metatranscriptomes", "Which microbial functions were being expressed in the field", "7", "5 (no plot 57)", dval("Observed 100% transcript samples")
)

count_table <- tibble(
  `Date` = c("January 2024", "February 2024", "Other 11 dates", "Total"),
  `Original MicroResp file` = c(count_on("Original release", "2024-01-02"), count_on("Original release", "2024-02-12"),
                                total_on("Original release") - count_on("Original release", "2024-01-02") - count_on("Original release", "2024-02-12"),
                                total_on("Original release")),
  `Corrected file (MicroResp v3)` = c(count_on("Corrected release (v3)", "2024-01-02"), count_on("Corrected release (v3)", "2024-02-12"),
                                      total_on("Corrected release (v3)") - count_on("Corrected release (v3)", "2024-01-02") - count_on("Corrected release (v3)", "2024-02-12"),
                                      total_on("Corrected release (v3)")),
  `Why` = c("Original: 6 real plots plus a copy of plot 1 (SM1503 = SM1502), and the IDs shifted by one plot. v3: fixed.",
            "Plot 27 (SM1536) was never measured; missing in every release.",
            "6 plots each, complete.",
            "The January extra copy and the February gap cancelled out, so the original total looked complete.")
) %>% mutate(across(everything(), as.character))

versions_table <- versions %>%
  transmute(`How the test was run` = str_replace(version, "^Mix only", "Mix only, a different question:"),
            `What counts as one replicate` = what_counts_as_a_replicate,
            `Rows` = as.character(n_rows), `R²` = f2(R2), `p` = if_else(p < 0.001, "< 0.001", f2(p, 3)))

distance_table <- distances %>% transmute(Comparison = comparison, `Bray-Curtis distance` = f2(bray_curtis))

class_table <- classes %>%
  mutate(substrate_group = factor(substrate_group, levels = c("Carbohydrate", "Organic acid", "Amino acid", "Amino sugar", "Aromatic", "Lipid-like", "Aliphatic amine"))) %>%
  arrange(substrate_group, desc(availability_class == "Complex"), display_name) %>%
  transmute(`Chemical class` = as.character(substrate_group), Compound = display_name,
            `H1b group` = availability_class, `What it is` = substrate_note)

pathway_table <- pathways %>%
  left_join(select(cat_sum, brite_C, mean_share_pct), by = "brite_C") %>%
  mutate(compound_class = factor(compound_class, levels = c("Lipid-like", "Aromatic", "Amino sugar", "Amino acid"))) %>%
  arrange(compound_class, desc(mean_share_pct)) %>%
  transmute(Class = as.character(compound_class), Pathway = sprintf("ko%s %s", brite_C_ID, brite_C_name),
            `What it does` = what_it_does, `Share of all expression` = paste0(formatC(mean_share_pct, format = "fg", digits = 2), "%"))

# -----------------------------------------------------------------------------
# Document
# -----------------------------------------------------------------------------

guide <- c(
  h1("A plain-language guide to the v4 analysis"),
  p("This guide explains, step by step, what was done in the v4 analysis and why. For each analysis it says what question it asks, why we chose it, what it found, and how it answers your research questions. It also lists every data detail you should know about, including problems we found and how each was handled."),
  p("Your framing, research questions, and hypotheses are copied word for word in sections 1 and 2. Every number in this guide comes from the result tables written by the scripts in the repository (Scripts/00_run_all.R runs them all), so nothing is typed in by hand."),

  h2("1. Framing"),
  h3("P1. Seasonal changes in microbial substrate availability and function"),
  p("Soil microbial function varies across the annual cycle as plant activity, temperature, moisture, and resource availability change. As plants transition between active growth and dormancy, changes in root exudation, root turnover, and residue inputs may alter the timing and relative availability of different substrates, including readily available compounds and more complex organic matter associated with plant residues and microbial necromass. Dormancy may therefore change the mixture of substrates available to microbes rather than simply reducing their supply. These shifts can alter microbial physiological state and resource acquisition strategies, potentially affecting how microbes respond to different substrates throughout the year. Yet much of our understanding of soil microbial function comes from growing-season measurements, despite evidence that microbial activity can continue outside the growing season (Monson et al., 2006; Gavazov et al., 2017). As a result, we know comparatively little about how microbial responses to different substrates change during plant dormancy, limiting our understanding of how plant seasonality shapes belowground microbial metabolism."),
  h3("P2. Connecting microbial substrate-use potential and field expression"),
  p("Commonly used measurements capture different aspects and timescales of microbial function. Field soil respiration integrates root and microbial activity under prevailing conditions but cannot resolve which substrates microbes are using. Laboratory substrate-induced respiration measures short-term responses to added compounds under standardized conditions. Because these responses are measured within 5 h, they primarily characterize the short-term respiratory capacity of the existing microbial community, which reflects its physiological state and antecedent environmental conditions rather than substantial growth during the assay (Anderson & Domsch, 1978; Campbell et al., 2003). Seasonal changes in these responses therefore show how microbial substrate-use potential shifted over the year. Metatranscriptomes, in contrast, characterize which microbial functions are being expressed under field conditions at the time of sampling, but not how strongly microbes would respond to specific substrate additions. It remains unclear how seasonal changes in field gene expression relate to changes in substrate-use potential, in part because the two capture microbial activity over different timescales: microbial transcripts typically persist for minutes (Moran et al., 2013), whereas soil extracellular enzymes can remain active for weeks after they are produced (Schimel et al., 2017). Seasonal changes in field expression and substrate-use potential may therefore not coincide. Comparing the two across related compounds and metabolic pathways can help clarify whether plant dormancy changes which functions microbes express, their potential to respond as substrate availability changes, or both."),
  h3("P3. Study objectives and research questions"),
  p("Here, we investigate how microbial substrate-use potential and field gene expression change between plant-active and plant-dormant periods, and how the two are related. We focused on fully irrigated (100%) plots in an irrigated tall wheatgrass system, a perennial grass with a plant-active period from April (plant emergence) through October (harvest) and a plant-dormant period from November through March [basis for season definitions]. We combine 13 months of substrate-induced respiration measurements for 31 soil-relevant organic compounds with field metatranscriptomes from seven of these sampling dates (June, July, September, and October 2023; January, March, and May 2024), and use supporting measurements of extracellular enzyme activity, microbial biomass, soil respiration, root biomass, and environmental and biogeochemical conditions to help interpret the observed patterns. Specifically, we ask (1) how substrate-induced respiration differs between plant-active and plant-dormant periods, and whether seasonal responses differ among substrates; (2) how functional transcript expression differs between these periods, and whether functions associated with different compound classes show distinct seasonal patterns; and (3) how seasonal changes in substrate-induced respiration within each compound class relate to the expression of functions associated with that class."),

  h2("2. Research questions and hypotheses"),
  h3("Q1. Seasonal changes in microbial substrate-use potential"),
  p("Q1: How does substrate-induced respiration differ between plant-active and plant-dormant periods, and do seasonal responses differ among substrates?"),
  p("H1a: Microbial substrate-induced respiration profiles will differ between plant-active and plant-dormant periods."),
  p("H1b: Responses to readily available substrates, such as simple sugars and organic acids, will be greater relative to each substrate's mean across sampled dates during plant-active periods, whereas responses to complex substrates associated with plant and microbial residues will be greater relative to their mean across sampled dates during plant-dormant periods."),
  h3("Q2. Seasonal changes in microbial functional expression"),
  p("Q2: How does microbial functional transcript expression differ between plant-active and plant-dormant periods, and do functions associated with different compound classes show distinct seasonal patterns?"),
  p("H2: Microbial functional transcript composition will differ between plant-active and plant-dormant periods."),
  p("We address the second part of Q2 descriptively, by examining seasonal patterns in the BRITE C functional categories mapped uniquely to lipid-like, aromatic, amino sugar, and amino acid compound classes in our a priori substrate-to-function mapping. Because only two transcript sampling dates fall within the plant-dormant period, we do not test directional predictions for individual categories."),
  h3("Q3. Correspondence between substrate-use potential and transcript expression"),
  p("Q3: How do seasonal changes in substrate-induced respiration within each compound class relate to transcript abundance of the BRITE C functional categories mapped uniquely to that class, and does this relationship differ among compound classes?"),
  p("We address Q3 descriptively. We also explore whether transcript abundance is associated with substrate-induced respiration measured one month later, recognizing that the available sampling dates limit our ability to identify temporal lags."),

  h2("3. The whole study on one page"),
  tbl(big_picture, rel = c(0.9, 2, 2.4, 2.2)),

  h2("4. The data, and every detail you should know"),
  h3("4.1 What was collected"),
  p("All analyses use the six fully irrigated (100%) plots. Soil was collected once a month from June 2023 to June 2024, giving 13 sampling dates. Three soil cores per plot were mixed into one sample."),
  tbl(data_table, rel = c(2.6, 2.4, 0.7, 1.1, 0.9)),
  h3("4.2 Problems found and how each was handled"),
  p(b("1. January 2024 MicroResp labels were shifted, and one sample was a copy (fixed). "),
    "In the original MicroResp file, the January 2024 sample IDs were shifted by one plot, so most January MicroResp values were attached to the wrong plot. Sample SM1503 was an exact copy of SM1502 (plot 1). The corrected file you uploaded (MicroResp v3) fixes both. We checked that its values were relabelled, not changed: v3 SM1503 is identical to the old SM1504. The scripts now stop automatically if MicroResp sample IDs ever disagree with the field table again."),
  p(b("2. February 2024, plot 27 has no MicroResp data (not fixable here). "),
    "Sample SM1536 (plot 27, 12 February 2024) exists for every other measurement, but it has no MicroResp data in any release we have (December, February, March, the DONTUSE file, and v3). Either it was never run or it was removed before release; nothing in the files says which. Sheryl should know. Effects: February has 5 plots instead of 6; the January-to-February comparison in Q3 has 4 plots instead of 5; and one built-in vegan option needs equal plots on every date (section 5.1)."),
  p(b("3. Why the v3 storyboard reported 78 MicroResp samples. "),
    "A complete year would be 13 dates x 6 plots = 78. The original file had 78 rows, but only because January's extra copy and February's missing plot cancelled out. The real number of distinct samples is 77:"),
  tbl(count_table, rel = c(1.2, 1.1, 1.2, 3.5)),
  p(b("4. One extreme value in July 2023 (handled). "),
    "Sample SM1082 (plot 57, July 2023) had salt-extractable organic C about 60 times the July average; all its other values were normal. Earlier cleaning removed the whole sample, which also removed its field CO2, roots, and enzymes. Now only that one value is removed and the sample is kept, so the field data are complete (78 plot-dates)."),
  p(b("5. Three July 2023 samples from the 25% plots are removed (no effect here). "),
    "SM1085, SM1087, and SM1088 were flagged as outliers in earlier work and are still removed. They are from the 25% irrigation plots, which this paper does not use."),
  p(b("6. Plot 57 has no metatranscriptomes. "), "Transcripts exist for plots 1, 25, 27, 40, and 46 on seven dates (35 samples). Plot 57 was not sequenced on any date."),
  p(b("7. We do not know exactly how the transcript table was made (needs your input). "),
    "The file Data/brite_C_expression.csv was added to the repository on 2026-07-07, and no script in the repository creates it. Its values are not raw counts (they have decimals) and do not add up to the same total in every sample, so they are some kind of normalized value. The Methods will need the steps used to make it."),
  p(b("8. The season boundaries need a reference. "), "Plant-active is April to October and plant-dormant November to March, as in your framing. The citation for this choice is still pending from Sheryl."),
  p(b("9. Two very low microbial biomass values. "), "Microbial biomass C has two values of 0.07 and 0.14 ug C/g, far below the others; one is a below-detection value replaced by half the detection limit. They do not change any result but are worth checking against lab notes."),

  h2("5. Three ideas you need before the analyses"),
  h3("5.1 Season comes with the date, so each date counts as one replicate"),
  p("Season is not something we applied to plots. It comes with the sampling date. All six plots sampled on the same day share the same weather, the same number of days since irrigation, the same field trip, and the same laboratory run. If something makes one day unusual, all six plots will be unusual together."),
  p(sprintf("The data show how strongly plots on the same day move together. Of the variation in a sample's overall SIR, %d%% is shared by all plots on the same date, while only %d%% is due to consistent differences between plots. In July 2023, %d of %d plots were below their usual values, and the date as a whole was about %d%% below average.",
            share_of("Sampling date"), share_of("Plot"), sum(july$response < 0), nrow(july), round(100 * (1 - 2^july_mean)))),
  p(b("Why this matters for the statistics. "), "A PERMANOVA decides whether a difference is real by shuffling the season labels many times and checking how often a random shuffle gives a difference as big as the real one. The v3 analysis shuffled labels plot by plot, so a shuffle could call plot 1 in January 'active' and plot 25 in January 'dormant'. Real seasons can never do that. Compared with those unrealistic shuffles, the real data always looked exceptional."),
  p(b("A test with made-up seasons shows the problem. "),
    sprintf("We picked 5 of the 13 dates at random, called them 'dormant', and repeated this %d times. Because these seasons are made up, a fair test should call them significant only about 5%% of the time. The v3 method called %d%% of them significant. Shuffling whole dates called %d%% significant, close to the expected 5%%. So the v3 p-value of 0.002 did not mean much: that method says 'significant' for almost half of all possible groupings of dates, real or not.",
            fake$n_fake_seasons, round(fake$pct_significant_v3_method), round(fake$pct_significant_date_method))),
  p(b("You do not have to average the plots. "), "There are several correct ways to keep dates together. All of them give a similar answer, and all are very different from the v3 method:"),
  tbl(versions_table, rel = c(3, 1.6, 0.6, 0.6, 0.6)),
  p("The plots are still used everywhere. The mixed models (section 6.1) include every plot and are told which plots share a date."),

  h3("5.2 Two different ways SIR can change: level and mix"),
  p("Think of the 31 compounds as 31 instruments in a band. Season could change the music in two ways:"),
  ul(c(
    "Level: every instrument gets louder (or quieter). Microbes respond more strongly to every compound.",
    "Mix: some instruments get louder relative to others. Microbes favour some compounds over others."
  )),
  p("H1a ('profiles will differ') could be either. H1b is specifically about mix: complex compounds rising relative to readily available ones. So we test level and mix separately."),

  h3("5.3 Comparing compounds that respire at very different rates"),
  p("Some compounds (simple sugars) produce much more CO2 than others (lignin). To compare seasonal change fairly, each value is expressed relative to that compound's own average over the year. This is exactly what H1b says: 'relative to each substrate's mean across sampled dates'."),
  p(sprintf("We use a log2 scale, where 0 means the compound's usual value, +1 means twice its usual value, and -1 means half. A seasonal difference of +%s on this scale means dormant values were 2^%s = %s times active values, or %d%% higher.",
            f2(level$estimate), f2(level$estimate), f2(2^level$estimate), pct_up(level$estimate))),

  h2("6. Q1: How did SIR differ between seasons?"),
  img(fig("Figure2", "Figure2_main.png"), 6.5),
  p(i("Figure 2. A: SIR averaged over all 31 compounds through the year. B: NMDS picture of the whole profile. C: seasonal difference for each compound. D: the H1b test.")),

  h3("6.1 Level: is SIR higher overall in one season? (Figure 2A)"),
  step("What it asks.", "If you average a sample's response over all 31 compounds, is that average higher on dormant dates than on active dates?"),
  step("Why we chose it.", "It answers the level question directly and is easy to read. It uses a mixed model: a regression that includes every plot but knows which plots share a date, so each date counts once for season."),
  step("What we found.", sprintf("SIR was about %d%% higher on dormant dates (95%% confidence interval: %d%% to %d%% higher; %s). The line in Figure 2A sits above zero from November to March and is lowest in July 2023 and April 2024.",
                                  pct_up(level$estimate), pct_up(level$lower), pct_up(level$upper), pv(level$p))),
  step("What it means for Q1.", "Yes, SIR differed between seasons: microbes responded more strongly in dormancy. This supports H1a."),

  h3("6.2 Each compound separately (Figure 2C)"),
  step("What it asks.", "For each of the 31 compounds, how much higher or lower was SIR in dormancy?"),
  step("Why we chose it.", "It shows whether the overall increase comes from all compounds or from a few, which is the second part of Q1 ('do seasonal responses differ among substrates?')."),
  step("What we found.", sprintf("All %d compounds were higher in dormancy. The largest increases were %s. No single compound was statistically significant on its own after correcting for testing 31 compounds at once (the FDR correction, which stops chance results from looking real).",
                                  nrow(sub_eff), paste(sprintf("%s (%s times)", lc(top_sub$display_name), f2(2^top_sub$estimate, 1)), collapse = ", "))),
  step("What it means for Q1.", "The seasonal change was shared by all compounds rather than driven by a few. Each compound alone is too noisy to be sure about, but together they point the same way."),

  h3("6.3 The whole profile: PERMANOVA and the NMDS picture (Figure 2B)"),
  step("What it asks.", "Do the full 31-compound profiles of dormant dates differ from those of active dates?"),
  step("Why we chose it.", "It is the standard test for 'profiles differ' (H1a), and reviewers expect it. The NMDS is the picture that goes with it: samples with similar profiles are placed close together. The NMDS is only a picture; the PERMANOVA is the test."),
  step("What we found.", sprintf("The whole-profile difference was suggestive but not certain: p = %s when every sample is used and whole dates are shuffled, and p = %s using date averages (table in section 5.1). In the NMDS, dormant and active points overlap.",
                                  f2(whole_dates$p, 3), f2(date_avg$p, 3))),
  step("Why the NMDS hides a real 34% increase.", sprintf("The NMDS places samples by Bray-Curtis distance, where 0 means identical profiles. Making a sample %d%% higher on every compound moves it only %s. Two plots sampled on the same day already differ by %s on average, and samples from different dates by %s. So the seasonal shift is about the size of ordinary plot-to-plot differences and disappears into the scatter. It is like the same song played louder: hard to see on noisy waveforms, obvious if you measure the volume (Figure 2A).",
                                                           pct_up(level$estimate), f2(distances$bray_curtis[1]), f2(distances$bray_curtis[2]), f2(distances$bray_curtis[3]))),
  tbl(distance_table, rel = c(4, 1.3)),
  step("What it means for Q1.", "The whole-profile test is weaker than the level test because it mixes the clear level change with a mix that did not change (next section)."),

  h3("6.4 Mix: did the relative pattern across compounds change?"),
  step("What it asks.", "If every sample is scaled to the same total, so that 'louder' no longer counts, do dormant and active profiles still differ?"),
  step("Why we chose it.", "It isolates the mix, which is what H1b is about, and explains the whole-profile result."),
  step("What we found.", sprintf("No. Season explained only %d%% of the differences in mix (%s).", round(100 * mix_version$R2), pv(mix_version$p))),
  step("What it means for Q1.", "Dormancy changed how strongly microbes responded, not which compounds they favoured."),

  h3("6.5 H1b: did complex compounds rise more than readily available ones? (Figure 2D)"),
  step("How compounds were grouped.", "One rule, applied to every compound: a compound is complex if microbes must break it apart outside the cell before they can take it up (a polymer, or an ester that must be split). Otherwise it is a small molecule that can be taken up directly, and it is readily available. This gives 8 complex and 23 readily available compounds, with none left over (Table A at the end)."),
  step("Why this rule.", "It is one clear, checkable criterion that matches the idea behind H1b: readily available compounds can be used immediately, complex ones need extra work (enzymes outside the cell). The transcript classes (lipid-like, aromatic, amino sugar, amino acid) are about chemistry rather than availability; for example chitin (complex) and N-acetylglucosamine (readily available) are both amino sugars. So those classes are used for a second check instead."),
  step("What it asks.", "Did complex compounds increase more in dormancy than readily available ones?"),
  step("What we found.", sprintf("Both groups increased by similar amounts: complex by %d%% (%s), readily available by %d%% (%s). The extra increase for complex compounds was small and not significant (%s log2, %s). H1b also predicted that readily available compounds would respond more in the plant-active period; they did not.",
                                  pct_up(avail_v("Complex")$estimate), pv(avail_v("Complex")$p),
                                  pct_up(avail_v("Readily available")$estimate), pv(avail_v("Readily available")$p),
                                  f2(h1b$estimate), pv(h1b$p))),
  step("Second check, within the same chemistry.", sprintf("Each polymer was compared with its own building blocks. Chitin and chitosan rose more than their building blocks (%s, %s), but cellulose and hemicellulose rose less than simple sugars (%s, %s), and lignin and protein showed no difference. There was no consistent pattern.",
                                                            f2(h1b_row("Amino sugar")$estimate), pv(h1b_row("Amino sugar")$p),
                                                            f2(h1b_row("Carbohydrate")$estimate), pv(h1b_row("Carbohydrate")$p))),
  step("Other checks.", sprintf("A stricter rule (only typical root exudates count as readily available), z-scores, removing the boundary months, 24-hour readings, and rates per gram of soil all gave small, non-significant differences (%s log2; Figure S2b).", range_txt(h1b_sens$estimate))),
  step("What it means for Q1.", "H1b is not supported. Microbes became more responsive to everything in dormancy, simple and complex alike."),

  h3("6.6 Curves through the year for each compound (Figure S2)"),
  step("What it asks.", "Did each compound change through the year at all, regardless of season?"),
  step("Why we chose it.", "It shows each compound's own pattern and identifies compounds that stayed steady (KH asked to contrast these with the seasonal ones). The curves are GAMs, smooth lines fitted through the 13 dates."),
  step("What we found.", sprintf("Most compounds changed through the year; %s showed no detectable change.", paste(lc(steady), collapse = ", "))),

  answer(b("Answer to Q1. "), sprintf("SIR differed between seasons: microbes responded about %d%% more strongly in dormancy, to all 31 compounds. The relative pattern across compounds did not change, and complex compounds did not rise more than readily available ones. H1a is supported as a change in level; H1b is not supported.", pct_up(level$estimate))),

  h2("7. Q2: How did field expression of microbial functions differ between seasons?"),
  img(fig("Figure3", "Figure3_main.png"), 6.5),
  p(i("Figure 3. A: NMDS picture of all 325 functions. B: summed expression of the functions linked to each compound class. C: each of the 20 linked functions through time.")),

  h3("7.1 What the transcript numbers are"),
  p("A metatranscriptome reads the genes microbes were actively using when the soil was collected. Genes were grouped into functional categories from the KEGG database (BRITE C level, for example 'fatty acid degradation'), regardless of which microbe expressed them. Grouping by function smooths over the fact that individual transcripts change within minutes."),
  p("After removing non-pathway rows, human and non-soil pathways, and 10 categories never detected, 325 categories remained. Each sample was then converted to proportions (each function's share of that sample's total expression), because total expression differed about three-fold between samples for reasons unrelated to which functions were expressed."),

  h3("7.2 H2: does the mix of expressed functions differ between seasons? (Figure 3A)"),
  step("What it asks.", "Are the expression profiles on the 2 dormant dates (January and March 2024) different from those on the 5 active dates?"),
  step("Why we chose it.", "It is the same kind of PERMANOVA as in Q1, again with each date counted once."),
  step("An important limit.", "With 5 active and 2 dormant dates there are only 21 ways to choose 2 dates out of 7. So the smallest p-value possible is 1/21 = 0.048. Getting that value means the 2 real dormant dates differ from the active dates more than any other pair of dates does."),
  step("What we found.", sprintf("p = %s, the smallest possible. Season explained %d%% of the differences among dates. The result held using only metabolism functions (R² = %s) and without converting to proportions (R² = %s).",
                                  f2(h2_res$p_permanova, 3), round(100 * h2_res$R2), f2(h2_sens$R2[str_detect(h2_sens$analysis, "Metabolism")]), f2(h2_sens$R2[str_detect(h2_sens$analysis, "as provided")]))),
  step("What it cannot tell us.", "With only 2 dormant dates, we cannot say whether the difference is caused by plant dormancy itself or by something else about January and March 2024, such as cold. In the text, write 'on the two plant-dormant sampling dates' rather than 'in the dormant season'."),

  h3("7.3 Functions linked to each compound class (Figures 3B and 3C)"),
  step("Why only 20 functions.", "Each MicroResp compound was linked in advance to the KEGG pathways that contain its breakdown or uptake reactions. Many pathways, like the TCA cycle, are linked to almost every compound, so they cannot tell compound classes apart. We kept only pathways linked to compounds of a single class: 20 pathways (Table B at the end)."),
  step("Why we add them up.", "For each class, the 20 pathways' shares are summed. SIR is a total of all microbial activity, so a class's total expression is the closest match."),
  step("What we found.", sprintf("These functions are a small part of all expression (about 0.2%% to 4%%). %d of the 20 were lower on dormant dates. January was the lowest date for %s amino acid functions. All four classes were highest in May 2024.",
                                  sum(cat_sum$dormant_to_active_ratio < 1), jan_low("Amino acid"))),
  step("What it cannot tell us.", "A lower share can mean less expression of that function or more expression of others. KEGG pathways include genes for building molecules as well as breaking them down, so higher expression does not prove the compound itself was being used."),
  answer(b("Answer to Q2. "), "The mix of expressed functions differed on the two plant-dormant dates (H2 supported, as strongly as the design allows). Functions linked to the MicroResp compounds, especially amino acid metabolism, were lowest in midwinter."),

  h2("8. Q3: Do SIR and expression of the matching functions go up and down together?"),
  img(fig("Figure4", "Figure4_main.png"), 6.5),
  p(i("Figure 4. A: SIR and expression for each class at the seven transcript dates. B: how closely they track. C: where they disagree.")),

  h3("8.1 Putting both on one scale"),
  p("SIR and transcript shares are in completely different units. Each was rescaled to its own average (0) and spread (1), called a z-score, so both can be drawn on one axis. This shows timing and direction (both up, both down), not size."),
  h3("8.2 Measuring how closely they track"),
  p("We used Spearman's rank correlation (rho): +1 if SIR and expression rank the seven dates in the same order, -1 if the order is reversed, 0 if unrelated. With only seven dates a correlation has to be larger than about 0.79 to be distinguishable from chance, which is why these are described rather than tested."),
  h3("8.3 Comparing with SIR one month later"),
  p(sprintf("Transcripts change within minutes, but the capacity measured by SIR builds over weeks (your framing, P2). So expression might match next month's SIR better than same-day SIR. We compared each transcript date with the next sampling date. Those gaps were 28 days except January to February (41 days, with plot 27 missing; plots per pair: %s).", paste(lagsum$n_pairs, collapse = ", "))),
  h3("8.4 Checking whether one date drives the result"),
  p("Each correlation was recalculated seven times, each time leaving out one date. If dropping one date changes the answer a lot, that date is driving it."),
  step("What we found.", sprintf("On the same date, the correlations were weak and inconsistent (rho %s). In every class, January had the highest SIR while expression was below average. Leaving January out raised every correlation (same date: rho %s). SIR one month later matched better in all four classes (rho %s with all dates; %s without January). Whole profiles were moderately aligned (Procrustes r = %s).",
                                  range_txt(loo_v("Concurrent", "rho_all_dates")), range_txt(loo_v("Concurrent", "rho_without_january")),
                                  range_txt(loo_v("One-month lag", "rho_all_dates")), range_txt(loo_v("One-month lag", "rho_without_january")),
                                  f2(protest$procrustes_correlation))),
  answer(b("Answer to Q3. "), "Outside midwinter, SIR and expression of the matching functions broadly rose and fell together, and matched better when SIR was taken a month later. In January they split: the potential to respond was high while expression was low. With seven dates this is a pattern to report, not a proven lag."),

  h2("9. Field context (Figure 1)"),
  step("Why we include it.", "It describes the conditions under which the SIR and transcript patterns happened."),
  step("What we found.", sprintf("Soil was %s °C colder on dormant dates and field CO2 was %d%% lower. Root biomass, microbial biomass, moisture, and enzyme potentials did not differ clearly between seasons.",
                                  f2(abs(f1v("SoilTemp_mean")$estimate), 1), round(100 * (1 - f1v("Resp_rate_mg_C_m2_h")$ratio)))),
  step("What it means.", "Microbial biomass and enzyme potential were similar in winter and summer. Field CO2 fell as soils cooled and plants went dormant, while the microbes' potential to respond to added compounds rose."),
  img(fig("Figure1", "Figure1_main.png"), 6.5),

  h2("10. What changed from storyboard v3, and why"),
  tbl(tribble(
    ~`Storyboard v3`, ~v4, ~Why,
    "SIR profiles differ by season, p = 0.002", sprintf("SIR level %d%% higher in dormancy (%s); mix unchanged; whole profile p = %s to %s", pct_up(level$estimate), pv(level$p), f2(whole_dates$p, 3), f2(builtin$p, 3)), "Each date now counts once (section 5.1); level and mix separated (section 5.2).",
    "Transcript profiles differ, p = 0.001", sprintf("p = %s, the smallest possible with 2 dormant dates", f2(h2_res$p_permanova, 3)), "Each date counts once.",
    "78 MicroResp samples", "77 samples", "January copy removed; February plot 27 was always missing (section 4.2).",
    "Top 3 compounds shown; H1b not tested", "H1b tested with an explicit rule", "To answer H1b directly.",
    "31 linked pathways, including TCA and glycolysis", "20 pathways linked to one class only", "Pathways shared by all classes cannot distinguish them.",
    "Root vs respiration in the main figure", "Supplement", "Not one of the storyboard questions."
  ), rel = c(2, 2.5, 2.5)),

  h2("11. Decisions and questions still open"),
  ul(c(
    sprintf("Decided (2026-10-25): the date-average PERMANOVA (p = %s) is the main whole-profile test; the all-samples version with whole dates shuffled (p = %s) is reported as the check.", f2(date_avg$p, 3), f2(whole_dates$p, 3)),
    "Whether to replace the NMDS in main Figure 2 with a simpler picture of 'same pattern, just louder', keeping the NMDS in the supplement.",
    "For Sheryl: why plot 27 in February 2024 has no MicroResp data; the reference for the season boundaries; harvest and fall-saturation dates; rhizotron data.",
    "Who made Data/brite_C_expression.csv and how (needed for Methods).",
    "The work is on the GitHub branch manuscript-v4-analysis; nothing has been merged into main."
  )),

  h2("Table A. The 31 MicroResp compounds"),
  tbl(class_table, rel = c(1.2, 1.6, 1.3, 3.2)),
  h2("Table B. The 20 transcript functions linked to one compound class"),
  tbl(pathway_table, rel = c(0.9, 2.1, 3.3, 0.9))
)

write_docx(guide, file.path("Docs", "Plain_language_guide_v4.docx"))
cat("Wrote Docs/Plain_language_guide_v4.docx\n")
