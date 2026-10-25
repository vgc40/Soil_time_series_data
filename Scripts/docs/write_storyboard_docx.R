# =============================================================================
# Write the v4 storyboard and analysis notes as Word documents (.docx)
# =============================================================================
#
# Every number in the documents is read from Results/manuscript_v4, so the
# text stays in sync with the analysis. Run after Scripts/00_run_all.R:
#   Rscript Scripts/docs/write_storyboard_docx.R
# The storyboard keeps the v3 sections; framing, questions, and hypotheses
# are reproduced exactly as written by VGC.

rm(list = ls())
suppressPackageStartupMessages(library(tidyverse))
source(file.path("Scripts", "docs", "docx_writer.R"))

out_dir <- "Docs"
dir.create(out_dir, showWarnings = FALSE)
res <- function(...) read_csv(file.path("Results", "manuscript_v4", ...), show_col_types = FALSE)
fig <- function(...) normalizePath(file.path("Figures", "manuscript_v4", ...), winslash = "/")

f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
pv <- function(p) if_else(p < 0.001, "p < 0.001", paste0("p = ", f2(p, 3)))
pct <- function(log2_value) sprintf("%+.0f%%", 100 * (2^log2_value - 1))
ci <- function(lo, hi, d = 2) sprintf("95%% CI %s to %s", f2(lo, d), f2(hi, d))

# -----------------------------------------------------------------------------
# Results used in the text
# -----------------------------------------------------------------------------

design <- res("sample_availability", "design_summary.csv")
dval <- function(item) design$value[design$item == item]
lagsum <- res("sample_availability", "lag_pair_summary.csv")

f1 <- res("Figure1", "Figure1_season_tests.csv")
f1v <- function(col) filter(f1, column == col)
enzyme_ratios <- f1$ratio[str_detect(f1$column, "Sum_")]

level <- res("Figure2", "Figure2A_overall_level_shift.csv")
var_share <- res("Figure2", "Figure2A_variance_shares_date_vs_plot.csv")
share_of <- function(src) round(100 * var_share$share[var_share$source == src])
level_dates <- res("Figure2", "Figure2A_overall_level_by_date.csv")
level_samples <- res("Figure2", "Figure2A_overall_level_by_sample.csv")
july <- filter(level_samples, Sampling_Date == as.Date("2023-07-05"))
july_mean <- level_dates$mean[level_dates$Sampling_Date == as.Date("2023-07-05")]
profile_tests <- res("Figure2", "Figure2B_season_PERMANOVA_profile_and_mix.csv")
whole <- filter(profile_tests, str_detect(analysis, "Whole"), level == "Date means (primary)")
mix <- filter(profile_tests, str_detect(analysis, "Mix"), level == "Date means (primary)")
sub_eff <- res("Figure2", "Figure2C_substrate_season_effects.csv")
h1b_tests <- res("Figure2", "Figure2D_H1b_tests.csv")
h1b_row <- function(cmp) filter(h1b_tests, comparison == cmp)
h1b <- h1b_row("All substrates: complex vs readily available")
h1b_perm <- res("Figure2", "Figure2D_H1b_permutation_check.csv")
by_avail <- res("Figure2", "Figure2D_season_effect_by_availability.csv")
h1b_sens <- res("Figure2", "FigureS2b_H1b_sensitivity.csv")
envfit <- res("Figure2", "Figure2_envfit_date_level.csv")
gam <- res("Figure2", "FigureS2_GAM_temporal_statistics.csv")

h2_all <- res("Figure3", "Figure3A_season_PERMANOVA_dispersion.csv")
h2_res <- filter(h2_all, level == "Date means (primary)")
h2_sens <- res("Figure3", "FigureS3_H2_sensitivity.csv") %>% filter(level == "Date means (primary)")
cat_sum <- res("Figure3", "Figure3B_unique_category_summary.csv")
share <- res("Figure3", "Figure3C_class_share_by_date.csv")
cat_dates <- res("Figure3", "Figure3B_unique_category_date_means.csv")
extremes <- cat_dates %>%
  group_by(compound_class, brite_C) %>%
  summarise(min_date = Sampling_Date[which.min(z)], max_date = Sampling_Date[which.max(z)],
            .groups = "drop")
jan_low <- function(cls) {
  e <- filter(extremes, compound_class == cls)
  sprintf("%d of %d", sum(e$min_date == as.Date("2024-01-02")), nrow(e))
}
lipid_low_2023 <- with(filter(extremes, compound_class == "Lipid-like"), sum(min_date < as.Date("2023-11-01")))
lipid_high_2024 <- with(filter(extremes, compound_class == "Lipid-like"), sum(max_date >= as.Date("2024-03-01")))

loo <- res("Figure4", "Figure4B_leave_one_date_out_summary.csv")
loo_v <- function(cmp, col) loo[[col]][loo$comparison == cmp]
div <- res("Figure4", "Figure4_class_values_by_date.csv") %>% filter(comparison == "Concurrent")
jan_div <- filter(div, Sampling_Date == as.Date("2024-01-02"))
protest <- res("Figure4", "FigureS4b_PROTEST.csv")

s1 <- res("FigureS1", "FigureS1_root_respiration_slopes.csv")
s1v <- function(primary, season) filter(s1, str_detect(analysis, "primary") == primary, Season == season)

classes <- read_csv(file.path("Data", "compound_class_mapping.csv"), show_col_types = FALSE)
pathways <- read_csv(file.path("Data", "compound_class_unique_BRITE_C.csv"), show_col_types = FALSE)

# Lower-case a substrate name for use mid-sentence (keeps D-, L-, N-, BSA, GABA).
lc <- function(x) if_else(str_detect(x, "^[A-Z][a-z]"), paste0(tolower(str_sub(x, 1, 1)), str_sub(x, 2)), x)
range_txt <- function(x, d = 2) sprintf("%s to %s", f2(min(x), d), f2(max(x), d))
n_pos <- sum(sub_eff$estimate > 0)
top_sub <- sub_eff %>% arrange(desc(estimate)) %>% slice(1:3)
steady <- sub_eff %>% filter(steady) %>% pull(display_name)
lower_dormant <- sum(cat_sum$dormant_to_active_ratio < 1)
share_range <- share %>% group_by(compound_class) %>% summarise(lo = min(mean), hi = max(mean))
gwc_env <- filter(envfit, variable == "Gravimetric water content")
lowest_dates <- level_dates %>% arrange(mean) %>% slice(1:2) %>% pull(Sampling_Date) %>% format("%B %Y")

# -----------------------------------------------------------------------------
# Tables
# -----------------------------------------------------------------------------

group_order <- c("Carbohydrate", "Organic acid", "Amino acid", "Amino sugar", "Aromatic",
                 "Lipid-like", "Aliphatic amine")
table1 <- classes %>%
  mutate(substrate_group = factor(substrate_group, levels = group_order)) %>%
  arrange(substrate_group, desc(availability_class == "Complex"), display_name) %>%
  transmute(
    `Chemical class` = as.character(substrate_group),
    Substrate = display_name,
    `H1b group` = availability_class,
    `What it is` = substrate_note,
    `Linked to transcripts (Q2, Q3)` = if_else(is.na(compound_class), "No", "Yes")
  )

table2 <- pathways %>%
  left_join(select(cat_sum, brite_C, mean_share_pct), by = "brite_C") %>%
  mutate(compound_class = factor(compound_class, levels = c("Lipid-like", "Aromatic", "Amino sugar", "Amino acid"))) %>%
  arrange(compound_class, desc(mean_share_pct)) %>%
  transmute(
    Class = as.character(compound_class),
    `KEGG pathway` = sprintf("ko%s %s", brite_C_ID, brite_C_name),
    `What the pathway does` = what_it_does,
    `Why it is linked to the class` = why_linked,
    `Share of all expression` = paste0(formatC(mean_share_pct, format = "fg", digits = 2), "%")
  )

glance <- tribble(
  ~Hypothesis, ~`What we asked the data`, ~Result, ~Verdict,
  "H1a", "Do substrate responses differ between plant-active and plant-dormant dates? Split into level (do all substrates respond more?) and mix (does the relative pattern change?).",
  sprintf("Level: %s in dormancy (%s). Mix: no change (R² = %s, %s). Whole profile: R² = %s, %s.",
          pct(level$estimate), pv(level$p), f2(mix$R2), pv(mix$p_permanova), f2(whole$R2), pv(whole$p_permanova)),
  "Supported for level: microbes responded more strongly to every substrate during dormancy. Not supported for mix.",
  "H1b", "Did complex substrates gain more in dormancy than readily available ones (and readily available ones more in the plant-active period)?",
  sprintf("Complex gained %s log2 more (%s). Readily available substrates also rose in dormancy (%s). Within classes the difference was inconsistent.",
          f2(h1b$estimate), pv(h1b$p), pct(by_avail$estimate[by_avail$availability == "Readily available"])),
  "Not supported.",
  "H2", "Does the composition of expressed functions differ between plant-active and plant-dormant dates?",
  sprintf("R² = %s, %s (smallest p possible with two dormant dates).", f2(h2_res$R2), pv(h2_res$p_permanova)),
  "Supported, as strongly as seven dates allow.",
  "Q2 (classes)", "How did the pathways linked to each compound class change? (descriptive)",
  sprintf("%d of 20 pathways had a lower share on dormant dates; %s amino acid pathways lowest in January; all four classes highest in May 2024.", lower_dormant, jan_low("Amino acid")),
  "Descriptive.",
  "Q3", "Do SIR and expression of the matching pathways rise and fall together? (descriptive)",
  sprintf("Same date: rho %s. SIR one month later: rho %s, higher in every class. Without January: same date %s.",
          range_txt(loo_v("Concurrent", "rho_all_dates")), range_txt(loo_v("One-month lag", "rho_all_dates")),
          range_txt(loo_v("Concurrent", "rho_without_january"))),
  "Descriptive: they diverge in midwinter and track each other better with a one-month offset."
)

glossary <- tribble(
  ~Term, ~Meaning,
  "SIR (substrate-induced respiration)", "CO2 released by a soil sample in the 5 h after a compound is added. A measure of how strongly the existing microbes can respond to that compound (potential), not what they do in the field.",
  "Carbon-normalized", "SIR divided by the amount of carbon added with the substrate, so substrates added at different carbon doses can be compared.",
  "log2 and fold change", "log2 values turn ratios into differences: +1 = twice as much, -1 = half, 0 = no change. A dormant-minus-active difference of 0.42 means dormant values were 2^0.42 = 1.34 times active values (34% higher).",
  "Annual mean (for each substrate)", "The average of the 13 date means for that substrate. Expressing each value relative to it removes the large differences in typical rate between substrates.",
  "z-score", "A value rescaled to its own average (0) and spread (1). Used when two very different measurements (SIR and transcripts) are compared on one axis: only timing and direction are comparable, not size.",
  "NMDS", "A picture that places samples so that similar samples are close together. It shows patterns but is not a test.",
  "Bray-Curtis", "The measure of how different two samples' profiles are, used for NMDS and PERMANOVA. It reflects differences in both overall size and relative pattern.",
  "PERMANOVA and R²", "The test behind the NMDS picture: how much of the difference among samples is explained by season (R², 0-1), and how often a random relabelling of dates gives as strong a result (p).",
  "Exact p-value", "With few dates, every possible way of labelling dates as active or dormant can be listed. The p-value is the share of those labellings that give a result at least as strong as the real one. With 5 active and 2 dormant transcript dates there are only 21 labellings, so p cannot be smaller than 1/21 = 0.048.",
  "Dispersion (spread)", "Whether samples in one season vary more among themselves than in the other. PERMANOVA can respond to spread as well as to a real shift, so spread is tested alongside it.",
  "Mixed model", "A regression that knows which samples share a date and which share a plot, so it does not count the six plots on one date as six independent observations of that season.",
  "FDR (false discovery rate)", "A correction for testing 31 substrates at once, so that chance results are not mistaken for real ones.",
  "GAM", "A smooth curve through time fitted to each substrate, used to ask whether a substrate changed through the year at all.",
  "BRITE C / KEGG pathway", "A functional category in the KEGG database (e.g., fatty acid degradation). Gene transcripts were grouped into these categories regardless of which microbe expressed them.",
  "Spearman rho", "A correlation based on ranks: +1 when two measurements rank the dates in the same order, -1 when the order is reversed, 0 when unrelated.",
  "Leave-one-date-out", "Repeating a calculation seven times, each time without one date, to see whether a single date drives the result.",
  "Procrustes", "Compares two whole-profile pictures (SIR and transcripts) to see whether samples sit in similar positions in both."
)

# -----------------------------------------------------------------------------
# Storyboard
# -----------------------------------------------------------------------------

storyboard <- c(
  h1("Storyboard v4: seasonal substrate-use potential and field gene expression"),
  p('<span class="note">', b("Note: "),
    "v4 was drafted with Claude (Anthropic) from the analysis in the GitHub repository (branch manuscript-v4-analysis); earlier versions used ChatGPT 5.6 Sol (enterprise). ",
    "Every number below is read directly from the result tables written by the scripts in Scripts/ (run order in Scripts/00_run_all.R). ",
    "The framing, research questions, and hypotheses are reproduced unchanged. Decisions, data changes, and open items are in the companion document Analysis_notes_v4.docx.</span>"),

  h2("Suggested Title"),
  p("Option 1: Plant dormancy strengthens soil microbial responses to added substrates while field expression of related functions dips in midwinter"),
  p("Option 2: Microbial substrate-use potential and field gene expression follow different seasonal paths under a perennial grass"),
  p(i("v3 title for reference: Soil microbial substrate-induced respiration and field transcript expression show seasonal patterns but limited correspondence.")),

  h2("Potential Journals"),
  p("SBB"),
  h2("Team members"),
  p("Vanessa Garayburu-Caruso, Josue Rodriguez-Ramos, Sheryl Bell, Evan Warburton, Sharon Zhao, Regan McDearis, David Degnan, Amy Zimmerman, and Kirsten Hofmockel"),
  h2("Acknowledgements"),
  p("Gregory, Kali Cornwell, Mischelle, other people that contributed? Trinidad (viromes; not in this paper)"),
  h2("GitHub link with code"),
  p("vgc40/Soil_time_series_data (https://github.com/vgc40/Soil_time_series_data); v4 analysis on branch manuscript-v4-analysis"),

  h2("Framing"),
  h3("P1. Seasonal changes in microbial substrate availability and function"),
  p("Soil microbial function varies across the annual cycle as plant activity, temperature, moisture, and resource availability change. As plants transition between active growth and dormancy, changes in root exudation, root turnover, and residue inputs may alter the timing and relative availability of different substrates, including readily available compounds and more complex organic matter associated with plant residues and microbial necromass. Dormancy may therefore change the mixture of substrates available to microbes rather than simply reducing their supply. These shifts can alter microbial physiological state and resource acquisition strategies, potentially affecting how microbes respond to different substrates throughout the year. Yet much of our understanding of soil microbial function comes from growing-season measurements, despite evidence that microbial activity can continue outside the growing season (Monson et al., 2006; Gavazov et al., 2017). As a result, we know comparatively little about how microbial responses to different substrates change during plant dormancy, limiting our understanding of how plant seasonality shapes belowground microbial metabolism."),
  h3("P2. Connecting microbial substrate-use potential and field expression"),
  p("Commonly used measurements capture different aspects and timescales of microbial function. Field soil respiration integrates root and microbial activity under prevailing conditions but cannot resolve which substrates microbes are using. Laboratory substrate-induced respiration measures short-term responses to added compounds under standardized conditions. Because these responses are measured within 5 h, they primarily characterize the short-term respiratory capacity of the existing microbial community, which reflects its physiological state and antecedent environmental conditions rather than substantial growth during the assay (Anderson & Domsch, 1978; Campbell et al., 2003). Seasonal changes in these responses therefore show how microbial substrate-use potential shifted over the year. Metatranscriptomes, in contrast, characterize which microbial functions are being expressed under field conditions at the time of sampling, but not how strongly microbes would respond to specific substrate additions. It remains unclear how seasonal changes in field gene expression relate to changes in substrate-use potential, in part because the two capture microbial activity over different timescales: microbial transcripts typically persist for minutes (Moran et al., 2013), whereas soil extracellular enzymes can remain active for weeks after they are produced (Schimel et al., 2017). Seasonal changes in field expression and substrate-use potential may therefore not coincide. Comparing the two across related compounds and metabolic pathways can help clarify whether plant dormancy changes which functions microbes express, their potential to respond as substrate availability changes, or both."),
  h3("P3. Study objectives and research questions"),
  p("Here, we investigate how microbial substrate-use potential and field gene expression change between plant-active and plant-dormant periods, and how the two are related. We focused on fully irrigated (100%) plots in an irrigated tall wheatgrass system, a perennial grass with a plant-active period from April (plant emergence) through October (harvest) and a plant-dormant period from November through March [basis for season definitions]. We combine 13 months of substrate-induced respiration measurements for 31 soil-relevant organic compounds with field metatranscriptomes from seven of these sampling dates (June, July, September, and October 2023; January, March, and May 2024), and use supporting measurements of extracellular enzyme activity, microbial biomass, soil respiration, root biomass, and environmental and biogeochemical conditions to help interpret the observed patterns. Specifically, we ask (1) how substrate-induced respiration differs between plant-active and plant-dormant periods, and whether seasonal responses differ among substrates; (2) how functional transcript expression differs between these periods, and whether functions associated with different compound classes show distinct seasonal patterns; and (3) how seasonal changes in substrate-induced respiration within each compound class relate to the expression of functions associated with that class."),

  h2("Research Questions & Corresponding Hypotheses"),
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

  h2("Results at a glance"),
  p("Each row links a hypothesis to what we asked of the data, what we found, and the verdict. The sections below explain how each result was reached."),
  tbl(glance, rel = c(0.7, 2.2, 2.6, 1.8)),

  h2("Approach"),
  p("This section explains what was measured, how each measurement was summarised, and why each analysis was chosen. Technical terms are defined in the glossary at the end."),

  h3("1. Samples and seasons"),
  p(sprintf("Soil was sampled monthly from June 2023 to June 2024 (13 dates) in the six fully irrigated (100%%) plots of the tall wheatgrass trial, three cores per plot pooled into one sample. Field, enzyme, and microbial biomass data exist for all %s plot-dates; MicroResp for %s (plot 27 was not assayed on 12 February 2024); metatranscriptomes for %s samples (five plots on seven dates; plot 57 was not sequenced). Dates from April (emergence) to October (harvest) are plant-active and dates from November to March plant-dormant. This gives 8 active and 5 dormant dates for MicroResp and 5 active and 2 dormant dates for transcripts.",
            dval("Observed 100% field records"), dval("Observed 100% MicroResp samples"), dval("Observed 100% transcript samples"))),

  h3("2. Why each sampling date, not each plot, counts as one replicate of season"),
  p("Season is not assigned to plots; it comes with the sampling date. All six plots sampled on the same day share that day's weather, the days since the last irrigation, the same field trip, and the same laboratory assay run. If something makes one date unusual, all six plots are unusual together."),
  p(sprintf("The data show this directly. Of the variation in a sample's overall SIR level, %d%% is shared by all plots on the same date, while only %d%% reflects consistent differences between plots (the remaining %d%% is plot-to-plot noise within a date). In July 2023, for example, %d of %d plots were below their substrate averages and the date as a whole was %s below average. Counting the six plots as six independent pieces of evidence would turn one unusual date into six apparent confirmations of a seasonal pattern. This problem is known as pseudoreplication (Hurlbert 1984).",
            share_of("Sampling date"), share_of("Plot"), share_of("Plot x date (residual)"),
            sum(july$response < 0), nrow(july), pct(july_mean) %>% str_remove("^-"))),
  p("All plots are still used. The mixed models include every plot and are told which plots share a date, so plots inform how variable the soil is on a given day, while dates carry the evidence about seasons. The multivariate tests compare date averages and report the plot-level version as a check."),

  h3("3. What MicroResp measures and how it was expressed"),
  p("MicroResp adds one compound to a small amount of soil (31 compounds plus water, three wells each, soil at 60% water-holding capacity, room temperature) and measures the CO2 released over 5 h. Because the assay is short and standardized, it measures the respiration potential of the microbes already present when the soil was collected, not what they were doing in the field. Rates are expressed per unit of carbon added."),
  p("Substrates differ greatly in their typical rates (simple sugars respire much faster than lignin). To compare seasonal change across substrates, as H1b specifies, each value was expressed relative to that substrate's own annual mean on a log2 scale: 0 means the substrate's usual rate, +1 twice its usual rate, -1 half. The seasonal difference for a substrate is then the dormant minus active value on this scale (a log2 fold change)."),

  h3("4. Two ways SIR can differ between seasons: level and mix"),
  p("A season could change SIR in two different ways. In a level change, microbes respond more (or less) to every substrate, like turning up the volume on all instruments. In a mix change, some substrates rise relative to others, like changing which instruments are loud. H1a ('profiles differ') could mean either, and H1b is specifically about the mix (complex substrates rising relative to readily available ones). We therefore tested them separately:"),
  ul(c(
    "Level: each sample's average response over all 31 substrates, compared between seasons with a mixed model (dates as replicates).",
    "Mix: each sample's 31 rates divided by its own total, so a sample that is uniformly higher has an identical profile; profiles then compared between seasons with PERMANOVA on date means.",
    "Whole profile (level and mix together): PERMANOVA on the raw rates, shown as an NMDS picture."
  )),

  h3("5. How substrates were grouped"),
  p("Each substrate has two labels (Table 1)."),
  p(b("H1b group (availability). "), "One explicit rule: a substrate is complex if microbes must break it apart outside the cell before they can take it up (a polymer, or an ester that must be split); otherwise it is a small molecule that can be taken up directly and is readily available. This gives 8 complex substrates (cellulose, hemicellulose, lignin, chitin, chitosan, the protein BSA, and the lipid esters methyl oleate and Tween 40) and 23 readily available ones. A single rule leaves no substrate unassigned. As a check, a stricter rule that counts only typical root exudates (sugars, organic acids, amino acids) as readily available gave the same answer (Fig. S2b)."),
  p(b("Chemical class. "), "Substrates were also grouped by chemistry, using the same classes as the transcript analysis (lipid-like, aromatic, amino sugar, amino acid) plus carbohydrates, organic acids, and one aliphatic amine. This matters for H1b: complex and readily available substrates also differ in chemistry, so we also compared each polymer with its own building blocks within a class (chitin and chitosan versus N-acetylglucosamine and glucosamine; lignin versus phenolic compounds; cellulose and hemicellulose versus simple sugars; protein versus amino acids). This holds chemistry constant and isolates availability."),

  h3("6. Transcripts"),
  p("Gene transcripts from the field samples were grouped into BRITE C functional categories (KEGG pathways), regardless of which microbe expressed them. Individual transcripts turn over in minutes, so grouping by function focuses on community-level investment rather than momentary signals of single genes (AZ comment). After removing hierarchy rows, human and non-soil pathways, and categories never detected, 325 categories remained. Each sample was converted to proportions (each category's share of that sample's total expression), because total expression differed about three-fold between samples for reasons unrelated to composition."),
  p("H2 was tested with PERMANOVA on date means. With 5 active and 2 dormant dates there are only 21 ways to choose two 'dormant' dates, so the smallest possible p-value is 1/21 = 0.048. A result at that value means the two dormant dates differ from the active dates more than any other pair of dates does, which is the strongest statement this design can make. It cannot separate plant dormancy from anything else unique to January and March 2024, such as cold."),

  h3("7. Linking substrates to transcripts"),
  p("An a priori table links each MicroResp substrate to the KEGG pathways that contain its breakdown or uptake reactions (Data/Microresp_pathways.csv). Many pathways, such as the TCA cycle or glycolysis, are linked to almost every substrate and cannot tell compound classes apart. We therefore used only pathways linked exclusively to substrates of one class: 20 pathways in total (Table 2). For each class, the expression of its pathways was summed, because a class's total investment is the closest analogue of SIR, which also sums all microbial activity. KEGG pathways contain both breakdown and building genes, so higher expression means more activity in that part of metabolism, not proof that the substrate itself was being used."),

  h3("8. Comparing SIR and transcripts (Q3)"),
  p(sprintf("For each compound class, SIR (mean over the class's substrates) and transcript abundance (summed class pathways) were rescaled to z-scores across the seven transcript dates, so the two can be drawn on one axis and compared for timing and direction. Spearman rank correlations were calculated between them (descriptive; with seven dates a correlation must exceed about 0.79 to be statistically distinguishable from zero). Because transcripts reflect the moment of sampling while the capacity measured by SIR builds over weeks (P2), transcripts were also compared with SIR measured at the next sampling date. Those intervals were 28 days except January to February (41 days; plot 27 missing; n plots per pair = %s). Each correlation was recalculated with each date left out in turn, to see whether one date drives it.", paste(lagsum$n_pairs, collapse = ", "))),

  h3("9. Supporting measurements"),
  p("Field soil respiration (chamber CO2 flux), soil temperature, gravimetric water content, root biomass hand-picked from the pooled 0-15 cm cores (not rhizosphere material), microbial biomass C and N, and potential activities of seven extracellular enzymes (summed as C-, N-, and P-acquiring). These describe the field conditions in which the SIR and transcript patterns arose."),

  h2("Key Figures"),

  h3("Figure 1. What were field conditions like in each season?"),
  img(fig("Figure1", "Figure1_main.png")),
  p('<span class="caption">', b("How to read Figure 1. "),
    "(A) Sampling design: each circle is a sampling date, with the number of plots sampled; the lower row marks the seven dates with metatranscriptomes. Green shading is the plant-active period, blue-grey the plant-dormant period; the blue bar is the drip-irrigation window. (B-I) Field variables through time. Grey points are plots; the black line and grey band are the date average and its 95% confidence interval; log axes are used where values span a wide range. Each subtitle gives the dormant-to-active ratio (or difference, for temperature) and its p-value, with dates as replicates.</span>"),
  p(b("What we found. ")),
  ul(c(
    sprintf("Soil temperature was %s °C lower on dormant dates.", f2(abs(f1v("SoilTemp_mean")$estimate), 1)),
    sprintf("Field soil respiration on dormant dates was %d%% of active-period rates (%s).", round(100 * f1v("Resp_rate_mg_C_m2_h")$ratio), pv(f1v("Resp_rate_mg_C_m2_h")$p)),
    sprintf("Root biomass (%s×), microbial biomass C (%s×), water content (%s×), and potential enzyme activities (%s-%s×) did not differ clearly between periods (all p > 0.15).",
            f2(f1v("Roots_gdw")$ratio), f2(f1v("MBC_ug_gdw_bgc")$ratio), f2(f1v("GWC")$ratio), f2(min(enzyme_ratios)), f2(max(enzyme_ratios)))
  )),
  p(b("What this means. "), "The fall in field CO2 during dormancy happened without a loss of microbial biomass or of the enzyme pool. The microbial community was still present in similar amounts in winter; what changed was the conditions it was working under (cold, no plant activity). This is the backdrop for the SIR results: any change in SIR is not simply because there were more or fewer microbes."),
  p(b("What this cannot tell us. "), "Field respiration combines root and microbial respiration and cannot separate them. Harvest and fall-saturation dates still need to be added to panel A."),

  h3("Figure 2. Q1: How did microbial responses to added substrates change with season?"),
  img(fig("Figure2", "Figure2_main.png")),
  p('<span class="caption">', b("How to read Figure 2. "),
    "(A) Level: the average SIR over all 31 substrates for each sample (grey points) and date (black line, 95% band), relative to each substrate's annual mean; above 0 means stronger responses than usual. (B) Whole profile: NMDS picture in which samples with similar 31-substrate profiles sit close together; small points are samples, large points are date averages (green active, blue-grey dormant). (C) Each substrate's dormant minus active difference with its 95% CI, grouped by chemical class and coloured by H1b group; right of 0 means higher in dormancy, +1 means twice as high; open symbols are not significant after correcting for 31 tests; 'steady' marks substrates with no detectable change through the year. (D) H1b: how much more complex substrates rose in dormancy than readily available ones, for all substrates and within each chemical class; right of 0 supports H1b.</span>"),
  p(b("What we found, step by step. ")),
  ul(c(
    sprintf("Panel A: average SIR sat above its annual level through the dormant months and was lowest in %s. With dates as replicates, SIR was %s higher on dormant dates (log2 %s, %s, %s).",
            paste(lowest_dates, collapse = " and "), str_remove(pct(level$estimate), "^\\+"), f2(level$estimate), ci(level$lower, level$upper), pv(level$p)),
    sprintf("Panel C: all %d of %d substrates had higher SIR in dormancy, but no single substrate was significant after correcting for 31 tests. The increase was shared by all substrates rather than driven by a few. The largest increases were for %s.",
            n_pos, nrow(sub_eff), paste(sprintf("%s (%s×)", lc(top_sub$display_name), f2(2^top_sub$estimate, 1)), collapse = ", ")),
    sprintf("Panel B: because the increase was shared, the relative mix of substrate responses did not change with season (mix-only test: R² = %s, %s). The whole-profile test, which combines level and mix, was intermediate (R² = %s, %s). In the NMDS, dormant and active samples overlap.",
            f2(mix$R2), pv(mix$p_permanova), f2(whole$R2), pv(whole$p_permanova)),
    sprintf("Panel D: complex substrates rose only slightly more than readily available ones (%s log2 extra, %s; permutation check %s), and readily available substrates did not respond more in the plant-active period; they also rose in dormancy (%s, %s). Within chemical classes the comparison was inconsistent: chitin and chitosan rose more than their monomers (%s, %s), but cellulose and hemicellulose rose less than simple sugars (%s, %s).",
            f2(h1b$estimate), pv(h1b$p), pv(h1b_perm$p_two_sided),
            pct(by_avail$estimate[by_avail$availability == "Readily available"]), pv(by_avail$p[by_avail$availability == "Readily available"]),
            f2(h1b_row("Amino sugar")$estimate), pv(h1b_row("Amino sugar")$p),
            f2(h1b_row("Carbohydrate")$estimate), pv(h1b_row("Carbohydrate")$p))
  )),
  p(b("Verdict. "), "H1a is supported as a level change: microbes responded more strongly to every substrate during dormancy. H1b is not supported: the change was in how strongly microbes responded, not in which substrates they favoured."),
  p(b("What this means. "), "In winter the soil community was primed to respond more strongly to any added carbon, simple or complex. Possible explanations, which this design cannot separate, include communities that are more carbon-limited when plants stop supplying root inputs and so respond more when carbon is added, better-preserved active cells or enzymes in cold soil, or a stronger response of winter soil to being warmed and wetted in the assay (a Birch-like effect). Water content at sampling was the field variable most aligned with the SIR ordination ",
    sprintf("(date-level envfit r² = %s, %s; exploratory).", f2(gwc_env$r2), pv(gwc_env$p))),
  p(b("What this cannot tell us. "), "SIR is measured at room temperature and 60% water-holding capacity, so it shows potential, not field rates. Thirteen dates give limited power; the level result (p = 0.045) should be reported with its confidence interval."),

  h3("Figure 3. Q2: How did field expression of microbial functions change with season?"),
  img(fig("Figure3", "Figure3_main.png")),
  p('<span class="caption">', b("How to read Figure 3. "),
    "(A) NMDS picture of all 325 functional categories: samples with similar expression profiles sit close together; small points are samples and large points date averages. (B) For each compound class, the summed share of all transcript expression in its uniquely linked pathways (Table 2); points are plots, the line the date average with 95% band. (C) Each of the 20 uniquely linked pathways at each date, coloured by how far it was above (red) or below (blue) its own average; this shows timing, not size, and the share of all expression is given in brackets.</span>"),
  p(b("What we found, step by step. ")),
  ul(c(
    sprintf("Panel A: the two dormant dates (January and March 2024) sit apart from the 2023 plant-active dates. Season explained %d%% of the differences among date averages (R² = %s, %s), the strongest result possible with two dormant dates; samples were equally spread in both seasons (spread %s). The result held using only metabolism categories (R² = %s) and without within-sample normalization (R² = %s).",
            round(100 * h2_res$R2), f2(h2_res$R2), pv(h2_res$p_permanova), pv(h2_res$p_dispersion),
            f2(h2_sens$R2[str_detect(h2_sens$analysis, "Metabolism")]), f2(h2_sens$R2[str_detect(h2_sens$analysis, "as provided")])),
    "May 2024 sits between the 2023 active dates and the dormant dates, so one year later expression had not returned to its June 2023 state.",
    sprintf("Panel B: the pathways linked to each class are a small part of all expression (lipid-like %s%%, aromatic %s%%, amino sugar %s%%, amino acid %s%% of the total). Amino acid pathways dropped sharply in January; all four classes reached their highest share in May 2024.",
            range_txt(share_range$lo[share_range$compound_class == "Lipid-like"] %>% c(share_range$hi[share_range$compound_class == "Lipid-like"])),
            range_txt(c(share_range$lo[share_range$compound_class == "Aromatic"], share_range$hi[share_range$compound_class == "Aromatic"])),
            range_txt(c(share_range$lo[share_range$compound_class == "Amino sugar"], share_range$hi[share_range$compound_class == "Amino sugar"])),
            range_txt(c(share_range$lo[share_range$compound_class == "Amino acid"], share_range$hi[share_range$compound_class == "Amino acid"]), 1)),
    sprintf("Panel C: %d of the 20 pathways had a lower share of expression on dormant dates. January was the lowest date for %s amino acid pathways and %s aromatic pathways. The lipid-like pathways followed a different course: %d of 4 were lowest in summer 2023 and %d of 4 highest in March-May 2024.",
            lower_dormant, jan_low("Amino acid"), jan_low("Aromatic"), lipid_low_2023, lipid_high_2024)
  )),
  p(b("Verdict. "), "H2 is supported: the composition of expressed functions differed on the two plant-dormant dates. Descriptively, expression of the substrate-related pathways, especially amino acid metabolism, dipped in midwinter."),
  p(b("What this means. "), "While the community's potential to respond to added substrates rose in winter (Figure 2), its field investment in the matching metabolic pathways fell, most clearly in January."),
  p(b("What this cannot tell us. "), "With two dormant dates we cannot say whether the difference reflects plant dormancy itself or other features of January and March 2024, such as cold. Shares are relative: a lower share can mean less expression of that pathway or more expression of others. Pathways include building as well as breakdown genes."),

  h3("Figure 4. Q3: Did SIR and expression of matching pathways rise and fall together?"),
  img(fig("Figure4", "Figure4_main.png")),
  p('<span class="caption">', b("How to read Figure 4. "),
    "(A) For each class, SIR (brown, mean of the class's substrates) and transcript abundance (blue, summed class pathways) at the seven transcript dates, both rescaled to their own average (0) and spread (1) so only timing and direction are comparable; rho is the rank correlation. (B) Rank correlations for same-date values (black) and for transcripts against SIR at the next sampling date (orange); dots use all dates, bars show the range when any one date is removed, and the X marks the value without January. (C) SIR minus transcript z-score at each date; brown means SIR was relatively higher than expression.</span>"),
  p(b("What we found, step by step. ")),
  ul(c(
    sprintf("Panel A: in every class, January had the highest SIR of the seven dates while expression of the matching pathways was below average (difference %s to %s in panel C).",
            f2(min(jan_div$divergence), 1), f2(max(jan_div$divergence), 1)),
    sprintf("Panel B: same-date correlations were weak and inconsistent (rho %s), none close to the 0.79 needed to be distinguishable from zero with seven dates.", range_txt(loo_v("Concurrent", "rho_all_dates"))),
    sprintf("Leaving dates out showed January was the main source of mismatch: without January, every correlation rose (same date rho %s).", range_txt(loo_v("Concurrent", "rho_without_january"))),
    sprintf("Comparing transcripts with SIR at the next sampling date gave higher correlations in all four classes (rho %s with all dates; %s without January).",
            range_txt(loo_v("One-month lag", "rho_all_dates")), range_txt(loo_v("One-month lag", "rho_without_january"))),
    sprintf("Whole profiles were moderately aligned (Procrustes r = %s, exploratory %s).", f2(protest$procrustes_correlation), pv(protest$protest_p))
  )),
  p(b("What this means. "), "Outside midwinter, substrate-use potential and expression of the matching pathways broadly rose and fell together, and they lined up better when transcripts were compared with SIR measured about a month later. This fits the different timescales in P2: expression changes first, capacity follows. In midwinter the two separated: the potential to respond stayed high while expression dropped. The two measurements therefore capture different, partly sequential, aspects of microbial function."),
  p(b("What this cannot tell us. "), "Seven dates are too few for formal tests, so the correlations are descriptive. The lag comparison uses unequal intervals, and the January to February pair is missing one plot. The pattern is a hypothesis to test with denser sampling, not an established lag."),

  h3("Supplementary figures"),
  p(b("Figure S1. Root biomass and field soil respiration. "),
    sprintf("Season-specific slopes from a mixed model with a date random intercept. All observations: plant-active slope %s (%s, %s); plant-dormant %s (%s). Excluding the three observations flagged by both a whole-dataset and a within-plot outlier screen (SM1621, SM1078, SM1080): plant-active %s (%s, %s). These points were kept in the primary analysis because no measurement error was found; they may be genuine hot spots or hot moments.",
            f2(s1v(TRUE, "Plant-active")$slope), ci(s1v(TRUE, "Plant-active")$lower, s1v(TRUE, "Plant-active")$upper), pv(s1v(TRUE, "Plant-active")$p),
            f2(s1v(TRUE, "Plant-dormant")$slope), pv(s1v(TRUE, "Plant-dormant")$p),
            f2(s1v(FALSE, "Plant-active")$slope), ci(s1v(FALSE, "Plant-active")$lower, s1v(FALSE, "Plant-active")$upper), pv(s1v(FALSE, "Plant-active")$p))),
  img(fig("FigureS1", "FigureS1_root_respiration.png"), 6),
  p(b("Figure S2. Each substrate through the year. "), sprintf("Smooth curves (GAMs) through each substrate's log2 SIR, grouped by chemical class and ordered by how much of the variation the curve explains. Substrates with no detectable change through the year: %s.", paste(lc(steady), collapse = ", "))),
  img(fig("Figure2", "FigureS2_substrate_GAM_trajectories.png"), 6.5),
  p(b("Figure S2b. Does the H1b answer depend on analysis choices? "),
    sprintf("The complex minus readily available difference under the stricter root-exudate rule, z-scores, without the transition dates, with 24-h readings, and per gram of soil. All estimates were small (%s log2) and none was significant.", range_txt(h1b_sens$estimate))),
  img(fig("Figure2", "FigureS2b_H1b_sensitivity.png"), 5.5),
  p(b("Figure S2c. Complex and readily available substrates through the year. "), "Both groups followed the same seasonal course, which is why the mix did not change."),
  img(fig("Figure2", "FigureS2c_availability_trajectories.png"), 5),
  p(b("Figure S3. Does the H2 answer depend on analysis choices? "), "Season R² for the primary analysis, metabolism-only categories, values without within-sample normalization, and without the October 2023 transition date."),
  img(fig("Figure3", "FigureS3_H2_sensitivity.png"), 5.5),
  p(b("Figure S4. Correspondence details. "), "(a) Class SIR against class transcript z-scores, same date and next date. (b) Procrustes residuals by date (higher = poorer agreement between the SIR and transcript pictures). (c) All correlations, at date and plot level and with an alternative way of summarising transcripts."),
  img(fig("Figure4", "FigureS4a_class_scatter.png"), 6.5),
  img(fig("Figure4", "FigureS4b_procrustes_by_date.png"), 4.5),
  img(fig("Figure4", "FigureS4c_all_correlations.png"), 5.5),

  h2("Take Home Messages"),
  ul(c(
    sprintf("Plant dormancy cut field soil respiration by about %d%% as soils cooled by %s °C, while microbial biomass, root biomass, and potential enzyme activities stayed similar.",
            round(100 * (1 - f1v("Resp_rate_mg_C_m2_h")$ratio)), f2(abs(f1v("SoilTemp_mean")$estimate), 0)),
    sprintf("At the same time, soil microbes responded about %s more strongly to every one of the 31 added substrates. Dormancy changed how strongly microbes responded, not which substrates they favoured: complex substrates did not gain clearly more than readily available ones.", pct(level$estimate) %>% str_remove("^\\+")),
    "Field expression of microbial functions differed on the two dormant dates, and expression of substrate-related pathways, especially amino acid metabolism, was lowest in January.",
    "Substrate-use potential and expression of the matching pathways broadly tracked each other except in midwinter, and tracked better when transcripts were compared with SIR a month later, suggesting that expression and potential describe different, possibly sequential, aspects of microbial function."
  )),

  h2("Table 1. MicroResp substrates and their groups"),
  p("H1b group follows one rule: complex = must be broken apart outside the cell before uptake (polymer or ester); readily available = small molecule taken up directly. Chemical classes are the same as those used for transcripts, extended to all substrates."),
  tbl(table1, rel = c(1.1, 1.5, 1.2, 3.0, 1.0)),

  h2("Table 2. Transcript pathways linked to each compound class"),
  p("Only pathways linked to substrates of a single class are used. Share of all expression is the average across the 35 samples."),
  tbl(table2, rel = c(0.9, 1.8, 2.6, 2.4, 0.9)),

  h2("Glossary"),
  tbl(glossary, rel = c(1.3, 4))
)

write_docx(storyboard, file.path(out_dir, "Storyboard_v4.docx"))

# -----------------------------------------------------------------------------
# Analysis notes
# -----------------------------------------------------------------------------

decisions <- tribble(
  ~`#`, ~Decision,
  "D1 (revised 2026-10-25)", "H1b groups by one rule: complex = polymers and esters that must be broken apart outside the cell (8); readily available = small molecules (23); no 'other'. Sensitivity: only root-exudate types count as readily available. Within-class polymer vs monomer comparisons added.",
  "D2", "Compound classes for Q2/Q3: lipid-like, aromatic, amino sugar, amino acid (including BSA and GABA); only BRITE C categories linked to a single class are used (20).",
  "D3", "Root biomass x respiration is a short SI analysis (Fig. S1).",
  "D4", "Top-3-by-GAM panels removed from the main figures; all 31 GAMs in Fig. S2.",
  "D5", "Transcripts as within-sample proportions over the 325 retained categories; values as provided as a sensitivity.",
  "D6", "SM1082 (July 2023, plot 57) kept with only SEOC set to NA.",
  "D7", "Seasons Apr-Oct / Nov-Mar; sensitivity without Oct 2023 and Apr 2024.",
  "D8", "Class aggregation: transcripts summed over a class's pathways; SIR as the mean of within-substrate z-scores.",
  "D9", "Fig. 3A uses all 325 categories; metabolism-only as sensitivity.",
  "D10", "Fig. 1 shows season shading and irrigation window; harvest/saturation markers pending.",
  "D11", "All root observations kept; sensitivity excludes the 3 points flagged by both screens.",
  "D12 (2026-10-25)", "H1a reported as two separate results, level and mix; NMDS drawn without connecting lines; Q3 includes a leave-one-date-out check."
)

review_answers <- c(
  b("Does SIR vary by season? ") %>% paste0(sprintf("Yes, in level: %s on dormant dates (%s), with all 31 substrates higher. Not in mix: the relative pattern across substrates did not change (%s). The whole-profile test (%s) mixes both, which is why it looked inconclusive. NMDS is a picture, not a test; dormant and active points overlap in it.", pct(level$estimate), pv(level$p), pv(mix$p_permanova), pv(whole$p_permanova))),
  b("Why plots cannot be independent replicates of season. ") %>% paste0(sprintf("Season comes with the date, and plots on a date share weather, irrigation timing, field trip, and assay run. %d%% of the variation in overall SIR level is shared by plots on the same date and only %d%% is due to plots. Counting plots as replicates (v3 approach) gives p = 0.001 for the whole profile on the same data; counting dates gives p = %s. Plots are still used in the mixed models, which account for shared dates.", share_of("Sampling date"), share_of("Plot"), f2(whole$p_permanova, 3))),
  b("Are Q2 and Q3 valid with two dormant dates? ") %>% paste0("H2 can be tested, but only in the sense 'the two dormant dates differ from the five active dates more than any other pair of dates'; p = 0.048 is the floor. It cannot separate plant dormancy from anything else about January and March 2024. The class part of Q2 and all of Q3 were already framed as descriptive, which suits seven dates. Q3 does not rely on the season split at all; it asks whether two measurements move together across dates, and the leave-one-date-out check shows January drives most of the mismatch. Suggested wording in results: 'on the two plant-dormant sampling dates' rather than 'in the dormant season'."),
  b("Why the H1b grouping changed. ") %>% paste0("v4 first used only substrates that clearly matched the H1b examples, leaving 14 as 'other'. The revised grouping uses one rule (must the substrate be broken apart outside the cell?), so every substrate is assigned. The transcript classes were not used for H1b directly because they mix polymers and monomers (chitin and NAG are both 'amino sugar'); instead the same classes are used to compare each polymer with its own building blocks."),
  b("Putting the branch online. ") %>% paste0("A pull request is not needed. Pushing the branch (git push -u origin manuscript-v4-analysis) puts it on GitHub separately from main; nothing is merged. The repository is private.")
)

changes_vs_v3 <- tribble(
  ~Result, ~`Storyboard v3`, ~`v4`, ~`Why it changed`,
  "SIR seasonality (H1a)", "PERMANOVA R² = 0.097, p = 0.002", sprintf("Level %s (%s); mix R² = %s (%s); whole profile R² = %s (%s)", pct(level$estimate), pv(level$p), f2(mix$R2), pv(mix$p_permanova), f2(whole$R2), pv(whole$p_permanova)), "Dates, not plots, as replicates; level and mix separated.",
  "Transcript seasonality (H2)", "R² = 0.234, p = 0.001", sprintf("R² = %s, %s (minimum possible)", f2(h2_res$R2), pv(h2_res$p_permanova)), "Date-level test; within-sample proportions.",
  "H1b", "Not tested (top 3 substrates shown)", sprintf("Complex vs readily available %s log2, %s; inconsistent within classes", f2(h1b$estimate), pv(h1b$p)), "Explicit availability rule and test.",
  "Linked pathways", "31 including TCA and glycolysis", "20 linked to a single class", "Pathways shared by many classes cannot distinguish them.",
  "SIR vs transcripts", "Top 3 substrates; Procrustes r = 0.335, p = 0.083", sprintf("Class level, same date and next date, leave-one-date-out; Procrustes r = %s", f2(protest$procrustes_correlation)), "Class-level design matching Q3; transcripts as proportions; MicroResp v3.",
  "Root vs respiration", "Main Figure 1B; 9 candidate outliers", "Supplement; 3 points flagged by two screens as sensitivity", "Not a storyboard question; SM1082 restored."
)

open_items <- c(
  b("BRITE C provenance (needs your input). ") %>% paste0("No code in the repository creates Data/brite_C_expression.csv (committed 2026-07-07). Values are non-integer and per-sample totals vary about three-fold, so they are not raw counts or TPM that sum to a constant. The upstream steps (read mapping, normalization, KO to BRITE roll-up) should be added under Scripts/upstream/ or cited in Methods."),
  b("Raw data are not on GitHub. ") %>% paste0("DataHub March/ is in .gitignore; Data/README.md lists every source file with an MD5 checksum."),
  b("For S. Bell. ") %>% paste0("Reference for the season definition; harvest, fall-saturation, and irrigation dates for Fig. 1; rhizotron data; any field concentrations of assay substrates (to compare with the MicroResp dose)."),
  b("Data package README errors. ") %>% paste0("Soil Phenotypes README collection-date line (2024-06-07 to 2024-07-03); EEA enzyme counts (six/four stated, five/two listed); EEA Collection_date = #N/A for SM1502."),
  b("Microbial biomass C. ") %>% paste0("Two very low values (0.07 and 0.14 ug C/g) after LOD/2 replacement; worth checking against lab records."),
  b("AI-use statement. ") %>% paste0("Needed for ChatGPT (v3) and Claude (v4); check SBB policy.")
)

notes <- c(
  h1("Analysis notes for storyboard v4"),
  p("Companion to Storyboard_v4.docx: answers to review questions, decisions, data changes, how results moved from v3, open items, and how to rerun everything."),
  h2("1. Answers to the 2026-10-25 review questions"),
  ul(review_answers),
  h2("2. Items that need your attention"),
  ul(open_items),
  h2("3. How results changed from storyboard v3"),
  tbl(changes_vs_v3, rel = c(1.2, 1.5, 2.2, 2)),
  h2("4. Data changes"),
  ul(c(
    "MicroResp v3 replaces the original release, in which January 2024 IDs SM1503-SM1514 were shifted by one plot (SM1503 duplicated SM1502). Values were relabelled, not changed; the availability audit now stops if MicroResp and field labels disagree.",
    "SM1082 restored with only SEOC set to NA; this changed the LOD/2 replacement for Olsen organic P (not used in any figure).",
    "Transcript filtering drops 10 categories never detected, leaving 325."
  )),
  h2("5. Decision log"),
  tbl(decisions, rel = c(1, 5)),
  h2("6. Reproducing the analysis"),
  p("From the repository root with the DataHub March folder in place: Rscript Scripts/00_run_all.R (about 5 minutes), then Rscript Scripts/docs/write_storyboard_docx.R. Data/README.md lists every input with its origin and checksum; package versions are in Results/manuscript_v4/sessionInfo.txt."),
  h2("7. Statistical methods in brief"),
  ul(c(
    "Exact date-level permutation: every assignment of season labels to dates is listed (1,287 for SIR, 21 for transcripts); p is the share giving a PERMANOVA pseudo-F at least as large as observed. Applied to date means (primary) and to samples with labels permuted by whole date (check). Dispersion tested the same way.",
    "Level: sample mean of log2 SIR relative to each substrate's annual mean ~ Season + (1|Date) + (1|Plot).",
    "Mix: PERMANOVA on profiles divided by each sample's total.",
    "H1b: response ~ Season x group + (1|Date) + (1|Date:group) + (1|Sample) + (1|Substrate:Season); the same model within each chemical class; permutation check with 99,999 random reassignments of the complex label.",
    "Q3: Spearman correlations of date means (n = 7) and plot-dates; leave-one-date-out; no p-values."
  ))
)

write_docx(notes, file.path(out_dir, "Analysis_notes_v4.docx"))
cat("Wrote", file.path(out_dir, c("Storyboard_v4.docx", "Analysis_notes_v4.docx")), sep = "\n")
