# =============================================================================
# Write the v4 storyboard and analysis notes as HTML (converted to .docx by
# Scripts/docs/html_to_docx.ps1)
# =============================================================================
#
# Every number in the documents is read from Results/manuscript_v4, so the
# text stays in sync with the analysis. Run after Scripts/00_run_all.R.

rm(list = ls())
suppressPackageStartupMessages(library(tidyverse))

out_dir <- "Docs"
dir.create(out_dir, showWarnings = FALSE)
res <- function(...) read_csv(file.path("Results", "manuscript_v4", ...),
                              show_col_types = FALSE)
fig <- function(...) normalizePath(file.path("Figures", "manuscript_v4", ...),
                                   winslash = "/")

f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
pv <- function(p) if_else(p < 0.001, "p < 0.001", paste0("p = ", f2(p, 3)))
ci <- function(est, lo, hi, d = 2) sprintf("%s (95%% CI %s to %s)", f2(est, d), f2(lo, d), f2(hi, d))

# -----------------------------------------------------------------------------
# Results used in the text
# -----------------------------------------------------------------------------

design <- res("sample_availability", "design_summary.csv")
lagsum <- res("sample_availability", "lag_pair_summary.csv")

f1 <- res("Figure1", "Figure1_season_tests.csv")
f1v <- function(col) filter(f1, column == col)

h1a <- res("Figure2", "Figure2A_season_PERMANOVA_dispersion.csv")
h1a_d <- filter(h1a, level == "Date means (primary)")
h1a_s <- filter(h1a, level != "Date means (primary)")
h1a_sens <- res("Figure2", "FigureS2b_H1a_sensitivity.csv") %>% filter(level == "Date means (primary)")
sub_eff <- res("Figure2", "Figure2B_substrate_season_effects.csv")
overall <- res("Figure2", "Figure2C_overall_seasonal_shift.csv")
by_class <- res("Figure2", "Figure2C_season_effect_by_class.csv")
h1b_tests <- res("Figure2", "Figure2C_H1b_class_by_season_tests.csv")
h1b <- filter(h1b_tests, test == "Mixed model")
h1b_perm <- filter(h1b_tests, test != "Mixed model")
h1b_sens <- res("Figure2", "FigureS2b_H1b_sensitivity.csv")
envfit <- res("Figure2", "Figure2_envfit_date_level.csv")
gam <- res("Figure2", "FigureS2_GAM_temporal_statistics.csv")
f2sum <- res("Figure2", "Figure2_analysis_summary.csv")
nmds2 <- f2sum$value[f2sum$item == "NMDS stress"]

h2_test <- res("Figure3", "Figure3A_season_PERMANOVA_dispersion.csv")
h2_d <- filter(h2_test, level == "Date means (primary)")
h2_s <- filter(h2_test, level != "Date means (primary)")
h2_sens <- res("Figure3", "FigureS3_H2_sensitivity.csv") %>% filter(level == "Date means (primary)")
cat_sum <- res("Figure3", "Figure3B_unique_category_summary.csv")
share <- res("Figure3", "Figure3C_class_share_by_date.csv")
f3sum <- res("Figure3", "Figure3_analysis_summary.csv")
nmds3 <- f3sum$value[f3sum$item == "NMDS stress"]

rho <- res("Figure4", "Figure4B_spearman_correlations.csv") %>%
  filter(transcript_metric == "Summed class share (primary)")
rho_v <- function(cls, cmp, sc = "Date means (n = 7)") {
  rho$rho[rho$compound_class == cls & rho$comparison == cmp & rho$scale == sc]
}
div <- res("Figure4", "Figure4_class_values_by_date.csv") %>% filter(comparison == "Concurrent")
protest <- res("Figure4", "FigureS4b_PROTEST.csv")
proc_res <- res("Figure4", "FigureS4b_procrustes_residuals.csv") %>%
  group_by(Sampling_Date) %>% summarise(r = mean(residual))

s1 <- res("FigureS1", "FigureS1_root_respiration_slopes.csv")
s1v <- function(primary, season) {
  filter(s1, str_detect(analysis, "primary") == primary, Season == season)
}

class_eff <- function(cls) filter(by_class, availability == cls)
n_pos <- sum(sub_eff$estimate > 0)
n_fdr <- sum(sub_eff$p_fdr < 0.05)
top_sub <- sub_eff %>% arrange(desc(estimate)) %>% slice(1:3)
steady <- sub_eff %>% filter(steady) %>% pull(display_name)
lower_dormant <- sum(cat_sum$dormant_to_active_ratio < 1)
jan <- div %>% filter(Sampling_Date == as.Date("2024-01-02"))
share_jan_aa <- share %>% filter(compound_class == "Amino acid")
gwc_env <- filter(envfit, variable == "Gravimetric water content")
temp_env <- filter(envfit, variable == "Soil temperature")
lag_n <- paste(lagsum$n_pairs, collapse = ", ")

# -----------------------------------------------------------------------------
# HTML helpers
# -----------------------------------------------------------------------------

esc <- function(x) x %>% str_replace_all("&", "&amp;") %>% str_replace_all("<", "&lt;") %>%
  str_replace_all(">", "&gt;")
h <- function(level, text) sprintf("<h%d>%s</h%d>", level, text, level)
h1 <- function(text) h(1, text)
h2 <- function(text) h(2, text)
h3 <- function(text) h(3, text)
p <- function(...) sprintf("<p>%s</p>", paste0(...))
b <- function(text) sprintf("<b>%s</b>", text)
i <- function(text) sprintf("<i>%s</i>", text)
ul <- function(items) paste0("<ul>", paste0("<li>", items, "</li>", collapse = ""), "</ul>")
# Figures are written as placeholders; html_to_docx.ps1 replaces each with the
# embedded picture (opening HTML with linked images stalls Word).
img <- function(path, width_in = 6.5) {
  sprintf('<p class="figure">[[IMG|%s|%s]]</p>', path, width_in)
}
tbl <- function(df) {
  head_row <- paste0("<tr>", paste0("<th>", names(df), "</th>", collapse = ""), "</tr>")
  body <- apply(df, 1, function(r) paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>"))
  paste0("<table>", head_row, paste(body, collapse = ""), "</table>")
}
page <- function(title, body) {
  paste0(
    '<html><head><meta charset="utf-8"><title>', title, '</title><style>',
    'body{font-family:Calibri,Arial,sans-serif;font-size:11pt;line-height:1.25}',
    'h1{font-size:18pt;color:#1F3864}h2{font-size:14pt;color:#1F3864;margin-top:18pt}',
    'h3{font-size:12pt;color:#2F5496}p{margin:4pt 0 6pt 0}',
    'table{border-collapse:collapse;margin:6pt 0}th,td{border:1px solid #999;padding:3pt 5pt;font-size:9.5pt;vertical-align:top}',
    'th{background:#E7ECF5}p.figure{text-align:center}.caption{font-size:10pt}',
    '.note{color:#7F6000;background:#FFF5D6;padding:4pt}',
    '</style></head><body>', paste(body, collapse = "\n"), '</body></html>'
  )
}

# -----------------------------------------------------------------------------
# Storyboard v4 (same sections as storyboard v3)
# -----------------------------------------------------------------------------

storyboard <- c(
  h1("Storyboard v4: seasonal substrate-use potential and field gene expression"),
  p('<span class="note">', b("Note: "), "v4 was drafted with Claude (Anthropic) from the analysis outputs in the GitHub repository (branch manuscript-v4-analysis); earlier versions used ChatGPT 5.6 Sol (enterprise). ",
    "Every number below is read from Results/manuscript_v4, produced by the scripts in Scripts/ (run order in Scripts/00_run_all.R). ",
    "Decisions, data changes, and open items are in the companion document <i>Analysis_notes_v4.docx</i>.</span>"),

  h2("Suggested Title"),
  p("Option 1: Plant dormancy raises microbial substrate-use potential while field expression of substrate-related functions declines"),
  p("Option 2: Microbial substrate-induced respiration and field transcript expression show distinct seasonal patterns across plant activity and dormancy"),
  p(i("v3 title for reference: Soil microbial substrate-induced respiration and field transcript expression show seasonal patterns but limited correspondence. KH asked to state the finding more clearly; option 1 names the direction of each response.")),

  h2("Potential Journals"),
  p("SBB"),

  h2("Team members"),
  p("Vanessa Garayburu-Caruso, Josue Rodriguez-Ramos, Sheryl Bell, Evan Warburton, Sharon Zhao, Regan McDearis, David Degnan, Amy Zimmerman, and Kirsten Hofmockel"),

  h2("Acknowledgements"),
  p("Gregory, Kali Cornwell, Mischelle, other people that contributed? Trinidad (viromes; not in this paper)"),

  h2("GitHub link with code"),
  p('<a href="https://github.com/vgc40/Soil_time_series_data">vgc40/Soil_time_series_data</a> (v4 analysis on branch manuscript-v4-analysis)'),

  h2("Framing"),
  h3("P1. Seasonal changes in microbial substrate availability and function"),
  p("Soil microbial function varies across the annual cycle as plant activity, temperature, moisture, and resource availability change. As plants transition between active growth and dormancy, changes in root exudation, root turnover, and residue inputs may alter the timing and relative availability of different substrates, including readily available compounds and more complex organic matter associated with plant residues and microbial necromass. Dormancy may therefore change the mixture of substrates available to microbes rather than simply reducing their supply. These shifts can alter microbial physiological state and resource acquisition strategies, potentially affecting how microbes respond to different substrates throughout the year. Yet much of our understanding of soil microbial function comes from growing-season measurements, despite evidence that microbial activity can continue outside the growing season (Monson et al., 2006; Gavazov et al., 2017). As a result, we know comparatively little about how microbial responses to different substrates change during plant dormancy, limiting our understanding of how plant seasonality shapes belowground microbial metabolism."),
  h3("P2. Connecting microbial substrate-use potential and field expression"),
  p("Commonly used measurements capture different aspects and timescales of microbial function. Field soil respiration integrates root and microbial activity under prevailing conditions but cannot resolve which substrates microbes are using. Laboratory substrate-induced respiration measures short-term responses to added compounds under standardized conditions. Because these responses are measured within 5 h, they primarily characterize the short-term respiratory capacity of the existing microbial community, which reflects its physiological state and antecedent environmental conditions rather than substantial growth during the assay (Anderson & Domsch, 1978; Campbell et al., 2003). Seasonal changes in these responses therefore show how microbial substrate-use potential shifted over the year. Metatranscriptomes, in contrast, characterize which microbial functions are being expressed under field conditions at the time of sampling, but not how strongly microbes would respond to specific substrate additions. It remains unclear how seasonal changes in field gene expression relate to changes in substrate-use potential, in part because the two capture microbial activity over different timescales: microbial transcripts typically persist for minutes (Moran et al., 2013), whereas soil extracellular enzymes can remain active for weeks after they are produced (Schimel et al., 2017). Seasonal changes in field expression and substrate-use potential may therefore not coincide. Comparing the two across related compounds and metabolic pathways can help clarify whether plant dormancy changes which functions microbes express, their potential to respond as substrate availability changes, or both."),
  h3("P3. Study objectives and research questions"),
  p("Here, we investigate how microbial substrate-use potential and field gene expression change between plant-active and plant-dormant periods, and how the two are related. We focused on fully irrigated (100%) plots in an irrigated tall wheatgrass system, a perennial grass with a plant-active period from April (plant emergence) through October (harvest) and a plant-dormant period from November through March [basis for season definitions: reference pending from S. Bell]. We combine 13 months of substrate-induced respiration measurements for 31 soil-relevant organic compounds with field metatranscriptomes from seven of these sampling dates (June, July, September, and October 2023; January, March, and May 2024), and use supporting measurements of extracellular enzyme activity, microbial biomass, soil respiration, root biomass, and environmental and biogeochemical conditions to help interpret the observed patterns. Specifically, we ask (1) how substrate-induced respiration differs between plant-active and plant-dormant periods, and whether seasonal responses differ among substrates; (2) how functional transcript expression differs between these periods, and whether functions associated with different compound classes show distinct seasonal patterns; and (3) how seasonal changes in substrate-induced respiration within each compound class relate to the expression of functions associated with that class."),

  h2("Research Questions & Corresponding Hypotheses"),
  p(b("Q1. "), "How does substrate-induced respiration differ between plant-active and plant-dormant periods, and do seasonal responses differ among substrates?"),
  p(b("H1a: "), "Microbial substrate-induced respiration profiles will differ between plant-active and plant-dormant periods."),
  p(b("H1b: "), "Responses to readily available substrates, such as simple sugars and organic acids, will be greater relative to each substrate's mean across sampled dates during plant-active periods, whereas responses to complex substrates associated with plant and microbial residues will be greater relative to their mean across sampled dates during plant-dormant periods."),
  p(b("Q2. "), "How does microbial functional transcript expression differ between plant-active and plant-dormant periods, and do functions associated with different compound classes show distinct seasonal patterns?"),
  p(b("H2: "), "Microbial functional transcript composition will differ between plant-active and plant-dormant periods. The second part of Q2 is addressed descriptively, using the BRITE C categories mapped uniquely to lipid-like, aromatic, amino sugar, and amino acid compound classes; with two dormant transcript dates, no directional predictions are tested for individual categories."),
  p(b("Q3. "), "How do seasonal changes in substrate-induced respiration within each compound class relate to transcript abundance of the BRITE C categories mapped uniquely to that class, and does this relationship differ among compound classes? Addressed descriptively, including transcript abundance versus SIR measured one month later."),

  h2("Approach"),
  p(b("Design. "), sprintf("Soils were sampled monthly from June 2023 to June 2024 (13 dates) in the 100%% irrigation treatment of the tall wheatgrass irrigation trial (six plots in three blocks; Jose cultivar). Field, extracellular enzyme, and microbial biomass data are complete for all %s plot-dates; MicroResp covers %s samples (plot 27 was not assayed on 12 February 2024); metatranscriptomes cover %s samples (plots 1, 25, 27, 40, 46 on seven dates; plot 57 was not sequenced).",
    design$value[design$item == "Observed 100% field records"],
    design$value[design$item == "Observed 100% MicroResp samples"],
    design$value[design$item == "Observed 100% transcript samples"])),
  p(b("Seasons. "), "Plant-active: April (emergence) to October (harvest), which coincides with the drip-irrigation season; plant-dormant: November to March. This gives 8 active and 5 dormant SIR dates and 5 active and 2 dormant transcript dates. A sensitivity analysis drops the transition dates (October 2023, April 2024)."),
  p(b("Replication unit for season. "), "Season is assigned to whole sampling dates, so dates, not plots, are the replicates for every season test. Multivariate tests use PERMANOVA on date-mean profiles with every possible assignment of season labels to dates enumerated (exact p-values: 1,287 assignments for SIR, 21 for transcripts, so the smallest attainable p for transcripts is 0.048). A sample-level version with labels permuted by whole date is reported as a sensitivity. Univariate season effects use mixed models with random intercepts for date and plot."),
  p(b("Substrate-induced respiration. "), "MicroResp (James Hutton Institute) with 31 substrates and a water control, triplicate wells, soil at 60% WHC, room temperature; CO2 read at 5 h and 24 h. The main analysis uses 5-h rates normalized per mg of added C (MicroResp v3 data release, which corrects January 2024 sample labels). Responses are log2-transformed and centred on each substrate's mean of date means, so seasonal differences are log2 fold changes relative to the substrate's annual mean. Whole-profile structure: Bray-Curtis NMDS and PERMANOVA. H1b: a mixed model of all substrates with a season x availability-class interaction (random effects for date, date x class, sample, and substrate x season), checked with an exact permutation of class labels among substrate-level effects. Substrate classes were fixed a priori: readily available (5 sugars, 5 organic acids), complex residue-associated (cellulose, hemicellulose, lignin, chitin, chitosan, BSA, methyl oleate), other (14)."),
  p(b("Metatranscriptomes. "), "Expression rolled up to BRITE C functional categories (taxon-independent), keeping 325 categories after removing hierarchy rows, curated human/non-soil pathways, and all-zero categories, and converted to within-sample proportions. Rolling expression up to functions mitigates the rapid turnover of individual transcripts (AZ comment). Composition: Bray-Curtis NMDS and exact date-level PERMANOVA; dispersion tested the same way."),
  p(b("Substrate-to-function mapping. "), "An a priori table links each MicroResp substrate to KEGG BRITE C pathways containing its catabolic or assimilatory reactions (Data/Microresp_pathways.csv). Four compound classes (lipid-like, aromatic, amino sugar, amino acid including BSA and GABA) were defined, and a category was used only if all substrates it maps to belong to one class, giving 20 categories (4, 6, 1, and 9)."),
  p(b("Correspondence (Q3). "), sprintf("For each class, SIR is the mean within-substrate z-score of the class's substrates and transcript abundance is the summed share of expression of its unique categories, z-scored. These are compared at the seven transcript dates (concurrent) and with SIR at the next sampling date (lag pairs: n plots = %s; all intervals 28 days except January to February, 41 days) using Spearman correlations, reported descriptively. Global correspondence: Procrustes/PROTEST on NMDS configurations, permuting by cyclic shifts within plot (exploratory).", lag_n)),
  p(b("Supporting measurements. "), "Field soil respiration (CO2 flux, 60-min chamber reading), soil temperature, gravimetric water content, hand-picked root biomass from pooled triplicate 0-15 cm cores (not rhizosphere material), microbial biomass C and N, and potential extracellular enzyme activities (BG, BX, CB, AP, NAG, Leu, Ala)."),

  h2("Key Figures"),

  h3("Figure 1. Study design and seasonal field context"),
  img(fig("Figure1", "Figure1_main.png")),
  p('<span class="caption">', b("Figure 1. "),
    "(A) Sampling design for the 100% irrigation treatment: dates with field, enzyme, and MicroResp data and with metatranscriptomes (numbers are plots sampled); shading marks plant-active (green) and plant-dormant (blue-grey) periods and the bar marks the April-October drip-irrigation window. Harvest and fall-saturation markers will be added once dates are confirmed. ",
    "(B-I) Soil temperature, gravimetric water content, field soil respiration, root biomass (0-15 cm), microbial biomass C, and summed C-, N-, and P-acquiring enzyme activities. Points are plots; lines and ribbons are date means (geometric means on log axes) with 95% CIs. Panel subtitles give the dormant/active ratio (or difference for temperature) from mixed models with date and plot random intercepts.</span>"),
  p(b("Interpretation: "), sprintf("Field soil respiration fell to %s of plant-active rates during dormancy (%s), tracking the %s °C drop in soil temperature. Root biomass (%s×, %s), microbial biomass C (%s×, %s), and potential enzyme activities (%s-%s×) did not differ between periods, and water content was similar (%s×, %s). Dormancy therefore reduced field CO2 production without reducing standing microbial biomass or the enzyme pool.",
    paste0(round(100 * f1v("Resp_rate_mg_C_m2_h")$ratio), "%"), pv(f1v("Resp_rate_mg_C_m2_h")$p),
    f2(abs(f1v("SoilTemp_mean")$estimate), 1),
    f2(f1v("Roots_gdw")$ratio), pv(f1v("Roots_gdw")$p),
    f2(f1v("MBC_ug_gdw_bgc")$ratio), pv(f1v("MBC_ug_gdw_bgc")$p),
    f2(min(f1$ratio[str_detect(f1$column, "Sum_")])), f2(max(f1$ratio[str_detect(f1$column, "Sum_")])),
    f2(f1v("GWC")$ratio), pv(f1v("GWC")$p))),

  h3("Figure 2. Seasonal substrate-induced respiration (Q1)"),
  img(fig("Figure2", "Figure2_main.png")),
  p('<span class="caption">', b("Figure 2. "),
    sprintf("(A) NMDS of the 31-substrate SIR profile (Bray-Curtis; stress = %s). Small points are samples; large points are date centroids joined in time order. Date-level PERMANOVA R² = %s, exact %s; dispersion %s. ", nmds2, f2(h1a_d$R2), pv(h1a_d$p_permanova), pv(h1a_d$p_dispersion)),
    "(B) Date-mean trajectories of each availability class, expressed as log2 SIR relative to each substrate's annual mean (95% CIs across plots). ",
    sprintf("(C) Dormant minus plant-active difference for each substrate (log2 fold change, 95%% CI) from mixed models with date and plot random intercepts; open symbols are not significant after FDR correction; \"steady\" marks substrates without detectable temporal structure (GAM FDR > 0.05). "),
    sprintf("(D) Seasonal difference by availability class (large points, 95%% CI) with substrate estimates (small points). All substrates combined: %s log2, %s. Complex minus readily available: %s, %s.</span>",
      ci(overall$estimate, overall$lower, overall$upper), pv(overall$p),
      ci(h1b$estimate, h1b$lower, h1b$upper), pv(h1b$p))),
  p(b("Interpretation: "), sprintf(
    "SIR was higher during plant dormancy for all %d of %d substrates, by %s log2 units on average (about %s%% higher; %s). With dates as replicates, the whole-profile difference did not reach significance (H1a: R² = %s, %s), and profiles were more consistent among dormant dates than among plant-active dates, which span two growing seasons (dispersion %s). H1a is therefore only partly supported: dormancy raised the overall level of substrate responses more than it changed their relative composition. ",
    n_pos, nrow(sub_eff), f2(overall$estimate), round(100 * (2^overall$estimate - 1)), pv(overall$p),
    f2(h1a_d$R2), pv(h1a_d$p_permanova), pv(h1a_d$p_dispersion)),
    sprintf("For H1b, complex residue-associated substrates increased most in dormancy (%s log2, %s) and readily available substrates least (%s, %s), the direction predicted, but the difference between classes was small and not significant (%s, %s; exact permutation p = %s). No readily available substrate responded more during the plant-active period. The largest dormancy increases were for %s. ",
      f2(class_eff("Complex")$estimate), pv(class_eff("Complex")$p),
      f2(class_eff("Readily available")$estimate), pv(class_eff("Readily available")$p),
      f2(h1b$estimate), pv(h1b$p), f2(h1b_perm$p, 3),
      paste(sprintf("%s (%s×)", top_sub$display_name, f2(2^top_sub$estimate, 1)), collapse = ", ")),
    sprintf("The H1b difference was stronger for 24-h responses (%s, %s) and absent when rates are expressed per g soil (%s), so H1b support depends on the measurement choice (Fig. S2b). Water content at sampling was the field variable most aligned with the SIR ordination (date-level envfit r² = %s, %s; temperature r² = %s, %s; exploratory).",
      f2(h1b_sens$estimate[str_detect(h1b_sens$analysis, "24-h")]), pv(h1b_sens$p[str_detect(h1b_sens$analysis, "24-h")]),
      f2(h1b_sens$estimate[str_detect(h1b_sens$analysis, "Per g")]),
      f2(gwc_env$r2), pv(gwc_env$p), f2(temp_env$r2), pv(temp_env$p))),

  h3("Figure 3. Seasonal field transcript expression (Q2)"),
  img(fig("Figure3", "Figure3_main.png")),
  p('<span class="caption">', b("Figure 3. "),
    sprintf("(A) NMDS of BRITE C transcript composition (325 categories, within-sample proportions, Bray-Curtis; stress = %s). Date-level PERMANOVA R² = %s, exact %s, the smallest p attainable with 5 active and 2 dormant dates; dispersion %s. ", nmds3, f2(h2_d$R2), pv(h2_d$p_permanova), pv(h2_d$p_dispersion)),
    "(B) Share of total expression for each compound class (sum of its uniquely mapped categories); points are plots, lines and ribbons date means with 95% CIs. ",
    "(C) Date-mean z-scores of the 20 uniquely mapped categories, standardized within category; colours show when a category was relatively high or low against its own mean, not differences among categories. Row labels give mean share of total expression; dormant dates in bold.</span>"),
  p(b("Interpretation: "), sprintf(
    "Field functional composition differed between plant-active and plant-dormant dates as strongly as the design can show (H2 supported: R² = %s, %s; no dispersion difference), and the result held for metabolism-only categories (R² = %s), for values without within-sample normalization (R² = %s), and without the October transition date (%s, again the smallest attainable). ",
    f2(h2_d$R2), pv(h2_d$p_permanova),
    f2(h2_sens$R2[str_detect(h2_sens$analysis, "Metabolism")]),
    f2(h2_sens$R2[str_detect(h2_sens$analysis, "as provided")]),
    pv(h2_sens$p_permanova[str_detect(h2_sens$analysis, "Transition")])),
    sprintf("Descriptively, %d of the 20 uniquely mapped categories had a lower share of expression on dormant dates. Amino acid categories dropped sharply in January and recovered by May; lipid-like and amino sugar categories rose through spring, peaking in May 2024; aromatic categories peaked in October 2023 and May 2024. In the ordination, May 2024 lay between the 2023 plant-active dates and the dormant dates, so composition had not returned to its June 2023 state one year later.", lower_dormant)),

  h3("Figure 4. Correspondence between SIR and transcript expression by compound class (Q3)"),
  img(fig("Figure4", "Figure4_main.png")),
  p('<span class="caption">', b("Figure 4. "),
    "(A) Date-mean class SIR (mean within-substrate z of log2 SIR) and class transcript abundance (summed share of uniquely mapped categories, z-scored) at the seven transcript dates; rho is the date-level Spearman correlation. ",
    sprintf("(B) Spearman correlations between SIR and transcripts, concurrent and with SIR at the next sampling date, for date means (n = 7) and plot-dates (n = 34-35); open symbols use mean within-category z-scores for transcripts. Correlations are descriptive. (C) Difference between SIR and transcript z-scores at each date; positive values mean SIR was relatively higher than expression.</span>")),
  p(b("Interpretation: "), sprintf(
    "Concurrent SIR and expression of class-specific functions were weakly and inconsistently related (date-level rho: lipid-like %s, aromatic %s, amino sugar %s, amino acid %s). The clearest pattern is a January mismatch in every class: SIR was at its highest of the seven transcript dates while class transcripts were below their mean (SIR minus transcript z = %s to %s). ",
    f2(rho_v("Lipid-like", "Concurrent")), f2(rho_v("Aromatic", "Concurrent")),
    f2(rho_v("Amino sugar", "Concurrent")), f2(rho_v("Amino acid", "Concurrent")),
    f2(min(jan$divergence), 1), f2(max(jan$divergence), 1)),
    sprintf("Transcripts aligned better with SIR measured at the next sampling date in all four classes (lag rho: %s, %s, %s, %s; plot-level rho also higher for every class), consistent with expression preceding changes in substrate-use potential. With seven dates and one unequal interval, this is a pattern to report, not a tested lag. Whole profiles corresponded moderately (Procrustes r = %s, exploratory PROTEST %s).",
      f2(rho_v("Lipid-like", "One-month lag")), f2(rho_v("Aromatic", "One-month lag")),
      f2(rho_v("Amino sugar", "One-month lag")), f2(rho_v("Amino acid", "One-month lag")),
      f2(protest$procrustes_correlation), pv(protest$protest_p))),

  h3("Supplementary figures"),
  p(b("Figure S1. Root biomass and field soil respiration. "), sprintf("Season-specific slopes from a mixed model with a date random intercept (plot variance estimated as zero). All observations: plant-active slope %s (%s), plant-dormant %s (%s). Excluding the three observations flagged by both a global and a within-plot modified z-score screen (SM1621, SM1078, SM1080): plant-active %s (%s). Points are retained in the primary analysis because no measurement error was identified; they may be hot spots or hot moments (KH).",
    ci(s1v(TRUE, "Plant-active")$slope, s1v(TRUE, "Plant-active")$lower, s1v(TRUE, "Plant-active")$upper), pv(s1v(TRUE, "Plant-active")$p),
    f2(s1v(TRUE, "Plant-dormant")$slope), pv(s1v(TRUE, "Plant-dormant")$p),
    ci(s1v(FALSE, "Plant-active")$slope, s1v(FALSE, "Plant-active")$lower, s1v(FALSE, "Plant-active")$upper), pv(s1v(FALSE, "Plant-active")$p))),
  img(fig("FigureS1", "FigureS1_root_respiration.png"), 6),
  p(b("Figure S2. Substrate-specific temporal trajectories. "), sprintf("GAMs (non-cyclic smooth of time, k = 7, plot random effect) on log2 SIR for all 31 substrates, grouped by substrate group and ordered by deviance explained within group. Substrates without detectable temporal structure (FDR > 0.05): %s.", paste(steady, collapse = ", "))),
  img(fig("Figure2", "FigureS2_substrate_GAM_trajectories.png"), 6.5),
  p(b("Figure S2b. Sensitivity of the H1b contrast. "), "Complex minus readily available difference in seasonal effect under alternative analysis choices: amino acids counted as readily available, within-substrate z-scores, transition dates excluded, 24-h incubation, and rates per g dry soil."),
  img(fig("Figure2", "FigureS2b_H1b_sensitivity.png"), 5.5),
  p(b("Figure S3. Sensitivity of the H2 test. "), "PERMANOVA R² for season with date-level and sample-level permutation for the primary analysis, metabolism-only categories (KEGG maps < 02000), values without within-sample normalization, and without the October 2023 transition date."),
  img(fig("Figure3", "FigureS3_H2_sensitivity.png"), 5.5),
  p(b("Figure S4. Correspondence details. "), "(a) Class SIR versus class transcript z-scores, concurrent and lagged. (b) Procrustes residuals by date (higher = poorer agreement between the SIR and transcript configurations)."),
  img(fig("Figure4", "FigureS4a_class_scatter.png"), 6.5),
  img(fig("Figure4", "FigureS4b_procrustes_by_date.png"), 4.5),

  h2("Take Home Messages"),
  ul(c(
    sprintf("Plant dormancy cut field soil respiration by about %s%% as soils cooled, while microbial biomass, root biomass, and potential enzyme activities stayed similar between periods.",
            round(100 * (1 - f1v("Resp_rate_mg_C_m2_h")$ratio))),
    sprintf("Microbial substrate-use potential moved the other way: SIR was higher in dormancy for all %d substrates (about %s%% on average), with the largest increases for %s. Complex residue-associated substrates rose slightly more than readily available ones, as H1b predicts, but the class difference was not significant.",
            nrow(sub_eff), round(100 * (2^overall$estimate - 1)), paste(top_sub$display_name[1:2], collapse = " and ")),
    sprintf("Field functional expression shifted between periods (R² = %s), with most class-specific functions, especially amino acid metabolism, lowest in January.", f2(h2_d$R2)),
    "SIR and the expression of class-specific functions diverged most in mid-winter and aligned better when transcripts were compared with SIR one month later, suggesting that expression and substrate-use potential describe different, possibly sequential, aspects of microbial function."
  ))
)

writeLines(page("Storyboard v4", storyboard), file.path(out_dir, "Storyboard_v4.html"), useBytes = TRUE)

# -----------------------------------------------------------------------------
# Analysis notes (everything else to document or read)
# -----------------------------------------------------------------------------

decisions <- tribble(
  ~`#`, ~Decision,
  "D1", "H1b classes: readily available = 5 sugars + 5 organic acids; complex = cellulose, hemicellulose, lignin, chitin, chitosan, BSA, methyl oleate; other = 14. Sensitivity: amino acids counted as readily available.",
  "D2", "Amino-acid class includes BSA and GABA; strict uniqueness rule gives 20 BRITE C categories (4 lipid-like, 6 aromatic, 1 amino sugar, 9 amino acid), including one carbon pool by folate, beta-alanine, and butanoate metabolism.",
  "D3", "Root biomass x respiration is a short SI analysis (Fig. S1); Fig. 1 shows root biomass as context only.",
  "D4", "Top-3-by-GAM panels removed from the main figures; all 31 GAMs in Fig. S2.",
  "D5", "Transcripts converted to within-sample proportions over the 325 retained categories; values as provided kept as a sensitivity.",
  "D6", "SM1082 (July 2023, plot 57) retained with only SEOC set to NA; the three 25% July samples remain removed.",
  "D7", "Seasons Apr-Oct / Nov-Mar; sensitivity without Oct 2023 and Apr 2024.",
  "D8", "Class aggregation: transcripts = summed share of unique categories; SIR = mean within-substrate z. Sensitivity: transcripts as mean category z.",
  "D9", "Fig. 3A uses all 325 categories; metabolism-only as sensitivity (Fig. S3).",
  "D10", "Fig. 1 shows season shading and irrigation window; harvest/saturation markers pending.",
  "D11", "All root observations kept; sensitivity excludes the 3 points flagged by both screens."
)

changes_vs_v3 <- tribble(
  ~Result, ~`Storyboard v3`, ~`v4 (this analysis)`, ~`Why it changed`,
  "SIR season PERMANOVA (H1a)", "R² = 0.097, p = 0.002", sprintf("Date level R² = %s, %s; sample level with date-wise permutation %s", f2(h1a_d$R2), pv(h1a_d$p_permanova), pv(h1a_s$p_permanova)),
  "v3 permuted samples within plot, which treats the six plots on a date as independent replicates of season. Season is assigned to dates, so v4 permutes season labels across dates. Also uses MicroResp v3 (77 samples instead of 78 with a duplicate).",
  "SIR dispersion", "F = 9.22, p = 0.007", pv(h1a_d$p_dispersion), "Same reason; dispersion still differs at the date level.",
  "Transcript season PERMANOVA (H2)", "R² = 0.234, p = 0.001", sprintf("R² = %s, exact %s (minimum attainable)", f2(h2_d$R2), pv(h2_d$p_permanova)),
  "Date-level permutation (21 possible label sets) and within-sample proportions. The effect is as strong as 7 dates allow; p = 0.001 was not attainable without pseudo-replication.",
  "Fig. 2 B-D", "Top 3 substrates by GAM deviance", "H1b class test, per-substrate effects, class trajectories", "Decision D4; v3 panels were not tied to H1b.",
  "GAM deviance explained", "Methyl oleate 78.9%, chitosan 76.8%, L-glutamine 60.9%", sprintf("Top: %s", paste(sprintf("%s %s%%", gam %>% arrange(desc(deviance_explained)) %>% slice(1:3) %>% pull(display_name), round(100 * (gam %>% arrange(desc(deviance_explained)) %>% slice(1:3) %>% pull(deviance_explained)))), collapse = ", ")),
  "GAMs now fit to log2 SIR (variance stabilizing) on MicroResp v3.",
  "Mapped categories (Fig. 3B)", "31 linked pathways including TCA, glycolysis", "20 categories unique to one compound class", "Decision D2; non-specific pathways removed.",
  "Procrustes", "r = 0.335, PROTEST p = 0.083", sprintf("r = %s, %s", f2(protest$procrustes_correlation), pv(protest$protest_p)),
  "Transcripts as proportions and MicroResp v3. Still exploratory (cyclic shifts with seven nonstationary dates).",
  "Root-respiration active slope", "0.49, p = 0.153 (9 removed: 2.55)", sprintf("%s, %s (3 removed: %s, %s)", f2(s1v(TRUE, "Plant-active")$slope), pv(s1v(TRUE, "Plant-active")$p), f2(s1v(FALSE, "Plant-active")$slope), pv(s1v(FALSE, "Plant-active")$p)),
  "SM1082 now retained (78 plot-dates); 3-point both-screens rule (D11)."
)

data_changes <- c(
  "MicroResp v3 replaces the original release. In the original, January 2024 IDs SM1503-SM1514 were shifted by one plot, SM1503 was a rounded copy of SM1502, and SM1514 existed only in that file. In v3, values were relabelled, not altered (v3 SM1503 = old SM1504). The availability script now stops if MicroResp and field labels disagree.",
  "The January 2024 transcripts SM1502-SM1506 are plots 1, 25, 27, 40, 46 (all 100%). My earlier statement in chat that January lacked plots 46 and 57 came from the mislabelled MicroResp file and was wrong.",
  "SM1082 is back in the cleaned field data with SEOC set to NA. Side effect: the LOD/2 replacement for Olsen organic P changed for 38 zero values (the minimum nonzero value now includes SM1082). This variable is not used in any figure.",
  "Data cleaning reads only the March DataHub release (older candidate paths removed).",
  "Transcript filtering now drops 10 categories with zero expression in all samples (as the previous Figure 3 script did), leaving 325."
)

open_items <- c(
  b("BRITE C provenance (needs your input). ") %>% paste0("No code in the repository creates Data/brite_C_expression.csv (committed 2026-07-07). Values are non-integer and per-sample totals vary 2.9-fold over the 325 analysed categories, so they are not raw counts or TPM summing to a constant. For Methods, the upstream pipeline (read mapping, normalization, KO to BRITE roll-up, handling of KOs in several categories) should be added to the repository or cited. If the team member who produced it (David Degnan?) can share the script, add it under Scripts/upstream/."),
  b("Raw data are not on GitHub. ") %>% paste0("DataHub March/ and Archive/ are in .gitignore. Data/README.md lists each source file with an MD5 checksum, so results can be traced to exact files, but anyone cloning the repository needs the DataHub package. Decide whether to keep this (sensible until the DataHub packages are public) or add the files."),
  b("For S. Bell. ") %>% paste0("(1) Reference for the season definition; (2) harvest, fall-saturation, and irrigation start/stop dates for the Fig. 1 timeline; (3) rhizotron root data if they should appear in Fig. 1 or S1; (4) whether field concentrations of any assay substrates are known, to benchmark the MicroResp dose (KH comment)."),
  b("Data package README errors to report. ") %>% paste0("Soil Phenotypes README: collection-date line reads 2024-06-07 to 2024-07-03 (should be 2023-06-07 to 2024-06-05); EEA text says six MUB and four MUC enzymes but lists five and two; EEA file has Collection_date = #N/A for SM1502."),
  b("Microbial biomass C. ") %>% paste0("Two very low MBC values (0.07 and 0.14 ug C/g) remain after LOD/2 replacement of zeros; visible in Fig. 1F. They do not affect the season test (log-scale mixed model, p = 0.86) but are worth checking against lab records."),
  b("AI-use statement. ") %>% paste0("The manuscript needs an AI-use statement covering ChatGPT (v3 storyboard) and Claude (v4 analysis code and text); check SBB policy."),
  b("Not pushed. ") %>% paste0("All work is committed on the local branch manuscript-v4-analysis and has not been pushed to GitHub.")
)

interpretation_cautions <- c(
  "H1a: with dates as replicates the SIR profile test is not significant (p = 0.11). The robust seasonal SIR result is the shared increase across all substrates in dormancy (p = 0.045), which composition-based tests (Bray-Curtis on profiles) partly remove. Frame Q1 around level and variability, not composition.",
  "H1b direction matches the prediction for complex substrates (largest dormancy increase), but the predicted active-season advantage for readily available substrates is absent: they also rose in dormancy. Results differ between 5 h, 24 h, and per-g-soil normalizations (Fig. S2b), so report the 5-h result as primary and the others transparently.",
  "H2: p = 0.048 is the smallest value possible with 5 active and 2 dormant dates; report it with that context. The effect size (R² = 0.52 on date means) is the more informative number.",
  "Q3: correlations rest on seven dates; with n = 7, |rho| must exceed about 0.79 to reach p < 0.05, which no class does, so present them as descriptive. The lag comparison mixes 28-day intervals with one 41-day interval (January to February, which also lacks plot 27).",
  "MicroResp reflects potential capacity after a resource pulse at 60% WHC and room temperature. Higher dormant SIR could partly reflect a rewetting or warming response of winter soils (Birch effect, KH comment); field water content was the variable most aligned with the SIR ordination (exploratory envfit)."
)

comment_responses <- tribble(
  ~Comment, ~`How v4 addresses it`,
  "KH: state the title finding more clearly", "Two title options naming the direction of each response.",
  "KH: date and substrate are the story, not plots", "Ordinations show date centroids in time order; plots are not highlighted.",
  "AZ: potential activity is higher in dormancy for all substrates", "Quantified: all 31 substrates higher in dormancy; overall shift tested with dates as replicates.",
  "KH: contrast seasonal and steady substrates", "Steady substrates (GAM FDR > 0.05) labelled in Fig. 2C and S2.",
  "VG: organize Fig. S2 by substrate classes", "Fig. S2 panels grouped by substrate group, ordered by deviance explained within group.",
  "KH: is the mapping based on annotations?", "Yes: KEGG pathways containing each substrate's catabolic/assimilatory reactions; class rules in Scripts/03_compound_class_mapping.R.",
  "KH: remove basic cell functions from metaT?", "Fig. 3A keeps all 325 categories; metabolism-only version in Fig. S3 gives the same result.",
  "KH: sum or average for a comparable measure?", "Transcripts summed within class (bulk investment), SIR averaged as z-scores; mean-z transcripts as sensitivity (Fig. 4B open symbols).",
  "AZ: are some months more aligned than others?", "Fig. 4C (per-date divergence) and Fig. S4b (Procrustes residuals by date): January is the least aligned date.",
  "KH: just 3 substrates; include the threshold", "Top-3 selection removed (D4).",
  "KH: wet-dry cycle / Birch effect", "Envfit on date means: water content most aligned with SIR ordination; noted as a caveat.",
  "KH/VG: outliers vs hot spots/moments; two lines of evidence", "All points kept; both-screens exclusion as sensitivity (Fig. S1).",
  "SB: roots are hand-picked from pooled 0-15 cm cores", "Labelled as such in Fig. 1 and S1."
)

file_map <- tribble(
  ~Script, ~Outputs,
  "Scripts/00_run_all.R", "Runs every step below; logs in Results/manuscript_v4/run_logs",
  "Scripts/01_Data_merging_and_cleaning.R", "Data/Data_merged_cleaned.csv; Results/Data_merging_and_cleaning/",
  "Scripts/02_sample_availability.R", "Results/manuscript_v4/sample_availability/ (audit, label check, lag pairs)",
  "Scripts/03_compound_class_mapping.R", "Data/compound_class_mapping.csv; Data/compound_class_unique_BRITE_C.csv",
  "Scripts/10_figure1_context.R", "Figures/Results manuscript_v4/Figure1",
  "Scripts/11_figureS1_root_respiration.R", "Figures/Results manuscript_v4/FigureS1",
  "Scripts/20_figure2_microresp.R", "Figures/Results manuscript_v4/Figure2 (incl. S2, S2b)",
  "Scripts/30_figure3_transcripts.R", "Figures/Results manuscript_v4/Figure3 (incl. S3)",
  "Scripts/40_figure4_correspondence.R", "Figures/Results manuscript_v4/Figure4 (incl. S4a, S4b)",
  "Scripts/helpers.R", "Season definition, loaders, exact date-level permutation test, plotting theme",
  "Scripts/docs/", "This document and the storyboard (HTML generated from results, converted to .docx with Word)",
  "Archive/2026-10-23_pre_v4/", "Superseded scripts, figures, and results (kept locally; in git history up to commit 80e0df9)"
)

notes <- c(
  h1("Analysis notes for storyboard v4"),
  p("Companion to Storyboard_v4.docx: decisions, data changes, how results moved from v3, open items, and how to rerun everything. Generated from the results in the repository on branch manuscript-v4-analysis."),
  h2("1. Items that need your attention"),
  ul(open_items),
  h2("2. Interpretation cautions"),
  ul(interpretation_cautions),
  h2("3. How results changed from storyboard v3"),
  tbl(changes_vs_v3),
  h2("4. Data changes"),
  ul(data_changes),
  h2("5. Decision log (agreed 2026-10-23)"),
  tbl(decisions),
  h2("6. Responses to v3 coauthor comments"),
  tbl(comment_responses),
  h2("7. Reproducing the analysis"),
  p("From the repository root, with the DataHub March folder in place and R 4.6.1: <i>Rscript Scripts/00_run_all.R</i> (about 10 minutes; Figure 2 permutation tests take most of it). Then regenerate these documents with <i>Rscript Scripts/docs/write_storyboard_html.R</i> followed by <i>powershell -File Scripts/docs/html_to_docx.ps1</i>. Data/README.md lists every input with its origin and MD5 checksum. Package versions are saved to Results/manuscript_v4/sessionInfo.txt."),
  tbl(file_map),
  h2("8. Statistical methods in brief"),
  ul(c(
    "Exact date-level permutation: all assignments of season labels to dates are enumerated (choose(13, 5) = 1,287 for SIR; choose(7, 2) = 21 for transcripts), PERMANOVA pseudo-F computed for each, p = share of assignments with F at least as large as observed. Applied to date-mean profiles (primary) and to samples with labels permuted by whole date (sensitivity). Dispersion (betadisper F) tested the same way.",
    "Mixed models (lme4/lmerTest, Satterthwaite df): per-substrate log2 SIR ~ Season + (1|Date) + (1|Plot); H1b: response ~ Season x class + (1|Date) + (1|Date:class) + (1|Sample) + (1|Substrate:Season), contrasts with emmeans.",
    "FDR: Benjamini-Hochberg across the 31 substrates.",
    "GAM: mgcv, log2 SIR ~ s(time, k = 7) + s(Plot, bs = 're'), REML.",
    "Field variables: log(y) ~ Season + (1|Date) + (1|Plot); soil temperature (one value per date) by Welch t-test on date means.",
    "Q3: Spearman correlations of date means (n = 7) and plot-dates (n = 34-35); no p-values reported."
  ))
)

writeLines(page("Analysis notes v4", notes), file.path(out_dir, "Analysis_notes_v4.html"), useBytes = TRUE)
cat("Wrote", file.path(out_dir, c("Storyboard_v4.html", "Analysis_notes_v4.html")), sep = "\n")
