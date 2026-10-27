# =============================================================================
# Write storyboard v5 and its companion notes as Word documents (.docx)
# =============================================================================
#
# Storyboard for the locked analysis plan v1 (analysis_plan_v1.md), in plain
# language. Framing, questions, and hypotheses are reproduced from
# VGC_Revised_Framing_Q_H.docx. Every number is read from Results/plan_v1
# (and, for the field-conditions figure, Results/manuscript_v4/Figure1).
# Run after Scripts/plan_v1/01-08:
#   Rscript Scripts/docs/write_storyboard_v5_docx.R

rm(list = ls())
suppressPackageStartupMessages(library(tidyverse))
source(file.path("Scripts", "docs", "docx_writer.R"))

out_dir <- "Docs"
res <- function(f) read_csv(file.path("Results", "plan_v1", f), show_col_types = FALSE)
fig <- function(f) normalizePath(file.path("Figures", "plan_v1", f), winslash = "/")

f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
pv <- function(p) if_else(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, format = "fg", digits = 2)))
ci <- function(lo, hi, d = 2) sprintf("95%% CI %s to %s", f2(lo, d), f2(hi, d))

# -----------------------------------------------------------------------------
# Results used in the text
# -----------------------------------------------------------------------------

counts <- res("cp1_sample_counts.csv") |> filter(Irrigation_trmt_pct == 100)
n_mr <- sum(counts$microresp)
n_tr <- sum(counts$transcripts)

season <- res("cp2_t1_t3_t6_season_tests.csv")
st <- function(test_name, norm = "per_C") filter(season, test == test_name, normalization == norm)
t1 <- st("1_H1a_simple_score")
t3 <- st("3_H1b_polymer_minus_simple")
tp <- st("polymer_score (for forest plot)")

roots <- res("cp2_t2_root_biomass_score.csv")
root_n <- res("cp2_t2_root_biomass_compound_counts.csv")
root_flags <- res("root_flags_active_100pct.csv") |> filter(root_flag)

coupling <- res("cp2_t4_t6_coupling.csv") |> filter(normalization == "per_C")
cp_all <- filter(coupling, dates == "all 13")
cp_dorm <- filter(coupling, dates == "dormant 5")

comp <- res("cp2_t5_compound_effects.csv")
simple_eff <- filter(comp, group == "Simple")
poly_eff <- filter(comp, group == "Polymer")
n_pos <- sum(comp$diff[comp$group != "Reference"] > 0)
ratio_simple <- simple_eff |> filter(ratio_induced_C > 0) |> pull(ratio_induced_C)
nice <- function(x) str_replace_all(x, c("_acid" = " acid", "_" = " ", "^D\\." = "", "^L\\." = "",
                                         "^N\\.Acetyl\\.D\\.glucosamine" = "N-acetylglucosamine",
                                         "Glucosamine hydrochloride" = "glucosamine")) |> tolower()
top_simple <- simple_eff |> arrange(desc(diff)) |> slice(1:3)
low_poly <- poly_eff |> arrange(diff) |> slice(1)

mixed <- res("cp2_mixed_models.csv")
mx_h1a <- filter(mixed, str_starts(model, "H1a"), term == "seasonDormant")
mx_int <- filter(mixed, term == "seasonDormant:groupPolymer")
window <- res("cp2_window_permutation.csv")
robust <- res("cp2_robustness.csv")
rb <- function(test_name, prefix) filter(robust, test == test_name, str_starts(variant, prefix))
loo1 <- robust |> filter(test == "1_simple_score", str_starts(variant, "e:"))
rep25 <- res("cp2_25pct_replication.csv")
enz <- res("cp2_t7_enzyme_season.csv")
enz_sub <- res("cp2_t7_enzyme_substrate.csv")
ratio8 <- res("cp2_t8_rate_ratio_season.csv")
field_cor <- res("cp2_explore_field_correlations.csv")
temp_simple <- filter(field_cor, score == "simple_score", variable == "SoilTemp_mean")
peaks <- res("cp2_explore_polymer_peaks.csv")

v4_f1 <- read_csv(file.path("Results", "manuscript_v4", "Figure1", "Figure1_season_tests.csv"),
                  show_col_types = FALSE)
f1v <- function(col) filter(v4_f1, column == col)

bridge <- res("cp4_bridge_q1.csv")
h2a <- res("cp4_h2a_expression_season.csv")
h2a_s <- filter(h2a, score == "Simple")
h2a_p <- filter(h2a, score == "Polymer")
h2a_d <- filter(h2a, score == "Polymer_minus_Simple")
h2b_sum <- res("cp4_h2b_summary_and_robustness.csv")
h2b_main <- h2b_sum[1, ]
h2b_perC <- filter(h2b_sum, str_starts(variant, "a"))
h2b_high <- filter(h2b_sum, str_starts(variant, "b"))
h2b_sens <- res("sens_q2_h2b_summary.csv") |> filter(str_starts(variant, "without"))
h2b_eff <- res("cp4_h2b_compound_effects.csv")
lag <- res("cp4_h2c_lag_advantage.csv")
step <- res("cp4_h2c_step_permutation_and_robustness.csv")
step_main <- step[1, ]
step_sens <- res("sens_q2_h2c_step_permutation.csv") |> filter(str_starts(variant, "without"))
lag_models <- res("cp4_lag_models_compare.csv")
lm_cnt <- lag_models |> group_by(pathway_step) |>
  summarise(n = n(), lag_better = sum(AIC_lag_minus_same < 0), .groups = "drop")
lag_ci_pos <- lag |> filter(ci_low > 0) |> pull(set)
lag_ci_neg <- lag |> filter(ci_high < 0) |> pull(set)
lag_pos <- lag |> group_by(pathway_step) |> summarise(n = n(), pos = sum(lag_advantage > 0))

gene_sets <- res("KO_gene_sets_APPROVED.csv")
set_status <- res("KO_gene_sets_RECOMMENDED_expression.csv")

sens_q1 <- res("sens_q1_test6_season.csv")
t6 <- st("6_H1a_simple_score", "per_MBC")

# -----------------------------------------------------------------------------
# Tables
# -----------------------------------------------------------------------------

compound_notes <- c(
  D.Fructose = "Fruit sugar; common in root exudates",
  D.Xylose = "Five-carbon sugar; building block of hemicellulose",
  D.Trehalose = "Two-glucose sugar that microbes make to survive stress",
  D.Cellobiose = "Two-glucose sugar released when cellulose is broken down",
  Citric_acid = "Organic acid released by roots",
  Malic_acid = "Organic acid released by roots",
  Oxalic_acid = "Organic acid released by roots and fungi",
  Succinic_acid = "Organic acid; also an everyday step in cell metabolism",
  Lactic_acid = "Organic acid made by fermenting microbes",
  Glutamic_acid = "Amino acid; common in root exudates",
  L.Arginine = "Amino acid rich in nitrogen",
  L.Glutamine = "Amino acid; a main carrier of nitrogen in cells",
  Glycine = "Smallest amino acid",
  L.Alanine = "Amino acid; common in root exudates",
  L.Leucine = "Amino acid with a branched carbon chain",
  Cellulose = "Main structural polymer of plant cell walls",
  Hemicellulose = "Plant cell-wall polymer of mixed sugars (mostly xylose)",
  Lignin = "Woody, ring-rich plant polymer; hard to break down",
  Chitin = "Polymer of fungal cell walls and insect shells",
  Chitosan = "Chitin with acetyl groups removed",
  D.Glucose = "Most common sugar; used as a reference only",
  N.Acetyl.D.glucosamine = "Building block of chitin (NAG)",
  Glucosamine_hydrochloride = "Building block of chitosan"
)

group_label <- c(Simple = "Simple compound", Polymer = "Polymer",
                 Reference = "Reference", Q2 = "Chitin and chitosan building block")

set_n <- set_status |> select(set, n_KOs, status)
table1 <- tibble(Substrate = names(compound_notes), `What it is` = unname(compound_notes)) |>
  left_join(comp |> select(Substrate, group), by = "Substrate") |>
  mutate(group = replace_na(group, "Q2")) |>
  left_join(set_n, by = c("Substrate" = "set")) |>
  transmute(
    Compound = str_replace_all(Substrate, c("_" = " ", "\\." = "-")) |>
      str_replace("^N-Acetyl-D-glucosamine", "N-acetyl-D-glucosamine") |>
      str_replace(" hydrochloride", ""),
    Group = group_label[group],
    `What it is` = `What it is`,
    `Gene set for Q2` = case_when(
      group == "Reference" ~ "Not used",
      is.na(status) ~ "None",
      status == "use in Q2" ~ sprintf("Yes (%d genes)", n_KOs),
      status == "no specific gene set" ~ "No specific genes",
      .default = "Genes barely expressed; not used"
    )
  )

enzyme_names <- c(
  CB = "Cellobiohydrolase: cuts cellobiose off cellulose",
  BG = "Beta-glucosidase: splits cellobiose into glucose",
  BX = "Beta-xylosidase: releases xylose from hemicellulose pieces",
  NAG = "N-acetylglucosaminidase: releases NAG from chitin pieces",
  Leu = "Leucine aminopeptidase: releases amino acids from proteins",
  Ala = "Alanine aminopeptidase: releases amino acids from proteins",
  AP = "Phosphatase: releases phosphate from organic matter"
)
table2 <- tibble(Enzyme = names(enzyme_names), `What it does` = unname(enzyme_names)) |>
  left_join(set_n, by = c("Enzyme" = "set")) |>
  mutate(`Genes in set` = as.character(n_KOs)) |>
  select(Enzyme, `What it does`, `Genes in set`)

glance <- tribble(
  ~Hypothesis, ~`What we asked`, ~`What we found`, ~Verdict,
  "H1a", "Do soils respond less to simple compounds in winter? In summer, do they respond more where there are more roots?",
  sprintf("Responses to simple compounds were higher in winter, not lower (+%s, %s). Roots showed no link (rho = %s).",
          f2(t1$diff), pv(t1$p_welch), f2(roots$rho[1])),
  "Not supported. The change went the opposite way.",
  "H1b", "Do polymer responses rise relative to simple-compound responses in winter?",
  sprintf("Polymers rose about as much as simple compounds (difference %s, %s).", f2(t3$diff), pv(t3$p_welch)),
  "Not supported.",
  "Priming / coupling", "Do the two groups move up and down together?",
  sprintf("Yes. Dates with high simple-compound responses also had high polymer responses (r = %s, %s).",
          f2(cp_all$pearson_r), pv(cp_all$pearson_p)),
  "Consistent with the data.",
  "H2a", "In winter, are genes for using simple compounds expressed less, and polymer genes relatively more? (descriptive)",
  sprintf("Simple-compound genes were lower on the two winter dates (%s). Polymer genes did not change (%s).",
          f2(h2a_s$diff), f2(h2a_p$diff)),
  "Partly: simple-compound genes fell as predicted; polymer genes did not rise.",
  "H2b", "For each compound, do gene expression and potential change in the same direction between seasons?",
  sprintf("%d of %d compounds matched; that is not more than chance (%s).",
          h2b_main$n_concordant, h2b_main$n_compounds, pv(h2b_main$binom_p)),
  "Not supported.",
  "H2c", "Does potential match last month's expression better than this month's, and more so for enzymes that work outside the cell?",
  sprintf("Outside-the-cell gene sets leaned toward a lag (difference %s, %s), but the result weakened in checks.",
          f2(step_main$diff), pv(step_main$p_perm)),
  "Not supported at the usual threshold; a hint worth testing again."
)

glossary <- tribble(
  ~Term, ~Meaning,
  "MicroResp (SIR)", "A lab test. A small amount of soil gets one compound, and we measure how much CO2 comes out in 5 hours. It shows how ready the microbes already in the soil are to use that compound. We call this 'potential'. It is not what they were doing in the field.",
  "Water control", "The same soil with only water added. Subtracting it removes the CO2 the soil would release anyway, so what is left is the response to the compound.",
  "Per unit of carbon", "Each response is divided by how much carbon was added with the compound, so compounds added in different amounts can be compared.",
  "Simple compound", "A small molecule (sugar, organic acid, amino acid) that microbes can take straight into the cell.",
  "Polymer", "A large molecule (cellulose, chitin, lignin) that microbes must first cut into pieces outside the cell, using enzymes they release.",
  "z-score / score", "Each compound's values rescaled so its own average is 0 and its typical spread is 1. This puts fast and slow compounds on one scale. A group score is the average z-score of the compounds in that group.",
  "Date mean", "The average of the five or six plots sampled on one day. We compare seasons using date means because all plots on one day share that day's weather and lab run.",
  "Difference (dormant − active)", "Winter average minus summer average. Above 0 means higher in winter.",
  "95% confidence interval (CI)", "The range of values that fits the data. If it does not include 0, the difference is unlikely to be chance.",
  "p-value", "How often a difference this large would appear by chance if seasons did not differ. Below 0.05 is the usual line for 'clear'.",
  "Correlation (r, rho)", "How closely two measurements rise and fall together: 1 = perfectly together, 0 = unrelated, −1 = opposite.",
  "Gene expression (transcripts)", "Which genes microbes were actively using in the field at the moment of sampling, read from RNA.",
  "Gene set", "The genes for one specific job, such as the enzymes that cut chitin, or the transporter and first step for taking up one sugar.",
  "Outside / inside the cell", "Enzymes that cut polymers work outside the cell and can last for days to weeks. Uptake and metabolism of simple compounds happen inside the cell, where molecules turn over quickly.",
  "Lag advantage", "How much better potential matches the previous month's gene expression than the same month's. Above 0 means the previous month fits better.",
  "Microbial biomass carbon (MBC)", "How much carbon is held in living microbes; a rough measure of how many microbes there are."
)

# -----------------------------------------------------------------------------
# Storyboard
# -----------------------------------------------------------------------------

storyboard <- c(
  h1("Storyboard v5: does winter shift what soil microbes are ready to eat, and do their genes show it?"),
  p('<span class="note">', b("Note: "),
    "v5 was drafted with Claude (Anthropic) from the locked analysis plan (analysis_plan_v1.md) on GitHub branch analysis-plan-v1. ",
    "Every number is read from the result tables in Results/plan_v1. The framing, questions, and hypotheses are reproduced from VGC_Revised_Framing_Q_H.docx. ",
    "This version is written in plain language. Checks, side results, and data issues are in the companion document Storyboard_v5_notes.docx.</span>"),

  h2("Suggested Title"),
  p("Option 1: In winter, soil microbes are more ready to use simple compounds, even as they express fewer genes to use them"),
  p("Option 2: Plant dormancy raises soil microbial readiness for simple and complex carbon alike, without a switch toward polymers"),
  p(i("v4 title for reference: Plant dormancy strengthens soil microbial responses to added substrates while field expression of related functions dips in midwinter.")),

  h2("Potential Journals"),
  p("SBB"),
  h2("Team members"),
  p("Vanessa Garayburu-Caruso, Josue Rodriguez-Ramos, Sheryl Bell, Evan Warburton, Sharon Zhao, Regan McDearis, David Degnan, Amy Zimmerman, and Kirsten Hofmockel"),
  h2("Acknowledgements"),
  p("Gregory, Kali Cornwell, Mischelle, other people that contributed? Trinidad (viromes; not in this paper)"),
  h2("GitHub link with code"),
  p("vgc40/Soil_time_series_data (https://github.com/vgc40/Soil_time_series_data); analysis on branch analysis-plan-v1"),

  h2("Framing"),
  h3("P1. Seasonal changes in microbial substrate availability and function"),
  p("Soil microbial function varies across the annual cycle as plant activity, temperature, moisture, and resource availability change. As plants transition between active growth and dormancy, changes in root exudation, root turnover, and residue inputs may alter the timing and relative availability of different substrates, including readily available compounds and more complex organic matter associated with plant residues and microbial necromass. Dormancy may therefore change the mixture of substrates available to microbes rather than simply reducing their supply. These shifts can alter microbial physiological state and resource acquisition strategies, potentially affecting how microbes respond to different substrates throughout the year. Yet much of our understanding of soil microbial function comes from growing-season measurements, despite evidence that microbial activity can continue outside the growing season. As a result, we know comparatively little about how microbial responses to different substrates change during plant dormancy, limiting our understanding of how plant seasonality shapes belowground microbial metabolism."),
  h3("P2. Connecting microbial substrate-use potential and field expression"),
  p("Commonly used measurements capture different aspects and timescales of microbial function. Field soil respiration integrates root and microbial activity under prevailing conditions but cannot resolve which substrates microbes are using. Laboratory substrate-induced respiration measures short-term responses to added compounds under standardized conditions. Because these responses are measured within 5 h, they primarily characterize the short-term respiratory capacity of the existing microbial community, which reflects its physiological state and antecedent environmental conditions rather than substantial growth during the assay. Seasonal changes in these responses therefore show how microbial substrate-use potential shifted over the year. Metatranscriptomes, in contrast, characterize which microbial functions are being expressed under field conditions at the time of sampling, but not how strongly microbes would respond to specific substrate additions. Relationships between transcripts and the processes they encode are often weak (Rocca et al. 2015), and timing may be one reason. Microbial transcripts typically persist for minutes, whereas soil extracellular enzymes can remain active for days to weeks after they are produced. This difference raises the possibility that substrate-use potential reflects expression from preceding weeks rather than the moment of sampling. If so, the effect should be strongest for the extracellular enzymes that break down polymers, since uptake and metabolism of simple compounds rely on intracellular machinery that turns over faster. Expression and potential could then shift in the same direction across seasons yet differ in monthly timing, with the lag depending on the pathway step. Whether such lags occur in field soils has, to our knowledge, not been tested. Comparing the two across related compounds and metabolic pathways can help clarify whether plant dormancy changes which functions microbes express, their potential to respond as substrate availability changes, or both."),
  h3("P3. Study objectives"),
  p("Here, we investigate how microbial substrate-use potential and field gene expression change between plant-active and plant-dormant periods, and how the two are related. We focused on fully irrigated (100%) plots in an irrigated tall wheatgrass system, a perennial grass with a plant-active period from April (plant emergence) through October (harvest) and a plant-dormant period from November through March. We combine 13 months of substrate-induced respiration measurements for 31 soil-relevant organic compounds with field metatranscriptomes from seven of these sampling dates (June, July, September, and October 2023; January, March, and May 2024), and use supporting measurements of extracellular enzyme activity, microbial biomass, soil respiration, root biomass, and environmental and biogeochemical conditions to help interpret the observed patterns."),

  h2("Research Questions & Corresponding Hypotheses"),
  h3("Q1. Does plant dormancy shift microbial substrate-use potential away from directly assimilable compounds and toward polymers that require extracellular breakdown?"),
  p(b("H1a, legacy of simple-compound supply. "), "Responses to simple compounds will be lower in dormant months than in active months. Within the active period, they will increase with root biomass. Mechanism: Root exudation maintains the microbial machinery for taking up and metabolizing these compounds, and that capacity declines when exudation stops (Zhalnina et al. 2018; Žifčáková et al. 2016)."),
  p(b("H1b, substrate switching. "), "Polymer responses will increase relative to simple-compound responses in dormancy. Mechanism: When labile carbon is scarce, investing in enzymes that break down polymers becomes worthwhile (Kaiser et al. 2010; Žifčáková et al. 2017; Allison & Vitousek 2005, to verify)."),
  p(b("Alternative, priming/coupling. "), "Simple-compound supply fuels polymer breakdown, so the two groups move together and their relative balance doesn't change."),
  h3("Q2. How do seasonal changes in field gene expression relate to changes in substrate-use potential?"),
  p(b("H2a, direction of expression (descriptive). "), "In dormancy, expression of genes for using simple compounds will decrease. Expression of polymer-degrading genes will increase relative to them. Mechanism: the same as H1. As exudation declines, microbes shift investment toward breaking down complex substrates (Žifčáková et al. 2016, 2017). Descriptive because only 2 of the 7 transcript dates fall in dormancy."),
  p(b("H2b, seasonal concordance. "), "For matched compounds, expression and potential will shift in the same direction between active and dormant periods. Alternative: discordance, meaning they shift in different directions, or only one of them shifts."),
  p(b("H2c, monthly timing. "), "Substrate-use potential will agree better with the previous month's expression than with same-date expression. The lag will be stronger for polymer breakdown (extracellular steps) than for uptake and metabolism of simple compounds (intracellular steps). Alternative: same-date tracking, with no lag. That's the expected result for intracellular steps."),

  h2("Results at a glance"),
  p("Each row is one hypothesis: what we asked, what the data showed, and the verdict. The sections below walk through each result."),
  tbl(glance, rel = c(0.8, 2.2, 2.6, 1.6)),

  h2("The story in brief"),
  p("Read this first. It is the train of thought of the paper in seven steps."),
  ul(c(
    sprintf("%s In winter the field went quiet: soils were %s °C colder and gave off about %d%% less CO2. But the microbes were still there in similar amounts (Figure 1).",
            b("1. Winter slows the field, not the microbes."),
            f2(abs(f1v("SoilTemp_mean")$estimate), 0), round(100 * (1 - f1v("Resp_rate_mg_C_m2_h")$ratio))),
    sprintf("%s We expected microbes to lose their readiness for simple compounds once roots stopped feeding them. The opposite happened: winter soils responded more strongly to simple compounds (Figure 2). Within the summer, soils with more roots did not respond more (Figure 3).",
            b("2. Microbes were more ready, not less.")),
    sprintf("%s Responses to polymers rose in winter too, by about the same amount. So microbes did not switch from simple to complex food. Simple and polymer responses rose and fell together through the year (Figure 4).",
            b("3. No switch to polymers.")),
    sprintf("%s In the field, genes for taking up and using simple compounds were expressed less on the two winter dates. Genes for cutting polymers did not change (Figure S5).",
            b("4. Genes told a different story.")),
    sprintf("%s Compound by compound, gene expression and lab readiness did not change in the same direction more often than chance (Figure 5).",
            b("5. Genes and readiness did not move together.")),
    sprintf("%s Readiness matched the previous month's genes somewhat better for enzymes that work outside the cell. This was a hint only; it did not hold up firmly in checks (Figures 6 and 7).",
            b("6. A possible time lag, not proven.")),
    sprintf("%s Winter readiness to use carbon is high across the board and is not explained by roots or by same-month gene expression. What drives it is the open question this paper leaves.",
            b("7. Take-away."))
  )),

  h2("Approach"),
  p("This section explains, in plain terms, what was measured and how. Technical terms are explained in the glossary at the end."),
  h3("1. Samples and seasons"),
  p(sprintf("Soil was sampled once a month from June 2023 to June 2024 (13 dates) in the six fully irrigated plots of a tall wheatgrass field. Dates from April to October count as the plant-active season (8 dates) and dates from November to March as the plant-dormant season (5 dates). The lab test (MicroResp) was run on %d samples; one plot was not tested in February 2024. Field gene expression was measured on %d samples: five plots on seven of the dates, two of them in winter (January and March 2024).", n_mr, n_tr)),
  h3("2. Why we compare dates, not plots"),
  p("Season comes with the sampling day. All plots sampled on one day share that day's weather, the same field trip, and the same lab run. So each day counts as one piece of evidence about season, and we compare the 5 winter days with the 8 summer days. Plots are averaged within each day first."),
  h3("3. Measuring readiness (Q1)"),
  p("For each sample and compound, we took the CO2 released in 5 hours, subtracted the CO2 from the same soil with only water, and divided by the amount of carbon added. That gives the response to the compound. Compounds respire at very different speeds, so each compound was rescaled to its own average (a z-score). We then averaged the 15 simple compounds into a simple score and the 5 polymers into a polymer score, for each sample and then each date."),
  h3("4. Measuring gene expression (Q2)"),
  p("For each compound, we listed the genes specific to using it, following a written rule. For polymers, these are the enzymes that cut them. For simple compounds, the transporter that brings them into the cell plus the first step that commits them to being used. Shared, general metabolism was left out. Genes for the seven measured soil enzymes formed seven more sets. Three compounds had no usable set: succinic acid has no specific genes, and the lignin and chitosan genes were barely expressed. Expression of each set was expressed as its share of all expression in that sample, on a log scale."),
  h3("5. Linking the two (Q2)"),
  p("We compared each compound's change between seasons in genes and in readiness (H2b). We also asked whether readiness matched gene expression from the same month or from the month before (H2c). In the plan, readiness for Q2 is divided by microbial biomass. As explained below, two samples with almost no measured biomass dominated that version, so we also show the version without that step."),

  h2("Key Figures"),

  h3("Figure 1. What were field conditions like in each season?"),
  img(normalizePath(file.path("Figures", "manuscript_v4", "Figure1", "Figure1_main.png"), winslash = "/")),
  p('<span class="caption">', b("How to read Figure 1. "),
    "From the v4 analysis, unchanged. (A) Sampling dates and which ones have gene data. (B-I) Field measurements through the year. Grey points are plots; the black line and band are the date average and its range. Shading marks the active (green) and dormant (blue-grey) seasons.</span>"),
  p(b("What it shows. ")),
  ul(c(
    sprintf("Soil was %s °C colder on winter dates.", f2(abs(f1v("SoilTemp_mean")$estimate), 1)),
    sprintf("Field CO2 release in winter was %d%% of the summer rate.", round(100 * f1v("Resp_rate_mg_C_m2_h")$ratio)),
    sprintf("Microbial biomass (%s times), root biomass (%s times), and soil enzymes did not differ clearly between seasons.",
            f2(f1v("MBC_ug_gdw_bgc")$ratio), f2(f1v("Roots_gdw")$ratio))
  )),
  p(b("Why it matters. "), "The microbes were present in winter in similar numbers. Any change in how they respond to added compounds is not simply because there were more or fewer of them."),

  h3("Figure 2. Q1, H1a: Were winter soils less ready to use simple compounds?"),
  img(fig("Fig1_season_simple_score.png")),
  p('<span class="caption">', b("How to read Figure 2. "),
    "Left panel: each coloured point is one sampling date; its height is the simple score for that date. Above 0 means a stronger response than that compound's yearly average. Black points and bars are each season's average and its 95% range. Right panel: the same, after dividing by microbial biomass; open circles are dates that include a sample with almost no measured biomass.</span>"),
  p(b("What we found. ")),
  ul(c(
    sprintf("Winter dates had higher simple scores than summer dates, by %s z-units (%s, %s; %d winter and %d summer dates).",
            f2(t1$diff), ci(t1$ci_low, t1$ci_high), pv(t1$p_welch), t1$n_dormant, t1$n_active),
    sprintf("In everyday units, most simple compounds released %s to %s times as much CO2 per unit carbon in winter as in summer.",
            f2(min(ratio_simple), 1), f2(max(ratio_simple), 1)),
    sprintf("The result held when we changed how the response was calculated, dropped any one compound, or used a statistical model that includes every plot (same size, %s, but %s because that model is more cautious).",
            f2(mx_h1a$estimate), pv(mx_h1a$p)),
    sprintf("It did not hold for the later reading (5 to 24 hours; difference %s, %s).", f2(rb("1_simple_score", "c:")$diff), pv(rb("1_simple_score", "c:")$p_welch)),
    sprintf("After dividing by microbial biomass, the difference shrank to %s (%s). That version is dominated by two samples with almost no measured biomass (see notes).",
            f2(t6$diff), pv(t6$p_welch))
  )),
  p(b("Verdict on H1a, first part. "), "Not supported. Readiness for simple compounds went up in winter, not down."),
  p(b("What this means. "), "Losing root inputs in winter did not reduce the microbes' readiness to use the compounds roots release. One possibility is that cold, carbon-starved microbes are 'hungry' and respond strongly to any carbon offered. Another is that the lab test, run warm and moist, wakes winter soil up more. This design cannot tell these apart."),

  h3("Figure 3. Q1, H1a: In summer, did soils with more roots respond more?"),
  img(fig("Fig2_simple_score_root_biomass.png"), 4),
  p('<span class="caption">', b("How to read Figure 3. "),
    "Each point is one summer sample. Across: root biomass (log scale). Up: simple score. Crosses are cores with unusually many roots. The dashed line uses all cores; the solid line leaves the crosses out.</span>"),
  p(b("What we found. ")),
  ul(c(
    sprintf("No link between roots and the simple score: rho = %s (%s, %d samples). Leaving out the %d root-dense cores did not change this (rho = %s, %s).",
            f2(roots$rho[1]), pv(roots$p[1]), roots$n[1], nrow(root_flags), f2(roots$rho[2]), pv(roots$p[2])),
    sprintf("Compound by compound, %d of 15 simple compounds had a positive link, about what chance alone would give.", root_n$n_positive[1])
  )),
  p(b("Verdict on H1a, second part. "), "Not supported. Roots in a core did not predict how ready its microbes were."),
  p(b("What this means. "), "Root mass in a soil core may not reflect how much roots are releasing at that moment, or exudates may affect a thin zone around roots that a bulk core dilutes."),

  h3("Figure 4. Q1, H1b and the alternative: Did microbes switch toward polymers in winter?"),
  img(fig("Fig3_forest_compound_effects.png"), 4),
  p('<span class="caption">', b("How to read Figure 4. "),
    "Each row is one compound. The dot is the winter-minus-summer difference; the line is its 95% range. Right of the dashed line means higher in winter. Blue: simple compounds. Orange: polymers. Grey: glucose (reference). The bottom rows summarize each group and the polymer-minus-simple difference.</span>"),
  p(b("What we found. ")),
  ul(c(
    sprintf("%d of 20 compounds were higher in winter; only %s was not. The largest increases among simple compounds were %s.",
            n_pos, nice(low_poly$Substrate), paste(nice(top_simple$Substrate), collapse = ", ")),
    sprintf("Polymers as a group rose by %s (%s), about as much as simple compounds (%s).", f2(tp$diff), pv(tp$p_welch), f2(t1$diff)),
    sprintf("So the polymer-minus-simple difference was small: %s (%s, %s). H1b predicted it would be clearly above 0.", f2(t3$diff), ci(t3$ci_low, t3$ci_high), pv(t3$p_welch)),
    sprintf("Simple and polymer scores rose and fell together across the 13 dates (r = %s, %s).", f2(cp_all$pearson_r), pv(cp_all$pearson_p))
  )),
  p(b("Verdict. "), "H1b is not supported: there was no switch toward polymers. The priming/coupling alternative fits the data: the two groups moved together."),
  p(b("What this means. "), "Winter changed how strongly microbes responded to carbon in general, not which kind they favoured. This matches the v4 storyboard, which found every one of the 31 compounds higher in winter."),
  p(b("What this cannot tell us. "), "Supporting the coupling alternative does not show that simple compounds cause polymer breakdown; it only shows the two rise and fall together."),

  h3("Figure 5. Q2, H2b: Did genes and readiness change in the same direction?"),
  img(fig("Fig4_H2b_quadrant.png")),
  p('<span class="caption">', b("How to read Figure 5. "),
    "Each point is one compound. Across: how much readiness changed from summer to winter. Up: how much its genes changed. Points in the top-right or bottom-left boxes changed the same way in both (a match). Left panel: readiness divided by microbial biomass, as in the plan. Right panel: readiness per unit carbon, without that step.</span>"),
  p(b("What we found. ")),
  ul(c(
    sprintf("Left panel: %d of %d compounds matched (%s); the correlation was weak (rho = %s, %s). Almost every compound sits at the same spot across, because one sample with almost no biomass drives all of them.",
            h2b_main$n_concordant, h2b_main$n_compounds, pv(h2b_main$binom_p), f2(h2b_main$spearman_rho), pv(h2b_main$spearman_p)),
    sprintf("Right panel: %d of %d matched (%s; rho = %s).", h2b_perC$n_concordant, h2b_perC$n_compounds, pv(h2b_perC$binom_p), f2(h2b_perC$spearman_rho)),
    sprintf("Without the two low-biomass samples: %d of %d matched (%s).", h2b_sens$n_concordant, h2b_sens$n_compounds, pv(h2b_sens$binom_p)),
    "A typical pattern (right panel): readiness for glycine, alanine, and leucine went up in winter while their genes went down."
  )),
  p(b("Verdict on H2b. "), "Not supported. Genes and readiness did not change in the same direction more often than chance."),
  p(b("What this means. "), "Of the outcomes listed in the hypotheses, the closest is 'only potential shifts': readiness rose in winter while the genes for many of the same compounds were expressed less. Readiness may reflect machinery built earlier, rather than what microbes were doing at the moment of sampling."),

  h3("Figure 6. Q2, H2c: Did readiness match last month's genes better than this month's?"),
  img(fig("Fig5_H2c_lag_advantage.png"), 4),
  p('<span class="caption">', b("How to read Figure 6. "),
    "Each row is one gene set. The dot is the lag advantage: how much better readiness matched the previous month's genes than the same month's. Right of 0 means the previous month fit better. Lines are 95% ranges. Top: enzymes that work outside the cell. Bottom: uptake and use inside the cell. Coloured vertical lines are each group's average.</span>"),
  p(b("What we found. ")),
  ul(c(
    sprintf("Outside the cell, %d of %d gene sets fit the previous month better; inside the cell, %d of %d.",
            lag_pos$pos[lag_pos$pathway_step == "extracellular"], lag_pos$n[lag_pos$pathway_step == "extracellular"],
            lag_pos$pos[lag_pos$pathway_step == "intracellular"], lag_pos$n[lag_pos$pathway_step == "intracellular"]),
    sprintf("On average, outside-the-cell sets had a larger lag advantage than inside-the-cell sets (by %s, %s). This is the direction H2c predicted, but it is above the usual 0.05 line.",
            f2(step_main$diff), pv(step_main$p_perm)),
    sprintf("Only two sets showed a clear lag on their own: the %s enzymes. One showed the reverse: the N-acetylglucosamine uptake set.",
            paste(lag_ci_pos, collapse = " and ")),
    sprintf("Without the two low-biomass samples, the group difference shrank to %s (%s).", f2(step_sens$diff), pv(step_sens$p_perm)),
    sprintf("A second check, sample by sample, pointed the other way: the previous month fit better for %d of %d inside-the-cell sets but only %d of %d outside-the-cell sets.",
            lm_cnt$lag_better[lm_cnt$pathway_step == "intracellular"], lm_cnt$n[lm_cnt$pathway_step == "intracellular"],
            lm_cnt$lag_better[lm_cnt$pathway_step == "extracellular"], lm_cnt$n[lm_cnt$pathway_step == "extracellular"])
  )),
  p(b("Verdict on H2c. "), "Not supported at the usual threshold. The main test leaned the predicted way, but the checks did not agree."),
  p(b("What this cannot tell us. "), "Seven dates give only seven pairs to correlate, and the months between samples were not all the same length. A lag would need denser sampling to test well."),

  h3("Figure 7. Example: the chitin pathway over the year"),
  img(fig("Fig6_example_chitin_NAG_pathway.png")),
  p('<span class="caption">', b("How to read Figure 7. "),
    "Every line is rescaled to its own average (0) so they can share one axis; only timing can be compared, not size. (A) Chitin-cutting genes, chitin readiness, and the soil enzyme that releases the chitin building block (NAG). (B) Genes for taking up and using NAG, and NAG readiness. Grey shading is winter.</span>"),
  p(b("What it shows. ")),
  ul(c(
    "Chitin-cutting genes were lowest in January. Chitin readiness (per unit carbon) peaked in December and NAG enzyme activity in January.",
    "The versions divided by biomass are flat except for two spikes (February and May 2024), each caused by one sample with almost no measured biomass."
  )),
  p(b("Why it matters. "), "It shows the general Q2 result in one pathway: winter readiness was high while the matching genes were low."),

  h3("Supplementary figures"),
  p(b("Figure S1. Soil enzymes by season (Q1, interpretation). "),
    sprintf("None of the four enzymes that cut cellulose, hemicellulose, or chitin pieces differed between seasons (winter-to-summer ratios %s to %s; all p > 0.3). Their activity did not track readiness for the matching polymer (rho %s to %s, all p > 0.28).",
            f2(min(enz$ratio)), f2(max(enz$ratio)), f2(min(enz_sub$rho)), f2(max(enz_sub$rho)))),
  img(fig("FigS1_enzymes.png")),
  p(b("Figure S2. Early versus later rates (Q1, interpretation). "),
    sprintf("The 5-to-24-hour rate divided by the first-5-hour rate did not differ between seasons for simple compounds (%s) or polymers (%s). So the winter increase was not about microbes growing more during the test.",
            pv(ratio8$p_welch[ratio8$group == "Simple"]), pv(ratio8$p_welch[ratio8$group == "Polymer"]))),
  img(fig("FigS2_rate_ratio_24h_5h.png"), 4),
  p(b("Figure S3. Simple and polymer scores through the year. "), "Both scores climbed through autumn, peaked in winter, and fell in spring, following the same curve."),
  img(fig("FigS3_score_time_series.png")),
  p(b("Figure S4. Samples with almost no measured biomass. "), "Shows why the biomass-divided versions are unreliable: two samples sit far above all others."),
  img(fig("FigS4_mbc_normalized_scores_flags.png")),
  p(b("Figure S5. Gene expression by season (Q2, H2a; descriptive). "),
    sprintf("Simple-compound genes were lower on the two winter dates (%s, %s). Polymer genes did not change (%s, %s). There are only two winter dates, so no p-values are given.",
            f2(h2a_s$diff), ci(h2a_s$ci_low, h2a_s$ci_high), f2(h2a_p$diff), ci(h2a_p$ci_low, h2a_p$ci_high))),
  img(fig("FigS5_H2a_expression_scores.png"), 4),
  p(b("Figure S6. Do the seven gene dates tell the same Q1 story as all 13 dates? "),
    sprintf("Yes. On the seven dates the simple score was again higher in winter (%s, %s), and the polymer-minus-simple difference was again small (%s, %s).",
            f2(bridge$diff[1]), pv(bridge$p_welch[1]), f2(bridge$diff[2]), pv(bridge$p_welch[2]))),
  img(fig("FigS6_bridge_Q1.png"), 4),
  p(b("Figure S7. Sample-by-sample lag check (Q2, H2c). "), "Points above the dashed line fit the previous month better. The fits were weak overall."),
  img(fig("FigS7_lag_models_plot_level.png"), 4),

  h2("Take Home Messages"),
  ul(c(
    sprintf("In winter, the field slowed down (%s °C colder, about %d%% less CO2), but soil microbes became more ready to use added carbon.",
            f2(abs(f1v("SoilTemp_mean")$estimate), 0), round(100 * (1 - f1v("Resp_rate_mg_C_m2_h")$ratio))),
    "They became more ready for simple compounds and polymers alike. There was no switch toward polymers; the two rose and fell together.",
    "Roots did not explain summer differences in readiness.",
    "In the field, genes for using simple compounds were expressed less in winter, the opposite of readiness. Genes and readiness did not change in step.",
    "There is a hint that readiness reflects gene expression from the month before, especially for enzymes that work outside the cell, but it is not firm."
  )),

  h2("Table 1. Compounds in the analysis"),
  p("The 15 simple compounds and 5 polymers make up the Q1 groups. Glucose is shown as a reference only. N-acetylglucosamine and glucosamine are used in Q2, paired with chitin and chitosan."),
  tbl(table1, rel = c(1.4, 1.4, 2.8, 1.4)),

  h2("Table 2. Soil enzymes and their gene sets"),
  tbl(table2, rel = c(0.8, 3.5, 1)),

  h2("Glossary"),
  tbl(glossary, rel = c(1.3, 4))
)

write_docx(storyboard, file.path(out_dir, "Storyboard_v5.docx"))

# -----------------------------------------------------------------------------
# Companion notes
# -----------------------------------------------------------------------------

rob_rows <- tribble(
  ~Check, ~`Simple score (test 1)`, ~`Polymer − simple (test 3)`,
  "Main result", sprintf("%s (%s)", f2(t1$diff), pv(t1$p_welch)), sprintf("%s (%s)", f2(t3$diff), pv(t3$p_welch)),
  "No water subtraction", sprintf("%s (%s)", f2(rb("1_simple_score", "a:")$diff), pv(rb("1_simple_score", "a:")$p_welch)),
  sprintf("%s (%s)", f2(rb("3_polymer_minus_simple", "a:")$diff), pv(rb("3_polymer_minus_simple", "a:")$p_welch)),
  "Per gram of soil instead of per unit carbon", sprintf("%s (%s)", f2(rb("1_simple_score", "b:")$diff), pv(rb("1_simple_score", "b:")$p_welch)),
  sprintf("%s (%s)", f2(rb("3_polymer_minus_simple", "b:")$diff), pv(rb("3_polymer_minus_simple", "b:")$p_welch)),
  "Later reading (5 to 24 h)", sprintf("%s (%s)", f2(rb("1_simple_score", "c:")$diff), pv(rb("1_simple_score", "c:")$p_welch)),
  sprintf("%s (%s)", f2(rb("3_polymer_minus_simple", "c:")$diff), pv(rb("3_polymer_minus_simple", "c:")$p_welch)),
  "Leave out one compound at a time (range)", sprintf("%s to %s", f2(min(loo1$diff)), f2(max(loo1$diff))),
  with(filter(robust, test == "3_polymer_minus_simple", str_starts(variant, "e:")), sprintf("%s to %s", f2(min(diff)), f2(max(diff)))),
  "Drop weak compounds (cellobiose, cellulose, chitin)", sprintf("%s (%s)", f2(rb("1_simple_score", "f:")$diff), pv(rb("1_simple_score", "f:")$p_welch)),
  sprintf("%s (%s)", f2(rb("3_polymer_minus_simple", "f:")$diff), pv(rb("3_polymer_minus_simple", "f:")$p_welch)),
  "25% irrigation plots", sprintf("%s (%s)", f2(rep25$diff[1]), pv(rep25$p_welch[1])), sprintf("%s (%s)", f2(rep25$diff[2]), pv(rep25$p_welch[2]))
)

questions <- c(
  b("Figure 1. ") %>% paste0("I reused the field-conditions figure from v4 (made with the earlier pipeline and its cleaned data table). Should it be redone with the plan v1 inputs, or is v4's version fine as context?"),
  b("Root-dense cores. ") %>% paste0(sprintf("Your framing document mentions three root-dense cores; the plan's rule flags %d (%s). Figure 3 marks all %d. Is that what you want?", nrow(root_flags), paste(root_flags$Sample_cat, collapse = ", "), nrow(root_flags))),
  b("Biomass-normalized potential. ") %>% paste0("The plan uses it as the main potential for Q2, but it is unreliable here (see section 2). The storyboard shows both versions. Do you want to make the per-carbon version the main one? That would be a change to the locked plan and would go in the deviation log."),
  b("Titles. ") %>% paste0("Both title options state the main Q1 result. Edit freely."),
  b("Citations. ") %>% paste0("The hypotheses keep your note 'Allison & Vitousek 2005, to verify'.")
)

notes <- c(
  h1("Notes for storyboard v5"),
  p("Companion to Storyboard_v5.docx. Results, checks, and data issues that support the story but do not belong in the main line. All numbers come from Results/plan_v1; the full record of decisions is Results/plan_v1/deviation_log.md."),

  h2("1. Questions for you"),
  ul(questions),

  h2("2. Why the biomass-divided results are not used in the main story"),
  p("The plan divides readiness by microbial biomass (MBC) for one Q1 check and for all of Q2. Four samples had MBC recorded as 0 (below the lab blank). Following an earlier project rule, these were set to half the lowest measured value (0.07 µg C per g). The lowest measured value was 0.14 (SM1539)."),
  p("Dividing by a number close to zero makes a huge result. SM1625 (May 2024, set to 0.07) and SM1539 (February 2024, 0.14) came out about 1,000 times higher than any other sample. Every biomass-divided score is therefore driven by these two samples (storyboard Figure S4)."),
  p(sprintf("We reran those analyses without them. The Q1 simple-score difference became %s (%s), and the H2b match count %d of %d (%s). But the next-lowest sample, SM1420 (December 2023, MBC 7.78), then took over the Q1 scores in the same way. Removing samples one at a time does not fix the problem, because dividing by biomass always magnifies the samples with the least biomass.",
            f2(sens_q1$diff[sens_q1$test == "6_H1a_simple_score" & str_starts(sens_q1$variant, "without")]),
            pv(sens_q1$p_welch[sens_q1$test == "6_H1a_simple_score" & str_starts(sens_q1$variant, "without")]),
            h2b_sens$n_concordant, h2b_sens$n_compounds, pv(h2b_sens$binom_p))),

  h2("3. Q1 checks"),
  p("Each cell is the winter-minus-summer difference (z-units) and its p-value, on date means. Full table with corrected p-values: cp2_robustness.csv."),
  tbl(rob_rows, rel = c(2.4, 1.6, 1.6)),
  ul(c(
    sprintf("Model with every plot (mixed model): winter effect on simple compounds %s (%s); the polymer-versus-simple interaction %s (%s). Same size as the date-level test, larger p because the model has few dates to work with.",
            f2(mx_h1a$estimate), pv(mx_h1a$p), f2(mx_int$estimate), pv(mx_int$p)),
    sprintf("Sliding-window check: of the %d possible five-month windows, the real winter window had the %s largest simple-score difference. With 13 dates this check cannot go below p = %s.",
            window$n_positions[1], if_else(window$rank_abs[1] == 1, "single", "not the"), f2(window$min_possible_p[1], 2)),
    sprintf("Soil temperature and the simple score were negatively linked across dates (r = %s, %s): colder dates, higher readiness. Exploratory.", f2(temp_simple$pearson_r), pv(temp_simple$pearson_p)),
    sprintf("Peak dates (exploratory): %s.", paste(sprintf("%s %s", tolower(peaks$Substrate), format(peaks$peak_date, "%b %Y")), collapse = "; "))
  )),

  h2("4. Q2 checks"),
  ul(c(
    sprintf("H2b with high-confidence genes only: %d of %d matched (%s). This is the only version at p < 0.05, and it is the one that drops the most genes.",
            h2b_high$n_concordant, h2b_high$n_compounds, pv(h2b_high$binom_p)),
    sprintf("H2c group difference across versions: main %s (%s); per unit carbon %s (%s); high-confidence genes %s (%s); without the two low-biomass samples %s (%s).",
            f2(step$diff[1]), pv(step$p_perm[1]), f2(step$diff[2]), pv(step$p_perm[2]), f2(step$diff[3]), pv(step$p_perm[3]),
            f2(step_sens$diff), pv(step_sens$p_perm)),
    "Sample-by-sample lag models: most fits explained very little (marginal R² mostly below 0.1), and in 21 of 52 fits the model found no difference between dates.",
    "Gene expression units are geTMM (length-corrected, normalized counts), added up for each gene family (KO). The script that did that adding up was not found."
  )),

  h2("5. Gene sets"),
  p("Rules (approved as D8): drop low-confidence genes; drop the branched-chain amino acid transporters from leucine (shared with two other amino acids and 92% of the signal); succinic acid has no specific genes; sets that are zero in 9 or more of 35 samples are not used (lignin, chitosan). The full list with reasons is Results/plan_v1/KO_gene_sets_APPROVED.csv."),

  h2("6. Data issues"),
  ul(c(
    "Carbon added per well (open, waiting on S. Bell): the file's per-carbon values for succinic acid, chitin, and methyl oleate do not match its own carbon column. The analysis uses the column. This does not affect z-scores, so the main Q1 and Q2 results are not affected unless chitin was weighed differently in each well.",
    "The 24 h reading is the rate from 5 to 24 hours (confirmed from the earlier manuscript methods).",
    "February 2024, plot 27 (SM1536) has no MicroResp measurement. Plot 57 has no gene data on any date."
  )),

  h2("7. What changed from storyboard v4"),
  ul(c(
    "New questions and hypotheses (VGC_Revised_Framing_Q_H.docx), set out in a locked plan before analysis.",
    "Responses now have the water control subtracted and are grouped into 15 simple compounds and 5 polymers (v4 used all 31 compounds in a readily-available/complex split).",
    "Q2 now uses compound-specific gene sets rather than KEGG pathway shares, and adds the time-lag test.",
    "The main Q1 finding is the same as v4: winter soils responded more strongly, across compound types."
  )),

  h2("8. Figure files"),
  p("Storyboard figure numbers differ from the plan's file names: Figure 1 = Figures/manuscript_v4/Figure1/Figure1_main; Figure 2 = Fig1_season_simple_score; Figure 3 = Fig2_simple_score_root_biomass; Figure 4 = Fig3_forest_compound_effects; Figure 5 = Fig4_H2b_quadrant; Figure 6 = Fig5_H2c_lag_advantage; Figure 7 = Fig6_example_chitin_NAG_pathway. Supplementary figures keep their plan numbers (FigS1 to FigS7 in Figures/plan_v1)."),

  h2("9. Reproducing"),
  p("From the repository root: run Scripts/plan_v1/01 to 08 in order, then Rscript Scripts/docs/write_storyboard_v5_docx.R.")
)

write_docx(notes, file.path(out_dir, "Storyboard_v5_notes.docx"))
cat("Wrote", file.path(out_dir, c("Storyboard_v5.docx", "Storyboard_v5_notes.docx")), sep = "\n")
