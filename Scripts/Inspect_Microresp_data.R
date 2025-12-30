# ==============================================================================
# PRIORITY ANALYSIS 1: Validate RF Enzyme Patterns with Substrate Utilization
# ==============================================================================

# Load required packages
library(tidyverse)
library(vegan)
library(ggplot2)
library(patchwork)
library(broom)
library(emmeans)

# ------------------------------------------------------------------------------
# 1.1 Load and prepare data
# ------------------------------------------------------------------------------

# Load MicroResp data
microresp <- read_csv("DataHub/Soil_Phenotypes/SoilSFA_A3.1_MicroResp.csv")
names(microresp)[1] = "Sample_cat"
# Load EEA data  
eea <- read_csv("DataHub/Soil_Phenotypes/SoilSFA_A3.1_EEA.csv")

# Load biogeochemistry data for season definitions
bgc <- read_csv("DataHub/Soil_Biogeochemistry/SoilSFA_A3.1_Soil_Biogeochemistry.csv")

# Define seasons based on irrigation period (from README)
# Growing season: April-October (irrigation application period)
# Non-growing season: November-March

bgc <- bgc %>%
  mutate(
    Collection_date = as.Date(Collection_date, format = "%m/%d/%y"),
    Month_num = lubridate::month(Collection_date),
    Season = case_when(
      Month_num >= 4 & Month_num <= 10 ~ "Growing",
      Month_num >= 11 | Month_num <= 3 ~ "Non-Growing",
      TRUE ~ NA_character_
    )
  )

# Join season info to EEA data
eea_season <- eea %>%
  left_join(bgc %>% select(Sample_cat, Season, Month_num), 
            by = "Sample_cat")

# ------------------------------------------------------------------------------
# 1.2 Prepare MicroResp data - filter and reshape
# ------------------------------------------------------------------------------

# Use 5hr timepoint and C-normalized data (per your methods)
microresp_analysis <- microresp %>%
  filter(
    Timepoint == 5,
    Normalization == "ngCO2-C_mgC_h",
    Substrate != "Water"  # Remove water control
  ) %>%
  # Join season info
  left_join(bgc %>% dplyr::select(Sample_cat, Season, Month_num, Irrigation_trmt_pct...7), 
            by = "Sample_cat")

# Define substrate functional categories based on your hypotheses
substrate_categories <- microresp_analysis %>%
  distinct(Substrate, Substrate_category, Substrate_formula, Substrate_NOSC) %>%
  mutate(
    Functional_group = case_when(
      # Labile substrates (root exudate-like)
      Substrate %in% c("D-Glucose", "D-Fructose", "D-Trehalose", 
                       "Citric_acid", "Malic_acid", "Succinic_acid", "Lactic_acid",
                       "L-Leucine", "L-Alanine", "Glycine", "L-Glutamine", "L-Arginine") ~ "Labile",
      # Polymeric/complex substrates
      Substrate %in% c("Cellulose", "Hemicellulose", "Chitin", "Chitosan", 
                       "Lignin", "BSA") ~ "Polymeric",
      # Intermediate complexity
      Substrate %in% c("D-Cellobiose", "D-Xylose", 
                       "N-Acetyl-D-glucosamine", "Glucosamine_hydrochloride") ~ "Intermediate",
      # Hydrophobic
      Substrate %in% c("Methyl_oleate", "Tween40") ~ "Hydrophobic",
      TRUE ~ "Other"
    ),
    # More specific groupings for enzyme validation
    Substrate_type = case_when(
      Substrate %in% c("L-Leucine", "L-Alanine", "Glycine", "L-Glutamine", 
                       "Glutamic_acid", "L-Arginine") ~ "Amino_acids",
      Substrate %in% c("Citric_acid", "Malic_acid", "Succinic_acid", 
                       "Lactic_acid", "Oxalic_acid") ~ "Organic_acids",
      Substrate %in% c("D-Glucose", "D-Fructose", "D-Trehalose") ~ "Simple_sugars",
      Substrate %in% c("Cellulose", "Hemicellulose", "D-Cellobiose", "D-Xylose") ~ "Cellulose_related",
      Substrate %in% c("Chitin", "Chitosan", "N-Acetyl-D-glucosamine", 
                       "Glucosamine_hydrochloride") ~ "Chitin_related",
      Substrate == "Lignin" ~ "Lignin",
      Substrate == "BSA" ~ "Protein",
      Substrate %in% c("Methyl_oleate", "Tween40") ~ "Lipid",
      TRUE ~ "Other"
    )
  )

# Add functional groups back to microresp data
microresp_analysis <- microresp_analysis %>%
  left_join(substrate_categories, by = c("Substrate", "Substrate_category", 
                                         "Substrate_formula", "Substrate_NOSC"))

# ------------------------------------------------------------------------------
# 1.3 Calculate mean SIR by substrate category and season
# ------------------------------------------------------------------------------

sir_by_category <- microresp_analysis %>%
  filter(!is.na(Season), Functional_group != "Other") %>%
  group_by(Season, Functional_group) %>%
  summarise(
    Mean_SIR = mean(SIR, na.rm = TRUE),
    SE_SIR = sd(SIR, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# Statistical test: Two-way ANOVA (season × substrate category)
sir_for_anova <- microresp_analysis %>%
  filter(!is.na(Season), Functional_group != "Other")

anova_model <- aov(SIR ~ Season * Functional_group, data = sir_for_anova)
anova_results <- summary(anova_model)
print(anova_results)

# Post-hoc tests if interaction is significant
posthoc <- emmeans(anova_model, pairwise ~ Season | Functional_group, adjust = "fdr")
print(posthoc)

# ------------------------------------------------------------------------------
# 1.4 Plot seasonal substrate utilization profiles (Figure for Analysis 1)
# ------------------------------------------------------------------------------

p1 <- ggplot(sir_by_category, aes(x = Functional_group, y = Mean_SIR, fill = Season)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = Mean_SIR - SE_SIR, ymax = Mean_SIR + SE_SIR),
                position = position_dodge(width = 0.8), width = 0.3) +
  scale_fill_manual(values = c("Growing" = "#2E7D32", "Non-Growing" = "#1565C0")) +
  labs(
    title = "Seasonal Substrate Utilization Profiles",
    subtitle = "H1: Elevated Leu (Growing) and CB (Non-Growing) correspond to substrate preferences",
    x = "Substrate Category",
    y = "Mean Substrate-Induced Respiration\n(ng CO₂-C / mg C added / hr)",
    fill = "Season"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.subtitle = element_text(size = 10, color = "gray30")
  )

print(p1)
ggsave("Fig_Q1_Seasonal_Substrate_Profiles.png", p1, width = 10, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 1.5 Enzyme-substrate correlation analysis
# ------------------------------------------------------------------------------

# Join EEA and MicroResp data
# Calculate mean SIR for key substrate types per sample
sir_wide <- microresp_analysis %>%
  group_by(Sample_cat, Season, Substrate_type, Irrigation_trmt_pct) %>%
  summarise(Mean_SIR = mean(SIR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Substrate_type, 
              values_from = Mean_SIR,
              names_prefix = "SIR_")

# Join with enzyme data
enzyme_substrate <- eea_season %>%
  left_join(sir_wide, by = c("Sample_cat", "Season")) %>%
  filter(!is.na(Season))  # Remove samples without season info

# ------------------------------------------------------------------------------
# 1.6 Key correlations: Leu vs Amino Acids; CB vs Cellulose
# ------------------------------------------------------------------------------

# Leucine aminopeptidase vs amino acid respiration
cor_leu_growing <- enzyme_substrate %>%
  filter(Season == "Growing", !is.na(SIR_Amino_acids)) %>%
  with(cor.test(Leu_nmol_gdw_h, SIR_Amino_acids, method = "spearman"))

cor_leu_nongrowing <- enzyme_substrate %>%
  filter(Season == "Non-Growing", !is.na(SIR_Amino_acids)) %>%
  with(cor.test(Leu_nmol_gdw_h, SIR_Amino_acids, method = "spearman"))

# Cellobiohydrolase vs cellulose-related respiration
cor_cb_growing <- enzyme_substrate %>%
  filter(Season == "Growing", !is.na(SIR_Cellulose_related)) %>%
  with(cor.test(CB_nmol_gdw_h, SIR_Cellulose_related, method = "spearman"))

cor_cb_nongrowing <- enzyme_substrate %>%
  filter(Season == "Non-Growing", !is.na(SIR_Cellulose_related)) %>%
  with(cor.test(CB_nmol_gdw_h, SIR_Cellulose_related, method = "spearman"))

# Print correlation results
cat("\n=== VALIDATION OF H1: Enzyme-Substrate Correlations ===\n")
cat("\n--- Leucine Aminopeptidase vs Amino Acid SIR ---\n")
cat("Growing Season: rho =", round(cor_leu_growing$estimate, 3), 
    ", p =", format.pval(cor_leu_growing$p.value, digits = 3), 
    ", n =", cor_leu_growing$parameter + 2, "\n")
cat("Non-Growing Season: rho =", round(cor_leu_nongrowing$estimate, 3), 
    ", p =", format.pval(cor_leu_nongrowing$p.value, digits = 3),
    ", n =", cor_leu_nongrowing$parameter + 2, "\n\n")

cat("--- Cellobiohydrolase vs Cellulose-Related SIR ---\n")
cat("Growing Season: rho =", round(cor_cb_growing$estimate, 3), 
    ", p =", format.pval(cor_cb_growing$p.value, digits = 3),
    ", n =", cor_cb_growing$parameter + 2, "\n")
cat("Non-Growing Season: rho =", round(cor_cb_nongrowing$estimate, 3), 
    ", p =", format.pval(cor_cb_nongrowing$p.value, digits = 3),
    ", n =", cor_cb_nongrowing$parameter + 2, "\n")

# Create correlation summary table
cor_summary <- data.frame(
  Enzyme = c("Leucine Aminopeptidase", "Leucine Aminopeptidase",
             "Cellobiohydrolase", "Cellobiohydrolase"),
  Substrate_Type = c("Amino Acids", "Amino Acids",
                     "Cellulose-Related", "Cellulose-Related"),
  Season = c("Growing", "Non-Growing", "Growing", "Non-Growing"),
  Spearman_rho = c(cor_leu_growing$estimate, cor_leu_nongrowing$estimate,
                   cor_cb_growing$estimate, cor_cb_nongrowing$estimate),
  P_value = c(cor_leu_growing$p.value, cor_leu_nongrowing$p.value,
              cor_cb_growing$p.value, cor_cb_nongrowing$p.value),
  N = c(cor_leu_growing$parameter + 2, cor_leu_nongrowing$parameter + 2,
        cor_cb_growing$parameter + 2, cor_cb_nongrowing$parameter + 2)
)

write_csv(cor_summary, "Table_Enzyme_Substrate_Correlations.csv")
print(cor_summary)

# ------------------------------------------------------------------------------
# 1.7 Scatterplots: Enzyme vs Substrate SIR
# ------------------------------------------------------------------------------

# Leu vs Amino Acids
p2 <- ggplot(enzyme_substrate %>% filter(!is.na(SIR_Amino_acids)), 
             aes(x = Leu_nmol_gdw_h, y = SIR_Amino_acids, color = Season)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = c("Growing" = "#2E7D32", "Non-Growing" = "#1565C0")) +
  facet_wrap(~ Season, scales = "free") +
  labs(
    title = "Leucine Aminopeptidase Activity vs Amino Acid Respiration",
    subtitle = paste0("Growing: ρ=", round(cor_leu_growing$estimate, 2), 
                      ", p=", format.pval(cor_leu_growing$p.value, digits = 2),
                      " | Non-Growing: ρ=", round(cor_leu_nongrowing$estimate, 2),
                      ", p=", format.pval(cor_leu_nongrowing$p.value, digits = 2)),
    x = "Leucine Aminopeptidase (nmol/g/h)",
    y = "Amino Acid SIR\n(ng CO₂-C / mg C / hr)",
    color = "Season"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(size = 9, color = "gray30")
  )

# CB vs Cellulose
p3 <- ggplot(enzyme_substrate %>% filter(!is.na(SIR_Cellulose_related)), 
             aes(x = CB_nmol_gdw_h, y = SIR_Cellulose_related, color = Season)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = c("Growing" = "#2E7D32", "Non-Growing" = "#1565C0")) +
  facet_wrap(~ Season, scales = "free") +
  labs(
    title = "Cellobiohydrolase Activity vs Cellulose-Related Respiration",
    subtitle = paste0("Growing: ρ=", round(cor_cb_growing$estimate, 2), 
                      ", p=", format.pval(cor_cb_growing$p.value, digits = 2),
                      " | Non-Growing: ρ=", round(cor_cb_nongrowing$estimate, 2),
                      ", p=", format.pval(cor_cb_nongrowing$p.value, digits = 2)),
    x = "Cellobiohydrolase (nmol/g/h)",
    y = "Cellulose-Related SIR\n(ng CO₂-C / mg C / hr)",
    color = "Season"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(size = 9, color = "gray30")
  )

# Combined figure
combined_correlations <- p2 / p3
print(combined_correlations)
ggsave("Fig_Q1_Enzyme_Substrate_Correlations.png", combined_correlations, 
       width = 10, height = 8, dpi = 300)

# ==============================================================================
# PRIORITY ANALYSIS 2: Substrate Chemistry & Metabolic Strategies
# ==============================================================================

# ------------------------------------------------------------------------------
# 2.1 NOSC relationship by season
# ------------------------------------------------------------------------------

# Calculate mean SIR per substrate per season
sir_by_substrate_season <- microresp_analysis %>%
  filter(!is.na(Substrate_NOSC), Substrate_NOSC != "na", !is.na(Season)) %>%
  mutate(Substrate_NOSC = as.numeric(Substrate_NOSC)) %>%
  group_by(Season, Substrate, Substrate_NOSC, Functional_group) %>%
  summarise(
    Mean_SIR = mean(SIR, na.rm = TRUE),
    SE_SIR = sd(SIR, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# Linear models: NOSC vs SIR by season
lm_nosc_growing <- sir_by_substrate_season %>%
  filter(Season == "Growing") %>%
  lm(Mean_SIR ~ Substrate_NOSC, data = .)

lm_nosc_nongrowing <- sir_by_substrate_season %>%
  filter(Season == "Non-Growing") %>%
  lm(Mean_SIR ~ Substrate_NOSC, data = .)

# Extract model summaries
summary_growing <- summary(lm_nosc_growing)
summary_nongrowing <- summary(lm_nosc_nongrowing)

cat("\n=== NOSC vs SIR Regression Results ===\n")
cat("Growing Season: R² =", round(summary_growing$r.squared, 3),
    ", slope =", round(coef(lm_nosc_growing)[2], 3),
    ", p =", format.pval(summary_growing$coefficients[2,4], digits = 3), "\n")
cat("Non-Growing Season: R² =", round(summary_nongrowing$r.squared, 3),
    ", slope =", round(coef(lm_nosc_nongrowing)[2], 3),
    ", p =", format.pval(summary_nongrowing$coefficients[2,4], digits = 3), "\n")

# Save regression results
nosc_regression_summary <- data.frame(
  Season = c("Growing", "Non-Growing"),
  R_squared = c(summary_growing$r.squared, summary_nongrowing$r.squared),
  Slope = c(coef(lm_nosc_growing)[2], coef(lm_nosc_nongrowing)[2]),
  Intercept = c(coef(lm_nosc_growing)[1], coef(lm_nosc_nongrowing)[1]),
  P_value = c(summary_growing$coefficients[2,4], summary_nongrowing$coefficients[2,4])
)
write_csv(nosc_regression_summary, "Table_NOSC_SIR_Regression.csv")

# ------------------------------------------------------------------------------
# 2.2 NOSC relationship plot
# ------------------------------------------------------------------------------

p4 <- ggplot(sir_by_substrate_season, 
             aes(x = Substrate_NOSC, y = Mean_SIR, color = Season, shape = Functional_group)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, aes(group = Season), alpha = 0.2) +
  scale_color_manual(values = c("Growing" = "#2E7D32", "Non-Growing" = "#1565C0")) +
  scale_shape_manual(values = c("Labile" = 16, "Polymeric" = 17, 
                                "Intermediate" = 15, "Hydrophobic" = 18)) +
  labs(
    title = "Substrate NOSC vs Substrate-Induced Respiration by Season",
    subtitle = paste0("Growing: R²=", round(summary_growing$r.squared, 2),
                      ", slope=", round(coef(lm_nosc_growing)[2], 1),
                      " | Non-Growing: R²=", round(summary_nongrowing$r.squared, 2),
                      ", slope=", round(coef(lm_nosc_nongrowing)[2], 1)),
    x = "Nominal Oxidation State of Carbon (NOSC)",
    y = "Mean Substrate-Induced Respiration\n(ng CO₂-C / mg C / hr)",
    color = "Season",
    shape = "Substrate Category"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    plot.subtitle = element_text(size = 9, color = "gray30")
  )

print(p4)
ggsave("Fig_Q1_NOSC_vs_SIR.png", p4, width = 10, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 2.3 Substrate category comparisons (detailed)
# ------------------------------------------------------------------------------

# Calculate mean SIR by substrate type (more detailed than functional group)
sir_by_type_season <- microresp_analysis %>%
  filter(!is.na(Season), Substrate_type != "Other") %>%
  group_by(Season, Substrate_type) %>%
  summarise(
    Mean_SIR = mean(SIR, na.rm = TRUE),
    SE_SIR = sd(SIR, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# Two-way ANOVA: Season × Substrate Type
sir_for_anova_detailed <- microresp_analysis %>%
  filter(!is.na(Season), Substrate_type != "Other")

anova_detailed <- aov(SIR ~ Season * Substrate_type, data = sir_for_anova_detailed)
anova_detailed_results <- summary(anova_detailed)
print(anova_detailed_results)

# Post-hoc comparisons
posthoc_detailed <- emmeans(anova_detailed, pairwise ~ Season | Substrate_type, 
                            adjust = "fdr")

# Extract and save pairwise comparisons
posthoc_contrasts <- as.data.frame(posthoc_detailed$contrasts)
write_csv(posthoc_contrasts, "Table_Posthoc_Season_Substrate_Type.csv")

# ------------------------------------------------------------------------------
# 2.4 Detailed substrate type plot
# ------------------------------------------------------------------------------

p5 <- ggplot(sir_by_type_season, aes(x = reorder(Substrate_type, Mean_SIR), 
                                     y = Mean_SIR, fill = Season)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = Mean_SIR - SE_SIR, ymax = Mean_SIR + SE_SIR),
                position = position_dodge(width = 0.8), width = 0.3) +
  scale_fill_manual(values = c("Growing" = "#2E7D32", "Non-Growing" = "#1565C0")) +
  coord_flip() +
  labs(
    title = "Substrate-Specific Respiration Patterns by Season",
    x = "Substrate Type",
    y = "Mean Substrate-Induced Respiration\n(ng CO₂-C / mg C / hr)",
    fill = "Season"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

print(p5)
ggsave("Fig_Q1_Substrate_Type_Details.png", p5, width = 10, height = 7, dpi = 300)

# ------------------------------------------------------------------------------
# 2.5 N-rich vs C-only substrate comparison
# ------------------------------------------------------------------------------

# Categorize substrates by N content
microresp_n_analysis <- microresp_analysis %>%
  filter(!is.na(Season)) %>%
  mutate(
    N_content = case_when(
      Substrate %in% c("L-Leucine", "L-Alanine", "Glycine", "L-Glutamine", 
                       "Glutamic_acid", "L-Arginine", "Chitin", "Chitosan",
                       "N-Acetyl-D-glucosamine", "Glucosamine_hydrochloride", 
                       "BSA", "Phenylethylamine", "GABA", "Diethylamine") ~ "N-containing",
      TRUE ~ "C-only"
    )
  )

# Calculate mean SIR by N content and season
sir_by_n_season <- microresp_n_analysis %>%
  group_by(Season, N_content) %>%
  summarise(
    Mean_SIR = mean(SIR, na.rm = TRUE),
    SE_SIR = sd(SIR, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# Statistical test
anova_n_content <- aov(SIR ~ Season * N_content, data = microresp_n_analysis)
summary(anova_n_content)

# Calculate ratio of N-containing to C-only for each sample
sir_n_ratio <- microresp_n_analysis %>%
  group_by(Sample_cat, Season, N_content) %>%
  summarise(Mean_SIR = mean(SIR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = N_content, values_from = Mean_SIR) %>%
  mutate(N_to_C_ratio = `N-containing` / `C-only`)

# Test if ratio differs by season
wilcox_n_ratio <- wilcox.test(N_to_C_ratio ~ Season, data = sir_n_ratio)

cat("\n=== N-rich vs C-only Substrate Preference ===\n")
cat("Wilcoxon test for N:C ratio difference by season:\n")
cat("W =", wilcox_n_ratio$statistic, ", p =", 
    format.pval(wilcox_n_ratio$p.value, digits = 3), "\n")

# Plot N-containing vs C-only comparison
p6 <- ggplot(sir_by_n_season, aes(x = N_content, y = Mean_SIR, fill = Season)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = Mean_SIR - SE_SIR, ymax = Mean_SIR + SE_SIR),
                position = position_dodge(width = 0.8), width = 0.3) +
  scale_fill_manual(values = c("Growing" = "#2E7D32", "Non-Growing" = "#1565C0")) +
  labs(
    title = "N-Containing vs C-Only Substrate Respiration",
    subtitle = "H1: Growing season microbes favor N-rich substrates during plant competition",
    x = "Substrate N Content",
    y = "Mean Substrate-Induced Respiration\n(ng CO₂-C / mg C / hr)",
    fill = "Season"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    plot.subtitle = element_text(size = 10, color = "gray30")
  )

print(p6)
ggsave("Fig_Q1_N_vs_C_Substrates.png", p6, width = 8, height = 6, dpi = 300)

# ==============================================================================
# PRIORITY ANALYSIS 3: GWC and Irrigation Treatment Effects (Q2)
# ==============================================================================

# ------------------------------------------------------------------------------
# 3.1 Irrigation treatment substrate profiles
# ------------------------------------------------------------------------------

# Calculate mean SIR by substrate and irrigation treatment
sir_by_irrigation <- microresp_analysis %>%
  filter(!is.na(Irrigation_trmt_pct), Substrate_type != "Other") %>%
  group_by(Irrigation_trmt_pct, Substrate_type) %>%
  summarise(
    Mean_SIR = mean(SIR, na.rm = TRUE),
    SE_SIR = sd(SIR, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(Irrigation_label = paste0(Irrigation_trmt_pct, "% Irrigation"))

# Heatmap data preparation
sir_heatmap_data <- microresp_analysis %>%
  filter(!is.na(Irrigation_trmt_pct)) %>%
  group_by(Irrigation_trmt_pct, Substrate) %>%
  summarise(Mean_SIR = mean(SIR, na.rm = TRUE), .groups = "drop") %>%
  mutate(Irrigation_label = paste0(Irrigation_trmt_pct, "% Irrigation"))

# Create heatmap
p7 <- ggplot(sir_heatmap_data, aes(x = Irrigation_label, y = Substrate, fill = Mean_SIR)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#0571b0", mid = "#f7f7f7", high = "#ca0020",
                       midpoint = median(sir_heatmap_data$Mean_SIR, na.rm = TRUE),
                       name = "Mean SIR\n(ng CO₂-C/mg C/hr)") +
  labs(
    title = "Substrate Utilization Profiles Across Irrigation Treatments",
    subtitle = "H2b: Chronic water limitation selects for different substrate preferences",
    x = "Irrigation Treatment",
    y = "Substrate"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.subtitle = element_text(size = 10, color = "gray30")
  )

print(p7)
ggsave("Fig_Q2_Irrigation_Substrate_Heatmap.png", p7, width = 8, height = 12, dpi = 300)

# ------------------------------------------------------------------------------
# 3.2 Polymer/hydrophobic substrate ratios (H2b test)
# ------------------------------------------------------------------------------

# Calculate polymer and labile SIR for each sample
sir_ratios <- microresp_analysis %>%
  filter(!is.na(Season), !is.na(Irrigation_trmt_pct)) %>%
  group_by(Sample_cat, Season, Irrigation_trmt_pct, Functional_group) %>%
  summarise(Mean_SIR = mean(SIR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Functional_group, values_from = Mean_SIR) %>%
  mutate(
    Polymer_to_Labile_ratio = (Polymeric + Hydrophobic) / Labile,
    Irrigation_label = paste0(Irrigation_trmt_pct, "% Irrigation")
  )

# Two-way ANOVA: Season × Irrigation on ratio
anova_ratio <- aov(Polymer_to_Labile_ratio ~ Season * Irrigation_trmt_pct, 
                   data = sir_ratios %>% filter(!is.na(Polymer_to_Labile_ratio)))
summary(anova_ratio)

# Post-hoc
posthoc_ratio <- emmeans(anova_ratio, pairwise ~ Irrigation_trmt_pct | Season, adjust = "fdr")
print(posthoc_ratio)

# Plot ratio
p8 <- ggplot(sir_ratios %>% filter(!is.na(Polymer_to_Labile_ratio)), 
             aes(x = Irrigation_label, y = Polymer_to_Labile_ratio, fill = Season)) +
  geom_boxplot(outlier.alpha = 0.5) +
  scale_fill_manual(values = c("Growing" = "#2E7D32", "Non-Growing" = "#1565C0")) +
  labs(
    title = "Polymer+Hydrophobic : Labile Substrate Ratio",
    subtitle = "H2b: Higher ratio in 25% irrigation = adaptation to diffusion-limited conditions",
    x = "Irrigation Treatment",
    y = "(Polymeric + Hydrophobic) / Labile SIR Ratio",
    fill = "Season"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    plot.subtitle = element_text(size = 10, color = "gray30")
  )

print(p8)
ggsave("Fig_Q2_Polymer_Labile_Ratio.png", p8, width = 10, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 3.3 GWC as continuous predictor of substrate preference
# ------------------------------------------------------------------------------

# Join GWC data to MicroResp
microresp_gwc <- microresp_analysis %>%
  left_join(bgc %>% select(Sample_cat, GWC), by = "Sample_cat") %>%
  filter(!is.na(GWC))

# Key substrates to plot
key_substrates <- c("Methyl_oleate", "Cellulose", "D-Glucose", "L-Leucine")

microresp_gwc_key <- microresp_gwc %>%
  filter(Substrate %in% key_substrates) %>%
  mutate(Irrigation_label = paste0(Irrigation_trmt_pct, "%"))

# GWC vs SIR for key substrates
p9 <- ggplot(microresp_gwc_key, aes(x = GWC, y = SIR, color = Irrigation_label)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
  facet_wrap(~ Substrate, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("25%" = "#d73027", "100%" = "#4575b4"),
                     name = "Irrigation Treatment") +
  labs(
    title = "GWC Effects on Key Substrate Respiration",
    subtitle = "H2a: GWC as continuous moisture predictor shows non-linear threshold effects",
    x = "Gravimetric Water Content (%)",
    y = "Substrate-Induced Respiration\n(ng CO₂-C / mg C / hr)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    plot.subtitle = element_text(size = 10, color = "gray30")
  )

print(p9)
ggsave("Fig_Q2_GWC_Substrate_Relationships.png", p9, width = 10, height = 8, dpi = 300)

# ------------------------------------------------------------------------------
# 3.4 Substrate solubility analysis (H2b mechanism test)
# ------------------------------------------------------------------------------

# Get solubility data and calculate irrigation difference
solubility_analysis <- microresp_analysis %>%
  filter(!is.na(Irrigation_trmt_pct), !is.na(Substrate_solubility),
         Substrate_solubility != "na", Substrate_solubility != "infinite",
         Substrate_solubility != "insoluble") %>%
  mutate(Solubility_numeric = as.numeric(Substrate_solubility)) %>%
  group_by(Substrate, Substrate_NOSC, Solubility_numeric, Irrigation_trmt_pct) %>%
  summarise(Mean_SIR = mean(SIR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Irrigation_trmt_pct, values_from = Mean_SIR,
              names_prefix = "Irr_") %>%
  mutate(SIR_diff_100_minus_25 = Irr_100 - Irr_25)

# Correlation between solubility and irrigation effect
cor_solubility <- cor.test(solubility_analysis$Solubility_numeric, 
                           solubility_analysis$SIR_diff_100_minus_25,
                           method = "spearman", use = "complete.obs")

cat("\n=== Solubility vs Irrigation Effect ===\n")
cat("Spearman correlation: rho =", round(cor_solubility$estimate, 3),
    ", p =", format.pval(cor_solubility$p.value, digits = 3), "\n")
cat("Interpretation: Positive rho = high solubility substrates show greater benefit from 100% irrigation\n")

# Plot solubility analysis
p10 <- ggplot(solubility_analysis, 
              aes(x = log10(Solubility_numeric + 0.1), y = SIR_diff_100_minus_25)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "#ca0020") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Substrate Solubility vs Irrigation Effect on Respiration",
    subtitle = paste0("H2b test: ρ=", round(cor_solubility$estimate, 2),
                      ", p=", format.pval(cor_solubility$p.value, digits = 2)),
    x = "Log10 Substrate Solubility (g/100g water)",
    y = "SIR Difference (100% - 25% Irrigation)\n(ng CO₂-C / mg C / hr)",
    caption = "Positive values = higher respiration under 100% irrigation"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.subtitle = element_text(size = 10, color = "gray30"))

print(p10)
ggsave("Fig_Q2_Solubility_vs_Irrigation_Effect.png", p10, width = 8, height = 6, dpi = 300)

# ==============================================================================
# SUMMARY STATISTICS OUTPUT
# ==============================================================================

# ==============================================================================
# SUMMARY STATISTICS OUTPUT (continued from where it cut off)
# ==============================================================================

# Create comprehensive summary table (line was incomplete)
summary_stats <- list(
  Enzyme_Substrate_Correlations = cor_summary,
  NOSC_Regression = nosc_regression_summary,
  Seasonal_ANOVA = broom::tidy(anova_model),
  Detailed_ANOVA = broom::tidy(anova_detailed),
  N_Content_ANOVA = broom::tidy(anova_n_content),
  Irrigation_Ratio_ANOVA = broom::tidy(anova_ratio)
)

# Save all summary statistics
write_csv(bind_rows(summary_stats, .id = "Analysis"), 
          "Summary_Statistics_All_Analyses.csv")

# ==============================================================================
# SESSION INFO
# ==============================================================================

# Save session info for reproducibility
writeLines(capture.output(sessionInfo()), "Session_Info.txt")

cat("\n=== Analysis Complete ===\n")
cat("All figures saved to working directory\n")
cat("All summary tables saved to working directory\n")