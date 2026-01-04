# ============================================================================
# MicroResp Analysis for A3.1 Prosser Time Series - UPDATED
# ============================================================================
rm(list = ls())
# Load required packages
library(tidyverse)
library(vegan)
library(indicspecies)
library(lme4)
library(lmerTest)
library(emmeans)
library(pheatmap)
library(viridis)
# Create output directories
dir.create("Figures", showWarnings = FALSE)
dir.create("Results", showWarnings = FALSE)
# Define custom season colors
season_colors <- c("Growing" = "#E69F00", "Non-Growing" = "#56B4E9")
# ============================================================================
# DATA PREPARATION
# ============================================================================
# Load MicroResp data
microresp <- read_csv("DataHub/Soil_Phenotypes/SoilSFA_A3.1_MicroResp.csv")
names(microresp)[1] = 'Sample_cat'
# Load biogeochemistry data
biogeo <- read_csv("Data/RF_final_cleaned_outliers_removed_zeros_replaced.csv")

# Filter to 5-hour timepoint and C-normalized data
microresp_5hr <- microresp %>%
  filter(Timepoint == 5,
         Normalization == "ngCO2-C_mgC_h") %>%
  filter(Substrate != "Water")

# Create season variable - WITH DATE PARSING FIX
microresp_5hr <- microresp_5hr %>%
  mutate(
    Collection_date = lubridate::ymd(Collection_date),  # ADDED FIX
    Month = lubridate::month(Collection_date),
    Season = case_when(
      Month %in% 4:10 ~ "Growing",
      TRUE ~ "Non-Growing"
    )
  )

# Define substrate categories
substrate_categories <- microresp_5hr %>%
  distinct(Substrate, Substrate_category) %>%
  mutate(
    Complexity = case_when(
      Substrate_category %in% c("carbohydrate", "amino acid A", "amino acid B",
                                "amino acid C", "amino acid D", "carboxylic acids") ~ "Labile",
      Substrate_category %in% c("polymer", "fatty acid ester", "plant wax ester") ~ "Complex",
      TRUE ~ "Other"
    )
  )
microresp_5hr <- microresp_5hr %>%
  left_join(substrate_categories, by = c("Substrate", "Substrate_category"))

# ============================================================================
# H1: SEASONAL SUBSTRATE SHIFTS
# ============================================================================
# ----------------------------------------------------------------------------
# 1.2 ENHANCED NMDS Ordination Analysis - ADDRESSING COMMENTS
# ----------------------------------------------------------------------------
# Create sample × substrate matrix
sir_matrix <- microresp_5hr %>%
  select(Sample_cat, Substrate, SIR) %>%
  pivot_wider(
    names_from = Substrate,
    values_from = SIR,
    values_fill = 0
  ) %>%
  column_to_rownames("Sample_cat")

# Run NMDS
set.seed(123)
nmds_result <- metaMDS(
  sir_matrix,
  distance = "bray",
  k = 2,
  trymax = 100
)

# Extract scores and add metadata (INCLUDING BOTH IRRIGATION TREATMENTS)
nmds_scores <- scores(nmds_result, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample_cat") %>%
  left_join(
    microresp_5hr %>% distinct(Sample_cat, Season, Irrigation_trmt_pct, Plot, Collection_date),
    by = "Sample_cat"
  ) %>%
  mutate(
    Month_label = format(Collection_date, "%b-%y")
  )

# Fit environmental vectors
env_fit <- envfit(
  nmds_result,
  sir_matrix,
  permutations = 999
)

# Extract ALL substrate scores with correlation coefficients
all_substrate_scores <- scores(env_fit, display = "vectors") %>%
  as.data.frame() %>%
  rownames_to_column("Substrate") %>%
  mutate(
    r = sqrt(NMDS1^2 + NMDS2^2)
  ) %>%
  arrange(desc(r))

# Print correlation coefficient summary
cat("\n=== Substrate Correlation Coefficient Summary ===\n")
cat("Range of correlation coefficients (r):\n")
cat("  Minimum:", round(min(all_substrate_scores$r), 3), "\n")
cat("  Maximum:", round(max(all_substrate_scores$r), 3), "\n")
cat("  Mean:", round(mean(all_substrate_scores$r), 3), "\n")
cat("  Median:", round(median(all_substrate_scores$r), 3), "\n")
cat("\nTop 10 substrates r range:",
    round(min(all_substrate_scores$r[1:10]), 3), "-",
    round(max(all_substrate_scores$r[1:10]), 3), "\n")
cat("Substrates 11-31 r range:",
    round(min(all_substrate_scores$r[11:31]), 3), "-",
    round(max(all_substrate_scores$r[11:31]), 3), "\n")

# Test for clear distinction (e.g., gap in correlation values)
r_diff_top10_to_11 <- all_substrate_scores$r[10] - all_substrate_scores$r[11]
cat("\nGap between 10th and 11th substrate:", round(r_diff_top10_to_11, 3), "\n")

# Save all correlation coefficients
write_csv(all_substrate_scores, "Results/substrate_envfit_correlations.csv")

# Select top 10 for visualization (scaled to 30% to avoid obscuring points)
substrate_scores_top10 <- all_substrate_scores %>%
  slice(1:10) %>%
  mutate(
    NMDS1_scaled = NMDS1 * 0.3,
    NMDS2_scaled = NMDS2 * 0.3
  )

# NMDS Plot - UPDATED with addressed comments
nmds_plot <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2,
                                     color = Season,
                                     shape = factor(Irrigation_trmt_pct))) +
  geom_point(size = 3, alpha = 0.7) +
  geom_segment(
    data = substrate_scores_top10,
    aes(x = 0, y = 0, xend = NMDS1_scaled, yend = NMDS2_scaled),
    arrow = arrow(length = unit(0.15, "cm")),
    color = "black",
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = substrate_scores_top10,
    aes(x = NMDS1_scaled * 1.15, y = NMDS2_scaled * 1.15, label = Substrate),
    size = 2.5,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = season_colors, name = "Season") +
  scale_shape_manual(values = c("25" = 16, "100" = 17),
                     name = "Irrigation (%)",
                     labels = c("25% FWC", "100% FWC")) +
  labs(
    title = "NMDS of MicroResp Substrate Utilization",
    subtitle = paste0("Stress = ", round(nmds_result$stress, 3),
                      " | Top 10 substrates: r = ",
                      round(min(substrate_scores_top10$r), 3), "-",
                      round(max(substrate_scores_top10$r), 3))
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )
print(nmds_plot)

# Save NMDS plot as PDF and PNG
ggsave("Figures/NMDS_MicroResp_substrates.pdf", nmds_plot, width = 10, height = 7)
ggsave("Figures/NMDS_MicroResp_substrates.png", nmds_plot, width = 10, height = 7, dpi = 300)

# PERMANOVA test (both irrigation treatments included)
metadata <- microresp_5hr %>%
  distinct(Sample_cat, Season, Irrigation_trmt_pct, Plot)
sir_matrix_ordered <- sir_matrix[metadata$Sample_cat, ]
permanova_result <- adonis2(
  sir_matrix_ordered ~ Season * Irrigation_trmt_pct,
  data = metadata,
  permutations = 999,
  method = "bray"
)
permanova_result

adonis2(sir_matrix_ordered ~ Season * Irrigation_trmt_pct,
        data=metadata, permutations=999, method="bray", by="term")

adonis2(sir_matrix_ordered ~ Season + Irrigation_trmt_pct,
        data=metadata, permutations=999, method="bray", by="term")

bd <- betadisper(vegdist(sir_matrix, method="bray"), metadata$Irrigation_trmt_pct)
permutest(bd, permutations=999)


# Save PERMANOVA results
capture.output(permanova_result, file = "Results/permanova_results.txt")

# ============================================================================
# H1: Additional Analyses
# ============================================================================
# 1.1 Substrate Category Comparison by Season
category_means <- microresp_5hr %>%
  filter(Complexity %in% c("Labile", "Complex")) %>%
  group_by(Sample_cat, Season, Plot, Irrigation_trmt_pct, Complexity) %>%
  summarise(Mean_SIR = mean(SIR, na.rm = TRUE), .groups = "drop")

model_h1 <- lmer(
  log1p(Mean_SIR) ~ Season * Complexity + (1|Plot),
  data = category_means
)
summary(model_h1)
anova(model_h1)
emmeans_h1 <- emmeans(model_h1, pairwise ~ Season | Complexity)

# Visualization
p_category <- ggplot(category_means, aes(x = Season, y = Mean_SIR, fill = Season)) +
  geom_boxplot() +
  facet_wrap(~Complexity) +
  scale_fill_manual(values = season_colors, name = "Season") +
  labs(
    title = "Substrate-Induced Respiration by Season and Complexity",
    y = "Mean SIR (ngCO2-C per mg C per h)",
    x = "Season"
  ) +
  theme_bw()
print(p_category)
ggsave("Figures/H1_substrate_category_by_season.pdf", p_category, width = 8, height = 5)
ggsave("Figures/H1_substrate_category_by_season.png", p_category, width = 8, height = 5, dpi = 300)

# 1.4 Substrate Diversity Indices
diversity_indices <- microresp_5hr %>%
  group_by(Sample_cat) %>%
  summarise(
    Shannon = diversity(SIR, index = "shannon"),
    Simpson = diversity(SIR, index = "simpson"),
    .groups = "drop"
  ) %>%
  left_join(
    microresp_5hr %>% distinct(Sample_cat, Season, Irrigation_trmt_pct, Plot),
    by = "Sample_cat"
  )

model_shannon <- lmer(Shannon ~ Season + (1|Plot), data = diversity_indices)
model_simpson <- lmer(Simpson ~ Season + (1|Plot), data = diversity_indices)

# Visualization
p_diversity <- diversity_indices %>%
  pivot_longer(cols = c(Shannon, Simpson), names_to = "Index", values_to = "Value") %>%
  ggplot(aes(x = Season, y = Value, fill = Season)) +
  geom_boxplot() +
  facet_wrap(~Index, scales = "free_y") +
  scale_fill_manual(values = season_colors, name = "Season") +
  theme_bw()
print(p_diversity)
ggsave("Figures/H1_substrate_diversity.pdf", p_diversity, width = 8, height = 5)
ggsave("Figures/H1_substrate_diversity.png", p_diversity, width = 8, height = 5, dpi = 300)

# ============================================================================
# H2b: MOISTURE-DRIVEN METABOLIC STRATEGY SHIFTS
# ============================================================================
# Calculate complex:labile ratio
complexity_ratio <- microresp_5hr %>%
  filter(Complexity %in% c("Labile", "Complex")) %>%
  group_by(Sample_cat, Season, Plot, Irrigation_trmt_pct, Complexity) %>%
  summarise(Mean_SIR = mean(SIR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Complexity, values_from = Mean_SIR) %>%
  mutate(Complex_Labile_Ratio = Complex / Labile)

model_h2b <- lmer(
  log1p(Complex_Labile_Ratio) ~ factor(Irrigation_trmt_pct) * Season + (1|Plot),
  data = complexity_ratio
)
summary(model_h2b)

# Visualization
p_complexity_ratio <- ggplot(complexity_ratio,
                             aes(x = factor(Irrigation_trmt_pct),
                                 y = Complex_Labile_Ratio,
                                 fill = Season)) +
  geom_boxplot() +
  scale_fill_manual(values = season_colors, name = "Season") +
  labs(
    title = "Complex:Labile Substrate Ratio by Irrigation and Season",
    x = "Irrigation Treatment (% FWC)",
    y = "Complex:Labile Ratio"
  ) +
  theme_bw()
print(p_complexity_ratio)
ggsave("Figures/H2b_complexity_ratio.pdf", p_complexity_ratio, width = 8, height = 5)
ggsave("Figures/H2b_complexity_ratio.png", p_complexity_ratio, width = 8, height = 5, dpi = 300)

# ============================================================================
# H3: WINTER METABOLIC POTENTIAL
# ============================================================================
# Merge field respiration data
microresp_resp <- microresp_5hr %>%
  left_join(
    biogeo %>% select(Sample_cat, Resp_rate_mg_C_m2_h),
    by = "Sample_cat"
  ) %>%
  filter(!is.na(Resp_rate_mg_C_m2_h))

threshold <- 10
microresp_resp <- microresp_resp %>%
  mutate(
    Respiration_Level = if_else(
      Resp_rate_mg_C_m2_h < threshold,
      "Low",
      "Moderate_High"
    )
  )

# Calculate total MicroResp per sample
total_microresp <- microresp_resp %>%
  group_by(Sample_cat, Season, Respiration_Level, Resp_rate_mg_C_m2_h, Plot) %>%
  summarise(Total_MicroResp = sum(SIR, na.rm = TRUE), .groups = "drop")

# Correlation plot
p_field_vs_microresp <- ggplot(total_microresp,
                               aes(x = Resp_rate_mg_C_m2_h,
                                   y = Total_MicroResp,
                                   color = Season)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = season_colors, name = "Season") +
  labs(
    title = "Field Respiration vs. Total MicroResp Capacity",
    x = "Field Respiration (mg C m⁻² h⁻¹)",
    y = "Total MicroResp SIR (sum across all substrates)"
  ) +
  theme_bw()
print(p_field_vs_microresp)
ggsave("Figures/H3_field_vs_microresp.pdf", p_field_vs_microresp, width = 8, height = 6)
ggsave("Figures/H3_field_vs_microresp.png", p_field_vs_microresp, width = 8, height = 6, dpi = 300)

# Potential:Realized ratio
ratio_data <- total_microresp %>%
  mutate(Potential_Realized_Ratio = Total_MicroResp / Resp_rate_mg_C_m2_h)

model_ratio <- lmer(
  log1p(Potential_Realized_Ratio) ~ Season + (1|Plot),
  data = ratio_data
)
summary(model_ratio)

p_ratio <- ggplot(ratio_data,
                  aes(x = Season, y = Potential_Realized_Ratio, fill = Season)) +
  geom_boxplot() +
  scale_fill_manual(values = season_colors, name = "Season") +
  scale_y_log10() +
  labs(
    title = "Potential vs. Realized Activity Ratio",
    y = "Total MicroResp / Field Respiration (log scale)",
    subtitle = "Higher ratios indicate greater metabolic potential relative to in situ activity"
  ) +
  theme_bw()
print(p_ratio)
ggsave("Figures/H3_potential_realized_ratio.pdf", p_ratio, width = 8, height = 6)
ggsave("Figures/H3_potential_realized_ratio.png", p_ratio, width = 8, height = 6, dpi = 300)

# ============================================================================
# SAVE ALL DATA OUTPUTS
# ============================================================================
write_csv(all_substrate_scores, "Results/substrate_envfit_all_correlations.csv")
write_csv(diversity_indices, "Results/substrate_diversity_indices.csv")
write_csv(complexity_ratio, "Results/complex_labile_ratios.csv")
write_csv(ratio_data, "Results/microresp_field_ratio.csv")

# Save model summaries
capture.output(
  summary(model_h1),
  summary(model_h2b),
  summary(model_ratio),
  file = "Results/model_summaries.txt"
)

cat("\n=== Analysis Complete ===\n")
cat("All figures saved to Figures/ as PDF and PNG\n")
cat("All results saved to Results/\n")