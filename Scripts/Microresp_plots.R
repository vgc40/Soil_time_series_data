# ============================================================================
# MicroResp Analysis for A3.1 Prosser Time Series
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

# ============================================================================
# DATA PREPARATION
# ============================================================================

# Load MicroResp data
microresp <- read_csv("DataHub/Soil_Phenotypes/SoilSFA_A3.1_MicroResp.csv")
names(microresp)[1] = 'Sample_cat'  # Rename first column to Sample_cat

# Load biogeochemistry data for field respiration rates
biogeo <- read_csv("Data/RF_final_cleaned_outliers_removed_zeros_replaced.csv")

# Filter to 5-hour timepoint and C-normalized data
microresp_5hr <- microresp %>%
  filter(Timepoint == 5,
         Normalization == "ngCO2-C_mgC_h") %>%
  # Remove water control
  filter(Substrate != "Water")

# Create season variable
microresp_5hr <- microresp_5hr %>%
  mutate(
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
# 1.1 Substrate Category Comparison by Season
# ----------------------------------------------------------------------------

# Calculate mean SIR per category per sample
category_means <- microresp_5hr %>%
  filter(Complexity %in% c("Labile", "Complex")) %>%
  group_by(Sample_cat, Season, Plot, Irrigation_trmt_pct, Complexity) %>%
  summarise(
    Mean_SIR = mean(SIR, na.rm = TRUE),
    .groups = "drop"
  )

# Mixed-effects model: SIR ~ Season * Complexity + (1|Plot)
model_h1 <- lmer(
  log1p(Mean_SIR) ~ Season * Complexity + (1|Plot),
  data = category_means
)

summary(model_h1)
anova(model_h1)

# Post-hoc comparisons
emmeans_h1 <- emmeans(model_h1, pairwise ~ Season | Complexity)
emmeans_h1

# Visualization
ggplot(category_means, aes(x = Season, y = Mean_SIR, fill = Complexity)) +
  geom_boxplot() +
  scale_fill_viridis_d() +
  labs(
    title = "Substrate-Induced Respiration by Season and Complexity",
    y = "Mean SIR (ngCO2-C per mg C per h)",
    x = "Season"
  ) +
  theme_bw()

# ----------------------------------------------------------------------------
# 1.2 Ordination Analysis (NMDS)
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

# Extract scores and add metadata
nmds_scores <- scores(nmds_result, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample_cat") %>%
  left_join(
    microresp_5hr %>% distinct(Sample_cat, Season, Irrigation_trmt_pct, Plot),
    by = "Sample_cat"
  )

# PERMANOVA test
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

# Fit environmental vectors
env_fit <- envfit(
  nmds_result,
  sir_matrix,
  permutations = 999
)

# Extract top 10 substrates by correlation
substrate_scores <- scores(env_fit, display = "vectors") %>%
  as.data.frame() %>%
  rownames_to_column("Substrate") %>%
  mutate(
    r = sqrt(NMDS1^2 + NMDS2^2)
  ) %>%
  arrange(desc(r)) %>%
  slice(1:10)

# Visualization
ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Season, shape = factor(Irrigation_trmt_pct))) +
  geom_point(size = 3, alpha = 0.7) +
  geom_segment(
    data = substrate_scores,
    aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "black",
    inherit.aes = FALSE
  ) +
  geom_text(
    data = substrate_scores,
    aes(x = NMDS1, y = NMDS2, label = Substrate),
    size = 3,
    inherit.aes = FALSE
  ) +
  scale_color_viridis_d() +
  labs(
    title = "NMDS of MicroResp Substrate Utilization",
    subtitle = paste0("Stress = ", round(nmds_result$stress, 3))
  ) +
  theme_bw()

# ----------------------------------------------------------------------------
# 1.3 Indicator Substrate Analysis
# ----------------------------------------------------------------------------

# Prepare data for indicator species analysis
sir_matrix_indicator <- sir_matrix[metadata$Sample_cat, ]

# Run indicator analysis
indicator_result <- multipatt(
  sir_matrix_indicator,
  metadata$Season,
  control = how(nperm = 999)
)

summary(indicator_result)

# Extract significant indicators
indicator_summary <- indicator_result$sign %>%
  rownames_to_column("Substrate") %>%
  filter(p.value < 0.05) %>%
  arrange(p.value)

print(indicator_summary)

# ----------------------------------------------------------------------------
# 1.4 Substrate Diversity Indices
# ----------------------------------------------------------------------------

# Calculate Shannon and Simpson diversity per sample
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

# Test seasonal differences
model_shannon <- lmer(
  Shannon ~ Season + (1|Plot),
  data = diversity_indices
)

model_simpson <- lmer(
  Simpson ~ Season + (1|Plot),
  data = diversity_indices
)

summary(model_shannon)
summary(model_simpson)

# Visualization
diversity_indices %>%
  pivot_longer(cols = c(Shannon, Simpson), names_to = "Index", values_to = "Value") %>%
  ggplot(aes(x = Season, y = Value, fill = Season)) +
  geom_boxplot() +
  facet_wrap(~Index, scales = "free_y") +
  scale_fill_viridis_d() +
  theme_bw()

# ============================================================================
# H2b: MOISTURE-DRIVEN METABOLIC STRATEGY SHIFTS
# ============================================================================

# ----------------------------------------------------------------------------
# 2.1 Irrigation Treatment Comparison
# ----------------------------------------------------------------------------

# Calculate complex:labile ratio per sample
complexity_ratio <- microresp_5hr %>%
  filter(Complexity %in% c("Labile", "Complex")) %>%
  group_by(Sample_cat, Season, Plot, Irrigation_trmt_pct, Complexity) %>%
  summarise(Mean_SIR = mean(SIR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Complexity, values_from = Mean_SIR) %>%
  mutate(
    Complex_Labile_Ratio = Complex / Labile
  )

# Mixed-effects model
model_h2b <- lmer(
  log1p(Complex_Labile_Ratio) ~ factor(Irrigation_trmt_pct) * Season + (1|Plot),
  data = complexity_ratio
)

summary(model_h2b)
anova(model_h2b)

# Post-hoc comparisons
emmeans_h2b <- emmeans(model_h2b, pairwise ~ Irrigation_trmt_pct | Season)
emmeans_h2b

# ----------------------------------------------------------------------------
# 2.2 Within-Season Moisture Effects
# ----------------------------------------------------------------------------

# Merge with GWC data
microresp_gwc <- microresp_5hr %>%
  left_join(
    biogeo %>% select(Sample_cat, GWC),
    by = "Sample_cat"
  ) %>%
  filter(!is.na(GWC))

# For each season, test GWC × Complexity interaction
model_growing <- microresp_gwc %>%
  filter(Season == "Growing", Complexity %in% c("Labile", "Complex")) %>%
  lmer(
    log1p(SIR) ~ GWC * Complexity + factor(Irrigation_trmt_pct) + (1|Plot),
    data = .
  )

model_nongrowing <- microresp_gwc %>%
  filter(Season == "Non-Growing", Complexity %in% c("Labile", "Complex")) %>%
  lmer(
    log1p(SIR) ~ GWC * Complexity + factor(Irrigation_trmt_pct) + (1|Plot),
    data = .
  )

summary(model_growing)
summary(model_nongrowing)

# ----------------------------------------------------------------------------
# 2.3 Substrate-Specific Moisture Responses
# ----------------------------------------------------------------------------

# Select key substrates (examples from indicator analysis or NMDS)
key_substrates <- c("D-Glucose", "Citric_acid", "Cellulose", "Methyl_oleate")

substrate_moisture_models <- lapply(key_substrates, function(sub) {
  data_sub <- microresp_gwc %>%
    filter(Substrate == sub)
  
  model <- lmer(
    log1p(SIR) ~ GWC * Season + factor(Irrigation_trmt_pct) + (1|Plot),
    data = data_sub
  )
  
  list(
    substrate = sub,
    model = model,
    summary = summary(model)
  )
})

names(substrate_moisture_models) <- key_substrates

# Print summaries
lapply(substrate_moisture_models, function(x) x$summary)

# ----------------------------------------------------------------------------
# 2.4 Heatmap Visualization
# ----------------------------------------------------------------------------

# Calculate mean SIR for each substrate × irrigation × season combination
heatmap_data <- microresp_5hr %>%
  group_by(Substrate, Irrigation_trmt_pct, Season) %>%
  summarise(Mean_SIR = mean(SIR, na.rm = TRUE), .groups = "drop") %>%
  unite("Condition", Irrigation_trmt_pct, Season, sep = "_") %>%
  pivot_wider(
    names_from = Condition,
    values_from = Mean_SIR
  ) %>%
  column_to_rownames("Substrate")

# Create heatmap
pheatmap(
  log1p(heatmap_data),
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  color = viridis(100),
  main = "Substrate Utilization Patterns Across Conditions",
  fontsize_row = 8
)

# ============================================================================
# H3: WINTER METABOLIC POTENTIAL
# ============================================================================

# ----------------------------------------------------------------------------
# 3.1 Define Low Respiration Threshold
# ----------------------------------------------------------------------------

# Merge field respiration data
microresp_resp <- microresp_5hr %>%
  left_join(
    biogeo %>% 
      select(Sample_cat, Resp_rate_mg_C_m2_h),
    by = "Sample_cat"
  ) %>%
  filter(!is.na(Resp_rate_mg_C_m2_h))

# Define threshold (< 10 mg C m⁻² h⁻¹)
threshold <- 10

# Classify samples
microresp_resp <- microresp_resp %>%
  mutate(
    Respiration_Level = if_else(
      Resp_rate_mg_C_m2_h < threshold,
      "Low",
      "Moderate_High"
    )
  )

# ----------------------------------------------------------------------------
# 3.2 Test for Maintained Capacity
# ----------------------------------------------------------------------------

# Calculate total MicroResp activity per sample
total_microresp <- microresp_resp %>%
  group_by(Sample_cat, Season, Respiration_Level, Resp_rate_mg_C_m2_h) %>%
  summarise(
    Total_MicroResp = sum(SIR, na.rm = TRUE),
    .groups = "drop"
  )

# Compare total MicroResp between low and moderate/high field respiration
# (within non-growing season)
nongrowing_capacity <- total_microresp %>%
  filter(Season == "Non-Growing")

t.test(
  Total_MicroResp ~ Respiration_Level,
  data = nongrowing_capacity
)

# Correlation between field respiration and MicroResp
cor.test(
  total_microresp$Resp_rate_mg_C_m2_h,
  total_microresp$Total_MicroResp
)

# Season-specific correlations
total_microresp %>%
  group_by(Season) %>%
  summarise(
    cor = cor(Resp_rate_mg_C_m2_h, Total_MicroResp, use = "complete.obs"),
    p_value = cor.test(Resp_rate_mg_C_m2_h, Total_MicroResp)$p.value
  )

# Visualization
ggplot(total_microresp, aes(x = Resp_rate_mg_C_m2_h, y = Total_MicroResp, color = Season)) +
  geom_point(size = 3, alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_viridis_d() +
  labs(
    title = "Field Respiration vs. Total MicroResp Capacity",
    x = "Field Respiration (mg C m⁻² h⁻¹)",
    y = "Total MicroResp SIR (sum across all substrates)"
  ) +
  theme_bw()

# ----------------------------------------------------------------------------
# 3.3 Substrate Breadth Analysis
# ----------------------------------------------------------------------------

# Calculate number of substrates with detectable respiration per sample
# (Define detection threshold, e.g., > 0 or > 10th percentile)
detection_threshold <- quantile(microresp_resp$SIR, 0.10, na.rm = TRUE)

substrate_breadth <- microresp_resp %>%
  group_by(Sample_cat, Season, Respiration_Level) %>%
  summarise(
    N_Substrates_Detected = sum(SIR > detection_threshold, na.rm = TRUE),
    .groups = "drop"
  )

# Compare substrate breadth: low winter vs. high summer
breadth_comparison <- substrate_breadth %>%
  filter(
    (Season == "Non-Growing" & Respiration_Level == "Low") |
      (Season == "Growing" & Respiration_Level == "Moderate_High")
  )

t.test(
  N_Substrates_Detected ~ Season,
  data = breadth_comparison
)

# Visualization
ggplot(substrate_breadth, aes(x = Respiration_Level, y = N_Substrates_Detected, fill = Season)) +
  geom_boxplot() +
  scale_fill_viridis_d() +
  labs(
    title = "Metabolic Substrate Breadth",
    x = "Field Respiration Level",
    y = "Number of Substrates with Detectable SIR"
  ) +
  theme_bw()

# ----------------------------------------------------------------------------
# 3.4 MicroResp:Field Respiration Ratio
# ----------------------------------------------------------------------------

# Calculate ratio
ratio_data <- total_microresp %>%
  mutate(
    Potential_Realized_Ratio = Total_MicroResp / Resp_rate_mg_C_m2_h
  )

# Test seasonal differences in ratio
model_ratio <- lmer(
  log1p(Potential_Realized_Ratio) ~ Season + (1|Plot),
  data = ratio_data %>%
    left_join(
      microresp_5hr %>% distinct(Sample_cat, Plot),
      by = "Sample_cat"
    )
)

summary(model_ratio)

# Visualization
ggplot(ratio_data, aes(x = Season, y = Potential_Realized_Ratio, fill = Season)) +
  geom_boxplot() +
  scale_fill_viridis_d() +
  scale_y_log10() +
  labs(
    title = "Potential vs. Realized Activity Ratio",
    y = "Total MicroResp / Field Respiration (log scale)",
    subtitle = "Higher ratios indicate greater metabolic potential relative to in situ activity"
  ) +
  theme_bw()

# ============================================================================
# EXPORT RESULTS
# ============================================================================

# Save key results tables
write_csv(indicator_summary, "indicator_substrates_by_season.csv")
write_csv(diversity_indices, "substrate_diversity_indices.csv")
write_csv(complexity_ratio, "complex_labile_ratios.csv")
write_csv(substrate_breadth, "substrate_breadth_analysis.csv")
write_csv(ratio_data, "microresp_field_ratio.csv")

# Save PERMANOVA results
capture.output(permanova_result, file = "permanova_results.txt")

# Save model summaries
capture.output(
  summary(model_h1),
  summary(model_h2b),
  summary(model_ratio),
  file = "model_summaries.txt"
)
