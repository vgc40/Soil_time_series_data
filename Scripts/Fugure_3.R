# ================================================================================
# SETUP AND DATA LOADING
# ================================================================================
rm(list = ls())
# Load required libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(vegan)
library(ggplot2)
library(cowplot)
library(corrplot)
library(ggrepel)

# Load your data
microresp <- read.csv("DataHub/Soil_Phenotypes/SoilSFA_A3.1_MicroResp.csv")
names(microresp)[1] = "Sample_cat"

pca_season_lut <- read.csv("Data/RF_cleaned_with_PCA_season.csv") %>%
  mutate(Collection_date = as.Date(Collection_date)) %>%
  select(Collection_date, Plot, Season) %>%  
  distinct()

# Define your color palettes
season_cols <- c("Growing" = "#E69F00", "Non-Growing" = "#56B4E9")
trt_cols    <- c("100" = "#5AAE61", "25" = "#9970AB")

# ================================================================================
# DATA PREPARATION
# ================================================================================

# Filter for 5-hour timepoint and C-normalized data
microresp_filtered <- microresp %>%
  filter(Timepoint == 5,
         Normalization == "ngCO2-C_mgC_h",
         Substrate != "Water") %>%
  mutate(
    Collection_date = as.Date(Collection_date, format = "%m/%d/%y"),
    Month = as.numeric(format(Collection_date, "%m"))
  ) %>%
  left_join(pca_season_lut, by = c("Collection_date", "Plot")) %>%
  mutate(
    Season = if_else(
      is.na(Season) & Month == 7, "Growing",
      if_else(is.na(Season) & Month == 1, "Non-Growing", as.character(Season))
    ),
    Season = factor(Season, levels = c("Growing", "Non-Growing"))
  )

# Create substrate classification using YOUR categorization
substrate_classes <- microresp_filtered %>%
  distinct(Substrate, Substrate_category) %>%
  mutate(
    Complexity = case_when(
      # Complex: polymers/proteins + hydrophobic/aromatic classes (lower bioavailability)
      Substrate_category %in% c("polymer", "protein",
                                "fatty acid ester", "plant wax ester", "phenolic") ~ "Complex",
      # Labile: small, generally soluble substrates (incl. low-MW carboxylic acids)
      Substrate_category %in% c("carbohydrate",
                                "amino acid A", "amino acid B", "amino acid C", "amino acid D",
                                "carboxylic acids", "amine", "signaling molecule") ~ "Labile",
      # Fallback: collapse any remaining categories into Labile (change if desired)
      TRUE ~ "Labile"
    )
  )

# Join complexity back to main data
microresp_filtered <- microresp_filtered %>%
  left_join(substrate_classes, by = c("Substrate", "Substrate_category"))

# ================================================================================
# CALCULATE COMPLEX/LABILE RATIOS
# ================================================================================

# Calculate mean SIR for Complex and Labile substrates per sample
complex_labile_ratios <- microresp_filtered %>%
  filter(Complexity %in% c("Complex", "Labile")) %>%
  group_by(Sample_cat, Season, Plot, Irrigation_trmt_pct, Complexity) %>%
  summarise(Mean_SIR = mean(SIR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    id_cols = c(Sample_cat, Season, Plot, Irrigation_trmt_pct),
    names_from = Complexity,
    values_from = Mean_SIR
  ) %>%
  mutate(
    Complex_Labile_Ratio = Complex / Labile
  )

# Display summary of the calculated ratios
print("Summary of calculated Complex/Labile ratios:")
summary(complex_labile_ratios$Complex_Labile_Ratio)

# Check for any missing values
print(paste("Number of samples with missing ratios:", 
            sum(is.na(complex_labile_ratios$Complex_Labile_Ratio))))

# ================================================================================
# ANALYSIS 2: COMPLEX/LABILE RATIO TEST
# ================================================================================

# Statistical test of complex/labile ratio across seasons
ratio_model <- lmer(Complex_Labile_Ratio ~ Season * Irrigation_trmt_pct + (1|Plot), 
                    data = complex_labile_ratios)

# Summary
#print("=" %R% 50)
print("COMPLEX/LABILE RATIO MODEL SUMMARY")
#print("=" %R% 50)
summary(ratio_model)

# ANOVA table
print("ANOVA Results:")
anova(ratio_model)

# Post-hoc comparisons
emmeans_ratio <- emmeans(ratio_model, pairwise ~ Season)
print("Estimated Marginal Means:")
print(emmeans_ratio)

# Visualization of ratio
p_ratio <- ggplot(complex_labile_ratios, 
                  aes(x = Season, y = Complex_Labile_Ratio, fill = Season)) +
  geom_violin(alpha = 0.5) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  facet_wrap(~Irrigation_trmt_pct, labeller = label_both) +
  scale_fill_manual(values = season_cols) +
  labs(y = "Complex/Labile Ratio",
       x = "Season") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Figures/complex_labile_ratio_season2.pdf", p_ratio, width = 10, height = 6)

# Effect size calculation
ratio_summary <- complex_labile_ratios %>%
  group_by(Season, Irrigation_trmt_pct) %>%
  summarise(
    Mean_Ratio = mean(Complex_Labile_Ratio, na.rm = TRUE),
    SD_Ratio = sd(Complex_Labile_Ratio, na.rm = TRUE),
    SE_Ratio = SD_Ratio / sqrt(n()),
    .groups = "drop"
  )

print("Summary statistics for Complex/Labile Ratio:")
print(ratio_summary)

# Save the calculated ratios for later use
write.csv(complex_labile_ratios, 
          "Results/complex_labile_ratios_calculated.csv", 
          row.names = FALSE)

# ================================================================================
# PANEL A: SUBSTRATE CATEGORY BY SEASON WITH WILCOXON TEST
# ================================================================================

# Calculate mean SIR by complexity category and season
category_summary <- microresp_filtered %>%
  filter(Complexity %in% c("Complex", "Labile")) %>%
  group_by(Sample_cat, Season, Complexity, Irrigation_trmt_pct) %>%
  summarise(Mean_SIR = mean(SIR, na.rm = TRUE), .groups = "drop")

# Perform Wilcoxon test for each complexity category
wilcox_results <- category_summary %>%
  group_by(Complexity) %>%
  summarise(
    p_value = wilcox.test(Mean_SIR ~ Season)$p.value,
    .groups = "drop"
  )

print("Wilcoxon test results:")
print(wilcox_results)

# Create Panel A
panel_a <- ggplot(category_summary, 
                  aes(x = Season, y = Mean_SIR, fill = Season)) +
  geom_violin(alpha = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.3, size = 1) +
  facet_wrap(~ Complexity, scales = "free_y") +
  scale_fill_manual(values = season_cols) +
  labs(
    title = "A",
    y = "Mean SIR (ngCO2-C per mg C per h)",
    x = "Season"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 12),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  # Add p-values from Wilcoxon tests
  geom_text(data = wilcox_results, 
            aes(x = 1.5, y = Inf, 
                label = paste0("p = ", round(p_value, 3))),
            inherit.aes = FALSE,
            vjust = 1.5, size = 3.5)

# ================================================================================
# PANEL B: NMDS WITH ENVFIT
# ================================================================================

# Create matrix for NMDS (samples x substrates)
sir_matrix <- microresp_filtered %>%
  select(Sample_cat, Substrate, SIR) %>%
  pivot_wider(names_from = Substrate, values_from = SIR, values_fill = 0) %>%
  column_to_rownames("Sample_cat") %>%
  as.matrix()

# Get metadata in same order
metadata <- microresp_filtered %>%
  select(Sample_cat, Season, Irrigation_trmt_pct) %>%
  distinct() %>%
  filter(Sample_cat %in% rownames(sir_matrix)) %>%
  arrange(match(Sample_cat, rownames(sir_matrix)))

# Ensure order matches
sir_matrix_ordered <- sir_matrix[metadata$Sample_cat, ]

# Run NMDS
set.seed(123)
nmds <- metaMDS(sir_matrix_ordered, distance = "bray", k = 2, trymax = 100)

# Perform envfit to identify important substrates
envfit_result <- envfit(nmds, sir_matrix_ordered, permutations = 999)

# Extract top substrates based on r-squared
envfit_scores <- as.data.frame(scores(envfit_result, display = "vectors"))
envfit_scores$Substrate <- rownames(envfit_scores)
envfit_scores$r_squared <- envfit_result$vectors$r
envfit_scores$p_value <- envfit_result$vectors$pvals

# Filter for significant substrates and select top 10
top_substrates <- envfit_scores %>%
  filter(p_value < 0.05) %>%
  arrange(desc(r_squared)) %>%
  head(10)

print("Top 10 substrates from envfit:")
print(top_substrates %>% select(Substrate, r_squared, p_value))

# CORRECTED: Extract NMDS scores using nmds$points
nmds_data <- data.frame(
  Sample_cat = rownames(nmds$points),
  NMDS1 = nmds$points[, 1],
  NMDS2 = nmds$points[, 2],
  Season = metadata$Season,
  Irrigation = factor(metadata$Irrigation_trmt_pct)
)

# Scale substrate vectors for plotting
arrow_scale <- 0.5  # Adjust this to make arrows appropriate length
top_substrates_scaled <- top_substrates %>%
  mutate(
    NMDS1_arrow = NMDS1 * arrow_scale,
    NMDS2_arrow = NMDS2 * arrow_scale
  )

# Create Panel B
panel_b <- ggplot(nmds_data, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Season, shape = Irrigation), 
             size = 3, alpha = 0.7) +
  scale_color_manual(values = season_cols) +
  scale_shape_manual(values = c("25" = 17, "100" = 16),
                     labels = c("25% FWC", "100% FWC")) +
  # Add substrate vectors
  geom_segment(data = top_substrates_scaled,
               aes(x = 0, y = 0, xend = NMDS1_arrow, yend = NMDS2_arrow),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "gray30", linewidth = 0.5) +
  geom_text_repel(data = top_substrates_scaled,
                  aes(x = NMDS1_arrow, y = NMDS2_arrow, label = Substrate),
                  size = 3, color = "gray30", 
                  max.overlaps = 20,
                  segment.color = NA) +
  labs(
    title = "B",
    subtitle = paste0("Stress = ", round(nmds$stress, 3), 
                      " | Top 10 substrates (envfit p < 0.05)"),
    color = "Season",
    shape = "Irrigation"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )


# ================================================================================
# PANEL C: HEATMAP OF EFFECT SIZES WITH MIXED EFFECTS MODELS
# ================================================================================

# Run mixed effects models for top substrates from envfit
top_substrate_names <- top_substrates$Substrate

# Initialize results dataframe
lmm_results <- data.frame()

for(sub in top_substrate_names) {
  
  # Filter data for this substrate
  sub_data <- microresp_filtered %>%
    filter(Substrate == sub) %>%
    drop_na(SIR)
  
  # Skip if insufficient data
  if(nrow(sub_data) < 10) {
    message(paste("Skipping", sub, "- insufficient data"))
    next
  }
  
  # Fit linear mixed model
  tryCatch({
    model <- lmer(SIR ~ Season * Irrigation_trmt_pct + (1|Plot), 
                  data = sub_data)
    
    # Extract results
    coef_summary <- summary(model)$coefficients
    
    # Check if SeasonNon-Growing exists
    if("SeasonNon-Growing" %in% rownames(coef_summary)) {
      
      lmm_results <- rbind(lmm_results, data.frame(
        Substrate = sub,
        Season_estimate = coef_summary["SeasonNon-Growing", "Estimate"],
        Season_SE = coef_summary["SeasonNon-Growing", "Std. Error"],
        Season_pvalue = coef_summary["SeasonNon-Growing", "Pr(>|t|)"]
      ))
    }
  }, error = function(e) {
    message(paste("Model failed for", sub, ":", e$message))
  })
}

# Apply FDR correction
lmm_results <- lmm_results %>%
  mutate(
    Season_pvalue_adj = p.adjust(Season_pvalue, method = "fdr"),
    Significant = Season_pvalue_adj < 0.05,
    Significance_marker = ifelse(Significant, "*", "")
  )

print("Mixed effects model results for top substrates:")
print(lmm_results)

# Create heatmap data
heatmap_data <- lmm_results %>%
  mutate(
    Substrate = factor(Substrate, 
                       levels = Substrate[order(Season_estimate)])
  )

# Create Panel C
panel_c <- ggplot(heatmap_data, 
                  aes(x = "Season Effect", y = Substrate)) +
  geom_tile(aes(fill = Season_estimate), 
            color = "white", linewidth = 0.8) +
  geom_text(aes(label = Significance_marker), 
            size = 8, vjust = 0.75, fontface = "bold") +
  scale_fill_gradient2(
    low = season_cols["Growing"], 
    mid = "white", 
    high = season_cols["Non-Growing"],
    midpoint = 0,
    name = "Effect Size\n(Non-Growing\n- Growing)"
  ) +
  labs(
    title = "C",
    subtitle = "* indicates FDR-adjusted p < 0.05",
    y = "Substrate"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0, size = 10),
    legend.position = "right",
    panel.grid = element_blank()
  )

# ================================================================================
# COMBINE PANELS INTO FINAL FIGURE
# ================================================================================

# Combine all three panels
final_figure <- plot_grid(
  panel_a,
  plot_grid(panel_b, panel_c, ncol = 2, rel_widths = c(1.2, 0.8)),
  ncol = 1,
  rel_heights = c(0.8, 1.2)
)

# Save figure
ggsave("Figures/Figure_3_MicroResp_3panel_figure2.pdf", 
       final_figure, 
       width = 14, 
       height = 12,
       units = "in")

ggsave("Figures/Figure_3_MicroResp_3panel_figure2.png", 
       final_figure, 
       width = 14, 
       height = 12,
       units = "in",
       dpi = 300)

print("Figure saved successfully!")

# ================================================================================
# SAVE STATISTICAL RESULTS
# ================================================================================

# Save Wilcoxon results
write.csv(wilcox_results, 
          "Results/wilcoxon_category_season.csv", 
          row.names = FALSE)

# Save envfit results
write.csv(top_substrates, 
          "Results/envfit_top10_substrates.csv", 
          row.names = FALSE)

# Save mixed model results
write.csv(lmm_results, 
          "Results/lmm_seasonal_effects.csv", 
          row.names = FALSE)

print("Statistical results saved!")

# ================================================================================
# OPTIONAL: CREATE INDIVIDUAL PANELS FOR PRESENTATIONS
# ================================================================================

ggsave("Figures/panel_A_categories2.pdf", panel_a, width = 8, height = 6)
ggsave("Figures/panel_B_NMD2S.pdf", panel_b, width = 10, height = 8)
ggsave("Figures/panel_C_heatmap2.pdf", panel_c, width = 6, height = 8)

print("Individual panels saved!")