# ==============================================================================
# PUBLICATION-READY THREE-PANEL FIGURE FOR H2B AND H3 ANALYSIS
# ==============================================================================
rm(list = ls())
# Load required libraries
library(tidyverse)
library(vegan)
library(patchwork)
library(pheatmap)
library(ggrepel)
library(gridExtra)
library(ggplotify)
library(lubridate)

# Define colors
season_cols <- c("Growing" = "#E69F00", "Non-Growing" = "#56B4E9")
trt_cols    <- c("100" = "#5AAE61", "25" = "#9970AB")

# Load your data
microresp <- read.csv("DataHub/Soil_Phenotypes/SoilSFA_A3.1_MicroResp.csv")
names(microresp)[1] = "Sample_cat"
biogeochem <- read.csv("Data/RF_final_cleaned_outliers_removed_zeros_replaced.csv")# Define growing season based on your specification


# Define growing season (April-October) and non-growing season (November-March)
define_season <- function(date) {
  month <- month(mdy(date))
  ifelse(month %in% 4:10, "Growing", "Non-Growing")
}

# Prepare MicroResp data (5hr, C-normalized)
microresp_prep <- microresp %>%
  filter(Timepoint == 5,
         Normalization == "ngCO2-C_mgC_h",
         Substrate != "Water") %>%
  mutate(Season = define_season(Collection_date),
         Season = factor(Season, levels = c("Growing", "Non-Growing")),
         Irrigation_trmt_pct = as.character(Irrigation_trmt_pct))

# Verify required columns exist
cat("\n=== Checking Required Columns ===\n")
cat("Sample_cat exists:", "Sample_cat" %in% colnames(microresp_prep), "\n")
cat("Irrigation_trmt_pct exists:", "Irrigation_trmt_pct" %in% colnames(microresp_prep), "\n")
cat("Substrate exists:", "Substrate" %in% colnames(microresp_prep), "\n")
cat("SIR exists:", "SIR" %in% colnames(microresp_prep), "\n\n")

# ==============================================================================
# PANEL A: NMDS PLOTS WITH ENVFIT (SEPARATE BY IRRIGATION)
# ==============================================================================

create_nmds_panel <- function(irrigation_level, data) {
  
  cat(paste0("\n=== Processing ", irrigation_level, "% Irrigation Treatment ===\n"))
  
  # Filter for specific irrigation treatment
  data_subset <- data %>%
    filter(Irrigation_trmt_pct == irrigation_level,
           !is.na(Season))
  
  cat("Number of samples:", length(unique(data_subset$Sample_cat)), "\n")
  
  # Create wide matrix for ordination
  comm_matrix <- data_subset %>%
    select(Sample_cat, Substrate, SIR) %>%
    pivot_wider(names_from = Substrate, 
                values_from = SIR,
                values_fill = 0) %>%
    column_to_rownames("Sample_cat") %>%
    as.matrix()
  
  # Remove columns with all zeros
  comm_matrix <- comm_matrix[, colSums(comm_matrix) > 0]
  
  cat("Number of substrates:", ncol(comm_matrix), "\n")
  
  # Get metadata for these samples
  sample_meta <- data_subset %>%
    distinct(Sample_cat, Season)
  
  # Run NMDS
  set.seed(123)
  nmds <- metaMDS(comm_matrix, distance = "bray", k = 2, trymax = 100)
  
  cat("NMDS converged. Stress:", round(nmds$stress, 4), "\n")
  
  # PERMANOVA
  permanova_result <- adonis2(comm_matrix ~ Season, 
                              data = sample_meta, 
                              permutations = 999,
                              method = "bray")
  
  # Envfit for substrates
  ef <- envfit(nmds, comm_matrix, permutations = 999)
  
  # Extract top 10 substrates by r2
  ef_scores <- scores(ef, display = "vectors")
  ef_df <- as.data.frame(ef_scores) %>%
    rownames_to_column("Substrate") %>%
    mutate(r2 = ef$vectors$r,
           pval = ef$vectors$pvals) %>%
    filter(pval < 0.05) %>%
    arrange(desc(r2)) %>%
    slice_head(n = 10)
  
  cat("Number of significant substrates (p < 0.05):", nrow(ef_df), "\n")
  
  # Scale arrows for better visualization
  arrow_scale <- 0.8
  ef_df <- ef_df %>%
    mutate(NMDS1 = NMDS1 * arrow_scale,
           NMDS2 = NMDS2 * arrow_scale)
  
  # NMDS scores
  nmds_scores <- as.data.frame(scores(nmds, display = "sites")) %>%
    rownames_to_column("Sample_cat") %>%
    left_join(sample_meta, by = "Sample_cat")
  
  # Create plot with FIXED subtitle formatting
  p <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Season)) +
    geom_point(size = 3.5, alpha = 0.8) +
    stat_ellipse(aes(group = Season), 
                 linetype = 2, 
                 linewidth = 0.8,
                 level = 0.95) +
    scale_color_manual(values = season_cols) +
    # Add envfit arrows
    geom_segment(data = ef_df,
                 aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
                 arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
                 color = "grey30", 
                 linewidth = 0.6,
                 alpha = 0.7,
                 inherit.aes = FALSE) +
    geom_text_repel(data = ef_df,
                    aes(x = NMDS1, y = NMDS2, label = Substrate),
                    color = "grey20", 
                    size = 3,
                    fontface = "italic",
                    box.padding = 0.5,
                    max.overlaps = 15,
                    inherit.aes = FALSE) +
    theme_bw() +
    labs(title = paste0(irrigation_level, "% Irrigation"),
         subtitle = paste0("Stress: ", round(nmds$stress, 3),
                           "\nPERMANOVA R² = ", round(permanova_result$R2[1], 3),
                           ", p = ", format.pval(permanova_result$`Pr(>F)`[1], digits = 3)),
         color = "Season") +
    theme(legend.position = "bottom",
          legend.title = element_text(face = "bold", size = 11),
          legend.text = element_text(size = 10),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          panel.grid.minor = element_blank(),
          axis.title = element_text(face = "bold", size = 11))
  
  return(list(plot = p, 
              envfit_df = ef_df,
              permanova = permanova_result,
              stress = nmds$stress))
}

# Create both NMDS plots
nmds_25 <- create_nmds_panel("25", microresp_prep)
nmds_100 <- create_nmds_panel("100", microresp_prep)

# Combine NMDS plots
panel_a <- (nmds_100$plot + nmds_25$plot) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

# ==============================================================================
# PANEL B: HEATMAP OF TOP ENVFIT SUBSTRATES
# ==============================================================================

# Get union of top substrates from both treatments
top_substrates <- unique(c(nmds_25$envfit_df$Substrate, 
                           nmds_100$envfit_df$Substrate))

cat("\n=== Top Substrates Selected for Heatmap ===\n")
cat(paste(top_substrates, collapse = ", "), "\n")
cat("Total number of unique substrates:", length(top_substrates), "\n\n")

# Prepare heatmap data - relative respiration within each sample
heatmap_data <- microresp_prep %>%
  filter(Substrate %in% top_substrates) %>%
  select(Sample_cat, Substrate, SIR, Irrigation_trmt_pct, Season) %>%
  # Calculate relative respiration within each sample
  group_by(Sample_cat) %>%
  mutate(Relative_SIR = 100 * SIR / sum(SIR, na.rm = TRUE)) %>%
  ungroup() %>%
  select(Sample_cat, Substrate, Relative_SIR, Irrigation_trmt_pct, Season) %>%
  pivot_wider(names_from = Substrate, 
              values_from = Relative_SIR,
              values_fill = 0)

# Create matrix for heatmap
heatmap_matrix <- heatmap_data %>%
  select(-Irrigation_trmt_pct, -Season) %>%
  column_to_rownames("Sample_cat") %>%
  as.matrix()

# Annotation for samples
sample_annotation <- heatmap_data %>%
  select(Sample_cat, Season, Irrigation_trmt_pct) %>%
  column_to_rownames("Sample_cat") %>%
  rename(Treatment = Irrigation_trmt_pct) %>%
  mutate(Treatment = factor(Treatment, levels = c("100", "25")))

# Annotation colors
ann_colors <- list(
  Season = season_cols,
  Treatment = trt_cols
)

# Create heatmap
set.seed(123)
heatmap_plot <- pheatmap(
  t(heatmap_matrix),  # Transpose so substrates are rows
  annotation_col = sample_annotation,
  annotation_colors = ann_colors,
  scale = "row",  # Z-score scaling across samples
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  show_colnames = FALSE,
  angle_col = 0,
  fontsize = 9,
  fontsize_row = 8,
  main = "Relative Respiration Patterns\n(Top Substrates)",
  border_color = NA,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  silent = TRUE
)

# Convert to ggplot object
panel_b <- as.ggplot(heatmap_plot)

# ==============================================================================
# PANEL C: POTENTIAL VS REALIZED ACTIVITY RATIO
# ==============================================================================

# Calculate ratios using your pre-calculated field respiration rate
ratio_data <- microresp_prep %>%
  group_by(Sample_cat) %>%
  summarize(Total_MicroResp = sum(SIR, na.rm = TRUE), .groups = "drop") %>%
  left_join(biogeochem %>% 
              select(Sample_cat, Resp_rate_mg_C_m2_h), 
            by = "Sample_cat") %>%
  left_join(microresp_prep %>% 
              distinct(Sample_cat, Season, Irrigation_trmt_pct, Collection_date), 
            by = "Sample_cat") %>%
  filter(!is.na(Resp_rate_mg_C_m2_h), 
         Resp_rate_mg_C_m2_h > 0,
         !is.na(Season)) %>%
  mutate(Ratio = Total_MicroResp / Resp_rate_mg_C_m2_h)

# Check if we have data
cat("\n=== Ratio Data Summary ===\n")
cat("Total samples with ratio data:", nrow(ratio_data), "\n")
print(table(ratio_data$Season))

# Wilcoxon test
wilcox_result <- wilcox.test(Ratio ~ Season, data = ratio_data)
wilcox_pval <- wilcox_result$p.value

# Format p-value
pval_label <- if(wilcox_pval < 0.001) {
  "italic(p) < 0.001"
} else {
  paste0("italic(p) == ", sprintf("%.3f", wilcox_pval))
}

# Calculate summary stats for annotations
ratio_summary <- ratio_data %>%
  group_by(Season) %>%
  summarize(
    median = median(Ratio, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

# Create boxplot
panel_c <- ggplot(ratio_data, aes(x = Season, y = Ratio, fill = Season)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, linewidth = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2.5) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_fill_manual(values = season_cols) +
  annotation_logticks(sides = "l") +
  labs(title = "Potential vs. Realized Activity",
       subtitle = parse(text = pval_label),
       y = expression(bold("Total MicroResp / Field Respiration")),
       x = NULL) +
  annotate("text", x = 1:2, 
           y = ratio_summary$median * 1.5,
           label = paste0("n = ", ratio_summary$n),
           size = 3.5, fontface = "italic") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.title.y = element_text(face = "bold", size = 11),
        axis.title.x = element_text(face = "bold", size = 11),
        axis.text = element_text(size = 10),
        panel.grid.minor = element_blank())

# ==============================================================================
# COMBINE ALL PANELS
# ==============================================================================

# Create final three-panel figure
final_figure <- (panel_a) / 
  (panel_b | panel_c) +
  plot_annotation(tag_levels = 'A',
                  tag_prefix = '(',
                  tag_suffix = ')') &
  theme(plot.tag = element_text(size = 14, face = "bold"))

# Adjust layout heights
final_figure <- final_figure + 
  plot_layout(heights = c(1, 1.2))

# Save figure
ggsave("Figure_H2b_H3_MicroResp_Analysis.pdf", 
       final_figure, 
       width = 14, 
       height = 12, 
       units = "in",
       dpi = 300)

ggsave("Figure_H2b_H3_MicroResp_Analysis.png", 
       final_figure, 
       width = 14, 
       height = 12, 
       units = "in",
       dpi = 300)

# ==============================================================================
# PRINT STATISTICAL RESULTS
# ==============================================================================

cat("\n================================================================================\n")
cat("STATISTICAL RESULTS FOR H2b AND H3 ANALYSIS\n")
cat("================================================================================\n")

cat("\n=== PANEL A: NMDS & PERMANOVA RESULTS ===\n")
cat("\n--- 100% Irrigation Treatment ---\n")
cat("NMDS Stress:", round(nmds_100$stress, 4), "\n\n")
print(nmds_100$permanova)
cat("\nTop 10 Substrates (by R²):\n")
print(nmds_100$envfit_df %>% 
        select(Substrate, r2, pval) %>% 
        mutate(r2 = round(r2, 3), pval = round(pval, 4)))

cat("\n--- 25% Irrigation Treatment ---\n")
cat("NMDS Stress:", round(nmds_25$stress, 4), "\n\n")
print(nmds_25$permanova)
cat("\nTop 10 Substrates (by R²):\n")
print(nmds_25$envfit_df %>% 
        select(Substrate, r2, pval) %>% 
        mutate(r2 = round(r2, 3), pval = round(pval, 4)))

cat("\n=== PANEL C: WILCOXON TEST RESULTS ===\n")
cat("Wilcoxon rank sum test\n\n")
print(wilcox_result)

cat("\n\nMedian Ratios by Season:\n")
print(ratio_summary)

cat("\n=== SAMPLE SIZES ===\n")
cat("\nSamples per Treatment × Season (for NMDS):\n")
print(microresp_prep %>%
        distinct(Sample_cat, Irrigation_trmt_pct, Season) %>%
        count(Irrigation_trmt_pct, Season) %>%
        pivot_wider(names_from = Season, values_from = n))

cat("\nSamples with Ratio Data (Panel C):\n")
print(table(ratio_data$Season))

cat("\n================================================================================\n")

# ==============================================================================
# OPTIONAL: EXPORT DATA FOR SUPPLEMENTAL MATERIALS
# ==============================================================================

# Export PERMANOVA results
permanova_summary <- bind_rows(
  nmds_100$permanova %>% 
    as.data.frame() %>% 
    rownames_to_column("Term") %>%
    mutate(Treatment = "100% Irrigation"),
  nmds_25$permanova %>% 
    as.data.frame() %>% 
    rownames_to_column("Term") %>%
    mutate(Treatment = "25% Irrigation")
)

write_csv(permanova_summary, "Table_S_PERMANOVA_results.csv")

# Export top substrates
top_substrates_summary <- bind_rows(
  nmds_100$envfit_df %>% mutate(Treatment = "100% Irrigation"),
  nmds_25$envfit_df %>% mutate(Treatment = "25% Irrigation")
) %>%
  select(Treatment, Substrate, r2, pval, NMDS1, NMDS2) %>%
  arrange(Treatment, desc(r2))

write_csv(top_substrates_summary, "Table_S_Top_Substrates_Envfit.csv")

# Export ratio data
write_csv(ratio_data, "Data_Potential_Realized_Ratios.csv")

cat("\n=== Supplemental files exported ===\n")
cat("- Table_S_PERMANOVA_results.csv\n")
cat("- Table_S_Top_Substrates_Envfit.csv\n")
cat("- Data_Potential_Realized_Ratios.csv\n")

# ==============================================================================
# ADDITIONAL ANALYSIS: TEST H2b SUBSTRATE COMPLEXITY HYPOTHESIS
# ==============================================================================

cat("\n=== ADDITIONAL ANALYSIS FOR H2b ===\n")
cat("Testing if substrate utilization patterns differ by complexity and treatment\n\n")

# Categorize substrates by complexity using your classification
substrate_complexity <- microresp_prep %>%
  distinct(Substrate, Substrate_category) %>%
  mutate(Complexity = case_when(
    # Complex: polymers/proteins + hydrophobic/aromatic classes
    Substrate_category %in% c("polymer", "protein", "fatty acid ester", 
                              "plant wax ester", "phenolic") ~ "Complex",
    # Labile: small, generally soluble substrates
    Substrate_category %in% c("carbohydrate", "amino acid A", "amino acid B", 
                              "amino acid C", "amino acid D", "carboxylic acids", 
                              "amine", "signaling molecule") ~ "Labile",
    # Handle control category
    Substrate_category == "control" ~ NA_character_,
    TRUE ~ "Labile"
  ))

# Calculate polymer/labile ratio
complexity_analysis <- microresp_prep %>%
  left_join(substrate_complexity, by = c("Substrate", "Substrate_category")) %>%
  filter(Complexity %in% c("Complex", "Labile"), !is.na(Complexity)) %>%
  group_by(Sample_cat, Irrigation_trmt_pct, Season, Complexity) %>%
  summarize(Total_SIR = sum(SIR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Complexity, values_from = Total_SIR, values_fill = 0) %>%
  mutate(Complex_Labile_Ratio = Complex / (Labile + 0.001))  # Add small constant to avoid division by zero

# Test interaction
complexity_model <- aov(Complex_Labile_Ratio ~ Irrigation_trmt_pct * Season, 
                        data = complexity_analysis)

cat("ANOVA: Complex/Labile Ratio ~ Treatment × Season\n")
print(summary(complexity_model))

cat("\n\nMean Complex/Labile Ratios:\n")
print(complexity_analysis %>%
        group_by(Irrigation_trmt_pct, Season) %>%
        summarize(Mean_Ratio = mean(Complex_Labile_Ratio, na.rm = TRUE),
                  SD_Ratio = sd(Complex_Labile_Ratio, na.rm = TRUE),
                  n = n(),
                  .groups = "drop"))

cat("\n================================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n")