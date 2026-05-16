# ==============================================================================
# PUBLICATION-READY FIGURE FOR H2B AND H3 ANALYSIS
# Layout (requested):
#   Row 1:  A = NMDS (100% irrigation)   |   B = NMDS (25% irrigation)
#   Row 2:  C = Heatmap                  |   D = Ratio boxplot
#
# Key design choices:
# - Do NOT force square aspect ratios (no coord_equal / coord_fixed)
# - Use ggplotify::as.ggplot(heatmap_plot$gtable) so the heatmap scales in patchwork
# - Give the heatmap more width than the boxplot
# - Clean borders and consistent margins
# - Automatic tags A–D
# - Save in Letter landscape for “sheet of paper” viewing
# ==============================================================================

rm(list = ls())

# ---- Libraries ----
library(tidyverse)
library(vegan)
library(patchwork)
library(pheatmap)
library(ggrepel)
library(ggplotify)
library(lubridate)
library(grid)

# ---- Colors ----
season_cols  <- c("Growing" = "#E69F00", "Non-Growing" = "#56B4E9")
trt_cols     <- c("100" = "#5AAE61", "25" = "#9970AB")
complex_cols <- c("Complex" = "#D95F02", "Labile" = "#1B9E77")

# ---- Data ----
microresp <- read.csv("DataHub/Soil_Phenotypes/SoilSFA_A3.1_MicroResp.csv")
names(microresp)[1] <- "Sample_cat"
#biogeochem <- read.csv("Data/RF_final_cleaned_outliers_removed_zeros_replaced.csv")
biogeochem <- read.csv("Data/RF_cleaned_with_PCA_season.csv")


pca_season_lut <- biogeochem %>%
  select(Sample_cat, Season) %>%
  distinct()

# ---- Prep microresp + join PCA season ----
microresp_prep <- microresp %>%
  filter(Timepoint == 5,
         Normalization == "ngCO2-C_mgC_h",
         Substrate != "Water") %>%
  mutate(
    Collection_date = as.Date(Collection_date, format = "%m/%d/%y"),
    Month = as.integer(format(Collection_date, "%m")),
    Irrigation_trmt_pct = as.character(Irrigation_trmt_pct)
  ) %>%
  left_join(pca_season_lut, by = "Sample_cat") %>%
  mutate(
    Season = if_else(
      is.na(Season) & Month == 7, "Growing",
      if_else(is.na(Season) & Month == 1, "Non-Growing", as.character(Season))
    ),
    Season = factor(Season, levels = c("Growing", "Non-Growing"))
  )


cat("\n=== Data Preparation Complete ===\n")
cat("Total rows:", nrow(microresp_prep), "\n")
cat("Unique samples:", length(unique(microresp_prep$Sample_cat)), "\n\n")

# ==============================================================================
# PANEL A/B: NMDS PLOTS WITH ENVFIT (SEPARATE BY IRRIGATION)
# ==============================================================================
create_nmds_panel <- function(irrigation_level, data) {
  
  cat(paste0("Processing ", irrigation_level, "% Irrigation Treatment\n"))
  
  data_subset <- data %>%
    filter(Irrigation_trmt_pct == irrigation_level,
           !is.na(Season))
  
  comm_matrix <- data_subset %>%
    select(Sample_cat, Substrate, SIR) %>%
    pivot_wider(names_from = Substrate, values_from = SIR, values_fill = 0) %>%
    column_to_rownames("Sample_cat") %>%
    as.matrix()
  
  comm_matrix <- comm_matrix[, colSums(comm_matrix) > 0, drop = FALSE]
  
  sample_meta <- data_subset %>%
    distinct(Sample_cat, Season)
  
  set.seed(123)
  nmds <- metaMDS(comm_matrix, distance = "bray", k = 2, trymax = 100)
  
  permanova_result <- adonis2(
    comm_matrix ~ Season,
    data = sample_meta,
    permutations = 999,
    method = "bray"
  )
  
  ef <- envfit(nmds, comm_matrix, permutations = 999)
  
  ef_scores <- scores(ef, display = "vectors")
  ef_df <- as.data.frame(ef_scores) %>%
    rownames_to_column("Substrate") %>%
    mutate(r2 = ef$vectors$r,
           pval = ef$vectors$pvals) %>%
    filter(pval < 0.05) %>%
    arrange(desc(r2)) %>%
    slice_head(n = 8)
  
  arrow_scale <- 0.8
  ef_df <- ef_df %>%
    mutate(NMDS1_arrow = NMDS1 * arrow_scale,
           NMDS2_arrow = NMDS2 * arrow_scale)
  
  nmds_scores <- as.data.frame(scores(nmds, display = "sites")) %>%
    rownames_to_column("Sample_cat") %>%
    left_join(sample_meta, by = "Sample_cat")
  
  # pull p-value for Season term
  perm_p <- permanova_result$`Pr(>F)`[1]
  perm_p_lab <- ifelse(is.na(perm_p), "NA",
                       ifelse(perm_p < 0.001, "<0.001", sprintf("%.3f", perm_p)))
  
  p <- ggplot(nmds_scores, aes(NMDS1, NMDS2)) +
    geom_point(aes(color = Season), size = 2.7, alpha = 0.8) +
    stat_ellipse(aes(group = Season, color = Season),
                 linetype = 2, linewidth = 0.7, level = 0.95) +
    scale_color_manual(values = season_cols) +
    geom_segment(
      data = ef_df,
      aes(x = 0, y = 0, xend = NMDS1_arrow, yend = NMDS2_arrow),
      arrow = arrow(length = unit(0.18, "cm")),
      color = "gray30",
      linewidth = 0.45
    ) +
    geom_text_repel(
      data = ef_df,
      aes(x = NMDS1_arrow, y = NMDS2_arrow, label = Substrate),
      color = "gray30",
      size = 2.8,
      max.overlaps = 30,
      segment.color = NA
    ) +
    theme_bw() +
    labs(
      title = paste0(irrigation_level, "% irrigation"),
      subtitle = paste0("Stress = ", round(nmds$stress, 3),
                        " | PERMANOVA p = ", perm_p_lab),
      color = "Season"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9),
      legend.position = "bottom",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 9),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.margin = margin(6, 6, 6, 6)
    )
  
  list(
    plot = p,
    envfit_df = ef_df,
    permanova = permanova_result,
    stress = nmds$stress
  )
}

nmds_100 <- create_nmds_panel("100", microresp_prep)
nmds_25  <- create_nmds_panel("25",  microresp_prep)

# ==============================================================================
# PANEL C: SUMMARIZED MEAN HEATMAP (Treatment × Season)
# ==============================================================================
top_substrates <- unique(c(nmds_100$envfit_df$Substrate, nmds_25$envfit_df$Substrate))

substrate_complexity <- microresp_prep %>%
  distinct(Substrate, Substrate_category) %>%
  mutate(Complexity = case_when(
    Substrate_category %in% c("polymer", "protein", "fatty acid ester",
                              "plant wax ester", "phenolic") ~ "Complex",
    Substrate_category %in% c("carbohydrate", "amino acid A", "amino acid B",
                              "amino acid C", "amino acid D", "carboxylic acids",
                              "amine", "signaling molecule") ~ "Labile",
    Substrate_category == "control" ~ NA_character_,
    TRUE ~ "Labile"
  ))

heatmap_summary <- microresp_prep %>%
  filter(Substrate %in% top_substrates) %>%
  left_join(substrate_complexity, by = c("Substrate", "Substrate_category")) %>%
  group_by(Sample_cat) %>%
  mutate(Relative_SIR = 100 * SIR / sum(SIR, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Irrigation_trmt_pct, Season, Substrate, Complexity) %>%
  summarize(
    Mean_Relative_SIR = mean(Relative_SIR, na.rm = TRUE),
    SE_Relative_SIR   = sd(Relative_SIR, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(Group = paste0(Irrigation_trmt_pct, "% - ", Season))

heatmap_matrix <- heatmap_summary %>%
  select(Substrate, Group, Mean_Relative_SIR) %>%
  pivot_wider(names_from = Group, values_from = Mean_Relative_SIR, values_fill = 0) %>%
  column_to_rownames("Substrate") %>%
  as.matrix()

column_order <- c("100% - Growing", "100% - Non-Growing",
                  "25% - Growing",  "25% - Non-Growing")
heatmap_matrix <- heatmap_matrix[, column_order, drop = FALSE]

col_annotation <- data.frame(
  Treatment = c("100", "100", "25", "25"),
  Season    = c("Growing", "Non-Growing", "Growing", "Non-Growing"),
  row.names = column_order
)

row_annotation <- substrate_complexity %>%
  filter(Substrate %in% rownames(heatmap_matrix)) %>%
  select(Substrate, Complexity) %>%
  distinct() %>%
  column_to_rownames("Substrate")

ann_colors <- list(
  Treatment  = trt_cols,
  Season     = season_cols,
  Complexity = complex_cols
)

set.seed(123)
heatmap_plot <- pheatmap(
  heatmap_matrix,
  annotation_col = col_annotation,
  annotation_row = row_annotation,
  annotation_colors = ann_colors,
  scale = "row",
  cluster_cols = FALSE,
  clustering_distance_rows = "correlation",
  clustering_method = "ward.D2",
   annotation_legend = FALSE,     # << big win: removes the chunky legend blocks
  show_colnames = FALSE,   # << hides “25% - Growing”, etc.
  treeheight_row = 35,           # smaller dendrogram height footprint
  treeheight_col = 10,           # small since we aren't clustering cols
  cellwidth = 26,                # << makes 4 columns readable in a 1/2 page panel
  cellheight = 16,               # << enlarge tiles; you only have ~12 rows
  fontsize = 10,
  fontsize_row = 9,
  fontsize_col = 10,
  angle_col = 0,
  border_color = "grey70",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(60),
  silent = TRUE
)

panel_C <- ggplotify::as.ggplot(heatmap_plot$gtable) +
  theme_void() +
  theme(
    plot.margin = margin(2, 2, 2, 2) # tighter margins to reclaim space
    #plot.background = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

# ==============================================================================
# PANEL D: POTENTIAL VS REALIZED ACTIVITY RATIO (Boxplot)
# ==============================================================================
ratio_data <- microresp_prep %>%
  group_by(Sample_cat) %>%
  summarize(Total_MicroResp = sum(SIR, na.rm = TRUE), .groups = "drop") %>%
  left_join(biogeochem %>% select(Sample_cat, Resp_rate_mg_C_m2_h), by = "Sample_cat") %>%
  left_join(microresp_prep %>% distinct(Sample_cat, Season, Irrigation_trmt_pct, Collection_date),
            by = "Sample_cat") %>%
  filter(!is.na(Resp_rate_mg_C_m2_h), Resp_rate_mg_C_m2_h > 0, !is.na(Season)) %>%
  mutate(Ratio = Total_MicroResp / Resp_rate_mg_C_m2_h)

wilcox_result <- wilcox.test(Ratio ~ Season, data = ratio_data)
wilcox_pval <- wilcox_result$p.value
pval_label <- if (wilcox_pval < 0.001) "p < 0.001" else paste0("p = ", sprintf("%.3f", wilcox_pval))

ratio_summary <- ratio_data %>%
  group_by(Season) %>%
  summarize(median = median(Ratio, na.rm = TRUE), n = n(), .groups = "drop")

panel_D <- ggplot(ratio_data, aes(x = Season, y = Ratio, fill = Season)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA, linewidth = 0.7) +
  geom_jitter(width = 0.18, alpha = 0.45, size = 2.2) +
  scale_fill_manual(values = season_cols) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  annotation_logticks(sides = "l") +
  labs(
    subtitle = paste("Wilcoxon:", pval_label),
    y = "Total MicroResp / Field Respiration (log10)",
    x = NULL
  ) +
  annotate(
    "text",
    x = 1:2,
    y = ratio_summary$median * 60,
    label = paste0("n = ", ratio_summary$n),
    size = 3.3
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(size = 10),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    plot.margin = margin(6, 6, 6, 6)
  )

# ==============================================================================
# COMBINE PANELS (requested structure)
# Row 1: A | B
# Row 2: C | D
# Give heatmap more width (C wider than D)
# ==============================================================================
panel_A <- nmds_100$plot
panel_B <- nmds_25$plot

final_figure <- (panel_A + panel_B) / (panel_C + panel_D) +
  plot_layout(
    heights = c(1, 1.05),
    widths  = c(1.45, 1)   # << more width to heatmap than boxplot in row 2
  ) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 14, face = "bold"),
    plot.tag.position = c(0.01, 0.99)
  )

# ==============================================================================
# SAVE FIGURE (sheet-of-paper friendly)
# Letter landscape: 11 x 8.5 in
# ==============================================================================
ggsave("Figures/Figure_H2b_H3_MicroResp_Analysis_paper_v4.pdf",
       final_figure, width = 11, height = 8.5, units = "in", dpi = 300)

ggsave("Figures/Figure_H2b_H3_MicroResp_Analysis_paper_v4.png",
       final_figure, width = 11, height = 8.5, units = "in", dpi = 300)

cat("\n=== Figure saved successfully ===\n")
cat("Figures/Figure_H2b_H3_MicroResp_Analysis_paper_v3.pdf\n")
cat("Figures/Figure_H2b_H3_MicroResp_Analysis_paper_v3.png\n")

# ==============================================================================
# PRINT STATISTICAL RESULTS
# ==============================================================================
cat("\n================================================================================\n")
cat("STATISTICAL RESULTS FOR H2b AND H3 ANALYSIS\n")
cat("================================================================================\n")

cat("\n=== PANEL A: NMDS 100% & PERMANOVA ===\n")
cat("NMDS Stress:", round(nmds_100$stress, 4), "\n\n")
print(nmds_100$permanova)
cat("\nTop 8 Substrates (by R²):\n")
print(nmds_100$envfit_df %>%
        select(Substrate, r2, pval) %>%
        mutate(r2 = round(r2, 3), pval = round(pval, 4)))

cat("\n=== PANEL B: NMDS 25% & PERMANOVA ===\n")
cat("NMDS Stress:", round(nmds_25$stress, 4), "\n\n")
print(nmds_25$permanova)
cat("\nTop 8 Substrates (by R²):\n")
print(nmds_25$envfit_df %>%
        select(Substrate, r2, pval) %>%
        mutate(r2 = round(r2, 3), pval = round(pval, 4)))

cat("\n=== PANEL C: HEATMAP SUMMARY STATISTICS ===\n")
cat("\nMean Relative Respiration by Treatment × Season:\n")
print(heatmap_summary %>%
        group_by(Irrigation_trmt_pct, Season) %>%
        summarize(
          Mean_Total = mean(Mean_Relative_SIR),
          SD_Total = sd(Mean_Relative_SIR),
          .groups = "drop"
        ))

cat("\n=== PANEL D: WILCOXON TEST RESULTS ===\n")
print(wilcox_result)
cat("\nMedian Ratios by Season:\n")
print(ratio_summary)

cat("\n=== SAMPLE SIZES ===\n")
cat("\nSamples per Treatment × Season (for NMDS):\n")
print(microresp_prep %>%
        distinct(Sample_cat, Irrigation_trmt_pct, Season) %>%
        count(Irrigation_trmt_pct, Season) %>%
        pivot_wider(names_from = Season, values_from = n))

cat("\n================================================================================\n")

# ==============================================================================
# EXPORT SUPPLEMENTAL DATA
# ==============================================================================
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
write_csv(permanova_summary, "Results/Table_S1_PERMANOVA_results.csv")

top_substrates_summary <- bind_rows(
  nmds_100$envfit_df %>% mutate(Treatment = "100% Irrigation"),
  nmds_25$envfit_df %>% mutate(Treatment = "25% Irrigation")
) %>%
  select(Treatment, Substrate, r2, pval, NMDS1, NMDS2) %>%
  arrange(Treatment, desc(r2))
write_csv(top_substrates_summary, "Results/Table_S2_Top_Substrates_Envfit.csv")

write_csv(heatmap_summary, "Results/Table_S3_Heatmap_Summary_Data.csv")
write_csv(ratio_data, "Results/Table_S4_Potential_Realized_Ratios.csv")

cat("\n=== Supplemental files exported ===\n")
cat("- Table_S1_PERMANOVA_results.csv\n")
cat("- Table_S2_Top_Substrates_Envfit.csv\n")
cat("- Table_S3_Heatmap_Summary_Data.csv\n")
cat("- Table_S4_Potential_Realized_Ratios.csv\n")

# ==============================================================================
# ADDITIONAL STATISTICAL TESTS FOR H2b (kept)
# ==============================================================================
cat("\n=== ADDITIONAL ANALYSIS FOR H2b ===\n")
cat("Testing substrate complexity differences between treatments\n\n")

complexity_test_data <- microresp_prep %>%
  left_join(substrate_complexity, by = c("Substrate", "Substrate_category")) %>%
  filter(Complexity %in% c("Complex", "Labile"), !is.na(Complexity)) %>%
  group_by(Sample_cat, Irrigation_trmt_pct, Season, Complexity) %>%
  summarize(Total_SIR = sum(SIR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Complexity, values_from = Total_SIR, values_fill = 0) %>%
  mutate(Complex_Labile_Ratio = Complex / (Labile + 0.001))

complexity_anova <- aov(Complex_Labile_Ratio ~ Irrigation_trmt_pct * Season,
                        data = complexity_test_data)

cat("ANOVA: Complex/Labile Ratio ~ Treatment × Season\n")
print(summary(complexity_anova))

cat("\nMean Complex/Labile Ratios:\n")
complexity_means <- complexity_test_data %>%
  group_by(Irrigation_trmt_pct, Season) %>%
  summarize(
    Mean_Ratio = mean(Complex_Labile_Ratio, na.rm = TRUE),
    SD_Ratio   = sd(Complex_Labile_Ratio, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
print(complexity_means)

write_csv(complexity_test_data, "Results/Table_S5_Complexity_Ratios.csv")

cat("\n================================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================================\n")
