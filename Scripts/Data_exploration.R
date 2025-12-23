# ============================================================================
# STEP 1: DATA EXPLORATION FOR RF ANALYSIS OF FIELD RESPIRATION
# Soils SFA A3.1 Prosser Time Series
# Research Question: Seasonal and moisture effects on respiration predictors
# ============================================================================

# Load required packages
rm(list = ls())
library(tidyverse)
library(lubridate)
library(corrplot)      # Correlation visualization
library(ggcorrplot)    # Alternative correlation plots
library(GGally)        # Pairs plots
library(naniar)        # Missing data visualization
library(VIM)           # Missing data patterns
library(caret)         # For findCorrelation function

# Set theme for consistent plotting
theme_set(theme_bw(base_size = 12))

# ----------------------------------------------------------------------------
# 1. LOAD AND PREPARE DATA
# ----------------------------------------------------------------------------

# Load biogeochemistry data
bgc_data <- read.csv("DataHub/Soil_Biogeochemistry/SoilSFA_A3.1_Soil_Biogeochemistry.csv") 

# Load EEA data
eea_data <- read.csv("DataHub/Soil_Phenotypes/SoilSFA_A3.1_EEA.csv") %>%
  rename(Irrigation_trmt_cat = Irrigaiton_trmt_cat)

# Parse dates and add temporal variables
bgc_data <- bgc_data %>%
  mutate(
    Collection_date = mdy(Collection_date),
    Month = month(Collection_date),
    Year = year(Collection_date),
    Month_name = month(Collection_date, label = TRUE, abbr = TRUE),
    
    # Define season based on irrigation period (April-October = Growing)
    Season = case_when(
      Month %in% c(4, 5, 6, 7, 8, 9, 10) ~ "Growing",
      Month %in% c(11, 12, 1, 2, 3) ~ "Non-Growing",
      TRUE ~ NA_character_
    ),
    Season = factor(Season, levels = c("Growing", "Non-Growing")),
    
    # Convert respiration to mg C m-2 h-1
    Resp_rate_mg_C_m2_h = CO2_C_60_m2 / 1000,
    
    # Create irrigation treatment factor
    Irrigation_trmt_pct = factor(Irrigation_trmt_pct, levels = c("25", "100"))
  )

# Merge with EEA data
rf_data <- bgc_data %>%
  left_join(
    eea_data %>% select(Sample_cat, Irrigation_trmt_cat, BG_nmol_gdw_h:MBN_ug_gdw),
    by = "Sample_cat",
    suffix = c("_bgc", "_eea")
  )

# Check merge success
cat("\n=== DATA MERGE CHECK ===\n")
cat("BGC samples:", nrow(bgc_data), "\n")
cat("EEA samples:", nrow(eea_data), "\n")
cat("Merged samples:", nrow(rf_data), "\n")
cat("Samples with respiration data:", sum(!is.na(rf_data$Resp_rate_mg_C_m2_h)), "\n")

# ----------------------------------------------------------------------------
# 2. DEFINE PREDICTOR VARIABLES (EXCLUDING RESIN STRIPS)
# ----------------------------------------------------------------------------

# Categorize predictors by biogeochemical function
soil_chemistry <- c("GWC", "pH", "SEOC_ug_gdw", "SETN_ug_gdw")

microbial_biomass <- c("MBC_ug_gdw_bgc", "MBN_ug_gdw_bgc")

extractable_P <- c("Olson_IP_ug_gdw", "Oslon_OP_ug_gdw", 
                   "MBP_IP_ug_gdw", "MBP_OP_ug_gdw")

# Individual enzyme activities (C-degrading)
enzymes_C <- c("BG_nmol_gdw_h", "BX_nmol_gdw_h", "CB_nmol_gdw_h")

# Individual enzyme activities (N-acquiring)
enzymes_N <- c("NAG_nmol_gdw_h", "Leu_nmol_gdw_h", "Ala_nmol_gdw_h")

# Individual enzyme activities (P-acquiring)
enzymes_P <- c("AP_nmol_gdw_h")

# Enzyme sums and ratios
enzyme_ratios <- c("Sum_Carbon", "Sum_Nitrogen", "Sum_Phosphorus",
                   "EEA_C_N", "EEA_C_P", "EEA_N_P")

# Environmental variables
environmental <- c("SoilTemp_mean", "SoilTemp_min", "SoilTemp_max", 
                   "Roots_gdw")

# ALL PREDICTORS (for initial assessment)
all_predictors <- c(
  soil_chemistry, 
  microbial_biomass, 
  extractable_P,
  enzymes_C, 
  enzymes_N, 
  enzymes_P,
  enzyme_ratios,
  environmental
)

# Response variable
response_var <- "Resp_rate_mg_C_m2_h"

# Create analysis dataset
rf_analysis <- rf_data %>%
  select(
    # Identifiers
    Sample_cat, Collection_date, Month, Month_name, Season, Year,
    Plot, Irrigation_trmt_pct, Irrigation_trmt_cat, Block, Rep,
    
    # All predictors
    all_of(all_predictors), 
    
    # Response
    all_of(response_var)
  ) %>%
  # Keep only samples with respiration data
  filter(!is.na(!!sym(response_var)))

cat("\n=== ANALYSIS DATASET ===\n")
cat("Total samples for RF analysis:", nrow(rf_analysis), "\n")
cat("Growing season samples:", sum(rf_analysis$Season == "Growing", na.rm = TRUE), "\n")
cat("Non-growing season samples:", sum(rf_analysis$Season == "Non-Growing", na.rm = TRUE), "\n")
cat("100% irrigation samples:", sum(rf_analysis$Irrigation_trmt_pct == "100", na.rm = TRUE), "\n")
cat("25% irrigation samples:", sum(rf_analysis$Irrigation_trmt_pct == "25", na.rm = TRUE), "\n")

# ----------------------------------------------------------------------------
# 3. MISSING DATA ASSESSMENT
# ----------------------------------------------------------------------------

cat("\n=== MISSING DATA SUMMARY ===\n\n")

# Calculate missing data for all predictors
missing_summary <- rf_analysis %>%
  select(all_of(all_predictors), all_of(response_var)) %>%
  summarise(across(everything(), 
                   list(
                     n_missing = ~sum(is.na(.)),
                     pct_missing = ~round(sum(is.na(.)) / n() * 100, 1)
                   ))) %>%
  pivot_longer(everything(), 
               names_to = c("Variable", ".value"),
               names_pattern = "(.+)_(n_missing|pct_missing)") %>%
  arrange(desc(pct_missing))

print(missing_summary, n = 30)

# Visualize missing data patterns
cat("\n=== GENERATING MISSING DATA VISUALIZATIONS ===\n")

# Missing data heatmap
p_missing_heat <- rf_analysis %>%
  select(all_of(all_predictors)) %>%
  vis_miss() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Missing Data Pattern Across All Predictors")

print(p_missing_heat)

# Missing data by month
missing_by_month <- rf_analysis %>%
  select(Month_name, Season, all_of(all_predictors)) %>%
  group_by(Month_name, Season) %>%
  summarise(
    across(all_of(all_predictors), ~sum(is.na(.))),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = -c(Month_name, Season),
               names_to = "Variable",
               values_to = "N_Missing")

p_missing_month <- ggplot(missing_by_month, 
                          aes(x = Month_name, y = Variable, fill = N_Missing)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", name = "N Missing") +
  facet_wrap(~Season, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Missing Data by Month",
       x = "Month", y = "Predictor Variable")

print(p_missing_month)

# ----------------------------------------------------------------------------
# 4. COMPLETE CASE ANALYSIS
# ----------------------------------------------------------------------------

cat("\n=== COMPLETE CASE ANALYSIS ===\n\n")

# Create complete case dataset
rf_complete <- rf_analysis %>%
  select(Sample_cat, Collection_date, Month_name, Season, 
         Irrigation_trmt_pct, Plot, Block, Rep,
         all_of(all_predictors), all_of(response_var)) %>%
  drop_na()

cat("Complete cases (all predictors):", nrow(rf_complete), "\n")
cat("Data retention:", round(nrow(rf_complete) / nrow(rf_analysis) * 100, 1), "%\n\n")

# Breakdown by season and irrigation
complete_summary <- rf_complete %>%
  group_by(Season, Irrigation_trmt_pct) %>%
  summarise(
    N_samples = n(),
    N_months = n_distinct(Month_name),
    Months = paste(unique(Month_name), collapse = ", "),
    .groups = "drop"
  )

cat("Complete cases by Season × Irrigation:\n")
print(complete_summary)

# Sample distribution across months
complete_by_month <- rf_complete %>%
  count(Month_name, Season, Irrigation_trmt_pct) %>%
  arrange(Month_name)

cat("\nComplete samples by month:\n")
print(complete_by_month, n = 30)

# ----------------------------------------------------------------------------
# 5. DATA DISTRIBUTION EXPLORATION
# ----------------------------------------------------------------------------

cat("\n=== EXPLORING DATA DISTRIBUTIONS ===\n")

# Summary statistics for complete dataset
summary_stats <- rf_complete %>%
  select(all_of(all_predictors), all_of(response_var)) %>%
  summarise(across(everything(), 
                   list(
                     min = ~min(., na.rm = TRUE),
                     q25 = ~quantile(., 0.25, na.rm = TRUE),
                     median = ~median(., na.rm = TRUE),
                     mean = ~mean(., na.rm = TRUE),
                     q75 = ~quantile(., 0.75, na.rm = TRUE),
                     max = ~max(., na.rm = TRUE),
                     sd = ~sd(., na.rm = TRUE)
                   ))) %>%
  pivot_longer(everything(),
               names_to = c("Variable", ".value"),
               names_pattern = "(.+)_(min|q25|median|mean|q75|max|sd)")

cat("\nSummary statistics saved. Displaying first 10 variables:\n")
print(head(summary_stats, 10))

# Visualize distributions - Response variable
p_resp_dist <- ggplot(rf_complete, aes(x = Resp_rate_mg_C_m2_h)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_density(aes(y = after_stat(count)), color = "red", linewidth = 1) +
  facet_grid(Season ~ Irrigation_trmt_pct, scales = "free") +
  labs(title = "Distribution of Respiration Rate",
       subtitle = "By Season and Irrigation Treatment",
       x = "Respiration Rate (mg C m⁻² h⁻¹)",
       y = "Count")

print(p_resp_dist)

# Boxplot of respiration by season and irrigation
p_resp_box <- ggplot(rf_complete, 
                     aes(x = Month_name, y = Resp_rate_mg_C_m2_h, 
                         fill = Irrigation_trmt_pct)) +
  geom_boxplot() +
  facet_wrap(~Season, scales = "free_x") +
  scale_fill_manual(values = c("25" = "#d95f02", "100" = "#1b9e77"),
                    name = "Irrigation",
                    labels = c("25%" = "25% (Drought)", 
                               "100%" = "100% (Control)")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Field Respiration Across Time",
       subtitle = "Comparing irrigation treatments by season",
       x = "Month", 
       y = "Respiration Rate (mg C m⁻² h⁻¹)")

print(p_resp_box)

# ----------------------------------------------------------------------------
# 6. CORRELATION ANALYSIS
# ----------------------------------------------------------------------------

cat("\n=== CORRELATION ANALYSIS ===\n")

# Calculate correlation matrix (complete cases only)
cor_matrix <- rf_complete %>%
  select(all_of(all_predictors), all_of(response_var)) %>%
  cor(use = "pairwise.complete.obs", method = "pearson")

# Identify highly correlated predictors (|r| > 0.8)
high_cor <- findCorrelation(cor_matrix, cutoff = 0.8, names = TRUE, 
                            verbose = FALSE)

cat("\nHighly correlated predictors (|r| > 0.8):\n")
if(length(high_cor) > 0) {
  print(high_cor)
} else {
  cat("No predictors with |r| > 0.8\n")
}

# Check correlation with response variable
cor_with_response <- cor_matrix[, response_var] %>%
  as.data.frame() %>%
  rownames_to_column("Variable") %>%
  rename(Correlation = 2) %>%
  filter(Variable != response_var) %>%
  arrange(desc(abs(Correlation)))

cat("\nTop 15 predictors correlated with respiration:\n")
print(head(cor_with_response, 15))

# Visualization: Full correlation matrix
p_cor_full <- ggcorrplot(cor_matrix, 
                         hc.order = TRUE,
                         type = "lower",
                         lab = FALSE,
                         outline.col = "white",
                         title = "Correlation Matrix: All Predictors") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8))

print(p_cor_full)

# Visualization: Correlation with respiration only
p_cor_resp <- cor_with_response %>%
  slice_head(n = 20) %>%
  mutate(Variable = fct_reorder(Variable, Correlation)) %>%
  ggplot(aes(x = Correlation, y = Variable, fill = Correlation)) +
  geom_col() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Top 20 Predictors Correlated with Respiration",
       x = "Pearson Correlation Coefficient",
       y = "Predictor Variable")

print(p_cor_resp)

# ----------------------------------------------------------------------------
# 7. EXAMINE SPECIFIC PROBLEMATIC CORRELATIONS
# ----------------------------------------------------------------------------

cat("\n=== EXAMINING MULTICOLLINEARITY ISSUES ===\n\n")

# Check correlation among enzyme ratios vs individual enzymes
enzyme_cor <- rf_complete %>%
  select(all_of(c(enzymes_C, enzymes_N, enzymes_P, enzyme_ratios))) %>%
  cor(use = "pairwise.complete.obs")

cat("Correlation between individual enzymes and their sums:\n")
# Sum_Carbon vs BG, BX, CB
cat("\nSum_Carbon correlations:\n")
print(round(enzyme_cor["Sum_Carbon", c("BG_nmol_gdw_h", "BX_nmol_gdw_h", "CB_nmol_gdw_h")], 3))

cat("\nSum_Nitrogen correlations:\n")
print(round(enzyme_cor["Sum_Nitrogen", c("NAG_nmol_gdw_h", "Leu_nmol_gdw_h", "Ala_nmol_gdw_h")], 3))

# Check soil temperature correlations (mean, min, max)
temp_cor <- rf_complete %>%
  select(starts_with("SoilTemp")) %>%
  cor(use = "pairwise.complete.obs")

cat("\n\nSoil temperature variable correlations:\n")
print(round(temp_cor, 3))

# Check MBC/MBN from two sources (bgc vs eea files)
if("MBC_ug_gdw_eea" %in% names(rf_complete) & "MBC_ug_gdw_bgc" %in% names(rf_complete)) {
  mb_comparison <- rf_complete %>%
    select(Sample_cat, MBC_ug_gdw_bgc, MBC_ug_gdw_eea, 
           MBN_ug_gdw_bgc, MBN_ug_gdw_eea) %>%
    drop_na()
  
  cat("\n\nMicrobial Biomass C comparison (BGC vs EEA files):\n")
  cat("Samples with both measurements:", nrow(mb_comparison), "\n")
  cat("Correlation MBC:", round(cor(mb_comparison$MBC_ug_gdw_bgc, 
                                    mb_comparison$MBC_ug_gdw_eea), 3), "\n")
  cat("Correlation MBN:", round(cor(mb_comparison$MBN_ug_gdw_bgc, 
                                    mb_comparison$MBN_ug_gdw_eea), 3), "\n")
}

# ----------------------------------------------------------------------------
# 8. IDENTIFY PREDICTOR GROUPS WITH HIGH INTERNAL CORRELATION
# ----------------------------------------------------------------------------

cat("\n=== INTERNAL CORRELATION WITHIN PREDICTOR GROUPS ===\n\n")

# Function to check within-group correlation
check_group_correlation <- function(data, vars, group_name) {
  if(length(vars) > 1) {
    cor_mat <- data %>%
      select(all_of(vars)) %>%
      cor(use = "pairwise.complete.obs")
    
    cat(group_name, ":\n")
    print(round(cor_mat, 2))
    cat("\n")
    
    # Find pairs with |r| > 0.8
    high_pairs <- which(abs(cor_mat) > 0.8 & abs(cor_mat) < 1, arr.ind = TRUE)
    if(nrow(high_pairs) > 0) {
      cat("  High correlations (|r| > 0.8):\n")
      for(i in 1:nrow(high_pairs)) {
        var1 <- rownames(cor_mat)[high_pairs[i, 1]]
        var2 <- colnames(cor_mat)[high_pairs[i, 2]]
        if(var1 < var2) {  # Avoid duplicates
          cat("    ", var1, "vs", var2, ": r =", 
              round(cor_mat[high_pairs[i, 1], high_pairs[i, 2]], 3), "\n")
        }
      }
    } else {
      cat("  No high correlations found\n")
    }
    cat("\n")
  }
}

# Check each group
check_group_correlation(rf_complete, soil_chemistry, "Soil Chemistry")
check_group_correlation(rf_complete, microbial_biomass, "Microbial Biomass")
check_group_correlation(rf_complete, extractable_P, "Extractable Phosphorus")
check_group_correlation(rf_complete, c(enzymes_C, enzymes_N, enzymes_P), 
                        "Individual Enzymes")
check_group_correlation(rf_complete, enzyme_ratios, "Enzyme Ratios")
check_group_correlation(rf_complete, environmental, "Environmental")

# ----------------------------------------------------------------------------
# 9. VISUALIZE KEY RELATIONSHIPS
# ----------------------------------------------------------------------------

cat("\n=== CREATING RELATIONSHIP VISUALIZATIONS ===\n")

# Scatterplot matrix for top correlated variables with respiration
top_vars <- cor_with_response %>%
  slice_head(n = 6) %>%
  pull(Variable)

p_pairs <- rf_complete %>%
  select(all_of(top_vars), Resp_rate_mg_C_m2_h, Season, Irrigation_trmt_pct) %>%
  ggpairs(
    columns = c(top_vars, "Resp_rate_mg_C_m2_h"),
    aes(color = Season, alpha = 0.5),
    upper = list(continuous = wrap("cor", size = 3)),
    lower = list(continuous = wrap("points", alpha = 0.3, size = 0.5)),
    title = "Relationships: Top 6 Predictors vs Respiration"
  ) +
  theme(strip.text = element_text(size = 8))

print(p_pairs)

# ----------------------------------------------------------------------------
# 10. SAVE RESULTS AND CREATE SUMMARY REPORT
# ----------------------------------------------------------------------------

# Save complete case dataset
write.csv(rf_complete, "RF_analysis_complete_cases.csv", row.names = FALSE)

# Save missing data summary
write.csv(missing_summary, "Missing_data_summary.csv", row.names = FALSE)

# Save correlation summary
write.csv(cor_with_response, "Correlation_with_respiration.csv", row.names = FALSE)

cat("\n=== FILES SAVED ===\n")
cat("1. RF_analysis_complete_cases.csv - Complete case dataset for RF\n")
cat("2. Missing_data_summary.csv - Missing data by variable\n")
cat("3. Correlation_with_respiration.csv - Correlations with response\n")

# Create summary report
cat("\n")
cat("==============================================================\n")
cat("        DATA COMPLETENESS SUMMARY FOR RF ANALYSIS            \n")
cat("==============================================================\n\n")

cat("STARTING DATASET:\n")
cat("  Total samples with respiration:", nrow(rf_analysis), "\n")
cat("  Time period: June 2023 - June 2024\n")
cat("  Number of predictors:", length(all_predictors), "\n\n")

cat("COMPLETE CASE DATASET:\n")
cat("  Complete samples:", nrow(rf_complete), "\n")
cat("  Data retention:", round(nrow(rf_complete) / nrow(rf_analysis) * 100, 1), "%\n")
cat("  Samples per predictor ratio:", 
    round(nrow(rf_complete) / length(all_predictors), 2), "\n\n")

cat("BY STRATIFICATION:\n")
for(i in 1:nrow(complete_summary)) {
  cat("  ", as.character(complete_summary$Season[i]), "×", 
      as.character(complete_summary$Irrigation_trmt_pct[i]), "%: ",
      complete_summary$N_samples[i], " samples (", 
      complete_summary$N_months[i], " months)\n", sep = "")
}

cat("\n==============================================================\n\n")

#######################################################################
# Check for outliers 
#######################################################################
# ============================================================================
# SYSTEMATIC OUTLIER DETECTION ACROSS ALL PREDICTORS
# ============================================================================

cat("\n=== OUTLIER DETECTION FOR ALL PREDICTORS ===\n\n")

# Function to detect outliers using multiple methods
detect_outliers <- function(data, variables) {
  
  outlier_results <- list()
  
  for(var in variables) {
    
    values <- data[[var]]
    
    # Method 1: Tukey's method (IQR-based)
    Q1 <- quantile(values, 0.25, na.rm = TRUE)
    Q3 <- quantile(values, 0.75, na.rm = TRUE)
    IQR <- Q3 - Q1
    lower_fence <- Q1 - 3 * IQR  # 3*IQR for extreme outliers
    upper_fence <- Q3 + 3 * IQR
    
    tukey_outliers <- which(values < lower_fence | values > upper_fence)
    
    # Method 2: Z-score (>3 or <-3 is extreme outlier)
    z_scores <- scale(values)
    zscore_outliers <- which(abs(z_scores) > 3)
    
    # Method 3: Statistical test (Grubbs test would require nortest package)
    # For now, just flag values >3 SD from mean
    
    if(length(tukey_outliers) > 0 | length(zscore_outliers) > 0) {
      outlier_results[[var]] <- list(
        tukey = tukey_outliers,
        zscore = zscore_outliers,
        n_tukey = length(tukey_outliers),
        n_zscore = length(zscore_outliers)
      )
    }
  }
  
  return(outlier_results)
}

# Detect outliers in all predictors
outliers_detected <- detect_outliers(rf_complete, all_predictors)

# Report outliers
cat("Variables with potential outliers:\n\n")

if(length(outliers_detected) > 0) {
  for(var_name in names(outliers_detected)) {
    cat("Variable:", var_name, "\n")
    cat("  Tukey method (3*IQR):", outliers_detected[[var_name]]$n_tukey, "outliers\n")
    cat("  Z-score method (|z|>3):", outliers_detected[[var_name]]$n_zscore, "outliers\n")
    
    # Show which samples
    if(outliers_detected[[var_name]]$n_tukey > 0) {
      outlier_samples <- rf_complete[outliers_detected[[var_name]]$tukey, 
                                     c("Sample_cat", "Month_name", "Season", 
                                       "Irrigation_trmt_pct", var_name)]
      cat("  Samples flagged by Tukey:\n")
      print(outlier_samples)
    }
    cat("\n")
  }
} else {
  cat("No extreme outliers detected in any predictors.\n")
}

# ============================================================================
# VISUALIZE OUTLIERS WITH BOXPLOTS
# ============================================================================

cat("\n=== CREATING OUTLIER VISUALIZATION ===\n")

# Create boxplots for all predictors to visually identify outliers
rf_complete_long <- rf_complete %>%
  select(Sample_cat, Season, Irrigation_trmt_pct, all_of(all_predictors)) %>%
  pivot_longer(cols = all_of(all_predictors),
               names_to = "Variable",
               values_to = "Value")

# Boxplots by variable
p_outliers <- ggplot(rf_complete_long, aes(x = Variable, y = Value)) +
  geom_boxplot(outlier.color = "red", outlier.size = 3) +
  facet_wrap(~Variable, scales = "free", ncol = 4) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(title = "Outlier Detection Across All Predictors",
       subtitle = "Red points indicate potential outliers (>1.5*IQR from quartiles)",
       x = "", y = "Value")

print(p_outliers)

# ============================================================================
# CHECK SPECIFIC SUSPICIOUS VALUES
# ============================================================================

cat("\n=== CHECKING FOR SUSPICIOUS VALUES ===\n\n")

# Check for zeros in key variables (often data quality issues)
zero_check <- rf_complete %>%
  select(Sample_cat, all_of(all_predictors)) %>%
  filter(if_any(all_of(all_predictors), ~ . == 0))

cat("Samples with zero values in predictors:\n")
cat("Number of samples:", nrow(zero_check), "\n")
if(nrow(zero_check) > 0) {
  # Count zeros by variable
  zero_counts <- zero_check %>%
    select(-Sample_cat) %>%
    summarise(across(everything(), ~sum(. == 0, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "Variable", values_to = "N_zeros") %>%
    filter(N_zeros > 0) %>%
    arrange(desc(N_zeros))
  
  print(zero_counts)
}

# Check for extreme values in each predictor
cat("\n\nExtreme values (>95th percentile) by variable:\n")
extreme_summary <- rf_complete %>%
  select(all_of(all_predictors)) %>%
  summarise(across(everything(),
                   list(
                     p95 = ~quantile(., 0.95, na.rm = TRUE),
                     max = ~max(., na.rm = TRUE),
                     n_above_p95 = ~sum(. > quantile(., 0.95, na.rm = TRUE))
                   ))) %>%
  pivot_longer(everything(),
               names_to = c("Variable", ".value"),
               names_pattern = "(.+)_(p95|max|n_above_p95)") %>%
  mutate(ratio_max_to_p95 = max / p95) %>%
  arrange(desc(ratio_max_to_p95))

print(extreme_summary, n = 20)

# ============================================================================
# DECISION RECOMMENDATION
# ============================================================================

cat("\n=== OUTLIER REMOVAL RECOMMENDATIONS ===\n\n")

# Identify samples to potentially remove
samples_to_investigate <- rf_complete %>%
  filter(
    SEOC_ug_gdw > 2000 |  # The extreme SEOC value
      MBP_IP_ug_gdw > 500 |  # Check for extreme P values
      MBP_OP_ug_gdw > 400    # Check for extreme P values
  ) %>%
  select(Sample_cat, Collection_date, Month_name, Plot, 
         SEOC_ug_gdw, SETN_ug_gdw, MBP_IP_ug_gdw, MBP_OP_ug_gdw,
         Resp_rate_mg_C_m2_h)

cat("Samples recommended for investigation/removal:\n")
print(samples_to_investigate)

cat("\n\nRECOMMENDATION:\n")
cat("1. REMOVE Sample SM1082 (SEOC = 2410) - clear measurement error\n")
cat("2. INVESTIGATE samples with extreme MBP values (if any)\n")
cat("3. CHECK samples with zero values for biological plausibility\n")
#################################################################
# ============================================================================
# IDENTIFYING ZERO VALUES AND OUTLIERS FOR REMOVAL
# ============================================================================

cat("\n=== REMOVING OUTLIER SAMPLES ===\n")

# First, remove the 4 flagged outlier samples from July 2023
outlier_samples <- c("SM1082", "SM1085", "SM1087", "SM1088")

rf_no_outliers <- rf_complete %>%
  filter(!Sample_cat %in% outlier_samples)

cat("Samples removed (July 2023 outliers):", length(outlier_samples), "\n")
cat("Remaining samples:", nrow(rf_no_outliers), "\n\n")

# Check how this affects stratification
outlier_impact <- rf_no_outliers %>%
  group_by(Season, Irrigation_trmt_pct) %>%
  summarise(N = n(), .groups = "drop")

cat("Stratification after outlier removal:\n")
print(outlier_impact)

# ============================================================================
# COMPREHENSIVE ZERO VALUE ANALYSIS
# ============================================================================

cat("\n\n=== ANALYZING ZERO VALUES ACROSS ALL PREDICTORS ===\n\n")

# Check for zeros in each predictor
zero_analysis <- rf_no_outliers %>%
  select(all_of(all_predictors)) %>%
  summarise(across(everything(), 
                   list(
                     n_zeros = ~sum(. == 0, na.rm = TRUE),
                     pct_zeros = ~round(sum(. == 0, na.rm = TRUE) / n() * 100, 1),
                     min_nonzero = ~min(.[. > 0], na.rm = TRUE)
                   ))) %>%
  pivot_longer(everything(),
               names_to = c("Variable", ".value"),
               names_pattern = "(.+)_(n_zeros|pct_zeros|min_nonzero)") %>%
  filter(n_zeros > 0) %>%
  arrange(desc(n_zeros))

cat("Variables with zero values:\n")
print(zero_analysis, n = 30)

# Total number of samples with at least one zero
samples_with_zeros <- rf_no_outliers %>%
  select(Sample_cat, all_of(all_predictors)) %>%
  filter(if_any(all_of(all_predictors), ~ . == 0))

cat("\n\nTotal samples containing at least one zero value:", nrow(samples_with_zeros), "\n")
cat("Percentage of dataset:", round(nrow(samples_with_zeros) / nrow(rf_no_outliers) * 100, 1), "%\n")

# Identify which samples have zeros in which variables
cat("\n\n=== SAMPLES WITH ZERO VALUES BY VARIABLE ===\n\n")

samples_with_zeros_detail <- rf_no_outliers %>%
  select(Sample_cat, Collection_date, Month_name, Season, 
         Irrigation_trmt_pct, all_of(all_predictors)) %>%
  pivot_longer(cols = all_of(all_predictors),
               names_to = "Variable",
               values_to = "Value") %>%
  filter(Value == 0) %>%
  arrange(Sample_cat, Variable)

# Count zeros per sample
zeros_per_sample <- samples_with_zeros_detail %>%
  count(Sample_cat, Month_name, Season, Irrigation_trmt_pct, name = "N_zeros") %>%
  arrange(desc(N_zeros))

cat("Samples with most zero values:\n")
print(head(zeros_per_sample, 15))

# Count zeros per variable
zeros_per_variable <- samples_with_zeros_detail %>%
  count(Variable, name = "N_samples_with_zero") %>%
  arrange(desc(N_samples_with_zero))

cat("\n\nVariables with most zero values:\n")
print(zeros_per_variable)

# ============================================================================
# VISUALIZE ZERO PATTERNS
# ============================================================================

cat("\n=== CREATING ZERO VALUE VISUALIZATIONS ===\n")

# Heatmap of zeros by sample and variable
p_zeros_heat <- samples_with_zeros_detail %>%
  mutate(has_zero = 1) %>%
  ggplot(aes(x = Variable, y = Sample_cat, fill = has_zero)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "darkred") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
        axis.text.y = element_text(size = 4)) +
  labs(title = "Heatmap of Zero Values",
       subtitle = "Samples (rows) × Variables (columns)",
       x = "Variable", y = "Sample")

print(p_zeros_heat)

# Zeros by month and season
zeros_by_month <- samples_with_zeros_detail %>%
  count(Month_name, Season, Variable) %>%
  ggplot(aes(x = Month_name, y = Variable, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), size = 2.5) +
  scale_fill_gradient(low = "white", high = "darkred", name = "N Zeros") +
  facet_wrap(~Season, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Zero Values by Month and Variable",
       x = "Month", y = "Variable")

print(zeros_by_month)

# ============================================================================
# BIOLOGICAL INTERPRETATION OF ZEROS
# ============================================================================

cat("\n\n=== INTERPRETATION OF ZERO VALUES ===\n\n")

cat("According to the README:\n")
cat("'Data from carbon, nitrogen and phosphorus analyses were checked for\n")
cat(" negative values which indicate that sample concentrations were lower\n")
cat(" than process blanks or controls. Negative values were converted to\n")
cat(" zeros for analyses but could also be considered below the detection limit.'\n\n")

# Which variables have zeros and are they biologically plausible?
cat("Variables with zeros - Interpretation:\n\n")

if("SEOC_ug_gdw" %in% zero_analysis$Variable) {
  cat("• SEOC_ug_gdw: Zeros indicate below detection limit (measurement artifact)\n")
}
if("SETN_ug_gdw" %in% zero_analysis$Variable) {
  cat("• SETN_ug_gdw: Zeros indicate below detection limit (measurement artifact)\n")
}
if("MBC_ug_gdw_bgc" %in% zero_analysis$Variable) {
  cat("• MBC_ug_gdw: Zeros indicate below detection limit (measurement artifact)\n")
}
if("MBN_ug_gdw_bgc" %in% zero_analysis$Variable) {
  cat("• MBN_ug_gdw: Zeros indicate below detection limit (measurement artifact)\n")
}
if("Oslon_OP_ug_gdw" %in% zero_analysis$Variable) {
  cat("• Oslon_OP_ug_gdw: Zeros may indicate below detection OR true absence\n")
}
if("Olson_IP_ug_gdw" %in% zero_analysis$Variable) {
  cat("• Olson_IP_ug_gdw: Zeros may indicate below detection OR true absence\n")
}
if("MBP_OP_ug_gdw" %in% zero_analysis$Variable) {
  cat("• MBP_OP_ug_gdw: Zeros may indicate below detection OR true absence\n")
}
if("MBP_IP_ug_gdw" %in% zero_analysis$Variable) {
  cat("• MBP_IP_ug_gdw: Zeros may indicate below detection OR true absence\n")
}

# ============================================================================
# DECISION POINT: HANDLING ZEROS
# ============================================================================

cat("\n\n=== OPTIONS FOR HANDLING ZERO VALUES ===\n\n")

cat("OPTION 1: Remove all samples with any zero values\n")
rf_no_zeros_complete <- rf_no_outliers %>%
  filter(if_all(all_of(all_predictors), ~ . > 0))

cat("  Samples retained:", nrow(rf_no_zeros_complete), 
    "out of", nrow(rf_no_outliers), 
    "(", round(nrow(rf_no_zeros_complete)/nrow(rf_no_outliers)*100, 1), "%)\n")

if(nrow(rf_no_zeros_complete) > 0) {
  strat_no_zeros <- rf_no_zeros_complete %>%
    group_by(Season, Irrigation_trmt_pct) %>%
    summarise(N = n(), .groups = "drop")
  print(strat_no_zeros)
}

cat("\n\nOPTION 2: Replace zeros with half the minimum detected value (LOD/2 approach)\n")
cat("  Common approach in analytical chemistry for values below detection\n")
cat("  Minimum non-zero values by variable:\n")

min_values <- zero_analysis %>%
  filter(n_zeros > 0) %>%
  select(Variable, min_nonzero) %>%
  mutate(LOD_half = round(min_nonzero / 2, 3))

print(min_values)

cat("\n\nOPTION 3: Replace zeros with a small constant (e.g., 0.001)\n")
cat("  Simple but may not reflect true detection limits\n")

cat("\n\nOPTION 4: Remove only variables with many zeros (if >X% zeros)\n")
cat("  Then handle remaining zeros with replacement\n")

# ============================================================================
# RECOMMENDED APPROACH
# ============================================================================

cat("\n\n=== RECOMMENDED APPROACH ===\n\n")

cat("Based on the README and typical analytical chemistry practices:\n\n")

cat("1. REMOVE outlier samples: SM1082, SM1085, SM1087, SM1088 (done above)\n\n")

cat("2. REPLACE zeros with LOD/2 (half minimum detected value) because:\n")
cat("   • README states zeros = 'below detection limit' not true biological zeros\n")
cat("   • This is standard practice for analytical chemistry data\n")
cat("   • Preserves sample size (critical for your small non-growing season groups)\n")
cat("   • More conservative than removing samples entirely\n\n")

cat("3. Variables to replace:\n")
cat("   • SEOC, SETN, MBC, MBN (analytical detection limits)\n")
cat("   • Olson_IP, Oslon_OP (may be true low values OR detection limits)\n")
cat("   • MBP_IP, MBP_OP (may be true low values OR detection limits)\n\n")

cat("This approach will retain maximum sample size while addressing\n")
cat("the measurement artifact issue documented in your README.\n")

# ============================================================================
# SAVE ZERO ANALYSIS RESULTS
# ============================================================================

write.csv(zero_analysis, "Zero_values_summary.csv", row.names = FALSE)
write.csv(zeros_per_sample, "Samples_with_zeros.csv", row.names = FALSE)
write.csv(samples_with_zeros_detail, "Zero_values_detail.csv", row.names = FALSE)

cat("\n✓ Files saved:\n")
cat("  - Zero_values_summary.csv\n")
cat("  - Samples_with_zeros.csv\n")
cat("  - Zero_values_detail.csv\n")

# ===== cleaning data =====

# ============================================================================
# REPLACE ZEROS WITH LOD/2 (Half Minimum Detected Value)
# ============================================================================

cat("\n=== REPLACING ZEROS WITH LOD/2 ===\n\n")

# Variables with zeros (from your analysis)
vars_with_zeros <- c("MBP_OP_ug_gdw", "Oslon_OP_ug_gdw", "Olson_IP_ug_gdw", 
                     "MBP_IP_ug_gdw", "SEOC_ug_gdw", "MBC_ug_gdw_bgc", 
                     "SETN_ug_gdw", "MBN_ug_gdw_bgc")

# Calculate minimum non-zero values for each variable
lod_values <- rf_no_outliers %>%
  select(all_of(vars_with_zeros)) %>%
  summarise(across(everything(), 
                   ~min(.[. > 0], na.rm = TRUE))) %>%
  pivot_longer(everything(), 
               names_to = "Variable", 
               values_to = "Min_detected") %>%
  mutate(LOD_half = Min_detected / 2)

cat("LOD/2 replacement values:\n")
print(lod_values)

# Replace zeros with LOD/2
rf_no_outliers_lod <- rf_no_outliers

for(var in vars_with_zeros) {
  # Get the LOD/2 value for this variable
  lod_replacement <- lod_values %>% 
    filter(Variable == var) %>% 
    pull(LOD_half)
  
  # Count zeros before replacement
  n_zeros_before <- sum(rf_no_outliers_lod[[var]] == 0, na.rm = TRUE)
  
  # Replace zeros with LOD/2
  rf_no_outliers_lod[[var]][rf_no_outliers_lod[[var]] == 0] <- lod_replacement
  
  # Confirm replacement
  n_zeros_after <- sum(rf_no_outliers_lod[[var]] == 0, na.rm = TRUE)
  
  cat(sprintf("  %s: %d zeros replaced with %.3f\n", 
              var, n_zeros_before, lod_replacement))
}

# Verify no zeros remain in these variables
cat("\n\nVerification - remaining zeros in treated variables:\n")
remaining_zeros <- rf_no_outliers_lod %>%
  select(all_of(vars_with_zeros)) %>%
  summarise(across(everything(), ~sum(. == 0, na.rm = TRUE)))

print(remaining_zeros)

# Final sample count
cat("\n\nFinal dataset after outlier removal and LOD replacement:\n")
cat("Total samples:", nrow(rf_no_outliers_lod), "\n")

final_strat_lod <- rf_no_outliers_lod %>%
  group_by(Season, Irrigation_trmt_pct) %>%
  summarise(N = n(), .groups = "drop")

print(final_strat_lod)

# Save the cleaned dataset
write.csv(rf_no_outliers_lod, 
          "RF_final_cleaned_outliers_removed_zeros_replaced.csv", 
          row.names = FALSE)

# ============================================================================
# CHECK AND REPLACE NEGATIVE RESPIRATION VALUES
# ============================================================================
cat("\n=== CHECKING FOR NEGATIVE RESPIRATION VALUES ===\n\n")

# Check for negative values in respiration
n_negative <- sum(rf_no_outliers_lod$Resp_rate_mg_C_m2_h < 0, na.rm = TRUE)
cat("Number of negative respiration values:", n_negative, "\n")

if(n_negative > 0) {
  # Show which samples have negative respiration
  negative_resp_samples <- rf_no_outliers_lod %>%
    filter(Resp_rate_mg_C_m2_h < 0) %>%
    select(Sample_cat, Collection_date, Month_name, Season, 
           Irrigation_trmt_pct, Resp_rate_mg_C_m2_h)
  
  cat("\nSamples with negative respiration:\n")
  print(negative_resp_samples)
  
  # Calculate minimum non-negative respiration value
  min_nonneg_resp <- min(rf_no_outliers_lod$Resp_rate_mg_C_m2_h[
    rf_no_outliers_lod$Resp_rate_mg_C_m2_h >= 0], 
    na.rm = TRUE)
  
  replacement_value <- min_nonneg_resp / 2
  
  cat("\n\nMinimum non-negative respiration:", round(min_nonneg_resp, 3), "mg C m⁻² h⁻¹\n")
  cat("Replacement value (half of minimum):", round(replacement_value, 3), "mg C m⁻² h⁻¹\n")
  
  # Replace negative values
  rf_no_outliers_lod$Resp_rate_mg_C_m2_h[
    rf_no_outliers_lod$Resp_rate_mg_C_m2_h < 0] <- replacement_value
  
  cat("\n✓ Negative respiration values replaced with", round(replacement_value, 3), "\n")
  
  # Verify no negative values remain
  n_negative_after <- sum(rf_no_outliers_lod$Resp_rate_mg_C_m2_h < 0, na.rm = TRUE)
  cat("Negative values remaining:", n_negative_after, "\n")
  
} else {
  cat("\n✓ No negative respiration values detected.\n")
}

# ============================================================================
# FINAL DATA SUMMARY AFTER ALL CLEANING STEPS
# ============================================================================
cat("\n\n=== FINAL CLEANED DATASET SUMMARY ===\n\n")
cat("Cleaning steps completed:\n")
cat("  1. Removed outlier samples (4 samples)\n")
cat("  2. Replaced zero values with LOD/2 (8 variables)\n")
cat("  3. Replaced negative respiration with half minimum (", n_negative, " values)\n\n")

cat("Final dataset:\n")
cat("  Total samples:", nrow(rf_no_outliers_lod), "\n")
cat("  Respiration range:", 
    round(min(rf_no_outliers_lod$Resp_rate_mg_C_m2_h, na.rm = TRUE), 3), 
    "to",
    round(max(rf_no_outliers_lod$Resp_rate_mg_C_m2_h, na.rm = TRUE), 3), 
    "mg C m⁻² h⁻¹\n\n")

# Show stratification
final_strat_complete <- rf_no_outliers_lod %>%
  group_by(Season, Irrigation_trmt_pct) %>%
  summarise(
    N = n(),
    Mean_Resp = round(mean(Resp_rate_mg_C_m2_h, na.rm = TRUE), 2),
    Min_Resp = round(min(Resp_rate_mg_C_m2_h, na.rm = TRUE), 2),
    Max_Resp = round(max(Resp_rate_mg_C_m2_h, na.rm = TRUE), 2),
    .groups = "drop"
  )

cat("Stratification with respiration statistics:\n")
print(final_strat_complete)

# ============================================================================
# SAVE FINAL CLEANED DATASET
# ============================================================================
write.csv(rf_no_outliers_lod,
          "RF_final_cleaned_outliers_removed_zeros_replaced.csv",
          row.names = FALSE)

cat("\n✓ Final cleaned dataset saved: RF_final_cleaned_outliers_removed_zeros_replaced.csv\n")
cat("  This file includes:\n")
cat("    • Outliers removed\n")
cat("    • Zeros replaced with LOD/2\n")
cat("    • Negative respiration replaced with half minimum\n\n")

# Optional: Create a summary document of all data cleaning steps
cleaning_summary <- data.frame(
  Step = c("Starting samples (complete cases)", 
           "After outlier removal",
           "After zero replacement",
           "After negative respiration replacement"),
  N_samples = c(nrow(rf_complete),
                nrow(rf_no_outliers),
                nrow(rf_no_outliers),
                nrow(rf_no_outliers_lod)),
  Notes = c("Complete cases from merged data",
            "Removed SM1082, SM1085, SM1087, SM1088",
            paste(length(vars_with_zeros), "variables treated"),
            paste(n_negative, "values replaced"))
)

write.csv(cleaning_summary, "Data_cleaning_summary.csv", row.names = FALSE)
cat("✓ Cleaning summary saved: Data_cleaning_summary.csv\n\n")

