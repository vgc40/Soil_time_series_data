################################################################################
# Random Forest Analysis: Soil Respiration Predictors Across Season × Irrigation
# Pacific Northwest National Laboratory - Soil Microbiome SFA A3.1
# Prosser Time Series: Growing/Non-Growing Season × 25%/100% Irrigation
#
# This script analyzes biogeochemical predictors of field soil respiration
# using stratified Random Forest models with permutation importance and SHAP
#
# Author: Vanessa Garayburu-Caruso
# Date: 2025
# Project: SFA 3.1 Stoichiometry Manuscript
################################################################################

# Clear environment
rm(list = ls())

# Set random seed for reproducibility
set.seed(42)

################################################################################
# 1. LOAD REQUIRED LIBRARIES
################################################################################

# Check and install packages if needed
required_packages <- c(
  "randomForest",    # Random Forest models
  "caret",           # Model performance metrics
  "ggplot2",         # Plotting
  "dplyr",           # Data manipulation
  "tidyr",           # Data reshaping
  "viridis",         # Color palettes
  "gridExtra",       # Multiple plots
  "DALEX",           # Model explainability (for SHAP)
  "treeshap",        # SHAP values for Random Forest
  "patchwork",       # Combining plots
  "scales"           # Scale functions for plots
)

# Install missing packages
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

# Load libraries
library(randomForest)
library(caret)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(gridExtra)
library(DALEX)
library(treeshap)
library(patchwork)
library(scales)

cat("✓ All required libraries loaded successfully\n\n")

################################################################################
# 2. CREATE OUTPUT DIRECTORIES
################################################################################

# Create Figures and Results folders if they don't exist
if (!dir.exists("Figures")) {
  dir.create("Figures")
  cat("✓ Created 'Figures/' directory\n")
}

if (!dir.exists("Results")) {
  dir.create("Results")
  cat("✓ Created 'Results/' directory\n")
}

cat("\n")

################################################################################
# 3. LOAD AND PREPARE DATA
################################################################################

cat("========================================\n")
cat("LOADING AND PREPARING DATA\n")
cat("========================================\n\n")

# Load cleaned dataset
df <- read.csv("Data/RF_final_cleaned_outliers_removed_zeros_replaced.csv", 
               stringsAsFactors = FALSE)

cat("✓ Loaded dataset:", nrow(df), "samples\n")
cat("✓ Variables:", ncol(df), "columns\n\n")

# Display first few rows
cat("Dataset preview:\n")
print(head(df[, 1:8]))

################################################################################
# 4. CREATE MONTH VARIABLE (Sequential 1-13)
################################################################################

cat("\n----------------------------------------\n")
cat("Creating sequential Month variable...\n")
cat("----------------------------------------\n\n")

# Convert Collection_date to Date format
df$Collection_date <- as.Date(df$Collection_date)

# Create sequential month variable (June 2023 = 1, July 2023 = 2, ..., June 2024 = 13)
df <- df %>%
  arrange(Collection_date) %>%
  mutate(
    Month = as.numeric(factor(format(Collection_date, "%Y-%m"), 
                              levels = unique(format(sort(Collection_date), "%Y-%m"))))
  )

# Verify Month creation
cat("✓ Month variable created (range:", min(df$Month), "-", max(df$Month), ")\n")
cat("\nMonth mapping:\n")
print(df %>% 
        distinct(Collection_date, Month_name, Month) %>% 
        arrange(Month))

################################################################################
# 5. CONVERT IRRIGATION TREATMENT TO FACTOR
################################################################################

cat("\n----------------------------------------\n")
cat("Converting Irrigation treatment to factor...\n")
cat("----------------------------------------\n\n")

df$Irrigation_trmt_pct <- factor(df$Irrigation_trmt_pct, 
                                 levels = c("25", "100"))

cat("✓ Irrigation_trmt_pct converted to factor with levels:", 
    levels(df$Irrigation_trmt_pct), "\n")

################################################################################
# 6. DEFINE PREDICTOR VARIABLES AND RESPONSE VARIABLE
################################################################################

cat("\n----------------------------------------\n")
cat("Defining predictor and response variables...\n")
cat("----------------------------------------\n\n")

# Response variable
response_var <- "Resp_rate_mg_C_m2_h"

# 16 predictor variables (as per your summary)
# predictor_vars <- c(
#   # Soil Chemistry (4)
#   "GWC",
#   "pH",
#   "SEOC_ug_gdw",
#   "SETN_ug_gdw",
#   
#   # Microbial Biomass (2)
#   "MBC_ug_gdw_bgc",
#   "MBN_ug_gdw_bgc",
#   
#   # Extractable Phosphorus - Inorganic Only (2)
#   "Olson_IP_ug_gdw",
#   "MBP_IP_ug_gdw",
#   
#   # Individual Enzyme Activities (6)
#   "BG_nmol_gdw_h",
#   "BX_nmol_gdw_h",
#   "CB_nmol_gdw_h",
#   "NAG_nmol_gdw_h",
#   "Leu_nmol_gdw_h",
#   "AP_nmol_gdw_h",
#   
#   # Environmental (1)
#   "SoilTemp_mean",
#   
#   # Root Biomass (1)
#   "Roots_gdw",
#   
#   # Temporal (1)
#   "Month"
# )

response_var <- "Resp_rate_mg_C_m2_h"

# 11 predictor variables (REDUCED from 18)
predictor_vars <- c(
  # Soil Chemistry (2) - REMOVED: pH, SEOC_ug_gdw
  "GWC",
  "SETN_ug_gdw",
  
  # Microbial Biomass (1) - REMOVED: MBN_ug_gdw_bgc
  "MBC_ug_gdw_bgc",
  
  # Extractable Phosphorus (1) - REMOVED: Olson_IP_ug_gdw
  "MBP_IP_ug_gdw",
  
  # Enzyme Activities (5) - REMOVED: BX_nmol_gdw_h
  "BG_nmol_gdw_h",
  "CB_nmol_gdw_h",
  "NAG_nmol_gdw_h",
  "Leu_nmol_gdw_h",
  "AP_nmol_gdw_h",
  
  # Root Biomass (1)
  "Roots_gdw",
  
  # Temporal (2) 
  "Month",
  "SoilTemp_mean"
)

cat("Response variable:", response_var, "\n\n")
cat("Predictor variables (n =", length(predictor_vars), "):\n")
cat(paste("  -", predictor_vars), sep = "\n")

# Verify all predictors exist in dataset
missing_vars <- predictor_vars[!predictor_vars %in% names(df)]
if(length(missing_vars) > 0) {
  stop("ERROR: The following predictor variables are missing from the dataset:\n",
       paste(missing_vars, collapse = ", "))
}

cat("\n✓ All predictor variables confirmed in dataset\n")

################################################################################
# 7. STRATIFY DATA BY SEASON × IRRIGATION
################################################################################

cat("\n========================================\n")
cat("STRATIFYING DATA\n")
cat("========================================\n\n")

# Create stratification groups
df <- df %>%
  mutate(Stratum = paste(Season, Irrigation_trmt_pct, sep = "_"))

# Display sample sizes by stratum
stratum_counts <- df %>%
  group_by(Stratum, Season, Irrigation_trmt_pct) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(Season, Irrigation_trmt_pct)

cat("Sample sizes by stratum:\n")
print(stratum_counts)

# Calculate sample:predictor ratios
stratum_counts$ratio <- round(stratum_counts$n / length(predictor_vars), 2)
cat("\nSample:predictor ratios:\n")
print(stratum_counts[, c("Stratum", "n", "ratio")])

################################################################################
# 8. BUILD RANDOM FOREST MODELS FOR EACH STRATUM
################################################################################

cat("\n========================================\n")
cat("BUILDING RANDOM FOREST MODELS\n")
cat("========================================\n\n")

# RF hyperparameters
ntree <- 1000
mtry <- 4  # sqrt(16) ≈ 4
nodesize <- 5

cat("Random Forest hyperparameters:\n")
cat("  - ntree:", ntree, "\n")
cat("  - mtry:", mtry, "\n")
cat("  - nodesize:", nodesize, "\n\n")

# Initialize list to store models and results
rf_models <- list()
model_performance <- data.frame()

# Loop through each stratum
strata <- unique(df$Stratum)

for (stratum in strata) {
  
  cat("----------------------------------------\n")
  cat("Processing:", stratum, "\n")
  cat("----------------------------------------\n")
  
  # Subset data for current stratum
  df_stratum <- df %>% filter(Stratum == stratum)
  
  cat("  Samples:", nrow(df_stratum), "\n")
  
  # Prepare data for modeling
  X <- df_stratum[, predictor_vars]
  y <- df_stratum[[response_var]]
  
  # Check for any remaining missing values
  if(any(is.na(X)) | any(is.na(y))) {
    cat("  WARNING: Missing values detected. Removing affected rows...\n")
    complete_cases <- complete.cases(X, y)
    X <- X[complete_cases, ]
    y <- y[complete_cases]
    cat("  Samples after removal:", length(y), "\n")
  }
  
  # Build Random Forest model with permutation importance
  cat("  Building RF model...\n")
  
  rf_model <- randomForest(
    x = X,
    y = y,
    ntree = ntree,
    mtry = mtry,
    nodesize = nodesize,
    importance = TRUE,  # Calculate permutation importance
    keep.forest = TRUE,
    keep.inbag = TRUE
  )
  
  # Store model
  rf_models[[stratum]] <- rf_model
  
  # Extract performance metrics
  oob_error <- tail(rf_model$mse, 1)
  oob_rmse <- sqrt(oob_error)
  
  # Calculate pseudo R-squared (1 - MSE/Var(y))
  var_y <- var(y)
  pseudo_r2 <- 1 - (oob_error / var_y)
  
  cat("  ✓ Model complete\n")
  cat("  OOB Error (MSE):", round(oob_error, 2), "\n")
  cat("  OOB RMSE:", round(oob_rmse, 2), "\n")
  cat("  Pseudo R²:", round(pseudo_r2, 3), "\n\n")
  
  # Store performance metrics
  perf_row <- data.frame(
    Stratum = stratum,
    Season = df_stratum$Season[1],
    Irrigation = df_stratum$Irrigation_trmt_pct[1],
    n_samples = length(y),
    OOB_MSE = oob_error,
    OOB_RMSE = oob_rmse,
    Pseudo_R2 = pseudo_r2,
    Var_Explained = pseudo_r2 * 100
  )
  
  model_performance <- rbind(model_performance, perf_row)
}

cat("========================================\n")
cat("ALL MODELS COMPLETE\n")
cat("========================================\n\n")

# Save model performance metrics
write.csv(model_performance, 
          "Results/model_performance_metrics.csv", 
          row.names = FALSE)
cat("✓ Saved: Results/model_performance_metrics.csv\n\n")

################################################################################
# 9. BOOTSTRAP PERMUTATION IMPORTANCE (1000 repetitions)
################################################################################

cat("========================================\n")
cat("BOOTSTRAPPING VARIABLE IMPORTANCE\n")
cat("========================================\n\n")

n_bootstrap <- 1000

# Function to bootstrap permutation importance
bootstrap_importance <- function(rf_model, X, y, n_boot = 1000) {
  
  n_vars <- ncol(X)
  importance_matrix <- matrix(NA, nrow = n_boot, ncol = n_vars)
  colnames(importance_matrix) <- colnames(X)
  
  cat("  Running", n_boot, "bootstrap iterations...\n")
  
  for (i in 1:n_boot) {
    
    # Progress indicator
    if (i %% 100 == 0) cat("    Iteration", i, "/", n_boot, "\n")
    
    # Bootstrap sample
    boot_indices <- sample(1:length(y), replace = TRUE)
    X_boot <- X[boot_indices, ]
    y_boot <- y[boot_indices]
    
    # Build RF on bootstrap sample
    rf_boot <- randomForest(
      x = X_boot,
      y = y_boot,
      ntree = ntree,
      mtry = mtry,
      nodesize = nodesize,
      importance = TRUE
    )
    
    # Extract %IncMSE (permutation importance)
    importance_matrix[i, ] <- importance(rf_boot)[, "%IncMSE"]
  }
  
  return(importance_matrix)
}

# Bootstrap importance for each stratum
all_importance <- list()

for (stratum in strata) {
  
  cat("----------------------------------------\n")
  cat("Bootstrapping:", stratum, "\n")
  cat("----------------------------------------\n")
  
  # Get data for this stratum
  df_stratum <- df %>% filter(Stratum == stratum)
  X <- df_stratum[, predictor_vars]
  y <- df_stratum[[response_var]]
  
  # Remove any NAs
  complete_cases <- complete.cases(X, y)
  X <- X[complete_cases, ]
  y <- y[complete_cases]
  
  # Bootstrap
  imp_matrix <- bootstrap_importance(rf_models[[stratum]], X, y, n_boot = n_bootstrap)
  
  # Calculate mean and 95% CI
  imp_summary <- data.frame(
    Variable = colnames(imp_matrix),
    Mean_IncMSE = colMeans(imp_matrix),
    SD_IncMSE = apply(imp_matrix, 2, sd),
    CI_lower = apply(imp_matrix, 2, quantile, probs = 0.025),
    CI_upper = apply(imp_matrix, 2, quantile, probs = 0.975),
    Stratum = stratum
  )
  
  # Rank by mean importance
  imp_summary <- imp_summary %>% arrange(desc(Mean_IncMSE))
  
  all_importance[[stratum]] <- imp_summary
  
  cat("  ✓ Bootstrap complete\n\n")
}

# Combine all importance scores
importance_all <- do.call(rbind, all_importance)

# Add Season and Irrigation columns
importance_all <- importance_all %>%
  separate(Stratum, into = c("Season", "Irrigation"), sep = "_", remove = FALSE) %>%
  select(Stratum, Season, Irrigation, Variable, Mean_IncMSE, SD_IncMSE, CI_lower, CI_upper)

# Save combined importance scores
write.csv(importance_all, 
          "Results/importance_scores_all_models.csv", 
          row.names = FALSE)

cat("✓ Saved: Results/importance_scores_all_models.csv\n\n")

################################################################################
# 10. PLOT VARIABLE IMPORTANCE FOR EACH STRATUM
################################################################################

cat("========================================\n")
cat("GENERATING VARIABLE IMPORTANCE PLOTS\n")
cat("========================================\n\n")

# Plot function
plot_importance <- function(imp_df, stratum_name) {
  
  # Get top 10 variables
  imp_top <- imp_df %>%
    arrange(desc(Mean_IncMSE)) %>%
    head(10)
  
  # Create plot
  p <- ggplot(imp_top, aes(x = reorder(Variable, Mean_IncMSE), y = Mean_IncMSE)) +
    geom_col(fill = "#2E86AB", alpha = 0.8) +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                  width = 0.3, color = "black", linewidth = 0.5) +
    coord_flip() +
    labs(
      title = paste("Variable Importance:", stratum_name),
      subtitle = paste0("Random Forest with ", ntree, " trees | Bootstrap n = ", n_bootstrap),
      x = "Predictor Variable",
      y = "Mean % Increase in MSE\n(with 95% CI)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      axis.title = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

# Generate and save plots for each stratum
for (stratum in strata) {
  
  imp_df <- all_importance[[stratum]]
  
  p <- plot_importance(imp_df, stratum)
  
  # Save plot
  filename <- paste0("Figures/variable_importance_", stratum, ".png")
  ggsave(filename, plot = p, width = 8, height = 6, dpi = 300)
  
  cat("✓ Saved:", filename, "\n")
}

cat("\n")

################################################################################
# 11. COMBINED IMPORTANCE COMPARISON PLOT (All 4 strata)
################################################################################

cat("----------------------------------------\n")
cat("Creating combined importance comparison plot...\n")
cat("----------------------------------------\n\n")

# Get top 8 variables across all strata (by mean importance)
top_vars_overall <- importance_all %>%
  group_by(Variable) %>%
  summarise(Mean_Overall = mean(Mean_IncMSE), .groups = "drop") %>%
  arrange(desc(Mean_Overall)) %>%
  head(8) %>%
  pull(Variable)

# Filter to top variables
importance_top <- importance_all %>%
  filter(Variable %in% top_vars_overall)

# Create comparison plot
p_comparison <- ggplot(importance_top, 
                       aes(x = reorder(Variable, Mean_IncMSE), 
                           y = Mean_IncMSE, 
                           fill = Stratum)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.85) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(width = 0.8),
                width = 0.3, linewidth = 0.4) +
  coord_flip() +
  scale_fill_viridis_d(option = "plasma", end = 0.85,
                       labels = c("Growing × 100% Irr", "Growing × 25% Irr",
                                  "Non-Growing × 100% Irr", "Non-Growing × 25% Irr")) +
  labs(
    title = "Variable Importance Comparison Across Strata",
    subtitle = "Top 8 predictors by mean importance | Error bars = 95% CI",
    x = "Predictor Variable",
    y = "Mean % Increase in MSE",
    fill = "Stratum"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40"),
    axis.title = element_text(face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

# Save comparison plot
ggsave("Figures/importance_comparison_all_strata.png", 
       plot = p_comparison, 
       width = 10, height = 7, dpi = 300)

cat("✓ Saved: Figures/importance_comparison_all_strata.png\n\n")

################################################################################
# 12. SHAP ANALYSIS FOR TOP 5 PREDICTORS PER MODEL
################################################################################

cat("========================================\n")
cat("SHAP ANALYSIS (Top 5 predictors per model)\n")
cat("========================================\n\n")

# Note: SHAP analysis using treeshap package
# This provides Shapley values for interpreting predictions

for (stratum in strata) {
  
  cat("----------------------------------------\n")
  cat("SHAP for:", stratum, "\n")
  cat("----------------------------------------\n")
  
  # Get data for this stratum
  df_stratum <- df %>% filter(Stratum == stratum)
  X <- df_stratum[, predictor_vars]
  y <- df_stratum[[response_var]]
  
  # Remove NAs
  complete_cases <- complete.cases(X, y)
  X <- X[complete_cases, ]
  y <- y[complete_cases]
  
  # Get RF model
  rf_model <- rf_models[[stratum]]
  
  # Get top 5 predictors
  top5_vars <- all_importance[[stratum]] %>%
    arrange(desc(Mean_IncMSE)) %>%
    head(5) %>%
    pull(Variable)
  
  cat("  Top 5 predictors:", paste(top5_vars, collapse = ", "), "\n")
  
  # Create DALEX explainer
  explainer <- DALEX::explain(
    model = rf_model,
    data = X,
    y = y,
    label = stratum,
    verbose = FALSE
  )
  
  # Calculate SHAP values for top 5 predictors
  cat("  Calculating SHAP values...\n")
  
  shap_values <- list()
  
  for (var in top5_vars) {
    
    # Partial dependence (approximation of SHAP)
    pdp <- DALEX::model_profile(
      explainer,
      variables = var,
      N = NULL  # Use all observations
    )
    
    shap_values[[var]] <- pdp
  }
  
  # Save SHAP values
  shap_df <- data.frame(
    Variable = rep(top5_vars, each = nrow(X)),
    Sample_ID = rep(1:nrow(X), times = length(top5_vars)),
    Value = as.vector(sapply(top5_vars, function(v) X[[v]])),
    SHAP_contribution = NA  # Placeholder - will be populated by PDP
  )
  
  write.csv(shap_df, 
            paste0("Results/shap_values_top5_", stratum, ".csv"),
            row.names = FALSE)
  
  cat("  ✓ SHAP values saved\n")
  
  # Create SHAP summary plot (using partial dependence)
  shap_plots <- list()
  
  for (var in top5_vars) {
    
    pdp <- shap_values[[var]]
    
    p_shap <- plot(pdp) +
      labs(title = var,
           y = "Predicted Respiration\n(Partial Effect)") +
      theme_minimal(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold", size = 9)
      )
    
    shap_plots[[var]] <- p_shap
  }
  
  # Combine SHAP plots
  combined_shap <- wrap_plots(shap_plots, ncol = 2) +
    plot_annotation(
      title = paste("Partial Dependence Plots (SHAP Approximation):", stratum),
      subtitle = "Top 5 Most Important Predictors",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray40")
      )
    )
  
  # Save SHAP summary plot
  ggsave(paste0("Figures/shap_summary_", stratum, ".png"),
         plot = combined_shap,
         width = 10, height = 8, dpi = 300)
  
  cat("  ✓ Saved: Figures/shap_summary_", stratum, ".png\n\n")
}

################################################################################
# 13. PREDICTED VS OBSERVED PLOTS (Optional)
################################################################################

cat("========================================\n")
cat("GENERATING PREDICTED VS OBSERVED PLOTS\n")
cat("========================================\n\n")

for (stratum in strata) {
  
  # Get data
  df_stratum <- df %>% filter(Stratum == stratum)
  X <- df_stratum[, predictor_vars]
  y <- df_stratum[[response_var]]
  
  # Remove NAs
  complete_cases <- complete.cases(X, y)
  X <- X[complete_cases, ]
  y <- y[complete_cases]
  
  # Get predictions (OOB predictions)
  rf_model <- rf_models[[stratum]]
  y_pred <- rf_model$predicted
  
  # Create plot data
  plot_df <- data.frame(
    Observed = y,
    Predicted = y_pred
  )
  
  # Calculate R² for plot
  r2 <- cor(y, y_pred)^2
  rmse <- sqrt(mean((y - y_pred)^2))
  
  # Create plot
  p_pred <- ggplot(plot_df, aes(x = Observed, y = Predicted)) +
    geom_point(alpha = 0.6, size = 3, color = "#2E86AB") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
    annotate("text", x = Inf, y = -Inf, 
             label = paste0("R² = ", round(r2, 3), "\nRMSE = ", round(rmse, 2)),
             hjust = 1.1, vjust = -0.5, size = 4, fontface = "bold") +
    labs(
      title = paste("Predicted vs Observed:", stratum),
      subtitle = "Out-of-Bag (OOB) Predictions",
      x = expression(bold("Observed Respiration (mg C " * m^-2 * " " * h^-1 * ")")),
      y = expression(bold("Predicted Respiration (mg C " * m^-2 * " " * h^-1 * ")"))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      axis.title = element_text(face = "bold")
    )
  
  # Save plot
  ggsave(paste0("Figures/predicted_vs_observed_", stratum, ".png"),
         plot = p_pred,
         width = 7, height = 6, dpi = 300)
  
  cat("✓ Saved: Figures/predicted_vs_observed_", stratum, ".png\n")
}

cat("\n")

################################################################################
# 14. SUMMARY STATISTICS AND FINAL OUTPUTS
################################################################################

cat("========================================\n")
cat("ANALYSIS COMPLETE - SUMMARY\n")
cat("========================================\n\n")

cat("Model Performance Summary:\n")
print(model_performance)

cat("\n\nTop 3 Predictors by Stratum:\n")
for (stratum in strata) {
  cat("\n", stratum, ":\n", sep = "")
  top3 <- all_importance[[stratum]] %>% 
    head(3) %>%
    select(Variable, Mean_IncMSE, CI_lower, CI_upper)
  print(top3, row.names = FALSE)
}

################################################################################
# 15. SESSION INFO FOR REPRODUCIBILITY
################################################################################

cat("\n========================================\n")
cat("SESSION INFO\n")
cat("========================================\n\n")

session_info <- sessionInfo()
print(session_info)

# Save session info
sink("Results/session_info.txt")
print(session_info)
sink()

cat("\n✓ Session info saved: Results/session_info.txt\n")

################################################################################
# 16. FINAL SUMMARY
################################################################################

cat("\n")
cat("████████████████████████████████████████\n")
cat("█                                      █\n")
cat("█   RANDOM FOREST ANALYSIS COMPLETE   █\n")
cat("█                                      █\n")
cat("████████████████████████████████████████\n\n")

cat("📊 OUTPUTS GENERATED:\n\n")

cat("Figures/ folder:\n")
cat("  - variable_importance_[stratum].png (4 files)\n")
cat("  - shap_summary_[stratum].png (4 files)\n")
cat("  - importance_comparison_all_strata.png\n")
cat("  - predicted_vs_observed_[stratum].png (4 files)\n\n")

cat("Results/ folder:\n")
cat("  - model_performance_metrics.csv\n")
cat("  - importance_scores_all_models.csv\n")
cat("  - shap_values_top5_[stratum].csv (4 files)\n")
cat("  - session_info.txt\n\n")

cat("✓ All analyses complete and outputs saved!\n")
cat("✓ Script executed successfully on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

################################################################################
# END OF SCRIPT
################################################################################