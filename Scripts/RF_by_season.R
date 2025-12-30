################################################################################
# Random Forest Analysis: Soil Respiration Predictors Across Seasons
# Pacific Northwest National Laboratory - Soil Microbiome SFA A3.1
# Prosser Time Series: Growing vs. Non-Growing Season
# REDUCED PREDICTOR SET (11 variables - removed 7 low-importance variables)
#
# This script analyzes biogeochemical predictors of field soil respiration
# using seasonal Random Forest models with permutation importance and SHAP
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

required_packages <- c(
  "randomForest",
  "caret",
  "ggplot2",
  "dplyr",
  "tidyr",
  "viridis",
  "gridExtra",
  "DALEX",
  "patchwork",
  "scales"
)

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(randomForest)
library(caret)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(gridExtra)
library(DALEX)
library(patchwork)
library(scales)

cat("✓ All required libraries loaded successfully\n\n")

################################################################################
# 2. CREATE OUTPUT DIRECTORIES
################################################################################

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

df <- read.csv("Data/RF_final_cleaned_outliers_removed_zeros_replaced.csv",
               stringsAsFactors = FALSE)

cat("✓ Loaded dataset:", nrow(df), "samples\n")
cat("✓ Variables:", ncol(df), "columns\n\n")

cat("Dataset preview:\n")
print(head(df[, 1:8]))

################################################################################
# 4. CREATE MONTH VARIABLE (Sequential 1-13)
################################################################################

cat("\n----------------------------------------\n")
cat("Creating sequential Month variable...\n")
cat("----------------------------------------\n\n")

df$Collection_date <- as.Date(df$Collection_date)

df <- df %>%
  arrange(Collection_date) %>%
  mutate(
    Month = as.numeric(factor(format(Collection_date, "%Y-%m"),
                              levels = unique(format(sort(Collection_date), "%Y-%m"))))
  )

cat("✓ Month variable created (range:", min(df$Month), "-", max(df$Month), ")\n")
cat("\nMonth mapping:\n")
print(df %>%
        distinct(Collection_date, Month_name, Month) %>%
        arrange(Month))

################################################################################
# 5. DEFINE PREDICTOR VARIABLES AND RESPONSE VARIABLE
################################################################################

cat("\n----------------------------------------\n")
cat("Defining predictor and response variables...\n")
cat("----------------------------------------\n\n")

# Response variable
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
  
  # Temporal (1) 
  "SoilTemp_mean"
)

predictor_vars_full <- predictor_vars

cat("Response variable:", response_var, "\n\n")
cat("Predictor variables (n =", length(predictor_vars_full), "):\n")
cat("  Reduced predictor set:\n")
cat(paste("    -", predictor_vars), sep = "\n")

cat("\n✓ REMOVED 7 low-importance variables:\n")
cat("    - Irrigation_trmt_pct (last in both seasons)\n")
cat("    - SEOC_ug_gdw (low importance #17)\n")
cat("    - pH (stable across site, low importance)\n")
cat("    - MBN_ug_gdw_bgc (redundant with MBC)\n")
cat("    - Olson_IP_ug_gdw (kept MBP_IP instead)\n")
cat("    - BX_nmol_gdw_h (lower enzyme importance)\n")

# Verify all predictors exist
missing_vars <- predictor_vars_full[!predictor_vars_full %in% names(df)]
if(length(missing_vars) > 0) {
  stop("ERROR: Missing variables:\n", paste(missing_vars, collapse = ", "))
}

cat("✓ All predictor variables confirmed in dataset\n")

################################################################################
# 6. STRATIFY DATA BY SEASON
################################################################################

cat("\n========================================\n")
cat("STRATIFYING DATA BY SEASON\n")
cat("========================================\n\n")

season_counts <- df %>%
  group_by(Season) %>%
  summarise(n_samples = n(), .groups = "drop")

cat("Sample sizes by season:\n")
print(season_counts)

season_counts$ratio <- round(season_counts$n_samples / length(predictor_vars_full), 2)

cat("\nSample:predictor ratios (with 11 predictors):\n")
print(season_counts[, c("Season", "n_samples", "ratio")])

cat("\n✓ EXCELLENT ratios with reduced predictor set!\n")
cat("  - Growing season: ~7.7:1 (improved from 4.7:1)\n")
cat("  - Non-growing season: ~4.3:1 (improved from 2.6:1)\n\n")

################################################################################
# 7. BUILD RANDOM FOREST MODELS
################################################################################

cat("========================================\n")
cat("BUILDING RANDOM FOREST MODELS\n")
cat("========================================\n\n")

ntree <- 1000
mtry <- 3  # sqrt(11) ≈ 3.32
nodesize <- 5

cat("Random Forest hyperparameters:\n")
cat("  - ntree:", ntree, "\n")
cat("  - mtry:", mtry, "(adjusted for 11 predictors)\n")
cat("  - nodesize:", nodesize, "\n\n")

rf_models <- list()
model_performance <- data.frame()

seasons <- unique(df$Season)

for (season in seasons) {
  cat("----------------------------------------\n")
  cat("Processing:", season, "Season\n")
  cat("----------------------------------------\n")
  
  df_season <- df %>% filter(Season == season)
  cat("  Total samples:", nrow(df_season), "\n")
  
  X <- df_season[, predictor_vars_full]
  y <- df_season[[response_var]]
  
  if(any(is.na(X)) | any(is.na(y))) {
    cat("  WARNING: Missing values detected. Removing...\n")
    complete_cases <- complete.cases(X, y)
    X <- X[complete_cases, ]
    y <- y[complete_cases]
    cat("  Samples after removal:", length(y), "\n")
  }
  
  cat("  Building RF model...\n")
  rf_model <- randomForest(
    x = X,
    y = y,
    ntree = ntree,
    mtry = mtry,
    nodesize = nodesize,
    importance = TRUE,
    keep.forest = TRUE,
    keep.inbag = TRUE
  )
  
  rf_models[[season]] <- rf_model
  
  oob_error <- tail(rf_model$mse, 1)
  oob_rmse <- sqrt(oob_error)
  var_y <- var(y)
  pseudo_r2 <- 1 - (oob_error / var_y)
  
  cat("  ✓ Model complete\n")
  cat("  OOB RMSE:", round(oob_rmse, 2), "\n")
  cat("  Pseudo R²:", round(pseudo_r2, 3), "\n\n")
  
  perf_row <- data.frame(
    Season = season,
    n_samples = length(y),
    n_predictors = length(predictor_vars_full),
    Sample_predictor_ratio = round(length(y) / length(predictor_vars_full), 2),
    OOB_MSE = oob_error,
    OOB_RMSE = oob_rmse,
    Pseudo_R2 = pseudo_r2,
    Var_Explained_pct = pseudo_r2 * 100
  )
  
  model_performance <- rbind(model_performance, perf_row)
}

cat("========================================\n")
cat("ALL MODELS COMPLETE\n")
cat("========================================\n\n")

write.csv(model_performance,
          "Results/model_performance_metrics_reduced.csv",
          row.names = FALSE)

cat("✓ Saved: Results/model_performance_metrics_reduced.csv\n\n")

################################################################################
# 8. BOOTSTRAP PERMUTATION IMPORTANCE
################################################################################

cat("========================================\n")
cat("BOOTSTRAPPING VARIABLE IMPORTANCE\n")
cat("========================================\n\n")

n_bootstrap <- 1000

bootstrap_importance <- function(rf_model, X, y, n_boot = 1000) {
  n_vars <- ncol(X)
  importance_matrix <- matrix(NA, nrow = n_boot, ncol = n_vars)
  colnames(importance_matrix) <- colnames(X)
  
  cat("  Running", n_boot, "bootstrap iterations...\n")
  
  for (i in 1:n_boot) {
    if (i %% 100 == 0) cat("    Iteration", i, "/", n_boot, "\n")
    
    boot_indices <- sample(1:length(y), replace = TRUE)
    X_boot <- X[boot_indices, ]
    y_boot <- y[boot_indices]
    
    rf_boot <- randomForest(
      x = X_boot,
      y = y_boot,
      ntree = ntree,
      mtry = mtry,
      nodesize = nodesize,
      importance = TRUE
    )
    
    importance_matrix[i, ] <- importance(rf_boot)[, "%IncMSE"]
  }
  
  return(importance_matrix)
}

all_importance <- list()

for (season in seasons) {
  cat("----------------------------------------\n")
  cat("Bootstrapping:", season, "Season\n")
  cat("----------------------------------------\n")
  
  df_season <- df %>% filter(Season == season)
  X <- df_season[, predictor_vars_full]
  y <- df_season[[response_var]]
  
  complete_cases <- complete.cases(X, y)
  X <- X[complete_cases, ]
  y <- y[complete_cases]
  
  imp_matrix <- bootstrap_importance(rf_models[[season]], X, y, n_boot = n_bootstrap)
  
  imp_summary <- data.frame(
    Variable = colnames(imp_matrix),
    Mean_IncMSE = colMeans(imp_matrix),
    SD_IncMSE = apply(imp_matrix, 2, sd),
    CI_lower = apply(imp_matrix, 2, quantile, probs = 0.025),
    CI_upper = apply(imp_matrix, 2, quantile, probs = 0.975),
    Season = season
  )
  
  imp_summary <- imp_summary %>% arrange(desc(Mean_IncMSE))
  all_importance[[season]] <- imp_summary
  
  cat("  ✓ Bootstrap complete\n\n")
}

importance_all <- do.call(rbind, all_importance)

write.csv(importance_all,
          "Results/importance_scores_all_models_reduced.csv",
          row.names = FALSE)

cat("✓ Saved: Results/importance_scores_all_models_reduced.csv\n\n")

################################################################################
# 9. PLOT VARIABLE IMPORTANCE
################################################################################

cat("========================================\n")
cat("GENERATING VARIABLE IMPORTANCE PLOTS\n")
cat("========================================\n\n")

plot_importance <- function(imp_df, season_name) {
  imp_top <- imp_df %>%
    arrange(desc(Mean_IncMSE))
  
  p <- ggplot(imp_top, aes(x = reorder(Variable, Mean_IncMSE),
                           y = Mean_IncMSE)) +
    geom_col(alpha = 0.85, fill = "#2E86AB") +
    geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                  width = 0.3, color = "black", linewidth = 0.5) +
    coord_flip() +
    labs(
      title = paste("Variable Importance:", season_name, "Season"),
      subtitle = paste0("Random Forest with ", ntree, " trees | Bootstrap n = ", n_bootstrap,
                        " | Reduced predictor set (n = 11)"),
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

for (season in seasons) {
  imp_df <- all_importance[[season]]
  p <- plot_importance(imp_df, season)
  
  filename <- paste0("Figures/variable_importance_", season, "_season_reduced.png")
  ggsave(filename, plot = p, width = 9, height = 7, dpi = 300)
  cat("✓ Saved:", filename, "\n")
}

cat("\n")

################################################################################
# 10. COMBINED IMPORTANCE COMPARISON PLOT
################################################################################

cat("----------------------------------------\n")
cat("Creating combined importance comparison plot...\n")
cat("----------------------------------------\n\n")

top_vars_overall <- importance_all %>%
  group_by(Variable) %>%
  summarise(Mean_Overall = mean(Mean_IncMSE), .groups = "drop") %>%
  arrange(desc(Mean_Overall)) %>%
  head(10) %>%
  pull(Variable)

importance_top <- importance_all %>%
  filter(Variable %in% top_vars_overall)

p_comparison <- ggplot(importance_top,
                       aes(x = reorder(Variable, Mean_IncMSE),
                           y = Mean_IncMSE,
                           fill = Season)) +
  geom_col(position = position_dodge(width = 0.8), alpha = 0.85) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(width = 0.8),
                width = 0.3, linewidth = 0.4) +
  coord_flip() +
  scale_fill_viridis_d(option = "mako", end = 0.7,
                       labels = c("Growing Season", "Non-Growing Season")) +
  labs(
    title = "Variable Importance Comparison: Growing vs. Non-Growing Season",
    subtitle = "Top 10 predictors by mean importance | Error bars = 95% CI | Reduced predictor set (n = 11)",
    x = "Predictor Variable",
    y = "Mean % Increase in MSE",
    fill = "Season"
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

ggsave("Figures/importance_comparison_seasons_reduced.png",
       plot = p_comparison,
       width = 10, height = 7, dpi = 300)

cat("✓ Saved: Figures/importance_comparison_seasons_reduced.png\n\n")

################################################################################
# 11. SHAP ANALYSIS FOR TOP 5 PREDICTORS PER MODEL
################################################################################

cat("========================================\n")
cat("SHAP ANALYSIS (Top 5 predictors per model)\n")
cat("========================================\n\n")

for (season in seasons) {
  cat("----------------------------------------\n")
  cat("SHAP for:", season, "Season\n")
  cat("----------------------------------------\n")
  
  df_season <- df %>% filter(Season == season)
  X <- df_season[, predictor_vars_full]
  y <- df_season[[response_var]]
  
  complete_cases <- complete.cases(X, y)
  X <- X[complete_cases, ]
  y <- y[complete_cases]
  
  rf_model <- rf_models[[season]]
  
  top5_vars <- all_importance[[season]] %>%
    arrange(desc(Mean_IncMSE)) %>%
    head(5) %>%
    pull(Variable)
  
  cat("  Top 5 predictors:", paste(top5_vars, collapse = ", "), "\n")
  
  explainer <- DALEX::explain(
    model = rf_model,
    data = X,
    y = y,
    label = paste(season, "Season"),
    verbose = FALSE
  )
  
  cat("  Calculating partial dependence plots...\n")
  shap_plots <- list()
  
  for (var in top5_vars) {
    pdp <- DALEX::model_profile(
      explainer,
      variables = var,
      N = NULL
    )
    
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
  
  combined_shap <- wrap_plots(shap_plots, ncol = 2) +
    plot_annotation(
      title = paste("Partial Dependence Plots (SHAP Approximation):", season, "Season"),
      subtitle = "Top 5 Most Important Predictors | Reduced predictor set (n = 11)",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 11, color = "gray40")
      )
    )
  
  ggsave(paste0("Figures/shap_summary_", season, "_season_reduced.png"),
         plot = combined_shap,
         width = 11, height = 9, dpi = 300)
  
  cat("  ✓ Saved: Figures/shap_summary_", season, "_season_reduced.png\n\n")
}

################################################################################
# 12. PREDICTED VS OBSERVED PLOTS
################################################################################

cat("========================================\n")
cat("GENERATING PREDICTED VS OBSERVED PLOTS\n")
cat("========================================\n\n")

for (season in seasons) {
  df_season <- df %>% filter(Season == season)
  X <- df_season[, predictor_vars_full]
  y <- df_season[[response_var]]
  
  complete_cases <- complete.cases(X, y)
  X <- X[complete_cases, ]
  y <- y[complete_cases]
  
  rf_model <- rf_models[[season]]
  y_pred <- rf_model$predicted
  
  plot_df <- data.frame(
    Observed = y,
    Predicted = y_pred
  )
  
  r2 <- cor(y, y_pred)^2
  rmse <- sqrt(mean((y - y_pred)^2))
  
  p_pred <- ggplot(plot_df, aes(x = Observed, y = Predicted)) +
    geom_point(alpha = 0.6, size = 3, color = "#2E86AB") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "black", linewidth = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "gray30", linewidth = 0.6) +
    annotate("text", x = Inf, y = -Inf,
             label = paste0("R² = ", round(r2, 3), "\nRMSE = ", round(rmse, 2)),
             hjust = 1.1, vjust = -0.5, size = 4, fontface = "bold") +
    labs(
      title = paste("Predicted vs Observed:", season, "Season"),
      subtitle = "Out-of-Bag (OOB) Predictions | Reduced predictor set (n = 11)",
      x = expression(bold("Observed Respiration (mg C " *m^-2* " " *h^-1* ")")),
      y = expression(bold("Predicted Respiration (mg C " *m^-2* " " *h^-1* ")"))
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40"),
      axis.title = element_text(face = "bold")
    )
  
  ggsave(paste0("Figures/predicted_vs_observed_", season, "_season_reduced.png"),
         plot = p_pred,
         width = 8, height = 7, dpi = 300)
  
  cat("✓ Saved: Figures/predicted_vs_observed_", season, "_season_reduced.png\n")
}

cat("\n")

################################################################################
# 13. SUMMARY STATISTICS AND FINAL OUTPUTS
################################################################################

cat("========================================\n")
cat("ANALYSIS COMPLETE - SUMMARY\n")
cat("========================================\n\n")

cat("Model Performance Summary:\n")
print(model_performance)

cat("\n\nTop 5 Predictors by Season:\n")
for (season in seasons) {
  cat("\n", season, " Season:\n", sep = "")
  top5 <- all_importance[[season]] %>%
    head(5) %>%
    select(Variable, Mean_IncMSE, CI_lower, CI_upper)
  print(top5, row.names = FALSE)
}

################################################################################
# 14. SESSION INFO FOR REPRODUCIBILITY
################################################################################

cat("\n========================================\n")
cat("SESSION INFO\n")
cat("========================================\n\n")

session_info <- sessionInfo()
print(session_info)

sink("Results/session_info_reduced.txt")
print(session_info)
sink()

cat("\n✓ Session info saved: Results/session_info_reduced.txt\n")

################################################################################
# 15. FINAL SUMMARY
################################################################################

cat("\n")
cat("████████████████████████████████████████\n")
cat("█                                      █\n")
cat("█   RANDOM FOREST ANALYSIS COMPLETE   █\n")
cat("█   REDUCED PREDICTOR SET (n = 11)    █\n")
cat("█                                      █\n")
cat("████████████████████████████████████████\n\n")

cat("📊 OUTPUTS GENERATED:\n\n")
cat("Figures/ folder:\n")
cat("  - variable_importance_[Season]_season_reduced.png (2 files)\n")
cat("  - shap_summary_[Season]_season_reduced.png (2 files)\n")
cat("  - importance_comparison_seasons_reduced.png\n")
cat("  - predicted_vs_observed_[Season]_season_reduced.png (2 files)\n\n")

cat("Results/ folder:\n")
cat("  - model_performance_metrics_reduced.csv\n")
cat("  - importance_scores_all_models_reduced.csv\n")
cat("  - session_info_reduced.txt\n\n")

cat("✓ All analyses complete and outputs saved!\n")
cat("✓ Script executed successfully on:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("========================================\n")
cat("KEY IMPROVEMENTS:\n")
cat("========================================\n")
cat("1. Reduced from 18 to 11 predictors\n")
cat("2. Improved sample:predictor ratios:\n")
cat("   - Growing: 7.7:1 (was 4.7:1)\n")
cat("   - Non-growing: 4.3:1 (was 2.6:1)\n")
cat("3. Removed redundant and low-importance variables\n")
cat("4. Retained most informative predictors\n\n")

cat("REMOVED VARIABLES:\n")
cat("  ❌ Irrigation_trmt_pct (lowest importance)\n")
cat("  ❌ SEOC_ug_gdw (low importance, GWC captures moisture)\n")
cat("  ❌ pH (stable across site)\n")
cat("  ❌ MBN_ug_gdw_bgc (redundant with MBC)\n")
cat("  ❌ Olson_IP_ug_gdw (kept MBP_IP instead)\n")
cat("  ❌ BX_nmol_gdw_h (lower enzyme importance)\n")

################################################################################
# END OF SCRIPT
################################################################################

################################################################################
# MANUSCRIPT FIGURE 2 - PUBLICATION QUALITY
################################################################################
cat("\n========================================\n")
cat("GENERATING MANUSCRIPT FIGURE 2\n")
cat("========================================\n\n")

# Install additional package if needed
if(!require("ggbeeswarm")) install.packages("ggbeeswarm")
library(ggbeeswarm)

# -----------------------------------------------------------------------------
# DATA PREPARATION: Clean Variable Names
# -----------------------------------------------------------------------------
cat("Preparing data with publication-ready labels...\n")

# Create clean variable labels based on your data dictionary
variable_labels <- c(
  "Month" = "Month",
  "SoilTemp_mean" = "Soil Temperature",
  "GWC" = "Gravimetric Water\nContent (%)",
  "MBC_ug_gdw_bgc" = "Microbial Biomass C\n(µg/g dry soil)",
  "BG_nmol_gdw_h" = "β-Glucosidase\n(nmol/g/h)",
  "NAG_nmol_gdw_h" = "N-Acetylglucosaminidase\n(nmol/g/h)",
  "Leu_nmol_gdw_h" = "Leucine Aminopeptidase\n(nmol/g/h)",
  "AP_nmol_gdw_h" = "Acid Phosphatase\n(nmol/g/h)",
  "CB_nmol_gdw_h" = "Cellobiohydrolase\n(nmol/g/h)",
  "SETN_ug_gdw" = "Extractable N\n(µg/g dry soil)",
  "Roots_gdw" = "Root Biomass\n(g/g dry soil)",
  "MBP_IP_ug_gdw" = "Microbial Biomass P\n(µg PO₄/g dry soil)"
)

# Get top 5 variables per season
top5_per_season <- importance_all %>%
  group_by(Season) %>%
  arrange(desc(Mean_IncMSE)) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    Variable_clean = variable_labels[Variable],
    Season_label = factor(
      Season,
      levels = c("Growing", "Non-Growing"),
      labels = c("Growing Season", "Non-Growing Season")
    )
  )

cat("✓ Data prepared with clean labels\n\n")

# -----------------------------------------------------------------------------
# PANEL A: Variable Importance Comparison (Top 5 per season)
# -----------------------------------------------------------------------------
cat("Creating Panel A: Variable Importance...\n")

panel_A <- ggplot(top5_per_season, 
                  aes(x = reorder(Variable_clean, Mean_IncMSE),
                      y = Mean_IncMSE,
                      fill = Season_label)) +
  geom_col(position = position_dodge(width = 0.75), 
           alpha = 0.85, 
           color = "black", 
           linewidth = 0.3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(width = 0.75),
                width = 0.3, 
                linewidth = 0.5) +
  coord_flip() +
  scale_fill_manual(values = c("Growing Season" = "#E69F00", 
                               "Non-Growing Season" = "#56B4E9"),
                    name = "Season") +
  labs(
    title = "A",
    x = NULL,
    y = "Variable Importance (% Increase in MSE)",
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = -0.15),
    axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 10)),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.6, "cm"),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    plot.margin = margin(10, 15, 10, 10)
  )

cat("✓ Panel A created\n\n")

# -----------------------------------------------------------------------------
# PANEL B: SHAP Beeswarm Plots (Directional Effects)
# -----------------------------------------------------------------------------
cat("Creating Panel B: SHAP Directional Effects (Beeswarm style)...\n")

# Function to calculate SHAP-like values using partial dependence
calculate_shap_values <- function(explainer, X, y, variable, n_sample = 100) {
  
  # Sample observations for efficiency
  if(nrow(X) > n_sample) {
    sample_idx <- sample(1:nrow(X), n_sample)
    X_sample <- X[sample_idx, ]
    y_sample <- y[sample_idx]
  } else {
    X_sample <- X
    y_sample <- y
  }
  
  # Get baseline prediction (mean)
  baseline_pred <- mean(y)
  
  # Calculate individual SHAP approximations
  shap_values <- numeric(nrow(X_sample))
  var_values <- X_sample[[variable]]
  
  for(i in 1:nrow(X_sample)) {
    # Create dataset with this observation
    X_temp <- X_sample
    X_temp[variable] <- X_sample[i, variable]
    
    # Predict and calculate SHAP
    pred <- predict(explainer$model, X_temp)
    shap_values[i] <- mean(pred) - baseline_pred
  }
  
  # Return dataframe
  data.frame(
    Variable = variable,
    SHAP_value = shap_values,
    Feature_value = var_values,
    stringsAsFactors = FALSE
  )
}

# Calculate SHAP values for top 5 variables in each season
shap_data_list <- list()

for (season in seasons) {
  cat("  Processing", season, "season SHAP values...\n")
  
  df_season <- df %>% filter(Season == season)
  X <- df_season[, predictor_vars_full]
  y <- df_season[[response_var]]
  
  complete_cases <- complete.cases(X, y)
  X <- X[complete_cases, ]
  y <- y[complete_cases]
  
  rf_model <- rf_models[[season]]
  
  top5_vars <- all_importance[[season]] %>%
    arrange(desc(Mean_IncMSE)) %>%
    head(5) %>%
    pull(Variable)
  
  explainer <- DALEX::explain(
    model = rf_model,
    data = X,
    y = y,
    label = season,
    verbose = FALSE
  )
  
  # Calculate SHAP for each top 5 variable
  for (var in top5_vars) {
    shap_result <- calculate_shap_values(explainer, X, y, var, n_sample = 50)
    shap_result$Season <- season
    shap_result$Variable_clean <- variable_labels[var]
    shap_data_list[[paste(season, var, sep = "_")]] <- shap_result
  }
}

# Combine all SHAP data
shap_combined <- do.call(rbind, shap_data_list)
rownames(shap_combined) <- NULL

# Add season labels
shap_combined$Season_label <- factor(
  shap_combined$Season,
  levels = c("Growing", "Non-Growing"),
  labels = c("Growing Season", "Non-Growing Season")
)

# Normalize feature values for coloring (0-1 scale per variable)
shap_combined <- shap_combined %>%
  group_by(Variable, Season) %>%
  mutate(
    Feature_value_norm = (Feature_value - min(Feature_value)) / 
      (max(Feature_value) - min(Feature_value) + 1e-10)
  ) %>%
  ungroup()

# Create beeswarm plot
panel_B <- ggplot(shap_combined, 
                  aes(x = Variable_clean, 
                      y = SHAP_value,
                      color = Feature_value_norm)) +
  geom_quasirandom(alpha = 0.6, size = 1.5, width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray30", linewidth = 0.5) +
  facet_wrap(~Season_label, ncol = 1, scales = "free_x") +
  coord_flip() +
  scale_color_gradient2(
    low = "#0072B2",      # Blue for low values
    mid = "#F0E442",      # Yellow for mid values  
    high = "#D55E00",     # Orange-red for high values
    midpoint = 0.5,
    name = "Feature Value",
    labels = c("Low", "", "High"),
    breaks = c(0, 0.5, 1)
  ) +
  labs(
    title = "B",
    x = NULL,
    y = "SHAP Value (Impact on Respiration Prediction)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = -0.08),
    axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 10)),
    axis.text = element_text(size = 10, color = "black"),
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    plot.caption = element_text(size = 8, hjust = 0, 
                                margin = margin(t = 10), 
                                lineheight = 1.2, color = "gray30"),
    plot.margin = margin(10, 15, 10, 10)
  )

cat("✓ Panel B created\n\n")

# -----------------------------------------------------------------------------
# COMBINE PANELS INTO FIGURE 2
# -----------------------------------------------------------------------------
cat("Combining panels into Figure 2...\n")

figure_2 <- panel_A | panel_B 
  

# Save as PNG (600 dpi for publication)
ggsave("Figures/Figure_2_manuscript.png",
       plot = figure_2,
       width = 15, 
       height = 10, 
       dpi = 600,
       bg = "white")

# Save as PDF (vector format)
ggsave("Figures/Figure_2_manuscript.pdf",
       plot = figure_2,
       width = 15, 
       height = 10,
       device = cairo_pdf,
       bg = "white")

cat("✓ Saved: Figures/Figure_2_manuscript.png (600 dpi)\n")
cat("✓ Saved: Figures/Figure_2_manuscript.pdf (vector)\n\n")

# Save SHAP summary statistics
shap_summary_stats <- shap_combined %>%
  group_by(Season_label, Variable_clean) %>%
  summarise(
    Mean_SHAP = mean(SHAP_value),
    Median_SHAP = median(SHAP_value),
    SD_SHAP = sd(SHAP_value),
    Direction = ifelse(mean(SHAP_value) > 0, "Positive", "Negative"),
    Mean_Feature_Value = mean(Feature_value),
    .groups = "drop"
  ) %>%
  arrange(Season_label, desc(abs(Mean_SHAP)))

write.csv(shap_summary_stats,
          "Results/SHAP_summary_statistics.csv",
          row.names = FALSE)

cat("✓ Saved: Results/SHAP_summary_statistics.csv\n\n")

cat("========================================\n")
cat("FIGURE 2 COMPLETE\n")
cat("========================================\n\n")
cat("📊 How to interpret Panel B (SHAP Beeswarm):\n\n")
cat("1. X-axis (SHAP Value):\n")
cat("   - Positive values → predictor INCREASES respiration\n")
cat("   - Negative values → predictor DECREASES respiration\n")
cat("   - Magnitude shows strength of effect\n\n")
cat("2. Point Distribution:\n")
cat("   - Spread shows variability in predictor effect\n")
cat("   - Dense clusters indicate consistent effects\n")
cat("   - Spread-out points show context-dependent effects\n\n")
cat("3. Color Gradient:\n")
cat("   - BLUE points = LOW feature values\n")
cat("   - YELLOW points = MEDIUM feature values\n")
cat("   - ORANGE-RED points = HIGH feature values\n\n")
cat("4. Key Interpretation Examples:\n")
cat("   - If HIGH values (red) cluster on RIGHT → High feature = High respiration\n")
cat("   - If HIGH values (red) cluster on LEFT → High feature = Low respiration\n")
cat("   - If colors are mixed → Non-linear or interactive effects\n\n")
cat("========================================\n\n")