rm(list = ls())
# Load required packages
library(tidyverse)
library(vegan)
library(randomForest)
library(caret)
library(ggplot2)
library(reshape2)

# Load your data
soil_data <- read.csv("Data/RF_final_cleaned_outliers_removed_zeros_replaced.csv")

# Define seasons
soil_data <- soil_data %>%
  mutate(
    Season = factor(Season, levels = c("Growing", "Non-Growing")),
    Month_name = factor(Month_name, levels = c("Jun", "Jul", "Aug", "Sep", 
                                               "Oct", "Nov", "Dec", "Jan", 
                                               "Feb", "Mar", "Apr", "May"))
  )

# Create seasonal subsets
growing_data <- soil_data %>% filter(Season == "Growing")
non_growing_data <- soil_data %>% filter(Season == "Non-Growing")

# EXPANDED predictor variables - adding back important variables
predictor_vars_expanded <- c(
  # Soil Chemistry (4) - ADDING BACK pH and SEOC
  "GWC",
  "pH",                    # Added back - can influence enzyme activity
  "SEOC_ug_gdw",          # Added back - available C substrate
  "SETN_ug_gdw",
  
  # Microbial Biomass (2) - ADDING BACK MBN
  "MBC_ug_gdw_bgc",
  "MBN_ug_gdw_bgc",       # Added back - for stoichiometry
  
  # Extractable Phosphorus (2) - ADDING BACK Olson_IP
  "Olson_IP_ug_gdw",      # Added back - available P
  "MBP_IP_ug_gdw",
  
  # Enzyme Activities (7) - ADDING BACK BX and Ala
  "BG_nmol_gdw_h",
  "BX_nmol_gdw_h",        # Added back - hemicellulose degradation
  "CB_nmol_gdw_h",
  "NAG_nmol_gdw_h",
  "Leu_nmol_gdw_h",
  "Ala_nmol_gdw_h",       # Added back - another N-acquiring enzyme
  "AP_nmol_gdw_h",
  
  # Enzyme Ratios (3) - NEW additions
  "EEA_C_N",              # Direct test of C vs N limitation
  "EEA_C_P",              # Direct test of C vs P limitation
  "Sum_Carbon",           # Summed C-acquiring enzymes
  
  # Root Biomass (1)
  "Roots_gdw",
  
  # Temporal (1)
  "SoilTemp_mean"
)

# Response variable
response_var <- "Resp_rate_mg_C_m2_h"

# Function to perform RDA with forward selection - UPDATED
perform_rda_expanded <- function(data, predictors, response, season_name) {
  
  # Remove rows with NA in response or predictors
  clean_data <- data %>%
    select(all_of(c(response, predictors))) %>%
    drop_na()
  
  cat("\n===========================================\n")
  cat("Season:", season_name, "\n")
  cat("Sample size:", nrow(clean_data), "\n")
  cat("Number of predictors:", length(predictors), "\n")
  cat("===========================================\n\n")
  
  # Prepare matrices
  Y <- clean_data[, response, drop = FALSE]
  X <- clean_data[, predictors]
  
  # Standardize predictors (z-scores)
  X_scaled <- scale(X)
  
  # Check for multicollinearity
  cat("Checking for high multicollinearity (VIF > 10):\n")
  # Simple correlation check
  cor_matrix <- cor(X_scaled, use = "complete.obs")
  high_cor <- which(abs(cor_matrix) > 0.9 & cor_matrix != 1, arr.ind = TRUE)
  if(nrow(high_cor) > 0) {
    cat("Warning: High correlations detected between:\n")
    for(i in 1:nrow(high_cor)) {
      cat("  ", rownames(cor_matrix)[high_cor[i,1]], "<->", 
          colnames(cor_matrix)[high_cor[i,2]], ": ", 
          round(cor_matrix[high_cor[i,1], high_cor[i,2]], 3), "\n")
    }
  } else {
    cat("No extreme multicollinearity detected (r > 0.9)\n")
  }
  cat("\n")
  
  # Initial RDA with all variables
  rda_full <- rda(Y ~ ., data = as.data.frame(X_scaled))
  
  cat("Full model R-squared:\n")
  print(RsquareAdj(rda_full))
  cat("\n")
  
  # Forward selection to identify significant predictors
  rda_null <- rda(Y ~ 1, data = as.data.frame(X_scaled))
  
  # Forward selection with permutation tests
  cat("Starting forward selection...\n\n")
  fwd_sel <- ordiR2step(rda_null, 
                        scope = formula(rda_full),
                        direction = "forward",
                        R2scope = TRUE,
                        permutations = how(nperm = 999))
  
  cat("\n--- FORWARD SELECTION SUMMARY ---\n")
  print(summary(fwd_sel))
  
  # Variance partitioning by functional groups
  cat("\n--- VARIANCE PARTITIONING ---\n")
  
  # Define functional groups
  c_enzymes <- c("BG_nmol_gdw_h", "BX_nmol_gdw_h", "CB_nmol_gdw_h", "Sum_Carbon")
  c_enzymes <- c_enzymes[c_enzymes %in% predictors]
  
  n_enzymes <- c("NAG_nmol_gdw_h", "Leu_nmol_gdw_h", "Ala_nmol_gdw_h")
  n_enzymes <- n_enzymes[n_enzymes %in% predictors]
  
  env_factors <- c("GWC", "SoilTemp_mean", "Roots_gdw", "pH")
  env_factors <- env_factors[env_factors %in% predictors]
  
  nutrients <- c("SETN_ug_gdw", "SEOC_ug_gdw", "MBC_ug_gdw_bgc", 
                 "MBN_ug_gdw_bgc", "MBP_IP_ug_gdw", "Olson_IP_ug_gdw")
  nutrients <- nutrients[nutrients %in% predictors]
  
  stoichiometry <- c("EEA_C_N", "EEA_C_P")
  stoichiometry <- stoichiometry[stoichiometry %in% predictors]
  
  # Try variance partitioning if we have enough variables in each group
  if(length(c_enzymes) >= 1 & length(n_enzymes) >= 1 & 
     length(env_factors) >= 1 & length(nutrients) >= 1) {
    
    cat("Performing variance partitioning with functional groups:\n")
    cat("  C-acquiring enzymes:", paste(c_enzymes, collapse=", "), "\n")
    cat("  N-acquiring enzymes:", paste(n_enzymes, collapse=", "), "\n")
    cat("  Environmental factors:", paste(env_factors, collapse=", "), "\n")
    cat("  Nutrient pools:", paste(nutrients, collapse=", "), "\n\n")
    
    varpart_result <- varpart(Y,
                              X_scaled[, c_enzymes, drop = FALSE],
                              X_scaled[, n_enzymes, drop = FALSE],
                              X_scaled[, env_factors, drop = FALSE],
                              X_scaled[, nutrients, drop = FALSE])
    
    print(varpart_result)
    
  } else {
    varpart_result <- NULL
    cat("Insufficient variables for variance partitioning\n")
  }
  
  # Calculate ANOVA for final model
  cat("\n--- ANOVA OF FINAL MODEL ---\n")
  final_anova <- anova.cca(fwd_sel, permutations = 999)
  print(final_anova)
  
  # Return results
  list(
    season = season_name,
    n_samples = nrow(clean_data),
    n_predictors = length(predictors),
    full_rda = rda_full,
    forward_selection = fwd_sel,
    variance_partition = varpart_result,
    final_anova = final_anova,
    selected_variables = labels(terms(fwd_sel))
  )
}

# Run expanded RDA for each season
cat("\n########## GROWING SEASON ##########\n")
rda_growing_expanded <- perform_rda_expanded(growing_data, 
                                             predictor_vars_expanded, 
                                             response_var, 
                                             "Growing Season")

cat("\n\n########## NON-GROWING SEASON ##########\n")
rda_non_growing_expanded <- perform_rda_expanded(non_growing_data, 
                                                 predictor_vars_expanded, 
                                                 response_var, 
                                                 "Non-Growing Season")

# ===================================================================
# VISUALIZATION OF RESULTS
# ===================================================================

# 1. Compare selected variables between seasons
cat("\n\n=== COMPARISON OF SELECTED VARIABLES ===\n")
cat("\nGrowing Season Selected Variables:\n")
print(rda_growing_expanded$selected_variables)

cat("\nNon-Growing Season Selected Variables:\n")
print(rda_non_growing_expanded$selected_variables)

# 2. Plot RDA biplots
par(mfrow = c(1, 2))
plot(rda_growing_expanded$forward_selection, 
     main = paste("Growing Season RDA\nR² =", 
                  round(RsquareAdj(rda_growing_expanded$forward_selection)$adj.r.squared, 3)))

plot(rda_non_growing_expanded$forward_selection, 
     main = paste("Non-Growing Season RDA\nR² =", 
                  round(RsquareAdj(rda_non_growing_expanded$forward_selection)$adj.r.squared, 3)))

# 3. Variance partitioning plots
par(mfrow = c(1, 2))
if(!is.null(rda_growing_expanded$variance_partition)) {
  plot(rda_growing_expanded$variance_partition, 
       digits = 2,
       main = "Variance Partitioning\nGrowing Season")
}

if(!is.null(rda_non_growing_expanded$variance_partition)) {
  plot(rda_non_growing_expanded$variance_partition, 
       digits = 2,
       main = "Variance Partitioning\nNon-Growing Season")
}

# 4. Create comparison table of R² values
comparison_df <- data.frame(
  Season = c("Growing", "Non-Growing"),
  R2_adjusted = c(
    RsquareAdj(rda_growing_expanded$forward_selection)$adj.r.squared,
    RsquareAdj(rda_non_growing_expanded$forward_selection)$adj.r.squared
  ),
  N_samples = c(rda_growing_expanded$n_samples, 
                rda_non_growing_expanded$n_samples),
  N_selected_vars = c(
    length(rda_growing_expanded$selected_variables),
    length(rda_non_growing_expanded$selected_variables)
  )
)

print(comparison_df)

# 5. Extract and compare variable importance
# Get the coefficients from each model
extract_coefficients <- function(rda_result) {
  coef_vals <- scores(rda_result$forward_selection, 
                      choices = 1, 
                      display = "bp", 
                      scaling = 0)
  
  coef_df <- data.frame(
    Variable = rownames(coef_vals),
    Coefficient = coef_vals[,1]
  ) %>%
    arrange(desc(abs(Coefficient)))
  
  return(coef_df)
}

growing_coefs <- extract_coefficients(rda_growing_expanded) %>%
  mutate(Season = "Growing")

non_growing_coefs <- extract_coefficients(rda_non_growing_expanded) %>%
  mutate(Season = "Non-Growing")

# Combine and visualize
library(ggplot2)
all_coefs <- bind_rows(growing_coefs, non_growing_coefs)

ggplot(all_coefs, aes(x = reorder(Variable, abs(Coefficient)), 
                      y = Coefficient, fill = Season)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(title = "RDA Variable Loadings by Season",
       x = "Variable",
       y = "RDA1 Coefficient") +
  theme_minimal() +
  scale_fill_manual(values = c("Growing" = "darkgreen", 
                               "Non-Growing" = "brown"))
