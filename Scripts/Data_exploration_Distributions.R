# Load required libraries
library(tidyverse)
library(moments)  # for skewness and kurtosis
library(gridExtra)  # for arranging plots
df <- read.csv("Data/RF_final_cleaned_outliers_removed_zeros_replaced.csv",
               stringsAsFactors = FALSE)
# Define predictor variables
predictor_vars <- c(
  "GWC",
  "SETN_ug_gdw",
  "MBC_ug_gdw_bgc",
  "MBP_IP_ug_gdw",
  "BG_nmol_gdw_h",
  "CB_nmol_gdw_h",
  "NAG_nmol_gdw_h",
  "Leu_nmol_gdw_h",
  "AP_nmol_gdw_h",
  "Roots_gdw",
  "SoilTemp_mean"
)

# Assuming your data is in 'df' - adjust accordingly
# df <- read.csv("your_data.csv")

# Function to calculate distribution statistics
calculate_distribution_stats <- function(data, vars) {
  stats_list <- list()
  
  for (var in vars) {
    x <- data[[var]][!is.na(data[[var]])]
    
    stats_list[[var]] <- data.frame(
      Variable = var,
      N = length(x),
      Mean = mean(x),
      Median = median(x),
      SD = sd(x),
      Min = min(x),
      Max = max(x),
      Range = max(x) - min(x),
      Skewness = skewness(x),
      Kurtosis = kurtosis(x),
      CV = sd(x) / mean(x) * 100  # Coefficient of Variation
    )
  }
  
  do.call(rbind, stats_list)
}

# Calculate statistics
dist_stats <- calculate_distribution_stats(df, predictor_vars)
print(dist_stats)

# Export statistics
write.csv(dist_stats, "Results/distribution_statistics.csv", row.names = FALSE)

# ============================================================
# VISUAL INSPECTION: Histograms and Q-Q plots
# ============================================================

# Create histogram and Q-Q plot for each variable
for (var in predictor_vars) {
  
  # Create a 2x2 plot layout
  par(mfrow = c(2, 2))
  
  # Get data (remove NAs)
  x <- df[[var]][!is.na(df[[var]])]
  
  # 1. Histogram of original data
  hist(x, main = paste("Histogram:", var), 
       xlab = var, col = "lightblue", breaks = 30)
  
  # 2. Q-Q plot of original data
  qqnorm(x, main = paste("Q-Q Plot:", var))
  qqline(x, col = "red")
  
  # 3. Histogram of log-transformed data (if all values > 0)
  if (all(x > 0)) {
    hist(log(x), main = paste("Histogram: log(", var, ")", sep = ""),
         xlab = paste("log(", var, ")", sep = ""), col = "lightgreen", breaks = 30)
    
    # 4. Q-Q plot of log-transformed data
    qqnorm(log(x), main = paste("Q-Q Plot: log(", var, ")", sep = ""))
    qqline(log(x), col = "red")
  } else {
    plot.new()
    text(0.5, 0.5, "Contains zeros or negative values\nLog transform not applicable", cex = 1.2)
    plot.new()
  }
  
  # Reset plot layout
  par(mfrow = c(1, 1))
  
  # Pause to view (optional - comment out if saving to PDF)
  # readline(prompt = "Press [enter] to continue to next variable")
}

# ============================================================
# ALTERNATIVE: Save all plots to a PDF
# ============================================================

pdf("Figures/distribution_diagnostics.pdf", width = 10, height = 8)

for (var in predictor_vars) {
  par(mfrow = c(2, 2))
  x <- df[[var]][!is.na(df[[var]])]
  
  hist(x, main = paste("Histogram:", var), 
       xlab = var, col = "lightblue", breaks = 30)
  
  qqnorm(x, main = paste("Q-Q Plot:", var))
  qqline(x, col = "red")
  
  if (all(x > 0)) {
    hist(log(x), main = paste("Histogram: log(", var, ")", sep = ""),
         xlab = paste("log(", var, ")", sep = ""), col = "lightgreen", breaks = 30)
    qqnorm(log(x), main = paste("Q-Q Plot: log(", var, ")", sep = ""))
    qqline(log(x), col = "red")
  }
}

dev.off()

# ============================================================
# SHAPIRO-WILK TEST FOR NORMALITY
# ============================================================

# Test normality (works best for n < 5000)
normality_tests <- data.frame()

for (var in predictor_vars) {
  x <- df[[var]][!is.na(df[[var]])]
  
  # Test original data
  if (length(x) < 5000) {
    shapiro_orig <- shapiro.test(x)
    
    # Test log-transformed data if applicable
    if (all(x > 0)) {
      shapiro_log <- shapiro.test(log(x))
      
      normality_tests <- rbind(normality_tests, data.frame(
        Variable = var,
        Original_p_value = shapiro_orig$p.value,
        Log_p_value = shapiro_log$p.value,
        Recommend_Log = ifelse(shapiro_log$p.value > shapiro_orig$p.value, "Yes", "No")
      ))
    } else {
      normality_tests <- rbind(normality_tests, data.frame(
        Variable = var,
        Original_p_value = shapiro_orig$p.value,
        Log_p_value = NA,
        Recommend_Log = "N/A (zeros/negatives)"
      ))
    }
  }
}

print(normality_tests)
write.csv(normality_tests, "Results/normality_tests.csv", row.names = FALSE)

# ============================================================
# TRANSFORMATION RECOMMENDATIONS
# ============================================================

transformation_recommendations <- data.frame()

for (var in predictor_vars) {
  x <- df[[var]][!is.na(df[[var]])]
  skew <- skewness(x)
  
  # Determine transformation based on skewness and presence of zeros
  has_zeros <- any(x <= 0)
  
  recommendation <- case_when(
    has_zeros & abs(skew) > 1 ~ "sqrt or cube root (has zeros/negatives)",
    !has_zeros & skew > 1 ~ "log or sqrt",
    !has_zeros & skew < -1 ~ "square or exponential",
    abs(skew) <= 1 ~ "none needed",
    TRUE ~ "examine visually"
  )
  
  transformation_recommendations <- rbind(transformation_recommendations, data.frame(
    Variable = var,
    Skewness = skew,
    Has_Zeros_Negatives = has_zeros,
    Recommendation = recommendation
  ))
}

print(transformation_recommendations)
write.csv(transformation_recommendations, "Results/transformation_recommendations.csv", row.names = FALSE)

# ============================================================
# COMPARISON PLOTS (ggplot2 version)
# ============================================================

library(ggplot2)

# Create density plots comparing original vs log-transformed
plot_list <- list()

for (var in predictor_vars) {
  x <- df[[var]][!is.na(df[[var]])]
  
  if (all(x > 0)) {
    plot_df <- data.frame(
      Original = x,
      Log_Transformed = log(x)
    ) %>%
      pivot_longer(cols = everything(), names_to = "Type", values_to = "Value")
    
    p <- ggplot(plot_df, aes(x = Value, fill = Type)) +
      geom_density(alpha = 0.5) +
      facet_wrap(~Type, scales = "free") +
      ggtitle(var) +
      theme_minimal()
    
    plot_list[[var]] <- p
  }
}

# Save all ggplots to PDF
pdf("Figures/density_comparisons.pdf", width = 10, height = 6)
for (p in plot_list) {
  print(p)
}
dev.off()

# ============================================================
# SUMMARY REPORT
# ============================================================

cat("\n=== DISTRIBUTION ANALYSIS SUMMARY ===\n\n")
cat("Files created:\n")
cat("  1. distribution_statistics.csv - Basic statistics\n")
cat("  2. normality_tests.csv - Shapiro-Wilk test results\n")
cat("  3. transformation_recommendations.csv - Suggested transformations\n")
cat("  4. distribution_diagnostics.pdf - Visual diagnostics\n")
cat("  5. density_comparisons.pdf - Density plot comparisons\n\n")

cat("Variables with high skewness (|skew| > 1):\n")
print(transformation_recommendations %>% 
        filter(abs(Skewness) > 1) %>% 
        select(Variable, Skewness, Recommendation))