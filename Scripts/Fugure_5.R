rm(list = ls())

# Load required libraries
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(broom.mixed)
library(performance)

# ---- Colors ----
season_cols  <- c("Growing" = "#E69F00", "Non-Growing" = "#56B4E9")
trt_cols     <- c("100" = "#5AAE61", "25" = "#9970AB")

# ---- Data ----
microresp <- read.csv("DataHub/Soil_Phenotypes/SoilSFA_A3.1_MicroResp.csv")
names(microresp)[1] <- "Sample_cat"
microresp$Substrate_NOSC <- as.numeric(microresp$Substrate_NOSC)

biogeochem <- read.csv("Data/RF_cleaned_with_PCA_season.csv")
pca_season_lut <- biogeochem %>%
  select(Sample_cat, Season) %>%
  distinct()

# ---- Prep microresp + join PCA season ----
microresp_data <- microresp %>%
  filter(Timepoint == 5,
         Normalization == "ngCO2-C_mgC_h",
         Substrate != "Water") %>%
  mutate(
    Collection_date = as.Date(Collection_date, format = "%m/%d/%y"),
    Month_num = as.integer(format(Collection_date, "%m")),
    Year = as.integer(format(Collection_date, "%y")),
    # Create distinct labels for Jun23 vs Jun24
    Month = case_when(
      Month_num == 6 & Year == 23 ~ "Jun23",
      Month_num == 6 & Year == 24 ~ "Jun24",
      TRUE ~ month.abb[Month_num]
    ),
    # Set month order with Jun23 and Jun24 at appropriate positions
    Month = factor(Month, 
                   levels = c("Jun23", "Jul", "Aug", "Sep", "Oct", "Nov", 
                              "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun24")),
    Irrigation_trmt_pct = as.character(Irrigation_trmt_pct)
  ) %>%
  left_join(pca_season_lut, by = "Sample_cat") %>%
  mutate(
    Season = if_else(
      is.na(Season) & Month_num == 7, "Growing",
      if_else(is.na(Season) & Month_num == 1, "Non-Growing", as.character(Season))
    ),
    Season = factor(Season, levels = c("Growing", "Non-Growing")),
    Moisture = factor(Irrigation_trmt_pct,
                      levels = c("25", "100"),
                      labels = c("25%", "100%"))
  ) %>%
  filter(!is.na(Substrate_NOSC),
         !is.na(SIR))

# ========================================
# Calculate median respiration for each NOSC value within each Month and Moisture
# ========================================

median_data <- microresp_data %>%
  group_by(Month, Moisture, Substrate_NOSC) %>%
  summarise(
    median_SIR = median(SIR, na.rm = TRUE),
    n_obs = n(),
    .groups = "drop"
  )

# ========================================
# MIXED EFFECTS MODELS
# ========================================

# Model with interaction effects
model_interaction <- lmer(SIR ~ Substrate_NOSC * Season * Moisture +
                            (1|Plot) + (1|Substrate),
                          data = microresp_data)

# Model with additive effects only
model_additive <- lmer(SIR ~ Substrate_NOSC + Season + Moisture +
                         (1|Plot) + (1|Substrate),
                       data = microresp_data)

# Compare models
cat("\n=== Model Comparison ===\n")
print(anova(model_additive, model_interaction))

# View model summary
cat("\n=== Interaction Model Summary ===\n")
print(summary(model_interaction))

# Extract model predictions for plotting
pred_data <- expand.grid(
  Substrate_NOSC = seq(min(microresp_data$Substrate_NOSC, na.rm = TRUE),
                       max(microresp_data$Substrate_NOSC, na.rm = TRUE),
                       length.out = 100),
  Season = c("Growing", "Non-Growing"),
  Moisture = c("25%", "100%"),
  Plot = unique(microresp_data$Plot)[1],
  Substrate = unique(microresp_data$Substrate)[1]
)

pred_data$predicted <- predict(model_interaction, newdata = pred_data, re.form = NA)

# ========================================
# FIGURE 1: Season facets with model predictions
# ========================================

p1 <- ggplot(microresp_data, aes(x = Substrate_NOSC, y = SIR, color = Moisture)) +
  geom_point(alpha = 0.6) +
  geom_line(data = pred_data,
            aes(x = Substrate_NOSC, y = predicted, color = Moisture),
            linewidth = 1.2) +
  facet_wrap(~ Season, ncol = 2) +
  scale_color_manual(values = c("25%" = "#9970AB", "100%" = "#5AAE61")) +
  labs(
    x = "NOSC",
    y = "Mean Carbon Normalized Respiration\n(ng CO2-C per mg substrate C per hour)",
    title = "Substrate-Induced Respiration vs NOSC by Season",
    color = "Irrigation Treatment"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank()
  )

# ========================================
# FIGURE 2: Monthly facets with median values
# ========================================

p_median <- ggplot(median_data, aes(x = Substrate_NOSC, y = median_SIR, color = Moisture)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  facet_wrap(~ Month, ncol = 4) +
  scale_color_manual(values = c("25%" = "#9970AB", "100%" = "#5AAE61")) +
  labs(
    x = "NOSC",
    y = "Median Carbon Normalized Respiration\n(ng CO2-C per mg substrate C per hour)",
    color = "Irrigation Treatment",
    title = "Median Substrate-Induced Respiration vs NOSC by Month"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 10, face = "bold")
  )

# Calculate R² and p-values for median data
median_monthly_stats <- median_data %>%
  group_by(Month, Moisture) %>%
  do(tidy(lm(median_SIR ~ Substrate_NOSC, data = .))) %>%
  filter(term == "Substrate_NOSC") %>%
  select(Month, Moisture, estimate, p.value)

median_monthly_r2 <- median_data %>%
  group_by(Month, Moisture) %>%
  do(glance(lm(median_SIR ~ Substrate_NOSC, data = .))) %>%
  select(Month, Moisture, r.squared)

median_monthly_stats <- left_join(median_monthly_stats, median_monthly_r2, 
                                  by = c("Month", "Moisture"))

# Get the axis limits for text positioning
y_max_median <- max(median_data$median_SIR, na.rm = TRUE)
x_min_median <- min(median_data$Substrate_NOSC, na.rm = TRUE)

# Add statistics to the median plot
p_median_with_stats <- p_median +
  geom_text(data = median_monthly_stats %>% filter(Moisture == "100%"),
            aes(x = x_min_median + 0.2, y = y_max_median * 0.95,
                label = paste0("R²=", sprintf("%.2f", r.squared), 
                               ", p=", sprintf("%.2e", p.value))),
            color = "#5AAE61", hjust = 0, vjust = 1, size = 2.5, 
            fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = median_monthly_stats %>% filter(Moisture == "25%"),
            aes(x = x_min_median + 0.2, y = y_max_median * 0.85,
                label = paste0("R²=", sprintf("%.2f", r.squared), 
                               ", p=", sprintf("%.2e", p.value))),
            color = "#9970AB", hjust = 0, vjust = 1, size = 2.5, 
            fontface = "bold", inherit.aes = FALSE)

# ========================================
# FIGURE 3: Substrate by NOSC
# ========================================

# Get unique substrate-NOSC combinations
substrate_nosc <- microresp_data %>%
  select(Substrate, Substrate_NOSC, Substrate_category) %>%
  distinct() %>%
  arrange(Substrate_category, Substrate_NOSC)

# Bar chart version
p_substrate_nosc <- ggplot(substrate_nosc, 
                           aes(x = reorder(Substrate, Substrate_NOSC), 
                               y = as.factor(Substrate_NOSC), 
                               fill = Substrate_category)) +
  geom_col(color = "black", size = 0.3) +
  coord_flip() +
  labs(
    x = "Substrate",
    y = "NOSC",
    fill = "Substrate Category",
    title = "NOSC Values by Substrate"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 8),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Boxplot version by category
p_category_nosc <- ggplot(substrate_nosc, 
                          aes(x = reorder(Substrate_category, Substrate_NOSC, median), 
                              y = Substrate_NOSC, 
                              fill = Substrate_category)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 2) +
  coord_flip() +
  labs(
    x = "Substrate Category",
    y = "NOSC",
    title = "NOSC Distribution by Substrate Category"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

# Print plots
print(p1)
print(p_median_with_stats)
print(p_substrate_nosc)
print(p_category_nosc)

# Save plots
ggsave("NOSC_respiration_by_season.png", p1, width = 10, height = 6, dpi = 300)
ggsave("NOSC_median_respiration_by_month.png", p_median_with_stats, 
       width = 14, height = 10, dpi = 300)
ggsave("NOSC_by_substrate.png", p_substrate_nosc, 
       width = 10, height = 12, dpi = 300)
ggsave("NOSC_by_category.png", p_category_nosc, 
       width = 10, height = 8, dpi = 300)

# Create faceted plot by month
p2 <- ggplot(microresp_data, aes(x = Substrate_NOSC, y = SIR, color = Moisture)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  facet_wrap(~ Month, ncol = 4) +
  scale_color_manual(values = c("25%" = "#9370DB", "100%" = "#2E8B57")) +
  labs(
    x = "NOSC",
    y = "Mean Carbon Normalized Respiration\n(ng CO2-C per mg substrate C per hour)",
    color = "Irrigation Treatment"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "gray90"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 10, face = "bold")
  )

# Calculate R² and p-values for each month and moisture level
monthly_stats <- microresp_data %>%
  group_by(Month, Moisture) %>%
  do(tidy(lm(SIR ~ Substrate_NOSC, data = .))) %>%
  filter(term == "Substrate_NOSC") %>%
  select(Month, Moisture, estimate, p.value)

monthly_r2 <- microresp_data %>%
  group_by(Month, Moisture) %>%
  do(glance(lm(SIR ~ Substrate_NOSC, data = .))) %>%
  select(Month, Moisture, r.squared)

monthly_stats <- left_join(monthly_stats, monthly_r2, by = c("Month", "Moisture"))

# Get the axis limits from the plot to position text better
y_max <- max(microresp_data$SIR, na.rm = TRUE)
x_min <- min(microresp_data$Substrate_NOSC, na.rm = TRUE)

# Add statistics to the plot with better positioning
p2_with_stats <- p2 +
  geom_text(data = monthly_stats %>% filter(Moisture == "100%"),
            aes(x = x_min + 0.2, y = y_max * 0.95,
                label = paste0("R²=", sprintf("%.2f", r.squared), 
                               ", p=", sprintf("%.2e", p.value))),
            color = "#2E8B57", hjust = 0, vjust = 1, size = 2.5, 
            fontface = "bold", inherit.aes = FALSE) +
  geom_text(data = monthly_stats %>% filter(Moisture == "25%"),
            aes(x = x_min + 0.2, y = y_max * 0.85,
                label = paste0("R²=", sprintf("%.2f", r.squared), 
                               ", p=", sprintf("%.2e", p.value))),
            color = "#9370DB", hjust = 0, vjust = 1, size = 2.5, 
            fontface = "bold", inherit.aes = FALSE)

# Print plots
print(p1)
print(p2_with_stats)
# Extract and display fixed effects from the interaction model
cat("\n=== Fixed Effects ===\n")
fixed_effects <- tidy(model_interaction, effects = "fixed")
print(fixed_effects)

# Calculate marginal and conditional R² for the mixed model
cat("\n=== Model R² Values ===\n")
r2_values <- r2(model_interaction)
print(r2_values)

# Print summary statistics
cat("\n=== NOSC Range by Substrate Category ===\n")
substrate_nosc %>%
  group_by(Substrate_category) %>%
  summarise(
    n_substrates = n(),
    min_NOSC = min(Substrate_NOSC),
    median_NOSC = median(Substrate_NOSC),
    max_NOSC = max(Substrate_NOSC)
  ) %>%
  print()

# Summary of observations per NOSC value
cat("\n=== Summary of observations per NOSC value ===\n")
median_data %>%
  group_by(Substrate_NOSC) %>%
  summarise(
    n_month_moisture_combos = n(),
    median_n_obs = median(n_obs),
    range_n_obs = paste(min(n_obs), "-", max(n_obs))
  ) %>%
  print()


# ==== Check NOSC =====

# Load required library
library(tidyverse)

# Function to parse chemical formula and calculate NOSC
calculate_NOSC <- function(formula) {
  if (is.na(formula) || formula == "na") {
    return(NA)
  }
  
  # Initialize element counts
  C <- 0
  H <- 0
  N <- 0
  O <- 0
  P <- 0
  S <- 0
  
  # Parse formula using regex to extract elements and their counts
  # Pattern matches element symbol followed by optional number
  pattern <- "([A-Z][a-z]?)(\\d*)"
  matches <- gregexpr(pattern, formula, perl = TRUE)
  elements <- regmatches(formula, matches)[[1]]
  
  for (element in elements) {
    # Extract element symbol and count
    parsed <- regmatches(element, regexec("([A-Z][a-z]?)(\\d*)", element))[[1]]
    if (length(parsed) >= 2) {
      symbol <- parsed[2]
      count <- ifelse(parsed[3] == "", 1, as.numeric(parsed[3]))
      
      # Assign to appropriate element
      if (symbol == "C") C <- count
      else if (symbol == "H") H <- count
      else if (symbol == "N") N <- count
      else if (symbol == "O") O <- count
      else if (symbol == "P") P <- count
      else if (symbol == "S") S <- count
    }
  }
  
  # Calculate NOSC using the formula: NOSC = 4 - {4C + H - 3N - 2O + 5P - 2S} / C
  if (C == 0) {
    return(NA)
  }
  
  NOSC <- 4 - ((4 * C + H - 3 * N - 2 * O + 5 * P - 2 * S) / C)
  
  return(NOSC)
}

# Read in the MicroResp data
microresp <- read.csv("DataHub/Soil_Phenotypes/SoilSFA_A3.1_MicroResp.csv")
names(microresp)[1] <- "Sample_cat"
microresp$Substrate_NOSC <- as.numeric(microresp$Substrate_NOSC)

# Get unique substrate-formula-NOSC combinations
substrate_comparison <- microresp %>%
  select(Substrate, Substrate_formula, Substrate_NOSC) %>%
  distinct() %>%
  mutate(
    Calculated_NOSC = sapply(Substrate_formula, calculate_NOSC),
    NOSC_difference = Substrate_NOSC - Calculated_NOSC,
    Percent_difference = ((Substrate_NOSC - Calculated_NOSC) / Substrate_NOSC) * 100
  ) %>%
  arrange(Substrate)

# View the comparison
print(substrate_comparison)

# Summary statistics
cat("\n=== Summary of NOSC Comparison ===\n")
cat("Number of substrates:", nrow(substrate_comparison), "\n")
cat("Number with matching NOSC (within 0.01):", 
    sum(abs(substrate_comparison$NOSC_difference) < 0.01, na.rm = TRUE), "\n")
cat("Mean absolute difference:", 
    mean(abs(substrate_comparison$NOSC_difference), na.rm = TRUE), "\n")
cat("Max absolute difference:", 
    max(abs(substrate_comparison$NOSC_difference), na.rm = TRUE), "\n")

# Identify substrates with large discrepancies
large_discrepancies <- substrate_comparison %>%
  filter(abs(NOSC_difference) > 0.1)

if (nrow(large_discrepancies) > 0) {
  cat("\n=== Substrates with NOSC differences > 0.1 ===\n")
  print(large_discrepancies)
}

# Create a comparison plot
ggplot(substrate_comparison, aes(x = Substrate_NOSC, y = Calculated_NOSC)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_text(data = large_discrepancies, 
            aes(label = Substrate), 
            hjust = -0.1, size = 2.5, check_overlap = TRUE) +
  labs(
    x = "Provided NOSC",
    y = "Calculated NOSC",
    title = "Comparison of Provided vs Calculated NOSC Values"
  ) +
  theme_bw()
