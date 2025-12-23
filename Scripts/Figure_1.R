# Load required libraries
rm(list = ls())
library(tidyverse)
library(patchwork)
library(ggpubr)

# Read data
data <- read.csv('DataHub/Soil_Biogeochemistry/SoilSFA_A3.1_Soil_Biogeochemistry.csv')


# Data cleaning and preparation
plot_data <- data %>%
  mutate(
    Collection_date = as.Date(Collection_date, format = "%m/%d/%y"),
    Month = format(Collection_date, "%b"),
    Year = format(Collection_date, "%y"),
    Month_Year = paste0(Month, Year),
    Irrigation_trmt_pct = factor(Irrigation_trmt_pct, levels = c(25, 100)),
    # Define growing season (April-October when irrigation applied)
    Season = case_when(
      lubridate::month(Collection_date) %in% 4:10 ~ "Growing",
      TRUE ~ "Non-Growing"
    ),
    Season = factor(Season, levels = c("Growing", "Non-Growing")),
    # Convert respiration to rate (mg C m-2 h-1)
    Respiration_rate = CO2_C_60_m2 / 1000
  ) %>%
  # Remove rows with missing respiration data
  filter(!is.na(Respiration_rate))

# Calculate mean and SE for ribbon plot
summary_data <- plot_data %>%
  group_by(Collection_date, Irrigation_trmt_pct) %>%
  summarise(
    mean_resp = mean(Respiration_rate, na.rm = TRUE),
    se_resp = sd(Respiration_rate, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    ymin = mean_resp - se_resp,
    ymax = mean_resp + se_resp
  )

# ===== PANEL A: Time Series with Line and Ribbon =====
# Using colors from your figure: 100% = Green, 25% = Purple
transition_date <- as.Date("2023-11-01")

p_timeseries <- ggplot(summary_data, 
                       aes(x = Collection_date, 
                           y = mean_resp, 
                           color = Irrigation_trmt_pct,
                           fill = Irrigation_trmt_pct,
                           group = Irrigation_trmt_pct)) +
  geom_ribbon(aes(ymin = ymin, ymax = ymax), alpha = 0.3, color = NA) +
  geom_line(linewidth = 1.2) +
  # Add vertical dashed line for season transition
  geom_vline(xintercept = as.numeric(transition_date), 
             linetype = "dashed", 
             color = "darkblue", 
             linewidth = 0.8) +
  # Add text annotation for the seasonal transition
  annotate("text", 
           x = transition_date, 
           y = max(summary_data$ymax, na.rm = TRUE) * 0.95,
           label = " ",
           angle = 90, 
           vjust = -0.5,
           size = 3.5, 
           color = "darkblue",
           fontface = "italic") +
  # Colors matching your figure
  scale_color_manual(values = c("25" = "#9970AB", "100" = "#5AAE61"),
                     labels = c("25% FWC", "100% FWC"),
                     name = "Treatment") +
  scale_fill_manual(values = c("25" = "#9970AB", "100" = "#5AAE61"),
                    labels = c("25% FWC", "100% FWC"),
                    name = "Treatment") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b'%y") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  labs(x = "",
       y = expression(CO[2]*" (mg C m"^-2*" h"^-1*")"),
       title = "") +
  theme_bw() +
  theme(
    legend.position = "right",  # Changed from c(0.85, 0.85) to "right"
    legend.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90")
  )
# ===== PANEL B: Boxplots with Wilcoxon Tests =====

# Wilcoxon test for each season (comparing 25% vs 100%)
wilcox_growing <- plot_data %>%
  filter(Season == "Growing") %>%
  wilcox.test(Respiration_rate ~ Irrigation_trmt_pct, data = .)

wilcox_nongrowing <- plot_data %>%
  filter(Season == "Non-Growing") %>%
  wilcox.test(Respiration_rate ~ Irrigation_trmt_pct, data = .)

print("\n=== Wilcoxon Test Results ===")
print(paste("Growing Season: W =", wilcox_growing$statistic, 
            ", p =", format.pval(wilcox_growing$p.value, digits = 3)))
print(paste("Non-Growing Season: W =", wilcox_nongrowing$statistic,
            ", p =", format.pval(wilcox_nongrowing$p.value, digits = 3)))

# Calculate y position for p-values (above each treatment group)
max_y <- max(plot_data$Respiration_rate, na.rm = TRUE)
y_position_growing <- max(plot_data$Respiration_rate[plot_data$Season == "Growing"], na.rm = TRUE) * 1.15
y_position_nongrowing <- max(plot_data$Respiration_rate[plot_data$Season == "Non-Growing"], na.rm = TRUE) * 1.15

# Create annotation dataframe with positions above each treatment
annotation_df <- data.frame(
  x = c(1, 2),  # Position above 25% and 100%
  y = c(y_position_growing, y_position_growing),
  label = c(
    paste0("Growing: p = ", format.pval(wilcox_growing$p.value, digits = 3)),
    paste0("Non-Growing: p = ", format.pval(wilcox_nongrowing$p.value, digits = 3))
  ),
  Irrigation_trmt_pct = factor(c(25, 100), levels = c(25, 100))
)

# Create boxplot
p_boxplot <- ggplot(plot_data, 
                    aes(x = Irrigation_trmt_pct, 
                        y = Respiration_rate, 
                        fill = Season)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.7, width = 0.6, position = position_dodge(0.8)) +
  geom_jitter(aes(group = Season), alpha = 0.3, size = 1.5, 
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8)) +
  # Using orange/blue colors for seasons
  scale_fill_manual(values = c("Growing" = "#E69F00", "Non-Growing" = "#56B4E9"),
                    name = "Season") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Irrigation Treatment (% FWC)",
       y = expression(CO[2]*" (mg C m"^-2*" h"^-1*")"),
       title = "") +
  # Add p-values above the plot
  annotate("text", 
           x = 1, 
           y = max_y * 1.1,
           label = paste0("Growing: p = ", format.pval(wilcox_growing$p.value, digits = 3)),
           size = 3.5, fontface = "italic", hjust = 0.5) +
  annotate("text", 
           x = 2, 
           y = max_y * 1.1,
           label = paste0("Non-Growing: p = ", format.pval(wilcox_nongrowing$p.value, digits = 3)),
           size = 3.5, fontface = "italic", hjust = 0.5) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )
# ===== COMBINE PANELS =====
combined_plot <- p_timeseries / p_boxplot +
  plot_layout(heights = c(1.2, 1.2)) +
  plot_annotation(
    tag_levels = 'A',
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )
  )

# Display plot
print(combined_plot)

# ===== SAVE OUTPUTS =====
ggsave("Figures/Field_Respiration_Combined.pdf", 
       combined_plot, 
       width = 10, 
       height = 10, 
       units = "in")

ggsave("Figures/Field_Respiration_Combined.png", 
       combined_plot, 
       width = 10, 
       height = 10, 
       units = "in",
       dpi = 300)

# Save statistical results
stat_results <- data.frame(
  Comparison = c("Growing Season", "Non-Growing Season"),
  Test = "Wilcoxon rank-sum test",
  W_statistic = c(wilcox_growing$statistic, wilcox_nongrowing$statistic),
  p_value = c(wilcox_growing$p.value, wilcox_nongrowing$p.value),
  Significant = c(wilcox_growing$p.value < 0.05, wilcox_nongrowing$p.value < 0.05)
)

write.csv(stat_results, "Results/Field_Respiration_Statistics.csv", row.names = FALSE)

