# ============================================================================
# 02 Multivariate Analysis
# Q1 Step 1: Multivariate temporal structure in MicroResp substrate use
# ============================================================================
#
# Objective:
#   Test whether MicroResp substrate-use profiles show structured temporal
#   dynamics across the annual cycle and whether those temporal shifts are
#   associated with environmental variables.
#
# Outputs:
#   Figures: Figures/
#   Tables:  Results/
#   Intermediate diagnostics: Results/Step1_diagnostics/

rm(list = ls())

# ----------------------------------------------------------------------------
# A. Load packages
# ----------------------------------------------------------------------------

required_packages <- c("vegan", "tidyverse", "viridis", "ggrepel", "patchwork")
missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

library(vegan)
library(tidyverse)
library(viridis)
library(ggrepel)
library(patchwork)

set.seed(42)

theme_set(theme_bw(base_size = 12))

# ----------------------------------------------------------------------------
# File paths and output directories
# ----------------------------------------------------------------------------

microresp_candidates <- c(
  "DataHub March/Soil_Phenotypes/SoilMicrobiomeSFA_prosser-time-series_MicroResp.csv",
  "DataHub/Soil_Phenotypes/SoilSFA_A3.1_MicroResp.csv",
  "SoilMicrobiomeSFA_prosser-time-series_MicroResp.csv"
)

# Use full clean dataset even if it has some NAs
env_candidates <- c(
  "Data/Data_merged_cleaned.csv"
)

microresp_file <- microresp_candidates[file.exists(microresp_candidates)][1]
env_file <- env_candidates[file.exists(env_candidates)][1]

if (is.na(microresp_file)) {
  stop(
    "Could not find MicroResp input. Checked: ",
    paste(microresp_candidates, collapse = "; ")
  )
}

if (is.na(env_file)) {
  stop(
    "Could not find cleaned environmental input. Checked: ",
    paste(env_candidates, collapse = "; ")
  )
}

dir.create("Figures", showWarnings = FALSE, recursive = TRUE)
dir.create("Results", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Step1_diagnostics", showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# B. Load and filter MicroResp data
# ----------------------------------------------------------------------------

microresp_raw <- read.csv(microresp_file, check.names = FALSE)

# The code below handles a common typo in the "Irrigation_trmt_cat" column name that was present
if ("Irrigaiton_trmt_cat" %in% names(microresp_raw) &&
    !"Irrigation_trmt_cat" %in% names(microresp_raw)) {
  microresp_raw <- microresp_raw %>%
    rename(Irrigation_trmt_cat = Irrigaiton_trmt_cat)
}

# Filter microresp data to 5h and normalized per C
microresp_filtered <- microresp_raw %>%
  mutate(
    Timepoint = suppressWarnings(as.numeric(Timepoint))
  ) %>%
  filter(
    Timepoint == 5,
    Normalization == "ngCO2-C_mgC_h",
    Substrate != "Water"
  ) %>%
  mutate(
    SIR = suppressWarnings(as.numeric(SIR))
  )

sir_parse_failures <- microresp_filtered %>%
  filter(is.na(SIR)) %>%
  select(Sample_cat, Substrate, SIR)

write_csv(
  sir_parse_failures,
  "Results/Step1_diagnostics/sir_parse_failures_after_filtering.csv"
)

n_unique_substrates <- n_distinct(microresp_filtered$Substrate)
unique_substrates <- sort(unique(microresp_filtered$Substrate))

cat("\n=== MICRORESP FILTER CHECK ===\n")
cat("Rows after Timepoint/Normalization/Substrate filtering:",
    nrow(microresp_filtered), "\n")
cat("Unique substrates remaining:", n_unique_substrates, "\n")
cat("Expected unique substrates: 31\n")
cat("Substrates:\n")
print(unique_substrates)

write_csv(
  tibble(Substrate = unique_substrates),
  "Results/Step1_diagnostics/unique_substrates_after_filtering.csv"
)

if (n_unique_substrates != 31) {
  warning(
    "Expected 31 unique substrates after filtering, but found ",
    n_unique_substrates,
    ". Review Results/Step1_diagnostics/unique_substrates_after_filtering.csv."
  )
}

# ----------------------------------------------------------------------------
# C. Build community matrix
# ----------------------------------------------------------------------------

sample_substrate_counts <- microresp_filtered %>%
  count(Sample_cat, Substrate, name = "n_records") %>%
  filter(n_records > 1)

write_csv(
  sample_substrate_counts,
  "Results/Step1_diagnostics/duplicate_sample_substrate_values.csv"
)

if (nrow(sample_substrate_counts) > 0) {
  cat("\nDuplicate sample-by-substrate values found; averaging duplicates.\n")
  cat("Duplicate combinations:", nrow(sample_substrate_counts), "\n")
} else {
  cat("\nNo duplicate sample-by-substrate values found.\n")
}

microresp_sample_substrate <- microresp_filtered %>%
  group_by(Sample_cat, Substrate) %>%
  summarise(
    SIR = if_else(all(is.na(SIR)), NA_real_, mean(SIR, na.rm = TRUE)),
    n_records = n(),
    .groups = "drop"
  )

community_wide <- microresp_sample_substrate %>%
  select(Sample_cat, Substrate, SIR) %>%
  pivot_wider(names_from = Substrate, values_from = SIR) %>%
  arrange(Sample_cat)

substrate_cols <- setdiff(names(community_wide), "Sample_cat")

initial_na_cells <- sum(is.na(community_wide[substrate_cols]))
initial_total_cells <- nrow(community_wide) * length(substrate_cols)
initial_na_pct <- initial_na_cells / initial_total_cells * 100

sample_missing_counts <- community_wide %>%
  mutate(n_missing_substrates = rowSums(is.na(across(all_of(substrate_cols))))) %>%
  select(Sample_cat, n_missing_substrates) %>%
  arrange(desc(n_missing_substrates), Sample_cat)

samples_removed_for_missing <- sample_missing_counts %>%
  filter(n_missing_substrates > 3)

write_csv(
  sample_missing_counts,
  "Results/Step1_diagnostics/sample_missing_substrate_counts.csv"
)
write_csv(
  samples_removed_for_missing,
  "Results/Step1_diagnostics/samples_removed_more_than_3_missing_substrates.csv"
)

community_wide <- community_wide %>%
  filter(!Sample_cat %in% samples_removed_for_missing$Sample_cat)

na_cells_after_sample_filter <- sum(is.na(community_wide[substrate_cols]))
total_cells_after_sample_filter <- nrow(community_wide) * length(substrate_cols)
na_pct_after_sample_filter <- na_cells_after_sample_filter /
  total_cells_after_sample_filter * 100

cat("\n=== COMMUNITY MATRIX MISSINGNESS ===\n")
cat("Initial NA cells:", initial_na_cells, "of", initial_total_cells,
    "(", round(initial_na_pct, 2), "% )\n")
cat("Samples removed for >3 missing substrates:",
    nrow(samples_removed_for_missing), "\n")
cat("NA cells after sample filtering:", na_cells_after_sample_filter, "of",
    total_cells_after_sample_filter,
    "(", round(na_pct_after_sample_filter, 2), "% )\n")

if (na_pct_after_sample_filter > 5) {
  stop(
    "More than 5% of community matrix cells are NA after removing samples ",
    "with >3 missing substrates. Review MicroResp missingness before NMDS."
  )
}

if (na_cells_after_sample_filter > 0) {
  community_wide[substrate_cols] <- community_wide[substrate_cols] %>%
    mutate(across(everything(), ~replace_na(., 0)))
  cat("Replaced", na_cells_after_sample_filter,
      "remaining NA cells with zero.\n")
} else {
  cat("No NA replacement was needed.\n")
}

nonzero_substrate_cols <- substrate_cols[
  colSums(as.matrix(community_wide[substrate_cols]), na.rm = TRUE) > 0
]

removed_all_zero_substrates <- setdiff(substrate_cols, nonzero_substrate_cols)

write_csv(
  tibble(Substrate = removed_all_zero_substrates),
  "Results/Step1_diagnostics/substrates_removed_all_zero.csv"
)

community_matrix <- community_wide %>%
  select(all_of(nonzero_substrate_cols)) %>%
  as.matrix()

rownames(community_matrix) <- community_wide$Sample_cat
mode(community_matrix) <- "numeric"

cat("\n=== FINAL MICRORESP MATRIX BEFORE ENVIRONMENTAL MATCHING ===\n")
cat("Samples:", nrow(community_matrix), "\n")
cat("Substrates:", ncol(community_matrix), "\n")

# ----------------------------------------------------------------------------
# D. Prepare environmental data and merge metadata
# ----------------------------------------------------------------------------

env_raw <- read.csv(env_file, check.names = FALSE)

required_env_cols <- c(
  "Sample_cat",
  "Collection_date",
  "Month",
  "Month_name",
  "Plot",
  "Irrigation_trmt_pct",
  "GWC",
  "pH",
  "SEOC_ug_gdw",
  "SETN_ug_gdw",
  "MBC_ug_gdw_bgc",
  "MBN_ug_gdw_bgc",
  "SoilTemp_mean",
  "Roots_gdw"
)

missing_env_cols <- setdiff(required_env_cols, names(env_raw))

if (length(missing_env_cols) > 0) {
  stop(
    "Environmental data are missing required columns: ",
    paste(missing_env_cols, collapse = ", ")
  )
}

env_subset <- env_raw %>%
  filter(Sample_cat %in% rownames(community_matrix))

community_without_env <- tibble(Sample_cat = rownames(community_matrix)) %>%
  anti_join(env_subset %>% select(Sample_cat), by = "Sample_cat")

env_without_community <- env_raw %>%
  anti_join(tibble(Sample_cat = rownames(community_matrix)), by = "Sample_cat") %>%
  select(Sample_cat)

write_csv(
  community_without_env,
  "Results/Step1_diagnostics/community_samples_without_environmental_data.csv"
)
write_csv(
  env_without_community,
  "Results/Step1_diagnostics/environmental_samples_without_community_data.csv"
)

if (nrow(community_without_env) > 0) {
  cat("\nDropping", nrow(community_without_env),
      "community samples without matching cleaned environmental data.\n")
}

common_samples <- intersect(rownames(community_matrix), env_subset$Sample_cat)
community_matrix <- community_matrix[common_samples, , drop = FALSE]

env_data <- env_subset %>%
  mutate(sample_order = match(Sample_cat, common_samples)) %>%
  arrange(sample_order) %>%
  select(-sample_order)

if (!identical(rownames(community_matrix), env_data$Sample_cat)) {
  stop("Community matrix rows and environmental metadata are not aligned.")
}

env_data <- env_data %>%
  mutate(
    Collection_date = as.Date(Collection_date, format = "%Y-%m-%d"),
    time_numeric = as.numeric(Collection_date - min(Collection_date, na.rm = TRUE)),
    year_month = floor_date(Collection_date, "month"),
    SEOC_SETN = if_else(SETN_ug_gdw == 0, NA_real_, SEOC_ug_gdw / SETN_ug_gdw),
    Irrigation_factor = factor(as.character(Irrigation_trmt_pct),
                               levels = c("25", "100"))
  )

month_lookup <- env_data %>%
  distinct(year_month) %>%
  arrange(year_month) %>%
  mutate(
    month_order = row_number(),
    month_label = format(year_month, "%b'%y")
  )

env_data <- env_data %>%
  left_join(month_lookup, by = "year_month") %>%
  mutate(
    month_factor = factor(month_order, levels = sort(unique(month_order))),
    Season_quarter = case_when(
      Month %in% c(12, 1, 2) ~ "Winter",
      Month %in% c(3, 4, 5) ~ "Spring",
      Month %in% c(6, 7, 8) ~ "Summer",
      Month %in% c(9, 10, 11) ~ "Fall",
      TRUE ~ NA_character_
    ),
    Season_quarter = factor(
      Season_quarter,
      levels = c("Winter", "Spring", "Summer", "Fall")
    )
  )

write_csv(
  month_lookup,
  "Results/Step1_diagnostics/month_order_lookup.csv"
)

cat("\n=== MATRIX AND ENVIRONMENTAL ALIGNMENT ===\n")
cat("Matched samples:", nrow(community_matrix), "\n")
cat("Substrates:", ncol(community_matrix), "\n")
cat("Month order lookup:\n")
print(month_lookup)

# ----------------------------------------------------------------------------
# E. NMDS ordination
# ----------------------------------------------------------------------------

set.seed(42)
nmds_k2 <- metaMDS(
  community_matrix,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  trace = FALSE
)

nmds_result <- nmds_k2
nmds_k <- 2
nmds_k3_stress <- NA_real_
nmds_k3_used <- FALSE

if (nmds_k2$stress > 0.20) {
  set.seed(42)
  nmds_k3 <- metaMDS(
    community_matrix,
    distance = "bray",
    k = 3,
    trymax = 200,
    autotransform = FALSE,
    trace = FALSE
  )

  nmds_k3_stress <- nmds_k3$stress
  stress_drop <- nmds_k2$stress - nmds_k3$stress

  if (stress_drop > 0.05) {
    nmds_result <- nmds_k3
    nmds_k <- 3
    nmds_k3_used <- TRUE
  }
}

nmds_stress <- nmds_result$stress
nmds_tries <- if ("tries" %in% names(nmds_result)) nmds_result$tries else NA_integer_
nmds_converged <- if ("converged" %in% names(nmds_result)) {
  nmds_result$converged
} else {
  !is.na(nmds_stress)
}

cat("\n=== NMDS SUMMARY ===\n")
cat("Selected k:", nmds_k, "\n")
cat("Final stress:", round(nmds_stress, 4), "\n")
cat("Converged:", nmds_converged, "\n")
cat("Tries needed:", nmds_tries, "\n")
if (!is.na(nmds_k3_stress)) {
  cat("k=3 stress:", round(nmds_k3_stress, 4), "\n")
  cat("k=3 selected for downstream analyses:", nmds_k3_used, "\n")
}

png("Figures/Step1_shepard_plot.png", width = 7, height = 6, units = "in", res = 300)
stressplot(nmds_result, main = paste0("Shepard plot, stress = ", round(nmds_stress, 3)))
dev.off()

nmds_scores <- scores(nmds_result, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample_cat") %>%
  left_join(env_data, by = "Sample_cat")

if (!all(c("NMDS1", "NMDS2") %in% names(nmds_scores))) {
  stop("NMDS site scores do not contain NMDS1 and NMDS2.")
}

# ----------------------------------------------------------------------------
# F. Environmental vector fitting
# ----------------------------------------------------------------------------

envfit_vars <- c(
  "SoilTemp_mean",
  "GWC",
  "SEOC_ug_gdw",
  "SETN_ug_gdw",
  "SEOC_SETN",
  "MBC_ug_gdw_bgc",
  "MBN_ug_gdw_bgc",
  "pH"
)

envfit_continuous_df <- env_data %>%
  select(all_of(envfit_vars)) %>%
  as.data.frame()

set.seed(42)
envfit_continuous <- envfit(
  nmds_result,
  envfit_continuous_df,
  permutations = 999,
  na.rm = TRUE
)

set.seed(42)
envfit_irrigation <- envfit(
  nmds_result,
  env_data %>% select(Irrigation_factor) %>% as.data.frame(),
  permutations = 999,
  na.rm = TRUE
)

envfit_scores <- scores(envfit_continuous, display = "vectors") %>%
  as.data.frame() %>%
  rownames_to_column("Variable")

envfit_results <- envfit_scores %>%
  mutate(
    R2 = envfit_continuous$vectors$r[Variable],
    p_value = envfit_continuous$vectors$pvals[Variable]
  ) %>%
  arrange(p_value)

if (!"NMDS2" %in% names(envfit_results)) {
  stop("Envfit vector scores do not contain NMDS2.")
}

envfit_significant <- envfit_results %>%
  filter(p_value < 0.05)

write_csv(envfit_results, "Results/Step1_envfit_results.csv")
write_csv(envfit_significant, "Results/Step1_envfit_significant.csv")

cat("\n=== ENVFIT CONTINUOUS SUMMARY ===\n")
print(envfit_continuous)
cat("\n=== ENVFIT IRRIGATION SUMMARY ===\n")
print(envfit_irrigation)

# ----------------------------------------------------------------------------
# G. PERMANOVA
# ----------------------------------------------------------------------------

set.seed(42)
permanova_continuous_time <- adonis2(
  community_matrix ~ time_numeric + Irrigation_factor,
  data = env_data,
  method = "bray",
  permutations = 999,
  by = "margin"
)

set.seed(42)
permanova_month_factor <- adonis2(
  community_matrix ~ month_factor + Irrigation_factor,
  data = env_data,
  method = "bray",
  permutations = 999,
  by = "margin"
)

permanova_continuous_df <- permanova_continuous_time %>%
  as.data.frame() %>%
  rownames_to_column("Term")

permanova_month_df <- permanova_month_factor %>%
  as.data.frame() %>%
  rownames_to_column("Term")

write_csv(permanova_continuous_df, "Results/Step1_PERMANOVA_continuous_time.csv")
write_csv(permanova_month_df, "Results/Step1_PERMANOVA_month_factor.csv")

cat("\n=== PERMANOVA: CONTINUOUS TIME ===\n")
print(permanova_continuous_df)
cat("\n=== PERMANOVA: MONTH FACTOR ===\n")
print(permanova_month_df)

# ----------------------------------------------------------------------------
# H. Dispersion homogeneity checks
# ----------------------------------------------------------------------------

bray_dist <- vegdist(community_matrix, method = "bray")

set.seed(42)
betadisper_irrigation <- betadisper(bray_dist, env_data$Irrigation_factor)
set.seed(42)
betadisper_irrigation_perm <- permutest(betadisper_irrigation, permutations = 999)

set.seed(42)
betadisper_season <- betadisper(bray_dist, env_data$Season_quarter)
set.seed(42)
betadisper_season_perm <- permutest(betadisper_season, permutations = 999)

betadisper_irrigation_df <- betadisper_irrigation_perm$tab %>%
  as.data.frame() %>%
  rownames_to_column("Term")

betadisper_season_df <- betadisper_season_perm$tab %>%
  as.data.frame() %>%
  rownames_to_column("Term")

write_csv(betadisper_irrigation_df, "Results/Step1_betadisper_irrigation.csv")
write_csv(betadisper_season_df, "Results/Step1_betadisper_season.csv")

betadisper_irrigation_p <- betadisper_irrigation_df %>%
  filter(Term == "Groups") %>%
  pull(`Pr(>F)`)

betadisper_season_p <- betadisper_season_df %>%
  filter(Term == "Groups") %>%
  pull(`Pr(>F)`)

cat("\n=== BETADISPER CHECKS ===\n")
cat("Irrigation dispersion p-value:", betadisper_irrigation_p, "\n")
cat("Season-quarter dispersion p-value:", betadisper_season_p, "\n")
if (!is.na(betadisper_irrigation_p) && betadisper_irrigation_p < 0.05) {
  cat("Caveat: irrigation dispersion is significantly heterogeneous.\n")
}
if (!is.na(betadisper_season_p) && betadisper_season_p < 0.05) {
  cat("Caveat: season-quarter dispersion is significantly heterogeneous.\n")
}

# ----------------------------------------------------------------------------
# I. Figure 1A: NMDS colored by sampling month
# ----------------------------------------------------------------------------

clean_vector_labels <- c(
  SoilTemp_mean = "Soil Temp",
  GWC = "GWC",
  SEOC_ug_gdw = "SEOC",
  SETN_ug_gdw = "SETN",
  SEOC_SETN = "SEOC:SETN",
  MBC_ug_gdw_bgc = "MBC",
  MBN_ug_gdw_bgc = "MBN",
  pH = "pH"
)

arrow_multiplier <- 0.6

envfit_plot_vectors <- envfit_significant %>%
  mutate(
    xend = NMDS1 * arrow_multiplier,
    yend = NMDS2 * arrow_multiplier,
    Label = recode(Variable, !!!clean_vector_labels)
  )

monthly_centroids <- nmds_scores %>%
  group_by(month_order, month_label) %>%
  summarise(
    NMDS1 = mean(NMDS1, na.rm = TRUE),
    NMDS2 = mean(NMDS2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(month_order)

month_breaks <- c(1, 3, 5, 7, 9, 11, 13)
month_break_labels <- month_lookup %>%
  filter(month_order %in% month_breaks) %>%
  arrange(month_order)

figure1a <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
  geom_path(
    data = monthly_centroids,
    aes(group = 1),
    linetype = "dashed",
    color = "black",
    linewidth = 0.7
  ) +
  geom_point(
    aes(color = month_order),
    size = 2.5,
    alpha = 0.8
  ) +
  geom_point(
    data = monthly_centroids,
    aes(fill = month_order),
    shape = 21,
    color = "black",
    size = 4,
    stroke = 0.9
  ) +
  geom_segment(
    data = envfit_plot_vectors,
    aes(x = 0, y = 0, xend = xend, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.18, "cm")),
    linewidth = 0.6,
    color = "black"
  ) +
  geom_text_repel(
    data = envfit_plot_vectors,
    aes(x = xend, y = yend, label = Label),
    inherit.aes = FALSE,
    size = 3.5,
    min.segment.length = 0
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = paste0("Stress = ", round(nmds_stress, 3)),
    hjust = 1.05,
    vjust = -0.7,
    size = 3.5
  ) +
  scale_color_viridis_c(
    name = "Month",
    option = "D",
    breaks = month_break_labels$month_order,
    labels = month_break_labels$month_label
  ) +
  scale_fill_viridis_c(option = "D", guide = "none") +
  labs(
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank()
  )

ggsave("Figures/Figure1A_NMDS_temporal.png", figure1a,
       width = 8, height = 6.5, dpi = 300)
ggsave("Figures/Figure1A_NMDS_temporal.pdf", figure1a,
       width = 8, height = 6.5)

# ----------------------------------------------------------------------------
# J. Figure 1B: NMDS colored by irrigation treatment
# ----------------------------------------------------------------------------

irrigation_row <- permanova_continuous_df %>%
  filter(Term == "Irrigation_factor")

irrigation_r2 <- irrigation_row %>% pull(R2)
irrigation_p <- irrigation_row %>% pull(`Pr(>F)`)

irrigation_annotation <- paste0(
  "PERMANOVA irrigation\n",
  "R2 = ", round(irrigation_r2, 3),
  ", p = ", signif(irrigation_p, 3)
)

figure1b <- ggplot(
  nmds_scores,
  aes(x = NMDS1, y = NMDS2, color = Irrigation_factor, shape = Irrigation_factor)
) +
  stat_ellipse(
    aes(group = Irrigation_factor),
    type = "t",
    level = 0.95,
    linewidth = 0.8
  ) +
  geom_point(size = 2.5, alpha = 0.8) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    label = irrigation_annotation,
    hjust = 1.05,
    vjust = -0.35,
    size = 3.5
  ) +
  scale_color_manual(
    name = "Irrigation",
    values = c("100" = "#0072B2", "25" = "#D55E00"),
    breaks = c("100", "25"),
    labels = c("100%", "25%")
  ) +
  scale_shape_manual(
    name = "Irrigation",
    values = c("100" = 16, "25" = 17),
    breaks = c("100", "25"),
    labels = c("100%", "25%")
  ) +
  labs(
    x = "NMDS1",
    y = "NMDS2"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank()
  )

ggsave("Figures/Figure1B_NMDS_irrigation.png", figure1b,
       width = 7, height = 6, dpi = 300)
ggsave("Figures/Figure1B_NMDS_irrigation.pdf", figure1b,
       width = 7, height = 6)

# ----------------------------------------------------------------------------
# K. Combined Figure 1
# ----------------------------------------------------------------------------

figure1_combined <- figure1a + figure1b +
  plot_annotation(tag_levels = "A")

ggsave("Figures/Figure1_combined.png", figure1_combined,
       width = 15, height = 6.5, dpi = 300)
ggsave("Figures/Figure1_combined.pdf", figure1_combined,
       width = 15, height = 6.5)

# ----------------------------------------------------------------------------
# L. Summary diagnostics table
# ----------------------------------------------------------------------------

time_row <- permanova_continuous_df %>%
  filter(Term == "time_numeric")

time_r2 <- time_row %>% pull(R2)
time_p <- time_row %>% pull(`Pr(>F)`)

significant_envfit_names <- envfit_significant %>%
  pull(Variable) %>%
  paste(collapse = "; ")

summary_diagnostics <- tibble(
  Metric = c(
    "NMDS stress",
    "Number of dimensions (k)",
    "Number of samples in community matrix",
    "Number of substrates in community matrix",
    "PERMANOVA R2 for time_numeric",
    "PERMANOVA p-value for time_numeric",
    "PERMANOVA R2 for Irrigation_factor",
    "PERMANOVA p-value for Irrigation_factor",
    "Number of significant envfit vectors (p < 0.05)",
    "Names of significant envfit vectors",
    "Betadisper p-value for irrigation",
    "Betadisper p-value for season (quarterly)"
  ),
  Value = c(
    as.character(round(nmds_stress, 4)),
    as.character(nmds_k),
    as.character(nrow(community_matrix)),
    as.character(ncol(community_matrix)),
    as.character(round(time_r2, 4)),
    as.character(signif(time_p, 4)),
    as.character(round(irrigation_r2, 4)),
    as.character(signif(irrigation_p, 4)),
    as.character(nrow(envfit_significant)),
    significant_envfit_names,
    as.character(signif(betadisper_irrigation_p, 4)),
    as.character(signif(betadisper_season_p, 4))
  )
)

write_csv(summary_diagnostics, "Results/Step1_summary_diagnostics.csv")

cat("\n=== STEP 1 SUMMARY DIAGNOSTICS ===\n")
print(summary_diagnostics)

# ----------------------------------------------------------------------------
# M. Save intermediate outputs for downstream steps
# ----------------------------------------------------------------------------

write.csv(
  community_matrix,
  "Results/Step1_community_matrix.csv",
  row.names = TRUE
)

write_csv(
  nmds_scores,
  "Results/Step1_nmds_scores_with_metadata.csv"
)

write_csv(
  env_data,
  "Results/Step1_env_data_prepared.csv"
)

write_csv(
  community_wide,
  "Results/Step1_diagnostics/community_matrix_before_environmental_matching.csv"
)

# ----------------------------------------------------------------------------
# Interpretation guidance
# ----------------------------------------------------------------------------

stress_quality <- case_when(
  nmds_stress < 0.15 ~ "good",
  nmds_stress <= 0.20 ~ "acceptable",
  TRUE ~ "poor"
)

envfit_text <- if (nrow(envfit_significant) > 0) {
  envfit_significant %>%
    mutate(
      Label = recode(Variable, !!!clean_vector_labels),
      Summary = paste0(Label, " (R2 = ", round(R2, 3),
                       ", p = ", signif(p_value, 3), ")")
    ) %>%
    pull(Summary) %>%
    paste(collapse = "; ")
} else {
  "No envfit vectors were significant at p < 0.05."
}

interpretation_lines <- c(
  paste0("NMDS stress was ", round(nmds_stress, 3), ", which is considered ",
         stress_quality, " for a ", nmds_k, "-dimensional solution."),
  paste0("The temporal NMDS figure should be interpreted by looking for an ordered ",
         "viridis color gradient and a coherent dashed centroid trajectory across months."),
  paste0("Significant environmental vectors: ", envfit_text),
  paste0("Irrigation PERMANOVA R2 was ", round(irrigation_r2, 3),
         " with p = ", signif(irrigation_p, 3),
         "; visual overlap is shown in Figure 1B."),
  paste0("Time_numeric PERMANOVA R2 was ", round(time_r2, 3),
         " with p = ", signif(time_p, 3),
         ", compared with irrigation R2 = ", round(irrigation_r2, 3), ".")
)

writeLines(interpretation_lines, "Results/Step1_interpretation_summary.txt")

cat("\n=== INTERPRETATION SUMMARY ===\n")
cat(paste(interpretation_lines, collapse = "\n"), "\n")

cat("\nQ1 Step 1 multivariate analysis complete.\n")
