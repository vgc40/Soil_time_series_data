# =============================================================================
# Figure 3: transcript seasonality and targeted MicroResp comparisons
# =============================================================================
# Linear workflow. Heatmap functions use abundance/prevalence, and Panel C
# gene sets are curated independently of temporal patterns.

rm(list = ls())

required_packages <- c(
  "tidyverse", "lubridate", "vegan", "permute", "pheatmap", "patchwork"
)
missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]
if (length(missing_packages) > 0) {
  stop("Install required packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(vegan)
  library(permute)
  library(pheatmap)
  library(patchwork)
})

set.seed(42)

season_colors <- c(
  "Plant-active, Apr-Oct" = "#2E8B57",
  "Plant-dormant, Nov-Mar" = "#5B6F8C"
)
figure_dir <- file.path("Figures", "manuscript_final", "Figure3")
result_dir <- file.path("Results", "manuscript_final", "Figure3")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

theme_manuscript <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.key = element_blank(),
    strip.text = element_text(face = "bold")
  )

# ---- Match 100% treatment transcript samples to field metadata -------------

field_data <- read_csv(
  file.path("Data", "Data_merged_cleaned.csv"),
  show_col_types = FALSE
) %>%
  filter(as.character(Irrigation_trmt_pct) == "100") %>%
  transmute(
    Sample_cat,
    Sampling_Date = as.Date(Collection_date),
    Plot_ID = factor(Plot),
    Block = factor(Block),
    Season = factor(
      if_else(
        month(Sampling_Date) %in% 4:10,
        "Plant-active, Apr-Oct",
        "Plant-dormant, Nov-Mar"
      ),
      levels = names(season_colors)
    )
  ) %>%
  distinct()

brite_c_wide <- read_csv(
  file.path("Data", "brite_C_expression.csv"),
  show_col_types = FALSE
) %>%
  rename(brite_C = 1) %>%
  filter(!is.na(brite_C), brite_C != "")

brite_d_wide <- read_csv(
  file.path("Data", "brite_D_expression.csv"),
  show_col_types = FALSE
) %>%
  rename(brite_D = 1) %>%
  filter(!is.na(brite_D), brite_D != "")

transcript_columns <- names(brite_c_wide)[-1]
transcript_metadata <- tibble(
  Sample_mt = transcript_columns,
  Sample_cat = str_remove(transcript_columns, "_MT$")
) %>%
  left_join(field_data, by = "Sample_cat") %>%
  arrange(Sampling_Date, Plot_ID)

unmatched_transcripts <- transcript_metadata %>%
  filter(is.na(Sampling_Date) | is.na(Plot_ID) | is.na(Season))
write_csv(
  unmatched_transcripts,
  file.path(result_dir, "Figure3_unmatched_transcript_samples.csv")
)

if (nrow(unmatched_transcripts) > 0) {
  stop("Some transcript samples do not match the 100% field metadata.")
}
if (!identical(names(brite_d_wide)[-1], transcript_columns)) {
  stop("BRITE C and BRITE D matrices have different sample columns.")
}
if (anyDuplicated(transcript_metadata$Sample_cat) > 0) {
  stop("Transcript sample identifiers are duplicated.")
}

# Apply the repository's curated non-soil/human pathway exclusions.
human_exclusion <- read.delim(
  file.path(
    "Data", "pathwaystoremove", "kegg_human_exclude_2026-07-13.txt"
  ),
  check.names = FALSE
) %>%
  transmute(pathway_code = str_extract(as.character(path_ko), "ko\\d{5}")) %>%
  filter(!is.na(pathway_code)) %>%
  distinct()

brite_c_long <- brite_c_wide %>%
  mutate(pathway_code = str_extract(brite_C, "ko\\d{5}")) %>%
  filter(
    !is.na(pathway_code),
    !pathway_code %in% human_exclusion$pathway_code
  ) %>%
  select(-pathway_code) %>%
  pivot_longer(
    cols = all_of(transcript_columns),
    names_to = "Sample_mt",
    values_to = "normalized_expression"
  ) %>%
  left_join(transcript_metadata, by = "Sample_mt") %>%
  mutate(normalized_expression = as.numeric(normalized_expression))

invalid_brite_c <- brite_c_long %>%
  filter(
    is.na(normalized_expression) |
      !is.finite(normalized_expression) |
      normalized_expression < 0
  )
if (nrow(invalid_brite_c) > 0) {
  write_csv(
    invalid_brite_c,
    file.path(result_dir, "Figure3_invalid_BRITE_C_values.csv")
  )
  stop("BRITE C expression contains invalid values.")
}

transcript_metadata %>%
  count(Sampling_Date, Season, name = "n_transcript_samples") %>%
  arrange(Sampling_Date) %>%
  write_csv(
    file.path(result_dir, "Figure3_transcript_samples_by_date.csv")
  )

# ---- Panel A: overall transcript seasonality ------------------------------

brite_c_sample_table <- brite_c_long %>%
  select(Sample_mt, brite_C, normalized_expression) %>%
  pivot_wider(
    names_from = brite_C,
    values_from = normalized_expression,
    values_fill = 0
  ) %>%
  arrange(match(Sample_mt, transcript_metadata$Sample_mt))

brite_c_matrix <- brite_c_sample_table %>%
  column_to_rownames("Sample_mt") %>%
  as.matrix()
brite_c_matrix <- brite_c_matrix[
  ,
  colSums(brite_c_matrix) > 0,
  drop = FALSE
]

if (any(rowSums(brite_c_matrix) == 0)) {
  stop("At least one transcript sample has zero BRITE C expression.")
}

transcript_metadata <- transcript_metadata %>%
  arrange(match(Sample_mt, rownames(brite_c_matrix)))

transcript_nmds <- metaMDS(
  brite_c_matrix,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  trace = FALSE
)

permutation_design <- how(
  nperm = 999,
  blocks = transcript_metadata$Plot_ID
)
season_permanova <- adonis2(
  brite_c_matrix ~ Season,
  data = transcript_metadata,
  permutations = permutation_design,
  method = "bray"
)
season_dispersion <- betadisper(
  vegdist(brite_c_matrix, method = "bray"),
  transcript_metadata$Season
)
season_dispersion_test <- permutest(
  season_dispersion,
  permutations = permutation_design
)

permanova_table <- season_permanova %>%
  as.data.frame() %>%
  rownames_to_column("term")
dispersion_table <- season_dispersion_test$tab %>%
  as.data.frame() %>%
  rownames_to_column("term")
write_csv(
  permanova_table,
  file.path(result_dir, "Figure3A_BRITE_C_season_PERMANOVA.csv")
)
write_csv(
  dispersion_table,
  file.path(result_dir, "Figure3A_BRITE_C_season_betadisper.csv")
)

nmds_scores <- scores(transcript_nmds, display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("Sample_mt") %>%
  as_tibble() %>%
  left_join(transcript_metadata, by = "Sample_mt")

nmds_centroids <- nmds_scores %>%
  group_by(Sampling_Date, Season) %>%
  summarise(
    NMDS1 = mean(NMDS1),
    NMDS2 = mean(NMDS2),
    .groups = "drop"
  ) %>%
  arrange(Sampling_Date)

permanova_f <- permanova_table %>%
  slice(1) %>%
  pull(F)
permanova_r2 <- permanova_table %>%
  slice(1) %>%
  pull(R2)
permanova_p <- permanova_table %>%
  slice(1) %>%
  pull(`Pr(>F)`)
dispersion_p <- dispersion_table %>%
  slice(1) %>%
  pull(`Pr(>F)`)

panel_a <- ggplot(nmds_scores, aes(NMDS1, NMDS2, color = Season)) +
  geom_path(
    data = nmds_centroids,
    aes(x = NMDS1, y = NMDS2, group = 1),
    inherit.aes = FALSE,
    color = "grey30",
    linewidth = 0.75,
    arrow = arrow(length = unit(0.10, "inches"), type = "closed")
  ) +
  geom_point(aes(shape = Plot_ID), size = 2.8, alpha = 0.85) +
  geom_point(
    data = nmds_centroids,
    aes(x = NMDS1, y = NMDS2, fill = Season),
    inherit.aes = FALSE,
    shape = 21,
    color = "black",
    size = 4
  ) +
  geom_text(
    data = nmds_centroids,
    aes(x = NMDS1, y = NMDS2, label = format(Sampling_Date, "%b")),
    inherit.aes = FALSE,
    nudge_y = 0.035,
    size = 3.2,
    check_overlap = TRUE
  ) +
  scale_color_manual(values = season_colors, drop = FALSE) +
  scale_fill_manual(values = season_colors, drop = FALSE) +
  labs(
    title = "A",
    subtitle = paste0(
      "BRITE C Bray-Curtis NMDS; stress = ",
      formatC(transcript_nmds$stress, digits = 3, format = "f"),
      "\nPERMANOVA R\u00B2 = ", formatC(permanova_r2, digits = 3, format = "f"),
      ", F = ", formatC(permanova_f, digits = 2, format = "f"),
      ", p = ", format.pval(permanova_p, digits = 2),
      "; dispersion p = ", format.pval(dispersion_p, digits = 2)
    ),
    x = "NMDS1",
    y = "NMDS2",
    color = "Season",
    fill = "Season",
    shape = "Plot"
  ) +
  theme_manuscript +
  theme(legend.position = "bottom") +
  guides(fill = "none", color = guide_legend(nrow = 2))

# ---- Panel B: compact BRITE C temporal heatmap ----------------------------

# Major functions are selected without temporal p-values: at least 80%
# prevalence, followed by the 40 highest mean-expression functions.
brite_c_ranking <- brite_c_long %>%
  group_by(brite_C) %>%
  summarise(
    mean_expression = mean(normalized_expression),
    detection_prevalence = mean(normalized_expression > 0),
    n_detected_samples = sum(normalized_expression > 0),
    .groups = "drop"
  ) %>%
  mutate(
    biologically_interpretable = !is.na(str_extract(brite_C, "ko\\d{5}")),
    prevalence_eligible = biologically_interpretable &
      detection_prevalence >= 0.80
  ) %>%
  arrange(desc(prevalence_eligible), desc(mean_expression)) %>%
  group_by(prevalence_eligible) %>%
  mutate(eligible_abundance_rank = row_number()) %>%
  ungroup() %>%
  mutate(
    included_in_heatmap = prevalence_eligible &
      eligible_abundance_rank <= 40
  )

write_csv(
  brite_c_ranking,
  file.path(result_dir, "Figure3B_BRITE_C_major_function_ranking.csv")
)

major_brite_c <- brite_c_ranking %>%
  filter(included_in_heatmap) %>%
  pull(brite_C)
if (length(major_brite_c) < 10) {
  stop("Too few BRITE C functions passed the heatmap criteria.")
}

brite_c_heatmap_long <- brite_c_long %>%
  filter(brite_C %in% major_brite_c) %>%
  group_by(brite_C) %>%
  mutate(
    function_mean = mean(normalized_expression),
    function_sd = sd(normalized_expression),
    expression_zscore = (normalized_expression - function_mean) / function_sd
  ) %>%
  ungroup() %>%
  group_by(brite_C, Sampling_Date, Season) %>%
  summarise(
    mean_expression_zscore = mean(expression_zscore),
    .groups = "drop"
  )

if (any(!is.finite(brite_c_heatmap_long$mean_expression_zscore))) {
  stop("A selected BRITE C function has zero variance.")
}

heatmap_dates <- sort(unique(transcript_metadata$Sampling_Date))
heatmap_date_labels <- format(heatmap_dates, "%b %Y")

brite_c_heatmap_table <- brite_c_heatmap_long %>%
  mutate(
    date_label = factor(
      format(Sampling_Date, "%b %Y"),
      levels = heatmap_date_labels
    ),
    function_label = str_remove(brite_C, "^\\d{5}\\s+") %>%
      str_remove("\\s+\\[PATH:ko\\d{5}\\]$") %>%
      str_wrap(width = 38)
  ) %>%
  select(function_label, date_label, mean_expression_zscore) %>%
  pivot_wider(
    names_from = date_label,
    values_from = mean_expression_zscore
  )

brite_c_heatmap_matrix <- brite_c_heatmap_table %>%
  column_to_rownames("function_label") %>%
  as.matrix()

row_clustering <- hclust(
  as.dist(
    1 - cor(
      t(brite_c_heatmap_matrix),
      use = "pairwise.complete.obs",
      method = "pearson"
    )
  ),
  method = "average"
)

heatmap_annotation <- transcript_metadata %>%
  distinct(Sampling_Date, Season) %>%
  arrange(Sampling_Date) %>%
  transmute(
    date_label = format(Sampling_Date, "%b %Y"),
    Season
  ) %>%
  column_to_rownames("date_label")

panel_b_heatmap <- pheatmap(
  brite_c_heatmap_matrix,
  cluster_rows = row_clustering,
  cluster_cols = FALSE,
  annotation_col = heatmap_annotation,
  annotation_colors = list(Season = season_colors),
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  breaks = seq(-2.5, 2.5, length.out = 102),
  fontsize = 8,
  fontsize_row = 6.7,
  fontsize_col = 8,
  border_color = NA,
  angle_col = 45,
  main = "B   Major BRITE C functions (expression z-score)",
  legend_breaks = c(-2, 0, 2),
  silent = TRUE
)
panel_b <- wrap_elements(full = panel_b_heatmap$gtable)

# ---- Panel C: targeted MicroResp-transcript comparisons -------------------

target_gene_table <- read_csv(
  file.path("Data", "merged_chitin_chitosan_methyl_oleate_gene_categories.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    target_ko = str_extract(str_to_upper(KO), "K\\d{5}"),
    evidence_text = str_to_lower(
      paste(Zip_category, Image_category, sep = " | ")
    ),
    target_family = case_when(
      str_detect(evidence_text, "methyl oleate|oleate") ~ "Methyl oleate",
      str_detect(
        evidence_text,
        "chitin|chitosan|glucosamine|glcn|glcna"
      ) ~ "Chitin / chitosan / N-acetylglucosamine",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(target_family), !is.na(target_ko)) %>%
  distinct(target_family, target_ko, .keep_all = TRUE)

target_family_levels <- c(
  "Chitin / chitosan / N-acetylglucosamine",
  "Methyl oleate"
)
target_gene_table <- target_gene_table %>%
  mutate(target_family = factor(target_family, levels = target_family_levels))

brite_d_target_long <- brite_d_wide %>%
  mutate(target_ko = str_extract(brite_D, "^K\\d{5}")) %>%
  inner_join(
    target_gene_table %>%
      select(target_family, target_ko),
    by = "target_ko"
  ) %>%
  pivot_longer(
    cols = all_of(transcript_columns),
    names_to = "Sample_mt",
    values_to = "normalized_expression"
  ) %>%
  mutate(normalized_expression = as.numeric(normalized_expression)) %>%
  group_by(target_family, target_ko, Sample_mt) %>%
  summarise(
    normalized_expression = sum(normalized_expression, na.rm = TRUE),
    .groups = "drop"
  )

target_ko_audit <- target_gene_table %>%
  select(
    target_family, target_ko, KEGG_Name_1, EC_number,
    Zip_category, Image_category, Zip_evidence, Image_evidence, Source
  ) %>%
  left_join(
    brite_d_target_long %>%
      group_by(target_family, target_ko) %>%
      summarise(
        present_in_brite_d = TRUE,
        detected_in_any_sample = any(normalized_expression > 0),
        .groups = "drop"
      ),
    by = c("target_family", "target_ko")
  ) %>%
  mutate(
    present_in_brite_d = replace_na(present_in_brite_d, FALSE),
    detected_in_any_sample = replace_na(detected_in_any_sample, FALSE),
    selection_basis = paste(
      "Curated biological annotation supplied before this analysis;",
      "temporal significance was not used."
    )
  )
write_csv(
  target_ko_audit,
  file.path(result_dir, "Figure3C_curated_transcript_gene_sets.csv")
)

target_expression <- brite_d_target_long %>%
  group_by(target_family, Sample_mt) %>%
  summarise(
    gene_set_expression = sum(normalized_expression),
    n_curated_kos_observed = n_distinct(target_ko),
    .groups = "drop"
  ) %>%
  group_by(target_family) %>%
  mutate(standardized_value = as.numeric(scale(gene_set_expression))) %>%
  ungroup() %>%
  left_join(transcript_metadata, by = "Sample_mt")

if (any(!is.finite(target_expression$standardized_value))) {
  stop("A curated transcript gene set has zero variance.")
}

microresp_substrate_map <- tribble(
  ~Substrate, ~target_family,
  "Chitin", "Chitin / chitosan / N-acetylglucosamine",
  "Chitosan", "Chitin / chitosan / N-acetylglucosamine",
  "N.Acetyl.D.glucosamine", "Chitin / chitosan / N-acetylglucosamine",
  "Glucosamine_hydrochloride", "Chitin / chitosan / N-acetylglucosamine",
  "Methyl_oleate", "Methyl oleate"
) %>%
  mutate(target_family = factor(target_family, levels = target_family_levels))

microresp_target_data <- read_csv(
  file.path(
    "DataHub March",
    "Soil_Phenotypes",
    "SoilMicrobiomeSFA_prosser-time-series_MicroResp.csv"
  ),
  col_types = cols(.default = col_character()),
  na = c("", "NA")
) %>%
  rename(Irrigation_trmt_cat = any_of("Irrigaiton_trmt_cat")) %>%
  mutate(
    Sampling_Date = mdy(Collection_date),
    Timepoint = suppressWarnings(as.numeric(Timepoint)),
    SIR = suppressWarnings(as.numeric(SIR))
  ) %>%
  filter(
    Irrigation_trmt_pct == "100",
    Timepoint == 5,
    Normalization == "ngCO2-C_mgC_h",
    Sample_cat %in% transcript_metadata$Sample_cat
  ) %>%
  inner_join(microresp_substrate_map, by = "Substrate") %>%
  left_join(
    transcript_metadata %>%
      select(Sample_cat, Plot_ID, Season),
    by = "Sample_cat"
  )

invalid_target_microresp <- microresp_target_data %>%
  filter(is.na(SIR) | !is.finite(SIR) | SIR < 0)
if (nrow(invalid_target_microresp) > 0) {
  write_csv(
    invalid_target_microresp,
    file.path(result_dir, "Figure3C_invalid_target_MicroResp_values.csv")
  )
  stop("Targeted MicroResp data contain invalid values.")
}

duplicate_target_microresp <- microresp_target_data %>%
  count(Sample_cat, Substrate, name = "record_count") %>%
  filter(record_count > 1)
if (nrow(duplicate_target_microresp) > 0) {
  stop("Targeted MicroResp sample-substrate records are duplicated.")
}

# Standardize before averaging so no high-magnitude substrate dominates.
microresp_family_response <- microresp_target_data %>%
  group_by(Substrate) %>%
  mutate(substrate_zscore = as.numeric(scale(SIR))) %>%
  ungroup() %>%
  group_by(target_family, Sample_cat, Sampling_Date, Plot_ID, Season) %>%
  summarise(
    standardized_value = mean(substrate_zscore),
    n_substrates = n_distinct(Substrate),
    .groups = "drop"
  )

missing_target_microresp <- transcript_metadata %>%
  select(Sample_cat, Sampling_Date, Plot_ID) %>%
  crossing(
    target_family = factor(
      target_family_levels,
      levels = target_family_levels
    )
  ) %>%
  anti_join(
    microresp_family_response,
    by = c("Sample_cat", "Sampling_Date", "Plot_ID", "target_family")
  )
write_csv(
  missing_target_microresp,
  file.path(result_dir, "Figure3C_missing_target_MicroResp_samples.csv")
)

target_timeseries <- bind_rows(
  microresp_family_response %>%
    transmute(
      target_family, Sample_cat, Sampling_Date, Plot_ID, Season,
      measure = "MicroResp response", standardized_value
    ),
  target_expression %>%
    transmute(
      target_family, Sample_cat, Sampling_Date, Plot_ID, Season,
      measure = "Transcript gene-set expression",
      standardized_value
    )
) %>%
  mutate(
    target_family = factor(target_family, levels = target_family_levels),
    measure = factor(
      measure,
      levels = c(
        "MicroResp response",
        "Transcript gene-set expression"
      )
    )
  )

target_timeseries_summary <- target_timeseries %>%
  group_by(target_family, measure, Sampling_Date, Season) %>%
  summarise(
    n = n(),
    mean_standardized_value = mean(standardized_value),
    se_standardized_value = sd(standardized_value) / sqrt(n),
    lower_se = mean_standardized_value - se_standardized_value,
    upper_se = mean_standardized_value + se_standardized_value,
    .groups = "drop"
  )
write_csv(
  target_timeseries,
  file.path(result_dir, "Figure3C_targeted_sample_level_standardized_values.csv")
)
write_csv(
  target_timeseries_summary,
  file.path(result_dir, "Figure3C_targeted_date_summary.csv")
)

# These correlations compare only the seven date-level means. They are
# descriptive checks of correspondence, not tests of mechanistic linkage.
target_date_pairs <- target_timeseries_summary %>%
  mutate(
    measure_key = recode(
      as.character(measure),
      "MicroResp response" = "microresp_z",
      "Transcript gene-set expression" = "transcript_z"
    )
  ) %>%
  select(
    target_family,
    Sampling_Date,
    measure_key,
    mean_standardized_value
  ) %>%
  pivot_wider(
    names_from = measure_key,
    values_from = mean_standardized_value
  )

target_correspondence <- target_date_pairs %>%
  group_by(target_family) %>%
  summarise(
    n_dates = sum(complete.cases(microresp_z, transcript_z)),
    pearson_r = cor(microresp_z, transcript_z, method = "pearson"),
    pearson_p = cor.test(
      microresp_z,
      transcript_z,
      method = "pearson"
    )$p.value,
    spearman_rho = cor(microresp_z, transcript_z, method = "spearman"),
    spearman_p = cor.test(
      microresp_z,
      transcript_z,
      method = "spearman",
      exact = FALSE
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    interpretation = paste(
      "Exploratory date-level correlation; descriptive only and not evidence",
      "of mechanistic linkage."
    )
  )

write_csv(
  target_correspondence,
  file.path(result_dir, "Figure3C_exploratory_date_level_correlations.csv")
)

correlation_caption <- target_correspondence %>%
  mutate(
    family_short = if_else(
      str_detect(as.character(target_family), "^Chitin"),
      "chitin family",
      "methyl oleate"
    ),
    correlation_statement = paste0(
      family_short,
      " Pearson r = ", formatC(pearson_r, digits = 2, format = "f"),
      " (p = ", format.pval(pearson_p, digits = 2),
      "); Spearman rho = ",
      formatC(spearman_rho, digits = 2, format = "f"),
      " (p = ", format.pval(spearman_p, digits = 2), ")"
    )
  ) %>%
  pull(correlation_statement) %>%
  paste(collapse = "; ")

figure_caption <- paste(
  "Panel B: the 40 functions were selected by abundance and prevalence,",
  "not temporal significance.",
  "Panel C: contrasting seasonal trajectories for pre-curated substrate",
  "families and transcript gene sets.",
  "Exploratory date-level correlations (n = 7 per family):",
  paste0(correlation_caption, "."),
  "These correlations are descriptive and do not establish mechanistic linkage."
) %>%
  str_wrap(width = 190)

write_lines(
  figure_caption,
  file.path(result_dir, "Figure3_caption.txt")
)

season_periods <- tribble(
  ~Season, ~xmin, ~xmax,
  "Plant-active, Apr-Oct", as.Date("2023-06-01"), as.Date("2023-10-31"),
  "Plant-dormant, Nov-Mar", as.Date("2023-11-01"), as.Date("2024-03-31"),
  "Plant-active, Apr-Oct", as.Date("2024-04-01"), as.Date("2024-06-30")
) %>%
  mutate(Season = factor(Season, levels = names(season_colors)))

panel_c <- ggplot(
  target_timeseries_summary,
  aes(Sampling_Date, mean_standardized_value)
) +
  geom_rect(
    data = season_periods,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Season),
    inherit.aes = FALSE,
    alpha = 0.16,
    color = NA
  ) +
  geom_ribbon(
    aes(x = Sampling_Date, ymin = lower_se, ymax = upper_se),
    inherit.aes = FALSE,
    fill = "grey45",
    alpha = 0.20
  ) +
  geom_hline(yintercept = 0, color = "grey60", linetype = 2) +
  geom_line(color = "black", linewidth = 0.8) +
  geom_point(
    aes(fill = Season),
    shape = 21,
    color = "black",
    size = 2.8
  ) +
  facet_grid(
    rows = vars(measure),
    cols = vars(target_family),
    scales = "free_y"
  ) +
  scale_fill_manual(values = season_colors, drop = FALSE) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b\n%Y") +
  labs(
    title = "C   Contrasting seasonal trajectories",
    subtitle = paste(
      "Biologically curated before temporal analysis; aligned date means",
      "and +/- SE are shown on separate response rows"
    ),
    x = NULL,
    y = "Standardized value (z-score)",
    fill = "Season"
  ) +
  theme_manuscript +
  theme(
    legend.position = "bottom",
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 9)
  ) +
  guides(fill = guide_legend(nrow = 1))

# ---- Save standalone panels and the combined figure -----------------------

ggsave(
  file.path(figure_dir, "Figure3A_transcript_BRITE_C_NMDS.png"),
  panel_a, width = 7.5, height = 6.2, dpi = 300
)
ggsave(
  file.path(figure_dir, "Figure3A_transcript_BRITE_C_NMDS.pdf"),
  panel_a, width = 7.5, height = 6.2
)
ggsave(
  file.path(figure_dir, "Figure3B_BRITE_C_temporal_heatmap.png"),
  panel_b, width = 9.3, height = 7.2, dpi = 300
)
ggsave(
  file.path(figure_dir, "Figure3B_BRITE_C_temporal_heatmap.pdf"),
  panel_b, width = 9.3, height = 7.2
)
ggsave(
  file.path(figure_dir, "Figure3C_targeted_MicroResp_transcript_series.png"),
  panel_c, width = 14.5, height = 6.4, dpi = 300
)
ggsave(
  file.path(figure_dir, "Figure3C_targeted_MicroResp_transcript_series.pdf"),
  panel_c, width = 14.5, height = 6.4
)

top_row <- panel_a + panel_b + plot_layout(widths = c(0.85, 1.35))
figure_3 <- top_row / panel_c +
  plot_layout(heights = c(1.0, 1.05)) +
  plot_annotation(
    caption = figure_caption,
    theme = theme(
      plot.caption = element_text(
        size = 8,
        hjust = 0,
        margin = margin(t = 8)
      )
    )
  )

ggsave(
  file.path(figure_dir, "Figure3_final.png"),
  figure_3, width = 16, height = 13, dpi = 300
)
ggsave(
  file.path(figure_dir, "Figure3_final.pdf"),
  figure_3, width = 16, height = 13
)

analysis_summary <- tibble(
  item = c(
    "Matched transcript samples",
    "Transcript plots",
    "Transcript sampling dates",
    "BRITE C functions in NMDS after exclusions",
    "BRITE C functions in compact heatmap",
    "NMDS stress",
    "PERMANOVA R2",
    "PERMANOVA F",
    "PERMANOVA p-value",
    "Dispersion p-value",
    "Targeted substrate families",
    "Curated KOs",
    "Curated KOs present in BRITE D",
    "Curated KOs detected in any transcript sample",
    "Chitin family exploratory date-level Pearson r",
    "Chitin family exploratory date-level Pearson p-value",
    "Chitin family exploratory date-level Spearman rho",
    "Chitin family exploratory date-level Spearman p-value",
    "Methyl oleate exploratory date-level Pearson r",
    "Methyl oleate exploratory date-level Pearson p-value",
    "Methyl oleate exploratory date-level Spearman rho",
    "Methyl oleate exploratory date-level Spearman p-value"
  ),
  value = c(
    nrow(transcript_metadata),
    n_distinct(transcript_metadata$Plot_ID),
    n_distinct(transcript_metadata$Sampling_Date),
    ncol(brite_c_matrix),
    length(major_brite_c),
    transcript_nmds$stress,
    permanova_r2,
    permanova_f,
    permanova_p,
    dispersion_p,
    n_distinct(target_gene_table$target_family),
    n_distinct(target_gene_table$target_ko),
    sum(target_ko_audit$present_in_brite_d),
    sum(target_ko_audit$detected_in_any_sample),
    target_correspondence$pearson_r[1],
    target_correspondence$pearson_p[1],
    target_correspondence$spearman_rho[1],
    target_correspondence$spearman_p[1],
    target_correspondence$pearson_r[2],
    target_correspondence$pearson_p[2],
    target_correspondence$spearman_rho[2],
    target_correspondence$spearman_p[2]
  )
)
write_csv(
  analysis_summary,
  file.path(result_dir, "Figure3_analysis_summary.csv")
)

cat("\nFigure 3 complete.\n")
print(analysis_summary, n = Inf)
