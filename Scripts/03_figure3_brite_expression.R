# =============================================================================
# Figure 3: field transcript seasonality and correspondence with strongest
# MicroResp temporal responses
# =============================================================================
#
# MAIN FIGURE
#
# A. Overall seasonal structure of field-expressed BRITE C functions
#
# B. Temporal expression of BRITE C pathways biologically linked a priori
#    to substrates tested with MicroResp
#
# C. Direct temporal comparison for the EXACT three substrates selected in
#    Figure 2 as having the strongest FDR-supported temporal MicroResp responses.
#
#    MicroResp and transcript expression are standardized SEPARATELY using
#    z-scores. This allows comparison of temporal direction and timing.
#    It does NOT make their absolute magnitudes quantitatively comparable.
#
# SUPPLEMENT
#
# Procrustes comparison of the complete MicroResp and BRITE C transcript
# multivariate structures.
#
# Scope: 100% irrigation treatment because metatranscriptomes are available
# only for that treatment.
# =============================================================================

rm(list = ls())


# =============================================================================
# Packages
# =============================================================================

required_packages <- c(
  "tidyverse",
  "lubridate",
  "vegan",
  "permute",
  "pheatmap",
  "patchwork"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  stop(
    "Install required packages: ",
    paste(missing_packages, collapse = ", ")
  )
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


# =============================================================================
# Plotting settings
# =============================================================================

season_colors <- c(
  "Plant-active, Apr-Oct" = "#2E8B57",
  "Plant-dormant, Nov-Mar" = "#5B6F8C"
)

figure_dir <- file.path(
  "Figures",
  "manuscript_final",
  "Figure3"
)

result_dir <- file.path(
  "Results",
  "manuscript_final",
  "Figure3"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  result_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

theme_manuscript <- theme_bw(
  base_size = 11
) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.key = element_blank(),
    strip.text = element_text(face = "bold")
  )


# =============================================================================
# Helper functions
# =============================================================================

clean_substrate_key <- function(x) {
  
  x %>%
    str_to_lower() %>%
    str_replace_all(
      "[^a-z0-9]",
      ""
    )
}


safe_zscore <- function(x) {
  
  x <- as.numeric(x)
  
  this_sd <- sd(
    x,
    na.rm = TRUE
  )
  
  if (
    !is.finite(this_sd) ||
    this_sd == 0
  ) {
    
    return(
      rep(
        NA_real_,
        length(x)
      )
    )
  }
  
  as.numeric(
    scale(x)
  )
}


# =============================================================================
# Field metadata for 100% irrigation treatment
# =============================================================================

field_data <- read_csv(
  file.path(
    "Data",
    "Data_merged_cleaned.csv"
  ),
  show_col_types = FALSE
) %>%
  
  filter(
    as.character(
      Irrigation_trmt_pct
    ) == "100"
  ) %>%
  
  transmute(
    Sample_cat,
    
    Sampling_Date =
      as.Date(
        Collection_date
      ),
    
    Plot_ID =
      factor(
        Plot
      ),
    
    Block =
      factor(
        Block
      ),
    
    Season =
      factor(
        if_else(
          month(
            Sampling_Date
          ) %in% 4:10,
          
          "Plant-active, Apr-Oct",
          
          "Plant-dormant, Nov-Mar"
        ),
        
        levels =
          names(
            season_colors
          )
      )
  ) %>%
  
  distinct()


# =============================================================================
# Load BRITE C and BRITE D transcript expression
# =============================================================================

brite_c_wide <- read_csv(
  file.path(
    "Data",
    "brite_C_expression.csv"
  ),
  show_col_types = FALSE
) %>%
  
  rename(
    brite_C = 1
  ) %>%
  
  filter(
    !is.na(
      brite_C
    ),
    brite_C != ""
  )


brite_d_wide <- read_csv(
  file.path(
    "Data",
    "brite_D_expression.csv"
  ),
  show_col_types = FALSE
) %>%
  
  rename(
    brite_D = 1
  ) %>%
  
  filter(
    !is.na(
      brite_D
    ),
    brite_D != ""
  )


transcript_columns <- names(
  brite_c_wide
)[-1]


transcript_metadata <- tibble(
  Sample_mt =
    transcript_columns,
  
  Sample_cat =
    str_remove(
      transcript_columns,
      "_MT$"
    )
) %>%
  
  left_join(
    field_data,
    by = "Sample_cat"
  ) %>%
  
  arrange(
    Sampling_Date,
    Plot_ID
  )


# =============================================================================
# Validate transcript metadata
# =============================================================================

unmatched_transcripts <- transcript_metadata %>%
  
  filter(
    is.na(
      Sampling_Date
    ) |
      is.na(
        Plot_ID
      ) |
      is.na(
        Season
      )
  )


write_csv(
  unmatched_transcripts,
  file.path(
    result_dir,
    "Figure3_unmatched_transcript_samples.csv"
  )
)


if (
  nrow(
    unmatched_transcripts
  ) > 0
) {
  
  stop(
    "Some transcript samples do not match the 100% field metadata."
  )
}


if (
  !identical(
    names(
      brite_d_wide
    )[-1],
    transcript_columns
  )
) {
  
  stop(
    "BRITE C and BRITE D matrices have different sample columns."
  )
}


if (
  anyDuplicated(
    transcript_metadata$Sample_cat
  ) > 0
) {
  
  stop(
    "Transcript sample identifiers are duplicated."
  )
}


# =============================================================================
# Remove curated non-soil/human pathways
# =============================================================================

human_exclusion <- read.delim(
  file.path(
    "Data",
    "pathwaystoremove",
    "kegg_human_exclude_2026-07-13.txt"
  ),
  check.names = FALSE
) %>%
  
  transmute(
    pathway_code =
      str_extract(
        as.character(
          path_ko
        ),
        "ko\\d{5}"
      )
  ) %>%
  
  filter(
    !is.na(
      pathway_code
    )
  ) %>%
  
  distinct()


# =============================================================================
# BRITE C long-format table
# =============================================================================

brite_c_long <- brite_c_wide %>%
  
  mutate(
    pathway_code =
      str_extract(
        brite_C,
        "ko\\d{5}"
      )
  ) %>%
  
  filter(
    !is.na(
      pathway_code
    ),
    
    !pathway_code %in%
      human_exclusion$pathway_code
  ) %>%
  
  select(
    -pathway_code
  ) %>%
  
  pivot_longer(
    cols =
      all_of(
        transcript_columns
      ),
    
    names_to =
      "Sample_mt",
    
    values_to =
      "normalized_expression"
  ) %>%
  
  left_join(
    transcript_metadata,
    by = "Sample_mt"
  ) %>%
  
  mutate(
    normalized_expression =
      as.numeric(
        normalized_expression
      )
  )


invalid_brite_c <- brite_c_long %>%
  
  filter(
    is.na(
      normalized_expression
    ) |
      !is.finite(
        normalized_expression
      ) |
      normalized_expression < 0
  )


if (
  nrow(
    invalid_brite_c
  ) > 0
) {
  
  write_csv(
    invalid_brite_c,
    file.path(
      result_dir,
      "Figure3_invalid_BRITE_C_values.csv"
    )
  )
  
  stop(
    "BRITE C expression contains invalid values."
  )
}


transcript_metadata %>%
  
  count(
    Sampling_Date,
    Season,
    name = "n_transcript_samples"
  ) %>%
  
  arrange(
    Sampling_Date
  ) %>%
  
  write_csv(
    file.path(
      result_dir,
      "Figure3_transcript_samples_by_date.csv"
    )
  )


# =============================================================================
# PANEL A
# Overall transcript functional seasonality
# =============================================================================

brite_c_sample_table <- brite_c_long %>%
  
  select(
    Sample_mt,
    brite_C,
    normalized_expression
  ) %>%
  
  pivot_wider(
    names_from =
      brite_C,
    
    values_from =
      normalized_expression,
    
    values_fill =
      0
  ) %>%
  
  arrange(
    match(
      Sample_mt,
      transcript_metadata$Sample_mt
    )
  )


brite_c_matrix <- brite_c_sample_table %>%
  
  column_to_rownames(
    "Sample_mt"
  ) %>%
  
  as.matrix()


brite_c_matrix <- brite_c_matrix[
  ,
  colSums(
    brite_c_matrix
  ) > 0,
  drop = FALSE
]


if (
  any(
    rowSums(
      brite_c_matrix
    ) == 0
  )
) {
  
  stop(
    "At least one transcript sample has zero BRITE C expression."
  )
}


transcript_metadata <- transcript_metadata %>%
  
  arrange(
    match(
      Sample_mt,
      rownames(
        brite_c_matrix
      )
    )
  )


# ---- NMDS -----------------------------------------------------------------

set.seed(42)

transcript_nmds <- metaMDS(
  brite_c_matrix,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  trace = FALSE
)


# ---- PERMANOVA ------------------------------------------------------------

permutation_design <- how(
  nperm = 999,
  blocks = transcript_metadata$Plot_ID
)


set.seed(42)

season_permanova <- adonis2(
  brite_c_matrix ~ Season,
  data = transcript_metadata,
  permutations = permutation_design,
  method = "bray"
)


# ---- Dispersion ------------------------------------------------------------

season_dispersion <- betadisper(
  vegdist(
    brite_c_matrix,
    method = "bray"
  ),
  transcript_metadata$Season
)


set.seed(42)

season_dispersion_test <- permutest(
  season_dispersion,
  permutations = permutation_design
)


permanova_table <- season_permanova %>%
  
  as.data.frame() %>%
  
  rownames_to_column(
    "term"
  )


dispersion_table <- season_dispersion_test$tab %>%
  
  as.data.frame() %>%
  
  rownames_to_column(
    "term"
  )


write_csv(
  permanova_table,
  file.path(
    result_dir,
    "Figure3A_BRITE_C_season_PERMANOVA.csv"
  )
)


write_csv(
  dispersion_table,
  file.path(
    result_dir,
    "Figure3A_BRITE_C_season_betadisper.csv"
  )
)


# ---- NMDS plotting scores -------------------------------------------------

nmds_scores <- scores(
  transcript_nmds,
  display = "sites"
) %>%
  
  as.data.frame() %>%
  
  rownames_to_column(
    "Sample_mt"
  ) %>%
  
  as_tibble() %>%
  
  left_join(
    transcript_metadata,
    by = "Sample_mt"
  )


permanova_r2 <- permanova_table %>%
  slice(1) %>%
  pull(R2)


permanova_p <- permanova_table %>%
  slice(1) %>%
  pull(`Pr(>F)`)


dispersion_p <- dispersion_table %>%
  slice(1) %>%
  pull(`Pr(>F)`)


panel_a_statistics <- paste0(
  "PERMANOVA: R\u00B2 = ",
  formatC(
    permanova_r2,
    digits = 3,
    format = "f"
  ),
  ", p = ",
  format.pval(
    permanova_p,
    digits = 2
  ),
  "\nDispersion: p = ",
  format.pval(
    dispersion_p,
    digits = 2
  ),
  "\nStress = ",
  formatC(
    transcript_nmds$stress,
    digits = 3,
    format = "f"
  )
)


panel_a <- ggplot(
  nmds_scores,
  aes(
    NMDS1,
    NMDS2,
    color = Season
  )
) +
  
  geom_point(
    aes(
      shape = Plot_ID
    ),
    size = 3,
    alpha = 0.85
  ) +
  
  annotate(
    "label",
    x = -Inf,
    y = Inf,
    label = panel_a_statistics,
    hjust = -0.05,
    vjust = 1.10,
    size = 3.4,
    lineheight = 1.05,
    linewidth = 0.25,
    fill = scales::alpha(
      "white",
      0.88
    )
  ) +
  
  scale_color_manual(
    values = season_colors,
    drop = FALSE
  ) +
  
  labs(
    title = "A",
    x = "NMDS1",
    y = "NMDS2",
    color = "Season",
    shape = "Plot"
  ) +
  
  theme_manuscript +
  
  theme(
    legend.position = "bottom"
  ) +
  
  guides(
    color = guide_legend(
      nrow = 2
    )
  )


# =============================================================================
# PANEL B
# MicroResp-linked BRITE C temporal heatmap
# =============================================================================


# =============================================================================
# Locate Microresp_pathways.csv
# =============================================================================

microresp_pathway_candidates <- c(
  
  "Microresp_pathways.csv",
  
  file.path(
    "Data",
    "Microresp_pathways.csv"
  ),
  
  file.path(
    "DataHub March",
    "Microresp_pathways.csv"
  ),
  
  paste0(
    "C:/Users/gara009/OneDrive - PNNL/",
    "Soil Microbiome SFA - A3.1 Rhizosphere and Moisture Effects/",
    "Prosser_Monthly/Manuscripts/A3-1-1_Prosser_time_biogeochem/",
    "Microresp_pathways.csv"
  )
)


existing_mapping_files <- microresp_pathway_candidates[
  file.exists(
    microresp_pathway_candidates
  )
]


if (
  length(
    existing_mapping_files
  ) == 0
) {
  
  stop(
    paste0(
      "Could not find Microresp_pathways.csv. Checked:\n",
      paste(
        microresp_pathway_candidates,
        collapse = "\n"
      )
    )
  )
}


microresp_pathway_file <- existing_mapping_files[1]


cat(
  "\nUsing MicroResp pathway mapping file:\n",
  microresp_pathway_file,
  "\n"
)


# =============================================================================
# Read MicroResp pathway mapping
# =============================================================================

microresp_pathway_mapping <- read_csv(
  microresp_pathway_file,
  show_col_types = FALSE
) %>%
  
  filter(
    !is.na(
      brite_C
    ),
    brite_C != ""
  ) %>%
  
  mutate(
    brite_C_ID =
      str_extract(
        brite_C,
        "^\\d{5}"
      ),
    
    brite_C_name =
      brite_C %>%
      str_remove(
        "^\\d{5}\\s+"
      ) %>%
      str_remove(
        "\\s+\\[PATH:ko\\d{5}\\]$"
      ),
    
    substrate_key =
      clean_substrate_key(
        Substrate
      )
  )


# =============================================================================
# Biological MicroResp group annotation only
# =============================================================================

mapping_with_groups <- microresp_pathway_mapping


pathway_group_counts <- mapping_with_groups %>%
  
  distinct(
    brite_C,
    Table_group,
    Substrate
  ) %>%
  
  count(
    brite_C,
    Table_group,
    name = "n_mapped_substrates"
  ) %>%
  
  group_by(
    brite_C
  ) %>%
  
  filter(
    n_mapped_substrates ==
      max(
        n_mapped_substrates
      )
  ) %>%
  
  summarise(
    Heatmap_group =
      if_else(
        n() == 1,
        first(
          Table_group
        ),
        "Mixed groups"
      ),
    .groups = "drop"
  )


# =============================================================================
# Transcript pathway coverage
# =============================================================================

pathway_expression_coverage <- brite_c_long %>%
  
  filter(
    brite_C %in%
      microresp_pathway_mapping$brite_C
  ) %>%
  
  group_by(
    brite_C
  ) %>%
  
  summarise(
    mean_expression =
      mean(
        normalized_expression
      ),
    
    detection_prevalence =
      mean(
        normalized_expression > 0
      ),
    
    .groups =
      "drop"
  )


# =============================================================================
# Date-level transcript expression and row z-score
# =============================================================================
#
# Each BRITE C pathway is standardized relative to its OWN temporal mean and SD.
#
# z = 0  -> approximately average expression for that pathway
# z > 0  -> above that pathway's annual average
# z < 0  -> below that pathway's annual average
#
# The heatmap therefore shows temporal pattern, not absolute transcript abundance.
# =============================================================================

pathway_date_expression <- brite_c_long %>%
  
  filter(
    brite_C %in%
      microresp_pathway_mapping$brite_C
  ) %>%
  
  group_by(
    brite_C,
    Sampling_Date,
    Season
  ) %>%
  
  summarise(
    mean_expression =
      mean(
        normalized_expression
      ),
    
    n_transcript_samples =
      n(),
    
    .groups =
      "drop"
  ) %>%
  
  group_by(
    brite_C
  ) %>%
  
  mutate(
    temporal_mean =
      mean(
        mean_expression
      ),
    
    temporal_sd =
      sd(
        mean_expression
      ),
    
    expression_zscore =
      if_else(
        is.finite(
          temporal_sd
        ) &
          temporal_sd > 0,
        
        (
          mean_expression -
            temporal_mean
        ) /
          temporal_sd,
        
        NA_real_
      )
  ) %>%
  
  ungroup()


# =============================================================================
# Determine which pathways can be shown
# =============================================================================

pathway_inclusion <- microresp_pathway_mapping %>%
  
  distinct(
    brite_C,
    brite_C_ID,
    brite_C_name
  ) %>%
  
  left_join(
    pathway_expression_coverage,
    by = "brite_C"
  ) %>%
  
  left_join(
    pathway_date_expression %>%
      distinct(
        brite_C,
        temporal_sd
      ),
    
    by = "brite_C"
  ) %>%
  
  mutate(
    represented_in_transcripts =
      !is.na(
        mean_expression
      ),
    
    sufficient_coverage =
      represented_in_transcripts &
      detection_prevalence > 0 &
      is.finite(
        temporal_sd
      ) &
      temporal_sd > 0,
    
    included_in_Figure3B =
      sufficient_coverage,
    
    exclusion_reason =
      case_when(
        !represented_in_transcripts ~
          "Not represented in transcript BRITE C data",
        
        detection_prevalence == 0 ~
          "No detected transcript expression",
        
        !is.finite(
          temporal_sd
        ) |
          temporal_sd == 0 ~
          "No temporal variance for row-wise standardization",
        
        TRUE ~
          "Included"
      )
  )


# =============================================================================
# Save exact mapping used in Panel B
# =============================================================================

figure3b_mapping_table <- mapping_with_groups %>%
  
  select(
    Substrate,
    Table_group,
    brite_C,
    brite_C_ID,
    brite_C_name
  ) %>%
  
  left_join(
    pathway_inclusion %>%
      select(
        brite_C,
        included_in_Figure3B,
        exclusion_reason,
        mean_expression,
        detection_prevalence
      ),
    
    by = "brite_C"
  ) %>%
  
  arrange(
    brite_C_ID,
    Table_group,
    Substrate
  )


write_csv(
  figure3b_mapping_table,
  file.path(
    result_dir,
    "Figure3B_MicroResp_to_BRITE_C_mapping_used.csv"
  )
)


included_pathways <- pathway_inclusion %>%
  
  filter(
    included_in_Figure3B
  ) %>%
  
  pull(
    brite_C
  )


if (
  length(
    included_pathways
  ) == 0
) {
  
  stop(
    "No mapped MicroResp-linked BRITE C pathways passed coverage checks."
  )
}


# =============================================================================
# Build heatmap matrix
# =============================================================================

brite_c_heatmap_table <- pathway_date_expression %>%
  
  filter(
    brite_C %in%
      included_pathways
  ) %>%
  
  left_join(
    microresp_pathway_mapping %>%
      distinct(
        brite_C,
        brite_C_name
      ),
    
    by = "brite_C"
  ) %>%
  
  mutate(
    date_label =
      factor(
        format(
          Sampling_Date,
          "%b %Y"
        ),
        
        levels =
          format(
            sort(
              unique(
                transcript_metadata$Sampling_Date
              )
            ),
            "%b %Y"
          )
      ),
    
    pathway_label =
      str_wrap(
        brite_C_name,
        width = 38
      )
  ) %>%
  
  select(
    pathway_label,
    date_label,
    expression_zscore
  ) %>%
  
  pivot_wider(
    names_from =
      date_label,
    
    values_from =
      expression_zscore
  )


brite_c_heatmap_matrix <- brite_c_heatmap_table %>%
  
  column_to_rownames(
    "pathway_label"
  ) %>%
  
  as.matrix()


# =============================================================================
# Cluster heatmap rows by transcript temporal pattern
# =============================================================================

row_cor <- cor(
  t(
    brite_c_heatmap_matrix
  ),
  use = "pairwise.complete.obs",
  method = "pearson"
)


row_cor[
  !is.finite(
    row_cor
  )
] <- 0


diag(
  row_cor
) <- 1


row_clustering <- hclust(
  as.dist(
    1 -
      row_cor
  ),
  method = "average"
)


# =============================================================================
# Heatmap annotations: biological MicroResp group only
# =============================================================================

pathway_annotations <- microresp_pathway_mapping %>%
  
  distinct(
    brite_C,
    brite_C_name
  ) %>%
  
  filter(
    brite_C %in%
      included_pathways
  ) %>%
  
  left_join(
    pathway_group_counts,
    by = "brite_C"
  ) %>%
  
  mutate(
    pathway_label =
      str_wrap(
        brite_C_name,
        width = 38
      )
  ) %>%
  
  select(
    pathway_label,
    
    MicroResp_group =
      Heatmap_group
  ) %>%
  
  column_to_rownames(
    "pathway_label"
  ) %>%
  
  as.data.frame()


heatmap_annotation <- transcript_metadata %>%
  
  distinct(
    Sampling_Date,
    Season
  ) %>%
  
  arrange(
    Sampling_Date
  ) %>%
  
  transmute(
    date_label =
      format(
        Sampling_Date,
        "%b %Y"
      ),
    
    Season
  ) %>%
  
  column_to_rownames(
    "date_label"
  )


microresp_group_levels <- sort(
  unique(
    pathway_annotations$MicroResp_group
  )
)


microresp_group_colors <- setNames(
  scales::hue_pal(
    l = 65,
    c = 90
  )(
    length(
      microresp_group_levels
    )
  ),
  microresp_group_levels
)


panel_b_heatmap <- pheatmap(
  brite_c_heatmap_matrix,
  
  cluster_rows =
    row_clustering,
  
  cluster_cols =
    FALSE,
  
  annotation_row =
    pathway_annotations,
  
  annotation_col =
    heatmap_annotation,
  
  annotation_names_row =
    FALSE,
  
  annotation_colors =
    list(
      Season =
        season_colors,
      
      MicroResp_group =
        microresp_group_colors
    ),
  
  color =
    colorRampPalette(
      c(
        "#2166AC",
        "white",
        "#B2182B"
      )
    )(
      101
    ),
  
  breaks =
    seq(
      -2.5,
      2.5,
      length.out = 102
    ),
  
  fontsize =
    8.5,
  
  fontsize_row =
    7.3,
  
  fontsize_col =
    8,
  
  border_color =
    NA,
  
  angle_col =
    45,
  
  main =
    paste(
      "B   MicroResp-linked field transcript pathways",
      "(transcript expression, row z-score)"
    ),
  
  legend_breaks =
    c(
      -2,
      0,
      2
    ),
  
  silent =
    TRUE
)


panel_b <- wrap_elements(
  full =
    panel_b_heatmap$gtable
)


# =============================================================================
# Load MicroResp data
# =============================================================================

microresp_data <- read_csv(
  file.path(
    "DataHub March",
    "Soil_Phenotypes",
    "SoilMicrobiomeSFA_prosser-time-series_MicroResp.csv"
  ),
  col_types = cols(
    .default =
      col_character()
  ),
  na = c(
    "",
    "NA"
  )
) %>%
  
  rename(
    Irrigation_trmt_cat =
      any_of(
        "Irrigaiton_trmt_cat"
      )
  ) %>%
  
  mutate(
    Sampling_Date =
      mdy(
        Collection_date
      ),
    
    Timepoint =
      suppressWarnings(
        as.numeric(
          Timepoint
        )
      ),
    
    SIR =
      suppressWarnings(
        as.numeric(
          SIR
        )
      ),
    
    substrate_key =
      clean_substrate_key(
        Substrate
      )
  ) %>%
  
  filter(
    Irrigation_trmt_pct ==
      "100",
    
    Timepoint ==
      5,
    
    Normalization ==
      "ngCO2-C_mgC_h",
    
    Substrate !=
      "Water"
  )


invalid_microresp <- microresp_data %>%
  
  filter(
    is.na(
      SIR
    ) |
      !is.finite(
        SIR
      ) |
      SIR < 0
  )


if (
  nrow(
    invalid_microresp
  ) > 0
) {
  
  write_csv(
    invalid_microresp,
    file.path(
      result_dir,
      "Figure3_invalid_MicroResp_values.csv"
    )
  )
  
  stop(
    "MicroResp data contain invalid values."
  )
}


# =============================================================================
# PANEL C
# Use exact top-three temporal substrates from Figure 2
# =============================================================================

figure2_ranking_file <- file.path(
  "Results",
  "manuscript_final",
  "Figure2",
  "Figure2_substrate_temporal_ranking.csv"
)


if (
  !file.exists(
    figure2_ranking_file
  )
) {
  
  stop(
    paste0(
      "Could not find the revised Figure 2 temporal ranking file:\n",
      figure2_ranking_file,
      "\nRun the revised Figure 2 script first."
    )
  )
}


figure2_temporal_ranking <- read_csv(
  figure2_ranking_file,
  show_col_types = FALSE
)


top3_table <- figure2_temporal_ranking %>%
  
  filter(
    displayed_in_main_figure
  ) %>%
  
  arrange(
    display_rank
  ) %>%
  
  slice_head(
    n = 3
  ) %>%
  
  transmute(
    Substrate,
    
    display_order =
      display_rank,
    
    temporal_rank,
    
    deviance_explained,
    
    p_adj_FDR,
    
    substrate_key =
      clean_substrate_key(
        Substrate
      ),
    
    display_name =
      str_replace_all(
        Substrate,
        "[._]",
        " "
      )
  )


if (
  nrow(
    top3_table
  ) != 3
) {
  
  stop(
    paste0(
      "Expected exactly three Figure 2 top temporal substrates, but found ",
      nrow(
        top3_table
      ),
      "."
    )
  )
}


cat(
  "\nTop three Figure 2 temporal MicroResp substrates used in Figure 3:\n"
)

print(
  top3_table %>%
    select(
      display_order,
      Substrate,
      deviance_explained,
      p_adj_FDR
    )
)


# =============================================================================
# Match top three substrates to BRITE C pathways
# =============================================================================

top3_pathway_mapping <- microresp_pathway_mapping %>%
  
  inner_join(
    top3_table %>%
      select(
        substrate_key,
        Figure2_Substrate =
          Substrate,
        display_name,
        display_order
      ),
    
    by =
      "substrate_key"
  ) %>%
  
  transmute(
    Substrate =
      Figure2_Substrate,
    
    display_name,
    
    display_order,
    
    brite_C,
    
    brite_C_ID,
    
    brite_C_name,
    
    Table_group
  ) %>%
  
  distinct()


missing_top3_mapping <- top3_table %>%
  
  anti_join(
    top3_pathway_mapping,
    by = "Substrate"
  )


write_csv(
  missing_top3_mapping,
  file.path(
    result_dir,
    "Figure3C_top3_missing_BRITE_C_mapping.csv"
  )
)


if (
  nrow(
    missing_top3_mapping
  ) > 0
) {
  
  stop(
    paste0(
      "At least one Figure 2 top-three substrate does not have a BRITE C ",
      "mapping in Microresp_pathways.csv. Review ",
      "Figure3C_top3_missing_BRITE_C_mapping.csv."
    )
  )
}


write_csv(
  top3_pathway_mapping,
  file.path(
    result_dir,
    "Figure3C_top3_MicroResp_BRITE_C_mapping.csv"
  )
)


# =============================================================================
# Determine matched transcript samples
# =============================================================================

matched_transcript_samples <- transcript_metadata %>%
  
  select(
    Sample_cat,
    Sampling_Date,
    Plot_ID,
    Season
  ) %>%
  
  distinct()


# =============================================================================
# MicroResp side of Panel C
# =============================================================================

microresp_top3 <- microresp_data %>%
  
  inner_join(
    top3_table %>%
      select(
        substrate_key,
        Figure2_Substrate =
          Substrate,
        display_name,
        display_order
      ),
    
    by =
      "substrate_key"
  ) %>%
  
  inner_join(
    matched_transcript_samples %>%
      select(
        Sample_cat,
        Transcript_Date =
          Sampling_Date,
        Transcript_Plot =
          Plot_ID,
        Transcript_Season =
          Season
      ),
    
    by =
      "Sample_cat"
  ) %>%
  
  filter(
    Sampling_Date ==
      Transcript_Date
  )


duplicate_top3_microresp <- microresp_top3 %>%
  
  count(
    Sample_cat,
    Figure2_Substrate,
    name = "record_count"
  ) %>%
  
  filter(
    record_count > 1
  )


if (
  nrow(
    duplicate_top3_microresp
  ) > 0
) {
  
  write_csv(
    duplicate_top3_microresp,
    file.path(
      result_dir,
      "Figure3C_duplicate_top3_MicroResp_records.csv"
    )
  )
  
  stop(
    "Duplicate matched MicroResp sample-substrate records detected."
  )
}


# =============================================================================
# Z-score MicroResp separately within substrate
# =============================================================================
#
# z = (SIR - substrate mean) / substrate SD
#
# This allows comparison of relative temporal change, not absolute magnitude.
# =============================================================================

microresp_top3_z <- microresp_top3 %>%
  
  group_by(
    Figure2_Substrate
  ) %>%
  
  mutate(
    standardized_value =
      safe_zscore(
        SIR
      )
  ) %>%
  
  ungroup()


if (
  any(
    !is.finite(
      microresp_top3_z$standardized_value
    )
  )
) {
  
  stop(
    "At least one top-three MicroResp substrate has zero or invalid variance."
  )
}


microresp_top3_date <- microresp_top3_z %>%
  
  group_by(
    Substrate =
      Figure2_Substrate,
    
    display_name,
    
    display_order,
    
    Sampling_Date,
    
    Season =
      Transcript_Season
  ) %>%
  
  summarise(
    n_samples =
      n(),
    
    mean_standardized_value =
      mean(
        standardized_value,
        na.rm = TRUE
      ),
    
    se_standardized_value =
      if_else(
        n() > 1,
        
        sd(
          standardized_value,
          na.rm = TRUE
        ) /
          sqrt(
            n()
          ),
        
        NA_real_
      ),
    
    .groups =
      "drop"
  ) %>%
  
  mutate(
    lower_se =
      mean_standardized_value -
      se_standardized_value,
    
    upper_se =
      mean_standardized_value +
      se_standardized_value,
    
    Measurement =
      "Potential substrate-response capacity (MicroResp)"
  )


# =============================================================================
# Transcript side of Panel C
# =============================================================================
#
# Each selected substrate is linked to one or more BRITE C pathways using the
# a priori Microresp_pathways.csv mapping.
#
# Each mapped transcript pathway is z-scored separately before pathways are
# averaged within a substrate/sample.
#
# This prevents high-abundance transcript pathways from dominating solely
# because their raw values are numerically larger.
# =============================================================================

transcript_top3_long <- brite_c_long %>%
  
  inner_join(
    top3_pathway_mapping %>%
      select(
        Substrate,
        display_name,
        display_order,
        brite_C,
        brite_C_name
      ),
    
    by =
      "brite_C"
  )


top3_transcript_pathway_audit <- top3_pathway_mapping %>%
  
  left_join(
    brite_c_long %>%
      
      group_by(
        brite_C
      ) %>%
      
      summarise(
        transcript_samples =
          n(),
        
        mean_expression =
          mean(
            normalized_expression
          ),
        
        expression_sd =
          sd(
            normalized_expression
          ),
        
        detected_any =
          any(
            normalized_expression > 0
          ),
        
        .groups =
          "drop"
      ),
    
    by =
      "brite_C"
  ) %>%
  
  mutate(
    usable_for_panel_c =
      !is.na(
        transcript_samples
      ) &
      is.finite(
        expression_sd
      ) &
      expression_sd > 0
  )


write_csv(
  top3_transcript_pathway_audit,
  file.path(
    result_dir,
    "Figure3C_top3_transcript_pathway_audit.csv"
  )
)


if (
  any(
    !top3_transcript_pathway_audit$usable_for_panel_c
  )
) {
  
  warning(
    paste0(
      "At least one mapped BRITE C pathway has zero variance or is absent. ",
      "It will be excluded from the standardized transcript trajectory. ",
      "See Figure3C_top3_transcript_pathway_audit.csv."
    )
  )
}


usable_top3_pathways <- top3_transcript_pathway_audit %>%
  
  filter(
    usable_for_panel_c
  ) %>%
  
  select(
    Substrate,
    brite_C
  )


transcript_top3_long <- transcript_top3_long %>%
  
  inner_join(
    usable_top3_pathways,
    by = c(
      "Substrate",
      "brite_C"
    )
  )


# =============================================================================
# Z-score each mapped transcript pathway independently
# =============================================================================

transcript_top3_z <- transcript_top3_long %>%
  
  group_by(
    Substrate,
    brite_C
  ) %>%
  
  mutate(
    pathway_zscore =
      safe_zscore(
        normalized_expression
      )
  ) %>%
  
  ungroup()


if (
  any(
    !is.finite(
      transcript_top3_z$pathway_zscore
    )
  )
) {
  
  stop(
    "Nonfinite transcript z-scores remain after transcript pathway filtering."
  )
}


# =============================================================================
# Average mapped transcript pathways within each sample
# =============================================================================

transcript_top3_sample <- transcript_top3_z %>%
  
  group_by(
    Substrate,
    display_name,
    display_order,
    Sample_cat,
    Sample_mt,
    Sampling_Date,
    Plot_ID,
    Season
  ) %>%
  
  summarise(
    standardized_value =
      mean(
        pathway_zscore,
        na.rm = TRUE
      ),
    
    n_mapped_pathways =
      n_distinct(
        brite_C
      ),
    
    .groups =
      "drop"
  )


transcript_top3_date <- transcript_top3_sample %>%
  
  group_by(
    Substrate,
    display_name,
    display_order,
    Sampling_Date,
    Season
  ) %>%
  
  summarise(
    n_samples =
      n(),
    
    mean_standardized_value =
      mean(
        standardized_value,
        na.rm = TRUE
      ),
    
    se_standardized_value =
      if_else(
        n() > 1,
        
        sd(
          standardized_value,
          na.rm = TRUE
        ) /
          sqrt(
            n()
          ),
        
        NA_real_
      ),
    
    .groups =
      "drop"
  ) %>%
  
  mutate(
    lower_se =
      mean_standardized_value -
      se_standardized_value,
    
    upper_se =
      mean_standardized_value +
      se_standardized_value,
    
    Measurement =
      "Related field transcript expression"
  )


# =============================================================================
# Combine Panel C measurements
# =============================================================================

target_timeseries_summary <- bind_rows(
  microresp_top3_date,
  transcript_top3_date
) %>%
  
  mutate(
    display_name =
      factor(
        display_name,
        
        levels =
          top3_table %>%
          arrange(
            display_order
          ) %>%
          pull(
            display_name
          )
      ),
    
    Measurement =
      factor(
        Measurement,
        
        levels =
          c(
            "Potential substrate-response capacity (MicroResp)",
            "Related field transcript expression"
          )
      )
  ) %>%
  
  arrange(
    display_order,
    Measurement,
    Sampling_Date
  )


write_csv(
  target_timeseries_summary,
  file.path(
    result_dir,
    "Figure3C_top3_targeted_date_summary.csv"
  )
)


# =============================================================================
# Save sample-level standardized data
# =============================================================================

microresp_sample_output <- microresp_top3_z %>%
  
  transmute(
    Substrate =
      Figure2_Substrate,
    
    display_name,
    
    display_order,
    
    Sample_cat,
    
    Sampling_Date,
    
    Plot_ID =
      Transcript_Plot,
    
    Season =
      Transcript_Season,
    
    Measurement =
      "Potential substrate-response capacity (MicroResp)",
    
    raw_value =
      SIR,
    
    standardized_value,
    
    n_mapped_pathways =
      NA_integer_
  )


transcript_sample_output <- transcript_top3_sample %>%
  
  transmute(
    Substrate,
    
    display_name,
    
    display_order,
    
    Sample_cat,
    
    Sampling_Date,
    
    Plot_ID,
    
    Season,
    
    Measurement =
      "Related field transcript expression",
    
    raw_value =
      NA_real_,
    
    standardized_value,
    
    n_mapped_pathways
  )


bind_rows(
  microresp_sample_output,
  transcript_sample_output
) %>%
  
  write_csv(
    file.path(
      result_dir,
      "Figure3C_top3_sample_level_standardized_values.csv"
    )
  )


# =============================================================================
# Seasonal shading for Panel C
# =============================================================================

season_periods <- tribble(
  ~Season, ~xmin, ~xmax,
  
  "Plant-active, Apr-Oct",
  as.Date(
    "2023-06-01"
  ),
  as.Date(
    "2023-10-31"
  ),
  
  "Plant-dormant, Nov-Mar",
  as.Date(
    "2023-11-01"
  ),
  as.Date(
    "2024-03-31"
  ),
  
  "Plant-active, Apr-Oct",
  as.Date(
    "2024-04-01"
  ),
  as.Date(
    "2024-06-30"
  )
) %>%
  
  mutate(
    Season =
      factor(
        Season,
        levels =
          names(
            season_colors
          )
      )
  )


# =============================================================================
# Panel C plot
# =============================================================================

panel_targeted <- ggplot(
  target_timeseries_summary,
  aes(
    Sampling_Date,
    mean_standardized_value
  )
) +
  
  geom_rect(
    data =
      season_periods,
    
    aes(
      xmin =
        xmin,
      
      xmax =
        xmax,
      
      ymin =
        -Inf,
      
      ymax =
        Inf,
      
      fill =
        Season
    ),
    
    inherit.aes =
      FALSE,
    
    alpha =
      0.14,
    
    color =
      NA
  ) +
  
  geom_ribbon(
    aes(
      ymin =
        lower_se,
      
      ymax =
        upper_se
    ),
    
    fill =
      "grey45",
    
    alpha =
      0.20
  ) +
  
  geom_hline(
    yintercept =
      0,
    
    color =
      "grey55",
    
    linetype =
      2,
    
    linewidth =
      0.5
  ) +
  
  geom_line(
    color =
      "black",
    
    linewidth =
      0.9
  ) +
  
  geom_point(
    aes(
      fill =
        Season
    ),
    
    shape =
      21,
    
    color =
      "black",
    
    size =
      2.8,
    
    stroke =
      0.5
  ) +
  
  facet_grid(
    rows =
      vars(
        Measurement
      ),
    
    cols =
      vars(
        display_name
      ),
    
    scales =
      "free_y"
  ) +
  
  scale_fill_manual(
    values =
      season_colors,
    
    drop =
      FALSE
  ) +
  
  scale_x_date(
    date_breaks =
      "2 months",
    
    date_labels =
      "%b\n%Y"
  ) +
  
  labs(
    title =
      "C   Strongest MicroResp temporal responses and related field transcript functions",
    
    subtitle =
      paste(
        "MicroResp and transcript expression were standardized separately;",
        "z-scores compare temporal direction and timing, not absolute magnitude"
      ),
    
    x =
      NULL,
    
    y =
      "Standardized value (z-score)",
    
    fill =
      "Season"
  ) +
  
  theme_manuscript +
  
  theme(
    legend.position =
      "bottom",
    
    strip.text.x =
      element_text(
        size = 10
      ),
    
    strip.text.y =
      element_text(
        size = 9
      )
  ) +
  
  guides(
    fill =
      guide_legend(
        nrow = 1
      )
  )


# =============================================================================
# Save Panel C
# =============================================================================

ggsave(
  file.path(
    figure_dir,
    "Figure3C_top3_MicroResp_transcript_trajectories.png"
  ),
  panel_targeted,
  width = 15,
  height = 6.7,
  dpi = 300
)


ggsave(
  file.path(
    figure_dir,
    "Figure3C_top3_MicroResp_transcript_trajectories.pdf"
  ),
  panel_targeted,
  width = 15,
  height = 6.7
)


# =============================================================================
# SUPPLEMENTARY ANALYSIS
# Global MicroResp-transcript Procrustes correspondence
# =============================================================================

procrustes_metadata <- transcript_metadata %>%
  
  filter(
    Sample_cat %in%
      microresp_data$Sample_cat
  ) %>%
  
  arrange(
    Plot_ID,
    Sampling_Date
  )


microresp_procrustes_data <- microresp_data %>%
  
  filter(
    Sample_cat %in%
      procrustes_metadata$Sample_cat
  )


duplicate_microresp <- microresp_procrustes_data %>%
  
  count(
    Sample_cat,
    Substrate,
    name =
      "record_count"
  ) %>%
  
  filter(
    record_count > 1
  )


if (
  nrow(
    duplicate_microresp
  ) > 0
) {
  
  stop(
    "Matched MicroResp sample-substrate records are duplicated."
  )
}


microresp_matrix <- microresp_procrustes_data %>%
  
  select(
    Sample_cat,
    Substrate,
    SIR
  ) %>%
  
  pivot_wider(
    names_from =
      Substrate,
    
    values_from =
      SIR
  ) %>%
  
  column_to_rownames(
    "Sample_cat"
  ) %>%
  
  as.matrix()


microresp_matrix <- microresp_matrix[
  procrustes_metadata$Sample_cat,
  ,
  drop = FALSE
]


transcript_procrustes_matrix <- brite_c_matrix[
  paste0(
    procrustes_metadata$Sample_cat,
    "_MT"
  ),
  ,
  drop = FALSE
]


rownames(
  transcript_procrustes_matrix
) <- procrustes_metadata$Sample_cat


if (
  any(
    is.na(
      microresp_matrix
    )
  ) |
  !identical(
    rownames(
      microresp_matrix
    ),
    rownames(
      transcript_procrustes_matrix
    )
  )
) {
  
  stop(
    "Procrustes matrices are incomplete or not aligned by sample."
  )
}


set.seed(43)

microresp_matched_nmds <- metaMDS(
  microresp_matrix,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  trace = FALSE
)


set.seed(44)

transcript_matched_nmds <- metaMDS(
  transcript_procrustes_matrix,
  distance = "bray",
  k = 2,
  trymax = 200,
  autotransform = FALSE,
  trace = FALSE
)


microresp_sites <- scores(
  microresp_matched_nmds,
  display = "sites"
)


transcript_sites <- scores(
  transcript_matched_nmds,
  display = "sites"
)


procrustes_fit <- procrustes(
  microresp_sites,
  transcript_sites,
  symmetric = TRUE
)


series_permutations <- how(
  nperm = 999,
  
  within =
    Within(
      type = "series",
      mirror = FALSE
    ),
  
  blocks =
    procrustes_metadata$Plot_ID
)


set.seed(45)

protest_result <- protest(
  microresp_sites,
  transcript_sites,
  permutations = series_permutations,
  symmetric = TRUE
)


permutation_description <- paste(
  "999 cyclic-shift permutations within Plot_ID; samples ordered",
  "chronologically within plot; mirror/reversal disabled. Exploratory because",
  "only seven nonstationary seasonal dates are available."
)


procrustes_scores <- tibble(
  Sample_cat =
    procrustes_metadata$Sample_cat,
  
  MicroResp1 =
    procrustes_fit$X[
      ,
      1
    ],
  
  MicroResp2 =
    procrustes_fit$X[
      ,
      2
    ],
  
  Transcript1 =
    procrustes_fit$Yrot[
      ,
      1
    ],
  
  Transcript2 =
    procrustes_fit$Yrot[
      ,
      2
    ]
) %>%
  
  left_join(
    procrustes_metadata %>%
      
      select(
        Sample_cat,
        Sampling_Date,
        Plot_ID,
        Season
      ),
    
    by =
      "Sample_cat"
  ) %>%
  
  mutate(
    residual_length =
      sqrt(
        (
          MicroResp1 -
            Transcript1
        )^2 +
          (
            MicroResp2 -
              Transcript2
          )^2
      )
  )


write_csv(
  procrustes_scores,
  file.path(
    result_dir,
    "FigureS3_sample_level_Procrustes_scores.csv"
  )
)


procrustes_results <- tibble(
  statistic = c(
    "Procrustes correlation",
    "PROTEST permutation p-value",
    "Procrustes residual sum of squares",
    "Matched samples",
    "Matched plots",
    "Matched sampling dates",
    "MicroResp matched NMDS stress",
    "Transcript matched NMDS stress"
  ),
  
  value = c(
    protest_result$t0,
    protest_result$signif,
    procrustes_fit$ss,
    nrow(
      procrustes_metadata
    ),
    n_distinct(
      procrustes_metadata$Plot_ID
    ),
    n_distinct(
      procrustes_metadata$Sampling_Date
    ),
    microresp_matched_nmds$stress,
    transcript_matched_nmds$stress
  ),
  
  permutation_design =
    permutation_description
)


write_csv(
  procrustes_results,
  file.path(
    result_dir,
    "FigureS3_PROTEST_results.csv"
  )
)


procrustes_points <- bind_rows(
  
  procrustes_scores %>%
    
    transmute(
      Sample_cat,
      Sampling_Date,
      Plot_ID,
      Season,
      
      Data_type =
        "Potential substrate-response capacity (MicroResp)",
      
      Axis1 =
        MicroResp1,
      
      Axis2 =
        MicroResp2
    ),
  
  procrustes_scores %>%
    
    transmute(
      Sample_cat,
      Sampling_Date,
      Plot_ID,
      Season,
      
      Data_type =
        "Field transcript expression",
      
      Axis1 =
        Transcript1,
      
      Axis2 =
        Transcript2
    )
) %>%
  
  mutate(
    Data_type =
      factor(
        Data_type,
        
        levels = c(
          "Potential substrate-response capacity (MicroResp)",
          "Field transcript expression"
        )
      )
  )


panel_procrustes_statistics <- paste0(
  "Procrustes r = ",
  formatC(
    protest_result$t0,
    digits = 3,
    format = "f"
  ),
  "\nExploratory PROTEST p = ",
  format.pval(
    protest_result$signif,
    digits = 2
  ),
  "\nn = ",
  nrow(
    procrustes_metadata
  ),
  " matched samples"
)


panel_procrustes <- ggplot() +
  
  geom_segment(
    data =
      procrustes_scores,
    
    aes(
      x =
        MicroResp1,
      
      y =
        MicroResp2,
      
      xend =
        Transcript1,
      
      yend =
        Transcript2,
      
      color =
        Season
    ),
    
    linewidth =
      0.55,
    
    alpha =
      0.60
  ) +
  
  geom_point(
    data =
      procrustes_points,
    
    aes(
      x =
        Axis1,
      
      y =
        Axis2,
      
      fill =
        Season,
      
      shape =
        Data_type
    ),
    
    color =
      "black",
    
    size =
      2.7,
    
    stroke =
      0.45
  ) +
  
  annotate(
    "label",
    
    x =
      -Inf,
    
    y =
      Inf,
    
    label =
      panel_procrustes_statistics,
    
    hjust =
      -0.05,
    
    vjust =
      1.10,
    
    size =
      3.3,
    
    lineheight =
      1.05,
    
    linewidth =
      0.25,
    
    fill =
      scales::alpha(
        "white",
        0.88
      )
  ) +
  
  scale_color_manual(
    values =
      season_colors,
    
    drop =
      FALSE
  ) +
  
  scale_fill_manual(
    values =
      season_colors,
    
    drop =
      FALSE
  ) +
  
  scale_shape_manual(
    values = c(
      "Potential substrate-response capacity (MicroResp)" =
        21,
      
      "Field transcript expression" =
        24
    )
  ) +
  
  coord_equal() +
  
  labs(
    title =
      "Global MicroResp-transcript correspondence",
    
    subtitle =
      "Shorter paired segments indicate stronger multivariate agreement",
    
    x =
      "Procrustes axis 1",
    
    y =
      "Procrustes axis 2",
    
    color =
      "Season",
    
    fill =
      "Season",
    
    shape =
      "Measurement"
  ) +
  
  theme_manuscript +
  
  theme(
    legend.position =
      "bottom"
  ) +
  
  guides(
    color =
      "none",
    
    fill =
      guide_legend(
        nrow = 2,
        override.aes =
          list(
            shape = 21
          )
      ),
    
    shape =
      guide_legend(
        nrow = 2
      )
  )


ggsave(
  file.path(
    figure_dir,
    "FigureS3_global_MicroResp_transcript_Procrustes.png"
  ),
  panel_procrustes,
  width = 8.2,
  height = 6.6,
  dpi = 300
)


ggsave(
  file.path(
    figure_dir,
    "FigureS3_global_MicroResp_transcript_Procrustes.pdf"
  ),
  panel_procrustes,
  width = 8.2,
  height = 6.6
)


# =============================================================================
# Save main Panels A and B
# =============================================================================

ggsave(
  file.path(
    figure_dir,
    "Figure3A_transcript_BRITE_C_NMDS.png"
  ),
  panel_a,
  width = 7.5,
  height = 6.2,
  dpi = 300
)


ggsave(
  file.path(
    figure_dir,
    "Figure3A_transcript_BRITE_C_NMDS.pdf"
  ),
  panel_a,
  width = 7.5,
  height = 6.2
)


ggsave(
  file.path(
    figure_dir,
    "Figure3B_MicroResp_linked_BRITE_C_heatmap.png"
  ),
  panel_b,
  width = 10.5,
  height = 7.4,
  dpi = 300
)


ggsave(
  file.path(
    figure_dir,
    "Figure3B_MicroResp_linked_BRITE_C_heatmap.pdf"
  ),
  panel_b,
  width = 10.5,
  height = 7.4
)


# =============================================================================
# FINAL MAIN FIGURE 3
# =============================================================================

top_row <- panel_a +
  panel_b +
  plot_layout(
    widths = c(
      0.82,
      1.38
    )
  )


figure_3 <- top_row /
  panel_targeted +
  
  plot_layout(
    heights = c(
      1.05,
      1.00
    )
  )


ggsave(
  file.path(
    figure_dir,
    "Figure3_final.png"
  ),
  figure_3,
  width = 16,
  height = 12.5,
  dpi = 300
)


ggsave(
  file.path(
    figure_dir,
    "Figure3_final.pdf"
  ),
  figure_3,
  width = 16,
  height = 12.5
)


# =============================================================================
# Main Figure 3 caption draft
# =============================================================================

top3_caption_names <- top3_table %>%
  
  arrange(
    display_order
  ) %>%
  
  pull(
    display_name
  ) %>%
  
  paste(
    collapse = ", "
  )


main_figure_caption <- paste0(
  
  "Figure 3. Seasonal field gene expression is structured but does not ",
  "necessarily track laboratory substrate-response capacity. ",
  
  "(A) Bray-Curtis NMDS of BRITE C transcript composition for the 100% ",
  "irrigation treatment. Field-expressed functional composition differed ",
  "between plant-active and plant-dormant samples ",
  
  "(PERMANOVA R2 = ",
  formatC(
    permanova_r2,
    digits = 3,
    format = "f"
  ),
  
  ", p = ",
  format.pval(
    permanova_p,
    digits = 2
  ),
  
  "; stress = ",
  formatC(
    transcript_nmds$stress,
    digits = 3,
    format = "f"
  ),
  
  "). Multivariate dispersion did not differ between seasons ",
  
  "(p = ",
  format.pval(
    dispersion_p,
    digits = 2
  ),
  
  "). ",
  
  "(B) Temporal expression of BRITE C pathways biologically linked a priori ",
  "to substrates included in the MicroResp assay. Rows show pathway-specific ",
  "z-scores calculated from date-mean transcript expression, so colors ",
  "represent expression relative to each pathway's own temporal mean rather ",
  "than absolute expression. Row annotations indicate the associated biological ",
  "MicroResp substrate group. ",
  
  "(C) Temporal comparison of the three strongest FDR-supported MicroResp ",
  "responses identified independently in Figure 2 (",
  top3_caption_names,
  ") with expression of their mapped BRITE C functions. MicroResp responses ",
  "were standardized separately within substrate and transcript expression ",
  "was standardized separately within mapped pathway before date-level ",
  "summarization. Z-scores therefore show relative temporal changes within ",
  "each measurement and should not be interpreted as equivalent absolute ",
  "magnitudes between MicroResp and transcript expression."
)


write_lines(
  str_wrap(
    main_figure_caption,
    width = 190
  ),
  file.path(
    result_dir,
    "Figure3_caption.txt"
  )
)


# =============================================================================
# Supplementary Procrustes caption
# =============================================================================

supplementary_procrustes_caption <- paste0(
  
  "Figure S3. Global correspondence between MicroResp substrate-response ",
  "profiles and field BRITE C transcript composition. Procrustes analysis ",
  "compared Bray-Curtis NMDS configurations for matched 100% irrigation ",
  "samples. Each segment connects the MicroResp and transcript positions for ",
  "the same sample after optimal rotation and scaling of the two ordinations; ",
  "shorter segments indicate greater agreement between the two multivariate ",
  "configurations. Procrustes correlation = ",
  formatC(
    protest_result$t0,
    digits = 3,
    format = "f"
  ),
  
  "; exploratory PROTEST p = ",
  format.pval(
    protest_result$signif,
    digits = 2
  ),
  
  "; n = ",
  nrow(
    procrustes_metadata
  ),
  
  " matched samples. Permutations were performed as cyclic shifts within ",
  "plots while preserving chronological sample order. Because only seven ",
  "nonstationary seasonal dates were available, the permutation p-value is ",
  "treated as exploratory."
)


write_lines(
  str_wrap(
    supplementary_procrustes_caption,
    width = 190
  ),
  file.path(
    result_dir,
    "FigureS3_Procrustes_caption.txt"
  )
)


# =============================================================================
# Analysis summary
# =============================================================================

analysis_summary <- tibble(
  item = c(
    "Matched transcript samples",
    "Transcript plots",
    "Transcript sampling dates",
    "BRITE C functions in Panel A NMDS",
    "Unique MicroResp-linked BRITE C pathways",
    "MicroResp-linked pathways included in Panel B",
    "Transcript NMDS stress",
    "PERMANOVA R2",
    "PERMANOVA p-value",
    "Dispersion p-value",
    "Figure 2 top temporal substrate 1",
    "Figure 2 top temporal substrate 2",
    "Figure 2 top temporal substrate 3",
    "Top-three substrate-pathway mappings",
    "Procrustes correlation",
    "PROTEST exploratory p-value",
    "Procrustes residual sum of squares",
    "Matched samples in Procrustes",
    "MicroResp matched NMDS stress",
    "Transcript matched NMDS stress"
  ),
  
  value = c(
    nrow(
      transcript_metadata
    ),
    
    n_distinct(
      transcript_metadata$Plot_ID
    ),
    
    n_distinct(
      transcript_metadata$Sampling_Date
    ),
    
    ncol(
      brite_c_matrix
    ),
    
    n_distinct(
      microresp_pathway_mapping$brite_C
    ),
    
    length(
      included_pathways
    ),
    
    transcript_nmds$stress,
    
    permanova_r2,
    
    permanova_p,
    
    dispersion_p,
    
    top3_table$Substrate[1],
    
    top3_table$Substrate[2],
    
    top3_table$Substrate[3],
    
    nrow(
      top3_pathway_mapping
    ),
    
    protest_result$t0,
    
    protest_result$signif,
    
    procrustes_fit$ss,
    
    nrow(
      procrustes_metadata
    ),
    
    microresp_matched_nmds$stress,
    
    transcript_matched_nmds$stress
  )
)


write_csv(
  analysis_summary,
  file.path(
    result_dir,
    "Figure3_analysis_summary.csv"
  )
)


# =============================================================================
# Console output
# =============================================================================

cat(
  "\nRevised Figure 3 complete.\n"
)


cat(
  "\nTop-three Figure 2 substrates used in Panel C:\n"
)


print(
  top3_table %>%
    
    select(
      display_order,
      Substrate,
      deviance_explained,
      p_adj_FDR
    )
)


cat(
  "\nBRITE C mappings used for those substrates:\n"
)


print(
  top3_pathway_mapping %>%
    
    arrange(
      display_order,
      brite_C_ID
    ) %>%
    
    select(
      Substrate,
      brite_C_ID,
      brite_C_name,
      Table_group
    )
)


cat(
  "\nAnalysis summary:\n"
)


print(
  analysis_summary,
  n = Inf
)