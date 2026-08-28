# =============================================================================
# Figure 2: seasonal MicroResp substrate-response capacity
# =============================================================================
#
# MicroResp was measured under standardized laboratory conditions at room
# temperature and 60% water-holding capacity. Results are interpreted as
# potential microbial substrate-induced respiration, not in situ activity.
#
# Scope: 100% irrigation, 5-hour carbon-normalized SIR responses.
#
# MAIN FIGURE:
#   A. NMDS of complete MicroResp substrate-response profile
#   B-D. Top 3 substrates with strongest temporal structure, ranked by
#        temporal GAM deviance explained among FDR-significant substrates
#
# SUPPLEMENT:
#   All substrate-specific temporal GAM trajectories
#   Fitted-trajectory clustering and heatmap
# =============================================================================

rm(list = ls())

required_packages <- c(
  "tidyverse",
  "lubridate",
  "vegan",
  "permute",
  "mgcv",
  "ggrepel",
  "cluster",
  "pheatmap",
  "RColorBrewer",
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
  library(mgcv)
  library(ggrepel)
  library(cluster)
  library(pheatmap)
  library(RColorBrewer)
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
  "Figure2"
)

result_dir <- file.path(
  "Results",
  "manuscript_final",
  "Figure2"
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
  base_size = 14
) +
  theme(
    panel.grid.minor = element_blank(),
    
    plot.title = element_text(
      face = "bold",
      size = 18
    ),
    
    plot.subtitle = element_text(
      size = 13
    ),
    
    axis.title = element_text(
      size = 15
    ),
    
    axis.text = element_text(
      size = 12
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 14
    ),
    
    legend.title = element_text(
      size = 13
    ),
    
    legend.text = element_text(
      size = 12
    ),
    
    legend.key = element_blank()
  )


# =============================================================================
# Load MicroResp data
# =============================================================================

microresp_file <- file.path(
  "DataHub March",
  "Soil_Phenotypes",
  "SoilMicrobiomeSFA_prosser-time-series_MicroResp.csv"
)

microresp_raw <- read_csv(
  microresp_file,
  col_types = cols(
    .default = col_character()
  ),
  na = c(
    "",
    "NA"
  )
)

microresp_data <- microresp_raw %>%
  
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
    
    Plot_ID =
      factor(
        Plot
      ),
    
    Block =
      factor(
        Block
      )
  ) %>%
  
  filter(
    Irrigation_trmt_pct == "100",
    Timepoint == 5,
    Normalization == "ngCO2-C_mgC_h",
    Substrate != "Water"
  ) %>%
  
  mutate(
    Season = factor(
      if_else(
        month(
          Sampling_Date
        ) %in% 4:10,
        "Plant-active, Apr-Oct",
        "Plant-dormant, Nov-Mar"
      ),
      levels = names(
        season_colors
      )
    ),
    
    time_days =
      as.numeric(
        Sampling_Date -
          min(
            Sampling_Date
          )
      ),
    
    sequential_month =
      dense_rank(
        Sampling_Date
      )
  ) %>%
  
  arrange(
    Sampling_Date,
    Plot_ID,
    Substrate
  )


# =============================================================================
# Basic data validation
# =============================================================================

invalid_sir <- microresp_data %>%
  filter(
    is.na(
      SIR
    ) |
      !is.finite(
        SIR
      ) |
      SIR < 0
  )

if (nrow(invalid_sir) > 0) {
  
  write_csv(
    invalid_sir,
    file.path(
      result_dir,
      "MicroResp_invalid_SIR_values.csv"
    )
  )
  
  stop(
    paste0(
      "The Bray-Curtis input contains missing, nonfinite, or negative SIR ",
      "values. Review the saved invalid-value table rather than imputing ",
      "silently."
    )
  )
}


duplicate_sample_substrate <- microresp_data %>%
  count(
    Sample_cat,
    Substrate,
    name = "record_count"
  ) %>%
  filter(
    record_count > 1
  )

if (nrow(duplicate_sample_substrate) > 0) {
  
  write_csv(
    duplicate_sample_substrate,
    file.path(
      result_dir,
      "MicroResp_duplicate_sample_substrate.csv"
    )
  )
  
  stop(
    "Duplicate Sample_cat by Substrate records require review."
  )
}


# =============================================================================
# Sample and substrate metadata
# =============================================================================

sample_metadata <- microresp_data %>%
  
  distinct(
    Sample_cat,
    Sampling_Date,
    Plot_ID,
    Block,
    Season,
    sequential_month,
    time_days
  ) %>%
  
  arrange(
    Sampling_Date,
    Plot_ID
  )


substrate_metadata <- microresp_data %>%
  
  distinct(
    Substrate,
    Substrate_category
  ) %>%
  
  arrange(
    Substrate
  )


# =============================================================================
# Bray-Curtis community matrix
# =============================================================================

community_wide <- microresp_data %>%
  
  select(
    Sample_cat,
    Substrate,
    SIR
  ) %>%
  
  pivot_wider(
    names_from = Substrate,
    values_from = SIR
  ) %>%
  
  right_join(
    sample_metadata,
    by = "Sample_cat"
  ) %>%
  
  arrange(
    Sampling_Date,
    Plot_ID
  )


substrate_columns <- substrate_metadata$Substrate


community_matrix <- community_wide %>%
  
  select(
    all_of(
      substrate_columns
    )
  ) %>%
  
  as.matrix()


rownames(
  community_matrix
) <- community_wide$Sample_cat

mode(
  community_matrix
) <- "numeric"


if (
  anyNA(
    community_matrix
  ) ||
  any(
    !is.finite(
      community_matrix
    )
  ) ||
  any(
    community_matrix < 0
  )
) {
  
  stop(
    "The final Bray-Curtis community matrix is not complete and nonnegative."
  )
}


write_csv(
  community_wide,
  file.path(
    result_dir,
    "Figure2_MicroResp_community_matrix_100pct.csv"
  )
)


# =============================================================================
# NMDS and repeated-plot seasonal tests
# =============================================================================

set.seed(42)

nmds_model <- metaMDS(
  community_matrix,
  distance = "bray",
  k = 2,
  trymax = 500,
  autotransform = FALSE,
  trace = FALSE
)


nmds_scores <- scores(
  nmds_model,
  display = "sites"
) %>%
  
  as.data.frame() %>%
  
  rownames_to_column(
    "Sample_cat"
  ) %>%
  
  as_tibble() %>%
  
  left_join(
    sample_metadata,
    by = "Sample_cat"
  )


# Restrict permutations within repeatedly sampled plots
permutation_design <- how(
  nperm = 999,
  blocks = sample_metadata$Plot_ID
)

set.seed(42)

season_permanova <- adonis2(
  community_matrix ~ Season,
  data = sample_metadata,
  method = "bray",
  permutations = permutation_design
) %>%
  
  as.data.frame() %>%
  
  rownames_to_column(
    "term"
  ) %>%
  
  as_tibble()


season_dispersion <- betadisper(
  vegdist(
    community_matrix,
    method = "bray"
  ),
  sample_metadata$Season
)

set.seed(42)

season_dispersion_test <- permutest(
  season_dispersion,
  permutations = permutation_design
)$tab %>%
  
  as.data.frame() %>%
  
  rownames_to_column(
    "term"
  ) %>%
  
  as_tibble()


write_csv(
  season_permanova,
  file.path(
    result_dir,
    "Figure2_NMDS_season_PERMANOVA.csv"
  )
)

write_csv(
  season_dispersion_test,
  file.path(
    result_dir,
    "Figure2_NMDS_season_betadisper.csv"
  )
)


# =============================================================================
# Extract NMDS statistics for plotting
# =============================================================================

permanova_row <- season_permanova %>%
  
  filter(
    term == "Season"
  ) %>%
  
  slice(
    1
  )


if (nrow(permanova_row) == 0) {
  
  permanova_row <- season_permanova %>%
    
    filter(
      !term %in% c(
        "Residual",
        "Total"
      )
    ) %>%
    
    slice(
      1
    )
}


dispersion_row <- season_dispersion_test %>%
  
  filter(
    !str_detect(
      str_to_lower(
        term
      ),
      "residual"
    )
  ) %>%
  
  slice(
    1
  )


permanova_R2 <-
  permanova_row[["R2"]][1]

permanova_F <-
  permanova_row[["F"]][1]

permanova_p <-
  permanova_row[["Pr(>F)"]][1]

dispersion_F <-
  dispersion_row[["F"]][1]

dispersion_p <-
  dispersion_row[["Pr(>F)"]][1]


format_p <- function(p) {
  
  if (
    length(p) == 0 ||
    is.na(p)
  ) {
    return(
      "NA"
    )
  }
  
  if (p < 0.001) {
    return(
      "< 0.001"
    )
  }
  
  paste0(
    "= ",
    sprintf(
      "%.3f",
      p
    )
  )
}


nmds_stats_label <- paste0(
  "PERMANOVA: R\u00B2 = ",
  sprintf(
    "%.3f",
    permanova_R2
  ),
  ", F = ",
  sprintf(
    "%.2f",
    permanova_F
  ),
  ", p ",
  format_p(
    permanova_p
  ),
  "\n",
  "Dispersion: F = ",
  sprintf(
    "%.2f",
    dispersion_F
  ),
  ", p ",
  format_p(
    dispersion_p
  )
)


# =============================================================================
# Environmental variables and envfit
# =============================================================================

field_environment <- read_csv(
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
    
    total_soil_respiration =
      Resp_rate_mg_C_m2_h,
    
    root_biomass =
      Roots_gdw,
    
    soil_temperature =
      SoilTemp_mean,
    
    gravimetric_water_content =
      GWC,
    
    microbial_biomass_C =
      MBC_ug_gdw_bgc
  )


envfit_candidates <- nmds_scores %>%
  
  left_join(
    field_environment,
    by = "Sample_cat"
  ) %>%
  
  select(
    Sample_cat,
    Plot_ID,
    NMDS1,
    NMDS2,
    soil_temperature,
    gravimetric_water_content,
    root_biomass,
    total_soil_respiration,
    microbial_biomass_C
  )


candidate_names <- c(
  "soil_temperature",
  "gravimetric_water_content",
  "root_biomass",
  "total_soil_respiration",
  "microbial_biomass_C"
)


envfit_complete <- envfit_candidates %>%
  
  drop_na(
    all_of(
      candidate_names
    )
  )


candidate_correlations <- cor(
  envfit_complete %>%
    
    select(
      all_of(
        candidate_names
      )
    ),
  
  use = "complete.obs"
)


selected_envfit_variables <- character()

envfit_screening <- tibble()


for (candidate_name in candidate_names) {
  
  if (
    length(
      selected_envfit_variables
    ) == 0
  ) {
    
    maximum_existing_correlation <- 0
    
  } else {
    
    maximum_existing_correlation <- max(
      abs(
        candidate_correlations[
          candidate_name,
          selected_envfit_variables,
          drop = TRUE
        ]
      )
    )
  }
  
  
  keep_candidate <-
    maximum_existing_correlation < 0.75
  
  
  envfit_screening <- bind_rows(
    envfit_screening,
    
    tibble(
      variable =
        candidate_name,
      
      maximum_absolute_correlation_with_retained =
        maximum_existing_correlation,
      
      retained =
        keep_candidate,
      
      decision =
        if_else(
          keep_candidate,
          "Retained for final envfit",
          "Excluded as redundant at absolute r >= 0.75"
        )
    )
  )
  
  
  if (keep_candidate) {
    
    selected_envfit_variables <- c(
      selected_envfit_variables,
      candidate_name
    )
  }
}


write_csv(
  envfit_screening,
  file.path(
    result_dir,
    "Figure2_envfit_variable_screening.csv"
  )
)


write_csv(
  as.data.frame(
    candidate_correlations
  ) %>%
    
    rownames_to_column(
      "variable"
    ),
  
  file.path(
    result_dir,
    "Figure2_envfit_candidate_correlations.csv"
  )
)


envfit_permutations <- how(
  nperm = 999,
  blocks = envfit_complete$Plot_ID
)


set.seed(42)

environment_fit <- envfit(
  envfit_complete %>%
    
    select(
      NMDS1,
      NMDS2
    ),
  
  envfit_complete %>%
    
    select(
      all_of(
        selected_envfit_variables
      )
    ),
  
  permutations =
    envfit_permutations
)


environment_vectors <- scores(
  environment_fit,
  display = "vectors"
) %>%
  
  as.data.frame() %>%
  
  rownames_to_column(
    "variable"
  ) %>%
  
  as_tibble() %>%
  
  mutate(
    R2 =
      environment_fit$vectors$r[
        variable
      ],
    
    p_value =
      environment_fit$vectors$pvals[
        variable
      ],
    
    p_adj_FDR =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  
  arrange(
    p_adj_FDR,
    desc(
      R2
    )
  )


write_csv(
  environment_vectors,
  file.path(
    result_dir,
    "Figure2_NMDS_envfit_full_results.csv"
  )
)


environment_plot_vectors <- environment_vectors %>%
  
  filter(
    p_adj_FDR < 0.05
  ) %>%
  
  slice_max(
    R2,
    n = 5,
    with_ties = FALSE
  )


if (
  nrow(
    environment_plot_vectors
  ) > 0
) {
  
  vector_multiplier <- 0.75 * min(
    max(
      abs(
        nmds_scores$NMDS1
      )
    ) /
      max(
        abs(
          environment_plot_vectors$NMDS1
        )
      ),
    
    max(
      abs(
        nmds_scores$NMDS2
      )
    ) /
      max(
        abs(
          environment_plot_vectors$NMDS2
        )
      )
  )
  
  
  environment_plot_vectors <- environment_plot_vectors %>%
    
    mutate(
      xend =
        NMDS1 *
        vector_multiplier,
      
      yend =
        NMDS2 *
        vector_multiplier,
      
      label =
        str_replace_all(
          variable,
          "_",
          " "
        )
    )
}


# =============================================================================
# Panel A: NMDS
# =============================================================================

panel_a <- ggplot(
  nmds_scores,
  aes(
    NMDS1,
    NMDS2
  )
) +
  
  geom_point(
    aes(
      color = Season,
      shape = Plot_ID
    ),
    size = 4.2,
    alpha = 0.90,
    stroke = 1
  )


if (
  nrow(
    environment_plot_vectors
  ) > 0
) {
  
  panel_a <- panel_a +
    
    geom_segment(
      data =
        environment_plot_vectors,
      
      aes(
        x = 0,
        y = 0,
        xend = xend,
        yend = yend
      ),
      
      inherit.aes =
        FALSE,
      
      arrow = arrow(
        length = unit(
          0.20,
          "cm"
        )
      ),
      
      color =
        "grey20",
      
      linewidth =
        0.9
    ) +
    
    geom_text_repel(
      data =
        environment_plot_vectors,
      
      aes(
        x = xend,
        y = yend,
        label = label
      ),
      
      inherit.aes =
        FALSE,
      
      color =
        "grey20",
      
      size =
        4.3,
      
      seed =
        42,
      
      box.padding =
        0.5,
      
      point.padding =
        0.4
    )
}


panel_a <- panel_a +
  
  annotate(
    "text",
    
    x = -Inf,
    y = Inf,
    
    hjust = -0.05,
    vjust = 1.10,
    
    label =
      nmds_stats_label,
    
    size =
      4.2
  ) +
  
  annotate(
    "text",
    
    x = Inf,
    y = -Inf,
    
    hjust = 1.05,
    vjust = -0.8,
    
    label = paste0(
      "Stress = ",
      sprintf(
        "%.3f",
        nmds_model$stress
      )
    ),
    
    size =
      4.2
  ) +
  
  scale_color_manual(
    values =
      season_colors,
    
    drop =
      FALSE
  ) +
  
  labs(
    title =
      "A",
    
    x =
      "NMDS1",
    
    y =
      "NMDS2"
  ) +
  
  coord_equal() +
  
  theme_manuscript +
  
  theme(
    legend.position =
      "none"
  )


# =============================================================================
# Temporal GAMs for all substrates
# =============================================================================

sir_long <- microresp_data %>%
  
  select(
    Sample_cat,
    Sampling_Date,
    Plot_ID,
    Block,
    Season,
    time_days,
    sequential_month,
    Substrate,
    Substrate_category,
    SIR
  )


prediction_dates <- sample_metadata %>%
  
  distinct(
    Sampling_Date,
    time_days
  ) %>%
  
  arrange(
    Sampling_Date
  )


dense_prediction_grid <- tibble(
  time_days = seq(
    min(
      sir_long$time_days
    ),
    
    max(
      sir_long$time_days
    ),
    
    length.out =
      250
  ),
  
  Sampling_Date =
    min(
      sir_long$Sampling_Date
    ) +
    round(
      time_days
    ),
  
  Plot_ID =
    factor(
      levels(
        sir_long$Plot_ID
      )[1],
      
      levels =
        levels(
          sir_long$Plot_ID
        )
    )
)


date_prediction_grid <- prediction_dates %>%
  
  mutate(
    Plot_ID =
      factor(
        levels(
          sir_long$Plot_ID
        )[1],
        
        levels =
          levels(
            sir_long$Plot_ID
          )
      )
  )


substrate_model_summary <-
  tibble()

substrate_dense_predictions <-
  tibble()

substrate_date_predictions <-
  tibble()


for (
  substrate_name in substrate_metadata$Substrate
) {
  
  substrate_data <- sir_long %>%
    
    filter(
      Substrate ==
        substrate_name
    )
  
  
  substrate_model <- gam(
    SIR ~
      s(
        time_days,
        k = 7
      ) +
      s(
        Plot_ID,
        bs = "re"
      ),
    
    data =
      substrate_data,
    
    method =
      "REML"
  )
  
  
  model_summary <-
    summary(
      substrate_model
    )
  
  
  smooth_row <- as.data.frame(
    model_summary$s.table
  ) %>%
    
    rownames_to_column(
      "term"
    ) %>%
    
    as_tibble() %>%
    
    filter(
      term ==
        "s(time_days)"
    )
  
  
  dense_prediction <- predict(
    substrate_model,
    
    newdata =
      dense_prediction_grid,
    
    type =
      "response",
    
    se.fit =
      TRUE,
    
    exclude =
      "s(Plot_ID)"
  )
  
  
  date_prediction <- predict(
    substrate_model,
    
    newdata =
      date_prediction_grid,
    
    type =
      "response",
    
    se.fit =
      TRUE,
    
    exclude =
      "s(Plot_ID)"
  )
  
  
  substrate_model_summary <- bind_rows(
    substrate_model_summary,
    
    tibble(
      Substrate =
        substrate_name,
      
      Substrate_category =
        first(
          substrate_data$Substrate_category
        ),
      
      n_observations =
        nrow(
          substrate_data
        ),
      
      n_plots =
        n_distinct(
          substrate_data$Plot_ID
        ),
      
      n_dates =
        n_distinct(
          substrate_data$Sampling_Date
        ),
      
      edf =
        smooth_row$edf,
      
      F_statistic =
        smooth_row$F,
      
      p_value =
        smooth_row$`p-value`,
      
      deviance_explained =
        model_summary$dev.expl
    )
  )
  
  
  substrate_dense_predictions <- bind_rows(
    substrate_dense_predictions,
    
    dense_prediction_grid %>%
      
      transmute(
        Substrate =
          substrate_name,
        
        Sampling_Date,
        
        time_days,
        
        fitted =
          as.numeric(
            dense_prediction$fit
          ),
        
        se =
          as.numeric(
            dense_prediction$se.fit
          ),
        
        lower_95 =
          fitted -
          1.96 *
          se,
        
        upper_95 =
          fitted +
          1.96 *
          se
      )
  )
  
  
  substrate_date_predictions <- bind_rows(
    substrate_date_predictions,
    
    date_prediction_grid %>%
      
      transmute(
        Substrate =
          substrate_name,
        
        Sampling_Date,
        
        time_days,
        
        fitted =
          as.numeric(
            date_prediction$fit
          ),
        
        se =
          as.numeric(
            date_prediction$se.fit
          ),
        
        lower_95 =
          fitted -
          1.96 *
          se,
        
        upper_95 =
          fitted +
          1.96 *
          se
      )
  )
}


# =============================================================================
# Multiple-testing correction and temporal ranking
# =============================================================================

substrate_model_summary <- substrate_model_summary %>%
  
  mutate(
    p_adj_FDR =
      p.adjust(
        p_value,
        method = "BH"
      )
  ) %>%
  
  arrange(
    desc(
      deviance_explained
    ),
    p_adj_FDR
  ) %>%
  
  mutate(
    temporal_rank =
      row_number()
  )


# =============================================================================
# Select top 3 temporal responses for main Figure 2
# =============================================================================
#
# Selection rule:
#
# 1. Temporal smooth must remain significant after FDR correction.
# 2. Among significant substrates, rank by temporal deviance explained.
# 3. Display the top three.
#
# Transcript behavior is NOT used in substrate selection.
# =============================================================================

top_temporal_substrates <- substrate_model_summary %>%
  
  filter(
    p_adj_FDR < 0.05
  ) %>%
  
  arrange(
    desc(
      deviance_explained
    ),
    p_adj_FDR
  ) %>%
  
  slice_head(
    n = 3
  ) %>%
  
  pull(
    Substrate
  )


if (
  length(
    top_temporal_substrates
  ) < 3
) {
  
  stop(
    "Fewer than three substrates are significant after FDR correction."
  )
}


substrate_selection_table <- substrate_model_summary %>%
  
  mutate(
    displayed_in_main_figure =
      Substrate %in%
      top_temporal_substrates,
    
    display_rank =
      match(
        Substrate,
        top_temporal_substrates
      ),
    
    selection_reason =
      case_when(
        displayed_in_main_figure ~
          paste0(
            "Top temporal response based on deviance explained among ",
            "FDR-significant substrate GAMs; display rank ",
            display_rank,
            "."
          ),
        
        TRUE ~
          paste0(
            "Not among the top three FDR-significant temporal responses ",
            "ranked by deviance explained."
          )
      )
  )


write_csv(
  substrate_selection_table,
  file.path(
    result_dir,
    "Figure2_substrate_temporal_ranking.csv"
  )
)


write_csv(
  substrate_dense_predictions,
  file.path(
    result_dir,
    "Figure2_substrate_fitted_dense_trajectories.csv"
  )
)


write_csv(
  substrate_date_predictions,
  file.path(
    result_dir,
    "Figure2_substrate_fitted_sampling_date_trajectories.csv"
  )
)


# =============================================================================
# Prepare top-3 data for main Figure 2
# =============================================================================

selected_plot_data <- sir_long %>%
  
  filter(
    Substrate %in%
      top_temporal_substrates
  ) %>%
  
  mutate(
    Substrate =
      factor(
        Substrate,
        levels =
          top_temporal_substrates
      ),
    
    display_label =
      str_replace_all(
        as.character(
          Substrate
        ),
        "[._]",
        " "
      )
  )


selected_predictions <- substrate_dense_predictions %>%
  
  filter(
    Substrate %in%
      top_temporal_substrates
  ) %>%
  
  mutate(
    Substrate =
      factor(
        Substrate,
        levels =
          top_temporal_substrates
      ),
    
    display_label =
      str_replace_all(
        as.character(
          Substrate
        ),
        "[._]",
        " "
      )
  )


# =============================================================================
# GAM statistics for displayed substrates
# =============================================================================

format_p_gam <- function(p) {
  
  case_when(
    is.na(
      p
    ) ~
      "NA",
    
    p < 0.001 ~
      "< 0.001",
    
    TRUE ~
      paste0(
        "= ",
        sprintf(
          "%.3f",
          p
        )
      )
  )
}


gam_plot_labels <- substrate_model_summary %>%
  
  filter(
    Substrate %in%
      top_temporal_substrates
  ) %>%
  
  mutate(
    Substrate =
      factor(
        Substrate,
        levels =
          top_temporal_substrates
      ),
    
    display_label =
      str_replace_all(
        as.character(
          Substrate
        ),
        "[._]",
        " "
      ),
    
    stat_label =
      paste0(
        "Deviance explained = ",
        sprintf(
          "%.1f",
          deviance_explained *
            100
        ),
        "%",
        "\nFDR p ",
        format_p_gam(
          p_adj_FDR
        )
      )
  )


# =============================================================================
# Seasonal shading
# =============================================================================

season_rectangles <- tribble(
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
# Panels B-D: top three temporal substrates
# =============================================================================

time_series_plot <- ggplot(
  selected_plot_data,
  
  aes(
    Sampling_Date,
    SIR
  )
) +
  
  geom_rect(
    data =
      season_rectangles,
    
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
      0.10,
    
    color =
      NA
  ) +
  
  geom_ribbon(
    data =
      selected_predictions,
    
    aes(
      x =
        Sampling_Date,
      
      ymin =
        lower_95,
      
      ymax =
        upper_95
    ),
    
    inherit.aes =
      FALSE,
    
    fill =
      "grey55",
    
    alpha =
      0.25
  ) +
  
  geom_line(
    data =
      selected_predictions,
    
    aes(
      x =
        Sampling_Date,
      
      y =
        fitted
    ),
    
    inherit.aes =
      FALSE,
    
    color =
      "black",
    
    linewidth =
      1.25
  ) +
  
  geom_point(
    aes(
      color =
        Season,
      
      shape =
        Plot_ID
    ),
    
    size =
      3,
    
    alpha =
      0.88,
    
    stroke =
      0.8
  ) +
  
  geom_text(
    data =
      gam_plot_labels,
    
    aes(
      x =
        Inf,
      
      y =
        Inf,
      
      label =
        stat_label
    ),
    
    inherit.aes =
      FALSE,
    
    hjust =
      1.05,
    
    vjust =
      1.05,
    
    size =
      3.7
  ) +
  
  facet_wrap(
    ~ display_label,
    
    nrow =
      1,
    
    scales =
      "free_y"
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
  
  scale_x_date(
    date_breaks =
      "3 months",
    
    date_labels =
      "%b\n%Y"
  ) +
  
  labs(
    title =
      "B-D",
    
    x =
      NULL,
    
    y =
      expression(
        "Potential SIR (ng CO"[2] *
          "-C " *
          mg^{-1} *
          " substrate-C " *
          h^{-1} *
          ")"
      ),
    
    color =
      "Season",
    
    shape =
      "Plot"
  ) +
  
  guides(
    fill =
      "none",
    
    color =
      guide_legend(
        order =
          1,
        
        override.aes =
          list(
            size = 4,
            alpha = 1
          )
      ),
    
    shape =
      guide_legend(
        order =
          2,
        
        nrow =
          1,
        
        override.aes =
          list(
            size = 4,
            alpha = 1
          )
      )
  ) +
  
  theme_manuscript +
  
  theme(
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1,
        size = 11
      ),
    
    strip.text =
      element_text(
        face = "bold",
        size = 15
      ),
    
    legend.position =
      "bottom"
  )


# =============================================================================
# Main Figure 2
# =============================================================================

figure_2 <- panel_a / time_series_plot +
  
  plot_layout(
    heights = c(
      1.05,
      1
    )
  )


ggsave(
  file.path(
    figure_dir,
    "Figure2_main.png"
  ),
  figure_2,
  width = 15,
  height = 10,
  dpi = 320
)


ggsave(
  file.path(
    figure_dir,
    "Figure2_main.pdf"
  ),
  figure_2,
  width = 15,
  height = 10
)


ggsave(
  file.path(
    figure_dir,
    "Figure2A_NMDS.png"
  ),
  panel_a,
  width = 8,
  height = 7,
  dpi = 320
)


ggsave(
  file.path(
    figure_dir,
    "Figure2A_NMDS.pdf"
  ),
  panel_a,
  width = 8,
  height = 7
)


ggsave(
  file.path(
    figure_dir,
    "Figure2B_D_top_temporal_substrates.png"
  ),
  time_series_plot,
  width = 15,
  height = 5.5,
  dpi = 320
)


ggsave(
  file.path(
    figure_dir,
    "Figure2B_D_top_temporal_substrates.pdf"
  ),
  time_series_plot,
  width = 15,
  height = 5.5
)


# =============================================================================
# SUPPLEMENTARY FIGURE
# All substrate-specific temporal trajectories
# =============================================================================

supp_plot_data <- sir_long %>%
  
  mutate(
    display_label =
      str_replace_all(
        Substrate,
        "[._]",
        " "
      )
  )


supp_predictions <- substrate_dense_predictions %>%
  
  mutate(
    display_label =
      str_replace_all(
        Substrate,
        "[._]",
        " "
      )
  )


supp_gam_labels <- substrate_model_summary %>%
  
  mutate(
    display_label =
      str_replace_all(
        Substrate,
        "[._]",
        " "
      ),
    
    stat_label =
      paste0(
        "Dev. = ",
        sprintf(
          "%.1f",
          deviance_explained *
            100
        ),
        "%",
        "\nFDR p ",
        format_p_gam(
          p_adj_FDR
        )
      )
  )


figure_s_temporal <- ggplot(
  supp_plot_data,
  
  aes(
    Sampling_Date,
    SIR
  )
) +
  
  geom_rect(
    data =
      season_rectangles,
    
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
      0.08,
    
    color =
      NA
  ) +
  
  geom_ribbon(
    data =
      supp_predictions,
    
    aes(
      x =
        Sampling_Date,
      
      ymin =
        lower_95,
      
      ymax =
        upper_95
    ),
    
    inherit.aes =
      FALSE,
    
    fill =
      "grey60",
    
    alpha =
      0.22
  ) +
  
  geom_line(
    data =
      supp_predictions,
    
    aes(
      x =
        Sampling_Date,
      
      y =
        fitted
    ),
    
    inherit.aes =
      FALSE,
    
    color =
      "black",
    
    linewidth =
      0.8
  ) +
  
  geom_point(
    aes(
      color =
        Season
    ),
    
    size =
      1.4,
    
    alpha =
      0.65
  ) +
  
  geom_text(
    data =
      supp_gam_labels,
    
    aes(
      x =
        Inf,
      
      y =
        Inf,
      
      label =
        stat_label
    ),
    
    inherit.aes =
      FALSE,
    
    hjust =
      1.03,
    
    vjust =
      1.05,
    
    size =
      2.4
  ) +
  
  facet_wrap(
    ~ display_label,
    
    ncol =
      4,
    
    scales =
      "free_y"
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
  
  scale_x_date(
    date_breaks =
      "6 months",
    
    date_labels =
      "%b\n%Y"
  ) +
  
  labs(
    x =
      NULL,
    
    y =
      expression(
        "Potential SIR (ng CO"[2] *
          "-C " *
          mg^{-1} *
          " substrate-C " *
          h^{-1} *
          ")"
      ),
    
    color =
      "Season"
  ) +
  
  guides(
    fill =
      "none"
  ) +
  
  theme_bw(
    base_size =
      10
  ) +
  
  theme(
    panel.grid.minor =
      element_blank(),
    
    strip.text =
      element_text(
        face = "bold",
        size = 9
      ),
    
    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1,
        size = 7
      ),
    
    axis.text.y =
      element_text(
        size = 7
      ),
    
    axis.title =
      element_text(
        size = 10
      ),
    
    legend.position =
      "bottom"
  )


ggsave(
  file.path(
    figure_dir,
    "FigureS_MicroResp_all_temporal_trajectories.png"
  ),
  figure_s_temporal,
  width = 14,
  height = 18,
  dpi = 320
)


ggsave(
  file.path(
    figure_dir,
    "FigureS_MicroResp_all_temporal_trajectories.pdf"
  ),
  figure_s_temporal,
  width = 14,
  height = 18
)


# =============================================================================
# Fitted-trajectory clustering
# =============================================================================

significant_substrates <- substrate_selection_table %>%
  
  filter(
    p_adj_FDR < 0.05
  ) %>%
  
  pull(
    Substrate
  )


trajectory_matrix <- substrate_date_predictions %>%
  
  filter(
    Substrate %in%
      significant_substrates
  ) %>%
  
  mutate(
    date_column =
      paste0(
        "date_",
        format(
          Sampling_Date,
          "%Y%m%d"
        )
      )
  ) %>%
  
  select(
    Substrate,
    date_column,
    fitted
  ) %>%
  
  pivot_wider(
    names_from =
      date_column,
    
    values_from =
      fitted
  ) %>%
  
  column_to_rownames(
    "Substrate"
  ) %>%
  
  as.matrix()


trajectory_sd <- apply(
  trajectory_matrix,
  1,
  sd
)


zero_variance_substrates <- tibble(
  Substrate =
    names(
      trajectory_sd
    ),
  
  fitted_trajectory_sd =
    trajectory_sd
) %>%
  
  filter(
    !is.finite(
      fitted_trajectory_sd
    ) |
      fitted_trajectory_sd ==
      0
  )


write_csv(
  zero_variance_substrates,
  file.path(
    result_dir,
    "Figure2_zero_variance_trajectories.csv"
  )
)


if (
  nrow(
    zero_variance_substrates
  ) > 0
) {
  
  trajectory_matrix <- trajectory_matrix[
    !rownames(
      trajectory_matrix
    ) %in%
      zero_variance_substrates$Substrate,
    ,
    drop = FALSE
  ]
}


zscore_matrix <- t(
  scale(
    t(
      trajectory_matrix
    )
  )
)


if (
  anyNA(
    zscore_matrix
  ) ||
  any(
    !is.finite(
      zscore_matrix
    )
  )
) {
  
  stop(
    "Z-scored fitted substrate trajectories contain nonfinite values."
  )
}


correlation_matrix <- cor(
  t(
    zscore_matrix
  ),
  method = "pearson"
)


correlation_validation <- tibble(
  check = c(
    "square",
    "symmetric",
    "finite",
    "diagonal_equals_one"
  ),
  
  passed = c(
    nrow(
      correlation_matrix
    ) ==
      ncol(
        correlation_matrix
      ),
    
    isTRUE(
      all.equal(
        correlation_matrix,
        t(
          correlation_matrix
        ),
        tolerance = 1e-10
      )
    ),
    
    all(
      is.finite(
        correlation_matrix
      )
    ),
    
    all(
      abs(
        diag(
          correlation_matrix
        ) -
          1
      ) <
        1e-10
    )
  )
)


write_csv(
  correlation_validation,
  file.path(
    result_dir,
    "Figure2_correlation_matrix_validation.csv"
  )
)


if (
  !all(
    correlation_validation$passed
  )
) {
  
  stop(
    "The substrate correlation matrix failed validation before as.dist()."
  )
}


correlation_distance_matrix <-
  1 -
  correlation_matrix


correlation_distance_matrix[
  correlation_distance_matrix <
    0
] <- 0


diag(
  correlation_distance_matrix
) <- 0


distance_object <- as.dist(
  correlation_distance_matrix
)


trajectory_hclust <- hclust(
  distance_object,
  method = "ward.D2"
)


k_max <- min(
  8,
  nrow(
    zscore_matrix
  ) -
    1
)


if (
  k_max <
  2
) {
  
  stop(
    "Fewer than three temporal substrate trajectories remain for clustering."
  )
}


silhouette_diagnostics <- tibble()


for (
  candidate_k in 2:k_max
) {
  
  candidate_clusters <- cutree(
    trajectory_hclust,
    k = candidate_k
  )
  
  
  candidate_silhouette <- silhouette(
    candidate_clusters,
    distance_object
  )
  
  
  silhouette_diagnostics <- bind_rows(
    silhouette_diagnostics,
    
    tibble(
      k =
        candidate_k,
      
      average_silhouette_width =
        mean(
          candidate_silhouette[
            ,
            "sil_width"
          ]
        ),
      
      minimum_silhouette_width =
        min(
          candidate_silhouette[
            ,
            "sil_width"
          ]
        ),
      
      n_negative_silhouette =
        sum(
          candidate_silhouette[
            ,
            "sil_width"
          ] <
            0
        )
    )
  )
}


selected_k <- silhouette_diagnostics %>%
  
  arrange(
    desc(
      average_silhouette_width
    ),
    k
  ) %>%
  
  slice_head(
    n = 1
  ) %>%
  
  pull(
    k
  )


silhouette_diagnostics <- silhouette_diagnostics %>%
  
  mutate(
    selected =
      k ==
      selected_k
  )


write_csv(
  silhouette_diagnostics,
  file.path(
    result_dir,
    "Figure2_clustering_silhouette_diagnostics.csv"
  )
)


raw_clusters <- cutree(
  trajectory_hclust,
  k = selected_k
)


cluster_assignments <- tibble(
  Substrate =
    names(
      raw_clusters
    ),
  
  raw_cluster =
    as.integer(
      raw_clusters
    ),
  
  peak_index =
    apply(
      zscore_matrix,
      1,
      which.max
    )
)


cluster_order <- cluster_assignments %>%
  
  group_by(
    raw_cluster
  ) %>%
  
  summarise(
    mean_peak_index =
      mean(
        peak_index
      ),
    
    .groups =
      "drop"
  ) %>%
  
  arrange(
    mean_peak_index,
    raw_cluster
  ) %>%
  
  mutate(
    Cluster_ID =
      paste0(
        "Cluster_",
        row_number()
      )
  )


cluster_assignments <- cluster_assignments %>%
  
  left_join(
    cluster_order,
    by = "raw_cluster"
  ) %>%
  
  left_join(
    substrate_metadata,
    by = "Substrate"
  ) %>%
  
  mutate(
    dendrogram_position =
      match(
        Substrate,
        
        trajectory_hclust$labels[
          trajectory_hclust$order
        ]
      )
  ) %>%
  
  arrange(
    Cluster_ID,
    dendrogram_position
  )


selected_silhouette <- silhouette(
  raw_clusters,
  distance_object
)


silhouette_table <- tibble(
  Substrate =
    names(
      raw_clusters
    ),
  
  silhouette_width =
    selected_silhouette[
      ,
      "sil_width"
    ]
)


cluster_assignments <- cluster_assignments %>%
  
  left_join(
    silhouette_table,
    by = "Substrate"
  )


write_csv(
  cluster_assignments %>%
    
    select(
      -dendrogram_position
    ),
  
  file.path(
    result_dir,
    "Figure2_substrate_cluster_assignments.csv"
  )
)


write_csv(
  as.data.frame(
    zscore_matrix
  ) %>%
    
    rownames_to_column(
      "Substrate"
    ),
  
  file.path(
    result_dir,
    "Figure2_zscored_fitted_trajectories.csv"
  )
)


write_csv(
  as.data.frame(
    correlation_matrix
  ) %>%
    
    rownames_to_column(
      "Substrate"
    ),
  
  file.path(
    result_dir,
    "Figure2_substrate_correlation_matrix.csv"
  )
)


# =============================================================================
# Clustering diagnostic plot
# =============================================================================

diagnostic_plot <- ggplot(
  silhouette_diagnostics,
  
  aes(
    k,
    average_silhouette_width
  )
) +
  
  geom_line(
    color =
      "grey30",
    
    linewidth =
      0.8
  ) +
  
  geom_point(
    aes(
      color =
        selected
    ),
    
    size =
      2.8
  ) +
  
  scale_color_manual(
    values = c(
      `FALSE` = "grey30",
      `TRUE` = "#D55E00"
    ),
    
    guide =
      "none"
  ) +
  
  scale_x_continuous(
    breaks =
      2:k_max
  ) +
  
  labs(
    title =
      "MicroResp fitted-trajectory clustering",
    
    subtitle =
      paste0(
        "Silhouette-selected k = ",
        selected_k
      ),
    
    x =
      "Number of clusters",
    
    y =
      "Average silhouette width"
  ) +
  
  theme_manuscript


ggsave(
  file.path(
    figure_dir,
    "FigureS2_clustering_silhouette_diagnostic.png"
  ),
  diagnostic_plot,
  width = 7,
  height = 5,
  dpi = 320
)


ggsave(
  file.path(
    figure_dir,
    "FigureS2_clustering_silhouette_diagnostic.pdf"
  ),
  diagnostic_plot,
  width = 7,
  height = 5
)


# =============================================================================
# Supplementary fitted-trajectory heatmap
# =============================================================================

ordered_substrates <-
  cluster_assignments$Substrate


heatmap_matrix <- zscore_matrix[
  ordered_substrates,
  ,
  drop = FALSE
]


colnames(
  heatmap_matrix
) <- format(
  prediction_dates$Sampling_Date,
  "%b'%y"
)


row_annotation <- cluster_assignments %>%
  
  select(
    Substrate,
    Cluster_ID,
    Substrate_category
  ) %>%
  
  column_to_rownames(
    "Substrate"
  )


cluster_palette <- setNames(
  hcl.colors(
    selected_k,
    "Dark 3"
  ),
  
  paste0(
    "Cluster_",
    seq_len(
      selected_k
    )
  )
)


category_levels <- sort(
  unique(
    substrate_metadata$Substrate_category
  )
)


category_palette <- setNames(
  hcl.colors(
    length(
      category_levels
    ),
    "Set 3"
  ),
  
  category_levels
)


annotation_colors <- list(
  Cluster_ID =
    cluster_palette,
  
  Substrate_category =
    category_palette
)


for (
  extension in c(
    "png",
    "pdf"
  )
) {
  
  pheatmap(
    heatmap_matrix,
    
    filename =
      file.path(
        figure_dir,
        
        paste0(
          "FigureS2_MicroResp_trajectory_heatmap.",
          extension
        )
      ),
    
    width =
      10,
    
    height =
      max(
        7,
        0.28 *
          nrow(
            heatmap_matrix
          )
      ),
    
    cluster_rows =
      FALSE,
    
    cluster_cols =
      FALSE,
    
    annotation_row =
      row_annotation,
    
    annotation_colors =
      annotation_colors,
    
    color =
      colorRampPalette(
        rev(
          brewer.pal(
            11,
            "RdBu"
          )
        )
      )(
        101
      ),
    
    breaks =
      seq(
        -3,
        3,
        length.out = 102
      ),
    
    border_color =
      NA,
    
    fontsize_row =
      8,
    
    fontsize_col =
      9,
    
    main =
      paste0(
        "100% irrigation fitted MicroResp trajectories; k = ",
        selected_k
      )
  )
}


# =============================================================================
# Analysis summary
# =============================================================================

analysis_summary <- tibble(
  item = c(
    "MicroResp observations",
    "MicroResp samples",
    "MicroResp plots",
    "MicroResp sampling dates",
    "Substrates",
    "NMDS stress",
    "Envfit samples with complete field data",
    "PERMANOVA permutation restriction",
    "Top temporal substrate 1",
    "Top temporal substrate 2",
    "Top temporal substrate 3",
    "FDR-significant temporal substrates",
    "Significant substrates clustered",
    "Zero-variance trajectories removed",
    "Silhouette-selected cluster number",
    "Two-cluster pattern supported"
  ),
  
  value = c(
    nrow(
      microresp_data
    ),
    
    n_distinct(
      microresp_data$Sample_cat
    ),
    
    n_distinct(
      microresp_data$Plot_ID
    ),
    
    n_distinct(
      microresp_data$Sampling_Date
    ),
    
    n_distinct(
      microresp_data$Substrate
    ),
    
    nmds_model$stress,
    
    nrow(
      envfit_complete
    ),
    
    "Permutations restricted within Plot_ID",
    
    top_temporal_substrates[1],
    
    top_temporal_substrates[2],
    
    top_temporal_substrates[3],
    
    sum(
      substrate_model_summary$p_adj_FDR <
        0.05,
      na.rm = TRUE
    ),
    
    nrow(
      zscore_matrix
    ),
    
    nrow(
      zero_variance_substrates
    ),
    
    selected_k,
    
    if_else(
      selected_k ==
        2,
      "Yes",
      "No"
    )
  )
)


write_csv(
  analysis_summary,
  file.path(
    result_dir,
    "Figure2_analysis_summary.csv"
  )
)


# =============================================================================
# Console output
# =============================================================================

cat(
  "\nFigure 2 analysis complete.\n"
)


cat(
  "\nAnalysis summary:\n"
)


print(
  analysis_summary
)


cat(
  "\nTop three temporal substrates displayed in Figure 2:\n"
)


print(
  substrate_selection_table %>%
    
    filter(
      displayed_in_main_figure
    ) %>%
    
    arrange(
      display_rank
    ) %>%
    
    select(
      display_rank,
      Substrate,
      Substrate_category,
      edf,
      F_statistic,
      deviance_explained,
      p_value,
      p_adj_FDR,
      selection_reason
    )
)


cat(
  "\nFull temporal ranking:\n"
)


print(
  substrate_selection_table %>%
    
    select(
      temporal_rank,
      Substrate,
      Substrate_category,
      edf,
      deviance_explained,
      p_adj_FDR,
      displayed_in_main_figure
    ) %>%
    
    arrange(
      temporal_rank
    )
)


cat(
  "\nClustering diagnostics:\n"
)


print(
  silhouette_diagnostics
)