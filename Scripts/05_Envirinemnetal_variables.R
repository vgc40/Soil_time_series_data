
############################################################
# Step 4: Environmental Attribution
# Q1 MicroResp temporal trajectory analysis
#
# Objective:
# Quantify how much variation in each temporal trend-group score
# is explained by measured environmental covariates.
############################################################

# -----------------------------
# 0. Setup
# -----------------------------

required_packages <- c(
  "tidyverse",
  "lubridate",
  "mgcv",
  "car",
  "broom",
  "patchwork",
  "RColorBrewer"
)

installed <- rownames(installed.packages())
for (pkg in required_packages) {
  if (!pkg %in% installed) {
    install.packages(pkg, dependencies = TRUE)
  }
}

library(tidyverse)
library(lubridate)
library(mgcv)
library(car)
library(broom)
library(patchwork)
library(RColorBrewer)

dir.create("Results", showWarnings = FALSE, recursive = TRUE)
dir.create("Figures", showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1. Input files
# -----------------------------


microresp_file <- "DataHub March/Soil_Phenotypes/SoilMicrobiomeSFA_prosser-time-series_MicroResp.csv"
env_file <- "Data/Data_merged_cleaned_complete_cases.csv"
cluster_file <- "Results/Step3_cluster_assignments_FINAL_correlation_wardD2.csv"

microresp_raw <- read_csv(microresp_file, show_col_types = FALSE)
env_raw <- read_csv(env_file, show_col_types = FALSE)
cluster_assignments <- read_csv(cluster_file, show_col_types = FALSE)

# Expected cluster columns:
# Substrate, Substrate_category, Cluster_ID, Silhouette_width

# -----------------------------
# 2. Prepare MicroResp data
# -----------------------------

microresp <- microresp_raw %>%
  filter(
    Timepoint == 5,
    Normalization == "ngCO2-C_mgC_h",
    Substrate != "Water"
  ) %>%
  mutate(
    Collection_date = lubridate::mdy(Collection_date),
    Plot = factor(Plot),
    Irrigation_factor = factor(
      Irrigation_trmt_pct,
      levels = c(25, 100),
      labels = c("25", "100")
    ),
    SIR = as.numeric(SIR)
  ) %>%
  filter(
    !is.na(Sample_cat),
    !is.na(Collection_date),
    !is.na(Plot),
    !is.na(Irrigation_factor),
    !is.na(Substrate),
    !is.na(SIR)
  )

first_date <- min(microresp$Collection_date, na.rm = TRUE)

microresp <- microresp %>%
  mutate(
    time_numeric = as.numeric(Collection_date - first_date)
  )

# -----------------------------
# 3. Compute per-sample trend-group scores
# -----------------------------
# For each sample and cluster:
# 1. log1p(SIR)
# 2. z-score within substrate across all samples
# 3. average z-scored values across substrates in the cluster

clustered_substrates <- cluster_assignments %>%
  select(Substrate, Substrate_category, Cluster_ID)

microresp_clustered <- microresp %>%
  inner_join(clustered_substrates, by = "Substrate")

if (nrow(microresp_clustered) == 0) {
  stop("No MicroResp rows matched Step 3 cluster assignments.")
}

substrate_scaling <- microresp_clustered %>%
  mutate(log1p_SIR = log1p(SIR)) %>%
  group_by(Substrate) %>%
  summarise(
    substrate_log_mean = mean(log1p_SIR, na.rm = TRUE),
    substrate_log_sd = sd(log1p_SIR, na.rm = TRUE),
    .groups = "drop"
  )

zero_sd_substrates <- substrate_scaling %>%
  filter(is.na(substrate_log_sd) | substrate_log_sd == 0) %>%
  pull(Substrate)

if (length(zero_sd_substrates) > 0) {
  warning(
    "Removing clustered substrates with zero or missing SD after log1p transform: ",
    paste(zero_sd_substrates, collapse = ", ")
  )
}

microresp_z <- microresp_clustered %>%
  mutate(log1p_SIR = log1p(SIR)) %>%
  left_join(substrate_scaling, by = "Substrate") %>%
  filter(
    !Substrate %in% zero_sd_substrates,
    !is.na(substrate_log_mean),
    !is.na(substrate_log_sd),
    substrate_log_sd > 0
  ) %>%
  mutate(
    substrate_z = (log1p_SIR - substrate_log_mean) / substrate_log_sd
  )

group_scores_long <- microresp_z %>%
  group_by(
    Sample_cat,
    Collection_date,
    time_numeric,
    Plot,
    Irrigation_factor,
    Irrigation_trmt_pct,
    Cluster_ID
  ) %>%
  summarise(
    group_score = mean(substrate_z, na.rm = TRUE),
    n_substrates_used = n_distinct(Substrate),
    .groups = "drop"
  )

write_csv(
  group_scores_long,
  "Results/Step4_group_scores_long.csv"
)

group_scores_wide <- group_scores_long %>%
  select(Sample_cat, Cluster_ID, group_score) %>%
  pivot_wider(
    names_from = Cluster_ID,
    values_from = group_score
  )

write_csv(
  group_scores_wide,
  "Results/Step4_group_scores_wide.csv"
)

# -----------------------------
# 4. Prepare environmental data
# -----------------------------

env <- env_raw %>%
  mutate(
    Collection_date = lubridate::ymd(Collection_date),
    Plot = factor(Plot),
    Irrigation_factor = factor(
      Irrigation_trmt_pct,
      levels = c(25, 100),
      labels = c("25", "100")
    ),
    SEOC_SETN = SEOC_ug_gdw / SETN_ug_gdw
  ) %>%
  select(
    Sample_cat,
    Collection_date,
    Plot,
    Irrigation_factor,
    Irrigation_trmt_pct,
    GWC,
    pH,
    SEOC_ug_gdw,
    SETN_ug_gdw,
    SEOC_SETN,
    MBC_ug_gdw_bgc,
    MBN_ug_gdw_bgc,
    SoilTemp_mean,
    Roots_gdw
  )

# -----------------------------
# 5. Merge group scores with environmental data
# -----------------------------

step4_data <- group_scores_long %>%
  left_join(
    env,
    by = c(
      "Sample_cat",
      "Collection_date",
      "Plot",
      "Irrigation_factor",
      "Irrigation_trmt_pct"
    )
  ) %>%
  filter(
    !is.na(group_score),
    !is.na(time_numeric),
    !is.na(SoilTemp_mean),
    !is.na(GWC),
    !is.na(SEOC_SETN),
    !is.na(MBC_ug_gdw_bgc),
    !is.na(Irrigation_factor),
    !is.na(Plot)
  )

write_csv(
  step4_data,
  "Results/Step4_modeling_dataset.csv"
)

if (nrow(step4_data) == 0) {
  stop("No complete cases available after merging group scores with environmental data.")
}

# -----------------------------
# 6. Helper functions
# -----------------------------

extract_model_stats <- function(model, cluster_id, model_name) {
  sm <- summary(model)
  
  tibble(
    Cluster_ID = cluster_id,
    Model = model_name,
    Deviance_explained = sm$dev.expl * 100,
    AIC = AIC(model),
    Adj_R2 = sm$r.sq
  )
}

extract_model_b_coefficients <- function(model, cluster_id) {
  coef_table <- summary(model)$p.table
  
  coef_df <- as.data.frame(coef_table) %>%
    rownames_to_column("Variable") %>%
    as_tibble()
  
  names(coef_df) <- c("Variable", "Coefficient", "SE", "t_value", "p_value")
  
  coef_df %>%
    filter(Variable != "(Intercept)") %>%
    mutate(
      Cluster_ID = cluster_id,
      CI_lower = Coefficient - 1.96 * SE,
      CI_upper = Coefficient + 1.96 * SE,
      Significant = if_else(p_value < 0.05, "Yes", "No")
    ) %>%
    select(
      Cluster_ID,
      Variable,
      Coefficient,
      SE,
      t_value,
      p_value,
      CI_lower,
      CI_upper,
      Significant
    )
}

compute_vif <- function(dat, cluster_id) {
  vif_lm <- lm(
    group_score ~ SoilTemp_mean + GWC + SEOC_SETN + MBC_ug_gdw_bgc + Irrigation_factor,
    data = dat
  )
  
  vif_raw <- car::vif(vif_lm)
  
  if (is.matrix(vif_raw)) {
    vif_df <- as.data.frame(vif_raw) %>%
      rownames_to_column("Variable") %>%
      as_tibble() %>%
      mutate(
        VIF = GVIF^(1 / (2 * Df))
      ) %>%
      select(Variable, VIF)
  } else {
    vif_df <- tibble(
      Variable = names(vif_raw),
      VIF = as.numeric(vif_raw)
    )
  }
  
  vif_df %>%
    mutate(
      Cluster_ID = cluster_id,
      VIF_flag = if_else(VIF >= 5, "High", "OK")
    ) %>%
    select(Cluster_ID, Variable, VIF, VIF_flag)
}

fit_cluster_models <- function(dat, cluster_id) {
  dat_cluster <- dat %>%
    filter(Cluster_ID == cluster_id) %>%
    mutate(
      Plot = factor(Plot),
      Irrigation_factor = factor(Irrigation_factor, levels = c("25", "100"))
    )
  
  n_complete <- nrow(dat_cluster)
  
  if (n_complete < 20) {
    warning("Cluster ", cluster_id, " has fewer than 20 complete observations.")
  }
  
  # Model A: time only
  model_a <- gam(
    group_score ~
      s(time_numeric, k = 7) +
      s(Plot, bs = "re"),
    data = dat_cluster,
    method = "REML"
  )
  
  # Model B: time + environment
  model_b <- gam(
    group_score ~
      s(time_numeric, k = 7) +
      SoilTemp_mean +
      GWC +
      SEOC_SETN +
      MBC_ug_gdw_bgc +
      Irrigation_factor +
      s(Plot, bs = "re"),
    data = dat_cluster,
    method = "REML"
  )
  
  # Model C: environment only
  model_c <- gam(
    group_score ~
      SoilTemp_mean +
      GWC +
      SEOC_SETN +
      MBC_ug_gdw_bgc +
      Irrigation_factor +
      s(Plot, bs = "re"),
    data = dat_cluster,
    method = "REML"
  )
  
  stats <- bind_rows(
    extract_model_stats(model_a, cluster_id, "A_time_only"),
    extract_model_stats(model_b, cluster_id, "B_time_environment"),
    extract_model_stats(model_c, cluster_id, "C_environment_only")
  ) %>%
    group_by(Cluster_ID) %>%
    mutate(
      AIC_Model_A = AIC[Model == "A_time_only"],
      delta_AIC_vs_A = AIC - AIC_Model_A
    ) %>%
    ungroup() %>%
    select(-AIC_Model_A)
  
  coefficients <- extract_model_b_coefficients(model_b, cluster_id)
  
  vif_df <- compute_vif(dat_cluster, cluster_id)
  
  list(
    model_a = model_a,
    model_b = model_b,
    model_c = model_c,
    stats = stats,
    coefficients = coefficients,
    vif = vif_df
  )
}

# -----------------------------
# 7. Fit models for each cluster
# -----------------------------

cluster_ids <- sort(unique(step4_data$Cluster_ID))

model_results <- vector("list", length(cluster_ids))
names(model_results) <- cluster_ids

for (cluster_id in cluster_ids) {
  message("Fitting Step 4 models for ", cluster_id)
  
  model_results[[cluster_id]] <- fit_cluster_models(
    dat = step4_data,
    cluster_id = cluster_id
  )
}

model_comparison_table <- map_dfr(model_results, "stats")
environmental_coefficients <- map_dfr(model_results, "coefficients")
vif_check <- map_dfr(model_results, "vif")

write_csv(
  model_comparison_table,
  "Results/Step4_model_comparison_table.csv"
)

write_csv(
  environmental_coefficients,
  "Results/Step4_environmental_coefficients.csv"
)

write_csv(
  vif_check,
  "Results/Step4_VIF_check.csv"
)

# -----------------------------
# 8. Optional diagnostic summaries
# -----------------------------

model_diagnostics <- map_dfr(
  names(model_results),
  function(cluster_id) {
    dat_cluster <- step4_data %>%
      filter(Cluster_ID == cluster_id)
    
    tibble(
      Cluster_ID = cluster_id,
      n_observations = nrow(dat_cluster),
      n_plots = n_distinct(dat_cluster$Plot),
      n_timepoints = n_distinct(dat_cluster$time_numeric),
      n_substrates_in_cluster = cluster_assignments %>%
        filter(Cluster_ID == cluster_id) %>%
        nrow()
    )
  }
)

write_csv(
  model_diagnostics,
  "Results/Step4_model_diagnostics.csv"
)

# -----------------------------
# 9. Figure 3A: Variance explained
# -----------------------------

variance_plot_data <- model_comparison_table %>%
  mutate(
    Model = factor(
      Model,
      levels = c(
        "A_time_only",
        "B_time_environment",
        "C_environment_only"
      ),
      labels = c(
        "A: Time only",
        "B: Time + environment",
        "C: Environment only"
      )
    ),
    Cluster_ID = factor(Cluster_ID, levels = cluster_ids)
  )

p_var <- ggplot(
  variance_plot_data,
  aes(
    x = Cluster_ID,
    y = Deviance_explained,
    fill = Model
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  labs(
    x = "Trend group",
    y = "Deviance explained (%)",
    fill = "Model"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(
  filename = "Figures/Figure3A_variance_explained.png",
  plot = p_var,
  width = 8,
  height = 5,
  dpi = 300
)

# -----------------------------
# 10. Figure 3B: Environmental coefficients
# -----------------------------
# Uses case_when() rather than recode() to avoid conflict with car::recode().

coef_plot_data <- environmental_coefficients %>%
  filter(
    Variable %in% c(
      "SoilTemp_mean",
      "GWC",
      "SEOC_SETN",
      "MBC_ug_gdw_bgc",
      "Irrigation_factor100"
    )
  ) %>%
  mutate(
    Variable = case_when(
      Variable == "SoilTemp_mean" ~ "Soil temperature",
      Variable == "GWC" ~ "Gravimetric water content",
      Variable == "SEOC_SETN" ~ "SEOC:SETN",
      Variable == "MBC_ug_gdw_bgc" ~ "Microbial biomass C",
      Variable == "Irrigation_factor100" ~ "Irrigation: 100% vs 25%",
      TRUE ~ Variable
    ),
    Variable = factor(
      Variable,
      levels = c(
        "Soil temperature",
        "Gravimetric water content",
        "SEOC:SETN",
        "Microbial biomass C",
        "Irrigation: 100% vs 25%"
      )
    ),
    Cluster_ID = factor(Cluster_ID, levels = cluster_ids)
  )

p_coef <- ggplot(
  coef_plot_data,
  aes(
    x = Coefficient,
    y = Variable,
    color = Significant
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_errorbarh(
    aes(xmin = CI_lower, xmax = CI_upper),
    height = 0.2
  ) +
  geom_point(size = 2.3) +
  facet_wrap(~ Cluster_ID, scales = "free_x") +
  labs(
    x = "Model B coefficient estimate ± 95% CI",
    y = NULL,
    color = "p < 0.05"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

ggsave(
  filename = "Figures/Figure3B_coefficient_plots.png",
  plot = p_coef,
  width = 10,
  height = 6,
  dpi = 300
)

# -----------------------------
# 11. Combined Figure 3
# -----------------------------

p_combined <- p_var / p_coef +
  plot_annotation(
    tag_levels = "A"
  )

ggsave(
  filename = "Figures/Figure3_combined.png",
  plot = p_combined,
  width = 10,
  height = 11,
  dpi = 300
)

ggsave(
  filename = "Figures/Figure3_combined.pdf",
  plot = p_combined,
  width = 10,
  height = 11
)

# -----------------------------
# 12. Console summary
# -----------------------------

cat("\nStep 4 complete.\n")
cat("Clusters modeled:", length(cluster_ids), "\n")
cat("Modeling rows:", nrow(step4_data), "\n\n")

cat("Model comparison summary:\n")
print(model_comparison_table)

cat("\nEnvironmental coefficients from Model B:\n")
print(environmental_coefficients)

cat("\nVIF check:\n")
print(vif_check)

high_vif <- vif_check %>%
  filter(VIF >= 5)

if (nrow(high_vif) > 0) {
  cat("\nWarning: VIF >= 5 detected for the following predictors:\n")
  print(high_vif)
} else {
  cat("\nNo VIF values >= 5 detected.\n")
}

cat("\nOutputs written to Results/ and Figures/.\n")