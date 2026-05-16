# ============================================================================
# 03 Substrate-Level Temporal GAMs
# Q1 Step 2: Individual MicroResp substrate temporal models
# ============================================================================
#
# Objective:
#   Fit one GAMM per MicroResp substrate to identify substrates with
#   significant nonlinear temporal trends across the annual cycle. These fitted
#   trajectories are saved for Step 3 trajectory clustering.

rm(list = ls())

# ----------------------------------------------------------------------------
# Required packages
# ----------------------------------------------------------------------------

required_packages <- c("mgcv", "tidyverse", "broom")
missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

library(mgcv)
library(tidyverse)
library(broom)

theme_set(theme_bw(base_size = 11))

dir.create("Figures", showWarnings = FALSE, recursive = TRUE)
dir.create("Results", showWarnings = FALSE, recursive = TRUE)
dir.create("Results/Step2_diagnostics", showWarnings = FALSE, recursive = TRUE)

# ----------------------------------------------------------------------------
# Input paths
# ----------------------------------------------------------------------------

community_matrix_file <- "Results/Step1_community_matrix.csv"
env_data_file <- "Results/Step1_env_data_prepared.csv"

microresp_candidates <- c(
  "DataHub March/Soil_Phenotypes/SoilMicrobiomeSFA_prosser-time-series_MicroResp.csv",
)

env_candidates <- c(
  "Data/Data_merged_cleaned.csv",
  "Data/Data_merged_cleaned_complete_cases.csv")

microresp_file <- microresp_candidates[file.exists(microresp_candidates)][1]
fallback_env_file <- env_candidates[file.exists(env_candidates)][1]

if (is.na(microresp_file)) {
  stop(
    "Could not find raw MicroResp file for substrate metadata. Checked: ",
    paste(microresp_candidates, collapse = "; ")
  )
}

# ----------------------------------------------------------------------------
# Substrate metadata
# ----------------------------------------------------------------------------

microresp_raw <- read.csv(microresp_file, check.names = FALSE)

if ("Irrigaiton_trmt_cat" %in% names(microresp_raw) &&
    !"Irrigation_trmt_cat" %in% names(microresp_raw)) {
  microresp_raw <- microresp_raw %>%
    rename(Irrigation_trmt_cat = Irrigaiton_trmt_cat)
}

substrate_metadata <- microresp_raw %>%
  mutate(Timepoint = suppressWarnings(as.numeric(Timepoint))) %>%
  filter(
    Timepoint == 5,
    Normalization == "ngCO2-C_mgC_h",
    Substrate != "Water"
  ) %>%
  distinct(Substrate, Substrate_category) %>%
  arrange(Substrate_category, Substrate)

write_csv(substrate_metadata, "Results/Step2_diagnostics/substrate_metadata.csv")

# ----------------------------------------------------------------------------
# A. Prepare data in long format
# ----------------------------------------------------------------------------

if (file.exists(community_matrix_file) && file.exists(env_data_file)) {
  cat("\nUsing Step 1 outputs:\n")
  cat("Community matrix:", community_matrix_file, "\n")
  cat("Environmental metadata:", env_data_file, "\n")

  community_matrix <- read.csv(
    community_matrix_file,
    row.names = 1,
    check.names = FALSE
  )

  env_data <- read.csv(env_data_file, check.names = FALSE)
} else {
  cat("\nStep 1 outputs were not found. Reconstructing from raw inputs.\n")

  if (is.na(fallback_env_file)) {
    stop(
      "Could not find environmental fallback file. Checked: ",
      paste(env_candidates, collapse = "; ")
    )
  }

  microresp_filtered <- microresp_raw %>%
    mutate(
      Timepoint = suppressWarnings(as.numeric(Timepoint)),
      SIR = suppressWarnings(as.numeric(SIR))
    ) %>%
    filter(
      Timepoint == 5,
      Normalization == "ngCO2-C_mgC_h",
      Substrate != "Water"
    )

  community_wide <- microresp_filtered %>%
    group_by(Sample_cat, Substrate) %>%
    summarise(SIR = mean(SIR, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Substrate, values_from = SIR) %>%
    arrange(Sample_cat)

  substrate_cols <- setdiff(names(community_wide), "Sample_cat")
  community_wide[substrate_cols] <- community_wide[substrate_cols] %>%
    mutate(across(everything(), ~replace_na(., 0)))

  community_matrix <- community_wide %>%
    select(all_of(substrate_cols)) %>%
    as.data.frame()

  rownames(community_matrix) <- community_wide$Sample_cat

  env_data <- read.csv(fallback_env_file, check.names = FALSE) %>%
    filter(Sample_cat %in% rownames(community_matrix)) %>%
    mutate(
      Collection_date = as.Date(Collection_date),
      time_numeric = as.numeric(Collection_date - min(Collection_date, na.rm = TRUE)),
      year_month = lubridate::floor_date(Collection_date, "month"),
      month_label = format(year_month, "%b'%y"),
      month_order = as.integer(factor(year_month, levels = sort(unique(year_month)))),
      Irrigation_factor = factor(as.character(Irrigation_trmt_pct),
                                 levels = c("25", "100"))
    )

  common_samples <- intersect(rownames(community_matrix), env_data$Sample_cat)
  community_matrix <- community_matrix[common_samples, , drop = FALSE]
  env_data <- env_data %>%
    mutate(sample_order = match(Sample_cat, common_samples)) %>%
    arrange(sample_order) %>%
    select(-sample_order)
}

community_long <- community_matrix %>%
  as.data.frame() %>%
  rownames_to_column("Sample_cat") %>%
  pivot_longer(
    cols = -Sample_cat,
    names_to = "Substrate",
    values_to = "SIR"
  ) %>%
  mutate(SIR = as.numeric(SIR))

model_data <- community_long %>%
  left_join(env_data, by = "Sample_cat") %>%
  left_join(substrate_metadata, by = "Substrate") %>%
  mutate(
    Plot = factor(Plot),
    Irrigation_factor = factor(as.character(Irrigation_factor),
                               levels = c("25", "100")),
    time_numeric = as.numeric(time_numeric)
  )

model_required_cols <- c("SIR", "time_numeric", "Plot", "Irrigation_factor")

model_na_summary <- model_data %>%
  summarise(across(all_of(model_required_cols), ~sum(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "Column", values_to = "N_missing")

write_csv(model_na_summary, "Results/Step2_diagnostics/model_input_na_summary.csv")

if (any(model_na_summary$N_missing > 0)) {
  stop("Missing values found in SIR, time_numeric, Plot, or Irrigation_factor.")
}

cat("\n=== STEP 2 MODEL DATA CHECK ===\n")
cat("Rows:", nrow(model_data), "\n")
cat("Unique samples:", n_distinct(model_data$Sample_cat), "\n")
cat("Unique substrates:", n_distinct(model_data$Substrate), "\n")
cat("Expected substrates: 31\n")
print(model_na_summary)

# ----------------------------------------------------------------------------
# B. Fit one GAMM per substrate
# ----------------------------------------------------------------------------

substrate_order <- substrate_metadata %>%
  pull(Substrate) %>%
  intersect(unique(model_data$Substrate))

if (length(substrate_order) != 31) {
  warning("Expected 31 substrates for modeling, but found ", length(substrate_order), ".")
}

gam_models <- list()
prediction_list <- list()
summary_rows <- list()

time_grid <- env_data %>%
  distinct(time_numeric, Collection_date, month_order, month_label) %>%
  arrange(time_numeric)

reference_irrigation <- factor("25", levels = c("25", "100"))
reference_plot <- levels(model_data$Plot)[1]

cat("\n=== FITTING GAMMS BY SUBSTRATE ===\n")

for (substrate_name in substrate_order) {
  substrate_data <- model_data %>%
    filter(Substrate == substrate_name) %>%
    arrange(time_numeric, Sample_cat)

  substrate_category <- substrate_data %>%
    distinct(Substrate_category) %>%
    pull(Substrate_category) %>%
    first()

  model_warning <- NA_character_
  model_error <- NA_character_

  fitted_model <- tryCatch(
    withCallingHandlers(
      gamm(
        SIR ~ s(time_numeric, k = 7) + Irrigation_factor,
        random = list(Plot = ~1),
        data = substrate_data,
        method = "REML"
      ),
      warning = function(w) {
        model_warning <<- conditionMessage(w)
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      model_error <<- conditionMessage(e)
      NULL
    }
  )

  if (is.null(fitted_model)) {
    summary_rows[[substrate_name]] <- tibble(
      Substrate = substrate_name,
      Substrate_category = substrate_category,
      EDF = NA_real_,
      F_statistic = NA_real_,
      p_value_smooth = NA_real_,
      Deviance_explained = NA_real_,
      Adj_R2 = NA_real_,
      Irrigation_coef = NA_real_,
      Irrigation_p = NA_real_,
      Scale = NA_real_,
      Converged = "No",
      Model_warning = model_warning,
      Model_error = model_error
    )

    cat("  FAILED:", substrate_name, "-", model_error, "\n")
    next
  }

  gam_models[[substrate_name]] <- fitted_model
  gam_summary <- summary(fitted_model$gam)

  smooth_table <- as.data.frame(gam_summary$s.table)
  smooth_row <- smooth_table[rownames(smooth_table) == "s(time_numeric)", , drop = FALSE]

  param_table <- as.data.frame(gam_summary$p.table)
  irrigation_row_name <- grep("^Irrigation_factor", rownames(param_table), value = TRUE)[1]

  irrigation_coef <- if (!is.na(irrigation_row_name)) {
    param_table[irrigation_row_name, "Estimate"]
  } else {
    NA_real_
  }

  irrigation_p <- if (!is.na(irrigation_row_name)) {
    param_table[irrigation_row_name, "Pr(>|t|)"]
  } else {
    NA_real_
  }

  fitted_values <- as.numeric(fitted(fitted_model$lme, level = 1))
  residual_sum_squares <- sum((substrate_data$SIR - fitted_values)^2, na.rm = TRUE)
  total_sum_squares <- sum((substrate_data$SIR - mean(substrate_data$SIR, na.rm = TRUE))^2,
                           na.rm = TRUE)
  deviance_explained <- if_else(
    total_sum_squares > 0,
    1 - residual_sum_squares / total_sum_squares,
    NA_real_
  )

  summary_rows[[substrate_name]] <- tibble(
    Substrate = substrate_name,
    Substrate_category = substrate_category,
    EDF = smooth_row$edf,
    F_statistic = smooth_row$F,
    p_value_smooth = smooth_row$`p-value`,
    Deviance_explained = deviance_explained,
    Adj_R2 = gam_summary$r.sq,
    Irrigation_coef = irrigation_coef,
    Irrigation_p = irrigation_p,
    Scale = gam_summary$scale,
    Converged = "Yes",
    Model_warning = model_warning,
    Model_error = model_error
  )

  prediction_data <- time_grid %>%
    mutate(
      Irrigation_factor = reference_irrigation,
      Plot = factor(reference_plot, levels = levels(model_data$Plot))
    )

  predicted <- predict(
    fitted_model$gam,
    newdata = prediction_data,
    type = "response",
    se.fit = TRUE
  )

  prediction_list[[substrate_name]] <- prediction_data %>%
    mutate(
      Substrate = substrate_name,
      Substrate_category = substrate_category,
      fitted_SIR = as.numeric(predicted$fit),
      SE = as.numeric(predicted$se.fit)
    ) %>%
    select(Substrate, Substrate_category, time_numeric, Collection_date,
           month_order, month_label, fitted_SIR, SE)

  cat("  Fit:", substrate_name, "\n")
}

gam_summary_table <- bind_rows(summary_rows) %>%
  mutate(
    p_adj_FDR = p.adjust(p_value_smooth, method = "BH"),
    Significant = if_else(!is.na(p_adj_FDR) & p_adj_FDR < 0.05, "Yes", "No")
  ) %>%
  select(
    Substrate,
    Substrate_category,
    EDF,
    F_statistic,
    p_value_smooth,
    p_adj_FDR,
    Significant,
    Deviance_explained,
    Adj_R2,
    Irrigation_coef,
    Irrigation_p,
    Scale,
    Converged,
    Model_warning,
    Model_error
  ) %>%
  arrange(p_adj_FDR, Substrate)

gam_trajectories <- bind_rows(prediction_list) %>%
  left_join(
    gam_summary_table %>% select(Substrate, p_adj_FDR, Significant),
    by = "Substrate"
  )

failed_models <- gam_summary_table %>%
  filter(Converged == "No")

write_csv(gam_summary_table, "Results/Step2_GAM_summary_table.csv")
write_csv(gam_trajectories, "Results/Step2_GAM_fitted_trajectories.csv")
write_csv(model_data, "Results/Step2_GAM_long_model_data.csv")
write_csv(failed_models, "Results/Step2_diagnostics/failed_gamm_models.csv")

# ----------------------------------------------------------------------------
# C-F. Multiple testing correction and summary statistics
# ----------------------------------------------------------------------------

significant_substrates <- gam_summary_table %>%
  filter(Significant == "Yes")

significant_irrigation <- gam_summary_table %>%
  filter(!is.na(Irrigation_p), Irrigation_p < 0.05)

cat("\n=== STEP 2 GAMM SUMMARY ===\n")
cat("Significant temporal trends (FDR < 0.05):",
    nrow(significant_substrates), "of", nrow(gam_summary_table), "\n")

if (nrow(significant_substrates) > 0) {
  cat("EDF range among significant substrates:",
      round(min(significant_substrates$EDF, na.rm = TRUE), 3), "to",
      round(max(significant_substrates$EDF, na.rm = TRUE), 3),
      "; median =", round(median(significant_substrates$EDF, na.rm = TRUE), 3),
      "\n")

  cat("Deviance explained range among significant substrates:",
      round(min(significant_substrates$Deviance_explained, na.rm = TRUE), 3), "to",
      round(max(significant_substrates$Deviance_explained, na.rm = TRUE), 3),
      "\n")
} else {
  cat("No significant temporal trends after FDR correction.\n")
}

cat("Substrates with significant irrigation effect (p < 0.05):",
    nrow(significant_irrigation), "\n")
cat("Models failed to converge:", nrow(failed_models), "\n")

step2_summary <- tibble(
  Metric = c(
    "Number of substrates modeled",
    "Significant temporal trends, FDR < 0.05",
    "Minimum EDF among significant substrates",
    "Median EDF among significant substrates",
    "Maximum EDF among significant substrates",
    "Minimum deviance explained among significant substrates",
    "Maximum deviance explained among significant substrates",
    "Substrates with irrigation p < 0.05",
    "Models failed to converge"
  ),
  Value = c(
    as.character(nrow(gam_summary_table)),
    as.character(nrow(significant_substrates)),
    if (nrow(significant_substrates) > 0) {
      as.character(round(min(significant_substrates$EDF, na.rm = TRUE), 4))
    } else {
      NA_character_
    },
    if (nrow(significant_substrates) > 0) {
      as.character(round(median(significant_substrates$EDF, na.rm = TRUE), 4))
    } else {
      NA_character_
    },
    if (nrow(significant_substrates) > 0) {
      as.character(round(max(significant_substrates$EDF, na.rm = TRUE), 4))
    } else {
      NA_character_
    },
    if (nrow(significant_substrates) > 0) {
      as.character(round(min(significant_substrates$Deviance_explained,
                             na.rm = TRUE), 4))
    } else {
      NA_character_
    },
    if (nrow(significant_substrates) > 0) {
      as.character(round(max(significant_substrates$Deviance_explained,
                             na.rm = TRUE), 4))
    } else {
      NA_character_
    },
    as.character(nrow(significant_irrigation)),
    as.character(nrow(failed_models))
  )
)

write_csv(step2_summary, "Results/Step2_GAM_summary_statistics.csv")

# ----------------------------------------------------------------------------
# G. Supplementary figure: 31-panel GAM curves
# ----------------------------------------------------------------------------

category_order <- substrate_metadata %>%
  distinct(Substrate_category) %>%
  arrange(Substrate_category) %>%
  pull(Substrate_category)

panel_order <- gam_summary_table %>%
  mutate(
    Substrate_category = factor(Substrate_category, levels = category_order),
    Significant_order = if_else(Significant == "Yes", 0, 1)
  ) %>%
  arrange(Substrate_category, Significant_order, p_adj_FDR, Substrate) %>%
  pull(Substrate)

plot_data <- model_data %>%
  left_join(
    gam_summary_table %>%
      select(Substrate, p_adj_FDR, Significant, EDF),
    by = "Substrate"
  ) %>%
  mutate(
    Substrate = factor(Substrate, levels = panel_order),
    Significant = factor(Significant, levels = c("Yes", "No"))
  )

trajectory_plot_data <- gam_trajectories %>%
  left_join(
    gam_summary_table %>% select(Substrate, EDF),
    by = "Substrate"
  ) %>%
  mutate(
    Substrate = factor(Substrate, levels = panel_order),
    Significant = factor(Significant, levels = c("Yes", "No")),
    lower = fitted_SIR - 1.96 * SE,
    upper = fitted_SIR + 1.96 * SE
  )

annotation_data <- gam_summary_table %>%
  mutate(
    Substrate = factor(Substrate, levels = panel_order),
    label = paste0(
      "EDF=", round(EDF, 2),
      ", FDR=", signif(p_adj_FDR, 2)
    )
  ) %>%
  left_join(
    plot_data %>%
      group_by(Substrate) %>%
      summarise(
        x = min(time_numeric, na.rm = TRUE),
        y = max(SIR, na.rm = TRUE),
        .groups = "drop"
      ),
    by = "Substrate"
  )

month_breaks <- time_grid$time_numeric
month_labels <- time_grid$month_label

non_significant_panel_background <- gam_summary_table %>%
  filter(Significant == "No") %>%
  mutate(Substrate = factor(Substrate, levels = panel_order))

all_gam_figure <- ggplot()

if (nrow(non_significant_panel_background) > 0) {
  all_gam_figure <- all_gam_figure +
    geom_rect(
      data = non_significant_panel_background,
      aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
      fill = "grey95",
      inherit.aes = FALSE
    )
}

all_gam_figure <- all_gam_figure +
  geom_point(
    data = plot_data,
    aes(x = time_numeric, y = SIR, color = Irrigation_factor),
    size = 0.8,
    alpha = 0.45
  ) +
  geom_ribbon(
    data = trajectory_plot_data,
    aes(x = time_numeric, ymin = lower, ymax = upper),
    fill = "grey55",
    alpha = 0.22
  ) +
  geom_line(
    data = trajectory_plot_data,
    aes(x = time_numeric, y = fitted_SIR, linetype = Significant),
    color = "black",
    linewidth = 0.55
  ) +
  geom_text(
    data = annotation_data,
    aes(x = x, y = y, label = label),
    hjust = 0,
    vjust = 1.1,
    size = 2.3,
    color = "black"
  ) +
  facet_wrap(~Substrate, scales = "free_y", ncol = 6) +
  scale_color_manual(
    name = "Irrigation",
    values = c("100" = "#0072B2", "25" = "#D55E00"),
    breaks = c("100", "25"),
    labels = c("100%", "25%")
  ) +
  scale_linetype_manual(
    values = c("Yes" = "solid", "No" = "dashed"),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = month_breaks,
    labels = month_labels
  ) +
  labs(
    x = "Month",
    y = "SIR (ng CO2-C / mg substrate-C / h)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    strip.text = element_text(size = 8),
    legend.position = "bottom"
  )

ggsave(
  "Figures/FigS_Step2_all_substrate_GAMs.png",
  all_gam_figure,
  width = 16,
  height = 14,
  dpi = 300
)

ggsave(
  "Figures/FigS_Step2_all_substrate_GAMs.pdf",
  all_gam_figure,
  width = 16,
  height = 14
)

# ----------------------------------------------------------------------------
# G. Summary significance figure
# ----------------------------------------------------------------------------

significance_plot_data <- gam_summary_table %>%
  mutate(
    p_adj_for_plot = pmax(p_adj_FDR, 1e-16),
    neg_log10_FDR = -log10(p_adj_for_plot),
    Substrate = fct_reorder(Substrate, neg_log10_FDR),
    Significant = factor(Significant, levels = c("No", "Yes"))
  )

significance_figure <- ggplot(
  significance_plot_data,
  aes(x = neg_log10_FDR, y = Substrate, fill = Substrate_category)
) +
  geom_col(width = 0.7, alpha = 0.85) +
  geom_vline(
    xintercept = -log10(0.05),
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  geom_point(
    aes(shape = Significant),
    x = 0.05,
    size = 2,
    color = "black"
  ) +
  scale_shape_manual(values = c("No" = 1, "Yes" = 16), name = "FDR < 0.05") +
  labs(
    x = "-log10(FDR-adjusted p-value; capped at 16)",
    y = NULL,
    fill = "Substrate category"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

ggsave(
  "Figures/Figure_Step2_GAM_significance.png",
  significance_figure,
  width = 9,
  height = 8,
  dpi = 300
)

ggsave(
  "Figures/Figure_Step2_GAM_significance.pdf",
  significance_figure,
  width = 9,
  height = 8
)

cat("\n=== FILES SAVED ===\n")
cat("Results/Step2_GAM_summary_table.csv\n")
cat("Results/Step2_GAM_fitted_trajectories.csv\n")
cat("Results/Step2_GAM_summary_statistics.csv\n")
cat("Figures/FigS_Step2_all_substrate_GAMs.png/pdf\n")
cat("Figures/Figure_Step2_GAM_significance.png/pdf\n")

cat("\nQ1 Step 2 substrate-level GAMM analysis complete.\n")
