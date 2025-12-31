################################################################################
# Random Forest Manuscript Analysis (NO Month) + LOPO CV + Bootstrap Importance
# + TRUE SHAP (fastshap) + PDPs + Figure 2 (matching your original colors/layout)
################################################################################

rm(list = ls())
set.seed(42)

# -----------------------------
# Libraries
# -----------------------------
pkgs <- c(
  "randomForest", "dplyr", "tidyr", "ggplot2", "lubridate",
  "purrr", "tibble", "scales",
  "fastshap", "ggbeeswarm", "patchwork", "pdp"
)
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# -----------------------------
# I/O
# -----------------------------
if (!dir.exists("Figures")) dir.create("Figures")
if (!dir.exists("Results")) dir.create("Results")

df <- read.csv("Data/RF_final_cleaned_outliers_removed_zeros_replaced.csv",
               stringsAsFactors = FALSE)

# -----------------------------
# Dates + Season definition (Growing = Apr–Oct; Non-growing = Nov–Mar)
# -----------------------------
df$Collection_date <- as.Date(df$Collection_date)

df <- df %>%
  mutate(
    month_num = lubridate::month(Collection_date),
    Season = case_when(
      month_num %in% 4:10 ~ "Growing",
      month_num %in% c(11, 12, 1, 2, 3) ~ "Non-Growing",
      TRUE ~ NA_character_
    ),
    Season = factor(Season, levels = c("Growing", "Non-Growing"))
  ) %>%
  filter(!is.na(Season))

# -----------------------------
# Response + predictors (NO Month)
# -----------------------------
response_var <- "Resp_rate_mg_C_m2_h"

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

stopifnot(all(c("Plot", response_var, predictor_vars) %in% names(df)))

# -----------------------------
# RF hyperparameters (match your original)
# -----------------------------
ntree <- 1000
mtry <- 3
nodesize <- 5

# -----------------------------
# Metrics helpers
# -----------------------------
rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
mae  <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)

# SSE-based R2 (can be negative; honest out-of-sample)
r2_sse <- function(obs, pred) {
  ss_res <- sum((obs - pred)^2, na.rm = TRUE)
  ss_tot <- sum((obs - mean(obs, na.rm = TRUE))^2, na.rm = TRUE)
  if (ss_tot == 0) return(NA_real_)
  1 - ss_res / ss_tot
}

# Correlation-based R2 (stable sometimes; can be NA if zero variance)
r2_cor <- function(obs, pred) {
  if (sd(obs, na.rm = TRUE) == 0 || sd(pred, na.rm = TRUE) == 0) return(NA_real_)
  cor(obs, pred, use = "complete.obs")^2
}

# -----------------------------
# 1) Leave-one-plot-out CV (primary validation)
# -----------------------------
cv_fold_metrics <- list()
cv_all_preds <- list()

for (season in levels(df$Season)) {
  
  df_s <- df %>% filter(Season == season) %>%
    drop_na(all_of(c(response_var, predictor_vars, "Plot")))
  
  plots <- sort(unique(df_s$Plot))
  
  for (p in plots) {
    train <- df_s %>% filter(Plot != p)
    test  <- df_s %>% filter(Plot == p)
    
    rf_fit <- randomForest(
      x = train[, predictor_vars],
      y = train[[response_var]],
      ntree = ntree,
      mtry = mtry,
      nodesize = nodesize,
      importance = TRUE
    )
    
    pred <- predict(rf_fit, newdata = test[, predictor_vars])
    
    cv_all_preds[[length(cv_all_preds) + 1]] <- tibble(
      Season = season,
      Plot = as.character(p),
      Collection_date = test$Collection_date,
      Observed = test[[response_var]],
      Predicted = as.numeric(pred)
    )
    
    cv_fold_metrics[[length(cv_fold_metrics) + 1]] <- tibble(
      Season = season,
      Model = "Mechanistic_RF_noMonth",
      group = as.character(p),  # <- consistent type prevents bind_rows errors
      n_test = nrow(test),
      RMSE = rmse(test[[response_var]], pred),
      MAE  = mae(test[[response_var]], pred),
      R2_SSE = r2_sse(test[[response_var]], pred),
      R2_COR = r2_cor(test[[response_var]], pred)
    )
  }
}

cv_results <- bind_rows(cv_fold_metrics)
cv_preds   <- bind_rows(cv_all_preds)

write.csv(cv_results, "Results/blocked_cv_fold_metrics_noMonth.csv", row.names = FALSE)
write.csv(cv_preds,   "Results/blocked_cv_predictions_noMonth.csv", row.names = FALSE)

cv_summary <- cv_results %>%
  group_by(Season, Model) %>%
  summarise(
    mean_RMSE = mean(RMSE, na.rm = TRUE),
    median_RMSE = median(RMSE, na.rm = TRUE),
    mean_MAE  = mean(MAE, na.rm = TRUE),
    median_R2_SSE = median(R2_SSE, na.rm = TRUE),
    median_R2_COR = median(R2_COR, na.rm = TRUE),
    frac_R2_SSE_pos = mean(R2_SSE > 0, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(cv_summary, "Results/blocked_cv_summary_noMonth.csv", row.names = FALSE)
print(cv_summary)

p_cv <- ggplot(cv_preds, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.65) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  facet_wrap(~Season, scales = "free") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Predicted vs Observed (Leave-one-plot-out CV)",
    subtitle = "Mechanistic RF (no Month)",
    x = "Observed respiration (mg C m^-2 h^-1)",
    y = "Predicted respiration (mg C m^-2 h^-1)"
  )
ggsave("Figures/pred_vs_obs_LOPO_noMonth.png", p_cv, width = 10, height = 5, dpi = 300)

# -----------------------------
# 2) Fit final seasonal RF models (for interpretation/SHAP/PDP)
# -----------------------------
rf_models <- list()
Xy_by_season <- list()

for (season in levels(df$Season)) {
  df_s <- df %>% filter(Season == season) %>%
    drop_na(all_of(c(response_var, predictor_vars)))
  
  X <- df_s[, predictor_vars, drop = FALSE]
  y <- df_s[[response_var]]
  
  rf_models[[season]] <- randomForest(
    x = X, y = y,
    ntree = ntree, mtry = mtry, nodesize = nodesize,
    importance = TRUE
  )
  
  Xy_by_season[[season]] <- list(df = df_s, X = X, y = y)
}

# -----------------------------
# 3) Plot-level bootstrap permutation importance (cluster bootstrap)
# -----------------------------
bootstrap_importance_by_plot <- function(df_s, predictors, response, n_boot = 1000) {
  plots <- unique(df_s$Plot)
  imp_mat <- matrix(NA_real_, nrow = n_boot, ncol = length(predictors))
  colnames(imp_mat) <- predictors
  
  for (b in seq_len(n_boot)) {
    boot_plots <- sample(plots, size = length(plots), replace = TRUE)
    boot_df <- purrr::map_dfr(boot_plots, ~df_s %>% filter(Plot == .x)) %>%
      drop_na(all_of(c(response, predictors)))
    
    rf_b <- randomForest(
      x = boot_df[, predictors],
      y = boot_df[[response]],
      ntree = ntree, mtry = mtry, nodesize = nodesize,
      importance = TRUE
    )
    
    imp_mat[b, ] <- importance(rf_b)[, "%IncMSE"]
    if (b %% 100 == 0) cat("bootstrap", b, "/", n_boot, "\n")
  }
  
  imp_mat
}

importance_list <- list()

for (season in levels(df$Season)) {
  df_s <- Xy_by_season[[season]]$df
  imp_mat <- bootstrap_importance_by_plot(df_s, predictor_vars, response_var, n_boot = 1000)
  
  imp_df <- tibble(
    Variable = colnames(imp_mat),
    Mean_IncMSE = colMeans(imp_mat, na.rm = TRUE),
    CI_lower = apply(imp_mat, 2, quantile, probs = 0.025, na.rm = TRUE),
    CI_upper = apply(imp_mat, 2, quantile, probs = 0.975, na.rm = TRUE),
    Season = season
  ) %>% arrange(desc(Mean_IncMSE))
  
  importance_list[[season]] <- imp_df
}

importance_all <- bind_rows(importance_list)
write.csv(importance_all, "Results/importance_bootstrap_noMonth.csv", row.names = FALSE)

top5 <- importance_all %>%
  group_by(Season) %>%
  slice_max(order_by = Mean_IncMSE, n = 5, with_ties = FALSE) %>%
  ungroup()

# -----------------------------
# 4) TRUE SHAP (fastshap)  <-- FIX: NO .parallel argument
# -----------------------------
pretty_labels <- c(
  "Month" = "Month",
  "SoilTemp_mean" = "Soil Temperature",
  "GWC" = "Gravimetric Water\nContent (%)",
  "MBC_ug_gdw_bgc" = "Microbial Biomass C\n(µg/g dry soil)",
  "BG_nmol_gdw_h" = "β-Glucosidase\n(nmol/g/h)",
  "NAG_nmol_gdw_h" = "N-Acetylglucosaminidase\n(nmol/g/h)",
  "Leu_nmol_gdw_h" = "Leucine Aminopeptidase\n(nmol/g/h)",
  "AP_nmol_gdw_h" = "Acid Phosphatase\n(nmol/g/h)",
  "CB_nmol_gdw_h" = "Cellobiohydrolase\n(nmol/g/h)",
  "SETN_ug_gdw" = "Extractable N\n(µg/g dry soil)",
  "Roots_gdw" = "Root Biomass\n(g/g dry soil)",
  "MBP_IP_ug_gdw" = "Microbial Biomass P\n(µg PO4/g dry soil)"
)

pred_fun <- function(object, newdata) predict(object, newdata = newdata)

shap_long_list <- list()

for (season in levels(df$Season)) {
  
  X <- Xy_by_season[[season]]$X
  rf_fit <- rf_models[[season]]
  
  vars_to_explain <- top5 %>% filter(Season == season) %>% pull(Variable)
  
  # fastshap: returns SHAP matrix for ALL predictors we pass in X
  # We'll compute for all predictors, then subset to top vars for plotting
  shap_mat <- fastshap::explain(
    object = rf_fit,
    X = X,
    pred_wrapper = pred_fun,
    nsim = 200,
    adjust = TRUE
    # <-- IMPORTANT: no .parallel argument here
  )
  
  shap_df <- as.data.frame(shap_mat)
  shap_df$row_id <- seq_len(nrow(shap_df))
  
  X2 <- as.data.frame(X)
  X2$row_id <- seq_len(nrow(X2))
  
  shap_long <- shap_df %>%
    select(row_id, all_of(vars_to_explain)) %>%
    pivot_longer(cols = all_of(vars_to_explain), names_to = "Variable", values_to = "SHAP_value") %>%
    left_join(
      X2 %>%
        select(row_id, all_of(vars_to_explain)) %>%
        pivot_longer(cols = all_of(vars_to_explain), names_to = "Variable", values_to = "Feature_value"),
      by = c("row_id", "Variable")
    ) %>%
    mutate(
      Season_label = factor(season, levels = c("Growing", "Non-Growing"),
                            labels = c("Growing Season", "Non-Growing Season")),
      Variable_clean = dplyr::recode(Variable, !!!pretty_labels)
    ) %>%
    group_by(Season_label, Variable) %>%
    mutate(
      Feature_value_norm = (Feature_value - min(Feature_value, na.rm = TRUE)) /
        (max(Feature_value, na.rm = TRUE) - min(Feature_value, na.rm = TRUE) + 1e-10)
    ) %>%
    ungroup()
  
  shap_long_list[[season]] <- shap_long
}

shap_plot_df <- bind_rows(shap_long_list)
write.csv(shap_plot_df, "Results/shap_long_noMonth.csv", row.names = FALSE)

# -----------------------------
# 5) PDPs for top 5 (supplemental)
# -----------------------------
pdp_out <- list()
for (season in levels(df$Season)) {
  rf_fit <- rf_models[[season]]
  X <- Xy_by_season[[season]]$X
  vars_to_plot <- top5 %>% filter(Season == season) %>% pull(Variable)
  
  for (v in vars_to_plot) {
    pd <- pdp::partial(rf_fit, pred.var = v, train = X,
                       grid.resolution = 25, progress = "none")
    pd$Season <- season
    pd$Variable <- v
    pdp_out[[paste(season, v, sep = "_")]] <- pd
  }
}
pdp_df <- bind_rows(pdp_out)
write.csv(pdp_df, "Results/pdp_noMonth.csv", row.names = FALSE)

# Optional PDP plot per season
for (season in levels(df$Season)) {
  pd_s <- pdp_df %>% filter(Season == season)
  
  # pdp::partial uses column name same as variable; reshape safely
  pd_long <- pd_s %>%
    pivot_longer(cols = -c(yhat, Season, Variable), names_to = "xname", values_to = "xval") %>%
    filter(xname == Variable)
  
  p_pdp <- ggplot(pd_long, aes(x = xval, y = yhat)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~Variable, scales = "free_x", ncol = 2) +
    theme_classic(base_size = 11) +
    labs(
      title = paste("Partial dependence plots:", season),
      x = NULL,
      y = "Partial dependence (mean prediction)"
    )
  
  ggsave(paste0("Figures/PDP_", season, "_noMonth.png"),
         p_pdp, width = 10, height = 7, dpi = 600)
}

# -----------------------------
# 6) Figure 2 (match your original script colors/layout)
# -----------------------------
season_colors <- c("Growing Season" = "#E69F00", "Non-Growing Season" = "#56B4E9")

panel_A <- top5 %>%
  mutate(
    Variable_clean = dplyr::recode(Variable, !!!pretty_labels),
    Season_label = factor(Season, levels = c("Growing", "Non-Growing"),
                          labels = c("Growing Season", "Non-Growing Season"))
  ) %>%
  ggplot(aes(x = reorder(Variable_clean, Mean_IncMSE),
             y = Mean_IncMSE,
             fill = Season_label)) +
  geom_col(position = position_dodge(width = 0.75),
           alpha = 0.85, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                position = position_dodge(width = 0.75),
                width = 0.3, linewidth = 0.5) +
  coord_flip() +
  scale_fill_manual(values = season_colors, name = "Season") +
  labs(title = "A", x = NULL, y = "Variable Importance (% Increase in MSE)") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = -0.15),
    axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 10)),
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 15, 10, 10)
  )

panel_B <- shap_plot_df %>%
  ggplot(aes(x = Variable_clean, y = SHAP_value, color = Feature_value_norm)) +
  ggbeeswarm::geom_quasirandom(alpha = 0.6, size = 1.5, width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.5) +
  facet_wrap(~Season_label, ncol = 1, scales = "free_y") +
  coord_flip() +
  scale_color_gradient2(
    low = "#0072B2",
    mid = "#F0E442",
    high = "#D55E00",
    midpoint = 0.5,
    name = "Feature Value",
    labels = c("Low", "", "High"),
    breaks = c(0, 0.5, 1)
  ) +
  labs(title = "B", x = NULL, y = "SHAP value (impact on prediction)") +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = -0.08),
    axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 10)),
    axis.text = element_text(size = 10, color = "black"),
    strip.background = element_rect(fill = "gray95", color = "black", linewidth = 0.5),
    strip.text = element_text(face = "bold", size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    axis.line = element_line(color = "black", linewidth = 0.5),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    plot.margin = margin(10, 15, 10, 10)
  )

figure_2 <- panel_A | panel_B

ggsave("Figures/Figure_2_manuscript_noMonth.png",
       plot = figure_2, width = 15, height = 10, dpi = 600, bg = "white")

# PDF export (vector). If cairo_pdf isn’t available, switch device="pdf"
ggsave("Figures/Figure_2_manuscript_noMonth.pdf",
       plot = figure_2, width = 15, height = 10, device = cairo_pdf, bg = "white")

message("DONE. Check Figures/ and Results/ directories.")
