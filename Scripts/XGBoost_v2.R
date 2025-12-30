################################################################################
# GPBoost mixed-effects ML (Windows-safe: numeric matrices)
# FIXES INCLUDED:
#   1) bagging not supported with GPBoost -> bagging_freq must be 0
#   2) predict() returns a LIST (with $response_mean) when gp_model is used
#      -> extract $response_mean instead of as.numeric(list)
#
# - Separate models by season
# - Random intercepts: Plot + Block
# - CV options: "date" (recommended) or "plot"
# - Held-out permutation importance (ΔRMSE)
# - SHAP-like contributions (tree part only): predcontrib=TRUE
# - Partial dependence plots (top 5 predictors)
################################################################################

rm(list = ls())
set.seed(42)

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(tidyr)
  library(ggplot2)
  library(gpboost)
})

# -----------------------------
# USER SETTINGS
# -----------------------------
DATA_PATH <- "Data/RF_final_cleaned_outliers_removed_zeros_replaced.csv"
USE_LOG1P <- TRUE

# Choose CV type:
#   "date" = leave-one-date-out within season (recommended for mechanism)
#   "plot" = leave-one-plot-out within season (new-plot generalization)
CV_TYPE <- "date"

NROUNDS <- 600
VERBOSE <- 0

# IMPORTANT: GPBoost does not support bagging -> bagging_freq must be 0
GB_PARAMS <- list(
  objective = "regression_l2",
  learning_rate = 0.03,
  num_leaves = 15L,
  max_depth = 5L,
  min_data_in_leaf = 8L,
  feature_fraction = 0.80,
  bagging_freq = 0L,     # REQUIRED
  lambda_l2 = 5.0
)

# permutation importance settings
N_PERM_PER_FOLD <- 1  # raise to 3-5 for smoother estimates

# PDP settings
PDP_POINTS <- 20

if (!dir.exists("Results")) dir.create("Results")
if (!dir.exists("Figures")) dir.create("Figures")

# -----------------------------
# Helpers
# -----------------------------
rmse_vec <- function(truth, pred) sqrt(mean((truth - pred)^2, na.rm = TRUE))

r2_sse <- function(truth, pred) {
  ss_res <- sum((truth - pred)^2, na.rm = TRUE)
  ss_tot <- sum((truth - mean(truth, na.rm = TRUE))^2, na.rm = TRUE)
  if (ss_tot == 0) return(NA_real_)
  1 - ss_res / ss_tot
}

make_folds <- function(dat, cv_type = c("date","plot")) {
  cv_type <- match.arg(cv_type)
  if (cv_type == "date") {
    groups <- sort(unique(dat$Collection_date))
    lapply(groups, function(g) list(test = which(dat$Collection_date == g), id = as.character(g)))
  } else {
    groups <- sort(unique(dat$Plot))
    lapply(groups, function(g) list(test = which(dat$Plot == g), id = as.character(g)))
  }
}

# Ensure numeric matrix with stable column names
to_num_matrix <- function(df, col_order) {
  m <- as.matrix(df[, col_order, drop = FALSE])
  storage.mode(m) <- "double"
  colnames(m) <- col_order
  m
}

# -----------------------------
# Fit + predict
# -----------------------------
fit_gpboost <- function(X_train_mat, y_train, group_train_mat, params, nrounds, verbose = 0) {
  gp_model <- GPModel(
    likelihood = "gaussian",
    group_data = group_train_mat
  )
  bst <- gpboost(
    data = X_train_mat,
    label = y_train,
    gp_model = gp_model,
    params = params,
    nrounds = nrounds,
    verbose = verbose
  )
  bst
}

# ---- CRITICAL FIX HERE ----
# When gp_model is present and predcontrib = FALSE, predict() returns a list
# containing $response_mean, not a numeric vector.
predict_gpboost_response <- function(bst, X_test_mat, group_test_mat) {
  pred <- predict(
    bst,
    data = X_test_mat,
    group_data_pred = group_test_mat,
    predcontrib = FALSE,
    predict_var = FALSE,
    pred_latent = FALSE
  )
  if (is.list(pred)) {
    if (!is.null(pred$response_mean)) return(as.numeric(pred$response_mean))
    # fallback (should rarely be needed)
    if (!is.null(pred$fixed_effect) && !is.null(pred$random_effect_mean)) {
      return(as.numeric(pred$fixed_effect + pred$random_effect_mean))
    }
    stop("predict() returned a list but no $response_mean found. Names: ",
         paste(names(pred), collapse = ", "))
  }
  as.numeric(pred)
}

# For SHAP-like contributions (tree part), predict() returns a matrix
# (and ignores gp_model/random effects when predcontrib=TRUE).
predict_gpboost_contrib <- function(bst, X_test_mat, group_test_mat) {
  pred <- predict(
    bst,
    data = X_test_mat,
    group_data_pred = group_test_mat,
    predcontrib = TRUE,
    predict_var = FALSE,
    pred_latent = FALSE
  )
  # Most often this is a numeric matrix already
  if (is.list(pred)) {
    # Extremely defensive fallback
    if (!is.null(pred$pred_contrib)) return(as.matrix(pred$pred_contrib))
    stop("predcontrib prediction returned a list unexpectedly. Names: ",
         paste(names(pred), collapse = ", "))
  }
  as.matrix(pred)
}

# -----------------------------
# Permutation importance (ΔRMSE on held-out)
# -----------------------------
perm_importance_fold <- function(bst, X_test_mat, y_test, group_test_mat, feature_names, n_perm = 1) {
  base_pred <- predict_gpboost_response(bst, X_test_mat, group_test_mat)
  base_rmse <- rmse_vec(y_test, base_pred)
  
  out <- vector("list", length(feature_names))
  names(out) <- feature_names
  
  for (v in feature_names) {
    deltas <- numeric(n_perm)
    j <- which(colnames(X_test_mat) == v)
    
    for (p in seq_len(n_perm)) {
      Xp <- X_test_mat
      Xp[, j] <- sample(Xp[, j], replace = FALSE)
      perm_pred <- predict_gpboost_response(bst, Xp, group_test_mat)
      deltas[p] <- rmse_vec(y_test, perm_pred) - base_rmse
    }
    
    out[[v]] <- tibble(Variable = v, importance = mean(deltas, na.rm = TRUE))
  }
  
  bind_rows(out) %>% mutate(base_rmse = base_rmse)
}

summarize_importance <- function(imp_long) {
  imp_long %>%
    group_by(Season, Variable) %>%
    summarise(
      mean_importance = mean(importance, na.rm = TRUE),
      median_importance = median(importance, na.rm = TRUE),
      lo = quantile(importance, 0.025, na.rm = TRUE),
      hi = quantile(importance, 0.975, na.rm = TRUE),
      prop_pos = mean(importance > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Season, desc(mean_importance))
}

# -----------------------------
# PDP
# -----------------------------
compute_pdp <- function(bst, X_ref_mat, group_ref_mat, var, grid_vals) {
  j <- which(colnames(X_ref_mat) == var)
  pd <- lapply(grid_vals, function(g) {
    Xg <- X_ref_mat
    Xg[, j] <- g
    pred_g <- predict_gpboost_response(bst, Xg, group_ref_mat)
    tibble(Variable = var, x = g, yhat = mean(pred_g, na.rm = TRUE))
  })
  bind_rows(pd)
}

# -----------------------------
# 1) Load data
# -----------------------------
df <- read.csv(DATA_PATH, stringsAsFactors = FALSE)
df$Collection_date <- as.Date(df$Collection_date)

df <- df %>%
  mutate(
    month_num = month(Collection_date),
    Season = case_when(
      month_num %in% 4:10 ~ "Growing",
      month_num %in% c(11, 12, 1, 2, 3) ~ "Non-Growing",
      TRUE ~ NA_character_
    ),
    Season = factor(Season, levels = c("Growing", "Non-Growing"))
  ) %>%
  filter(!is.na(Season))

if (!"Plot" %in% names(df)) stop("Missing Plot column")
df$Plot <- as.character(df$Plot)

if (!"Block" %in% names(df)) df$Block <- "B1"
df$Block <- as.character(df$Block)

response_raw <- "Resp_rate_mg_C_m2_h"
predictors <- c(
  "GWC","SETN_ug_gdw","MBC_ug_gdw_bgc","MBP_IP_ug_gdw",
  "BG_nmol_gdw_h","CB_nmol_gdw_h","NAG_nmol_gdw_h","Leu_nmol_gdw_h",
  "AP_nmol_gdw_h","Roots_gdw","SoilTemp_mean"
)

need_cols <- c("Season","Plot","Block","Collection_date", response_raw, predictors)
miss <- setdiff(need_cols, names(df))
if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

df <- df %>%
  mutate(y = if (USE_LOG1P) log1p(.data[[response_raw]]) else .data[[response_raw]])

# -----------------------------
# 2) Run per season
# -----------------------------
all_preds <- list()
all_metrics <- list()
all_imp_long <- list()
all_shap_long <- list()
all_pdp <- list()

for (season in levels(df$Season)) {
  
  message("\n==============================")
  message("Season: ", season, " | CV_TYPE: ", CV_TYPE)
  message("==============================")
  
  dat <- df %>%
    filter(Season == season) %>%
    select(Season, Collection_date, Plot, Block, y, all_of(predictors)) %>%
    drop_na()
  
  # Build X as numeric matrix with fixed column order
  X_df <- dat[, predictors, drop = FALSE]
  X_mat <- to_num_matrix(X_df, predictors)
  y_vec <- dat$y
  group_mat <- as.matrix(dat[, c("Plot","Block")])
  
  folds <- make_folds(dat, cv_type = CV_TYPE)
  pred_vec <- rep(NA_real_, nrow(dat))
  
  imp_fold_list <- list()
  shap_fold_list <- list()
  
  for (k in seq_along(folds)) {
    test_idx <- folds[[k]]$test
    train_idx <- setdiff(seq_len(nrow(dat)), test_idx)
    
    X_train <- X_mat[train_idx, , drop = FALSE]
    y_train <- y_vec[train_idx]
    g_train <- group_mat[train_idx, , drop = FALSE]
    
    X_test <- X_mat[test_idx, , drop = FALSE]
    y_test <- y_vec[test_idx]
    g_test <- group_mat[test_idx, , drop = FALSE]
    
    bst <- fit_gpboost(X_train, y_train, g_train, GB_PARAMS, NROUNDS, verbose = VERBOSE)
    
    pred_test <- predict_gpboost_response(bst, X_test, g_test)
    pred_vec[test_idx] <- pred_test
    
    # permutation importance (held-out)
    imp_k <- perm_importance_fold(bst, X_test, y_test, g_test, predictors, n_perm = N_PERM_PER_FOLD) %>%
      mutate(Season = season, fold = k, fold_id = folds[[k]]$id)
    imp_fold_list[[k]] <- imp_k
    
    # SHAP-like contributions (tree part ONLY)
    shap_mat <- predict_gpboost_contrib(bst, X_test, g_test)
    shap_df <- as.data.frame(shap_mat)
    shap_df$.row <- test_idx
    
    shap_long <- shap_df %>%
      pivot_longer(cols = -c(.row), names_to = "Feature", values_to = "Contribution") %>%
      mutate(Season = season, fold = k, fold_id = folds[[k]]$id)
    shap_fold_list[[k]] <- shap_long
  }
  
  preds <- dat %>% mutate(.pred = pred_vec)
  all_preds[[season]] <- preds
  all_metrics[[season]] <- tibble(
    Season = season,
    CV_TYPE = CV_TYPE,
    n = nrow(preds),
    n_plots = length(unique(preds$Plot)),
    n_dates = length(unique(preds$Collection_date)),
    pooled_RMSE = rmse_vec(preds$y, preds$.pred),
    pooled_R2_SSE = r2_sse(preds$y, preds$.pred)
  )
  
  imp_long <- bind_rows(imp_fold_list)
  all_imp_long[[season]] <- imp_long
  all_shap_long[[season]] <- bind_rows(shap_fold_list)
  
  # PDP from model trained on all data (season)
  bst_full <- fit_gpboost(X_mat, y_vec, group_mat, GB_PARAMS, NROUNDS, verbose = VERBOSE)
  
  imp_rank <- imp_long %>%
    group_by(Variable) %>%
    summarise(mean_importance = mean(importance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_importance))
  
  top_vars <- head(imp_rank$Variable, 5)
  pdp_list <- list()
  
  for (v in top_vars) {
    grid_vals <- as.numeric(quantile(X_mat[, v], probs = seq(0.05, 0.95, length.out = PDP_POINTS), na.rm = TRUE))
    pdp_list[[v]] <- compute_pdp(bst_full, X_mat, group_mat, v, grid_vals) %>%
      mutate(Season = season)
  }
  all_pdp[[season]] <- bind_rows(pdp_list)
}

preds_all <- bind_rows(all_preds)
metrics_all <- bind_rows(all_metrics)
imp_long_all <- bind_rows(all_imp_long)
imp_sum_all <- summarize_importance(imp_long_all)
shap_long_all <- bind_rows(all_shap_long)
pdp_all <- bind_rows(all_pdp)

# -----------------------------
# Save outputs
# -----------------------------
write.csv(metrics_all, "Results/gpboost_cv_metrics_by_season.csv", row.names = FALSE)
write.csv(preds_all, "Results/gpboost_cv_predictions_by_season.csv", row.names = FALSE)

write.csv(imp_long_all, "Results/gpboost_perm_importance_long.csv", row.names = FALSE)
write.csv(imp_sum_all, "Results/gpboost_perm_importance_summary.csv", row.names = FALSE)

write.csv(shap_long_all, "Results/gpboost_shap_predcontrib_long.csv", row.names = FALSE)
write.csv(pdp_all, "Results/gpboost_pdp_by_season.csv", row.names = FALSE)

# -----------------------------
# Quick plots
# -----------------------------
p_pred <- ggplot(preds_all, aes(x = y, y = .pred)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  facet_wrap(~Season, scales = "free") +
  theme_minimal(base_size = 12) +
  labs(
    title = paste0("GPBoost CV predictions by season (CV_TYPE=", CV_TYPE, ")"),
    subtitle = if (USE_LOG1P) "Outcome = log1p(respiration)" else "Outcome = respiration (raw)",
    x = "Observed y",
    y = "Predicted y"
  )
ggsave("Figures/gpboost_pred_vs_obs_by_season.png", p_pred, width = 10, height = 5, dpi = 300)

imp_plot <- imp_sum_all %>%
  group_by(Season) %>%
  slice_max(order_by = mean_importance, n = 10, with_ties = FALSE) %>%
  ungroup()

p_imp <- ggplot(imp_plot, aes(x = reorder(Variable, mean_importance), y = mean_importance)) +
  geom_col() +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
  coord_flip() +
  facet_wrap(~Season, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(
    title = "GPBoost held-out permutation importance (ΔRMSE)",
    subtitle = "Error bars = 95% CI across CV folds",
    x = NULL,
    y = "Mean ΔRMSE when permuted"
  )
ggsave("Figures/gpboost_perm_importance_top10.png", p_imp, width = 10, height = 6, dpi = 300)

p_pdp <- ggplot(pdp_all, aes(x = x, y = yhat)) +
  geom_line() +
  facet_grid(Season ~ Variable, scales = "free_x") +
  theme_minimal(base_size = 11) +
  labs(
    title = "Partial dependence (GPBoost)",
    subtitle = "Averaged predictions from model trained on all season data",
    x = "Predictor value",
    y = "Mean predicted y"
  )
ggsave("Figures/gpboost_pdp_top5.png", p_pdp, width = 14, height = 6, dpi = 300)

message("\nDONE. Outputs in Results/ and Figures/.")
################################################################################
