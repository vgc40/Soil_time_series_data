################################################################################
# XGBoost ML (separate models by season) with Plot-blocked CV (LOPO-style)
# + held-out permutation importance (RMSE increase)
#
# This version fixes your recipe errors by using recipe(y ~ ., data = dat)
# so Plot is present and can be assigned an "id" role.
#
# Outputs:
#   Results/xgb_best_params_by_season.csv
#   Results/xgb_cv_summary_by_season.csv
#   Results/xgb_cv_predictions_by_season.csv
#   Results/xgb_perm_importance_cv_summary.csv
#   Figures/xgb_pred_vs_obs_cv_by_season.png
################################################################################

rm(list = ls())
set.seed(42)

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(tidyr)
  library(ggplot2)
  library(tidymodels)
  library(xgboost)
})

if (!dir.exists("Results")) dir.create("Results")
if (!dir.exists("Figures")) dir.create("Figures")

# -----------------------------
# 1) Load data
# -----------------------------
df <- read.csv("Data/RF_final_cleaned_outliers_removed_zeros_replaced.csv",
               stringsAsFactors = FALSE)

df$Collection_date <- as.Date(df$Collection_date)

# -----------------------------
# 2) Define Season (Growing = Apr–Oct; Non-growing = Nov–Mar)
# -----------------------------
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
# 3) Standardize Plot column name (assumes you have Plot already)
# -----------------------------
if (!"Plot" %in% names(df)) {
  stop("Column 'Plot' not found in df. Your columns are:\n", paste(names(df), collapse = ", "))
}
df$Plot <- factor(df$Plot)

# -----------------------------
# 4) Outcome + predictors (NO Month)
# -----------------------------
response_raw <- "Resp_rate_mg_C_m2_h"

predictors <- c(
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

need_cols <- c(response_raw, predictors, "Plot", "Season")
missing_cols <- setdiff(need_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns:\n", paste(missing_cols, collapse = ", "))
}

# Use log1p to reduce skew/heteroskedasticity
USE_LOG1P <- TRUE
df <- df %>%
  mutate(y = if (USE_LOG1P) log1p(.data[[response_raw]]) else .data[[response_raw]])

# Honest R2_SSE (can be negative; use CV predictions)
r2_sse <- function(truth, pred) {
  ss_res <- sum((truth - pred)^2, na.rm = TRUE)
  ss_tot <- sum((truth - mean(truth, na.rm = TRUE))^2, na.rm = TRUE)
  if (ss_tot == 0) return(NA_real_)
  1 - ss_res / ss_tot
}

# -----------------------------
# 5) XGBoost model (tuned)
# -----------------------------
xgb_spec <-
  boost_tree(
    trees = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    loss_reduction = tune(),
    min_n = tune(),
    sample_size = tune(),
    mtry = tune()
  ) %>%
  set_engine("xgboost") %>%
  set_mode("regression")

# -----------------------------
# 6) Recipe + Plot-blocked CV folds
# -----------------------------
make_recipe <- function(dat) {
  # NOTE: y ~ . ensures Plot is part of the recipe data, so update_role works
  recipe(y ~ ., data = dat) %>%
    update_role(Plot, new_role = "id") %>%     # Plot used only for grouping, not prediction
    step_impute_median(all_numeric_predictors()) %>%
    step_zv(all_predictors())                  # remove any zero-variance predictors in folds
}

make_folds <- function(dat) {
  n_plots <- length(unique(dat$Plot))
  if (n_plots < 3) stop("Need at least 3 plots for grouped CV; found: ", n_plots)
  group_vfold_cv(dat, group = Plot, v = n_plots)  # LOPO-style
}

# Tuning grid
grid <- grid_latin_hypercube(
  trees(range = c(300L, 1500L)),
  tree_depth(range = c(2L, 8L)),
  learn_rate(range = c(-4, -1)),     # log10
  loss_reduction(range = c(-2, 2)),  # log10
  min_n(range = c(2L, 25L)),
  sample_prop(range = c(0.6, 1.0)),
  finalize(mtry(), tibble(mtry = length(predictors))),
  size = 30
)

metric_set_basic <- metric_set(rmse, mae)

# -----------------------------
# 7) Held-out permutation importance on assessment folds (RMSE increase)
# -----------------------------
perm_importance_cv <- function(wf_best, folds, predictors) {
  out <- list()
  
  for (i in seq_along(folds$splits)) {
    sp <- folds$splits[[i]]
    train <- analysis(sp)
    test  <- assessment(sp)
    
    fitted <- fit(wf_best, data = train)
    
    base_pred <- predict(fitted, new_data = test) %>% bind_cols(test %>% select(y))
    base_rmse <- rmse(base_pred, truth = y, estimate = .pred)$.estimate
    
    for (v in predictors) {
      test_perm <- test
      test_perm[[v]] <- sample(test_perm[[v]], replace = FALSE)
      
      perm_pred <- predict(fitted, new_data = test_perm) %>% bind_cols(test_perm %>% select(y))
      perm_rmse <- rmse(perm_pred, truth = y, estimate = .pred)$.estimate
      
      out[[length(out) + 1]] <- tibble(
        fold = i,
        Variable = v,
        base_rmse = base_rmse,
        perm_rmse = perm_rmse,
        importance = perm_rmse - base_rmse
      )
    }
  }
  
  bind_rows(out)
}

# -----------------------------
# 8) Run models separately by season
# -----------------------------
all_best <- list()
all_cv_summary <- list()
all_preds <- list()
all_permimp <- list()

for (season in levels(df$Season)) {
  
  message("\n==============================")
  message("Season: ", season)
  message("==============================")
  
  dat <- df %>%
    filter(Season == season) %>%
    select(y, Plot, all_of(predictors)) %>%
    drop_na()
  
  # Sanity check
  message("Rows: ", nrow(dat), " | Plots: ", length(unique(dat$Plot)))
  
  folds <- make_folds(dat)
  rec <- make_recipe(dat)
  
  # optional check: recipe preps successfully
  invisible(prep(rec, training = dat))
  
  wf <- workflow() %>% add_recipe(rec) %>% add_model(xgb_spec)
  
  # Tune
  tune_res <- tune_grid(
    wf,
    resamples = folds,
    grid = grid,
    metrics = metric_set_basic,
    control = control_grid(save_pred = TRUE, verbose = FALSE)
  )
  
  # If tuning fails, print notes and stop
  if (inherits(tune_res, "try-error")) {
    stop("tune_grid failed.")
  }
  
  best_params <- select_best(tune_res, metric = "rmse")
  all_best[[season]] <- best_params %>% mutate(Season = season)
  
  wf_best <- finalize_workflow(wf, best_params)
  
  # Best-model CV predictions (clean)
  fit_res <- fit_resamples(
    wf_best,
    resamples = folds,
    metrics = metric_set_basic,
    control = control_resamples(save_pred = TRUE)
  )
  
  preds <- collect_predictions(fit_res) %>%
    mutate(Season = season) %>%
    left_join(dat %>% mutate(.row = row_number()) %>% select(.row, Plot), by = ".row")
  
  pooled_r2 <- r2_sse(preds$y, preds$.pred)
  
  cv_metrics <- collect_metrics(fit_res) %>%
    mutate(Season = season)
  
  cv_summary <- cv_metrics %>%
    select(Season, .metric, mean, std_err) %>%
    pivot_wider(names_from = .metric, values_from = c(mean, std_err)) %>%
    mutate(pooled_R2_SSE = pooled_r2)
  
  all_cv_summary[[season]] <- cv_summary
  all_preds[[season]] <- preds
  
  # Held-out permutation importance
  message("Computing held-out permutation importance...")
  pim <- perm_importance_cv(wf_best, folds, predictors) %>%
    mutate(Season = season)
  
  all_permimp[[season]] <- pim
}

best_all <- bind_rows(all_best)
cv_summary_all <- bind_rows(all_cv_summary)
preds_all <- bind_rows(all_preds)
permimp_all <- bind_rows(all_permimp)

# Save outputs
write.csv(best_all, "Results/xgb_best_params_by_season.csv", row.names = FALSE)
write.csv(cv_summary_all, "Results/xgb_cv_summary_by_season.csv", row.names = FALSE)
write.csv(preds_all, "Results/xgb_cv_predictions_by_season.csv", row.names = FALSE)

perm_summary <- permimp_all %>%
  group_by(Season, Variable) %>%
  summarise(
    mean_importance = mean(importance, na.rm = TRUE),
    median_importance = median(importance, na.rm = TRUE),
    lo = quantile(importance, 0.025, na.rm = TRUE),
    hi = quantile(importance, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Season, desc(mean_importance))

write.csv(perm_summary, "Results/xgb_perm_importance_cv_summary.csv", row.names = FALSE)

# Plot predicted vs observed
p <- ggplot(preds_all, aes(x = y, y = .pred)) +
  geom_point(alpha = 0.55) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  facet_wrap(~Season, scales = "free") +
  theme_minimal(base_size = 12) +
  labs(
    title = "XGBoost (best tuned) CV predictions by season (Plot-blocked)",
    subtitle = if (USE_LOG1P) "Outcome = log1p(respiration)" else "Outcome = respiration (raw)",
    x = "Observed",
    y = "Predicted"
  )

ggsave("Figures/xgb_pred_vs_obs_cv_by_season.png", p, width = 10, height = 5, dpi = 300)

message("\nDONE. Check Results/ and Figures/.")
################################################################################
