################################################################################
# Permutation-null test for RF permutation importance (cluster-aware)
# - Seasons: Growing (Apr–Oct), Non-growing (Nov–Mar)
# - Primary predictors: mechanistic, NO Month
# - Null: shuffle response WITHIN plot (preserves plot clustering)
# Outputs:
#   Results/null_importance_summary.csv
#   Results/null_importance_raw.csv
#   Figures/null_importance_plot_[Season].png
################################################################################

rm(list = ls())
set.seed(42)

suppressPackageStartupMessages({
  library(randomForest)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(purrr)
})

# -----------------------------
# 1) Load data + season
# -----------------------------
df <- read.csv("Data/RF_final_cleaned_outliers_removed_zeros_replaced.csv",
               stringsAsFactors = FALSE)

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

# -----------------------------
# 2) Variables (NO Month)
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

# RF params (match your manuscript model)
ntree <- 1000
mtry <- 3
nodesize <- 5

if (!dir.exists("Results")) dir.create("Results")
if (!dir.exists("Figures")) dir.create("Figures")

# -----------------------------
# 3) Helper: fit and get %IncMSE
# -----------------------------
get_incMSE <- function(dat, predictors, response) {
  dat <- dat %>% drop_na(all_of(c(response, predictors)))
  rf <- randomForest(
    x = dat[, predictors, drop = FALSE],
    y = dat[[response]],
    ntree = ntree, mtry = mtry, nodesize = nodesize,
    importance = TRUE
  )
  imp <- importance(rf)[, "%IncMSE"]
  as.numeric(imp)
}

# -----------------------------
# 4) Null generator: shuffle y within Plot
# -----------------------------
shuffle_within_plot <- function(dat, response, plot_col = "Plot") {
  dat %>%
    group_by(.data[[plot_col]]) %>%
    mutate("{response}" := sample(.data[[response]], size = n(), replace = FALSE)) %>%
    ungroup()
}

# -----------------------------
# 5) Run observed + null
# -----------------------------
n_perm <- 500   # start with 200-500; increase to 1000 if you want smoother null

all_null <- list()
all_obs  <- list()

for (season in levels(df$Season)) {
  
  message("Season: ", season)
  
  d <- df %>% filter(Season == season) %>%
    drop_na(all_of(c("Plot", response_var, predictor_vars)))
  
  # Observed importance
  obs <- get_incMSE(d, predictor_vars, response_var)
  obs_df <- tibble(
    Season = season,
    Variable = predictor_vars,
    Observed_IncMSE = obs
  )
  all_obs[[season]] <- obs_df
  
  # Null importances
  null_mat <- matrix(NA_real_, nrow = n_perm, ncol = length(predictor_vars))
  colnames(null_mat) <- predictor_vars
  
  for (b in seq_len(n_perm)) {
    d_null <- shuffle_within_plot(d, response_var, plot_col = "Plot")
    null_mat[b, ] <- get_incMSE(d_null, predictor_vars, response_var)
    if (b %% 50 == 0) message("  perm ", b, "/", n_perm)
  }
  
  null_long <- as.data.frame(null_mat) %>%
    mutate(iter = row_number()) %>%
    pivot_longer(cols = all_of(predictor_vars), names_to = "Variable", values_to = "Null_IncMSE") %>%
    mutate(Season = season)
  
  all_null[[season]] <- null_long
}

obs_all  <- bind_rows(all_obs)
null_all <- bind_rows(all_null)

write.csv(null_all, "Results/null_importance_raw.csv", row.names = FALSE)

# -----------------------------
# 6) Summarize: compare observed to null
# -----------------------------
null_summary <- null_all %>%
  group_by(Season, Variable) %>%
  summarise(
    Null_mean = mean(Null_IncMSE, na.rm = TRUE),
    Null_p95  = quantile(Null_IncMSE, 0.95, na.rm = TRUE),
    Null_p99  = quantile(Null_IncMSE, 0.99, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(obs_all, by = c("Season","Variable")) %>%
  mutate(
    Exceeds_p95 = Observed_IncMSE > Null_p95,
    Exceeds_p99 = Observed_IncMSE > Null_p99,
    z_score = (Observed_IncMSE - Null_mean) / sd(null_all$Null_IncMSE, na.rm = TRUE) # rough global scale
  ) %>%
  arrange(Season, desc(Observed_IncMSE))

write.csv(null_summary, "Results/null_importance_summary.csv", row.names = FALSE)

print(null_summary %>% group_by(Season) %>% summarise(
  n_vars = n(),
  n_exceed_p95 = sum(Exceeds_p95, na.rm = TRUE),
  n_exceed_p99 = sum(Exceeds_p99, na.rm = TRUE),
  .groups = "drop"
))

# -----------------------------
# 7) Plot: observed vs null distribution (top variables)
# -----------------------------
for (season in levels(df$Season)) {
  
  top_vars <- obs_all %>%
    filter(Season == season) %>%
    arrange(desc(Observed_IncMSE)) %>%
    slice_head(n = 6) %>%
    pull(Variable)
  
  pdat_null <- null_all %>% filter(Season == season, Variable %in% top_vars)
  pdat_obs  <- obs_all  %>% filter(Season == season, Variable %in% top_vars)
  
  p <- ggplot(pdat_null, aes(x = Variable, y = Null_IncMSE)) +
    geom_violin(trim = TRUE, alpha = 0.35) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.5) +
    geom_point(data = pdat_obs, aes(y = Observed_IncMSE), size = 3) +
    coord_flip() +
    theme_minimal(base_size = 12) +
    labs(
      title = paste("Permutation-null importance test:", season),
      subtitle = "Violin/box = null (%IncMSE after shuffling respiration within plot); dot = observed",
      x = NULL, y = "% Increase in MSE (Permutation Importance)"
    )
  
  ggsave(paste0("Figures/null_importance_plot_", season, ".png"),
         p, width = 10, height = 6, dpi = 300)
}

message("DONE. See Results/null_importance_summary.csv and Figures/null_importance_plot_*.png")
################################################################################
