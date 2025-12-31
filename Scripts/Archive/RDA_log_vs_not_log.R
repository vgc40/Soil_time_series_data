################################################################################
# pRDA Sensitivity Grid (FIXED):
#   Resp: raw vs log1p
#   Predictors: raw vs log1p(skewed-only, auto)
#
# partial RDA: X ~ resp + Condition(Plot)
# blocked permutations by Plot
#
# Runs: All, Growing, Non-Growing
# Exports:
#   Results/pRDA_grid_model_stats.csv
#   Results/pRDA_grid_top5_arrows_*.csv
#   Results/pRDA_grid_transform_table_*.csv
################################################################################

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(permute)
})

# ----------------------------
# 1) Load data
# ----------------------------
df <- read.csv("Data/RF_final_cleaned_outliers_removed_zeros_replaced.csv",
               stringsAsFactors = FALSE)

df$Collection_date <- as.Date(df$Collection_date)
df$Plot <- as.factor(df$Plot)

df <- df %>%
  mutate(
    Season = factor(Season, levels = c("Growing", "Non-Growing")),
    Irrigation_trmt_pct = as.factor(Irrigation_trmt_pct),
    Resp = Resp_rate_mg_C_m2_h,
    Resp_log = log1p(Resp_rate_mg_C_m2_h)
  )

# ----------------------------
# 2) Predictor subset (11 vars)
# ----------------------------
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

need <- c("Season","Plot","Irrigation_trmt_pct","Resp","Resp_log", predictor_vars)
miss <- setdiff(need, names(df))
if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

dir.create("Results", showWarnings = FALSE)

# ----------------------------
# 3) Blocked permutations within Plot
# ----------------------------
perm_ctrl <- function(plot_factor, nperm = 999) {
  how(nperm = nperm, blocks = plot_factor)
}

# ----------------------------
# 4) Skewness helper (no extra packages)
# ----------------------------
skewness_fp <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 3) return(NA_real_)
  m <- mean(x)
  s <- sd(x)
  if (s == 0) return(0)
  sum(((x - m) / s)^3) / n
}

# ----------------------------
# 5) Predictor transform: raw vs log1p(skewed-only)
# ----------------------------
transform_predictors <- function(X, method = c("raw", "log1p_auto"),
                                 skew_thresh = 1.0) {
  method <- match.arg(method)
  
  X2 <- as.data.frame(X)
  tab <- tibble(
    Variable = names(X2),
    min = map_dbl(X2, ~ suppressWarnings(min(.x, na.rm = TRUE))),
    skew = map_dbl(X2, ~ skewness_fp(.x)),
    transform = "raw"
  )
  
  if (method == "log1p_auto") {
    do_log <- with(tab, (min >= 0) & is.finite(skew) & (skew > skew_thresh))
    tab$transform[do_log] <- "log1p"
    for (v in tab$Variable[do_log]) {
      X2[[v]] <- log1p(X2[[v]])
    }
  }
  
  list(X_trans = X2, transform_table = tab)
}

# ----------------------------
# 6) Prep data (drop NA, transform X, then scale)
# ----------------------------
prep_data <- function(dat, vars, resp_col, pred_transform) {
  dat2 <- dat %>%
    select(Season, Plot, Irrigation_trmt_pct, all_of(resp_col), all_of(vars)) %>%
    drop_na()
  
  X  <- dat2 %>% select(all_of(vars))
  tx <- transform_predictors(X, method = pred_transform, skew_thresh = 1.0)
  
  X_trans <- tx$X_trans
  Xs <- scale(X_trans)
  
  list(dat = dat2, X_raw = X, X_trans = X_trans, Xs = Xs, transform_table = tx$transform_table)
}

# ----------------------------
# 7) Summarize pRDA model (effect size + blocked p-values)
# ----------------------------
summarize_prda <- function(mod, dat_meta, dataset_label, resp_label, pred_label, nperm = 999) {
  ctrl <- perm_ctrl(dat_meta$Plot, nperm)
  
  a_overall <- anova.cca(mod, permutations = ctrl)
  a_term    <- anova.cca(mod, by = "term", permutations = ctrl)
  
  adj_r2 <- RsquareAdj(mod)$adj.r.squared
  frac_constrained <- mod$CCA$tot.chi / mod$tot.chi
  
  eig <- mod$CCA$eig
  rda1_frac_of_constrained <- ifelse(length(eig) > 0, eig[1] / sum(eig), NA_real_)
  
  tibble(
    Dataset = dataset_label,
    Resp_scale = resp_label,
    Pred_transform = pred_label,
    n = nrow(dat_meta),
    n_plots = n_distinct(dat_meta$Plot),
    Adj_R2 = adj_r2,
    Constrained_fraction = frac_constrained,
    RDA1_fraction_of_constrained = rda1_frac_of_constrained,
    F_overall = a_overall$F[1],
    p_overall = a_overall$`Pr(>F)`[1],
    F_resp = a_term$F[1],
    p_resp = a_term$`Pr(>F)`[1]
  )
}

# ----------------------------
# 8) Arrow table (top N by length) + envfit significance (blocked)
# ----------------------------
arrow_table <- function(mod, Xs, dat_meta, top_n = 5, nperm = 999) {
  ctrl <- perm_ctrl(dat_meta$Plot, nperm)
  
  sites <- scores(mod, display = "sites", scaling = 2)
  if (!("RDA1" %in% colnames(sites))) stop("No RDA1 in site scores (unexpected).")
  rda1 <- sites[, "RDA1"]
  
  sp <- as.data.frame(scores(mod, display = "species", scaling = 2))
  sp$Variable <- rownames(sp)
  
  # common single-constraint view: RDA1 vs PC1
  if (!("RDA1" %in% names(sp))) sp$RDA1 <- 0
  if (!("PC1"  %in% names(sp))) sp$PC1  <- 0
  
  sp <- sp %>%
    mutate(
      arrow_len = sqrt(RDA1^2 + PC1^2),
      angle_deg = atan2(PC1, RDA1) * 180 / pi,
      cor_with_RDA1 = map_dbl(Variable, ~ suppressWarnings(cor(Xs[, .x], rda1))),
      aligned_with_pos_RDA1 = RDA1 > 0
    ) %>%
    arrange(desc(arrow_len)) %>%
    slice_head(n = top_n)
  
  ef <- envfit(mod, Xs, permutations = ctrl)
  
  ef_tbl <- tibble(
    Variable = names(ef$vectors$pvals),
    envfit_r2 = as.numeric(ef$vectors$r^2),
    envfit_p  = as.numeric(ef$vectors$pvals)
  ) %>%
    mutate(envfit_p_fdr = p.adjust(envfit_p, method = "fdr"))
  
  sp %>% left_join(ef_tbl, by = "Variable")
}

# ----------------------------
# 9) One run (dataset × resp_scale × pred_transform)  **FIXED**
# ----------------------------
run_one <- function(dat, dataset_label, resp_col, resp_label,
                    pred_transform, pred_label, vars,
                    top_n = 5, nperm = 999) {
  
  pre <- prep_data(dat, vars, resp_col = resp_col, pred_transform = pred_transform)
  dat2 <- pre$dat
  
  # FIX: create a real column used by vegan formula
  dat2$resp_use <- dat2[[resp_col]]
  
  # partial RDA
  mod <- rda(pre$Xs ~ resp_use + Condition(Plot), data = dat2)
  
  stats  <- summarize_prda(mod, dat2, dataset_label, resp_label, pred_label, nperm = nperm)
  arrows <- arrow_table(mod, pre$Xs, dat2, top_n = top_n, nperm = nperm)
  
  list(model = mod, stats = stats, arrows = arrows, transform_table = pre$transform_table)
}

# ----------------------------
# 10) Run grid
# ----------------------------
datasets <- list(
  All = df,
  Growing = df %>% filter(Season == "Growing"),
  `Non-Growing` = df %>% filter(Season == "Non-Growing")
)

resp_scales <- tribble(
  ~resp_col,  ~resp_label,
  "Resp",     "raw",
  "Resp_log", "log1p"
)

pred_transforms <- tribble(
  ~pred_transform, ~pred_label,
  "raw",           "rawX",
  "log1p_auto",    "log1pX_auto"
)

all_stats <- list()
all_arrows <- list()

for (dname in names(datasets)) {
  for (i in seq_len(nrow(resp_scales))) {
    for (j in seq_len(nrow(pred_transforms))) {
      
      resp_col <- resp_scales$resp_col[i]
      resp_lab <- resp_scales$resp_label[i]
      
      pred_tr  <- pred_transforms$pred_transform[j]
      pred_lab <- pred_transforms$pred_label[j]
      
      cat("\n====================================\n")
      cat("Dataset:", dname, "| Resp:", resp_lab, "| X:", pred_lab, "\n")
      cat("====================================\n")
      
      out <- run_one(
        dat = datasets[[dname]],
        dataset_label = dname,
        resp_col = resp_col,
        resp_label = resp_lab,
        pred_transform = pred_tr,
        pred_label = pred_lab,
        vars = predictor_vars,
        top_n = 5,
        nperm = 999
      )
      
      key <- paste(dname, resp_lab, pred_lab, sep = "__")
      all_stats[[key]] <- out$stats
      
      arrows_out <- out$arrows %>%
        mutate(Dataset = dname, Resp_scale = resp_lab, Pred_transform = pred_lab) %>%
        select(Dataset, Resp_scale, Pred_transform,
               Variable, arrow_len, aligned_with_pos_RDA1, angle_deg, cor_with_RDA1,
               envfit_r2, envfit_p, envfit_p_fdr)
      
      all_arrows[[key]] <- arrows_out
      
      write.csv(arrows_out,
                paste0("Results/pRDA_grid_top5_arrows_", gsub("\\s+","_",dname),
                       "_", resp_lab, "_", pred_lab, ".csv"),
                row.names = FALSE)
      
      tt <- out$transform_table %>%
        mutate(Dataset = dname, Resp_scale = resp_lab, Pred_transform = pred_lab)
      
      write.csv(tt,
                paste0("Results/pRDA_grid_transform_table_", gsub("\\s+","_",dname),
                       "_", resp_lab, "_", pred_lab, ".csv"),
                row.names = FALSE)
    }
  }
}

model_stats <- bind_rows(all_stats)
write.csv(model_stats, "Results/pRDA_grid_model_stats.csv", row.names = FALSE)

arrows_combined <- bind_rows(all_arrows)
write.csv(arrows_combined, "Results/pRDA_grid_top5_arrows_COMBINED.csv", row.names = FALSE)

cat("\nSaved:\n",
    "  Results/pRDA_grid_model_stats.csv\n",
    "  Results/pRDA_grid_top5_arrows_*_raw|log1p_rawX|log1pX_auto.csv\n",
    "  Results/pRDA_grid_top5_arrows_COMBINED.csv\n",
    "  Results/pRDA_grid_transform_table_*.csv\n")
################################################################################
