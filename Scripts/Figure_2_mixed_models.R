################################################################################
# Publication workflow: Respiration drivers by season (Growing vs Non-Growing)
# - Univariate Resp response (monthly repeated measures by Plot)
# - Mixed model: Resp_log ~ biogeochem + time + treatment + (1|Plot)
# - Optional AR(1) within Plot (auto-selected by AIC and Phi)
# - Ranks predictors per season by drop-one ΔAIC + LRT + FDR
# - Fits a combined Season interaction model (optional) to test differences
# - Saves a publication-style figure:
#   A: Resp vs GWC by Season (scatter + loess)
#   B: Growing: top drivers forest plot (standardized effects + 95% CI)
#   C: Non-growing: top drivers forest plot (standardized effects + 95% CI)
################################################################################

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(tidyverse)
  library(nlme)
  library(patchwork)
  library(broom.mixed)
})

# ----------------------------
# 0) Settings
# ----------------------------
infile <- "RF_final_cleaned_outliers_removed_zeros_replaced.csv"
if (!file.exists(infile) && file.exists(file.path("Data", infile))) {
  infile <- file.path("Data", infile)
}
if (!file.exists(infile)) stop("Could not find input file: ", infile)

out_fig_dir <- "Figures"
out_res_dir <- "Results"
dir.create(out_fig_dir, showWarnings = FALSE)
dir.create(out_res_dir, showWarnings = FALSE)

# Your biogeochem predictor set
preds <- c(
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

top_n_to_plot <- 5

trt_cols <- c("100" = "#5AAE61", "25" = "#9970AB")
season_cols <- c("Growing" = "#E69F00", "Non-Growing" = "#56B4E9")

# ----------------------------
# 1) Load & prep data
# ----------------------------
df <- read.csv(infile, stringsAsFactors = FALSE)

resp_col   <- "Resp_rate_mg_C_m2_h"
date_col   <- "Collection_date"
plot_col   <- "Plot"
season_col <- "Season"
trt_col    <- "Irrigation_trmt_pct"

need <- c(date_col, plot_col, season_col, trt_col, resp_col, preds)
miss <- setdiff(need, names(df))
if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

df <- df %>%
  mutate(
    Collection_date = as.Date(.data[[date_col]]),
    Plot            = factor(.data[[plot_col]]),
    Season          = factor(.data[[season_col]], levels = c("Growing", "Non-Growing")),
    Irrigation      = factor(.data[[trt_col]]),
    Resp            = .data[[resp_col]],
    Resp_log        = log1p(Resp),
    time_num        = as.numeric(Collection_date)
  ) %>%
  select(Collection_date, time_num, Season, Plot, Irrigation, Resp, Resp_log, all_of(preds)) %>%
  drop_na() %>%
  arrange(Plot, time_num)

# Center time globally
df <- df %>% mutate(time_c = time_num - mean(time_num))

# Global z-scoring of predictors (stable comparability across seasons)
scale_params <- df %>%
  summarise(across(all_of(preds),
                   list(mu = ~mean(.x), sd = ~sd(.x)),
                   .names = "{.col}__{.fn}"))

zscore_global <- function(x, mu, sd) {
  if (!is.finite(sd) || sd == 0) return(rep(0, length(x)))
  (x - mu) / sd
}

for (v in preds) {
  mu <- scale_params[[paste0(v, "__mu")]]
  sd <- scale_params[[paste0(v, "__sd")]]
  df[[paste0(v, "_z")]] <- zscore_global(df[[v]], mu, sd)
}

preds_z <- paste0(preds, "_z")

dat_grow <- df %>% filter(Season == "Growing") %>% arrange(Plot, time_c)
dat_nong <- df %>% filter(Season == "Non-Growing") %>% arrange(Plot, time_c)

# ----------------------------
# 2) Model fitting utilities (NO update() used)
# ----------------------------
fixed_formula <- as.formula(
  paste("Resp_log ~ Irrigation + time_c +", paste(preds_z, collapse = " + "))
)

fit_lme <- function(data, fixed, use_ar1 = FALSE, method = "ML") {
  if (use_ar1) {
    lme(
      fixed = fixed,
      random = ~ 1 | Plot,
      correlation = corAR1(form = ~ time_c | Plot),
      data = data,
      method = method,
      na.action = na.omit,
      control = lmeControl(msMaxIter = 200, msVerbose = FALSE)
    )
  } else {
    lme(
      fixed = fixed,
      random = ~ 1 | Plot,
      data = data,
      method = method,
      na.action = na.omit,
      control = lmeControl(msMaxIter = 200, msVerbose = FALSE)
    )
  }
}

get_phi <- function(m) {
  if (is.null(m$modelStruct$corStruct)) return(NA_real_)
  tryCatch(as.numeric(coef(m$modelStruct$corStruct, unconstrained = FALSE)), error = function(e) NA_real_)
}

choose_best_model <- function(m_noar, m_ar1) {
  # If AR1 failed, choose no-AR
  if (inherits(m_ar1, "try-error")) return(list(best = m_noar, which = "NO_AR1"))
  
  # If no-AR failed (rare), choose AR1
  if (inherits(m_noar, "try-error")) return(list(best = m_ar1, which = "AR1"))
  
  a0 <- AIC(m_noar)
  a1 <- AIC(m_ar1)
  phi <- get_phi(m_ar1)
  phi0ish <- is.finite(phi) && abs(phi) < 0.05
  
  # Prefer simpler if AR1 not better or phi ~ 0
  if (a0 <= a1 || phi0ish) {
    list(best = m_noar, which = "NO_AR1")
  } else {
    list(best = m_ar1, which = "AR1")
  }
}

# ----------------------------
# 3) Fit season models (ML for ranking; REML for coefficients)
# ----------------------------
m_grow_noar <- try(fit_lme(dat_grow, fixed_formula, use_ar1 = FALSE, method = "ML"), silent = TRUE)
m_grow_ar1  <- try(fit_lme(dat_grow, fixed_formula, use_ar1 = TRUE,  method = "ML"), silent = TRUE)
sel_grow <- choose_best_model(m_grow_noar, m_grow_ar1)
m_grow_ml <- sel_grow$best
m_grow_reml <- try(update(m_grow_ml, method = "REML"), silent = TRUE)  # update() ok here because data is global

m_nong_noar <- try(fit_lme(dat_nong, fixed_formula, use_ar1 = FALSE, method = "ML"), silent = TRUE)
m_nong_ar1  <- try(fit_lme(dat_nong, fixed_formula, use_ar1 = TRUE,  method = "ML"), silent = TRUE)
sel_nong <- choose_best_model(m_nong_noar, m_nong_ar1)
m_nong_ml <- sel_nong$best
m_nong_reml <- try(update(m_nong_ml, method = "REML"), silent = TRUE)

cat("\nChosen model (Growing):", sel_grow$which, " AIC=", AIC(m_grow_ml), " Phi=", get_phi(tryCatch(m_grow_ar1, error=function(e) NULL)), "\n")
cat("Chosen model (Non-Growing):", sel_nong$which, " AIC=", AIC(m_nong_ml), " Phi=", get_phi(tryCatch(m_nong_ar1, error=function(e) NULL)), "\n")

# Export model fit summary
model_fit_summary <- function(m, label) {
  tibble(
    Season = label,
    n = nobs(m),
    n_plots = length(unique(m$data$Plot)),
    AIC = AIC(m),
    BIC = BIC(m),
    logLik = as.numeric(logLik(m)),
    AR1_phi = get_phi(m),
    rand_sd_intercept = as.numeric(VarCorr(m)[1, "StdDev"]),
    resid_sd = as.numeric(VarCorr(m)[2, "StdDev"])
  )
}
fit_summ <- bind_rows(
  model_fit_summary(m_grow_ml, "Growing"),
  model_fit_summary(m_nong_ml, "Non-Growing")
)
write.csv(fit_summ, file.path(out_res_dir, "season_model_fit_summary.csv"), row.names = FALSE)

# ----------------------------
# 4) Predictor ranking: drop-one refits (NO update())
# ----------------------------
drop_one_rank <- function(m_full_ml, data, preds_z, use_ar1) {
  full_aic <- AIC(m_full_ml)
  out <- purrr::map_dfr(preds_z, function(v) {
    fixed_drop <- as.formula(paste("Resp_log ~ Irrigation + time_c +",
                                   paste(setdiff(preds_z, v), collapse = " + ")))
    m_drop <- try(fit_lme(data, fixed_drop, use_ar1 = use_ar1, method = "ML"), silent = TRUE)
    if (inherits(m_drop, "try-error")) {
      return(tibble(
        Variable = gsub("_z$", "", v),
        Variable_z = v,
        deltaAIC = NA_real_,
        LRT = NA_real_,
        p = NA_real_
      ))
    }
    LRT <- 2 * (logLik(m_full_ml) - logLik(m_drop))
    tibble(
      Variable = gsub("_z$", "", v),
      Variable_z = v,
      deltaAIC = AIC(m_drop) - full_aic,
      LRT = as.numeric(LRT),
      p = pchisq(as.numeric(LRT), df = 1, lower.tail = FALSE)
    )
  }) %>%
    mutate(p_fdr = p.adjust(p, method = "fdr")) %>%
    arrange(desc(deltaAIC))
  out
}

use_ar1_grow <- identical(sel_grow$which, "AR1")
use_ar1_nong <- identical(sel_nong$which, "AR1")

imp_grow <- drop_one_rank(m_grow_ml, dat_grow, preds_z, use_ar1 = use_ar1_grow)
imp_nong <- drop_one_rank(m_nong_ml, dat_nong, preds_z, use_ar1 = use_ar1_nong)

write.csv(imp_grow, file.path(out_res_dir, "predictor_importance_Growing.csv"), row.names = FALSE)
write.csv(imp_nong, file.path(out_res_dir, "predictor_importance_NonGrowing.csv"), row.names = FALSE)

cat("\nTop drivers (Growing) by ΔAIC:\n")
print(imp_grow %>% slice_head(n = 10))
cat("\nTop drivers (Non-Growing) by ΔAIC:\n")
print(imp_nong %>% slice_head(n = 10))

# ----------------------------
# 5) Optional combined interaction model (tests if effects differ by season)
# ----------------------------
fixed_int <- as.formula(
  paste("Resp_log ~ Season*(Irrigation + time_c +", paste(preds_z, collapse = " + "), ")")
)

m_all_noar <- try(fit_lme(df, fixed_int, use_ar1 = FALSE, method = "ML"), silent = TRUE)

if (!inherits(m_all_noar, "try-error")) {
  an_all <- anova(m_all_noar)
  write.csv(as.data.frame(an_all), file.path(out_res_dir, "interaction_model_anova.csv"))
  cat("\nInteraction model ANOVA saved: Results/interaction_model_anova.csv\n")
} else {
  warning("Combined interaction model failed; skipping export.")
}

# ----------------------------
# 6) Publication figure panels
# ----------------------------

# Panel A: Resp vs GWC by Season (scatter + loess), colored by irrigation
panel_A <- df %>%
  ggplot(aes(x = GWC, y = Resp)) +
  geom_point(aes(fill = Irrigation), shape = 21, color = "black", size = 2.8, alpha = 0.9) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 0.8, color = "black") +
  facet_wrap(~Season, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = trt_cols, name = "Irrigation (%)") +
  theme_bw(base_size = 12) +
  labs(
    title = "A",
    x = "Gravimetric Water Content (GWC)",
    y = expression(Respiration~(mg~C~m^{-2}~h^{-1}))
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = c(0.02, 0.02),
    legend.justification = c(0, 0),
    legend.key = element_rect(fill = NA),
    legend.title = element_text(face = "bold")
  )

# Forest plot for top drivers from REML model (standardized effects)
coef_forest <- function(m_reml, imp_tbl, season_label, panel_title, top_n = 5) {
  
  top_vars_z <- imp_tbl %>% slice_head(n = top_n) %>% pull(Variable_z)
  
  coef_tbl <- broom.mixed::tidy(m_reml, effects = "fixed", conf.int = TRUE) %>%
    filter(term %in% top_vars_z) %>%
    mutate(
      Variable = gsub("_z$", "", term),
      Variable = factor(Variable, levels = rev(gsub("_z$", "", top_vars_z)))
    )
  
  ann <- paste0(season_label, "\n", "n = ", nobs(m_reml))
  
  ggplot(coef_tbl, aes(x = estimate, y = Variable)) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6) +
    geom_pointrange(aes(xmin = conf.low, xmax = conf.high), linewidth = 0.7) +
    theme_bw(base_size = 12) +
    labs(
      title = panel_title,
      x = "Standardized effect (95% CI)",
      y = NULL
    ) +
    annotate("text", x = Inf, y = Inf, label = ann, hjust = 1.05, vjust = 1.2, size = 3.3) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 10)
    )
}

# If REML refit failed for any reason, fall back to ML for plotting coefficients
if (inherits(m_grow_reml, "try-error")) m_grow_reml <- m_grow_ml
if (inherits(m_nong_reml, "try-error")) m_nong_reml <- m_nong_ml

panel_B <- coef_forest(m_grow_reml, imp_grow, "Growing", "B", top_n = top_n_to_plot)
panel_C <- coef_forest(m_nong_reml, imp_nong, "Non-Growing", "C", top_n = top_n_to_plot)

fig <- panel_A | (panel_B / panel_C) + plot_layout(widths = c(2, 1))

png_file <- file.path(out_fig_dir, "Figure_ABC_RespGWC_TopDrivers_MixedModels.png")
pdf_file <- file.path(out_fig_dir, "Figure_ABC_RespGWC_TopDrivers_MixedModels.pdf")

ggsave(png_file, fig, width = 15, height = 10, dpi = 300, bg = "white")
ggsave(pdf_file, fig, width = 15, height = 10, bg = "white")

cat("\nSaved figure:\n  ", png_file, "\n  ", pdf_file, "\n", sep = "")

# ----------------------------
# 7) Export coefficient tables for reporting
# ----------------------------
coef_grow <- broom.mixed::tidy(m_grow_reml, effects = "fixed", conf.int = TRUE) %>%
  mutate(Season = "Growing")
coef_nong <- broom.mixed::tidy(m_nong_reml, effects = "fixed", conf.int = TRUE) %>%
  mutate(Season = "Non-Growing")

coef_all <- bind_rows(coef_grow, coef_nong) %>%
  mutate(term_clean = case_when(
    term == "(Intercept)" ~ "(Intercept)",
    term == "time_c" ~ "time_c",
    term == "Irrigation100" ~ "Irrigation100",
    TRUE ~ gsub("_z$", "", term)
  ))

write.csv(coef_all, file.path(out_res_dir, "season_model_coefficients_REML.csv"), row.names = FALSE)

cat("\nDone. Outputs written to:\n",
    "  - Results/season_model_fit_summary.csv\n",
    "  - Results/predictor_importance_Growing.csv\n",
    "  - Results/predictor_importance_NonGrowing.csv\n",
    "  - Results/season_model_coefficients_REML.csv\n",
    "  - Figures/Figure_ABC_RespGWC_TopDrivers_MixedModels.(png|pdf)\n",
    sep = "")
################################################################################
# END
################################################################################
