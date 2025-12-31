################################################################################
# Partial RDA (subset vars) with AUTO log1p transform for skewed predictors
# + blocked permutation tests + figure: Panel A Resp~GWC + Panels B/C pRDA
################################################################################

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(permute)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

# ----------------------------
# 1) Load + prep
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

season_cols <- c("Growing" = "#E69F00", "Non-Growing" = "#56B4E9")
trt_cols    <- c("100" = "#5AAE61", "25" = "#9970AB")  # change if desired

# ----------------------------
# 2) Predictor subset
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

need <- c("Season","Plot","Irrigation_trmt_pct","Resp","Resp_log","GWC", predictor_vars)
miss <- setdiff(need, names(df))
if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

# ----------------------------
# 3) Helpers
# ----------------------------

# simple skewness (no extra packages)
skewness1 <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) return(NA_real_)
  m <- mean(x)
  s <- sd(x)
  if (s == 0) return(0)
  mean(((x - m) / s)^3)
}

# Blocked permutations within Plot
perm_ctrl <- function(plot_factor, nperm = 999) {
  how(nperm = nperm, blocks = plot_factor)
}

# Decide which predictors to log1p (skewed + nonnegative)
choose_log_vars <- function(X, skew_thresh = 1.0) {
  vars <- colnames(X)
  out <- tibble(
    Variable = vars,
    min_val = map_dbl(vars, ~ suppressWarnings(min(X[[.x]], na.rm = TRUE))),
    skew_pos = map_dbl(vars, ~ {
      v <- X[[.x]]
      v <- v[is.finite(v) & v > 0]
      if (length(v) < 5) return(NA_real_)
      skewness1(v)
    })
  ) %>%
    mutate(
      can_log = min_val >= 0,
      do_log1p = can_log & is.finite(skew_pos) & (skew_pos > skew_thresh)
    )
  out
}

# Prep data: drop NA, auto log1p skewed vars, scale X
prep_data <- function(dat, vars, skew_thresh = 1.0, verbose = TRUE) {
  dat2 <- dat %>%
    select(Season, Plot, Irrigation_trmt_pct, Resp, Resp_log, all_of(vars)) %>%
    drop_na()
  
  X <- dat2 %>% select(all_of(vars))
  
  log_tbl <- choose_log_vars(X, skew_thresh = skew_thresh)
  log_vars <- log_tbl %>% filter(do_log1p) %>% pull(Variable)
  
  if (verbose) {
    message("\nAuto log1p transform (skew_thresh = ", skew_thresh, "):")
    if (length(log_vars) == 0) message("  (none)")
    if (length(log_vars) > 0) message("  ", paste(log_vars, collapse = ", "))
  }
  
  # apply log1p to selected vars
  if (length(log_vars) > 0) {
    X <- X %>% mutate(across(all_of(log_vars), ~ log1p(pmax(.x, 0))))
  }
  
  Xs <- scale(X)
  
  list(dat = dat2, X = X, Xs = Xs, log_table = log_tbl, log_vars = log_vars)
}

fit_prda <- function(pre) {
  rda(pre$Xs ~ Resp_log + Condition(Plot), data = pre$dat)
}

summarize_prda <- function(mod, dat_meta, label, nperm = 999) {
  ctrl <- perm_ctrl(dat_meta$Plot, nperm)
  a_overall <- anova.cca(mod, permutations = ctrl)
  a_term    <- anova.cca(mod, by = "term", permutations = ctrl)
  
  tibble(
    Dataset = label,
    n = nrow(dat_meta),
    n_plots = n_distinct(dat_meta$Plot),
    Adj_R2 = RsquareAdj(mod)$adj.r.squared,
    Constrained_fraction = mod$CCA$tot.chi / mod$tot.chi,
    F_overall = unname(a_overall$F[1]),
    p_overall = unname(a_overall$`Pr(>F)`[1]),
    F_resp = unname(a_term$F[1]),
    p_resp = unname(a_term$`Pr(>F)`[1])
  )
}

get_axes <- function(site_scores) {
  cn <- colnames(site_scores)
  if (all(c("RDA1","RDA2") %in% cn)) return(c("RDA1","RDA2"))
  if (all(c("RDA1","PC1")  %in% cn)) return(c("RDA1","PC1"))
  cn[1:2]
}

arrow_table <- function(mod, Xs, dat_meta, top_n = 5, nperm = 999) {
  ctrl <- perm_ctrl(dat_meta$Plot, nperm)
  
  sites <- scores(mod, display = "sites", scaling = 2)
  rda1 <- sites[, "RDA1"]
  
  sp <- as.data.frame(scores(mod, display = "species", scaling = 2))
  sp$Variable <- rownames(sp)
  
  # guarantee coordinates exist
  if (!("RDA1" %in% names(sp))) sp$RDA1 <- 0
  if (!("PC1"  %in% names(sp))) sp$PC1  <- 0
  
  sp <- sp %>%
    mutate(
      arrow_len = sqrt(RDA1^2 + PC1^2),
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

plot_prda_panel <- function(pre, mod, stats_tbl, panel_title = "B",
                            top_n = 5, color_by = c("treatment","none")) {
  color_by <- match.arg(color_by)
  
  sites <- as.data.frame(scores(mod, display = "sites", scaling = 2))
  sites$Season <- pre$dat$Season
  sites$Irrigation <- pre$dat$Irrigation_trmt_pct
  
  ax <- get_axes(as.matrix(sites))
  xax <- ax[1]; yax <- ax[2]
  
  p <- ggplot(sites, aes(x = .data[[xax]], y = .data[[yax]]))
  
  if (color_by == "treatment") {
    p <- p +
      geom_point(aes(fill = Irrigation), shape = 21, color = "black",
                 size = 3, alpha = 0.9) +
      scale_fill_manual(values = trt_cols, name = "Irrigation (%)")
  } else {
    p <- p + geom_point(shape = 21, fill = "grey75", color = "black",
                        size = 3, alpha = 0.9)
  }
  
  sp <- as.data.frame(scores(mod, display = "species", scaling = 2))
  sp$Variable <- rownames(sp)
  if (!(xax %in% names(sp))) sp[[xax]] <- 0
  if (!(yax %in% names(sp))) sp[[yax]] <- 0
  
  sp2 <- sp %>%
    mutate(arrow_len = sqrt(.data[[xax]]^2 + .data[[yax]]^2)) %>%
    arrange(desc(arrow_len)) %>%
    slice_head(n = top_n)
  
  st <- stats_tbl[1, ]
  stat_label <- paste0(
    "Adj R\u00B2 = ", sprintf("%.3f", st$Adj_R2), "\n",
    "p(Resp) = ", sprintf("%.3f", st$p_resp), "\n",
    "n = ", st$n
  )
  
  p +
    geom_segment(
      data = sp2,
      aes(x = 0, y = 0, xend = .data[[xax]], yend = .data[[yax]]),
      inherit.aes = FALSE,
      arrow = arrow(length = unit(0.18, "cm")),
      linewidth = 0.7,
      color = "black",
      alpha = 0.9
    ) +
    ggrepel::geom_text_repel(
      data = sp2,
      aes(x = .data[[xax]], y = .data[[yax]], label = Variable),
      inherit.aes = FALSE,
      size = 3.2,
      max.overlaps = 50
    ) +
    annotate("text", x = Inf, y = Inf, label = stat_label,
             hjust = 1.05, vjust = 1.1, size = 3.4) +
    coord_equal() +
    theme_bw(base_size = 12) +
    labs(title = panel_title, x = xax, y = yax) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "none")
}

# ----------------------------
# 4) Panel A: Resp ~ GWC
# ----------------------------
panel_A <- df %>%
  select(Season, Irrigation_trmt_pct, Plot, Resp, GWC) %>%
  drop_na() %>%
  ggplot(aes(x = GWC, y = Resp)) +
  geom_point(aes(fill = Irrigation_trmt_pct),
             shape = 21, color = "black",
             size = 3, alpha = 0.9) +
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

# ----------------------------
# 5) Run models (auto log1p skewed predictors)
# ----------------------------
dir.create("Results", showWarnings = FALSE)
dir.create("Figures", showWarnings = FALSE)

run_one <- function(dat, label, vars, top_n = 5, nperm = 999, skew_thresh = 1.0) {
  pre <- prep_data(dat, vars, skew_thresh = skew_thresh, verbose = TRUE)
  mod <- fit_prda(pre)
  stats <- summarize_prda(mod, pre$dat, label, nperm = nperm)
  arrows <- arrow_table(mod, pre$Xs, pre$dat, top_n = top_n, nperm = nperm)
  list(pre = pre, mod = mod, stats = stats, arrows = arrows)
}

top_n <- 5
nperm <- 999
skew_thresh <- 1.0   # try 0.8–1.5 if you want more/less logging

out_all    <- run_one(df, "All", predictor_vars, top_n = top_n, nperm = nperm, skew_thresh = skew_thresh)
out_grow   <- run_one(df %>% filter(Season == "Growing"),
                      "Growing", predictor_vars, top_n = top_n, nperm = nperm, skew_thresh = skew_thresh)
out_nongro <- run_one(df %>% filter(Season == "Non-Growing"),
                      "Non-Growing", predictor_vars, top_n = top_n, nperm = nperm, skew_thresh = skew_thresh)

# Export model stats + arrow tables
write.csv(out_all$stats,    "Results/pRDA_subset_log_auto_stats_All.csv", row.names = FALSE)
write.csv(out_grow$stats,   "Results/pRDA_subset_log_auto_stats_Growing.csv", row.names = FALSE)
write.csv(out_nongro$stats, "Results/pRDA_subset_log_auto_stats_Non-Growing.csv", row.names = FALSE)

write.csv(out_all$arrows,    paste0("Results/pRDA_subset_log_auto_top", top_n, "_arrows_All.csv"), row.names = FALSE)
write.csv(out_grow$arrows,   paste0("Results/pRDA_subset_log_auto_top", top_n, "_arrows_Growing.csv"), row.names = FALSE)
write.csv(out_nongro$arrows, paste0("Results/pRDA_subset_log_auto_top", top_n, "_arrows_Non-Growing.csv"), row.names = FALSE)

all_stats <- bind_rows(out_all$stats, out_grow$stats, out_nongro$stats)
write.csv(all_stats, "Results/pRDA_subset_log_auto_model_stats_all.csv", row.names = FALSE)

# Also export which variables got log1p’d (so you can report it)
write.csv(out_all$pre$log_table,    "Results/pRDA_subset_log_auto_transform_table_All.csv", row.names = FALSE)
write.csv(out_grow$pre$log_table,   "Results/pRDA_subset_log_auto_transform_table_Growing.csv", row.names = FALSE)
write.csv(out_nongro$pre$log_table, "Results/pRDA_subset_log_auto_transform_table_Non-Growing.csv", row.names = FALSE)

# ----------------------------
# 6) Panels B/C
# ----------------------------
panel_B <- plot_prda_panel(out_grow$pre, out_grow$mod, out_grow$stats,
                           panel_title = "B", top_n = top_n, color_by = "treatment")

panel_C <- plot_prda_panel(out_nongro$pre, out_nongro$mod, out_nongro$stats,
                           panel_title = "C", top_n = top_n, color_by = "treatment")

# Combine layout
fig <- panel_A | (panel_B / panel_C) + plot_layout(widths = c(2, 1))

ggsave("Figures/Figure_ABC_RespGWC_pRDA_subset_log_auto_top5.png",
       fig, width = 15, height = 10, dpi = 300, bg = "white")
ggsave("Figures/Figure_ABC_RespGWC_pRDA_subset_log_auto_top5.pdf",
       fig, width = 15, height = 10, bg = "white")

cat("\nSaved figure:\n",
    "  Figures/Figure_ABC_RespGWC_pRDA_subset_log_auto_top5.png\n",
    "  Figures/Figure_ABC_RespGWC_pRDA_subset_log_auto_top5.pdf\n\n")

cat("Saved results:\n",
    "  Results/pRDA_subset_log_auto_model_stats_all.csv\n",
    "  Results/pRDA_subset_log_auto_stats_*.csv\n",
    "  Results/pRDA_subset_log_auto_top5_arrows_*.csv\n",
    "  Results/pRDA_subset_log_auto_transform_table_*.csv\n\n")

print(all_stats)
