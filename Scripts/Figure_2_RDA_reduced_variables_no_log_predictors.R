################################################################################
# Manuscript Figure (Subset pRDA):
# Panel A: Resp ~ GWC (by Season; legend once inside)
# Panel B: Growing partial RDA (Xs ~ Resp_log | Plot), top 5 arrows by |cor with RDA1|
# Panel C: Non-growing partial RDA (Xs ~ Resp_log | Plot), top 5 arrows by |cor with RDA1|
# Layout: A big on left; B top-right; C bottom-right
# Notes:
# - Predictors are the 11-variable subset
# - X scaled; Resp_log = log1p(resp)
# - Permutation tests blocked by Plot (repeated measures)
################################################################################

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(permute)
  library(ggrepel)
  library(patchwork)
})

# ----------------------------
# 1) Load data + prep
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

# Season colors you requested (used if you ever color by season)
season_cols <- c("Growing" = "#E69F00", "Non-Growing" = "#56B4E9")

# Treatment colors for Panel A and optional point fill in B/C
trt_cols <- c("100" = "#5AAE61", "25" = "#9970AB")

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

# ----------------------------
# 3) Helpers: prep + blocked perms
# ----------------------------
prep_data <- function(dat, vars) {
  dat2 <- dat %>%
    select(Season, Plot, Irrigation_trmt_pct, Resp, Resp_log, all_of(vars)) %>%
    drop_na()
  
  X  <- dat2 %>% select(all_of(vars))
  Xs <- scale(X)
  
  list(dat = dat2, X = X, Xs = Xs)
}

perm_ctrl <- function(plot_factor, nperm = 999) {
  how(nperm = nperm, blocks = plot_factor)
}

summarize_prda <- function(mod, dat_meta, label, nperm = 999) {
  ctrl <- perm_ctrl(dat_meta$Plot, nperm)
  
  a_overall <- anova.cca(mod, permutations = ctrl)
  a_term    <- anova.cca(mod, by = "term", permutations = ctrl)
  
  adj_r2 <- RsquareAdj(mod)$adj.r.squared
  frac_constrained <- mod$CCA$tot.chi / mod$tot.chi
  
  tibble(
    Dataset = label,
    n = nrow(dat_meta),
    n_plots = n_distinct(dat_meta$Plot),
    Adj_R2 = adj_r2,
    Constrained_fraction = frac_constrained,
    F_overall = a_overall$F[1],
    p_overall = a_overall$`Pr(>F)`[1],
    F_resp = a_term$F[1],
    p_resp = a_term$`Pr(>F)`[1]
  )
}

# ----------------------------
# 4) Panel A: Resp ~ GWC (faceted by Season)
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
# 5) Core plotting function: partial RDA by season
#    - arrows = top N by |cor_with_RDA1| (Resp relevance)
#    - flip RDA1 so +RDA1 aligns with higher Resp_log
#    - annotate Adj R2, p(Resp), n
# ----------------------------
plot_partial_rda_subset <- function(dat, vars, panel_title = "B",
                                    top_arrows = 5,
                                    color_by = c("treatment","season","none"),
                                    nperm = 999) {
  
  color_by <- match.arg(color_by)
  
  pre  <- prep_data(dat, vars)
  dat2 <- pre$dat
  Xs   <- pre$Xs
  
  # Fit partial RDA
  mod <- rda(Xs ~ Resp_log + Condition(Plot), data = dat2)
  
  # Stats (blocked perms)
  st <- summarize_prda(mod, dat2, label = panel_title, nperm = nperm)
  
  # Site scores (scaling 2)
  sites <- as.data.frame(scores(mod, display = "sites", scaling = 2))
  if (!("RDA1" %in% names(sites))) stop("Expected RDA1 in site scores.")
  sites$Resp_log <- dat2$Resp_log
  sites$Season <- dat2$Season
  sites$Irrigation <- dat2$Irrigation_trmt_pct
  
  # Determine 2 axes for plotting (typically RDA1 + PC1)
  axes <- intersect(c("RDA1","RDA2","PC1","PC2"), names(sites))[1:2]
  if (length(axes) < 2) stop("Could not find two ordination axes to plot.")
  x_axis <- axes[1]
  y_axis <- axes[2]
  
  # Species scores (arrows)
  sp <- as.data.frame(scores(mod, display = "species", scaling = 2))
  sp$Variable <- rownames(sp)
  
  # Ensure axes exist in species too
  if (!(x_axis %in% names(sp))) sp[[x_axis]] <- 0
  if (!(y_axis %in% names(sp))) sp[[y_axis]] <- 0
  
  # ---- Anchor axis direction for interpretability ----
  # Flip RDA1 so +RDA1 corresponds to higher Resp_log
  if ("RDA1" %in% names(sites)) {
    cor_rda1_resp <- suppressWarnings(cor(sites$RDA1, sites$Resp_log))
    if (is.finite(cor_rda1_resp) && cor_rda1_resp < 0) {
      sites$RDA1 <- -sites$RDA1
      if ("RDA1" %in% names(sp)) sp$RDA1 <- -sp$RDA1
    }
  }
  
  # ---- Choose arrows by Resp relevance ----
  # Correlation between each variable (Xs) and RDA1 site scores
  rda1 <- sites$RDA1
  sp$cor_with_RDA1 <- map_dbl(sp$Variable, ~ suppressWarnings(cor(Xs[, .x], rda1)))
  sp$abs_cor_with_RDA1 <- abs(sp$cor_with_RDA1)
  
  # For transparency (optional): total arrow length in plotted axes
  sp$arrow_len_plotspace <- sqrt(sp[[x_axis]]^2 + sp[[y_axis]]^2)
  
  sp_top <- sp %>%
    arrange(desc(abs_cor_with_RDA1)) %>%
    slice_head(n = top_arrows)
  
  # --- Add envfit stats for the displayed arrows (blocked perms) ---
  ctrl <- perm_ctrl(dat2$Plot, nperm)
  ef <- envfit(mod, Xs, permutations = ctrl)
  ef_tbl <- tibble(
    Variable = names(ef$vectors$pvals),
    envfit_r2 = as.numeric(ef$vectors$r^2),
    envfit_p  = as.numeric(ef$vectors$pvals)
  ) %>% mutate(envfit_p_fdr = p.adjust(envfit_p, method = "fdr"))
  
  sp_top <- sp_top %>% left_join(ef_tbl, by = "Variable")
  
  # Points layer
  p <- ggplot(sites, aes(x = .data[[x_axis]], y = .data[[y_axis]]))
  
  if (color_by == "treatment") {
    p <- p +
      geom_point(aes(fill = Irrigation), shape = 21, color = "black",
                 size = 3, alpha = 0.9) +
      scale_fill_manual(values = trt_cols)
  } else if (color_by == "season") {
    p <- p +
      geom_point(aes(fill = Season), shape = 21, color = "black",
                 size = 3, alpha = 0.9) +
      scale_fill_manual(values = season_cols)
  } else {
    p <- p +
      geom_point(shape = 21, fill = "grey75", color = "black",
                 size = 3, alpha = 0.9)
  }
  
  # Annotation string
  ann <- paste0(
    "Adj R² = ", sprintf("%.3f", st$Adj_R2), "\n",
    "p(Resp) = ", sprintf("%.3f", st$p_resp), "\n",
    "n = ", st$n
  )
  
  # Plot with arrows + labels
  p +
    geom_segment(
      data = sp_top,
      aes(x = 0, y = 0, xend = .data[[x_axis]], yend = .data[[y_axis]]),
      inherit.aes = FALSE,
      arrow = arrow(length = unit(0.18, "cm")),
      linewidth = 0.75,
      color = "black"
    ) +
    ggrepel::geom_text_repel(
      data = sp_top,
      aes(x = .data[[x_axis]], y = .data[[y_axis]], label = Variable),
      inherit.aes = FALSE,
      size = 3.2,
      max.overlaps = 50
    ) +
    annotate("text", x = Inf, y = Inf, label = ann,
             hjust = 1.08, vjust = 1.15, size = 3.5) +
    coord_equal() +
    theme_bw(base_size = 12) +
    labs(title = panel_title, x = x_axis, y = y_axis) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
}

# ----------------------------
# 6) Make Panels B/C using subset vars
# ----------------------------
panel_B <- plot_partial_rda_subset(
  df %>% filter(Season == "Growing"),
  vars = predictor_vars,
  panel_title = "B",
  top_arrows = 5,
  color_by = "treatment",
  nperm = 999
)

panel_C <- plot_partial_rda_subset(
  df %>% filter(Season == "Non-Growing"),
  vars = predictor_vars,
  panel_title = "C",
  top_arrows = 5,
  color_by = "treatment",
  nperm = 999
)

# ----------------------------
# 7) Combine layout: A big left, B/C stacked right
# ----------------------------
fig <- panel_A | (panel_B / panel_C) +
  plot_layout(widths = c(2, 1))

dir.create("Figures", showWarnings = FALSE)

ggsave("Figures/Figure_ABC_RespGWC_pRDA_subset_top5.png",
       fig, width = 15, height = 10, dpi = 300, bg = "white")

ggsave("Figures/Figure_ABC_RespGWC_pRDA_subset_top5.pdf",
       fig, width = 15, height = 10, bg = "white")

cat("Saved:\n",
    "  Figures/Figure_ABC_RespGWC_pRDA_subset_top5.png\n",
    "  Figures/Figure_ABC_RespGWC_pRDA_subset_top5.pdf\n")

################################################################################
# 8) Export stats + top 5 arrow tables (Resp-relevant ranking)
################################################################################
dir.create("Results", showWarnings = FALSE)

export_subset_stats_and_arrows <- function(dat, label, vars, top_n = 5, nperm = 999) {
  
  pre  <- prep_data(dat, vars)
  dat2 <- pre$dat
  Xs   <- pre$Xs
  
  mod <- rda(Xs ~ Resp_log + Condition(Plot), data = dat2)
  st  <- summarize_prda(mod, dat2, label, nperm = nperm)
  
  # Compute arrows ranked by |cor_with_RDA1| (Resp relevance)
  sites <- scores(mod, display = "sites", scaling = 2)
  rda1 <- sites[, "RDA1"]
  
  # Flip direction so +RDA1 aligns with higher Resp_log
  cor_rda1_resp <- suppressWarnings(cor(rda1, dat2$Resp_log))
  if (is.finite(cor_rda1_resp) && cor_rda1_resp < 0) rda1 <- -rda1
  
  sp <- as.data.frame(scores(mod, display = "species", scaling = 2))
  sp$Variable <- rownames(sp)
  if (!("RDA1" %in% names(sp))) sp$RDA1 <- 0
  if (!("PC1"  %in% names(sp))) sp$PC1  <- 0
  
  sp <- sp %>%
    mutate(
      cor_with_RDA1 = map_dbl(Variable, ~ suppressWarnings(cor(Xs[, .x], rda1))),
      abs_cor_with_RDA1 = abs(cor_with_RDA1),
      arrow_len_total = sqrt(RDA1^2 + PC1^2)
    ) %>%
    arrange(desc(abs_cor_with_RDA1)) %>%
    slice_head(n = top_n)
  
  ctrl <- perm_ctrl(dat2$Plot, nperm)
  ef <- envfit(mod, Xs, permutations = ctrl)
  
  ef_tbl <- tibble(
    Variable = names(ef$vectors$pvals),
    envfit_r2 = as.numeric(ef$vectors$r^2),
    envfit_p  = as.numeric(ef$vectors$pvals)
  ) %>%
    mutate(envfit_p_fdr = p.adjust(envfit_p, method = "fdr"))
  
  arrows <- sp %>% left_join(ef_tbl, by = "Variable")
  
  write.csv(st, file = paste0("Results/pRDA_subset_stats_", label, ".csv"), row.names = FALSE)
  write.csv(arrows, file = paste0("Results/pRDA_subset_top", top_n, "_arrows_", label, ".csv"), row.names = FALSE)
  
  list(model = mod, stats = st, arrows = arrows)
}

out_all    <- export_subset_stats_and_arrows(df, "All", predictor_vars, top_n = 5, nperm = 999)
out_grow   <- export_subset_stats_and_arrows(df %>% filter(Season=="Growing"), "Growing", predictor_vars, top_n = 5, nperm = 999)
out_nongro <- export_subset_stats_and_arrows(df %>% filter(Season=="Non-Growing"), "Non-Growing", predictor_vars, top_n = 5, nperm = 999)

all_stats <- bind_rows(out_all$stats, out_grow$stats, out_nongro$stats)
write.csv(all_stats, "Results/pRDA_subset_model_stats_all.csv", row.names = FALSE)

cat("\n=== Subset model stats (All/Growing/Non-Growing) ===\n")
print(all_stats)

cat("\n=== Top arrows (Growing) ===\n")
print(out_grow$arrows)

cat("\n=== Top arrows (Non-Growing) ===\n")
print(out_nongro$arrows)

################################################################################
# END
################################################################################
