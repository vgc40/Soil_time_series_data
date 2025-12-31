################################################################################
# FINAL: Manuscript Figure (Subset pRDA) + Stats Export
# Panel A: Resp ~ GWC (by Season; legend once inside)
# Panel B: Growing pRDA (Xs ~ Resp | Plot), top 5 arrows by |cor(X, RDA1)|
# Panel C: Non-growing pRDA (Xs ~ Resp | Plot), top 5 arrows by |cor(X, RDA1)|
# - Predictors = 11-variable subset
# - X scaled (raw X then scale)
# - Blocked permutations by Plot
# - Exports model stats + top5 arrows tables + combined stats table
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

season_cols <- c("Growing" = "#E69F00", "Non-Growing" = "#56B4E9")
trt_cols    <- c("100" = "#5AAE61", "25" = "#9970AB")

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

# Build top-N arrow table ranked by |cor(X, RDA1)| + envfit stats (blocked)
arrow_table <- function(mod, Xs, dat_meta, top_n = 5, nperm = 999) {
  
  ctrl <- perm_ctrl(dat_meta$Plot, nperm)
  
  # Site scores (scaling 2)
  sites <- scores(mod, display = "sites", scaling = 2)
  if (!("RDA1" %in% colnames(sites))) stop("No RDA1 in site scores.")
  rda1 <- sites[, "RDA1"]
  
  # Flip so +RDA1 corresponds to higher Resp (stored in dat_meta$Resp_used)
  cor_rda1_resp <- suppressWarnings(cor(rda1, dat_meta$Resp_used))
  if (is.finite(cor_rda1_resp) && cor_rda1_resp < 0) rda1 <- -rda1
  
  # Species (biplot) scores
  sp <- as.data.frame(scores(mod, display = "species", scaling = 2))
  sp$Variable <- rownames(sp)
  
  # Many 1-constraint models show RDA1 vs PC1
  if (!("RDA1" %in% names(sp))) sp$RDA1 <- 0
  if (!("PC1"  %in% names(sp))) sp$PC1  <- 0
  
  sp <- sp %>%
    mutate(
      cor_with_RDA1 = map_dbl(Variable, ~ suppressWarnings(cor(Xs[, .x], rda1))),
      abs_cor_with_RDA1 = abs(cor_with_RDA1),
      arrow_len_plot = sqrt(RDA1^2 + PC1^2)
    ) %>%
    arrange(desc(abs_cor_with_RDA1)) %>%
    slice_head(n = top_n)
  
  # envfit significance for each variable (blocked perms)
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
# 5) pRDA plot helper (also returns stats + arrows)
# ----------------------------
plot_prda_and_stats <- function(dat, vars, panel_title,
                                resp_col = c("Resp","Resp_log"),
                                top_arrows = 5,
                                color_by = c("treatment","season","none"),
                                nperm = 999) {
  
  resp_col <- match.arg(resp_col)
  color_by <- match.arg(color_by)
  
  pre  <- prep_data(dat, vars)
  dat2 <- pre$dat
  Xs   <- pre$Xs
  
  # store which Resp we used for flipping/arrow stats
  dat2$Resp_used <- dat2[[resp_col]]
  
  # Fit model
  form <- as.formula(paste0("Xs ~ ", resp_col, " + Condition(Plot)"))
  mod <- rda(form, data = dat2)
  
  # Stats + top arrows
  st     <- summarize_prda(mod, dat2, panel_title, nperm = nperm)
  arrows <- arrow_table(mod, Xs, dat2, top_n = top_arrows, nperm = nperm)
  
  # Plot scores
  sites <- as.data.frame(scores(mod, display = "sites", scaling = 2))
  if (!("RDA1" %in% names(sites))) stop("Expected RDA1 in site scores.")
  sites$Season    <- dat2$Season
  sites$Irrigation <- dat2$Irrigation_trmt_pct
  sites$Resp_used <- dat2$Resp_used
  
  # Force axes to RDA1 vs PC1 when available
  x_axis <- "RDA1"
  y_axis <- if ("PC1" %in% names(sites)) "PC1" else if ("RDA2" %in% names(sites)) "RDA2" else stop("No PC1 or RDA2.")
  
  # Species scores for arrows (ensure axes exist)
  sp <- as.data.frame(scores(mod, display = "species", scaling = 2))
  sp$Variable <- rownames(sp)
  if (!(x_axis %in% names(sp))) sp[[x_axis]] <- 0
  if (!(y_axis %in% names(sp))) sp[[y_axis]] <- 0
  
  # Keep only the arrows we selected (by Variable)
  sp_top <- sp %>% filter(Variable %in% arrows$Variable)
  
  # Flip direction so +RDA1 = higher Resp_used
  cor_rda1_resp <- suppressWarnings(cor(sites$RDA1, sites$Resp_used))
  if (is.finite(cor_rda1_resp) && cor_rda1_resp < 0) {
    sites$RDA1 <- -sites$RDA1
    sp_top$RDA1 <- -sp_top$RDA1
  }
  
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
  
  ann <- paste0(
    "Adj R² = ", sprintf("%.3f", st$Adj_R2), "\n",
    "p(Resp) = ", sprintf("%.3f", st$p_resp), "\n",
    "n = ", st$n
  )
  
  p <- p +
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
  
  list(plot = p, model = mod, stats = st, arrows = arrows)
}

# ----------------------------
# 6) Run seasonal pRDA panels (use raw Resp for B/C)
# ----------------------------
out_B <- plot_prda_and_stats(
  df %>% filter(Season == "Growing"),
  vars = predictor_vars,
  panel_title = "B",
  resp_col = "Resp",
  top_arrows = 5,
  color_by = "treatment",
  nperm = 999
)

out_C <- plot_prda_and_stats(
  df %>% filter(Season == "Non-Growing"),
  vars = predictor_vars,
  panel_title = "C",
  resp_col = "Resp",
  top_arrows = 5,
  color_by = "treatment",
  nperm = 999
)

panel_B <- out_B$plot
panel_C <- out_C$plot

# ----------------------------
# 7) Combine + save figure
# ----------------------------
fig <- panel_A | (panel_B / panel_C) +
  plot_layout(widths = c(2, 1))

dir.create("Figures", showWarnings = FALSE)

ggsave("Figures/Figure_ABC_RespGWC_pRDA_subset_top5.png",
       fig, width = 15, height = 10, dpi = 300, bg = "white")

ggsave("Figures/Figure_ABC_RespGWC_pRDA_subset_top5.pdf",
       fig, width = 15, height = 10, bg = "white")

cat("Saved figure:\n",
    "  Figures/Figure_ABC_RespGWC_pRDA_subset_top5.png\n",
    "  Figures/Figure_ABC_RespGWC_pRDA_subset_top5.pdf\n")

# ----------------------------
# 8) Export stats + arrows for All/Growing/Non-Growing
# ----------------------------
dir.create("Results", showWarnings = FALSE)

# Also export "All" stats/arrows (not plotted)
out_all <- plot_prda_and_stats(
  df,
  vars = predictor_vars,
  panel_title = "All",
  resp_col = "Resp",     # keep raw Resp for comparability; you can switch to Resp_log if desired
  top_arrows = 5,
  color_by = "none",
  nperm = 999
)

# Write individual files
write.csv(out_all$stats, "Results/pRDA_subset_stats_All.csv", row.names = FALSE)
write.csv(out_B$stats,   "Results/pRDA_subset_stats_Growing.csv", row.names = FALSE)
write.csv(out_C$stats,   "Results/pRDA_subset_stats_Non-Growing.csv", row.names = FALSE)

write.csv(out_all$arrows, "Results/pRDA_subset_top5_arrows_All.csv", row.names = FALSE)
write.csv(out_B$arrows,   "Results/pRDA_subset_top5_arrows_Growing.csv", row.names = FALSE)
write.csv(out_C$arrows,   "Results/pRDA_subset_top5_arrows_Non-Growing.csv", row.names = FALSE)

# Combined model stats table
all_stats <- bind_rows(out_all$stats, out_B$stats, out_C$stats)
write.csv(all_stats, "Results/pRDA_subset_model_stats_all.csv", row.names = FALSE)

cat("\nSaved results:\n",
    "  Results/pRDA_subset_stats_All.csv\n",
    "  Results/pRDA_subset_stats_Growing.csv\n",
    "  Results/pRDA_subset_stats_Non-Growing.csv\n",
    "  Results/pRDA_subset_top5_arrows_All.csv\n",
    "  Results/pRDA_subset_top5_arrows_Growing.csv\n",
    "  Results/pRDA_subset_top5_arrows_Non-Growing.csv\n",
    "  Results/pRDA_subset_model_stats_all.csv\n")

print(all_stats)
