################################################################################
# Manuscript Figure: Panel A (Resp ~ GWC) + Panel B/C (partial RDA by season)
# Layout: A large on LEFT; B (top-right) + C (bottom-right)
# Arrows: top 5 (NOW selected by |cor with Resp-axis (RDA1)|, not total arrow length)
# Legend: ONLY once, INSIDE Panel A
# Theme: theme_bw(); square scaling for ordinations via coord_equal()
################################################################################

rm(list = ls())
set.seed(123)

suppressPackageStartupMessages({
  library(tidyverse)
  library(vegan)
  library(ggrepel)
  library(patchwork)
  library(scales)   # alpha()
  library(permute)
})

# ----------------------------
# 1) Load + basic prep
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

# Requested season palette (kept for consistency if needed)
season_cols <- c("Growing" = "#E69F00", "Non-Growing" = "#56B4E9")

# Treatment palette for Panel A
trt_cols <- c("100" = "#5AAE61", "25" = "#9970AB")

# ----------------------------
# 2) Biogeochem variables (RDA matrix)
# ----------------------------
bio_vars <- c(
  "GWC","pH","SEOC_ug_gdw","SETN_ug_gdw","MBC_ug_gdw_bgc","MBN_ug_gdw_bgc",
  "Olson_IP_ug_gdw","MBP_IP_ug_gdw","BG_nmol_gdw_h","BX_nmol_gdw_h","CB_nmol_gdw_h",
  "NAG_nmol_gdw_h","Leu_nmol_gdw_h","Ala_nmol_gdw_h","AP_nmol_gdw_h","EEA_C_N",
  "EEA_C_P","Sum_Carbon","Roots_gdw","SoilTemp_mean"
)

need <- c("Season","Plot","Irrigation_trmt_pct","Resp_log","GWC", bio_vars)
miss <- setdiff(need, names(df))
if (length(miss) > 0) stop("Missing columns: ", paste(miss, collapse = ", "))

# ----------------------------
# 3) Helper: prep (drop NA, scale X)
# ----------------------------
prep_data <- function(dat) {
  dat2 <- dat %>%
    select(Season, Plot, Irrigation_trmt_pct, Resp_log, all_of(bio_vars)) %>%
    drop_na()
  
  X <- dat2 %>% select(all_of(bio_vars))
  Xs <- scale(X)
  
  list(dat = dat2, Xs = Xs)
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
    # Legend INSIDE panel A (bottom-left)
    legend.position = c(0.02, 0.02),
    legend.justification = c(0, 0),
    legend.key = element_rect(fill = NA),
    legend.title = element_text(face = "bold")
  )

# ----------------------------
# 5) Permutation control: blocked within Plot
# ----------------------------
perm_ctrl <- function(plot_factor, nperm = 999) {
  how(nperm = nperm, blocks = plot_factor)
}

# ----------------------------
# 6) Summarize pRDA model + blocked permutation tests
# ----------------------------
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
# 7) Arrow table: now ranks by |cor with RDA1| (Resp axis)
# ----------------------------
arrow_table <- function(mod, Xs, dat_meta, top_n = 5, nperm = 999) {
  
  ctrl <- perm_ctrl(dat_meta$Plot, nperm)
  
  # Site scores: constrained axis is RDA1
  sites <- scores(mod, display = "sites", scaling = 2)
  rda1  <- sites[, "RDA1"]
  
  # Species scores = biplot arrows
  sp <- as.data.frame(scores(mod, display = "species", scaling = 2))
  sp$Variable <- rownames(sp)
  
  if (!("RDA1" %in% names(sp))) sp$RDA1 <- 0
  if (!("PC1"  %in% names(sp))) sp$PC1  <- 0
  
  # Correlation of each biogeochem variable with Resp-axis
  sp <- sp %>%
    mutate(
      cor_with_RDA1 = map_dbl(Variable, ~ cor(Xs[, .x], rda1)),
      abs_cor_with_RDA1 = abs(cor_with_RDA1),
      arrow_len_total = sqrt(RDA1^2 + PC1^2),
      angle_deg = atan2(PC1, RDA1) * 180 / pi
    ) %>%
    arrange(desc(abs_cor_with_RDA1)) %>%
    slice_head(n = top_n)
  
  # envfit significance for each variable vector (blocked perms)
  ef <- envfit(mod, Xs, permutations = ctrl)
  ef_r2 <- ef$vectors$r^2
  ef_p  <- ef$vectors$pvals
  
  ef_tbl <- tibble(
    Variable = names(ef_p),
    envfit_r2 = as.numeric(ef_r2),
    envfit_p = as.numeric(ef_p)
  ) %>%
    mutate(envfit_p_fdr = p.adjust(envfit_p, method = "fdr"))
  
  sp %>%
    left_join(ef_tbl, by = "Variable") %>%
    mutate(
      aligned_with_pos_RDA1 = RDA1 > 0
    )
}

# ----------------------------
# 8) Helper: partial RDA plot with top arrows
#    (Top 5 by |cor_with_RDA1|; axis flipped so +RDA1 = higher Resp_log)
# ----------------------------
plot_partial_rda <- function(dat, panel_title, top_arrows = 5,
                             color_by = c("none","treatment","season"),
                             nperm = 999) {
  
  color_by <- match.arg(color_by)
  
  pre <- prep_data(dat)
  dat2 <- pre$dat
  Xs <- pre$Xs
  
  # partial RDA: X ~ Resp | Plot
  mod <- rda(Xs ~ Resp_log + Condition(Plot), data = dat2)
  
  # Stats for annotation (blocked permutations)
  stats <- summarize_prda(mod, dat2, label = panel_title, nperm = nperm)
  
  # Site scores
  sites_raw <- as.data.frame(scores(mod, display = "sites", choices = 1:2, scaling = 2))
  sites_raw$Season <- dat2$Season
  sites_raw$Irrigation <- dat2$Irrigation_trmt_pct
  sites_raw$Resp_log <- dat2$Resp_log
  
  # Identify axes used (typically RDA1 + PC1)
  site_axes <- intersect(c("RDA1","RDA2","PC1","PC2"), names(sites_raw))[1:2]
  if (length(site_axes) < 2) stop("Could not find two axes for plotting.")
  x_axis <- site_axes[1]
  y_axis <- site_axes[2]
  
  # Species scores
  sp_raw <- as.data.frame(scores(mod, display = "species", choices = 1:2, scaling = 2))
  sp_raw$Variable <- rownames(sp_raw)
  
  # Ensure both axes exist in species too
  if (!(x_axis %in% names(sp_raw))) sp_raw[[x_axis]] <- 0
  if (!(y_axis %in% names(sp_raw))) sp_raw[[y_axis]] <- 0
  
  # ---- Make axis direction consistent: +RDA1 means higher Resp_log ----
  # Only do this if x-axis is RDA1
  if (x_axis == "RDA1") {
    cor_rda1_resp <- cor(sites_raw$RDA1, sites_raw$Resp_log)
    if (is.finite(cor_rda1_resp) && cor_rda1_resp < 0) {
      sites_raw$RDA1 <- -sites_raw$RDA1
      if ("RDA1" %in% names(sp_raw)) sp_raw$RDA1 <- -sp_raw$RDA1
    }
  }
  
  # ---- Select top arrows by Resp-relevance: |cor_with_RDA1| ----
  # (Compute correlation using the potentially flipped RDA1 site scores)
  rda1_sites <- if ("RDA1" %in% names(sites_raw)) sites_raw$RDA1 else NULL
  if (is.null(rda1_sites)) stop("RDA1 site scores not found; cannot rank arrows by Resp-axis.")
  
  sp_raw$cor_with_RDA1 <- map_dbl(sp_raw$Variable, ~ cor(Xs[, .x], rda1_sites))
  sp_raw$abs_cor_with_RDA1 <- abs(sp_raw$cor_with_RDA1)
  
  # For transparency, also compute total arrow length in this plot space
  sp_raw$arrow_len_total <- sqrt(sp_raw[[x_axis]]^2 + sp_raw[[y_axis]]^2)
  
  sp_top <- sp_raw %>%
    arrange(desc(abs_cor_with_RDA1)) %>%
    slice_head(n = top_arrows)
  
  # Build plot
  p <- ggplot(sites_raw, aes(x = .data[[x_axis]], y = .data[[y_axis]]))
  
  if (color_by == "treatment") {
    p <- p +
      geom_point(aes(fill = Irrigation), shape = 21, color = "black", size = 3, alpha = 0.9) +
      scale_fill_manual(values = trt_cols, name = "Irrigation (%)")
  } else if (color_by == "season") {
    p <- p +
      geom_point(aes(fill = Season), shape = 21, color = "black", size = 3, alpha = 0.9) +
      scale_fill_manual(values = season_cols, name = "Season")
  } else {
    p <- p +
      geom_point(shape = 21, fill = "grey75", color = "black", size = 3, alpha = 0.9)
  }
  
  # Annotation text (model-level)
  ann <- paste0(
    "Adj R² = ", sprintf("%.3f", stats$Adj_R2), "\n",
    "p(Resp) = ", sprintf("%.3f", stats$p_resp), "\n",
    "n = ", stats$n
  )
  
  p +
    geom_segment(
      data = sp_top,
      aes(x = 0, y = 0, xend = .data[[x_axis]], yend = .data[[y_axis]]),
      inherit.aes = FALSE,
      arrow = arrow(length = unit(0.18, "cm")),
      linewidth = 0.7,
      color = "black",
      alpha = 0.9
    ) +
    ggrepel::geom_text_repel(
      data = sp_top,
      aes(x = .data[[x_axis]], y = .data[[y_axis]], label = Variable),
      inherit.aes = FALSE,
      size = 3.2,
      max.overlaps = 50
    ) +
    annotate("text",
             x = Inf, y = Inf,
             label = ann,
             hjust = 1.08, vjust = 1.2,
             size = 3.5) +
    coord_equal() +
    theme_bw(base_size = 12) +
    labs(
      title = panel_title,
      x = x_axis,
      y = y_axis
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none"  # keep legend only in Panel A
    )
}

# ----------------------------
# 9) Panels B & C (partial RDA by season)
#    (Color by treatment optional; you chose treatment)
# ----------------------------
panel_B <- plot_partial_rda(
  df %>% filter(Season == "Growing"),
  panel_title = "B",
  top_arrows = 5,
  color_by = "treatment",
  nperm = 999
)

panel_C <- plot_partial_rda(
  df %>% filter(Season == "Non-Growing"),
  panel_title = "C",
  top_arrows = 5,
  color_by = "treatment",
  nperm = 999
)

# ----------------------------
# 10) Combine: A on left (bigger), B/C stacked on right
# ----------------------------
fig <- panel_A | (panel_B / panel_C) +
  plot_layout(widths = c(2, 1))

dir.create("Figures", showWarnings = FALSE)

ggsave("Figures/Figure_ABC_RespGWC_pRDA_top5.png",
       fig, width = 15, height = 10, dpi = 300, bg = "white")

ggsave("Figures/Figure_ABC_RespGWC_pRDA_top5.pdf",
       fig, width = 15, height = 10, bg = "white")

cat("Saved:\n  Figures/Figure_ABC_RespGWC_pRDA_top5.png\n  Figures/Figure_ABC_RespGWC_pRDA_top5.pdf\n")

################################################################################
# ==== Stats exports (blocked perms; top 5 by |cor_with_RDA1|) ====
################################################################################
dir.create("Results", showWarnings = FALSE)

run_and_export <- function(dat, label, top_n = 5, nperm = 999) {
  
  pre <- prep_data(dat)
  mod <- rda(pre$Xs ~ Resp_log + Condition(Plot), data = pre$dat)
  
  stats <- summarize_prda(mod, pre$dat, label, nperm = nperm)
  arrows <- arrow_table(mod, pre$Xs, pre$dat, top_n = top_n, nperm = nperm)
  
  write.csv(stats,  file = paste0("Results/pRDA_stats_", gsub("\\s+","_",label), ".csv"),
            row.names = FALSE)
  write.csv(arrows, file = paste0("Results/pRDA_top", top_n, "_arrows_", gsub("\\s+","_",label), ".csv"),
            row.names = FALSE)
  
  list(model = mod, stats = stats, arrows = arrows)
}

out_all    <- run_and_export(df, "All", top_n = 5)
out_grow   <- run_and_export(df %>% filter(Season=="Growing"), "Growing", top_n = 5)
out_nongro <- run_and_export(df %>% filter(Season=="Non-Growing"), "Non-Growing", top_n = 5)

all_stats <- bind_rows(out_all$stats, out_grow$stats, out_nongro$stats)
write.csv(all_stats, "Results/pRDA_model_stats_all.csv", row.names = FALSE)

print(all_stats)
print(out_grow$arrows)
print(out_nongro$arrows)

################################################################################
# END
################################################################################
