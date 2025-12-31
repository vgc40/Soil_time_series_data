################################################################################
# Manuscript Figure: Panel A (Resp ~ GWC) + Panel B/C (partial RDA by season)
# Layout: A large on LEFT; B (top-right) + C (bottom-right)
# Arrows: top 5 only
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

# Treatment palette for Panel A (set as you like)
# (Pick distinct colors; here I match your season palette for visual harmony.)
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
  facet_wrap(~Season,scales = 'free_y', nrow = 2) +
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
# 5) Helper: partial RDA plot with top arrows (TOP 5)
# ----------------------------
plot_partial_rda <- function(dat, panel_title, top_arrows = 5, color_by = c("none","treatment","season")) {
  
  color_by <- match.arg(color_by)
  
  pre <- prep_data(dat)
  dat2 <- pre$dat
  Xs <- pre$Xs
  
  # partial RDA: X ~ Resp | Plot
  mod <- rda(Xs ~ Resp_log + Condition(Plot), data = dat2)
  
  # Site scores (often RDA1 + PC1 for single constraint)
  sites <- as.data.frame(scores(mod, display = "sites", choices = 1:2, scaling = 2)) %>%
    rename_with(make.names)
  sites$Season <- dat2$Season
  sites$Irrigation <- dat2$Irrigation_trmt_pct
  
  # Species scores (arrows)
  sp <- as.data.frame(scores(mod, display = "species", choices = 1:2, scaling = 2)) %>%
    rename_with(make.names)
  sp$var <- rownames(sp)
  
  # Axes present in sites
  site_axes <- intersect(c("RDA1","RDA2","PC1","PC2"), names(sites))[1:2]
  if (length(site_axes) < 2) stop("Could not find two axes for plotting.")
  x_axis <- site_axes[1]
  y_axis <- site_axes[2]
  
  # Ensure species has same axes
  if (!(x_axis %in% names(sp))) sp[[x_axis]] <- 0
  if (!(y_axis %in% names(sp))) sp[[y_axis]] <- 0
  
  # Arrow length -> pick top 5
  sp$arrow_len <- sqrt(sp[[x_axis]]^2 + sp[[y_axis]]^2)
  sp_top <- sp %>% arrange(desc(arrow_len)) %>% slice_head(n = top_arrows)
  
  p <- ggplot(sites, aes(x = .data[[x_axis]], y = .data[[y_axis]]))
  
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
      aes(x = .data[[x_axis]], y = .data[[y_axis]], label = var),
      inherit.aes = FALSE,
      size = 3.2,
      max.overlaps = 50
    ) +
    coord_equal() +
    theme_bw(base_size = 12) +
    labs(title = panel_title, x = x_axis, y = y_axis) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none"  # keep legend only in Panel A
    )
}

# ----------------------------
# 6) Panels B & C (partial RDA by season) with TOP 5 arrows
#    (Optional: color_by="treatment" if you want to see irrigation separation)
# ----------------------------
panel_B <- plot_partial_rda(
  df %>% filter(Season == "Growing"),
  panel_title = "B",
  top_arrows = 5,
  color_by = "treatment"
)

panel_C <- plot_partial_rda(
  df %>% filter(Season == "Non-Growing"),
  panel_title = "C",
  top_arrows = 5,
  color_by = "treatment"
)

# ----------------------------
# 7) Combine: A on left (bigger), B/C stacked on right
# ----------------------------
fig <- panel_A | (panel_B / panel_C) +
  plot_layout(widths = c(2, 1))

ggsave("Figures/Figure_ABC_RespGWC_pRDA_top5.png",
       fig, width = 15, height = 10, dpi = 300, bg = "white")

ggsave("Figures/Figure_ABC_RespGWC_pRDA_top5.pdf",
       fig, width = 15, height = 10, bg = "white")

cat("Saved:\n  Figures/Figure_ABC_RespGWC_pRDA_top5.png\n  Figures/Figure_ABC_RespGWC_pRDA_top5.pdf\n")
################################################################################
