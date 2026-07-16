############################################################
# Targeted MicroResp-to-BRITE C temporal analysis
#
# Required first step:
#   Run Scripts/09_BRITE_GAM_random_plot_all_pathways.R
#
# This script maps Data/Microresp_pathways.csv to the BRITE C GAM results by
# exact ko##### pathway code. It retains every matched target pathway whether
# or not the all-pathway temporal test was significant, and also calculates a
# target-set FDR across the unique matched BRITE C pathways.
############################################################

rm(list = ls())

required_packages <- c("tidyverse", "readr")
missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org",
    dependencies = TRUE
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
})

theme_set(theme_bw(base_size = 11))

# ---- Paths ------------------------------------------------------------------

all_gam_dir <- file.path(
  "Results",
  "BRITE_temporal_gam_random_plot_all_pathways"
)

microresp_file <- file.path("Data", "Microresp_pathways.csv")
results_dir <- file.path("Results", "MicroResp_targeted_BRITE_C_gam")
figures_dir <- file.path("Figures", "MicroResp_targeted_BRITE_C_gam")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

required_inputs <- c(
  microresp_file,
  file.path(all_gam_dir, "brite_C_GAM_summary_table.csv"),
  file.path(all_gam_dir, "brite_C_GAM_fitted_trajectories.csv"),
  file.path(all_gam_dir, "brite_C_GAM_long_model_data.csv")
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs) > 0) {
  stop(
    "Missing required input(s): ",
    paste(missing_inputs, collapse = ", "),
    ". Run Scripts/09_BRITE_GAM_random_plot_all_pathways.R first."
  )
}

# ---- Read and standardize mappings -----------------------------------------

microresp_map <- read_csv(microresp_file, show_col_types = FALSE) %>%
  mutate(
    Substrate = as.character(Substrate),
    brite_C = as.character(brite_C),
    pathway_code = str_extract(brite_C, "ko\\d{5}"),
    pathway_numeric_code = str_remove(pathway_code, "^ko")
  ) %>%
  filter(!is.na(Substrate), Substrate != "") %>%
  distinct()

if (anyNA(microresp_map$pathway_code)) {
  warning(
    "Some MicroResp mapping rows lack a ko##### BRITE C code. ",
    "They will be retained in the mapping audit but cannot be matched."
  )
}

gam_summary <- read_csv(
  file.path(all_gam_dir, "brite_C_GAM_summary_table.csv"),
  show_col_types = FALSE
) %>%
  mutate(pathway_code = str_extract(brite_class, "ko\\d{5}"))

gam_predictions <- read_csv(
  file.path(all_gam_dir, "brite_C_GAM_fitted_trajectories.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    Collection_date = as.Date(Collection_date),
    pathway_code = str_extract(brite_class, "ko\\d{5}")
  )

gam_model_data <- read_csv(
  file.path(all_gam_dir, "brite_C_GAM_long_model_data.csv"),
  show_col_types = FALSE
) %>%
  mutate(
    Collection_date = as.Date(Collection_date),
    pathway_code = str_extract(brite_class, "ko\\d{5}")
  )

# ---- Mapping audit ----------------------------------------------------------

mapping_audit <- microresp_map %>%
  left_join(
    gam_summary %>%
      select(
        pathway_code,
        matched_brite_class = brite_class,
        fit_status,
        p_value,
        p_adj_fdr_all_brite_C = p_adj_fdr,
        significant_time_all_brite_C = significant_time,
        deviance_explained_full,
        time_variance_explained
      ),
    by = "pathway_code"
  ) %>%
  mutate(
    mapping_status = case_when(
      is.na(pathway_code) ~ "No ko code in MicroResp mapping",
      is.na(matched_brite_class) ~ "No match in filtered BRITE C data",
      fit_status == "Converged" ~ "Matched; model converged",
      TRUE ~ paste0("Matched; model status: ", fit_status)
    )
  ) %>%
  arrange(Substrate, pathway_code)

write_csv(
  mapping_audit,
  file.path(results_dir, "MicroResp_to_BRITE_C_mapping_audit.csv")
)

coverage_summary <- mapping_audit %>%
  summarise(
    n_mapping_rows = n(),
    n_substrates = n_distinct(Substrate),
    n_unique_requested_pathways = n_distinct(pathway_code, na.rm = TRUE),
    n_mapping_rows_matched = sum(!is.na(matched_brite_class)),
    n_unique_pathways_matched = n_distinct(
      pathway_code[!is.na(matched_brite_class)],
      na.rm = TRUE
    ),
    n_unique_converged_pathways = n_distinct(
      pathway_code[fit_status == "Converged"],
      na.rm = TRUE
    )
  )

write_csv(
  coverage_summary,
  file.path(results_dir, "MicroResp_BRITE_C_mapping_coverage.csv")
)

# ---- Unique target-pathway inference ---------------------------------------

target_pathways <- microresp_map %>%
  filter(!is.na(pathway_code)) %>%
  group_by(pathway_code, brite_C) %>%
  summarise(
    n_linked_substrates = n_distinct(Substrate),
    linked_substrates = paste(sort(unique(Substrate)), collapse = "; "),
    substrate_categories = paste(
      sort(unique(Substrate_category)),
      collapse = "; "
    ),
    table_groups = paste(sort(unique(Table_group)), collapse = "; "),
    .groups = "drop"
  ) %>%
  left_join(gam_summary, by = "pathway_code") %>%
  mutate(
    p_value_target_set = if_else(
      fit_status == "Converged",
      p_value,
      NA_real_
    ),
    p_adj_fdr_microresp_targets = p.adjust(
      p_value_target_set,
      method = "BH"
    ),
    significant_time_microresp_targets = if_else(
      !is.na(p_adj_fdr_microresp_targets) &
        p_adj_fdr_microresp_targets < 0.05,
      "Yes",
      "No"
    )
  ) %>%
  select(-p_value_target_set) %>%
  arrange(p_adj_fdr_microresp_targets, p_value, pathway_code)

write_csv(
  target_pathways,
  file.path(results_dir, "MicroResp_target_BRITE_C_GAM_summary.csv")
)

write_csv(
  filter(target_pathways, significant_time_microresp_targets == "Yes"),
  file.path(
    results_dir,
    "MicroResp_target_BRITE_C_significant_temporal_pathways.csv"
  )
)

# ---- Targeted model data and trajectories ----------------------------------

matched_codes <- target_pathways %>%
  filter(!is.na(brite_class)) %>%
  pull(pathway_code) %>%
  unique()

target_model_data <- gam_model_data %>%
  filter(pathway_code %in% matched_codes) %>%
  left_join(
    target_pathways %>%
      select(
        pathway_code,
        linked_substrates,
        significant_time_microresp_targets,
        p_adj_fdr_microresp_targets
      ),
    by = "pathway_code"
  )

target_predictions <- gam_predictions %>%
  filter(pathway_code %in% matched_codes) %>%
  left_join(
    target_pathways %>%
      select(
        pathway_code,
        linked_substrates,
        significant_time_microresp_targets,
        p_adj_fdr_microresp_targets
      ),
    by = "pathway_code"
  )

write_csv(
  target_model_data,
  file.path(results_dir, "MicroResp_target_BRITE_C_long_model_data.csv")
)

write_csv(
  target_predictions,
  file.path(results_dir, "MicroResp_target_BRITE_C_fitted_trajectories.csv")
)

# A substrate-expanded table preserves the explicit many-to-many mapping for
# later comparison with the actual MicroResp substrate measurements.
substrate_pathway_trajectories <- microresp_map %>%
  select(
    Substrate,
    Substrate_category,
    Table_group,
    pathway_code,
    mapped_brite_C = brite_C
  ) %>%
  inner_join(target_predictions, by = "pathway_code")

write_csv(
  substrate_pathway_trajectories,
  file.path(
    results_dir,
    "MicroResp_substrate_to_BRITE_C_fitted_trajectories.csv"
  )
)

# ---- Figures ----------------------------------------------------------------

pathway_labels <- target_pathways %>%
  transmute(
    pathway_code,
    plot_label = paste0(
      str_remove(brite_C, "\\s+\\[PATH:ko\\d{5}\\]$"),
      "\nMicroResp: ",
      linked_substrates
    )
  )

plot_points <- target_model_data %>%
  left_join(pathway_labels, by = "pathway_code")

plot_lines <- target_predictions %>%
  left_join(pathway_labels, by = "pathway_code")

n_target_pathways <- n_distinct(plot_lines$pathway_code)
plot_height <- max(8, ceiling(n_target_pathways / 3) * 2.8)

all_target_plot <- ggplot(plot_points, aes(Collection_date, Expression)) +
  geom_point(alpha = 0.4, size = 1, color = "grey35") +
  geom_ribbon(
    data = plot_lines,
    aes(Collection_date, ymin = lower, ymax = upper),
    inherit.aes = FALSE,
    fill = "#9ecae1",
    alpha = 0.35
  ) +
  geom_line(
    data = plot_lines,
    aes(Collection_date, fitted),
    inherit.aes = FALSE,
    color = "#08519c",
    linewidth = 0.9
  ) +
  facet_wrap(~ plot_label, scales = "free_y", ncol = 3) +
  labs(
    title = "BRITE C pathways linked to MicroResp substrates",
    subtitle = paste0(
      "All matched pathways are shown, regardless of temporal significance. ",
      "Lines exclude Plot random effects."
    ),
    x = "Collection date",
    y = "BRITE C expression"
  ) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 7)
  )

ggsave(
  file.path(figures_dir, "MicroResp_target_BRITE_C_all_temporal_GAMs.png"),
  all_target_plot,
  width = 14,
  height = plot_height,
  dpi = 300,
  limitsize = FALSE
)

ggsave(
  file.path(figures_dir, "MicroResp_target_BRITE_C_all_temporal_GAMs.pdf"),
  all_target_plot,
  width = 14,
  height = plot_height,
  limitsize = FALSE
)

significant_codes <- target_pathways %>%
  filter(significant_time_microresp_targets == "Yes") %>%
  pull(pathway_code)

if (length(significant_codes) > 0) {
  significant_plot <- ggplot(
    filter(plot_points, pathway_code %in% significant_codes),
    aes(Collection_date, Expression)
  ) +
    geom_point(alpha = 0.4, size = 1, color = "grey35") +
    geom_ribbon(
      data = filter(plot_lines, pathway_code %in% significant_codes),
      aes(Collection_date, ymin = lower, ymax = upper),
      inherit.aes = FALSE,
      fill = "#a1d99b",
      alpha = 0.35
    ) +
    geom_line(
      data = filter(plot_lines, pathway_code %in% significant_codes),
      aes(Collection_date, fitted),
      inherit.aes = FALSE,
      color = "#238b45",
      linewidth = 0.9
    ) +
    facet_wrap(~ plot_label, scales = "free_y", ncol = 3) +
    labs(
      title = "Temporally significant MicroResp-linked BRITE C pathways",
      subtitle = "Significance uses FDR across the targeted unique BRITE C pathway set",
      x = "Collection date",
      y = "BRITE C expression"
    ) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 7)
    )

  significant_height <- max(7, ceiling(length(significant_codes) / 3) * 2.8)

  ggsave(
    file.path(
      figures_dir,
      "MicroResp_target_BRITE_C_significant_temporal_GAMs.png"
    ),
    significant_plot,
    width = 14,
    height = significant_height,
    dpi = 300,
    limitsize = FALSE
  )

  ggsave(
    file.path(
      figures_dir,
      "MicroResp_target_BRITE_C_significant_temporal_GAMs.pdf"
    ),
    significant_plot,
    width = 14,
    height = significant_height,
    limitsize = FALSE
  )
}

# Heatmap of pathway-specific fitted shape. Each pathway is standardized across
# dates so the heatmap compares temporal shape rather than expression magnitude.
heatmap_data <- target_predictions %>%
  group_by(pathway_code, brite_class) %>%
  mutate(
    fitted_sd = sd(fitted, na.rm = TRUE),
    fitted_z = if_else(
      is.finite(fitted_sd) & fitted_sd > 0,
      (fitted - mean(fitted, na.rm = TRUE)) / fitted_sd,
      0
    )
  ) %>%
  ungroup() %>%
  mutate(
    pathway_label = str_remove(
      brite_class,
      "\\s+\\[PATH:ko\\d{5}\\]$"
    ),
    pathway_label = fct_reorder(pathway_label, p_adj_fdr_microresp_targets),
    date_label = format(Collection_date, "%b %Y")
  )

trajectory_heatmap <- ggplot(
  heatmap_data,
  aes(date_label, pathway_label, fill = fitted_z)
) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    name = "Fitted\nz-score"
  ) +
  labs(
    title = "Temporal shapes of MicroResp-linked BRITE C pathways",
    subtitle = "Each pathway is standardized across dates",
    x = NULL,
    y = NULL
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 8)
  )

ggsave(
  file.path(figures_dir, "MicroResp_target_BRITE_C_trajectory_heatmap.png"),
  trajectory_heatmap,
  width = 10,
  height = max(7, 0.3 * n_target_pathways),
  dpi = 300,
  limitsize = FALSE
)

ggsave(
  file.path(figures_dir, "MicroResp_target_BRITE_C_trajectory_heatmap.pdf"),
  trajectory_heatmap,
  width = 10,
  height = max(7, 0.3 * n_target_pathways),
  limitsize = FALSE
)

run_manifest <- tibble(
  mapping_file = microresp_file,
  all_pathway_gam_directory = all_gam_dir,
  match_key = "Exact ko##### pathway code extracted from brite_C",
  inclusion_rule = "All matched MicroResp-linked BRITE C pathways",
  targeted_fdr_scope = "Unique matched, converged BRITE C pathways",
  results_dir = results_dir,
  figures_dir = figures_dir,
  run_date = as.character(Sys.Date())
)

write_csv(run_manifest, file.path(results_dir, "run_manifest.csv"))

cat("\nMicroResp-targeted BRITE C analysis complete.\n")
print(coverage_summary)
cat(
  "Target pathways significant after targeted FDR: ",
  sum(target_pathways$significant_time_microresp_targets == "Yes"),
  "\n",
  sep = ""
)
cat("Results: ", results_dir, "\n", sep = "")
cat("Figures: ", figures_dir, "\n", sep = "")

