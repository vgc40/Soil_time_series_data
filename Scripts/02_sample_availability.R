# =============================================================================
# Sample availability and matching audit for the final manuscript figures
# =============================================================================
#
# This script does not modify source data. It inventories the 100% irrigation
# samples used to connect field respiration, root biomass, MicroResp, and
# metatranscriptomes. Inputs: cleaned field table (Scripts/01), MicroResp v3
# (DataHub March), and the sample columns of Data/brite_C_expression.csv.

rm(list = ls())

required_packages <- c("tidyverse", "lubridate")
missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

result_dir <- file.path("Results", "manuscript_v4", "sample_availability")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

field_file <- file.path("Data", "Data_merged_cleaned.csv")
microresp_file <- file.path(
  "DataHub March",
  "Soil_Phenotypes",
  "SoilMicrobiomeSFA_prosser-time-series_MicroResp_v3.csv"
)
transcript_file <- file.path("Data", "brite_C_expression.csv")

required_files <- c(field_file, microresp_file, transcript_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop("Missing required input files: ", paste(missing_files, collapse = "; "))
}

# ---- Field respiration and root biomass ------------------------------------

field_100 <- read_csv(field_file, show_col_types = FALSE) %>%
  filter(as.character(Irrigation_trmt_pct) == "100") %>%
  mutate(
    Sampling_Date = as.Date(Collection_date),
    Plot_ID = as.character(Plot),
    Block = as.character(Block),
    Season = if_else(
      month(Sampling_Date) %in% 4:10,
      "Plant-active, Apr-Oct",
      "Plant-dormant, Nov-Mar"
    )
  ) %>%
  arrange(Sampling_Date, Plot_ID)

field_duplicate_records <- field_100 %>%
  count(Sample_cat, Sampling_Date, Plot_ID, name = "record_count") %>%
  filter(record_count > 1)

field_sample_table <- field_100 %>%
  distinct(
    Sample_cat,
    Sampling_Date,
    Plot_ID,
    Block,
    Season,
    Roots_gdw,
    Resp_rate_mg_C_m2_h
  )

field_dates <- field_sample_table %>%
  distinct(Sampling_Date) %>%
  arrange(Sampling_Date)

field_plots <- field_sample_table %>%
  distinct(Plot_ID, Block) %>%
  arrange(as.numeric(Plot_ID))

expected_field_design <- crossing(
  field_dates,
  field_plots
)

missing_expected_field_records <- expected_field_design %>%
  anti_join(
    field_sample_table %>% select(Sampling_Date, Plot_ID),
    by = c("Sampling_Date", "Plot_ID")
  )

# ---- MicroResp --------------------------------------------------------------

microresp_100 <- read_csv(microresp_file, show_col_types = FALSE, na = c("", "NA", "na")) %>%
  rename(
    Irrigation_trmt_cat = any_of("Irrigaiton_trmt_cat")
  ) %>%
  mutate(
    Sampling_Date = mdy(Collection_date),
    Timepoint = suppressWarnings(as.numeric(Timepoint)),
    SIR = suppressWarnings(as.numeric(SIR)),
    Plot_ID = as.character(Plot),
    Block = as.character(Block)
  ) %>%
  filter(
    as.character(Irrigation_trmt_pct) == "100",
    Timepoint == 5,
    Normalization == "ngCO2-C_mgC_h",
    Substrate != "Water"
  )

microresp_sample_table <- microresp_100 %>%
  distinct(Sample_cat, Sampling_Date, Plot_ID, Block) %>%
  mutate(
    Season = if_else(
      month(Sampling_Date) %in% 4:10,
      "Plant-active, Apr-Oct",
      "Plant-dormant, Nov-Mar"
    )
  ) %>%
  arrange(Sampling_Date, Plot_ID)

microresp_duplicate_measurements <- microresp_100 %>%
  count(Sample_cat, Substrate, name = "record_count") %>%
  filter(record_count > 1)

# The original MicroResp release had January 2024 Sample_cat labels shifted by
# one plot. Stop if any MicroResp sample disagrees with the field table on
# plot or date, so a mislabelled release cannot silently enter the analysis.
microresp_label_conflicts <- microresp_sample_table %>%
  inner_join(
    field_sample_table %>% select(Sample_cat, field_date = Sampling_Date,
                                  field_plot = Plot_ID),
    by = "Sample_cat"
  ) %>%
  filter(Sampling_Date != field_date | Plot_ID != field_plot)

write_csv(
  microresp_label_conflicts,
  file.path(result_dir, "microresp_label_conflicts_with_field_table.csv")
)

if (nrow(microresp_label_conflicts) > 0) {
  stop("MicroResp Sample_cat labels disagree with the field table; see ",
       "microresp_label_conflicts_with_field_table.csv")
}

missing_expected_microresp_records <- expected_field_design %>%
  anti_join(
    microresp_sample_table %>% select(Sampling_Date, Plot_ID),
    by = c("Sampling_Date", "Plot_ID")
  )

# ---- Metatranscriptomes -----------------------------------------------------

transcript_columns <- names(
  read_csv(transcript_file, n_max = 0, show_col_types = FALSE)
)

transcript_sample_table <- tibble(
  Sample_mt = transcript_columns[-1],
  Sample_cat = str_remove(Sample_mt, "_MT$")
) %>%
  left_join(
    field_sample_table %>%
      select(Sample_cat, Sampling_Date, Plot_ID, Block, Season),
    by = "Sample_cat"
  ) %>%
  arrange(Sampling_Date, Plot_ID)

# ---- One-month lag pairs (Q3) ----------------------------------------------
# Transcript samples paired with MicroResp from the same plot at the next
# sampling date.

sampling_dates <- sort(unique(microresp_sample_table$Sampling_Date))

lag_pairs <- transcript_sample_table %>%
  filter(!is.na(Sampling_Date)) %>%
  mutate(
    next_date = sampling_dates[match(Sampling_Date, sampling_dates) + 1]
  ) %>%
  left_join(
    microresp_sample_table %>%
      select(next_date = Sampling_Date, Plot_ID, next_Sample_cat = Sample_cat),
    by = c("next_date", "Plot_ID")
  )

lag_pair_summary <- lag_pairs %>%
  group_by(transcript_date = Sampling_Date, next_date) %>%
  summarise(
    days_between = as.numeric(first(next_date) - first(transcript_date)),
    n_transcript_plots = n(),
    n_pairs = sum(!is.na(next_Sample_cat)),
    plots_without_next_microresp = paste(Plot_ID[is.na(next_Sample_cat)],
                                         collapse = "; "),
    .groups = "drop"
  )

write_csv(lag_pairs, file.path(result_dir, "lag_pairs_transcript_to_next_microresp.csv"))
write_csv(lag_pair_summary, file.path(result_dir, "lag_pair_summary.csv"))

# ---- Counts by date ---------------------------------------------------------

all_sampling_dates <- tibble(
  Sampling_Date = sort(
    unique(c(
      field_sample_table$Sampling_Date,
      microresp_sample_table$Sampling_Date,
      transcript_sample_table$Sampling_Date
    ))
  )
)

respiration_counts <- field_sample_table %>%
  filter(!is.na(Resp_rate_mg_C_m2_h)) %>%
  count(Sampling_Date, name = "total_field_soil_respiration")

root_counts <- field_sample_table %>%
  filter(!is.na(Roots_gdw)) %>%
  count(Sampling_Date, name = "root_biomass")

microresp_counts <- microresp_sample_table %>%
  count(Sampling_Date, name = "microresp")

transcript_counts <- transcript_sample_table %>%
  filter(!is.na(Sampling_Date)) %>%
  count(Sampling_Date, name = "metatranscriptomes")

matched_root_respiration_counts <- field_sample_table %>%
  filter(
    !is.na(Roots_gdw),
    !is.na(Resp_rate_mg_C_m2_h)
  ) %>%
  count(Sampling_Date, name = "matched_root_and_field_respiration")

matched_microresp_transcript_counts <- microresp_sample_table %>%
  inner_join(
    transcript_sample_table %>%
      filter(!is.na(Sampling_Date)) %>%
      select(Sample_cat, transcript_date = Sampling_Date),
    by = "Sample_cat"
  ) %>%
  filter(Sampling_Date == transcript_date) %>%
  count(Sampling_Date, name = "matched_microresp_and_transcripts")

availability_by_date <- all_sampling_dates %>%
  left_join(respiration_counts, by = "Sampling_Date") %>%
  left_join(root_counts, by = "Sampling_Date") %>%
  left_join(microresp_counts, by = "Sampling_Date") %>%
  left_join(transcript_counts, by = "Sampling_Date") %>%
  left_join(matched_root_respiration_counts, by = "Sampling_Date") %>%
  left_join(matched_microresp_transcript_counts, by = "Sampling_Date") %>%
  mutate(
    across(-Sampling_Date, ~ replace_na(.x, 0L)),
    Season = if_else(
      month(Sampling_Date) %in% 4:10,
      "Plant-active, Apr-Oct",
      "Plant-dormant, Nov-Mar"
    )
  ) %>%
  relocate(Season, .after = Sampling_Date)

availability_long <- availability_by_date %>%
  pivot_longer(
    cols = -c(Sampling_Date, Season),
    names_to = "data_type",
    values_to = "n_observations"
  )

availability_summary <- availability_long %>%
  group_by(data_type) %>%
  summarise(
    n_observations = sum(n_observations),
    n_dates_with_data = sum(n_observations > 0),
    first_date = min(Sampling_Date[n_observations > 0]),
    last_date = max(Sampling_Date[n_observations > 0]),
    .groups = "drop"
  )

# ---- Missing and unmatched samples -----------------------------------------

matching_issues <- bind_rows(
  missing_expected_field_records %>%
    transmute(
      issue_type = "Missing expected field record",
      Sample_cat = NA_character_,
      Sampling_Date,
      Plot_ID,
      Block,
      details = "No field record for an expected 100% irrigation plot-date."
    ),
  field_sample_table %>%
    filter(is.na(Resp_rate_mg_C_m2_h)) %>%
    transmute(
      issue_type = "Missing field respiration",
      Sample_cat,
      Sampling_Date,
      Plot_ID,
      Block,
      details = "Field record exists, but total soil respiration is missing."
    ),
  field_sample_table %>%
    filter(is.na(Roots_gdw)) %>%
    transmute(
      issue_type = "Missing root biomass",
      Sample_cat,
      Sampling_Date,
      Plot_ID,
      Block,
      details = "Field record exists, but root biomass is missing."
    ),
  missing_expected_microresp_records %>%
    transmute(
      issue_type = "Missing expected MicroResp sample",
      Sample_cat = NA_character_,
      Sampling_Date,
      Plot_ID,
      Block,
      details = "No 5-hour carbon-normalized MicroResp sample for an expected 100% irrigation plot-date."
    ),
  microresp_sample_table %>%
    anti_join(field_sample_table %>% select(Sample_cat), by = "Sample_cat") %>%
    transmute(
      issue_type = "MicroResp sample unmatched to field table",
      Sample_cat,
      Sampling_Date,
      Plot_ID,
      Block,
      details = "MicroResp sample has no matching Sample_cat in the cleaned field table."
    ),
  field_sample_table %>%
    anti_join(microresp_sample_table %>% select(Sample_cat), by = "Sample_cat") %>%
    transmute(
      issue_type = "Field sample unmatched to MicroResp",
      Sample_cat,
      Sampling_Date,
      Plot_ID,
      Block,
      details = "Field sample has no matching 5-hour carbon-normalized MicroResp sample."
    ),
  transcript_sample_table %>%
    filter(is.na(Sampling_Date)) %>%
    transmute(
      issue_type = "Transcript sample unmatched to field metadata",
      Sample_cat,
      Sampling_Date,
      Plot_ID,
      Block,
      details = "Transcript sample name could not be matched to the cleaned field table."
    ),
  transcript_sample_table %>%
    filter(!is.na(Sampling_Date)) %>%
    anti_join(microresp_sample_table %>% select(Sample_cat), by = "Sample_cat") %>%
    transmute(
      issue_type = "Transcript sample unmatched to MicroResp",
      Sample_cat,
      Sampling_Date,
      Plot_ID,
      Block,
      details = "Transcript sample has no matching 5-hour carbon-normalized MicroResp sample."
    ),
  microresp_sample_table %>%
    anti_join(transcript_sample_table %>% select(Sample_cat), by = "Sample_cat") %>%
    transmute(
      issue_type = "MicroResp sample without transcript",
      Sample_cat,
      Sampling_Date,
      Plot_ID,
      Block,
      details = "No metatranscriptome was collected for this MicroResp sample; this is reported as availability, not assumed to be an error."
    )
) %>%
  arrange(issue_type, Sampling_Date, Plot_ID)

design_summary <- tibble(
  item = c(
    "Blocks in repository metadata",
    "Plots across both treatments",
    "Plots in 100% treatment",
    "Expected monthly sampling points",
    "Expected 100% plot-date records",
    "Observed 100% field records",
    "Observed 100% MicroResp samples",
    "Observed 100% transcript samples",
    "Observed transcript sampling dates",
    "Field duplicate sample records",
    "MicroResp duplicate sample-substrate measurements"
  ),
  value = c(
    n_distinct(field_100$Block),
    12,
    n_distinct(field_sample_table$Plot_ID),
    n_distinct(field_sample_table$Sampling_Date),
    nrow(expected_field_design),
    nrow(field_sample_table),
    nrow(microresp_sample_table),
    nrow(transcript_sample_table),
    n_distinct(transcript_sample_table$Sampling_Date, na.rm = TRUE),
    nrow(field_duplicate_records),
    nrow(microresp_duplicate_measurements)
  )
)

write_csv(
  availability_by_date,
  file.path(result_dir, "sample_availability_by_date.csv")
)
write_csv(
  availability_long,
  file.path(result_dir, "sample_availability_long.csv")
)
write_csv(
  availability_summary,
  file.path(result_dir, "sample_availability_summary.csv")
)
write_csv(
  matching_issues,
  file.path(result_dir, "missing_and_unmatched_samples.csv")
)
write_csv(
  design_summary,
  file.path(result_dir, "design_summary.csv")
)
write_csv(
  field_duplicate_records,
  file.path(result_dir, "field_duplicate_records.csv")
)
write_csv(
  microresp_duplicate_measurements,
  file.path(result_dir, "microresp_duplicate_measurements.csv")
)

cat("\nSample availability audit complete.\n")
print(availability_by_date)
cat("\nDesign summary:\n")
print(design_summary)
cat("\nMatching issues by type:\n")
print(count(matching_issues, issue_type, sort = TRUE))
