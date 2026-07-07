suppressPackageStartupMessages({
  library(tidyverse)
})
rm(list = ls())

# -----------------------------
# File paths
# -----------------------------
expr_file  <- "C:/Users/gara009/OneDrive - PNNL/Soil Microbiome SFA - A3.1 Rhizosphere and Moisture Effects/Prosser_Monthly/Manuscripts/A3-1-1_Prosser_time_biogeochem/Transcriptomic_abundance/geTMM_normalized_expression_noZeroes.csv"
meta_file  <- "Data/Data_merged_cleaned.csv"
annotation_file <- "C:/Users/gara009/OneDrive - PNNL/Soil Microbiome SFA - A3.1 Rhizosphere and Moisture Effects/Prosser_Monthly/Manuscripts/A3-1-1_Prosser_time_biogeochem/Annotations/2.5kbSizeOut.emapper.annotations.kegg_annotated.csv"

# -----------------------------
# Read files
# -----------------------------
expr  <- read.csv(expr_file, check.names = FALSE)
names(expr)[1] <- "Query"  # rename first column 
meta  <- read.csv(meta_file, check.names = FALSE)
annotation <- read.csv(annotation_file, check.names = FALSE)

samples = colnames(expr)[-1]  # get sample names from expression data
#samples = gsub("_MT", "", samples)  # remove "MT" from sample names

# merge annotations with expression 

data = merge(expr, annotation, by = "Query",  all = TRUE)

# Now aggregate the data by brite level

# Group by brite D and then add the group rows for the column samples and produce a final dataframe of brite d per sample 
brite_D = data %>%
  group_by(brite_D) %>%
  summarise(across(all_of(samples), sum, na.rm = TRUE)) %>%
  ungroup()
write.csv(brite_D, "Data/brite_D_expression.csv", row.names = FALSE)

# Group by brite C 
brite_C = data %>%
  group_by(brite_C) %>%
  summarise(across(all_of(samples), sum, na.rm = TRUE)) %>%
  ungroup()
write.csv(brite_C, "Data/brite_C_expression.csv", row.names = FALSE)


# Group by gene "Query"
gene = data %>%
  group_by(Query) %>%
  summarise(across(all_of(samples), sum, na.rm = TRUE)) %>%
  ungroup()
write.csv(gene, "Data/Query_expression.csv", row.names = FALSE)

