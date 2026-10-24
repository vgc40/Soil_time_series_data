# =============================================================================
# 03 A priori substrate classes and uniquely mapped BRITE C categories
# =============================================================================
#
# Writes the two classification tables every figure script reads:
#
#   Data/compound_class_mapping.csv
#     One row per MicroResp substrate with
#       - availability_class (H1b): Readily available / Complex / Other
#       - compound_class (Q2, Q3): Lipid-like / Aromatic / Amino sugar /
#         Amino acid, or NA when the substrate is not in one of these classes
#
#   Data/compound_class_unique_BRITE_C.csv
#     BRITE C categories that Data/Microresp_pathways.csv links to substrates
#     of exactly one compound class.
#
# Class assignments are decisions D1 and D2 (2026-10-23) and are fixed here,
# before looking at any seasonal result. Data/Microresp_pathways.csv is the
# a priori substrate-to-BRITE C mapping (KEGG pathways containing the
# catabolic or assimilatory reactions for each substrate).

rm(list = ls())
source(file.path("Scripts", "helpers.R"))

# -----------------------------------------------------------------------------
# Substrate classes
# -----------------------------------------------------------------------------

readily_available <- c(
  "D.Glucose", "D.Fructose", "D.Trehalose", "D.Xylose", "D.Cellobiose",
  "Citric_acid", "Lactic_acid", "Malic_acid", "Oxalic_acid", "Succinic_acid"
)

complex_residue <- c(
  "Cellulose", "Hemicellulose", "Lignin", "Chitin", "Chitosan", "BSA",
  "Methyl_oleate"
)

amino_acids <- c(
  "Glutamic_acid", "L.Arginine", "L.Glutamine", "Glycine", "L.Alanine",
  "L.Leucine"
)

compound_classes <- list(
  "Lipid-like" = c("Methyl_oleate", "Tween40"),
  "Aromatic" = c("Salicylic_acid", "Vanillic_acid", "Lignin", "Phenylethylamine"),
  "Amino sugar" = c("Chitin", "Chitosan", "Glucosamine_hydrochloride",
                    "N.Acetyl.D.glucosamine"),
  "Amino acid" = c(amino_acids, "BSA", "GABA")
)

availability_rationale <- c(
  "Readily available" = "Low-molecular-weight sugar or organic acid typical of root exudates.",
  "Complex" = "Polymer or ester associated with plant or microbial residues.",
  "Other" = "Not assigned a priori to either H1b group."
)

pathways <- read_csv(input_files$microresp_pathways, show_col_types = FALSE)

substrates <- pathways %>%
  distinct(Substrate, Substrate_category, Table_group)

compound_lookup <- enframe(compound_classes, name = "compound_class",
                           value = "Substrate") %>%
  unnest(Substrate)

stopifnot(!anyDuplicated(compound_lookup$Substrate))
stopifnot(all(c(readily_available, complex_residue, compound_lookup$Substrate) %in%
                substrates$Substrate))

substrate_classes <- substrates %>%
  left_join(compound_lookup, by = "Substrate") %>%
  mutate(
    availability_class = case_when(
      Substrate %in% readily_available ~ "Readily available",
      Substrate %in% complex_residue ~ "Complex",
      TRUE ~ "Other"
    ),
    availability_class_sensitivity = if_else(
      Substrate %in% amino_acids, "Readily available", availability_class
    ),
    availability_rationale = availability_rationale[availability_class],
    display_name = Substrate %>%
      str_replace_all("_", " ") %>%
      str_replace("^D\\.", "D-") %>%
      str_replace("^L\\.", "L-") %>%
      str_replace("N\\.Acetyl\\.D\\.glucosamine", "N-acetyl-D-glucosamine") %>%
      str_replace("Glucosamine hydrochloride", "Glucosamine-HCl"),
    source = "Decisions D1/D2 (2026-10-23); Table_group from Data/Microresp_pathways.csv"
  ) %>%
  arrange(availability_class, compound_class, Substrate)

stopifnot(nrow(substrate_classes) == 31)

# -----------------------------------------------------------------------------
# BRITE C categories unique to one compound class
# -----------------------------------------------------------------------------

brite_links <- pathways %>%
  distinct(Substrate, brite_C) %>%
  left_join(compound_lookup, by = "Substrate")

# A category is unique when every substrate it is linked to belongs to the
# same compound class (substrates outside the four classes count as another
# class, so sharing with them disqualifies a category).
unique_brite_c <- brite_links %>%
  mutate(class_or_other = coalesce(compound_class, "Outside the four classes")) %>%
  group_by(brite_C) %>%
  summarise(
    n_classes = n_distinct(class_or_other),
    compound_class = first(class_or_other),
    linked_substrates = paste(sort(unique(Substrate)), collapse = "; "),
    .groups = "drop"
  ) %>%
  filter(n_classes == 1, compound_class != "Outside the four classes") %>%
  mutate(
    brite_C_ID = str_extract(brite_C, "^\\d{5}"),
    brite_C_name = str_remove(brite_C, "^\\d{5}\\s+") %>%
      str_remove("\\s*\\[PATH:ko\\d{5}\\]"),
    compound_class = factor(compound_class, levels = names(compound_classes))
  ) %>%
  select(compound_class, brite_C_ID, brite_C_name, brite_C, linked_substrates) %>%
  arrange(compound_class, brite_C_ID)

expected_counts <- c("Lipid-like" = 4, "Aromatic" = 6, "Amino sugar" = 1,
                     "Amino acid" = 9)
observed_counts <- table(unique_brite_c$compound_class)[names(expected_counts)]

if (!isTRUE(all(as.integer(observed_counts) == expected_counts))) {
  print(unique_brite_c)
  stop("Unique BRITE C counts differ from the 4/6/1/9 agreed in decision D2.")
}

write_csv(substrate_classes, input_files$substrate_classes)
write_csv(unique_brite_c, input_files$unique_brite_c)

cat("Substrate classes written for", nrow(substrate_classes), "substrates\n")
print(count(substrate_classes, availability_class))
cat("\nUnique BRITE C categories:", nrow(unique_brite_c), "\n")
print(count(unique_brite_c, compound_class))
