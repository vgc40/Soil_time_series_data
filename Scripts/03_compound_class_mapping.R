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

# Availability (H1b), decision D1 revised 2026-10-25. One rule: a substrate is
# "complex" when microbes must break it apart outside the cell (a polymer, or
# an ester that must be hydrolysed) before the pieces can be taken up;
# otherwise it is a small molecule taken up directly ("readily available").
complex_residue <- c(
  "Cellulose", "Hemicellulose", "Lignin", "Chitin", "Chitosan", "BSA",
  "Methyl_oleate", "Tween40"
)

# Sensitivity: only typical root-exudate compounds count as readily available
# (sugars, organic acids, amino acids); other small molecules are excluded.
amino_acids <- c(
  "Glutamic_acid", "L.Arginine", "L.Glutamine", "Glycine", "L.Alanine",
  "L.Leucine", "GABA"
)
root_exudate_types <- c(
  "D.Glucose", "D.Fructose", "D.Trehalose", "D.Xylose", "D.Cellobiose",
  "Citric_acid", "Lactic_acid", "Malic_acid", "Oxalic_acid", "Succinic_acid",
  amino_acids
)

# Compound classes. The first four are the classes linked to transcripts
# (Q2, Q3, decision D2); the rest extend the same idea to every substrate so
# SIR can be summarised by chemistry as well as by availability.
substrate_groups <- list(
  "Lipid-like" = c("Methyl_oleate", "Tween40"),
  "Aromatic" = c("Lignin", "Salicylic_acid", "Vanillic_acid", "Phenylethylamine"),
  "Amino sugar" = c("Chitin", "Chitosan", "Glucosamine_hydrochloride",
                    "N.Acetyl.D.glucosamine"),
  "Amino acid" = c(amino_acids, "BSA"),
  "Carbohydrate" = c("Cellulose", "Hemicellulose", "D.Cellobiose", "D.Glucose",
                     "D.Fructose", "D.Trehalose", "D.Xylose"),
  "Organic acid" = c("Citric_acid", "Lactic_acid", "Malic_acid", "Oxalic_acid",
                     "Succinic_acid"),
  "Aliphatic amine" = "Diethylamine"
)
compound_classes <- substrate_groups[c("Lipid-like", "Aromatic", "Amino sugar", "Amino acid")]

# Within-class H1b check: the polymer and its building blocks in one class.
# Lipid-like has no monomer in the assay (oleic acid was not tested).
polymer_monomer_classes <- c("Aromatic", "Amino sugar", "Amino acid", "Carbohydrate")

substrate_notes <- c(
  "Cellulose" = "Plant cell-wall polysaccharide (glucose polymer).",
  "Hemicellulose" = "Plant cell-wall polysaccharide (mainly xylose polymer).",
  "Lignin" = "Plant cell-wall aromatic polymer.",
  "Chitin" = "Fungal cell-wall and arthropod polymer of N-acetylglucosamine.",
  "Chitosan" = "Deacetylated chitin (glucosamine polymer).",
  "BSA" = "Protein (bovine serum albumin); must be cut into amino acids/peptides.",
  "Methyl_oleate" = "Fatty acid methyl ester; stands in for plant wax esters.",
  "Tween40" = "Polysorbate fatty acid ester; stands in for lipid esters.",
  "D.Cellobiose" = "Disaccharide released from cellulose.",
  "D.Glucose" = "Simple sugar; common root exudate.",
  "D.Fructose" = "Simple sugar; common root exudate.",
  "D.Trehalose" = "Disaccharide; microbial and plant stress sugar.",
  "D.Xylose" = "Simple sugar released from hemicellulose.",
  "Citric_acid" = "Organic acid; common root exudate.",
  "Lactic_acid" = "Organic acid; microbial fermentation product.",
  "Malic_acid" = "Organic acid; common root exudate.",
  "Oxalic_acid" = "Organic acid; root and fungal exudate.",
  "Succinic_acid" = "Organic acid; central metabolite and root exudate.",
  "Glutamic_acid" = "Amino acid; common root exudate.",
  "L.Arginine" = "Amino acid.",
  "L.Glutamine" = "Amino acid.",
  "Glycine" = "Amino acid.",
  "L.Alanine" = "Amino acid.",
  "L.Leucine" = "Amino acid.",
  "GABA" = "Non-protein amino acid (gamma-aminobutyric acid); signalling molecule.",
  "Glucosamine_hydrochloride" = "Amino sugar monomer (building block of chitosan).",
  "N.Acetyl.D.glucosamine" = "Amino sugar monomer (building block of chitin).",
  "Salicylic_acid" = "Plant phenolic acid; signalling compound.",
  "Vanillic_acid" = "Phenolic acid released from lignin.",
  "Phenylethylamine" = "Aromatic amine.",
  "Diethylamine" = "Small aliphatic amine."
)

pathways <- read_csv(input_files$microresp_pathways, show_col_types = FALSE)

substrates <- pathways %>%
  distinct(Substrate, Substrate_category, Table_group)

compound_lookup <- enframe(compound_classes, name = "compound_class",
                           value = "Substrate") %>%
  unnest(Substrate)
group_lookup <- enframe(substrate_groups, name = "substrate_group",
                        value = "Substrate") %>%
  unnest(Substrate)

stopifnot(!anyDuplicated(group_lookup$Substrate))
stopifnot(setequal(group_lookup$Substrate, substrates$Substrate))
stopifnot(all(c(complex_residue, names(substrate_notes)) %in% substrates$Substrate))

substrate_classes <- substrates %>%
  left_join(compound_lookup, by = "Substrate") %>%
  left_join(group_lookup, by = "Substrate") %>%
  mutate(
    availability_class = if_else(Substrate %in% complex_residue,
                                 "Complex", "Readily available"),
    availability_class_sensitivity = case_when(
      Substrate %in% complex_residue ~ "Complex",
      Substrate %in% root_exudate_types ~ "Readily available",
      TRUE ~ "Excluded"
    ),
    availability_rationale = if_else(
      availability_class == "Complex",
      "Polymer or ester: must be broken apart outside the cell before uptake.",
      "Small molecule: can be taken up directly."
    ),
    polymer_monomer_class = if_else(substrate_group %in% polymer_monomer_classes,
                                    substrate_group, NA_character_),
    substrate_note = substrate_notes[Substrate],
    display_name = Substrate %>%
      str_replace_all("_", " ") %>%
      str_replace("^D\\.", "D-") %>%
      str_replace("^L\\.", "L-") %>%
      str_replace("N\\.Acetyl\\.D\\.glucosamine", "N-acetyl-D-glucosamine") %>%
      str_replace("Glucosamine hydrochloride", "Glucosamine-HCl"),
    source = "Decisions D1 (revised 2026-10-25) and D2; Table_group from Data/Microresp_pathways.csv"
  ) %>%
  arrange(availability_class, substrate_group, Substrate)

stopifnot(sum(substrate_classes$availability_class == "Complex") == 8)

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

# Plain-language summary of each KEGG pathway and why it is linked to the
# class's substrates. KEGG pathways hold genes for both breakdown and
# synthesis, so expression indicates activity of the pathway, not use of the
# assay substrate itself.
pathway_descriptions <- tribble(
  ~brite_C_ID, ~what_it_does, ~why_linked,
  "00071", "Beta-oxidation: breaks fatty acids down to acetyl-CoA for energy.", "Fatty acids are released when methyl oleate and Tween 40 esters are split.",
  "00073", "Builds and modifies plant cuticle and root-barrier lipids (cutin, suberin, waxes).", "Methyl oleate stands in for plant wax esters.",
  "00561", "Makes and breaks down glycerol-based lipids (e.g., triacylglycerols) with lipases.", "Lipases and esterases split lipid esters such as methyl oleate and Tween 40.",
  "00564", "Makes and turns over membrane phospholipids.", "Tween 40 is a fatty acid ester of a polyol, handled by related ester-cleaving enzymes.",
  "00350", "Breaks down tyrosine and related phenolic rings (e.g., via homogentisate).", "Shared ring-cleavage steps for lignin-derived and amine-linked aromatics.",
  "00360", "Breaks down phenylalanine and phenylethylamine (phenylacetate route).", "Direct route for phenylethylamine; shared steps with lignin-derived aromatics.",
  "00362", "Opens and degrades benzene rings of benzoate and related aromatics.", "Central route for aromatic acids such as vanillic acid and lignin fragments.",
  "00626", "Degrades fused-ring aromatics; salicylate is a key intermediate.", "Salicylic acid is processed through this pathway.",
  "00627", "Degrades aminobenzoates and related substituted aromatics.", "Shared steps for phenolic acids and lignin fragments.",
  "00940", "Makes and modifies phenylpropanoids, the building blocks of lignin.", "Linked to lignin and its phenylpropanoid units.",
  "00520", "Converts amino sugars (glucosamine, N-acetylglucosamine) and builds cell-wall sugars; includes chitin and chitosan breakdown products.", "Chitin, chitosan, glucosamine and NAG all enter metabolism here.",
  "00220", "Makes arginine and links it to glutamate and urea-cycle intermediates.", "Arginine and glutamate metabolism.",
  "00250", "Interconverts alanine, aspartate and glutamate; includes the GABA shunt.", "Alanine, glutamate, glutamine and GABA; amino acids released from protein.",
  "00260", "Interconverts glycine, serine and threonine.", "Glycine; amino acids released from protein.",
  "00270", "Metabolises the sulfur amino acids cysteine and methionine.", "Amino acids released from protein (BSA).",
  "00280", "Breaks down the branched-chain amino acids valine, leucine and isoleucine.", "Leucine; amino acids released from protein.",
  "00330", "Interconverts arginine, proline and glutamate; includes polyamines.", "Arginine, glutamate and GABA; amino acids released from protein.",
  "00410", "Metabolises beta-alanine and related compounds.", "Linked to GABA metabolism.",
  "00650", "Metabolises butanoate (butyrate); the GABA shunt feeds succinate through here.", "GABA is broken down through the GABA shunt.",
  "00670", "Moves one-carbon units using folate; glycine cleavage feeds this pool.", "Glycine is cleaved into one-carbon units here."
)

unique_brite_c <- unique_brite_c %>%
  left_join(pathway_descriptions, by = "brite_C_ID")
stopifnot(!anyNA(unique_brite_c$what_it_does))

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
