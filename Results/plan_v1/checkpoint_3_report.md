# Checkpoint 3 report: draft KO gene sets (plan §4)

Script: `Scripts/plan_v1/05_ko_gene_sets.R`. Source: KEGG names in
`Data/brite_D_expression.csv` (9,022 KOs, 35 samples), starting from the
curated chitin/chitosan list. **Draft for approval. No Q2 analysis has been
run.**

## Files

| File | Contents |
|---|---|
| `KO_gene_sets_DRAFT.csv` | One row per set × KO: set, set type, group, pathway step (extracellular / intracellular), KO, KEGG name, step, reason for inclusion, confidence, other sets sharing the KO |
| `KO_gene_sets_DRAFT_expression.csv` | Summed expression per set across the 35 samples (all KOs, and high-confidence KOs only) |
| `KO_gene_sets_DRAFT_ko_expression.csv` | Mean expression and number of nonzero samples for every KO in the draft |
| `KO_gene_sets_DRAFT_notes.csv` | Compounds proposed as "no specific gene set" |

Expression values are geTMM (gene-length-corrected TMM counts per million)
summed per KO; see the deviation log, resolution of O2. Values below are the
file's values summed over the KOs in a set. "Share" is the set's sum as parts
per million of the sample's total over all KOs.

## How the draft was built

- **Polymers (extracellular):** enzymes that cleave the polymer. Chitin and
  chitosan start from the curated list. `nagZ` (K01207) from that list is
  placed in the NAG enzyme set, not the chitin set, because it acts on
  chitobiose and oligomers rather than the polymer. The NagE PTS, NagK, NagA,
  GspK, and NagB entries from the curated list go to the NAG and glucosamine
  sets.
- **Simple compounds, NAG, glucosamine (intracellular):** a compound-specific
  transporter plus the first committed catabolic step. Glycolysis, the TCA
  cycle, the pentose phosphate pathway, and general amino acid metabolism are
  excluded. For example, malate dehydrogenase and succinate dehydrogenase are
  out, as are branched-chain aminotransferase and the shared BCKDH complex
  for leucine. Biosynthetic enzymes (trehalose synthases, glutamine
  synthetase, NADP-glutamate dehydrogenase) are out.
- **Enzyme sets:** KOs carrying the EC number of each assay (CB 3.2.1.91,
  BG 3.2.1.21, BX 3.2.1.37, NAG 3.2.1.52, Leu 3.4.11.1, Ala 3.4.11.2 and
  3.4.11.14, AP 3.1.3.2). All are labelled extracellular, as the plan says,
  although some (BglB, PepN, NagZ) are mainly intracellular in bacteria.
- **Confidence:** high = specific to the compound and step; medium = the
  right step but shared with a close relative, or of uncertain direction;
  low = a candidate that is likely shared or indirect.
- **Glucose** has no set; the plan does not request one.

## Summary of the draft sets

| Set | Group | Step | KOs (high) | Mean sum | Share (ppm) | Samples with sum 0 | Mean sum, high only |
|---|---|---|---|---|---|---|---|
| Cellulose | Polymer | extracellular | 3 (2) | 601 | 758 | 0 | 598 |
| Hemicellulose | Polymer | extracellular | 8 (2) | 200 | 252 | 0 | 122 |
| Lignin | Polymer | extracellular | 2 (0) | **2.8** | 3.2 | **13** | — |
| Chitin | Polymer | extracellular | 7 (4) | 258 | 309 | 0 | 164 |
| Chitosan | Polymer | extracellular | 2 (1) | **1.6** | 2.1 | **20** | 1.4 |
| D.Fructose | Simple | intracellular | 11 (7) | 1,526 | 1,676 | 0 | 1,468 |
| D.Xylose | Simple | intracellular | 6 (5) | 1,071 | 1,312 | 0 | 1,053 |
| D.Trehalose | Simple | intracellular | 10 (7) | 111 | 136 | 0 | 24 |
| D.Cellobiose | Simple | intracellular | 9 (4) | 246 | 295 | 0 | 234 |
| Citric_acid | Simple | intracellular | 8 (5) | 448 | 528 | 0 | 41 |
| Malic_acid | Simple | intracellular | 4 (1) | 462 | 538 | 0 | **0** |
| Oxalic_acid | Simple | intracellular | 5 (4) | 65 | 76 | 0 | 60 |
| Succinic_acid | Simple | intracellular | 5 (0) | 232 | 268 | 1 | — |
| Lactic_acid | Simple | intracellular | 9 (7) | 65 | 74 | 0 | 41 |
| Glutamic_acid | Simple | intracellular | 13 (6) | 1,323 | 1,556 | 0 | 166 |
| L.Arginine | Simple | intracellular | 13 (11) | 381 | 437 | 0 | 322 |
| L.Glutamine | Simple | intracellular | 9 (4) | 956 | 1,183 | 0 | 860 |
| Glycine | Simple | intracellular | 7 (5) | 1,035 | 1,191 | 0 | 1,029 |
| L.Alanine | Simple | intracellular | 5 (2) | 565 | 677 | 0 | 355 |
| L.Leucine | Simple | intracellular | 10 (3) | **11,635** | **14,459** | 0 | 23 |
| N.Acetyl.D.glucosamine | Q2 | intracellular | 7 (7) | 37 | 43 | 2 | 37 |
| Glucosamine_hydrochloride | Q2 | intracellular | 6 (1) | 67 | 80 | 0 | 0.7 |
| CB | Enzyme | extracellular | 1 (1) | 114 | 147 | 0 | 114 |
| BG | Enzyme | extracellular | 3 (2) | 261 | 311 | 0 | 219 |
| BX | Enzyme | extracellular | 3 (3) | **5.6** | 6.8 | **3** | 5.6 |
| NAG | Enzyme | extracellular | 3 (1) | 587 | 743 | 0 | 42 |
| Leu | Enzyme | extracellular | 3 (1) | 130 | 152 | 0 | 64 |
| Ala | Enzyme | extracellular | 2 (1) | 91 | 107 | 0 | 66 |
| AP | Enzyme | extracellular | 7 (3) | 171 | 207 | 0 | 33 |

20 of the 176 draft KOs are zero in all 35 samples (listed in
`KO_gene_sets_DRAFT_ko_expression.csv`), including CHIB (K20547), citA
(K03288), citF (K01643), maeN (K11616), and three of the acid phosphatases.

## Decisions needed

1. **Succinic acid.** Only shared C4-dicarboxylate transporters qualify; the
   first catabolic step is succinate dehydrogenase (TCA, excluded). Proposed:
   mark as **"no specific gene set"**.
2. **Lignin.** The file has no laccase, lignin peroxidase, or manganese
   peroxidase KOs. The draft has only a DyP-type peroxidase (medium) and EfeB
   (low). Mean sum 2.8; zero in 13 of 35 samples. Options: keep as a
   near-zero set, or mark "no specific gene set".
3. **Chitosan.** Chitosanase plus exo-glucosaminidase; mean sum 1.6, zero in
   20 of 35 samples. Same choice as lignin.
4. **BX enzyme.** Three beta-xylosidases; mean sum 5.6, zero in 3 samples.
   Endo-xylanases are in the hemicellulose set, not here.
5. **Leucine is dominated by one shared transporter.** `livK` (K01999,
   branched-chain amino acid binding protein, shared with Ile and Val) has a
   mean of 10,719, which is 92% of the set and about 1.3% of all expression.
   The leucine-specific steps (Ivd, MccA/B) sum to 23. Options: keep the
   liv transporter, drop it, or use the high-confidence KOs only.
6. **Other sets carried by medium or low KOs:**
   - Glutamic acid: `gltI` (glutamate/aspartate binding protein, medium)
     is 837 of 1,323.
   - Citric acid: the putative TctABC transporter (medium) is 407 of 448.
   - Malic acid: the only high-confidence KO (MaeN) is zero; the set is DctA
     (low, shared with succinate) and malic enzymes (medium).
   - NAG enzyme: the eukaryotic hexosaminidase HEXA_B (medium) is 545 of 587;
     bacterial NagZ is 42.
   - Trehalose: ThuEFG (medium, shared with maltose) is 86 of 111.
   - Glucosamine: the only high-confidence KO (GspK) averages 0.7; the set is
     mostly NagB (shared with the NAG pathway).
7. **AP enzyme.** The assay buffer is pH 7.91. Alkaline phosphatases PhoA and
   PhoD are included at low confidence. Keep or drop?
8. **KOs shared between sets:** K19668 (cellulose and CB), K01256 (Leu and
   Ala), K11103 (malate and succinate), K03310 and K11737 (glycine and
   alanine). Shared KOs make those sets' values correlated in H2b and H2c.
   Keep shared KOs, or assign each to one set?
9. **Extracellular label for enzyme sets.** The plan labels all assayed
   enzymes extracellular. Several KOs are mainly intracellular in bacteria
   (BglB, PepN, NagZ, LAP3). The draft follows the plan.

## Open items carried forward

- **O1:** resolved, 5–24 h rate.
- **O2:** resolved, geTMM summed per KO (generating code not found).
- **O3:** mg C per well for succinic acid, chitin, methyl oleate.

---

## Addendum: recommended gene sets (D8)

At VGC's request, Claude applied its recommended choice for each decision
above. Rules and reasons are in the deviation log (D8). **Still a proposal;
Q2 does not start until approved.**

### Choice for each decision

| # | Decision | Recommended choice |
|---|---|---|
| 1 | Succinic acid | No specific gene set |
| 2 | Lignin | Insufficient expression (near-zero rule); not used in Q2 |
| 3 | Chitosan | Insufficient expression (near-zero rule); not used in Q2 |
| 4 | BX enzyme | Keep (zero in 3 of 35 samples, under the rule's limit of 9); flagged as low expression |
| 5 | Leucine | Drop shared branched-chain transporters; keep leucine dehydrogenase, Ivd, MccA, MccB |
| 6 | Sets carried by medium or low KOs | Drop all low-confidence KOs; keep medium |
| 7 | AP alkaline phosphatases | Dropped (low confidence, by rule 1) |
| 8 | Shared KOs | Keep in each set they belong to |
| 9 | Extracellular label for enzymes | Keep, as the plan states |

### Recommended sets

| Set | Group | Step | KOs (high) | Mean sum | Samples with sum 0 | Status |
|---|---|---|---|---|---|---|
| Cellulose | Polymer | extracellular | 2 (2) | 598 | 0 | use in Q2 |
| Hemicellulose | Polymer | extracellular | 8 (2) | 200 | 0 | use in Q2 |
| Chitin | Polymer | extracellular | 6 (4) | 250 | 0 | use in Q2 |
| Lignin | Polymer | extracellular | 1 (0) | 0.1 | 33 | insufficient expression |
| Chitosan | Polymer | extracellular | 2 (1) | 1.6 | 20 | insufficient expression |
| D.Fructose | Simple | intracellular | 11 (7) | 1,526 | 0 | use in Q2 |
| D.Xylose | Simple | intracellular | 6 (5) | 1,071 | 0 | use in Q2 |
| D.Trehalose | Simple | intracellular | 10 (7) | 111 | 0 | use in Q2 |
| D.Cellobiose | Simple | intracellular | 9 (4) | 246 | 0 | use in Q2 |
| Citric_acid | Simple | intracellular | 8 (5) | 448 | 0 | use in Q2 |
| Malic_acid | Simple | intracellular | 3 (1) | 243 | 0 | use in Q2 |
| Oxalic_acid | Simple | intracellular | 5 (4) | 65 | 0 | use in Q2 |
| Succinic_acid | Simple | intracellular | 0 | — | — | no specific gene set |
| Lactic_acid | Simple | intracellular | 9 (7) | 65 | 0 | use in Q2 |
| Glutamic_acid | Simple | intracellular | 10 (6) | 1,132 | 0 | use in Q2 |
| L.Arginine | Simple | intracellular | 13 (11) | 381 | 0 | use in Q2 |
| L.Glutamine | Simple | intracellular | 9 (4) | 956 | 0 | use in Q2 |
| Glycine | Simple | intracellular | 7 (5) | 1,035 | 0 | use in Q2 |
| L.Alanine | Simple | intracellular | 3 (2) | 361 | 0 | use in Q2 |
| L.Leucine | Simple | intracellular | 4 (3) | 87 | 0 | use in Q2 |
| N.Acetyl.D.glucosamine | Q2 | intracellular | 7 (7) | 37 | 2 | use in Q2 |
| Glucosamine_hydrochloride | Q2 | intracellular | 2 (1) | 4.3 | 3 | use in Q2 |
| AP | Enzyme | extracellular | 5 (3) | 36 | 0 | use in Q2 |
| Ala | Enzyme | extracellular | 2 (1) | 91 | 0 | use in Q2 |
| BG | Enzyme | extracellular | 3 (2) | 261 | 0 | use in Q2 |
| BX | Enzyme | extracellular | 3 (3) | 5.6 | 3 | use in Q2 |
| CB | Enzyme | extracellular | 1 (1) | 114 | 0 | use in Q2 |
| Leu | Enzyme | extracellular | 3 (1) | 130 | 0 | use in Q2 |
| NAG | Enzyme | extracellular | 3 (1) | 587 | 0 | use in Q2 |

### What this means for Q2 (counts only)

- **Compound sets used in Q2:** 19 (3 polymers, 14 simple compounds, NAG,
  glucosamine). These are the compounds in H2b and the SIR sets in H2c.
- **Pathway-step groups for H2c (test 5):** extracellular = 3 polymers + 7
  enzymes = 10 sets; intracellular = 16 sets.
- **H2a polymer expression score** rests on 3 sets (cellulose,
  hemicellulose, chitin).
- **Lowest-expression sets kept:** BX (mean 5.6), glucosamine (4.3), NAG
  compound (37); their centred log ratios will lean on the zero replacement
  in §5.1.
- AP drops from a mean of 171 to 36 without the alkaline phosphatases.
