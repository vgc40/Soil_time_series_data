# =============================================================================
# 05. Draft KO gene sets (plan §4, Checkpoint 3)
# =============================================================================
#
# Drafts one gene set per compound (polymers, simple compounds, NAG,
# glucosamine) and per assayed enzyme from the KEGG names in
# brite_D_expression.csv, following the inclusion rule in plan §4. Writes the
# draft table and the summed expression of each set. No Q2 analysis is run
# until the sets are approved.

source(file.path("Scripts", "plan_v1", "00_setup.R"))

# -----------------------------------------------------------------------------
# Draft gene sets
# -----------------------------------------------------------------------------
# step: what the KO does in the set. reason: why it is included. Confidence:
# high = specific to the compound and the step; medium = specific step but
# shared with a close relative or of uncertain direction; low = included as a
# candidate, likely shared or indirect.

ext <- "extracellular"
int <- "intracellular"

draft <- tribble(
  ~set, ~set_type, ~group, ~pathway_step, ~KO, ~step, ~reason, ~confidence,
  # ---- Polymers ------------------------------------------------------------
  "Cellulose", "compound", "Polymer", ext, "K01179", "depolymerization", "endo-1,4-beta-glucanase cleaves cellulose chains", "high",
  "Cellulose", "compound", "Polymer", ext, "K19668", "depolymerization", "cellobiohydrolase releases cellobiose from cellulose ends (also in CB set)", "high",
  "Cellulose", "compound", "Polymer", ext, "K20542", "depolymerization", "endoglucanase BcsZ; encoded in bacterial cellulose synthesis operon, periplasmic", "low",
  "Hemicellulose", "compound", "Polymer", ext, "K01181", "depolymerization", "endo-1,4-beta-xylanase cleaves the xylan backbone", "high",
  "Hemicellulose", "compound", "Polymer", ext, "K15924", "depolymerization", "glucuronoarabinoxylan endo-1,4-beta-xylanase", "high",
  "Hemicellulose", "compound", "Polymer", ext, "K01218", "depolymerization", "endo-1,4-beta-mannanase; mannan-type hemicellulose", "medium",
  "Hemicellulose", "compound", "Polymer", ext, "K19355", "depolymerization", "endo-1,4-beta-mannanase (MAN)", "medium",
  "Hemicellulose", "compound", "Polymer", ext, "K01209", "side-chain removal", "alpha-L-arabinofuranosidase removes arabinose side chains of arabinoxylan", "medium",
  "Hemicellulose", "compound", "Polymer", ext, "K15921", "side-chain removal", "arabinoxylan arabinofuranohydrolase", "medium",
  "Hemicellulose", "compound", "Polymer", ext, "K20844", "side-chain removal", "non-reducing end alpha-L-arabinofuranosidase", "medium",
  "Hemicellulose", "compound", "Polymer", ext, "K05972", "side-chain removal", "acetylxylan esterase removes acetyl groups from xylan", "medium",
  "Lignin", "compound", "Polymer", ext, "K15733", "oxidative depolymerization", "dye-decolorizing (DyP-type) peroxidase; bacterial lignin oxidation", "medium",
  "Lignin", "compound", "Polymer", ext, "K16301", "oxidative depolymerization", "EfeB, DyP-type peroxidase; primary role in iron uptake", "low",
  "Chitin", "compound", "Polymer", ext, "K01183", "depolymerization", "chitinase (curated list)", "high",
  "Chitin", "compound", "Polymer", ext, "K13381", "depolymerization", "bifunctional chitinase/lysozyme ChiA (curated list)", "high",
  "Chitin", "compound", "Polymer", ext, "K20547", "depolymerization", "basic endochitinase B (curated list)", "high",
  "Chitin", "compound", "Polymer", ext, "K21606", "depolymerization", "endo-chitodextinase (curated list)", "high",
  "Chitin", "compound", "Polymer", ext, "K03791", "depolymerization", "putative chitinase (curated list)", "medium",
  "Chitin", "compound", "Polymer", ext, "K03933", "depolymerization", "chitin-binding protein CbpD, lytic polysaccharide monooxygenase (AA10)", "medium",
  "Chitin", "compound", "Polymer", ext, "K01452", "deacetylation", "chitin deacetylase converts chitin to chitosan; modifies rather than cleaves (curated list)", "low",
  "Chitosan", "compound", "Polymer", ext, "K01233", "depolymerization", "chitosanase (curated list)", "high",
  "Chitosan", "compound", "Polymer", ext, "K15855", "depolymerization", "exo-1,4-beta-D-glucosaminidase releases glucosamine from chitosan", "medium",
  # ---- Simple compounds ----------------------------------------------------
  "D.Fructose", "compound", "Simple", int, "K10552", "transport", "fructose ABC transporter, substrate-binding FrcB", "high",
  "D.Fructose", "compound", "Simple", int, "K10553", "transport", "fructose ABC transporter, permease FrcC", "high",
  "D.Fructose", "compound", "Simple", int, "K10554", "transport", "fructose ABC transporter, ATPase FrcA", "high",
  "D.Fructose", "compound", "Simple", int, "K02770", "transport", "fructose PTS EIIBC/EIIC FruA", "high",
  "D.Fructose", "compound", "Simple", int, "K02769", "transport", "fructose PTS EIIB", "high",
  "D.Fructose", "compound", "Simple", int, "K02768", "transport", "fructose PTS EIIA FruB", "high",
  "D.Fructose", "compound", "Simple", int, "K11183", "transport", "fructose PTS multiphosphoryl transfer protein FruB/Fpr", "medium",
  "D.Fructose", "compound", "Simple", int, "K11195", "transport", "fructose PTS EIIB LevE (lev family also takes mannose)", "medium",
  "D.Fructose", "compound", "Simple", int, "K02771", "transport", "fructose PTS EIID LevG (lev family also takes mannose)", "medium",
  "D.Fructose", "compound", "Simple", int, "K00882", "first committed step", "1-phosphofructokinase phosphorylates PTS-derived fructose-1-P", "high",
  "D.Fructose", "compound", "Simple", int, "K00847", "first committed step", "fructokinase; also acts on fructose from sucrose", "medium",
  "D.Xylose", "compound", "Simple", int, "K10543", "transport", "D-xylose ABC transporter, substrate-binding XylF", "high",
  "D.Xylose", "compound", "Simple", int, "K10544", "transport", "D-xylose ABC transporter, permease XylH", "high",
  "D.Xylose", "compound", "Simple", int, "K10545", "transport", "D-xylose ABC transporter, ATPase XylG", "high",
  "D.Xylose", "compound", "Simple", int, "K08138", "transport", "xylose:H+ symporter XylE", "high",
  "D.Xylose", "compound", "Simple", int, "K01805", "first committed step", "xylose isomerase XylA", "high",
  "D.Xylose", "compound", "Simple", int, "K00854", "second step", "xylulokinase XylB; xylulose can also come from other pentoses", "medium",
  "D.Trehalose", "compound", "Simple", int, "K02819", "transport", "trehalose PTS EIIBC TreB", "high",
  "D.Trehalose", "compound", "Simple", int, "K17311", "transport", "trehalose ABC transporter, substrate-binding TreS", "high",
  "D.Trehalose", "compound", "Simple", int, "K17312", "transport", "trehalose ABC transporter, permease TreT", "high",
  "D.Trehalose", "compound", "Simple", int, "K17314", "transport", "trehalose ABC transporter, ATPase TreV", "high",
  "D.Trehalose", "compound", "Simple", int, "K10236", "transport", "trehalose/maltose ABC transporter, substrate-binding ThuE; shared with maltose", "medium",
  "D.Trehalose", "compound", "Simple", int, "K10237", "transport", "trehalose/maltose ABC transporter, permease ThuF/SugA; shared with maltose", "medium",
  "D.Trehalose", "compound", "Simple", int, "K10238", "transport", "trehalose/maltose ABC transporter, permease ThuG/SugB; shared with maltose", "medium",
  "D.Trehalose", "compound", "Simple", int, "K01226", "first committed step", "trehalose-6-phosphate hydrolase TreC (after PTS uptake)", "high",
  "D.Trehalose", "compound", "Simple", int, "K01194", "first committed step", "alpha,alpha-trehalase TreA/TreF", "high",
  "D.Trehalose", "compound", "Simple", int, "K05342", "first committed step", "trehalose phosphorylase", "high",
  "D.Cellobiose", "compound", "Simple", int, "K10240", "transport", "cellobiose ABC transporter, substrate-binding CebE", "high",
  "D.Cellobiose", "compound", "Simple", int, "K10241", "transport", "cellobiose ABC transporter, permease CebF", "high",
  "D.Cellobiose", "compound", "Simple", int, "K10242", "transport", "cellobiose ABC transporter, permease CebG", "high",
  "D.Cellobiose", "compound", "Simple", int, "K02759", "transport", "cellobiose PTS EIIA CelC/ChbA; also chitobiose", "medium",
  "D.Cellobiose", "compound", "Simple", int, "K02760", "transport", "cellobiose PTS EIIB CelA/ChbB; also chitobiose", "medium",
  "D.Cellobiose", "compound", "Simple", int, "K02761", "transport", "cellobiose PTS EIIC CelB/ChbC; also chitobiose", "medium",
  "D.Cellobiose", "compound", "Simple", int, "K00702", "first committed step", "cellobiose phosphorylase", "high",
  "D.Cellobiose", "compound", "Simple", int, "K01222", "first committed step", "6-phospho-beta-glucosidase CelF (after PTS uptake); other aryl-glucosides too", "medium",
  "D.Cellobiose", "compound", "Simple", int, "K01223", "first committed step", "6-phospho-beta-glucosidase BglA (after PTS uptake); other aryl-glucosides too", "medium",
  "Citric_acid", "compound", "Simple", int, "K03288", "transport", "citrate/tricarballylate:H+ symporter CitA/TcuC", "high",
  "Citric_acid", "compound", "Simple", int, "K03300", "transport", "citrate-metal:H+ symporter, CitMHS family", "high",
  "Citric_acid", "compound", "Simple", int, "K09477", "transport", "citrate:succinate antiporter CitT", "high",
  "Citric_acid", "compound", "Simple", int, "K07793", "transport", "tripartite tricarboxylate transporter TctA (putative)", "medium",
  "Citric_acid", "compound", "Simple", int, "K07794", "transport", "tripartite tricarboxylate transporter TctB (putative)", "medium",
  "Citric_acid", "compound", "Simple", int, "K07795", "transport", "tripartite tricarboxylate transporter TctC (putative)", "medium",
  "Citric_acid", "compound", "Simple", int, "K01643", "first committed step", "citrate lyase alpha subunit CitF (citrate cleavage outside the TCA cycle)", "high",
  "Citric_acid", "compound", "Simple", int, "K01644", "first committed step", "citrate lyase beta subunit CitE", "high",
  "Malic_acid", "compound", "Simple", int, "K11616", "transport", "malate:Na+ symporter MaeN", "high",
  "Malic_acid", "compound", "Simple", int, "K11103", "transport", "aerobic C4-dicarboxylate transporter DctA; also succinate, fumarate (also in succinate set)", "low",
  "Malic_acid", "compound", "Simple", int, "K00027", "first committed step", "NAD malic enzyme MaeA; anaplerotic node, malate decarboxylation", "medium",
  "Malic_acid", "compound", "Simple", int, "K00029", "first committed step", "NADP malic enzyme MaeB; anaplerotic node, malate decarboxylation", "medium",
  "Oxalic_acid", "compound", "Simple", int, "K08177", "transport", "oxalate/formate antiporter OxlT", "high",
  "Oxalic_acid", "compound", "Simple", int, "K07749", "first committed step", "formyl-CoA transferase Frc activates oxalate to oxalyl-CoA", "high",
  "Oxalic_acid", "compound", "Simple", int, "K18702", "first committed step", "CoA:oxalate CoA-transferase UctC", "medium",
  "Oxalic_acid", "compound", "Simple", int, "K01577", "second step", "oxalyl-CoA decarboxylase Oxc", "high",
  "Oxalic_acid", "compound", "Simple", int, "K01569", "first committed step", "oxalate decarboxylase OxdD", "high",
  "Succinic_acid", "compound", "Simple", int, "K11103", "transport", "aerobic C4-dicarboxylate transporter DctA; also malate, fumarate", "low",
  "Succinic_acid", "compound", "Simple", int, "K07791", "transport", "anaerobic C4-dicarboxylate transporter DcuA", "low",
  "Succinic_acid", "compound", "Simple", int, "K07792", "transport", "anaerobic C4-dicarboxylate transporter DcuB", "low",
  "Succinic_acid", "compound", "Simple", int, "K03326", "transport", "C4-dicarboxylate transporter DcuC/DcuD", "low",
  "Succinic_acid", "compound", "Simple", int, "K11106", "transport", "L-tartrate/succinate antiporter TtdT", "low",
  "Lactic_acid", "compound", "Simple", int, "K00427", "transport", "L-lactate permease LldP/LctP", "high",
  "Lactic_acid", "compound", "Simple", int, "K03303", "transport", "lactate permease LutP/LctP", "high",
  "Lactic_acid", "compound", "Simple", int, "K08178", "transport", "lactate transporter JEN (fungal)", "medium",
  "Lactic_acid", "compound", "Simple", int, "K00101", "first committed step", "L-lactate dehydrogenase (cytochrome) LldD, oxidizes lactate to pyruvate", "high",
  "Lactic_acid", "compound", "Simple", int, "K18928", "first committed step", "L-lactate dehydrogenase complex LldE (LutA)", "high",
  "Lactic_acid", "compound", "Simple", int, "K18929", "first committed step", "L-lactate dehydrogenase complex LldF (LutB)", "high",
  "Lactic_acid", "compound", "Simple", int, "K00782", "first committed step", "L-lactate dehydrogenase complex LldG (LutC)", "high",
  "Lactic_acid", "compound", "Simple", int, "K03777", "first committed step", "D-lactate dehydrogenase (quinone) Dld", "high",
  "Lactic_acid", "compound", "Simple", int, "K00102", "first committed step", "D-lactate dehydrogenase (cytochrome) LDHD", "medium",
  "Glutamic_acid", "compound", "Simple", int, "K11102", "transport", "proton/glutamate symporter GltP", "high",
  "Glutamic_acid", "compound", "Simple", int, "K03312", "transport", "glutamate:Na+ symporter GltS", "high",
  "Glutamic_acid", "compound", "Simple", int, "K10005", "transport", "glutamate ABC transporter, substrate-binding GluB", "high",
  "Glutamic_acid", "compound", "Simple", int, "K10006", "transport", "glutamate ABC transporter, permease GluC", "high",
  "Glutamic_acid", "compound", "Simple", int, "K10007", "transport", "glutamate ABC transporter, permease GluD", "high",
  "Glutamic_acid", "compound", "Simple", int, "K10008", "transport", "glutamate ABC transporter, ATPase GluA", "high",
  "Glutamic_acid", "compound", "Simple", int, "K10001", "transport", "glutamate/aspartate ABC transporter GltI; also aspartate", "medium",
  "Glutamic_acid", "compound", "Simple", int, "K10002", "transport", "glutamate/aspartate ABC transporter GltK; also aspartate", "medium",
  "Glutamic_acid", "compound", "Simple", int, "K10003", "transport", "glutamate/aspartate ABC transporter GltJ; also aspartate", "medium",
  "Glutamic_acid", "compound", "Simple", int, "K10004", "transport", "glutamate/aspartate ABC transporter GltL; also aspartate", "medium",
  "Glutamic_acid", "compound", "Simple", int, "K00260", "first committed step", "catabolic NAD glutamate dehydrogenase RocG/GudB; GDH is also central N metabolism", "low",
  "Glutamic_acid", "compound", "Simple", int, "K15371", "first committed step", "NAD glutamate dehydrogenase GDH2; also central N metabolism", "low",
  "Glutamic_acid", "compound", "Simple", int, "K01580", "first committed step", "glutamate decarboxylase GadB (to GABA); acid-stress role", "low",
  "L.Arginine", "compound", "Simple", int, "K09996", "transport", "arginine ABC transporter, substrate-binding ArtJ", "high",
  "L.Arginine", "compound", "Simple", int, "K09997", "transport", "arginine ABC transporter, substrate-binding ArtI", "high",
  "L.Arginine", "compound", "Simple", int, "K09999", "transport", "arginine ABC transporter, permease ArtQ", "high",
  "L.Arginine", "compound", "Simple", int, "K10022", "transport", "arginine/ornithine ABC transporter, substrate-binding AotJ", "high",
  "L.Arginine", "compound", "Simple", int, "K10023", "transport", "arginine/ornithine ABC transporter, permease AotM", "high",
  "L.Arginine", "compound", "Simple", int, "K10024", "transport", "arginine/ornithine ABC transporter, permease AotQ", "high",
  "L.Arginine", "compound", "Simple", int, "K10025", "transport", "arginine/ornithine ABC transporter, ATPase AotP", "high",
  "L.Arginine", "compound", "Simple", int, "K02205", "transport", "arginine/ornithine permease RocE/RocC", "high",
  "L.Arginine", "compound", "Simple", int, "K10013", "transport", "lysine/arginine/ornithine ABC transporter ArgT; also lysine", "medium",
  "L.Arginine", "compound", "Simple", int, "K03758", "transport", "arginine:ornithine antiporter ArcD; also lysine permease", "medium",
  "L.Arginine", "compound", "Simple", int, "K01478", "first committed step", "arginine deiminase ArcA", "high",
  "L.Arginine", "compound", "Simple", int, "K01476", "first committed step", "arginase RocF", "high",
  "L.Arginine", "compound", "Simple", int, "K00673", "first committed step", "arginine N-succinyltransferase AstA", "high",
  "L.Glutamine", "compound", "Simple", int, "K10036", "transport", "glutamine ABC transporter, substrate-binding GlnH", "high",
  "L.Glutamine", "compound", "Simple", int, "K10037", "transport", "glutamine ABC transporter, permease GlnP", "high",
  "L.Glutamine", "compound", "Simple", int, "K10038", "transport", "glutamine ABC transporter, ATPase GlnQ", "high",
  "L.Glutamine", "compound", "Simple", int, "K11626", "transport", "putative sodium/glutamine symporter GlnT", "medium",
  "L.Glutamine", "compound", "Simple", int, "K10039", "transport", "aspartate/glutamate/glutamine ABC transporter Peb1A; also Asp, Glu", "medium",
  "L.Glutamine", "compound", "Simple", int, "K10040", "transport", "aspartate/glutamate/glutamine ABC transporter Peb1B; also Asp, Glu", "medium",
  "L.Glutamine", "compound", "Simple", int, "K10041", "transport", "aspartate/glutamate/glutamine ABC transporter Peb1C; also Asp, Glu", "medium",
  "L.Glutamine", "compound", "Simple", int, "K01425", "first committed step", "glutaminase GlsA", "high",
  "L.Glutamine", "compound", "Simple", int, "K05597", "first committed step", "glutaminase-asparaginase AnsB; also asparagine", "medium",
  "Glycine", "compound", "Simple", int, "K11737", "transport", "D-serine/D-alanine/glycine transporter CycA; also D-Ser, D-Ala (also in alanine set)", "medium",
  "Glycine", "compound", "Simple", int, "K03310", "transport", "alanine or glycine:cation symporter, AGCS family (also in alanine set)", "medium",
  "Glycine", "compound", "Simple", int, "K00281", "first committed step", "glycine cleavage system P protein GcvP", "high",
  "Glycine", "compound", "Simple", int, "K00282", "first committed step", "glycine cleavage system P protein subunit 1 GcvPA", "high",
  "Glycine", "compound", "Simple", int, "K00283", "first committed step", "glycine cleavage system P protein subunit 2 GcvPB", "high",
  "Glycine", "compound", "Simple", int, "K00605", "first committed step", "glycine cleavage system T protein GcvT", "high",
  "Glycine", "compound", "Simple", int, "K02437", "first committed step", "glycine cleavage system H protein GcvH", "high",
  "L.Alanine", "compound", "Simple", int, "K03310", "transport", "alanine or glycine:cation symporter, AGCS family (also in glycine set)", "medium",
  "L.Alanine", "compound", "Simple", int, "K11737", "transport", "CycA, D-alanine/D-serine/glycine transporter; L-alanine uptake uncertain (also in glycine set)", "low",
  "L.Alanine", "compound", "Simple", int, "K00259", "first committed step", "alanine dehydrogenase Ald", "high",
  "L.Alanine", "compound", "Simple", int, "K19244", "first committed step", "alanine dehydrogenase Ala", "high",
  "L.Alanine", "compound", "Simple", int, "K00285", "first committed step", "D-amino-acid dehydrogenase DadA; L-alanine route via racemase", "low",
  "L.Leucine", "compound", "Simple", int, "K01999", "transport", "branched-chain amino acid ABC transporter LivK; shared with Ile, Val", "medium",
  "L.Leucine", "compound", "Simple", int, "K01997", "transport", "branched-chain amino acid ABC transporter LivH; shared with Ile, Val", "medium",
  "L.Leucine", "compound", "Simple", int, "K01998", "transport", "branched-chain amino acid ABC transporter LivM; shared with Ile, Val", "medium",
  "L.Leucine", "compound", "Simple", int, "K01995", "transport", "branched-chain amino acid ABC transporter LivG; shared with Ile, Val", "medium",
  "L.Leucine", "compound", "Simple", int, "K01996", "transport", "branched-chain amino acid ABC transporter LivF; shared with Ile, Val", "medium",
  "L.Leucine", "compound", "Simple", int, "K03311", "transport", "branched-chain amino acid:cation transporter BrnQ (LIVCS); shared with Ile, Val", "medium",
  "L.Leucine", "compound", "Simple", int, "K00263", "first committed step", "leucine dehydrogenase; also Ile, Val", "medium",
  "L.Leucine", "compound", "Simple", int, "K00253", "first leucine-specific step", "isovaleryl-CoA dehydrogenase Ivd", "high",
  "L.Leucine", "compound", "Simple", int, "K01968", "leucine-specific step", "3-methylcrotonyl-CoA carboxylase alpha MccA", "high",
  "L.Leucine", "compound", "Simple", int, "K01969", "leucine-specific step", "3-methylcrotonyl-CoA carboxylase beta MccB", "high",
  # ---- Held for Q2: NAG and glucosamine ------------------------------------
  "N.Acetyl.D.glucosamine", "compound", "Q2", int, "K02803", "transport", "N-acetylglucosamine PTS EIIB NagE (curated list)", "high",
  "N.Acetyl.D.glucosamine", "compound", "Q2", int, "K02804", "transport", "N-acetylglucosamine PTS EIICB(A) NagE (curated list)", "high",
  "N.Acetyl.D.glucosamine", "compound", "Q2", int, "K10200", "transport", "N-acetylglucosamine ABC transporter, substrate-binding NgcE", "high",
  "N.Acetyl.D.glucosamine", "compound", "Q2", int, "K10201", "transport", "N-acetylglucosamine ABC transporter, permease NgcF", "high",
  "N.Acetyl.D.glucosamine", "compound", "Q2", int, "K10202", "transport", "N-acetylglucosamine ABC transporter, permease NgcG", "high",
  "N.Acetyl.D.glucosamine", "compound", "Q2", int, "K00884", "first committed step", "N-acetylglucosamine kinase NagK (curated list)", "high",
  "N.Acetyl.D.glucosamine", "compound", "Q2", int, "K01443", "second step", "N-acetylglucosamine-6-phosphate deacetylase NagA (curated list)", "high",
  "Glucosamine_hydrochloride", "compound", "Q2", int, "K18676", "first committed step", "glucosamine kinase GspK (curated list)", "high",
  "Glucosamine_hydrochloride", "compound", "Q2", int, "K02564", "second step", "glucosamine-6-phosphate deaminase NagB; also downstream of NAG (curated list)", "medium",
  "Glucosamine_hydrochloride", "compound", "Q2", int, "K02793", "transport", "mannose PTS EIIA ManX; carries glucosamine but mainly mannose/glucose", "low",
  "Glucosamine_hydrochloride", "compound", "Q2", int, "K02794", "transport", "mannose PTS EIIAB ManX; carries glucosamine but mainly mannose/glucose", "low",
  "Glucosamine_hydrochloride", "compound", "Q2", int, "K02795", "transport", "mannose PTS EIIC ManY; carries glucosamine but mainly mannose/glucose", "low",
  "Glucosamine_hydrochloride", "compound", "Q2", int, "K02796", "transport", "mannose PTS EIID ManZ; carries glucosamine but mainly mannose/glucose", "low",
  # ---- Assayed enzymes (all labelled extracellular, plan §4) ---------------
  "CB", "enzyme", "Enzyme", ext, "K19668", "assayed activity", "cellulose 1,4-beta-cellobiosidase, EC 3.2.1.91 (also in cellulose set)", "high",
  "BG", "enzyme", "Enzyme", ext, "K01188", "assayed activity", "beta-glucosidase, EC 3.2.1.21", "high",
  "BG", "enzyme", "Enzyme", ext, "K05349", "assayed activity", "beta-glucosidase BglX, EC 3.2.1.21 (periplasmic)", "high",
  "BG", "enzyme", "Enzyme", ext, "K05350", "assayed activity", "beta-glucosidase BglB, EC 3.2.1.21 (mostly cytoplasmic)", "medium",
  "BX", "enzyme", "Enzyme", ext, "K01198", "assayed activity", "xylan 1,4-beta-xylosidase XynB, EC 3.2.1.37", "high",
  "BX", "enzyme", "Enzyme", ext, "K15920", "assayed activity", "xylan 1,4-beta-xylosidase XYL4, EC 3.2.1.37", "high",
  "BX", "enzyme", "Enzyme", ext, "K22268", "assayed activity", "xylan 1,4-beta-xylosidase XylA, EC 3.2.1.37", "high",
  "NAG", "enzyme", "Enzyme", ext, "K01207", "assayed activity", "beta-N-acetylhexosaminidase NagZ, EC 3.2.1.52 (curated list)", "high",
  "NAG", "enzyme", "Enzyme", ext, "K12373", "assayed activity", "hexosaminidase HEXA_B, EC 3.2.1.52 (eukaryotic)", "medium",
  "NAG", "enzyme", "Enzyme", ext, "K14459", "assayed activity", "hexosaminidase HEX, EC 3.2.1.52 (fungal)", "medium",
  "Leu", "enzyme", "Enzyme", ext, "K01255", "assayed activity", "leucyl aminopeptidase PepA, EC 3.4.11.1", "high",
  "Leu", "enzyme", "Enzyme", ext, "K11142", "assayed activity", "cytosol aminopeptidase LAP3, EC 3.4.11.1", "medium",
  "Leu", "enzyme", "Enzyme", ext, "K01256", "assayed activity", "aminopeptidase N PepN; also hydrolyses Leu-AMC (also in Ala set)", "medium",
  "Ala", "enzyme", "Enzyme", ext, "K01256", "assayed activity", "aminopeptidase N PepN (alanyl aminopeptidase), EC 3.4.11.2 (also in Leu set)", "high",
  "Ala", "enzyme", "Enzyme", ext, "K08776", "assayed activity", "puromycin-sensitive aminopeptidase, EC 3.4.11.14", "medium",
  "AP", "enzyme", "Enzyme", ext, "K09474", "assayed activity", "acid phosphatase class A PhoN, EC 3.1.3.2", "high",
  "AP", "enzyme", "Enzyme", ext, "K03788", "assayed activity", "acid phosphatase class B AphA, EC 3.1.3.2", "high",
  "AP", "enzyme", "Enzyme", ext, "K01093", "assayed activity", "4-phytase / acid phosphatase AppA, EC 3.1.3.26 3.1.3.2", "high",
  "AP", "enzyme", "Enzyme", ext, "K14379", "assayed activity", "tartrate-resistant acid phosphatase ACP5, EC 3.1.3.2 (eukaryotic)", "medium",
  "AP", "enzyme", "Enzyme", ext, "K14410", "assayed activity", "lysosomal acid phosphatase ACP2, EC 3.1.3.2 (eukaryotic)", "medium",
  "AP", "enzyme", "Enzyme", ext, "K01077", "assayed activity", "alkaline phosphatase PhoA, EC 3.1.3.1; assay buffer pH 7.91", "low",
  "AP", "enzyme", "Enzyme", ext, "K01113", "assayed activity", "alkaline phosphatase PhoD, EC 3.1.3.1; assay buffer pH 7.91", "low"
)

# Compounds the plan asks for that have no specific gene set
no_set <- tribble(
  ~set, ~set_type, ~group, ~note,
  "Succinic_acid", "compound", "Simple",
  "Transporters only (all shared C4-dicarboxylate carriers, low confidence); the first catabolic step is succinate dehydrogenase (TCA cycle, excluded). Candidate for 'no specific gene set'."
)

# -----------------------------------------------------------------------------
# Attach KEGG names and check every KO is in the expression file
# -----------------------------------------------------------------------------

expr <- read_csv(input_files$brite_d, show_col_types = FALSE)
names(expr)[1] <- "brite_D"
expr <- expr |>
  filter(!is.na(brite_D), brite_D != "") |>
  mutate(KO = str_extract(brite_D, "^K\\d{5}"),
         KEGG_name = str_trim(str_remove(brite_D, "^K\\d{5}\\s+")))
sample_cols <- str_subset(names(expr), "_MT$")

missing_kos <- setdiff(draft$KO, expr$KO)
if (length(missing_kos) > 0) stop("KOs not in expression file: ",
                                  paste(missing_kos, collapse = ", "))

draft_out <- draft |>
  left_join(expr |> select(KO, KEGG_name), by = "KO") |>
  mutate(n_sets = n_distinct(set), .by = KO,
         shared_with = map2_chr(KO, set, \(k, s) {
           paste(setdiff(draft$set[draft$KO == k], s), collapse = ", ")
         })) |>
  select(set, set_type, group, pathway_step, KO, KEGG_name, step, reason,
         confidence, shared_with)
write_csv(draft_out, file.path(dir_results, "KO_gene_sets_DRAFT.csv"))
write_csv(no_set, file.path(dir_results, "KO_gene_sets_DRAFT_notes.csv"))

# -----------------------------------------------------------------------------
# Summed expression of each gene set (units unconfirmed, O2)
# -----------------------------------------------------------------------------

expr_long <- expr |>
  select(KO, all_of(sample_cols)) |>
  pivot_longer(-KO, names_to = "sample", values_to = "value")

sample_totals <- expr_long |>
  summarise(total = sum(value), .by = sample)

ko_summary <- expr_long |>
  filter(KO %in% draft$KO) |>
  summarise(mean_value = mean(value), n_samples_nonzero = sum(value > 0),
            .by = KO)

set_sums <- draft |>
  distinct(set, set_type, group, pathway_step, KO, confidence) |>
  inner_join(expr_long, by = "KO", relationship = "many-to-many") |>
  summarise(value = sum(value), .by = c(set, set_type, group, pathway_step,
                                         sample))

summarise_sets <- function(d) {
  d |>
    left_join(sample_totals, by = "sample") |>
    summarise(mean_sum = mean(value), median_sum = median(value),
              min_sum = min(value), max_sum = max(value),
              n_samples_zero = sum(value == 0),
              mean_share_ppm = mean(value / total) * 1e6,
              .by = c(set, set_type, group, pathway_step))
}

set_expression <- summarise_sets(set_sums) |>
  left_join(draft |>
              summarise(n_KOs = n(), n_high = sum(confidence == "high"),
                        .by = set),
            by = "set")

high_sums <- draft |>
  filter(confidence == "high") |>
  distinct(set, set_type, group, pathway_step, KO) |>
  inner_join(expr_long, by = "KO", relationship = "many-to-many") |>
  summarise(value = sum(value), .by = c(set, set_type, group, pathway_step,
                                         sample)) |>
  summarise_sets() |>
  select(set, mean_sum_high_only = mean_sum, n_samples_zero_high_only = n_samples_zero)

set_expression <- set_expression |>
  left_join(high_sums, by = "set") |>
  arrange(factor(set_type, c("compound", "enzyme")),
          factor(group, c("Polymer", "Simple", "Q2", "Enzyme")), set)
write_csv(set_expression, file.path(dir_results, "KO_gene_sets_DRAFT_expression.csv"))

ko_detail <- draft_out |>
  left_join(ko_summary, by = "KO") |>
  select(set, KO, KEGG_name, confidence, mean_value, n_samples_nonzero)
write_csv(ko_detail, file.path(dir_results, "KO_gene_sets_DRAFT_ko_expression.csv"))

options(width = 200)
cat("Samples:", length(sample_cols), " KOs in file:", n_distinct(expr$KO),
    " KOs in draft:", n_distinct(draft$KO), "\n")
print(set_expression, n = Inf)
print(ko_detail |> filter(n_samples_nonzero == 0) |> select(set, KO, KEGG_name), n = Inf)
