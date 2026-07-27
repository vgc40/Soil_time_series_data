import fs from "node:fs/promises";
import { FileBlob, PresentationFile } from "@oai/artifact-tool";

const root = "C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/Soil_time_series_data";
const work = `${root}/.codex_work/factor1_deck`;
const starter = `${work}/template-starter.pptx`;
const output = `${root}/MEFISTO_Factor1_gene_expression_BRITE_clustering_results.pptx`;

const presentation = await PresentationFile.importPptx(await FileBlob.load(starter));
const sourceInspect = await fs.readFile(`${work}/template-inspect/template-inspect.ndjson`, "utf8");
const sourceRecordList = sourceInspect.split(/\r?\n/).filter(Boolean)
  .map((line) => JSON.parse(line)).filter((record) => record.id);
const sourceRecords = new Map(sourceRecordList.map((record) => [record.id, record]));
const currentSnapshot = await presentation.inspect({
  kind: "slide,textbox,shape,image,notes",
  maxChars: 200000,
});
const currentRecords = currentSnapshot.ndjson.split(/\r?\n/).filter(Boolean)
  .map((line) => JSON.parse(line));

function getElement(outputSlide, sourceId) {
  const sourceRecord = sourceRecords.get(sourceId);
  if (!sourceRecord) throw new Error(`Missing source record: ${sourceId}`);
  const targetKind = sourceRecord.kind === "textbox" ? "textbox" : sourceRecord.kind;
  let match = currentRecords.find((record) => {
    if (record.slide !== outputSlide || record.kind !== targetKind) return false;
    if (sourceRecord.name && record.name === sourceRecord.name) return true;
    if (!sourceRecord.bbox || !record.bbox) return false;
    return sourceRecord.bbox.every((value, index) => Math.abs(value - record.bbox[index]) < 1);
  });
  if (!match) {
    const sourcePeers = sourceRecordList.filter((record) =>
      record.slide === sourceRecord.slide && record.kind === sourceRecord.kind);
    const ordinal = sourcePeers.findIndex((record) => record.id === sourceId);
    const currentPeers = currentRecords.filter((record) =>
      record.slide === outputSlide && record.kind === targetKind);
    match = currentPeers[ordinal];
  }
  if (!match) throw new Error(`No current match for ${sourceId} on slide ${outputSlide}`);
  return presentation.resolve(match.id);
}

function setText(slide, id, value) {
  getElement(slide, id).text = value;
}

async function replaceImage(slide, id, path, alt) {
  const image = getElement(slide, id);
  const bytes = await fs.readFile(path);
  const frame = image.frame;
  const geometry = image.geometry;
  const rotation = image.rotation;
  image.replace({
    blob: bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength),
    contentType: "image/png",
    alt,
    fit: "contain",
  });
  image.frame = frame;
  image.geometry = geometry;
  image.rotation = rotation;
}

function setNotes(slideNumber, lines) {
  presentation.slides.items[slideNumber - 1].speakerNotes.setText(
    `[Sources]\n${lines.join("\n")}`,
  );
}

// Slide 1
setText(1, "sh/ts7md4r2", "SOILS SFA | PROSSER TIME SERIES");
setText(1, "sh/547294r6", "MEFISTO Factor 1 genes contain three seasonal expression programs");
setText(1, "sh/7qp4be9c", "Factor 1 gene selection, BRITE aggregation, temporal GAMs, MicroResp overlap, and trajectory clustering");
setText(1, "sh/65g3298r", "100% irrigation metatranscriptomes | Seven dates | July 2026");
setNotes(1, [
  "Scripts/13_Factor1_expression_MicroResp_BRITE_C_D.R",
  "SFA3-1-1_VGC_Storyboard_v2.docx",
]);

// Slide 2
setText(2, "sh/t8byxkn2", "FACTOR 1 FOLLOW-ON ANALYSIS");
setText(2, "sh/9072xkry", "The workflow asks which Factor 1 transcripts align with MicroResp-linked functions");
setText(2, "sh/b29kza94", "2");
setText(2, "sh/a10jqpsj", "MEFISTO Factor 1 features were filtered from the original expression matrix, annotated to KEGG, aggregated to BRITE C pathways and BRITE D KOs, modeled through time, intersected with the MicroResp map, and clustered by fitted trajectory shape.");
setText(2, "sh/r65knqtk", "770k");
setText(2, "sh/q5wjelsz", "unique Factor 1 transcript features");
setText(2, "sh/mtwrmxg7", "516");
setText(2, "sh/nu58f2hs", "unique BRITE C pathways represented");
setText(2, "sh/ahwrqhgj", "6,789");
setText(2, "sh/bip8jmho", "unique BRITE D KOs represented");
setText(2, "sh/cnupgny5", "GAM: Expression ~ s(time, k = 6) + s(Plot, random effect). Clustering: row z-score, 1 - Pearson correlation, Ward.D2.");
setText(2, "sh/do3q9szq", "Source: factor1_gene_summary.csv and run_manifest.csv");
setNotes(2, [
  "Results/Factor1_BRITE_MicroResp/factor1_gene_summary.csv",
  "Results/Factor1_BRITE_MicroResp/run_manifest.csv",
  "Scripts/13_Factor1_expression_MicroResp_BRITE_C_D.R",
]);

// Slide 3
setText(3, "sh/obq90bml", "FACTOR 1 FUNCTIONAL BREADTH");
setText(3, "sh/cza94vmx", "Factor 1 genes span broad pathway- and KO-level functional space");
setText(3, "sh/298ryl4v", "3");
setText(3, "sh/n6ls3alk", "BRITE C");
setText(3, "sh/n2l4fq98", "BRITE D");
setText(3, "sh/m1c3mlsn", "The most gene-rich annotations include central metabolism, biosynthesis, transport, energy conservation, and information-processing functions.");
setText(3, "sh/943mhgre", "Bars count Factor 1 genes per annotation; they do not show loading sign or temporal significance.");
await replaceImage(3, "im/3q5onqdc", `${root}/Figures/Factor1_BRITE_MicroResp/factor1_brite_C_gene_counts.png`, "Top BRITE C pathways by number of Factor 1 genes");
await replaceImage(3, "im/gne5cvel", `${root}/Figures/Factor1_BRITE_MicroResp/factor1_brite_D_gene_counts.png`, "Top BRITE D KOs by number of Factor 1 genes");
setNotes(3, [
  "Results/Factor1_BRITE_MicroResp/factor1_brite_C_summary.csv",
  "Results/Factor1_BRITE_MicroResp/factor1_brite_D_summary.csv",
  "Figures/Factor1_BRITE_MicroResp/factor1_brite_C_gene_counts.png",
  "Figures/Factor1_BRITE_MicroResp/factor1_brite_D_gene_counts.png",
]);

// Slide 4
setText(4, "sh/obq90bml", "FACTOR 1 TEMPORAL GAMS");
setText(4, "sh/cza94vmx", "Temporal structure is widespread, but not universal, across Factor 1 functions");
setText(4, "sh/298ryl4v", "4");
setText(4, "sh/n6ls3alk", "BRITE C");
setText(4, "sh/n2l4fq98", "BRITE D");
setText(4, "sh/m1c3mlsn", "Of 348 pathway models, 337 converged and 145 were time-significant. Of 6,338 KO models, 5,784 converged and 925 were time-significant after FDR correction.");
setText(4, "sh/943mhgre", "Panels show the FDR-significant subset; fitted curves exclude the plot random effect.");
await replaceImage(4, "im/3q5onqdc", `${root}/Figures/Factor1_BRITE_MicroResp/temporal_gam/brite_C_significant_temporal_GAMs.png`, "Significant Factor 1 BRITE C temporal GAMs");
await replaceImage(4, "im/gne5cvel", `${root}/Figures/Factor1_BRITE_MicroResp/temporal_gam/brite_D_significant_temporal_GAMs.png`, "Significant Factor 1 BRITE D temporal GAMs");
setNotes(4, [
  "Results/Factor1_BRITE_MicroResp/temporal_gam/Factor1_GAM_run_summary.csv",
  "Figures/Factor1_BRITE_MicroResp/temporal_gam/brite_C_significant_temporal_GAMs.png",
  "Figures/Factor1_BRITE_MicroResp/temporal_gam/brite_D_significant_temporal_GAMs.png",
]);

// Slide 5
setText(5, "sh/obq90bml", "MICRORESP OVERLAP AND CLUSTERING");
setText(5, "sh/cza94vmx", "The MicroResp-linked Factor 1 subset supports three expression clusters");
setText(5, "sh/298ryl4v", "5");
setText(5, "sh/n6ls3alk", "BRITE C");
setText(5, "sh/n2l4fq98", "BRITE D");
setText(5, "sh/m1c3mlsn", "MicroResp mapping matched 31/516 Factor 1 pathways and 1,113/6,789 Factor 1 KOs. After model matching, 31 pathways and 985 KOs were clustered; both elbows selected k = 3.");
setText(5, "sh/943mhgre", "Target-set FDR significance: 20 BRITE C pathways and 230 BRITE D KOs.");
await replaceImage(5, "im/3q5onqdc", `${root}/Figures/Factor1_BRITE_MicroResp/clustering/brite_C/brite_C_Factor1_elbow_diagnostic.png`, "Factor 1 MicroResp-linked BRITE C elbow diagnostic");
await replaceImage(5, "im/gne5cvel", `${root}/Figures/Factor1_BRITE_MicroResp/clustering/brite_D/brite_D_Factor1_elbow_diagnostic.png`, "Factor 1 MicroResp-linked BRITE D elbow diagnostic");
setNotes(5, [
  "Results/Factor1_BRITE_MicroResp/factor1_microresp_overlap_summary.csv",
  "Results/Factor1_BRITE_MicroResp/clustering/Factor1_MicroResp_clustering_run_summary.csv",
  "Results/Factor1_BRITE_MicroResp/clustering/brite_C_Factor1_MicroResp_elbow_diagnostics.csv",
  "Results/Factor1_BRITE_MicroResp/clustering/brite_D_Factor1_MicroResp_elbow_diagnostics.csv",
]);

// Slide 6
setText(6, "sh/doj29oba", "FACTOR 1 | BRITE C");
setText(6, "sh/1cj2d8b6", "Pathway clusters separate early-high, rebound, and late-activation modes");
setText(6, "sh/rm1k7yt4", "6");
setText(6, "sh/3ihk3et8", "Cluster 1 | 12 pathways: highest in Jun-Jul, minimum in Jan");
setText(6, "sh/kbm987y5", "Cluster 2 | 13 pathways: early-high, Jan minimum, strong May rebound");
setText(6, "sh/i94r6xgz", "Cluster 3 | 6 pathways: low through Oct, then rises to a May maximum");
setText(6, "sh/x8vaxsfe", "Connection to Factor 1");
setText(6, "sh/a5kr2xg3", "Cluster 3 most closely resembles the winter/early-spring activation described in the storyboard. Clusters 1 and 2 show that Factor 1 membership also contains early-season and rebound programs.");
setText(6, "sh/v6tsv2xo", "Source: brite_C Factor 1 cluster assignments and average trajectories");
await replaceImage(6, "im/b6xofahg", `${root}/Figures/Factor1_BRITE_MicroResp/clustering/brite_C/brite_C_Factor1_cluster_average_trajectories.png`, "Factor 1 MicroResp-linked BRITE C cluster-average trajectories");
setNotes(6, [
  "Results/Factor1_BRITE_MicroResp/clustering/brite_C_Factor1_MicroResp_cluster_assignments.csv",
  "Results/Factor1_BRITE_MicroResp/clustering/brite_C_Factor1_MicroResp_cluster_average_trajectories.csv",
]);

// Slide 7
setText(7, "sh/kzmdova1", "FACTOR 1 | BRITE D");
setText(7, "sh/cb2tkvap", "KO-level clusters resolve a large sustained late-activation program");
setText(7, "sh/ydkbm5sv", "7");
setText(7, "sh/eh0ba1sr", "Cluster 1 | 477 KOs: highest in Jun-Jul, sustained decline afterward");
setText(7, "sh/baxkjqlk", "Cluster 2 | 245 KOs: early-high, Jan minimum, strong May rebound");
setText(7, "sh/xcz2l03q", "Cluster 3 | 263 KOs: low in summer, increases from Jan through May");
setText(7, "sh/25ojml43", "Functional examples");
setText(7, "sh/p8fih03e", "Late/rebound groups include nitrate and nitrite reduction, storage-carbon, amino-acid, sugar-use, and N-oxide functions. Early-high genes include urease, polyamine, TCA, and aromatic-use functions.");
setText(7, "sh/o761ov2t", "Source: brite_D Factor 1 cluster assignments and average trajectories");
await replaceImage(7, "im/8rqxsni1", `${root}/Figures/Factor1_BRITE_MicroResp/clustering/brite_D/brite_D_Factor1_cluster_average_trajectories.png`, "Factor 1 MicroResp-linked BRITE D cluster-average trajectories");
setNotes(7, [
  "Results/Factor1_BRITE_MicroResp/clustering/brite_D_Factor1_MicroResp_cluster_assignments.csv",
  "Results/Factor1_BRITE_MicroResp/clustering/brite_D_Factor1_MicroResp_cluster_average_trajectories.csv",
]);

// Slide 8
setText(8, "sh/obq90bml", "FACTOR 1 TRAJECTORY HEATMAPS");
setText(8, "sh/cza94vmx", "The three timing modes are coherent across pathway and KO resolutions");
setText(8, "sh/298ryl4v", "8");
setText(8, "sh/n6ls3alk", "BRITE C");
setText(8, "sh/n2l4fq98", "BRITE D");
setText(8, "sh/m1c3mlsn", "BRITE C displays all 31 mapped pathways. BRITE D displays the top 15 KOs per cluster by time-only variance explained, while all 985 KOs remain in the clustering.");
setText(8, "sh/943mhgre", "Rows are z-scored; colors compare temporal shape, not expression magnitude.");
await replaceImage(8, "im/3q5onqdc", `${root}/Figures/Factor1_BRITE_MicroResp/clustering/brite_C/brite_C_Factor1_MicroResp_heatmap.png`, "Factor 1 MicroResp-linked BRITE C heatmap");
await replaceImage(8, "im/gne5cvel", `${root}/Figures/Factor1_BRITE_MicroResp/clustering/brite_D/brite_D_Factor1_MicroResp_heatmap.png`, "Factor 1 MicroResp-linked BRITE D heatmap");
setNotes(8, [
  "Results/Factor1_BRITE_MicroResp/clustering/brite_C_Factor1_MicroResp_cluster_assignments.csv",
  "Results/Factor1_BRITE_MicroResp/clustering/brite_D_Factor1_MicroResp_heatmap_feature_selection.csv",
]);

// Slide 9
setText(9, "sh/t8byxkn2", "STORYBOARD SYNTHESIS");
setText(9, "sh/9072xkry", "Factor 1 supports partial coordination, with an important sign-specific caveat");
setText(9, "sh/b29kza94", "9");
setText(9, "sh/a10jqpsj", "The late-activation clusters support the storyboard's interpretation that part of the Factor 1 transcript signal tracks the winter/early-spring MicroResp phenotype. Other clusters reveal early-high and rebound programs within the same Factor 1 feature set.");
setText(9, "sh/r65knqtk", "Support");
setText(9, "sh/q5wjelsz", "263 KOs and 6 pathways rise from winter toward May");
setText(9, "sh/mtwrmxg7", "Caveat");
setText(9, "sh/nu58f2hs", "This analysis combines positive and negative Factor 1 loadings");
setText(9, "sh/ahwrqhgj", "Next");
setText(9, "sh/bip8jmho", "Repeat clustering separately by loading sign");
setText(9, "sh/cnupgny5", "Factor 1 membership is broader than the sign-specific reconstructed trajectory. Separating loading direction would provide the cleanest direct test of the storyboard's negative-weight winter signal.");
setText(9, "sh/do3q9szq", "Source: Factor 1 cluster outputs and manuscript storyboard");
setNotes(9, [
  "SFA3-1-1_VGC_Storyboard_v2.docx",
  "Data/factor1_all_omic_genes.csv",
  "Results/Factor1_BRITE_MicroResp/clustering/Factor1_MicroResp_clustering_run_summary.csv",
  "Scripts/13_Factor1_expression_MicroResp_BRITE_C_D.R",
]);

await fs.mkdir(`${work}/final-render`, { recursive: true });
for (let i = 0; i < presentation.slides.items.length; i += 1) {
  const slide = presentation.slides.items[i];
  const png = await presentation.export({ slide, format: "png", scale: 1 });
  await fs.writeFile(`${work}/final-render/slide-${i + 1}.png`, new Uint8Array(await png.arrayBuffer()));
  const layout = await slide.export({ format: "layout" });
  await fs.writeFile(`${work}/final-render/slide-${i + 1}.layout.json`, await layout.text());
}
const montage = await presentation.export({ format: "webp", montage: true, scale: 1 });
await fs.writeFile(`${work}/final-render/montage.webp`, new Uint8Array(await montage.arrayBuffer()));
const pptx = await PresentationFile.exportPptx(presentation);
await pptx.save(output);
console.log(output);
