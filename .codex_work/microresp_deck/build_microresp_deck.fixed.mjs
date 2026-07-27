import fs from "node:fs/promises";
import { FileBlob, PresentationFile } from "@oai/artifact-tool";

const root = "C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/Soil_time_series_data";
const work = `${root}/.codex_work/microresp_deck`;
const starter = `${work}/template-starter.pptx`;
const output = `${root}/MicroResp_BRITE_C_D_targeted_temporal_results.pptx`;

const presentation = await PresentationFile.importPptx(await FileBlob.load(starter));
const inspectText = await fs.readFile(`${work}/template-starter.pptx.inspect.ndjson`, "utf8");
const originalRecords = new Map(
  inspectText.split(/\\r?\\n/).filter(Boolean).map((line) => JSON.parse(line))
    .filter((record) => record.id).map((record) => [record.id, record]),
);
function resolveOriginal(id) {
  const record = originalRecords.get(id);
  if (!record) throw new Error(`Missing starter inspection record: ${id}`);
  const slide = presentation.slides.items[record.slide - 1];
  if (record.kind === "notes") return slide.speakerNotes;
  if (record.kind === "textbox" || record.kind === "shape") {
    const shape = slide.shapes.items.find((item) => item.name === record.name);
    if (!shape) throw new Error(`Missing shape ${record.name} on slide ${record.slide}`);
    return shape;
  }
  if (record.kind === "image") {
    const image = slide.images.items.find((item) => item.name === record.name);
    if (!image) throw new Error(`Missing image ${record.name} on slide ${record.slide}`);
    return image;
  }
  throw new Error(`Unsupported record kind for ${id}: ${record.kind}`);
}

function setText(id, value) {
  resolveOriginal(id).text = value;
}

async function replaceImage(id, path, alt) {
  const image = resolveOriginal(id);
  const bytes = await fs.readFile(path);
  const oldFrame = image.frame;
  const oldCrop = image.crop;
  const oldFit = image.fit;
  const oldGeometry = image.geometry;
  const oldRotation = image.rotation;
  image.replace({
    blob: bytes.buffer.slice(bytes.byteOffset, bytes.byteOffset + bytes.byteLength),
    contentType: "image/png",
    alt,
    ...(oldFit ? { fit: oldFit } : {}),
  });
  image.frame = oldFrame;
  image.crop = oldCrop;
  image.geometry = oldGeometry;
  image.rotation = oldRotation;
}

function setNotes(id, lines) {
  resolveOriginal(id).setText(`[Sources]\n${lines.join("\n")}`);
}

// Slide 1
setText("sh/kjad0fet", "SOILS SFA | PROSSER TIME SERIES");
setText("sh/wfqtwvex", "MicroResp-linked expression resolves three seasonal programs");
setText("sh/yh8by5wn", "Targeted BRITE C pathways and BRITE D KOs: GAM review, temporal clustering, and storyboard synthesis");
setText("sh/jihc7ad8", "100% irrigation metatranscriptomes | Seven dates | July 2026");
setNotes("nt/tbtajt", [
  "Storyboard: SFA3-1-1_VGC_Storyboard_v2.docx",
  "Analysis: Scripts/12_MicroResp_targeted_BRITE_C_D_GAM_and_clustering.R",
]);

// Slide 2
setText("sh/hgzm9crq", "MICRORESP-LINKED EXPRESSION");
setText("sh/hsrq14vy", "Filtering narrows the transcriptome to functions tied to assayed substrates");
setText("sh/72p8vudw", "2");
setText("sh/m1gr2pwb", "The explicit MicroResp-to-KEGG map was joined to previously fitted random-plot GAMs. All mapped, converged trajectories were clustered after within-feature z-scoring; targeted-set FDR identifies the strongest time-varying subset.");
setText("sh/zy58rux0", "31");
setText("sh/yxw7ypwf", "BRITE C pathways mapped and clustered");
setText("sh/a58n6x03", "985");
setText("sh/bq14f218", "BRITE D KOs clustered from 1,312 requested");
setText("sh/218n2d0r", "k = 3");
setText("sh/n2h4bi1c", "automatic elbow at both annotation levels");
setText("sh/8zu5gri1", "Significant within the targeted set: 20 BRITE C pathways and 230 BRITE D KOs. Distance: 1 - Pearson correlation; linkage: Ward.D2.");
setText("sh/9036pwzm", "Source: Results/MicroResp_targeted_BRITE_C_D/*_run_summary.csv");
setNotes("nt/cvatf6", [
  "Results/MicroResp_targeted_BRITE_C_D/MicroResp_BRITE_C_D_run_summary.csv",
  "Data/MicroResp_all_BRITE_C_D_complete.csv",
  "Scripts/12_MicroResp_targeted_BRITE_C_D_GAM_and_clustering.R",
]);

// Slide 3
setText("sh/sn6tgrut", "MICRORESP-LINKED EXPRESSION");
setText("sh/grqtkbup", "Both functional resolutions support three temporal expression clusters");
setText("sh/6lobehcn", "3");
setText("sh/jixs3mdw", "BRITE C");
setText("sh/ripczip4", "BRITE D");
setText("sh/6hgvqd8j", "Within-cluster dissimilarity drops sharply through k = 3 at both levels. Higher k values continue to subdivide trajectories but add progressively smaller gains.");
setText("sh/tk7u187a", "Source: MicroResp-targeted BRITE C/D elbow diagnostics; automatic selection, no manual override.");
await replaceImage("im/bidgb2xo", `${root}/Figures/MicroResp_targeted_BRITE_C_D/brite_C/brite_C_elbow_diagnostic.png`, "BRITE C elbow diagnostic");
await replaceImage("im/0n6xo7e5", `${root}/Figures/MicroResp_targeted_BRITE_C_D/brite_D/brite_D_elbow_diagnostic.png`, "BRITE D elbow diagnostic");
setNotes("nt/c72spz", [
  "Results/MicroResp_targeted_BRITE_C_D/brite_C_MicroResp_elbow_diagnostics.csv",
  "Results/MicroResp_targeted_BRITE_C_D/brite_D_MicroResp_elbow_diagnostics.csv",
]);

// Slide 4
setText("sh/q5knmdcb", "MICRORESP-LINKED BRITE C");
setText("sh/210nidcf", "Twenty of 31 mapped pathways show targeted-set temporal significance");
setText("sh/c725onuh", "4");
setText("sh/sby5cjud", "Cluster 1 | 12 pathways: biosynthesis, glycan use, and aromatic degradation");
setText("sh/dovaxkf6", "Cluster 2 | 13 pathways: TCA, N cycling, amino-acid and fatty-acid metabolism");
setText("sh/rmd8vax0", "Cluster 3 | 6 pathways: glycolysis, pyruvate, lipid, butanoate functions");
setText("sh/0r69szeh", "GAM evidence");
setText("sh/zqxszuxw", "The mapped pathway set contains broad, nonlinear time trends rather than a single synchronized response. Strong signals include arginine biosynthesis, N metabolism, amino-sugar metabolism, and aromatic-compound degradation.");
setText("sh/ep4rqpwr", "Source: brite_C_MicroResp_GAM_statistics.csv and targeted GAM panels");
await replaceImage("im/e1cf61sb", `${root}/Figures/MicroResp_targeted_BRITE_C_D/brite_C/brite_C_targeted_GAMs_page_01.png`, "Example GAM panels for MicroResp-linked BRITE C pathways");
setNotes("nt/9ws2ke", [
  "Results/MicroResp_targeted_BRITE_C_D/brite_C_MicroResp_GAM_statistics.csv",
  "Figures/MicroResp_targeted_BRITE_C_D/brite_C/brite_C_targeted_GAMs_page_01.png",
]);

// Slide 5
setText("sh/yhsrqtwr", "MICRORESP-LINKED BRITE C");
setText("sh/m5cbudwn", "Pathway clusters separate broad winter activity from shorter seasonal modes");
setText("sh/kja9s3ex", "5");
setText("sh/5gnaxsvm", "Cluster 1 | 12 pathways");
setText("sh/1gv6xsvm", "Pentose phosphate, glycan degradation, amino-acid biosynthesis, benzoate/aminobenzoate degradation, and folate-linked C1 metabolism.");
setText("sh/etknmdcv", "Cluster 2 | 13 pathways");
setText("sh/dsbmtsva", "TCA and glyoxylate cycles, N metabolism, amino-sugar metabolism, fatty-acid degradation, and several amino-acid degradation routes.");
setText("sh/qpknidcz", "Cluster 3 | 6 pathways");
setText("sh/xof65cve", "Glycolysis, pyruvate and butanoate metabolism, glycerolipid and glycerophospholipid metabolism, plus naphthalene degradation.");
setText("sh/wn65w7u9", "The heatmap includes every mapped, converged BRITE C pathway. Rows are z-scored, so color represents trajectory shape rather than absolute expression.");
setText("sh/wrepsv25", "Source: brite_C_MicroResp_cluster_assignments.csv and brite_C_MicroResp_heatmap.png");
await replaceImage("im/ypgrip0j", `${root}/Figures/MicroResp_targeted_BRITE_C_D/brite_C/brite_C_MicroResp_heatmap.png`, "Heatmap of all mapped BRITE C pathway trajectories");
setNotes("nt/opm0y2", [
  "Results/MicroResp_targeted_BRITE_C_D/brite_C_MicroResp_cluster_assignments.csv",
  "Figures/MicroResp_targeted_BRITE_C_D/brite_C/brite_C_MicroResp_heatmap.png",
]);

// Slide 6
setText("sh/lcn6toni", "MICRORESP-LINKED BRITE D");
setText("sh/tgn6xone", "KO-level GAMs reveal strong but heterogeneous seasonal expression");
setText("sh/ja5ore5c", "6");
setText("sh/36po3ilw", "Cluster 1 | 477 KOs (48%): includes urease, polyamine, TCA, and aromatic-use genes");
setText("sh/0z2h4vq9", "Cluster 2 | 245 KOs (25%): includes nitrate reduction and storage-carbon functions");
setText("sh/edkz2l8j", "Cluster 3 | 263 KOs (27%): includes sugar, amino-acid, and N-oxide functions");
setText("sh/dcbitg7y", "Resolution matters");
setText("sh/q90zyl87", "Of 1,312 mapped KOs requested, 985 converged and were clustered; 230 were temporally significant within the targeted set. KO diversity shows that pathway-level averages can conceal opposing constituent trajectories.");
setText("sh/ra90rqps", "Source: brite_D_MicroResp_GAM_statistics.csv and targeted GAM panels");
await replaceImage("im/z6hwfi10", `${root}/Figures/MicroResp_targeted_BRITE_C_D/brite_D/brite_D_targeted_GAMs_page_01.png`, "Example GAM panels for MicroResp-linked BRITE D KOs");
setNotes("nt/ac20vx", [
  "Results/MicroResp_targeted_BRITE_C_D/brite_D_MicroResp_GAM_statistics.csv",
  "Figures/MicroResp_targeted_BRITE_C_D/brite_D/brite_D_targeted_GAMs_page_01.png",
]);

// Slide 7
setText("sh/md8j2dgr", "MICRORESP-LINKED BRITE D");
setText("sh/y98jydgf", "The strongest KO trajectories connect C use with nitrogen and stress metabolism");
setText("sh/8fa14nyh", "7");
setText("sh/pkf21ozy", "Cluster 1 | 477 KOs");
setText("sh/p4vut0ji", "Representative high-variance features include ureA/ureB, ornithine decarboxylase, aconitate hydratase, PRPP synthetase, and aromatic-degradation enzymes.");
setText("sh/214vylkr", "Cluster 2 | 245 KOs");
setText("sh/x8fex0je", "Representative features include narG/narH, nirK, glgC, glycerol-3-phosphate dehydrogenase, and amino-acid biosynthesis functions.");
setText("sh/a5ov25k3", "Cluster 3 | 263 KOs");
setText("sh/lwze1kjq", "Representative features include alcohol dehydrogenase, arabinose isomerase, chitobiose phosphorylase, asparagine synthetase, and norB.");
setText("sh/kv6d8f25", "For readability, the BRITE D heatmap shows the top 15 features per cluster ranked by time-only variance explained; all 985 converged KOs remain in the clustering.");
setText("sh/kfex4jah", "Source: brite_D_MicroResp_heatmap_feature_selection.csv");
await replaceImage("im/a18rml8r", `${root}/Figures/MicroResp_targeted_BRITE_C_D/brite_D/brite_D_MicroResp_heatmap.png`, "Selected BRITE D KO heatmap, top 15 per cluster");
setNotes("nt/owshb2", [
  "Results/MicroResp_targeted_BRITE_C_D/brite_D_MicroResp_heatmap_feature_selection.csv",
  "Results/MicroResp_targeted_BRITE_C_D/brite_D_MicroResp_cluster_assignments.csv",
]);

// Slide 8
setText("sh/f69oj2t8", "STORYBOARD SYNTHESIS");
setText("sh/9wfe1cri", "Targeted expression partly supports the winter phenotype, not a one-to-one match");
setText("sh/byxw3m9o", "8");
setText("sh/wz6dcrq9", "Storyboard claim");
setText("sh/y18feh8z", "Targeted BRITE C");
setText("sh/hsr2pof2", "Targeted BRITE D");
setText("sh/5wr2t8vy", "Broad winter/early-spring phenotype");
setText("sh/ju9kryds", "Three modes; N, TCA, amino-sugar, and fatty-acid pathways represented");
setText("sh/tkbmx8va", "Three modes; nitrate, urease, C-use, amino-acid, and stress KOs represented");
setText("sh/edw36ton", "Sharp Sep-Oct fall pulse");
setText("sh/o3ylcn6p", "No single dominant pathway cluster uniquely reproduces the pulse");
setText("sh/2hg3adoz", "No single dominant KO cluster uniquely reproduces the pulse");
setText("sh/v218nqhs", "Cross-scale conclusion");
setText("sh/h4jqp0zy", "Partial coordination across expressed pathways, consistent with the storyboard");
setText("sh/j618raho", "Heterogeneity at KO resolution reinforces assay-specific variation");
setText("sh/zah8f6hk", "Interpretation: seasonal substrate-use capacity and expression overlap, but transient phenotypes may require denser sampling or remain assay-specific.");
setNotes("nt/f9oh6p", [
  "SFA3-1-1_VGC_Storyboard_v2.docx",
  "Results/MicroResp_targeted_BRITE_C_D/brite_C_MicroResp_cluster_assignments.csv",
  "Results/MicroResp_targeted_BRITE_C_D/brite_D_MicroResp_cluster_assignments.csv",
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
