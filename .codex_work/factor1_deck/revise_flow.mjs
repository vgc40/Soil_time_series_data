import fs from "node:fs/promises";
import { FileBlob, PresentationFile } from "@oai/artifact-tool";

const root = "C:/Users/gara009/OneDrive - PNNL/Documents/GitHub/Soil_time_series_data";
const input = `${root}/MEFISTO_Factor1_gene_expression_BRITE_clustering_results.pptx`;
const output = `${root}/MEFISTO_Factor1_gene_expression_BRITE_clustering_results_revised.pptx`;
const work = `${root}/.codex_work/factor1_deck/revised-render`;

const presentation = await PresentationFile.importPptx(await FileBlob.load(input));

// Remove the dense GAM mosaic and the heatmap comparison.
presentation.slides.items[7].delete();
presentation.slides.items[3].delete();

function shape(slideNumber, name) {
  const slide = presentation.slides.items[slideNumber - 1];
  const result = slide.shapes.items.find((item) => item.name === name);
  if (!result) throw new Error(`Missing ${name} on slide ${slideNumber}`);
  return result;
}

// Reframe the elbow slide as the bridge between filtering and trend clusters.
shape(4, "Rectangle 5").text = "FILTERING AND CLUSTER SELECTION";
shape(4, "Rectangle 1").text = "MicroResp-linked Factor 1 functions resolve into three temporal modes";
shape(4, "Rectangle 10").text = "The MicroResp map matched 31 Factor 1 pathways and 1,113 Factor 1 KOs. After model matching, 31 pathways and 985 KOs were clustered; both elbow diagnostics selected k = 3.";
shape(4, "Rectangle 11").text = "Target-set FDR significance: 20 BRITE C pathways and 230 BRITE D KOs.";

// Make the two central trend slides explicit and parallel.
shape(5, "Rectangle 5").text = "BRITE C CLUSTER TRENDS";
shape(5, "Rectangle 1").text = "BRITE C pathways form three distinct temporal expression trajectories";
shape(5, "Rectangle 14").text = "Storyboard connection";
shape(5, "Rectangle 15").text = "Cluster 3 rises from winter through May and most closely resembles the storyboard's winter/early-spring Factor 1 signal. Clusters 1 and 2 capture early-high and rebound behavior.";

shape(6, "Rectangle 5").text = "BRITE D CLUSTER TRENDS";
shape(6, "Rectangle 1").text = "BRITE D KOs independently recover the same three timing modes";
shape(6, "Rectangle 14").text = "Storyboard connection";
shape(6, "Rectangle 15").text = "Cluster 3 contains 263 KOs that rise from winter through May. KO-level resolution also reveals a large early-high group and a strong May-rebound group.";

// Renumber visible slide markers after removal.
for (let slideNumber = 2; slideNumber <= presentation.slides.items.length; slideNumber += 1) {
  shape(slideNumber, "Rectangle 3").text = String(slideNumber);
}

await fs.mkdir(work, { recursive: true });
for (let i = 0; i < presentation.slides.items.length; i += 1) {
  const slide = presentation.slides.items[i];
  const png = await presentation.export({ slide, format: "png", scale: 1 });
  await fs.writeFile(`${work}/slide-${i + 1}.png`, new Uint8Array(await png.arrayBuffer()));
  const layout = await slide.export({ format: "layout" });
  await fs.writeFile(`${work}/slide-${i + 1}.layout.json`, await layout.text());
}
const montage = await presentation.export({ format: "webp", montage: true, scale: 1 });
await fs.writeFile(`${work}/montage.webp`, new Uint8Array(await montage.arrayBuffer()));

const pptx = await PresentationFile.exportPptx(presentation);
await pptx.save(output);
console.log(output);
