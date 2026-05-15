# LongViewSC

**Interactive visualisation of gene and isoform expression in long-read single-cell data**

Developed by the [Clark Laboratory](https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/clark-lab), University of Melbourne

🌐 [Open App Online](https://longviewsc.researchsoftware.unimelb.edu.au/LongViewSC/) &nbsp;|&nbsp; 📖 [FLAMESv2 Data Prep Tutorial](https://sefi196.github.io/FLAMESv2_LR_sc_tutorial/)

---

## Workflow

| Step | Action |
|------|--------|
| **1. Load your data** | Upload a Seurat object (`.rds`, `.qs`, or `.qs2`) with gene and isoform expression data, or try the built-in demo dataset. |
| **2. Configure** | Type a gene name in the sidebar. Select your isoform assay, dimensional reduction, and metadata grouping column. Optionally load a GTF for transcript structure. |
| **3. Run Analysis** | Press **Run Analysis** to generate all plots. The top 4 isoforms are pre-selected — use the coloured chips in the sidebar to toggle isoforms across every tab. |
| **4. Explore & export** | Navigate the tabs to explore cell-type maps, isoform statistics, expression plots, transcript structures, heatmaps, and proportions. Download individual plots or generate a full HTML report. |

---

## Visualisation Tabs

### 1 — Overview
A high-level view of your dataset: a UMAP coloured by your chosen grouping variable, a gene-level feature plot, and a violin plot comparing expression across groups.

- Select a **dimensional reduction** (e.g. UMAP, PCA) in the sidebar.
- Choose a **metadata column** to colour cells by (cluster, cell type, condition).
- Hit **Run Analysis** to generate all plots — the same button updates every tab.

> **Quick start:** Click **Demo Data** on the landing page to load the built-in dataset, then press Run Analysis.

---

### 2 — Isoform Statistics
A searchable, sortable table listing every isoform detected for your selected gene, with expression and detection statistics. Use this to identify the most relevant isoforms before diving into other plots.

- Isoforms are named `ENST00000XXXXX.X-GENE` (e.g. `ENST00000224237.9-VIM`).
- A tick (✓) marks isoforms currently selected in the sidebar chip panel.
- Clicking a row also selects/deselects that isoform.

> **Chip selector:** Coloured chips in the sidebar toggle which isoforms appear across the Expression, Heatmap, Proportions and Trajectory tabs.

---

### 3 — Isoform Expression
Single-cell resolution plots for each selected isoform: a feature plot projected onto the isoform-level UMAP embedding and a violin plot across groups.

- **Feature plots** use the isoform-level UMAP (`umap_iso`) if present, otherwise the gene-level UMAP.
- **Violin plots** show expression per isoform, grouped by your chosen metadata column.
- Only isoforms selected via the chip panel are shown.

> **Important:** Set the **Isoform Assay** field in the sidebar to the correct assay name (e.g. `iso`) before clicking Run Analysis.

---

### 4 — Transcript Structure
Visualises the exon-intron architecture of each selected isoform from your uploaded GTF file, enabling direct structural comparison.

- Exons are drawn as filled blocks; introns as directional arrows showing strand.
- Upload a **.gtf** file in the sidebar — if omitted this tab is empty.
- GTF transcript IDs must match the ENST IDs in your isoform assay.

> **Naming:** Isoforms must follow the format `ENST00000XXXXX.X-GENENAME` for automatic GTF matching.

---

### 5 — Pseudobulk Heatmap
Aggregates isoform counts per group, normalises to log(CPM+1), and renders an interactive clustered heatmap — ideal for comparing relative isoform usage across conditions or cell types at a glance.

- Rows are isoforms; columns are groups; colour intensity = log(CPM+1).
- Rows are clustered by expression similarity.
- The heatmap is fully interactive — hover for exact values, zoom, pan, and download via the camera icon.

---

### 6 — Isoform Proportions
Faceted pie charts — one per group — showing each selected isoform as a proportion of total gene expression. Ideal for detecting isoform switching between cell types or conditions.

- Selected isoforms are shown as distinct colours; all others are grouped as **Other isoforms** (grey).
- Percentages are displayed for slices > 5%.
- Set **Min counts** in the sidebar to hide groups with low coverage.

> **Isoform switching:** Compare the dominant (largest) slice across groups — a shift in colour indicates a switch in the preferred isoform.

---

### 7 — Expression Trajectory
Plots single-cell normalised isoform expression across an ordered sequence of groups (e.g. developmental stages), with a loess smooth overlaid per isoform panel.

- **Drag and drop** the group labels in the sidebar to define the left-to-right ordering.
- Each isoform gets its own panel with an independent y-axis.
- A loess smooth with a shaded confidence ribbon highlights overall expression trends per isoform.

> **Pairing tip:** Use Proportions to see *which* isoform dominates, and Trajectory to see *when* it changes along the axis.

---

### 8 — Compare Two Isoforms
Directly compare the expression of any two isoforms (or genes) at single-cell resolution across three complementary views:

- **Blend:** a UMAP blend — one isoform in red, the other in blue, co-expressing cells in purple. Works across genes or within a gene.
- **Scatter:** each cell is a dot — X axis = isoform 1 expression, Y axis = isoform 2 expression. Includes a linear fit, R², and optional marginal distributions.
- **Feature Map:** a categorical UMAP classifying every cell as **Both** (co-expressed), **Only X**, **Only Y**, or **Neither** — ideal for mutual-exclusivity analysis.

> **Tip:** Use the **Show cells** filter to restrict comparisons to specific cell groups.

---

## Input Data Format

### Seurat Object (`.rds` / `.qs` / `.qs2`)
- Must contain a gene assay (e.g. `RNA`) and an isoform assay (e.g. `iso`).
- Isoform IDs must follow `ENST00000XXXXX.X-GENENAME` — the gene suffix must match gene assay names (e.g. `ENST00000544301.7-VIM` ↔ `VIM`).
- Requires at least one dimensional reduction (UMAP, PCA).
- Metadata must include a column for cell grouping.

See the [FLAMESv2 Tutorial](https://sefi196.github.io/FLAMESv2_LR_sc_tutorial/) for generating compatible objects.

### GTF Annotation (`.gtf`) — Optional
- Standard GTF format with a `transcript_id` attribute.
- Transcript IDs must match the ENST IDs in your isoform assay.
- Only required for the **Transcript Structure** tab.
- Files up to 5 GB supported locally.

### Performance Tips
- **Online:** Keep uploads under 200 MB for best responsiveness.
- **Local:** Use `qs::qsave()` format for large objects — up to 10x faster loading than `.rds`.
- **qs / qs2:** Both formats are supported; the app auto-detects which library to use.

---

## Installation

### Option 1 — Use Online
Access LongViewSC directly in your browser — no setup required.

```
https://longviewsc.researchsoftware.unimelb.edu.au/LongViewSC/
```

Recommended for files < 200 MB.

### Option 2 — Local Installation
For large files (> 200 MB) or offline use. Requires conda or miniconda.

```bash
git clone https://github.com/Sefi196/LongViewSC.git
cd LongViewSC
conda env create -f environment.yml
conda activate LongViewSC_env
```

Install ggtranscript (not on conda):

```bash
Rscript -e "devtools::install_github('dzhang32/ggtranscript')"
```

Launch the app:

```bash
Rscript -e "shiny::runApp('app.R', launch.browser=TRUE)"
```

Runs entirely locally. No data is sent online.

---

## Contact & Citation

Developed by **Sefi Prawer** — Clark Laboratory, University of Melbourne.
📧 [sefi.prawer@unimelb.edu.au](mailto:sefi.prawer@unimelb.edu.au)

If you use LongViewSC please cite our work.
