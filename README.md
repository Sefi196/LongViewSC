# LongViewSC

**Interactive visualisation of gene and isoform expression in long-read single-cell data.**

LongViewSC is a Shiny application that lets researchers explore gene and isoform expression from long-read scRNA-seq data without writing any code. Upload a Seurat object, select a gene, and navigate six visualisation tabs.

## Try it online

**<https://longviewsc.researchsoftware.unimelb.edu.au/LongViewSC/>**

> Files under ~200 MB load quickly online. For larger datasets use the local installation below.

## Visualisation tabs

| Tab | What it shows |
|-----|---------------|
| **Overview** | UMAP coloured by group, gene feature plot, violin plot |
| **Isoform Statistics** | Searchable table of all detected isoforms for a gene |
| **Isoform Expression** | Per-isoform feature plots and violin plots |
| **Transcript Structure** | Exon-intron architecture from GTF (requires GTF upload) |
| **Pseudobulk Heatmap** | Interactive clustered heatmap of log(CPM+1) isoform expression |
| **Proportions** | Faceted pie charts showing isoform fractions per group |
| **Trajectory** | Expression across an ordered sequence of groups with loess smooth |

## Input data

### Seurat Object (`.rds` / `.qs` / `.qs2`)
- Must contain a **gene assay** (e.g. `RNA`) and an **isoform assay** (e.g. `iso`).
- Isoform features must follow the naming convention: `ENST00000XXXXX.X-GENENAME`  
  (e.g. `ENST00000544301.7-VIM`).
- Requires at least one dimensional reduction (UMAP, PCA, etc.).
- Metadata must include a column for cell grouping (cluster, cell type, condition).

### GTF File (`.gtf`) — optional
- Standard GTF with a `transcript_id` attribute matching your isoform assay IDs.
- Only needed for the **Transcript Structure** tab.

> See the [FLAMESv2 Long-Read Single-Cell Tutorial](https://sefi196.github.io/FLAMESv2_LR_sc_tutorial/) for instructions on generating compatible files.

## Local installation

### 1. Install conda or miniconda

See <https://www.anaconda.com/docs/getting-started/miniconda/install>.

### 2. Clone the repository

```bash
git clone https://github.com/Sefi196/LongViewSC.git
cd LongViewSC
```

### 3. Create and activate the conda environment

```bash
conda env create -f environment.yml
conda activate LongViewSC_env
```

This installs all R dependencies except `ggtranscript`, which must be installed from GitHub.

### 4. Install ggtranscript

```r
Rscript -e "devtools::install_github('dzhang32/ggtranscript')"
```

### 5. Launch the app

```bash
Rscript -e "shiny::runApp('app.R', launch.browser=TRUE)"
```

The app opens in your default browser and runs entirely locally — no data is sent online.

## R package dependencies

| Package | Source |
|---------|--------|
| shiny, shinyjs, shinydashboard | CRAN |
| Seurat | CRAN |
| ggplot2, patchwork, scales | CRAN |
| plotly, heatmaply | CRAN |
| DT | CRAN |
| tibble, dplyr, tidyr | CRAN |
| markdown, sortable | CRAN |
| qs, qs2 | CRAN |
| rtracklayer | Bioconductor |
| ggtranscript | GitHub (dzhang32/ggtranscript) |

## Citation

If you use LongViewSC in your research, please cite our work: *(citation placeholder — to be updated on publication)*.

## Contact & support

- Email: [sefi.prawer@unimelb.edu.au](mailto:sefi.prawer@unimelb.edu.au)
- GitHub issues: <https://github.com/Sefi196/LongViewSC/issues>
- Lab: [Clark Laboratory, University of Melbourne](https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/clark-lab)

---

**Developed by Sefi Prawer at the Clark Laboratory, University of Melbourne.**
