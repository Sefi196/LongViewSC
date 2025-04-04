<h1 style="text-align: center;">

Isoform Explorer - Welcome Guide

</h1>

<h1 style="text-align: center;">

## Overview

</h2>

Isoform Explorer is a Shiny app designed to help researchers visualize isoform expression from single-cell sequencing data. This tool enables easy comparison of isoform-level expression across different conditions and clusters using heatmaps and transcript structure visualizations.

<https://sefi196.github.io/FLAMESv2_LR_sc_tutorial/>

## What You Can Do

-   **Pseudobulk Heatmap**: Aggregate and visualize isoform expression across clusters or conditions.\
-   **Transcript Structure Plot**: View exon-intron structures of selected isoforms.\
-   **Consistent Isoform Ordering**: Ensure alignment between heatmaps and transcript structure plots.\
-   **Dynamic Selection**: Choose specific isoforms or all isoforms of a gene for analysis.

## How to Use the App

1.  **Select a Dataset**
    -   The app loads a pre-processed Seurat object containing isoform-level expression data.\
    -   You can choose between different assays and clustering strategies.
2.  **Choose Isoforms to Plot**
    -   Search for a gene or manually select isoforms of interest.\
    -   The heatmap and transcript plots will update accordingly.
3.  **Adjust Plot Settings**
    -   Customize the heatmap dendrogram and ordering.\
    -   Toggle display options for better visualization.
4.  **Interpret Results**
    -   The heatmap shows **pseudobulk log-transformed CPM values** across clusters.\
    -   The transcript plot provides a **structural view of isoform architecture**.

💡 *Tip: Click on different isoforms to explore expression differences across clusters!*

# Single Cell Isoform Expression Viewer

## Overview

The Single Cell Isoform Expression Viewer is an interactive **Shiny** application designed to visualize isoform expression from single-cell RNA sequencing (scRNA-seq) data. This tool supports **Seurat objects** and allows users to explore gene and isoform-level expression patterns with a variety of plots, including feature plots, violin plots, dot plots, and pseudobulk heatmaps.

## Features

-   **Upload and process Seurat objects (.rds) and GTF annotation files (.gtf).**
-   **Visualize gene expression** using feature plots and violin plots.
-   **Explore isoform expression** at the single-cell level.
-   **Analyze transcript structure** using isoform-specific plots.
-   **Generate heatmaps** to compare isoform expression across conditions.

## Installation and Requirements

The application requires the following R packages:

``` r
library(shiny)
library(Seurat)
library(ggplot2)
library(ggtranscript)
library(rtracklayer)
library(heatmaply)
library(tibble)
library(markdown)
library(shinyjs)
library(shinydashboard)
```

To run the app, ensure these dependencies are installed using:

``` r
install.packages(c("shiny", "ggplot2", "heatmaply", "tibble", "markdown", "shinyjs", "shinydashboard"))
BiocManager::install(c("Seurat", "ggtranscript", "rtracklayer"))
```

## Running the Application

Launch the app by running the following command in R:

``` r
shiny::runApp("path/to/app")
```

Ensure that the required source files (`plot_gene_transcripts.R` and `plot_pseudobulk_heatmap.R`) are available in the working directory.

## File Formats

### Input Files

-   **Seurat Object (.rds)**: A serialized Seurat object containing expression data, metadata, and dimensional reductions.
-   **GTF File (.gtf)**: A gene annotation file used to map isoforms to gene symbols.

### Output Interpretation

-   **Gene Feature Plot**: Displays the spatial distribution of gene expression in the selected dimensional reduction (e.g., UMAP, tSNE).
-   **Cell Type Plot**: Clusters of cells based on metadata annotations.
-   **Violin Plot**: Compares gene expression across metadata groups.
-   **Isoform Feature Plot**: Visualizes the expression of isoforms.
-   **Dot Plot**: Summarizes isoform expression patterns across groups.
-   **Transcript Structure Plot**: Shows exon/intron structures of different isoforms.
-   **Pseudobulk Heatmap**: Provides an overview of isoform expression across conditions.

## Usage

1.  **Upload Data**
    -   Load a Seurat `.rds` file.
    -   Upload a `.gtf` file for annotation.
2.  **Select Parameters**
    -   Choose a gene to analyze.
    -   Select a dimensional reduction method.
    -   Specify the isoform assay.
    -   Define the metadata column to group cells.
    -   Set the number of isoforms to display.
3.  **Generate Plots**
    -   Click **GO** to visualize selected features.
4.  **Interpret Results**
    -   Explore expression patterns and isoform distributions.

## Contact and Support

For further inquiries or support, please contact [**genomeprot\@outlook.com**](mailto:genomeprot@outlook.com) or visit our [Clark Laboratory website](https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/clark-lab).

------------------------------------------------------------------------

**Developed by Sefi Prawer in the Clark Laboratory, University of Melbourne**
