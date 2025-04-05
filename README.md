<h1 style="text-align: center;">

LongViewSC: Interactive Visualization of Gene and Isoform Expression in Single-Cell Data

</h1>

<h1 style="text-align: center;">

## Overview

</h2>

**LongViewSC** is an interactive Shiny app designed to visualize both gene and isoform expression in Long read single-cell data. This user-friendly tool is built to provide accessible and intuitive visualization for researchers, regardless of their coding experience. This user-friendly interface allows for intuitive exploration of your data with various visualization tools such as FeaturePlots, heatmaps and transcript structure views. With these tools researchers can easily compare gene and isoform expression across different conditions, clusters and cell types.

If you find this application helpful, please cite our work: {PLACEHOLDER}.

<https://sefi196.github.io/FLAMESv2_LR_sc_tutorial/>

## File Formats

### Input Files

-   **Seurat Object (.rds)**: A serialized Seurat object containing expression data, metadata, and dimensional reductions.
-   **GTF File (.gtf)**: A gene annotation file used to map isoforms to gene symbols.

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

## Overview

The Single Cell Isoform Expression Viewer is an interactive **Shiny** application designed to visualize isoform expression from single-cell RNA sequencing (scRNA-seq) data. This tool supports **Seurat objects** and allows users to explore gene and isoform-level expression patterns with a variety of plots, including feature plots, violin plots, dot plots, and pseudobulk heatmaps.

## Features

-   **Upload and process Seurat objects (.rds) and GTF annotation files (.gtf).**
-   **Visualize gene expression** using feature plots and violin plots.
-   **Explore isoform expression** at the single-cell level.
-   **Analyze transcript structure** using isoform-specific plots.
-   **Generate heatmaps** to compare isoform expression across conditions.

## Installation and Requirements

### **Option 1: Access LongViewSC online** (currently under development)  

<https://sefi196.shinyapps.io/sc_expresstion_view/>

For optimal performance, we recommend that users upload files smaller than 200MB. For those who wish to explore larger files, please refer to **Option 2** for alternative installation options.

### **Option 2 : Local installation** 

For users exploring large file \> 200MB or many files sequentially uploads may be slow. To ovecome this, a local intalltion is reccomeded.

The application has dependencies that we have provided as a conda environment file. to install and run the app:

1.  **Install conda or miniconda**

2.  **Clone this repo:**

``` r
git clone https://github.com/Sefi196/LongViewSC.git 
cd LongViewSC
```

3.  **Run the installing script**

``` r
cd LongViewSC
chmod +x install.sh
./install.sh
```

4.  **Launch the app by running the following command in R:**

``` r
shiny::runApp("/path/to/your/ShinyApp")
```

This will launch the app in a web browser.

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
    -   Specify the isoform assay. **Don't forget to do this!**
    -   Define the metadata column to group cells.
    -   Set the number of isoforms to display.
3.  **Generate Plots**
    -   Click **GO** to visualize selected features.
4.  **Interpret Results**
    -   Explore expression patterns and isoform structure
5.  **Download your data**

## Contact and Support

For further inquiries or support, please contact sefi.prawer\@unimelb.edu.au or leave a comment on the github page

------------------------------------------------------------------------

**Developed by Sefi Prawer in the Clark Laboratory, University of Melbourne**
