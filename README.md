<h1 style="text-align: center;">

LongViewSC: Interactive Visualization of Gene and Isoform Expression in Single-Cell Data

</h1>

<h1 style="text-align: center;">

## Overview

</h2>

**LongViewSC** is an interactive Shiny app designed to visualize both gene and isoform expression in Long read single-cell data. This user-friendly tool is built to provide accessible and intuitive visualization for researchers, regardless of their coding experience and allow for intuitive exploration of your data. The application provides various visualization tools to help researchers easily compare gene and isoform expression across different conditions, clusters and cell types.

### What You Can Do

-   **Visualize gene expression** using feature plots and violin plots.
-   **Explore isoform expression** at the single-cell level.
-   **Analyze transcript structure** View exon-intron structures of selected isoforms.
-   **Pseudobulk Heatmap**: Aggregate and visualize isoform expression across clusters or conditions.
-   **Dynamic Selection**: Choose specific isoforms or all isoforms of a gene for analysis.

**If you find this application helpful, please cite our work: {PLACEHOLDER}.**

## How to Use the App

1.  **Upload Data**
    -   Load a Seurat `.rds` file (required).
    -   Upload a `.gtf` file for annotation (optional).
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
    -   (Under development)

## File Formats

### Input Files

-   A **Seurat Object (.rds)**: A Seurat object containing gene and isoform expression data, metadata, and dimensional reductions.
    -   The Seurat object **must contain** a gene and isoform assay.
    -   **Isoforms** must be labeled in the following format: `isoformID-geneID` (e.g., `ENST00000544301.7-VIM`).
-   **GTF File (.gtf)**: A gene annotation file.

📌 **Note:** If you are unsure about the format of your sample files, please refer to this comprehensive tutorial on generating these files from gene and count matrices: [FLAMESv2 Long-Read Single-Cell Tutorial](https://sefi196.github.io/FLAMESv2_LR_sc_tutorial/)

## Installation and Requirements

### **Option 1: Access LongViewSC online** (currently under development)

<http://portal.unimelb-lrscview.cloud.edu.au:3838/LongViewSC/>

For optimal performance, we recommend that users upload files smaller than 200MB. For those who wish to explore larger files, please refer to **Option 2** for alternative installation options.

### **Option 2 : Local installation**

For users exploring large file \> 200MB or many files sequentially. In this case uploads may be slow and to overcome this, a local installation is recommended. The application has dependencies that we have provided as a conda environment file.

1.  **Install conda or miniconda**

Installation of conda or miniconda will depend on your system. see <https://www.anaconda.com/docs/getting-started/miniconda/install> for details.

2.  **Clone this repo:**

``` bash
git clone https://github.com/Sefi196/LongViewSC.git 
cd LongViewSC
```

3.  **Create and laod the conda environment**

``` bash
conda env create -f environment.yml
conda activate LongViewSC_env
```

4.  **Install ggtranscript in R**

``` r
#open R
R
# Ensure devtools has been installed by runing 
library("devtools")

#If it has not been loaded install devtools
install.packages("devtools")
#install ggtranscript
devtools::install_github("dzhang32/ggtranscript")

# exit R
q()
```

5.  **Launch the app by running the following command in R:**

``` r
Rscript -e "shiny::runApp('<path_to_app.R>', launch.browser = TRUE)"
```

This will launch the app in your default web browser.

📌 **Note:** Although the app opens in your browser, it runs **entirely locally** on your machine. No internet connection is required, and no data is sent online.

## Contact and Support

For further inquiries or support, please contact sefi.prawer\@unimelb.edu.au or leave a comment on the github page

------------------------------------------------------------------------

**Developed by Sefi Prawer in the Clark Laboratory, University of Melbourne**
