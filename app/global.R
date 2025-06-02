suppressPackageStartupMessages({
  suppressMessages({
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
  library(patchwork) 
  library(DT)
  })
})

options(shiny.maxRequestSize = 5000 * 1024^2)  # Increases limit to 50GB

# read in demo data 
# 1) Read the demo Seurat object once, at startup
seurat_obj_demo <- readRDS("demo_data/Day55_tutorial_gene_and_isoform_seurat.rds")

# 2) Read the demo GTF once, at startup
gtf_obj_demo <- rtracklayer::import("demo_data/demo_isoform_annotated.gtf") %>% as_tibble()

## Plotting and utils scripts 
source("R/plot_gene_transcripts.R")
source("R/plot_pseudobulk_heatmap.R")
source("R/utils.R")
