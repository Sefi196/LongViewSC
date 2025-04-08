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
  })
})

options(shiny.maxRequestSize = 2000 * 1024^2)  # Increases limit to 20GB

source("R/plot_gene_transcripts.R")
source("R/plot_pseudobulk_heatmap.R")
source("R/utils.R")
