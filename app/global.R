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
    library(qs)
    library(qs2)
    library(sortable)   # drag-to-order UI for trajectory tab
  })
})

options(shiny.maxRequestSize = 5000 * 1024^2)  # Increases limit to 50GB

# Helper: load a Seurat object from .rds, .qs, or .qs2 based on file extension.
# `original_name` should be the user-facing filename (e.g. input$seurat_file$name)
# because Shiny stores uploads with a tmp path that has no useful extension.
load_seurat_object <- function(path, original_name = NULL) {
  ext_source <- if (!is.null(original_name)) original_name else path
  ext <- tolower(tools::file_ext(ext_source))
  
  if (ext == "rds") {
    readRDS(path)
  } else if (ext %in% c("qs", "qs2")) {
    # Try qs::qread first; if it fails with a format error, fall back to qs2.
    # Some pipelines save .qs files using the qs2 format (and vice versa).
    tryCatch(
      qs::qread(path),
      error = function(e) {
        if (grepl("format not detected|not a qs file|invalid|magic", conditionMessage(e), ignore.case = TRUE)) {
          tryCatch(
            qs2::qs_read(path),
            error = function(e2) {
              stop("Could not read file as qs or qs2 format. Original error: ", conditionMessage(e))
            }
          )
        } else {
          stop(e)
        }
      }
    )
  } else {
    stop(paste0(
      "Unsupported file format: '.", ext, "'. ",
      "Please upload a .rds, .qs, or .qs2 Seurat object."
    ))
  }
}

# read in demo data 
# 1) Read the demo Seurat object once, at startup
seurat_obj_demo <- readRDS("demo_data/Day55_tutorial_gene_and_isoform_seurat.rds")

# 2) Read the demo GTF once, at startup
gtf_obj_demo <- rtracklayer::import("demo_data/demo_isoform_annotated.gtf") %>% as_tibble()

## Plotting and utils scripts 
source("R/plot_gene_transcripts.R")
source("R/plot_pseudobulk_heatmap.R")
source("R/utils.R")
source("R/plot_isoform_pie.R")
source("R/plot_expression_trajectory.R")