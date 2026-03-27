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
    library(rmarkdown)
    library(shinyjs)
    library(shinydashboard)
    library(patchwork)
    library(DT)
    library(qs)
    library(qs2)
    library(sortable)      # drag-to-order UI for trajectory tab
    library(plotly)        # Sankey and interactive plots
    library(shinyWidgets)  # material toggle switches
  })
})

options(shiny.maxRequestSize = 5000 * 1024^2)  # Increases limit to 50GB

# Global ggplot2 defaults — larger, bolder, sans-serif text for all plots
theme_update(
  text         = element_text(family = "sans"),
  plot.title   = element_text(size = 16, face = "bold", family = "sans"),
  axis.title   = element_text(size = 14, face = "bold", family = "sans"),
  axis.text    = element_text(size = 13, family = "sans"),
  strip.text   = element_text(size = 13, face = "bold", family = "sans"),
  legend.title = element_text(size = 13, face = "bold", family = "sans"),
  legend.text  = element_text(size = 12, family = "sans")
)

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

# Discover Clark Lab datasets from the clark_data/ directory.
# Each subdirectory should contain one Seurat file (.rds/.qs/.qs2) and
# optionally one GTF file (.gtf). The folder name becomes the dataset label.
clark_data_dir <- "clark_data"
clark_datasets <- list()
if (dir.exists(clark_data_dir)) {
  dirs <- list.dirs(clark_data_dir, full.names = TRUE, recursive = FALSE)
  for (d in dirs) {
    dataset_name <- basename(d)
    seurat_files <- list.files(d, pattern = "\\.(rds|qs|qs2)$", full.names = TRUE, ignore.case = TRUE)
    gtf_files    <- list.files(d, pattern = "\\.gtf$",           full.names = TRUE, ignore.case = TRUE)
    if (length(seurat_files) > 0) {
      clark_datasets[[dataset_name]] <- list(
        seurat_path = seurat_files[1],
        gtf_path    = if (length(gtf_files) > 0) gtf_files[1] else NULL
      )
    }
  }
}

# Clark Lab dataset password (set CLARK_PASSWORD env var to override)
clark_password <- Sys.getenv("CLARK_PASSWORD", unset = "clark2024")

## Plotting and utils scripts
source("R/plot_gene_transcripts.R")
source("R/plot_pseudobulk_heatmap.R")
source("R/utils.R")
source("R/plot_isoform_pie.R")
source("R/plot_expression_trajectory.R")
# --- Usage logging ---
suppressPackageStartupMessages({
  library(DBI)
  library(RSQLite)
})
APP_VERSION <- "dev"
source("R/usage_logger.R")
