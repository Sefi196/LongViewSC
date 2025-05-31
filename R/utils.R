# A helper to register a downloadHandler for a given reactive plot and ID
# Helper function
downloadModalServer <- function(id, plot_expr, prefix) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input[[paste0("download_", id)]], {
      showModal(modalDialog(
        title = paste("Download:", gsub("_", " ", prefix)),
        
        # Use inline-flex and smaller styling for numeric inputs
        
        # Stack Width and Height vertically
        div(style = "display: flex; flex-direction: column; gap: 8px; margin-bottom: 10px;",
            numericInput(ns(paste0("width_", id)), "Width (in)", 8, min = 4, max = 20, width = "180px"),
            numericInput(ns(paste0("height_", id)), "Height (in)", 6, min = 4, max = 20, width = "180px")
        ),
        # Compact download buttons
        div(style = "display: flex; flex-wrap: wrap; gap: 6px; justify-content: space-between;",
            downloadButton(ns(paste0("download_", id, "_pdf")), "PDF", style = "padding: 6px 10px; font-size: 85%;"),
            downloadButton(ns(paste0("download_", id, "_png")), "PNG", style = "padding: 6px 10px; font-size: 85%;"),
            #downloadButton(ns(paste0("download_", id, "_tiff")), "TIFF", style = "padding: 6px 10px; font-size: 85%;"),
            downloadButton(ns(paste0("download_", id, "_jpeg")), "JPEG", style = "padding: 6px 10px; font-size: 85%;")
        ),
        
        easyClose = TRUE,
        footer = NULL,
        size = "s"  # smallest modal size available
      ))
    })
    
    # Download handlers
    output[[paste0("download_", id, "_pdf")]] <- downloadHandler(
      filename = function() paste0(prefix, "_", Sys.Date(), ".pdf"),
      content = function(file) {  
        ggsave(file, plot = plot_expr(), width = input[[paste0("width_", id)]], height = input[[paste0("height_", id)]], device = "pdf")
      }
    )
    
    output[[paste0("download_", id, "_png")]] <- downloadHandler(
      filename = function() paste0(prefix, "_", Sys.Date(), ".png"),
      content = function(file) {  
        ggsave(file, plot = plot_expr(), width = input[[paste0("width_", id)]], height = input[[paste0("height_", id)]], device = "png", dpi = 300)
      }
    )
    
    #output[[paste0("download_", id, "_tiff")]] <- downloadHandler(
     # filename = function() paste0(prefix, "_", Sys.Date(), ".tiff"),
      #content = function(file) {
      #  ggsave(file, plot = plot_expr(), width = input[[paste0("width_", id)]], height = input[[paste0("height_", id)]], device = "tiff", dpi = 300)
      #}
    #)
    
    output[[paste0("download_", id, "_jpeg")]] <- downloadHandler(
      filename = function() paste0(prefix, "_", Sys.Date(), ".jpeg"),
      content = function(file) {  
        ggsave(file, plot = plot_expr(), width = input[[paste0("width_", id)]], height = input[[paste0("height_", id)]], device = "jpeg", dpi = 300)
      }
    )
  })
}

#download for heatmaply plots
downloadPlotlyModalServer <- function(id, plot_expr, prefix) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input[[paste0("download_", id)]], {
      showModal(modalDialog(
        title = paste("Download:", gsub("_", " ", prefix)),
        div(style = "margin-bottom: 10px;",
            helpText("This is an interactive Plotly plot. You can download an HTML file here or as a PNG using the camera icon in the top-right corner of the plot.")
        ),
        downloadButton(ns(paste0("download_", id, "_html")), "Download HTML",
                       style = "padding: 6px 10px; font-size: 85%;"),
        easyClose = TRUE,
        footer = NULL,
        size = "s"
      ))
    })
    
    output[[paste0("download_", id, "_html")]] <- downloadHandler(
      filename = function() paste0(prefix, "_", Sys.Date(), ".html"),
      content = function(file) {  
        htmlwidgets::saveWidget(plot_expr(), file = file, selfcontained = TRUE)
      }
    )
  })
}

