#server side

server <- function(input, output, session) {
  #helper function for dyanmic seletcion of assay options 
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  #code for GPT helping with reset fucntionalltiy of back button 
  # 1) Keep track of “last reset” so that plots only
  #    render when input$GO > lastReset()
  # -----------------------------------------------
  
  lastReset <- reactiveVal(0)
  
  
  #### Main logic ### 
  # Hide everything except landing page initially
  shinyjs::hide("instructionsPage")
  shinyjs::hide("mainUI")
  
  # When 'Start Analysis' is clicked, hide landing page and show main UI
  observeEvent(input$startBtn, {
    shinyjs::hide("landingPage")
    shinyjs::hide("instructionsPage")  # Ensure instructions page is hidden
    shinyjs::show("mainUI")
    shinyjs::show("seurat_file")
    shinyjs::show("gtf")
    shinyjs::show("resetBtn")
  })
  
  # When the mailbox button is clicked, show/hide the dropdown menu
  observeEvent(input$mailbox, {
    toggle("dropdownMenu")
  })
  
  # When 'Instructions' is clicked, hide landing page and show instructions
  observeEvent(input$instructionsBtn, {
    shinyjs::hide("landingPage")
    shinyjs::hide("mainUI")  # Ensure main UI is hidden
    shinyjs::show("instructionsPage")
  })
  
  # Add functionality for back button from instructions to landing page
  # Back button logic (conditional on demo_mode)
  observeEvent(input$backToLanding, {
    # Always hide the mainUI and show the landing page
    shinyjs::hide("mainUI")
    shinyjs::show("landingPage")
    shinyjs::hide("instructionsPage")
    
    # If we are in demo_mode, then drop demo data and reset everything
    if (app_state$demo_mode) {
      # 1) Clear the stored Seurat + GTF
      seurat_obj(NULL)
      gtf(NULL)
      app_state$demo_mode <- FALSE
      
      # 2) Reset all the form inputs (so that feature/selectInput/etc. go blank)
      shinyjs::reset("analysisForm")
      
      # 3) Remove any warning message
      output$isoforms_warning <- renderUI({ NULL })
      
      # 4) Force‐blank every plot by bumping lastReset()
      lastReset(input$GO)
    }
  })
  
  
  
  observeEvent(input$resetBtn, {
    # 1) Record the current GO count so that all existing plots disappear
    lastReset(input$GO)
    
    # 2) Clear every selector inside #analysisForm
    shinyjs::reset("analysisForm")
    
    # 3) Clear any custom warning
    output$isoforms_warning <- renderUI({ NULL })
    
    shinyjs::reset("seurat_file")
    shinyjs::reset("gtf")
  })
  
  # Demo Button functionality
  observeEvent(input$DemoBtn, {
    app_state$demo_mode <- TRUE  # <--- Track that demo data was loaded
    
    shinyjs::hide("seurat_file")
    shinyjs::hide("gtf")
    shinyjs::hide("resetBtn")
    
    # Show spinner
    shinyjs::show("spinner")  # Assuming you have a div with id="spinner"
    # Hide spinner after loading data
    
    # Update the reactive value with the demo Seurat object
    seurat_obj(seurat_obj_demo)
    
    # Dynamically update dropdown menus for reductions, assays, and features after demo Seurat object is loaded
    observe({
      req(seurat_obj())  # Ensure the Seurat object is available before updating UI
      print("Updating UI based on Demo Seurat object.")  # Debugging print
      
      # Get available reductions, assays, and group_by
      reductions <- names(seurat_obj()@reductions)
      isoform_assay <- names(seurat_obj()@assays)
      group_by <- colnames(seurat_obj()@meta.data)[!sapply(seurat_obj()@meta.data, is.numeric)]
      
      # Logic to select defaults
      default_isoform_assay <- isoform_assay[grepl("^iso(form)?$", isoform_assay, ignore.case = TRUE)][1]
      default_group_by <- group_by[grepl("seurat_clusters|cell.*type|annotation", group_by, ignore.case = TRUE)][1]
      default_reduction <- reductions[grepl("umap|tsne|iso", reductions, ignore.case = TRUE)][1]
      
      # Update UI elements
      updateSelectInput(session, "reduction", choices = reductions, selected = default_reduction %||% reductions[1])
      updateSelectInput(session, "isoform_assay", choices = isoform_assay, selected = default_isoform_assay %||% isoform_assay[1])
      updateSelectInput(session, "group_by", choices = group_by, selected = default_group_by %||% group_by[1])
      updateSelectizeInput(session, "feature", choices = rownames(seurat_obj()), server = TRUE)
    })
    
    # Update the reactive value with the new GTF object
    gtf(gtf_obj_demo)
    
    # Show notification when demo data is loaded
    showNotification("Demo data loaded successfully! Ready for analysis.", type = "message", duration = 5)
    
    # Show the main UI after demo data is loaded
    shinyjs::hide("landingPage")
    shinyjs::hide("spinner")
    shinyjs::show("mainUI")
  })
  
  #set the state of the app to control the demo mode resetting the input objects
  app_state <- reactiveValues(
    demo_mode = FALSE
  )
  
  #### Main logic ####
  # should wrap this up in a separate script 
  seurat_obj <- reactiveVal(NULL)
  gtf <- reactiveVal(NULL)
  
  # Handle Seurat object file upload
  observeEvent(input$seurat_file, {
    req(input$seurat_file)  # Ensure a file is uploaded
    
    shinyjs::show("spinner")
    obj <- readRDS(input$seurat_file$datapath)  # Load the Seurat object
    
    #set app status 
    app_state$demo_mode <- FALSE
    # Update the reactive value with the new Seurat object
    seurat_obj(obj)
    
    # Dynamically update dropdown menus for reductions, assays, and features after Seurat object is loaded
    observe({
      req(seurat_obj())  # Ensure the Seurat object is available before updating UI
      print("Updating UI based on Seurat object.")  # Debugging print
      
      # Get available reductions, assays, and group_by
      reductions <- names(seurat_obj()@reductions)
      isoform_assay <- names(seurat_obj()@assays)
      group_by <- colnames(seurat_obj()@meta.data)[!sapply(seurat_obj()@meta.data, is.numeric)]
      
      print(paste("Reductions available:", paste(reductions, collapse = ", ")))
      print(paste("Isoform assays available:", paste(isoform_assay, collapse = ", ")))
      print(paste("Group_by columns:", paste(group_by, collapse = ", ")))
      
      # Logic to select defaults
      default_isoform_assay <- isoform_assay[grepl("^iso(form)?$", isoform_assay, ignore.case = TRUE)][1]
      default_group_by <- group_by[grepl("seurat_clusters|cell.*type|annotation", group_by, ignore.case = TRUE)][1]
      default_reduction <- reductions[grepl("umap|tsne|iso", reductions, ignore.case = TRUE)][1]
      
      # Update UI elements
      updateSelectInput(session, "reduction", choices = reductions, selected = default_reduction %||% reductions[1])
      updateSelectInput(session, "isoform_assay", choices = isoform_assay, selected = default_isoform_assay %||% isoform_assay[1])
      updateSelectInput(session, "group_by", choices = group_by, selected = default_group_by %||% group_by[1])
      updateSelectizeInput(session, "feature", choices = rownames(seurat_obj()), server = TRUE)
      
      shinyjs::hide("spinner")
      
      if (!is.null(default_isoform_assay)) {
        showNotification(paste("Auto-selected assay:", default_isoform_assay), type = "message", duration = 4)
      }
      
      if (!is.null(default_group_by)) {
        showNotification(paste("Auto-selected metadata column:", default_group_by), type = "message", duration = 4)
      }
    })
  })
  
  # Handle GTF file upload
  observeEvent(input$gtf, {
    req(input$gtf)  # Ensure a file is uploaded
    shinyjs::show("spinner")
    gtf_path <- input$gtf$datapath  # Get the uploaded file path
    
    # Import GTF file using rtracklayer
    gtf_obj <- rtracklayer::import(gtf_path) %>% as_tibble()
    
    # Update the reactive value with the new GTF object
    gtf(gtf_obj)
    
    shinyjs::hide("spinner")
    #set app status
    app_state$demo_mode <- FALSE
  })
  
  
  # Reactive expression triggered by "GO" button (updates everything)
  # Extract features related to isoforms
  # NEW: reactive version—auto‐invalidates when input$feature or input$GO is missing
  filtered_data <- eventReactive(input$GO, {
    req(seurat_obj(), input$feature)
    
    assay_name     <- input$isoform_assay
    isoform_feats  <- rownames(seurat_obj()@assays[[assay_name]]@features)
    joined         <- JoinLayers(seurat_obj())
    expr_matrix    <- GetAssayData(joined, assay = assay_name, slot = "data")
    
    matching_feats <- grep(
      paste0("(^|-|\\b)", input$feature, "($|\\b)"),
      isoform_feats, value = TRUE
    )
    subset_expr    <- expr_matrix[matching_feats, , drop = FALSE]
    total_expr     <- Matrix::rowSums(subset_expr)
    matching_feats <- names(sort(total_expr, decreasing = TRUE))
    
    list(
      celltype_plot     = DimPlot(seurat_obj(), reduction = input$reduction, group.by = input$group_by),
      feature_plot_gene = FeaturePlot(seurat_obj(), features = input$feature, reduction = input$reduction),
      vln_plot          = VlnPlot(seurat_obj(), features = input$feature, group.by = input$group_by),
      isoform_features  = matching_feats
    )
  })
  
  # Reactive function to get the top N isoform features
  # provide warning if too many isoforms selected 
  isoform_features_to_plot <- eventReactive(input$GO, {
    req(filtered_data(), input$number_of_isoforms)
    
    all_isoforms <- filtered_data()$isoform_features
    gene_name    <- input$feature
    num_isoforms <- length(all_isoforms)
    
    if (input$number_of_isoforms > num_isoforms) {
      msg <- paste(gene_name, "has", num_isoforms, "isoforms.")
      output$isoforms_warning <- renderUI({
        HTML(
          paste0(
            "<p style='color:red; font-size:20px; text-align:center;'>",
            msg,
            "</p>"
          )
        )
      })
    } else {
      output$isoforms_warning <- renderUI({ NULL })
    }
    
    head(all_isoforms, input$number_of_isoforms)
  })
  
  # Reactive function for the isoform Feature Plot
  isoform_plot <- reactive({
    req(isoform_features_to_plot())   # auto‐invalidates if isoform_features_to_plot() is NULL
    plots <- FeaturePlot(
      seurat_obj(),
      features  = isoform_features_to_plot(),
      reduction = input$reduction,
      order     = TRUE
    )
    plots <- lapply(plots, function(pl) pl + theme(plot.title = element_text(size = 11)))
    wrap_plots(plots = plots, ncol = 4)
  })
  
  dotplot_isoform <- reactive({
    req(isoform_features_to_plot())
    DotPlot(
      seurat_obj(),
      features = isoform_features_to_plot(),
      assay    = input$isoform_assay,
      group.by = input$group_by
    ) + theme(axis.text.x = element_text(angle = 80, hjust = 1))
  })
  
  
  # Reactive expression to trigger the heatmap plot when the button is clicked
  reactive_heatmap <- reactive({
    req(isoform_features_to_plot())
    plot_pseudobulk_heatmap(
      seurat_obj       = seurat_obj(),
      group.by         = input$group_by,
      isoforms_to_plot = isoform_features_to_plot(),
      isoform_assay    = input$isoform_assay
    )
  })
  
  #### Render the main plots (with the new req) ####
  
  # Gene Feature Plot
  output$feature_plot_gene <- renderPlot({
    req(input$GO > lastReset(), filtered_data())
    filtered_data()$feature_plot
  })
  
  # Violin Plot
  output$vln_plot <- renderPlot({
    req(input$GO > lastReset(), filtered_data())
    filtered_data()$vln_plot
  })
  
  # Cell Type DimPlot
  output$celltype_plot <- renderPlot({
    req(input$GO > lastReset(), filtered_data())
    filtered_data()$celltype_plot
  })
  
  # Isoform Feature Plot
  output$feature_plot_iso <- renderPlot({
    req(input$GO > lastReset(), isoform_plot())
    isoform_plot()
  })
  
  # Dot Plot for Isoforms
  output$dot_plot_iso <- renderPlot({
    req(input$GO > lastReset(), dotplot_isoform())
    dotplot_isoform()
  })
  
  # Transcript Structure (from GTF)
  output$transcript_plot <- renderPlot({
    req(input$GO > lastReset(), gtf(), isoform_features_to_plot())
    if (is.null(gtf()) || nrow(gtf()) == 0) {
      showModal(modalDialog(
        title = "Missing GTF File",
        "Please provide a valid GTF to build the isoform stack.",
        easyClose = TRUE,
        footer    = NULL
      ))
    } else {
      plot_gene_transcripts(
        isoforms_to_plot = isoform_features_to_plot(),
        gtf              = gtf()
      )
    }
  })
  
  # Pseudobulk Heatmap (Plotly)
  output$heatmap_plot <- renderPlotly({
    req(input$GO > lastReset(), reactive_heatmap())
    reactive_heatmap()
  })
  

  # —— download handlers must be guarded ——  
  
  downloadModalServer(
    "celltype_plot",
    reactive({
      req(req(input$GO > lastReset()))
      filtered_data()$celltype_plot
    }),
    "Cell type_plot"
  )
  
  downloadModalServer(
    "feature_plot_gene",
    reactive({
      req(req(input$GO > lastReset()))
      filtered_data()$feature_plot
    }),
    "Feature plot - Gene"
  )
  
  downloadModalServer(
    "vln_plot",
    reactive({
      req(req(input$GO > lastReset()))
      filtered_data()$vln_plot
    }),
    "Violin plot"
  )
  
  downloadModalServer(
    "feature_plot_isoform",
    reactive({
      req(req(input$GO > lastReset()))
      isoform_plot()
    }),
    "Feature plot - Isoform"
  )
  
  downloadModalServer(
    "dot_plot_isoform",
    reactive({
      req(req(input$GO > lastReset()))
      dotplot_isoform()
    }),
    "Dot plot - Isoform"
  )
  
  downloadModalServer(
    "Isoform_TranscriptStructure",
    reactive({
      req(req(input$GO > lastReset()))
      plot_gene_transcripts(
        isoforms_to_plot = isoform_features_to_plot(),
        gtf = gtf()
      )
    }),
    "Transcript Structure"
  )
  
  downloadPlotlyModalServer(
    "pseudobulk_heatmap",
    reactive({
      req(req(input$GO > lastReset()))
      reactive_heatmap()
    }),
    "Pseudobulk heatmap"
  )
  
}