#server side

server <- function(input, output, session) {
  #helper function for dyanmic seletcion of assay options 
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  #code for GPT helping with reset fucntionalltiy of back button 
  # 1) Keep track of “last reset” so that plots only
  #    render when input$GO > lastReset()
  # -----------------------------------------------
  
  lastReset <- reactiveVal(0)
  
  # 1) A reactiveValues container to hold the actual isoform names that are checked.
  #    Initialize it empty; we'll populate it after filtered_data() runs.
  rv <- reactiveValues(
    selected_isoform_names = character(0)
  )
  
  # Always reset chip selection to top 4 whenever a new gene is loaded
  observeEvent(filtered_data(), {
    feats <- isoform_features_to_plot()
    if (length(feats) > 0) {
      rv$selected_isoform_names <- feats[seq_len(min(4, length(feats)))]
    }
  }, ignoreInit = TRUE)
  
  # Selection is now chip-only — isoform_table_rows_selected observer removed
  
  # Chip click — toggle the clicked isoform in/out of selected set
  observeEvent(input$chip_clicked, {
    id  <- input$chip_clicked
    cur <- rv$selected_isoform_names
    if (id %in% cur) {
      rv$selected_isoform_names <- cur[cur != id]
    } else {
      rv$selected_isoform_names <- c(cur, id)
    }
  })
  
  
  #### Main logic ### 
  # Hide everything except landing page initially
  shinyjs::hide("instructionsPage")
  shinyjs::hide("mainUI")
  
  # When 'Start Analysis' is clicked, hide landing page and show main UI
  observeEvent(input$startBtn, {
    shinyjs::hide("landingPage")
    shinyjs::hide("instructionsPage")  # Ensure instructions page is hidden
    shinyjs::show("mainUI")
    shinyjs::show("file_inputs_panel")
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
      join_and_store(NULL)
      gtf(NULL)
      app_state$demo_mode <- FALSE
      
      # 2) Reset all the form inputs (so that feature/selectInput/etc. go blank)
      shinyjs::reset("analysisForm")
      
      # 3) Remove any warning message
      #output$isoforms_warning <- renderUI({ NULL })
      
      # 4) Force‐blank every plot by bumping lastReset()
      lastReset(input$GO)
    }
  })
  
  
  
  observeEvent(input$resetBtn, {
    lastReset(input$GO)
    rv$selected_isoform_names <- character(0)
    shinyjs::reset("analysisForm")
    output$isoforms_warning <- renderUI({ NULL })
    shinyjs::reset("seurat_file")
    shinyjs::reset("gtf")
    shinyjs::show("file_inputs_panel")
  })
  
  # Demo Button functionality
  observeEvent(input$DemoBtn, {
    app_state$demo_mode <- TRUE  # <--- Track that demo data was loaded
    
    shinyjs::hide("file_inputs_panel")
    shinyjs::hide("resetBtn")
    
    # Show spinner
    shinyjs::show("spinner")  # Assuming you have a div with id="spinner"
    # Hide spinner after loading data
    
    # Load demo Seurat + GTF, join layers once, then update UI
    seurat_obj(seurat_obj_demo)
    join_and_store(seurat_obj_demo)
    gtf(gtf_obj_demo)
    update_seurat_ui()
    
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
  seurat_obj    <- reactiveVal(NULL)
  seurat_joined <- reactiveVal(NULL)   # set once at load, never recomputed reactively
  gtf           <- reactiveVal(NULL)
  
  # Call this everywhere seurat_obj is set — JoinLayers once and store
  join_and_store <- function(obj) {
    if (is.null(obj)) { seurat_joined(NULL); return() }
    for (assay_name in names(obj@assays)) {
      tryCatch(
        { obj <- JoinLayers(obj, assay = assay_name) },
        error = function(e) {
          message("Could not join layers for assay '", assay_name, "': ", conditionMessage(e))
        }
      )
    }
    seurat_joined(obj)
  }
  
  # ── Helper: update all dropdowns after any Seurat object is loaded ──────────
  update_seurat_ui <- function() {
    req(seurat_obj())
    reductions    <- names(seurat_obj()@reductions)
    isoform_assay <- names(seurat_obj()@assays)
    group_by      <- colnames(seurat_obj()@meta.data)[!sapply(seurat_obj()@meta.data, is.numeric)]
    
    default_isoform_assay <- isoform_assay[grepl("^iso(form)?$", isoform_assay, ignore.case = TRUE)][1]
    default_group_by      <- group_by[grepl("seurat_clusters|cell.*type|annotation", group_by, ignore.case = TRUE)][1]
    default_reduction     <- reductions[grepl("umap|tsne|iso", reductions, ignore.case = TRUE)][1]
    
    updateSelectInput(session, "reduction",     choices = reductions,    selected = default_reduction     %||% reductions[1])
    updateSelectInput(session, "isoform_assay", choices = isoform_assay, selected = default_isoform_assay %||% isoform_assay[1])
    updateSelectInput(session, "group_by",      choices = group_by,      selected = default_group_by      %||% group_by[1])
    
    updateSelectizeInput(session, "feature", choices = rownames(seurat_obj()), server = TRUE)
    shinyjs::hide("spinner")
    
    if (!is.null(default_isoform_assay))
      showNotification(paste("Auto-selected assay:", default_isoform_assay), type = "message", duration = 4)
    if (!is.null(default_group_by))
      showNotification(paste("Auto-selected metadata column:", default_group_by), type = "message", duration = 4)
  }
  
  # Handle Seurat object file upload (local)
  observeEvent(input$seurat_file, {
    req(input$seurat_file)  # Ensure a file is uploaded
    
    shinyjs::show("spinner")
    obj <- tryCatch(
      load_seurat_object(input$seurat_file$datapath, input$seurat_file$name),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 8)
        NULL
      }
    )
    if (is.null(obj)) return()
    app_state$demo_mode <- FALSE
    rv$selected_isoform_names <- character(0)
    lastReset(input$GO)
    seurat_obj(obj)
    join_and_store(obj)
    update_seurat_ui()
  })
  
  # Handle Seurat object load from HPC path
  observeEvent(input$load_hpc_seurat, {
    path <- trimws(input$seurat_hpc_path)
    if (nchar(path) == 0) {
      showNotification("Please enter an HPC file path.", type = "warning", duration = 5)
      return()
    }
    if (!file.exists(path)) {
      showNotification(paste("File not found:", path), type = "error", duration = 8)
      return()
    }
    shinyjs::show("spinner")
    obj <- tryCatch(
      load_seurat_object(path, path),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 8)
        NULL
      }
    )
    if (is.null(obj)) { shinyjs::hide("spinner"); return() }
    app_state$demo_mode <- FALSE
    rv$selected_isoform_names <- character(0)
    lastReset(input$GO)
    seurat_obj(obj)
    join_and_store(obj)
    showNotification(paste("Loaded:", basename(path)), type = "message", duration = 5)
    update_seurat_ui()
  })
  
  # Handle GTF file upload (local)
  observeEvent(input$gtf, {
    req(input$gtf)
    shinyjs::show("spinner")
    gtf_obj <- tryCatch(
      rtracklayer::import(input$gtf$datapath) %>% as_tibble(),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 8)
        NULL
      }
    )
    if (is.null(gtf_obj)) { shinyjs::hide("spinner"); return() }
    gtf(gtf_obj)
    shinyjs::hide("spinner")
    app_state$demo_mode <- FALSE
  })
  
  # Handle GTF load from HPC path
  observeEvent(input$load_hpc_gtf, {
    path <- trimws(input$gtf_hpc_path)
    if (nchar(path) == 0) {
      showNotification("Please enter an HPC file path.", type = "warning", duration = 5)
      return()
    }
    if (!file.exists(path)) {
      showNotification(paste("File not found:", path), type = "error", duration = 8)
      return()
    }
    shinyjs::show("spinner")
    gtf_obj <- tryCatch(
      rtracklayer::import(path) %>% as_tibble(),
      error = function(e) {
        showNotification(conditionMessage(e), type = "error", duration = 8)
        NULL
      }
    )
    if (is.null(gtf_obj)) { shinyjs::hide("spinner"); return() }
    gtf(gtf_obj)
    shinyjs::hide("spinner")
    app_state$demo_mode <- FALSE
    showNotification(paste("Loaded GTF:", basename(path)), type = "message", duration = 5)
  })
  
  # Reactive expression triggered by "GO" button (updates everything)
  # Extract features related to isoforms
  # NEW: reactive version—auto‐invalidates when input$feature or input$GO is missing
  filtered_data <- eventReactive(input$GO, {
    req(seurat_obj(), input$feature)
    
    assay_name     <- input$isoform_assay
    isoform_feats  <- rownames(seurat_obj()@assays[[assay_name]]@features)
    expr_matrix    <- GetAssayData(seurat_joined(), assay = assay_name, layer = "data")
    
    matching_feats <- grep(
      paste0("(^|-|\\b)", input$feature, "($|\\b)"),
      isoform_feats, value = TRUE
    )
    subset_expr    <- expr_matrix[matching_feats, , drop = FALSE]
    total_expr     <- Matrix::rowSums(subset_expr)
    matching_feats <- names(sort(total_expr, decreasing = TRUE))
    
    list(
      celltype_plot     = DimPlot(seurat_joined(), reduction = input$reduction, group.by = input$group_by),
      feature_plot_gene = FeaturePlot(seurat_joined(), features = input$feature, reduction = input$reduction),
      vln_plot          = VlnPlot(seurat_joined(), features = input$feature, group.by = input$group_by),
      isoform_features  = matching_feats,
      expr_matrix       = expr_matrix   # <-- cached so isoform_table doesn't re-fetch
    )
  }, ignoreNULL = FALSE)  # <---- This makes it “fire once at launch” even though GO=0
  
  # Reactive function to get the top N isoform features
  #isoform_features_to_plot <- eventReactive(input$GO, {
  # req(filtered_data(), input$number_of_isoforms)
  
  #all_isoforms <- filtered_data()$isoform_features
  #gene_name    <- input$feature
  #num_isoforms <- length(all_isoforms)
  
  #if (input$number_of_isoforms > num_isoforms) {
  #  msg <- paste(gene_name, "has", num_isoforms, "isoforms.")
  #  output$isoforms_warning <- renderUI({
  #    HTML(
  #      paste0(
  #        "<p style='color:red; font-size:20px; text-align:center;'>",
  #        msg,
  #        "</p>"
  #      )
  #    )
  #  })
  #} else {
  #  output$isoforms_warning <- renderUI({ NULL })
  #}
  
  #ead(all_isoforms, input$number_of_isoforms)
  #})
  
  # Derived directly from filtered_data — no need for a second eventReactive
  isoform_features_to_plot <- reactive({
    req(filtered_data())
    filtered_data()$isoform_features
  })
  
  
  ## select the isofroms required to plot 
  #selected_isoforms <- reactive({
  #  req(isoform_features_to_plot())
  #  all_feats <- isoform_features_to_plot()
  #  sel_idx   <- rv$selected_rows
  
  # if (is.null(sel_idx) || length(sel_idx) == 0) {
  #  character(0)
  #} else {
  # Make sure we only index up to length(all_feats):
  # valid_idx <- sel_idx[sel_idx <= length(all_feats)]
  #  all_feats[valid_idx]
  #}
  #})
  
  isoform_plot <- reactive({
    req(input$GO > lastReset(), selected_isoforms(), seurat_joined())
    plots <- FeaturePlot(
      seurat_joined(),
      features  = selected_isoforms(),
      reduction = input$reduction,
      order     = TRUE
    )
    lapply(plots, function(pl) pl + theme(plot.title = element_text(size = 11))) %>%
      wrap_plots(ncol = 4)
  })
  
  dotplot_isoform <- reactive({
    req(input$GO > lastReset(), selected_isoforms(), seurat_joined())
    DotPlot(
      seurat_joined(),
      features = selected_isoforms(),
      assay    = input$isoform_assay,
      group.by = input$group_by
    ) + theme(axis.text.x = element_text(angle = 80, hjust = 1))
  })
  
  # Reactive function for the isoform Feature Plot
  #isoform_plot <- eventReactive(input$GO, {
  #req(selected_isoforms())
  #plots <- FeaturePlot(
  #seurat_obj(),
  #features  = selected_isoforms(),
  #reduction = input$reduction,
  # order     = TRUE
  #)
  #plots <- lapply(plots, function(pl) pl + theme(plot.title = element_text(size = 11)))
  # wrap_plots(plots = plots, ncol = 4)
  #})
  
  #dotplot_isoform <-  eventReactive(input$GO, {
  #   req(selected_isoforms())
  #    DotPlot(
  #seurat_obj(),
  # features = selected_isoforms(),
  #  assay    = input$isoform_assay,
  #   group.by = input$group_by
  #  ) + theme(axis.text.x = element_text(angle = 80, hjust = 1))
  #})
  
  
  # Reactive expression to trigger the heatmap plot when the button is clicked
  reactive_heatmap <- reactive({
    req(input$GO > lastReset(), selected_isoforms(), seurat_joined())
    plot_pseudobulk_heatmap(
      seurat_obj       = seurat_joined(),
      group.by         = input$group_by,
      isoforms_to_plot = selected_isoforms(),
      isoform_assay    = input$isoform_assay
    )
  })
  
  
  # ============================
  # Render the DataTable itself
  # ── Isoform stats table (read-only — selection via chips in sidebar) ─────
  output$isoform_table <- DT::renderDataTable({
    req(filtered_data(), length(selected_isoforms()) > 0)
    feats    <- selected_isoforms()   # only chip-selected isoforms
    expr_mat <- filtered_data()$expr_matrix
    expr_sub <- expr_mat[feats, , drop = FALSE]
    n_cells  <- ncol(expr_mat)        # denominator = all cells
    
    total_expr  <- as.numeric(Matrix::rowSums(expr_sub))
    num_cells   <- as.integer(Matrix::rowSums(expr_sub > 0))
    mean_expr   <- round(total_expr / n_cells, 4)
    max_expr    <- round(as.numeric(Matrix::rowMaxs(expr_sub)), 3)
    pct_cells   <- round(num_cells / n_cells * 100, 1)
    short_label <- sub("-[^-]+$", "", feats)
    
    df_iso <- data.frame(
      `Transcript ID`      = short_label,
      `Mean Norm. Expr.`   = mean_expr,
      `Max Norm. Expr.`    = max_expr,
      `% Cells Expressing` = pct_cells,
      `Cells (n)`          = num_cells,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    
    DT::datatable(
      df_iso,
      rownames  = FALSE,
      selection = "none",
      options   = list(
        dom            = 't',
        scrollY        = '280px',
        scrollCollapse = TRUE,
        paging         = FALSE,
        scrollX        = TRUE,
        order          = list(list(1, 'desc'))
      )
    ) %>%
      DT::formatRound(columns = c("Mean Norm. Expr.", "Max Norm. Expr."), digits = 3) %>%
      DT::formatStyle(
        "% Cells Expressing",
        background         = DT::styleColorBar(c(0, 100), '#d4eaf7'),
        backgroundSize     = '98% 70%',
        backgroundRepeat   = 'no-repeat',
        backgroundPosition = 'center'
      )
  }, server = FALSE)
  
  
  selected_isoforms <- reactive({
    req(isoform_features_to_plot())
    ranked <- isoform_features_to_plot()
    ranked[ranked %in% rv$selected_isoform_names]
  })
  
  # Shared colour map: rank 1 = hue 1 — used by chips, transcript plot, pie
  isoform_colour_map <- reactive({
    feats <- isoform_features_to_plot()
    req(length(feats) > 0)
    setNames(scales::hue_pal()(length(feats)), feats)
  })
  
  
  #### Render the main plots (with the new req) ####
  
  # 1) Gene Feature Plot (from filtered_data()$feature_plot_gene)
  output$feature_plot_gene <- renderPlot({
    req(input$GO > lastReset(), filtered_data())
    filtered_data()$feature_plot_gene   # <--- note “feature_plot_gene” matches the list name
  })
  
  # 2) Violin Plot
  output$vln_plot <- renderPlot({
    req(input$GO > lastReset(), filtered_data())
    filtered_data()$vln_plot
  })
  
  # 3) Cell Type DimPlot
  output$celltype_plot <- renderPlot({
    req(input$GO > lastReset(), filtered_data())
    filtered_data()$celltype_plot
  })
  
  # 4) Isoform Feature Plot
  output$feature_plot_iso <- renderPlot({
    req(isoform_plot())
    isoform_plot()
  })
  
  output$dot_plot_iso <- renderPlot({
    req(dotplot_isoform())
    dotplot_isoform()
  })
  
  # Transcript Structure (from GTF)
  output$transcript_plot <- renderPlot({
    req(input$GO > lastReset(), gtf(), selected_isoforms(), seurat_joined())
    if (is.null(gtf()) || nrow(gtf()) == 0) {
      showModal(modalDialog(
        title = "Missing GTF File",
        "Please provide a valid GTF to build the isoform stack.",
        easyClose = TRUE,
        footer    = NULL
      ))
    } else {
      sel     <- selected_isoforms()
      col_map <- isoform_colour_map()
      # GTF uses transcript IDs with version suffix but without -GENENAME
      # e.g. feature name "ENST00000380394.9-RPS6" -> GTF id "ENST00000380394.9"
      gtf_keys <- sub("-[^-]+$", "", sel)       # strip -GENENAME, keep .version
      fill_map <- stats::setNames(col_map[sel], gtf_keys)
      p <- plot_gene_transcripts(
        isoforms_to_plot = sel,
        gtf              = gtf()
      )
      p + ggplot2::scale_fill_manual(values = fill_map, na.value = "grey80")
    }
  })
  
  # 6) Pseudobulk Heatmap (Plotly)
  output$heatmap_plot <- renderPlotly({
    req(reactive_heatmap())
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
  
  # ── Isoform Proportions Pie Chart ──────────────────────────────────────────
  pie_plot <- reactive({
    req(input$GO > lastReset(), seurat_joined(), input$feature,
        input$isoform_assay, input$group_by, selected_isoforms())
    # Use the same colour map as chips — key by stripped transcript ID
    col_map        <- isoform_colour_map()
    stripped_names <- sub("-[^-]+$", "", names(col_map))  # strip -GENENAME
    stripped_names <- sub("[.][^.]*$", "", stripped_names) # strip version
    names(col_map) <- stripped_names
    plotIsoformPieFromSeurat(
      seurat_obj        = seurat_joined(),
      gene              = input$feature,
      selected_isoforms = selected_isoforms(),
      ncol              = input$pie_ncol %||% 4L,
      assay             = input$isoform_assay,
      layer             = "counts",
      cell_type_col     = input$group_by,
      min_counts        = input$pie_min_counts,
      colour_map        = col_map
    )
  })
  
  output$pie_plot <- renderPlot({
    req(pie_plot())
    pie_plot()
  })
  
  downloadModalServer(
    "isoform_pie",
    reactive({
      req(input$GO > lastReset())
      pie_plot()
    }),
    "Isoform Proportions Pie"
  )
  
  # ── Expression Trajectory ────────────────────────────────────────────────────
  
  # Step 1: populate the include selectize when the column changes
  # Populate trajectory include list from the sidebar group_by column
  observeEvent(input$group_by, {
    req(seurat_obj(), input$group_by)
    lvls <- sort(unique(as.character(seurat_obj()@meta.data[[input$group_by]])))
    updateSelectizeInput(session, "traj_include",
                         choices  = lvls,
                         selected = lvls)
  })
  
  # Step 2: render the rank_list from whatever is currently selected in traj_include
  output$traj_order_ui <- renderUI({
    req(input$traj_include)
    selected <- as.character(input$traj_include)
    sortable::rank_list(
      text     = NULL,
      labels   = selected,
      input_id = "traj_cell_order",
      options  = sortable::sortable_options(animation = 150)
    )
  })
  
  trajectory_plot <- reactive({
    req(input$GO > lastReset(), seurat_joined(), selected_isoforms(), input$group_by)
    ordered_lvls <- as.character(unlist(input$traj_cell_order))
    validate(
      need(length(ordered_lvls) >= 2,
           "Please select at least 2 cell types and press GO.")
    )
    plotExpressionTrajectory(
      seurat_obj     = seurat_joined(),
      isoforms       = selected_isoforms(),
      cell_type_col  = input$group_by,
      ordered_levels = ordered_lvls,
      assay          = input$isoform_assay,
      layer          = "data"
    )
  })
  
  output$trajectory_plot <- renderPlot({
    req(trajectory_plot())
    trajectory_plot()
  })
  
  # ── Isoform chip panel ────────────────────────────────────────────────────
  output$isoform_chips_panel <- renderUI({
    feats <- isoform_features_to_plot()
    if (is.null(feats) || length(feats) == 0) return(NULL)
    
    # Expression totals for the bar widths (reuse cached matrix)
    expr_mat   <- filtered_data()$expr_matrix
    expr_sub   <- expr_mat[feats, , drop = FALSE]
    totals     <- as.numeric(Matrix::rowSums(expr_sub))
    max_total  <- max(totals, 1)
    pct_widths <- round(totals / max_total * 100)
    
    # Colour palette — same hues as the pie chart, grey for unselected
    n_feats  <- length(feats)
    hues     <- scales::hue_pal()(n_feats)
    selected <- rv$selected_isoform_names
    
    # Short label: strip trailing -GENENAME
    short_labels <- sub("-[^-]+$", "", feats)
    
    # Build one chip div per isoform
    chips <- lapply(seq_along(feats), function(i) {
      feat      <- feats[i]
      is_sel    <- feat %in% selected
      bg_colour <- if (is_sel) hues[i] else "#f0f0f0"
      chip_cls  <- if (is_sel) "iso-chip chip-selected" else "iso-chip"
      
      tags$div(
        class   = chip_cls,
        style   = paste0("background:", bg_colour, "; border-color:", bg_colour, ";"),
        onclick = sprintf(
          "Shiny.setInputValue('chip_clicked', '%s', {priority: 'event'});",
          feat
        ),
        tags$span(class = "chip-label",
                  title = feat,          # full ID on hover
                  short_labels[i]),
        tags$div(class = "chip-bar-wrap",
                 tags$div(class = "chip-bar",
                          style = paste0("width:", pct_widths[i], "%;")))
      )
    })
    
    n_sel <- length(selected)
    tagList(
      tags$div(class = "chip-panel-title",
               paste0("Isoforms (", n_sel, " of ", n_feats, " selected)")),
      tags$div(class = "chip-panel-hint", "Click to toggle • Bar = relative expression"),
      tags$div(style = "display: flex; flex-wrap: wrap;", chips)
    )
  })
  
  downloadModalServer(
    "trajectory_plot",
    reactive({
      req(input$GO > lastReset())
      trajectory_plot()
    }),
    "Expression Trajectory"
  )
  
}