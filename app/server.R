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
    feats <- filtered_data()$isoform_features
    if (length(feats) == 0) return()
    # Preserve the user's current selection if it intersects with the new
    # feature set (same gene re-run). Only reset to top-4 if the gene changed
    # or no previously-selected isoforms exist in the new feature set.
    prev_sel <- rv$selected_isoform_names
    keep     <- prev_sel[prev_sel %in% feats]
    if (length(keep) > 0) {
      rv$selected_isoform_names <- keep
    } else {
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
  shinyjs::hide("analysisForm")
  
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
    
    # Demo: wipe everything. Own-data: keep analysis, just go back.
    if (app_state$demo_mode) {
      seurat_obj(NULL)
      join_and_store(NULL)
      gtf(NULL)
      app_state$demo_mode <- FALSE
      rv$selected_isoform_names <- character(0)
      lastReset(input$GO)
      shinyjs::reset("analysisForm")
      shinyjs::hide("analysisForm")
      shinyjs::show("file_inputs_panel")
    }
  })
  
  
  
  observeEvent(input$resetBtn, {
    # Full clean slate: clear data, plots, chips, form
    seurat_obj(NULL)
    join_and_store(NULL)
    gtf(NULL)
    rv$selected_isoform_names <- character(0)
    lastReset(input$GO)
    shinyjs::reset("analysisForm")
    shinyjs::hide("analysisForm")
    output$isoforms_warning <- renderUI({ NULL })
    shinyjs::reset("seurat_file")
    shinyjs::reset("gtf")
    shinyjs::show("file_inputs_panel")
    # traj_include_ui clears automatically via lastReset guard in renderUI
  })
  
  # Demo Button functionality
  observeEvent(input$DemoBtn, {
    app_state$demo_mode <- TRUE
    # Clear chips and plot state without resetting the form
    # (reset() would undo the dropdown values set by update_seurat_ui)
    rv$selected_isoform_names <- character(0)
    lastReset(input$GO)
    
    shinyjs::hide("file_inputs_panel")
    shinyjs::hide("resetBtn")
    shinyjs::show("spinner")
    
    seurat_obj(seurat_obj_demo)
    join_and_store(seurat_obj_demo)
    gtf(gtf_obj_demo)
    update_seurat_ui()
    # server=TRUE registers choices lazily — selected= is ignored at call time.
    # shinyjs::delay waits 300ms for the selectize widget to initialise client-side,
    # then sends the selection. This is more reliable than onFlushed for server-mode.
    # Delay 300ms for selectize to init, set VIM, then after another 400ms
    # auto-click GO so the top 4 isoforms are pre-selected for the user.
    shinyjs::delay(500, {
      updateSelectizeInput(session, "feature", selected = "VIM")
      shinyjs::delay(600, {
        shinyjs::click("GO")
        # Extra delay: if chips are still empty after GO fires
        # (selectize may not have delivered input$feature in time),
        # force-populate from whatever filtered_data computed.
        shinyjs::delay(800, {
          feats <- isolate(filtered_data()$isoform_features)
          if (length(feats) > 0 && length(isolate(rv$selected_isoform_names)) == 0) {
            rv$selected_isoform_names <- feats[seq_len(min(4, length(feats)))]
          }
        })
      })
    })
    
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
  
  # ── Sidebar collapse via shinyjs ─────────────────────────────────────────────
  observeEvent(input$sidebarToggle_click, {
    shinyjs::toggle(id = "sidebar_panel_wrap", anim = TRUE, animType = "slide", time = 0.3)
  }, ignoreInit = TRUE)
  
  # ── Helper: update all dropdowns after any Seurat object is loaded ──────────
  update_seurat_ui <- function() {
    req(seurat_obj())
    reductions    <- names(seurat_obj()@reductions)
    isoform_assay <- names(seurat_obj()@assays)
    group_by      <- colnames(seurat_obj()@meta.data)[!vapply(seurat_obj()@meta.data, function(x) is.numeric(x) || is.list(x), logical(1))]
    
    default_isoform_assay <- isoform_assay[grepl("^iso(form)?$", isoform_assay, ignore.case = TRUE)][1]
    default_group_by      <- group_by[grepl("seurat_clusters|cell.*type|annotation", group_by, ignore.case = TRUE)][1]
    default_reduction     <- reductions[grepl("umap|tsne|iso", reductions, ignore.case = TRUE)][1]
    
    updateSelectInput(session, "reduction",     choices = reductions,    selected = default_reduction     %||% reductions[1])
    updateSelectInput(session, "isoform_assay", choices = isoform_assay, selected = default_isoform_assay %||% isoform_assay[1])
    updateSelectInput(session, "group_by",      choices = group_by,      selected = default_group_by      %||% group_by[1])
    
    updateSelectizeInput(session, "feature", choices = rownames(seurat_obj()), server = TRUE)
    shinyjs::show("analysisForm")
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
    
    # Ensure gene-level assay is default so FeaturePlot/VlnPlot find the gene immediately
    seu_g <- seurat_joined()
    gene_assay_name <- setdiff(names(seu_g@assays), assay_name)[1]
    if (!is.na(gene_assay_name)) Seurat::DefaultAssay(seu_g) <- gene_assay_name

    list(
      celltype_plot     = DimPlot(seu_g, reduction = input$reduction, group.by = input$group_by),
      feature_plot_gene = FeaturePlot(seu_g, features = input$feature, reduction = input$reduction),
      vln_plot          = VlnPlot(seu_g, features = input$feature, group.by = input$group_by),
      isoform_features  = matching_feats,
      expr_matrix       = expr_matrix,
      feature           = input$feature
    )
  }, ignoreNULL = FALSE)
  
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
  
  # Snapshot only the settings at GO press (not isoforms).
  # isoform_plot and dotplot_isoform read selected_isoforms() live
  # so chip toggles update them, but changing reduction/assay/group_by
  # does not re-render until GO is pressed again.
  iso_settings <- eventReactive(input$GO, {
    list(
      reduction = input$reduction,
      assay     = input$isoform_assay,
      group_by  = input$group_by
    )
  }, ignoreNULL = FALSE)
  
  isoform_plot <- reactive({
    req(input$GO > lastReset(), iso_settings(), seurat_joined())
    req(length(selected_isoforms()) > 0)
    plots <- FeaturePlot(
      seurat_joined(),
      features  = selected_isoforms(),
      reduction = iso_settings()$reduction,
      order     = TRUE
    )
    lapply(plots, function(pl) pl + theme(plot.title = element_text(size = 11))) %>%
      wrap_plots(ncol = 4)
  })
  
  dotplot_isoform <- reactive({
    req(input$GO > lastReset(), iso_settings(), seurat_joined())
    req(length(selected_isoforms()) > 0)
    DotPlot(
      seurat_joined(),
      features = selected_isoforms(),
      assay    = iso_settings()$assay,
      group.by = iso_settings()$group_by
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
    req(input$GO > lastReset(), seurat_joined(), iso_settings())
    req(length(selected_isoforms()) > 0)
    plot_pseudobulk_heatmap(
      seurat_obj       = seurat_joined(),
      group.by         = iso_settings()$group_by,
      isoforms_to_plot = selected_isoforms(),
      isoform_assay    = iso_settings()$assay
    )
  })
  
  
  # ============================
  # Render the DataTable itself
  # ── Isoform stats table (read-only — selection via chips in sidebar) ─────
  output$isoform_table <- DT::renderDataTable({
    req(input$GO > lastReset(), filtered_data(), isoform_features_to_plot())
    feats    <- isoform_features_to_plot()   # all isoforms for the gene
    expr_mat <- filtered_data()$expr_matrix
    expr_sub <- expr_mat[feats, , drop = FALSE]
    n_cells  <- ncol(expr_mat)
    
    total_expr  <- as.numeric(Matrix::rowSums(expr_sub))
    num_cells   <- as.integer(Matrix::rowSums(expr_sub > 0))
    pct_cells   <- round(num_cells / n_cells * 100, 1)
    short_label <- sub("-[^-]+$", "", feats)
    is_selected <- feats %in% rv$selected_isoform_names
    # Median expression only among cells that actually express the isoform
    expr_dense  <- as.matrix(expr_sub)
    median_expr <- round(apply(expr_dense, 1, function(x) {
      vals <- x[x > 0]
      if (length(vals) == 0) 0 else median(vals)
    }), 2)
    
    df_iso <- data.frame(
      `Selected`                = ifelse(is_selected, "✔", ""),
      `Transcript ID`           = short_label,
      `Total Norm. Counts`      = round(total_expr, 1),
      `% Cells Expressing`      = pct_cells,
      `Cells (n)`               = num_cells,
      `Median Expr. (non-zero)` = median_expr,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    
    DT::datatable(
      df_iso,
      rownames   = FALSE,
      selection  = "none",
      extensions = "Buttons",
      options    = list(
        dom            = 'Bt',
        buttons        = list(list(extend = "csv", text = "Download CSV",
                                   filename = "isoform_statistics")),
        scrollY        = '380px',
        scrollCollapse = TRUE,
        paging         = FALSE,
        scrollX        = TRUE,
        order          = list(list(2, 'desc'))
      )
    ) %>%
      DT::formatRound(columns = c("Total Norm. Counts", "Median Expr. (non-zero)"), digits = 1) %>%
      DT::formatStyle(
        "% Cells Expressing",
        background         = DT::styleColorBar(c(0, 100), '#d4eaf7'),
        backgroundSize     = '98% 70%',
        backgroundRepeat   = 'no-repeat',
        backgroundPosition = 'center'
      ) %>%
      DT::formatStyle(
        "Selected",
        color      = "#2e7d32",
        fontWeight = "bold",
        fontSize   = "16px",
        textAlign  = "center"
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
  
  output$feature_plot_gene <- renderPlot({
    req(input$GO > lastReset(), filtered_data())
    filtered_data()$feature_plot_gene
  }) |> bindCache(filtered_data()$feature, iso_settings()$reduction, lastReset())
  output$vln_plot <- renderPlot({
    req(input$GO > lastReset(), filtered_data())
    filtered_data()$vln_plot
  }) |> bindCache(filtered_data()$feature, iso_settings()$group_by, lastReset())
  output$celltype_plot <- renderPlot({
    req(input$GO > lastReset(), filtered_data())
    filtered_data()$celltype_plot
  }) |> bindCache(iso_settings()$reduction, iso_settings()$group_by, filtered_data()$feature, lastReset())
  
  # 4) Isoform Feature Plot
  output$feature_plot_iso <- renderPlot({
    req(isoform_plot())
    isoform_plot()
  }) |> bindCache(selected_isoforms(), iso_settings()$reduction, lastReset())
  
  output$dot_plot_iso <- renderPlot({
    req(dotplot_isoform())
    dotplot_isoform()
  }) |> bindCache(selected_isoforms(), iso_settings()$assay, iso_settings()$group_by, lastReset())
  
  # Transcript Structure (from GTF)
  output$transcript_plot <- renderPlot({
    req(input$GO > lastReset(), selected_isoforms(), seurat_joined())
    if (is.null(gtf()) || length(gtf()) == 0 || (is.data.frame(gtf()) && nrow(gtf()) == 0)) {
      # Return a friendly placeholder ggplot
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0.5, y = 0.6, size = 6, colour = "#888",
                          label = "No GTF file loaded") +
        ggplot2::annotate("text", x = 0.5, y = 0.45, size = 4, colour = "#aaa",
                          label = "Upload a GTF file (or load from HPC) to view transcript structures.") +
        ggplot2::theme_void() +
        ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)
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
  }) |> bindCache(selected_isoforms(), iso_settings()$group_by, iso_settings()$assay, lastReset())
  
  
  
  # —— download handlers must be guarded ——  
  
  downloadModalServer(
    "celltype_plot",
    reactive({ req(input$GO > lastReset()); filtered_data()$celltype_plot }),
    "Cell type_plot"
  )
  downloadModalServer(
    "feature_plot_gene",
    reactive({ req(input$GO > lastReset()); filtered_data()$feature_plot_gene }),
    "Feature plot - Gene"
  )
  downloadModalServer(
    "vln_plot",
    reactive({ req(input$GO > lastReset()); filtered_data()$vln_plot }),
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
    req(input$GO > lastReset(), seurat_joined(), iso_settings(), filtered_data())
    req(length(selected_isoforms()) > 0)
    # Use the same colour map as chips — key by stripped transcript ID
    col_map        <- isoform_colour_map()
    stripped_names <- sub("-[^-]+$", "", names(col_map))  # strip -GENENAME
    stripped_names <- sub("[.][^.]*$", "", stripped_names) # strip version
    names(col_map) <- stripped_names
    plotIsoformPieFromSeurat(
      seurat_obj        = seurat_joined(),
      gene              = filtered_data()$feature,
      selected_isoforms = selected_isoforms(),
      ncol              = input$pie_ncol %||% 4L,
      assay             = iso_settings()$assay,
      layer             = "counts",
      cell_type_col     = iso_settings()$group_by,
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
  
  # Render the selectize as uiOutput so it is fully destroyed on reset
  # (updateSelectizeInput is unreliable for clearing server-side selectize).
  output$traj_include_ui <- renderUI({
    if (is.null(seurat_joined()) || input$GO <= lastReset()) return(NULL)
    req(iso_settings())
    col  <- iso_settings()$group_by
    lvls <- sort(unique(as.character(seurat_joined()@meta.data[[col]])))
    selectizeInput("traj_include", label = NULL,
                   choices  = lvls,
                   selected = lvls,
                   multiple = TRUE,
                   options  = list(placeholder = "Select cell types…"))
  })
  
  # Render drag list once GO has fired (iso_settings snapshotted).
  # Falls back to all valid levels if traj_include selectize hasn't round-tripped yet.
  output$traj_order_ui <- renderUI({
    req(input$GO > lastReset(), filtered_data(), seurat_joined(), iso_settings())
    col      <- iso_settings()$group_by
    all_lvls <- sort(unique(as.character(seurat_joined()@meta.data[[col]])))
    selected <- if (length(input$traj_include) > 0) as.character(input$traj_include) else all_lvls
    sortable::rank_list(
      text     = NULL,
      labels   = selected,
      input_id = "traj_cell_order",
      options  = sortable::sortable_options(animation = 150)
    )
  })
  
  # reactive() so chip toggles still update it; iso_settings() snapshots
  # group_by and assay at GO-press so live input changes have no effect.
  trajectory_plot <- reactive({
    req(input$GO > lastReset(), seurat_joined(), iso_settings())
    req(length(selected_isoforms()) > 0)
    # Filter to only levels that exist in the current column —
    # traj_cell_order may briefly hold stale labels from the previous
    # metadata column while traj_order_ui re-renders, causing combine_vars.
    # traj_cell_order is populated client-side after the rank list renders,
    # so on the first GO it may be empty. Fall back to all valid levels so
    # the plot loads immediately without requiring a second GO press.
    valid_col_lvls <- sort(unique(as.character(seurat_joined()@meta.data[[iso_settings()$group_by]])))
    ordered_lvls   <- as.character(unlist(input$traj_cell_order))
    ordered_lvls   <- ordered_lvls[ordered_lvls %in% valid_col_lvls]
    if (length(ordered_lvls) < 2) ordered_lvls <- valid_col_lvls
    validate(
      need(length(ordered_lvls) >= 2,
           "Please select at least 2 cell types and press GO.")
    )
    plotExpressionTrajectory(
      seurat_obj     = seurat_joined(),
      isoforms       = selected_isoforms(),
      cell_type_col  = iso_settings()$group_by,
      ordered_levels = ordered_lvls,
      assay          = iso_settings()$assay,
      layer          = "data"
    )
  })
  
  output$trajectory_plot <- renderPlot({
    req(trajectory_plot())
    trajectory_plot()
  })
  
  trajectory_gene_plot <- reactive({
    req(input$GO > lastReset(), seurat_joined(), iso_settings(), filtered_data())
    # traj_cell_order is populated client-side after the rank list renders,
    # so on the first GO it may be empty. Fall back to all valid levels so
    # the plot loads immediately without requiring a second GO press.
    valid_col_lvls <- sort(unique(as.character(seurat_joined()@meta.data[[iso_settings()$group_by]])))
    ordered_lvls   <- as.character(unlist(input$traj_cell_order))
    ordered_lvls   <- ordered_lvls[ordered_lvls %in% valid_col_lvls]
    if (length(ordered_lvls) < 2) ordered_lvls <- valid_col_lvls
    validate(need(length(ordered_lvls) >= 2, "Please select at least 2 cell types and set order below."))
    # Gene from filtered_data snapshot — not live input$feature
    gene <- filtered_data()$feature
    req(gene)
    default_assay <- DefaultAssay(seurat_joined())
    plotExpressionTrajectory(
      seurat_obj     = seurat_joined(),
      isoforms       = gene,
      cell_type_col  = iso_settings()$group_by,
      ordered_levels = ordered_lvls,
      assay          = default_assay,
      layer          = "data"
    )
  })
  
  output$trajectory_gene_plot <- renderPlot({
    req(trajectory_gene_plot())
    trajectory_gene_plot()
  })
  
  # ── Isoform chip panel ────────────────────────────────────────────────────
  output$isoform_chips_panel <- renderUI({
    # Return nothing if no data loaded or no GO has run since last reset
    if (is.null(seurat_obj()) || input$GO <= lastReset()) return(NULL)
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
      tags$div(style = "display: flex; align-items: center; gap: 6px; margin-bottom: 4px;",
               tags$span(class = "chip-panel-title",
                         paste0("Isoforms (", n_sel, " of ", n_feats, " selected)")),
               tags$span(
                 class = "chip-info-icon",
                 `data-toggle`   = "popover",
                 `data-trigger`  = "hover focus",
                 `data-placement`= "right",
                 `data-content`  = "Click to select or deselect an isoform for plotting. The bar shows each isoform's total expression as a fraction of the most highly expressed isoform for this gene.",
                 `data-html`     = "false",
                 tabindex = "0",
                 "i"
               )
      ),
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
  
  downloadModalServer(
    "trajectory_gene_plot",
    reactive({
      req(input$GO > lastReset())
      trajectory_gene_plot()
    }),
    "Gene Expression Trajectory"
  )
  
}