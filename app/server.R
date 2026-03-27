#server side

server <- function(input, output, session) {

  # --- Anonymous usage logging ---
  .session_start_time <- Sys.time()
  .session_log_id     <- log_session_start(APP_VERSION)
  session$onSessionEnded(function() {
    duration_secs <- as.numeric(difftime(Sys.time(), .session_start_time, units = "secs"))
    ram_mb <- tryCatch(
      as.numeric(trimws(system(paste("ps -o rss= -p", Sys.getpid()), intern = TRUE))) / 1024,
      error = function(e) NA_real_
    )
    log_session_end(.session_log_id, duration_secs, ram_mb)
  })


  #helper function for dyanmic seletcion of assay options 
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  #code for GPT helping with reset fucntionalltiy of back button 
  # 1) Keep track of “last reset” so that plots only
  #    render when input$GO > lastReset()
  # -----------------------------------------------
  
  lastReset  <- reactiveVal(0)
  gene_cache <- reactiveVal(list())   # in-session cache: key -> filtered_data result
  chip_help_dismissed <- reactiveVal(FALSE)

  observeEvent(input$dismiss_chip_help, { chip_help_dismissed(TRUE) })

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

  # Select all / Deselect all chip buttons
  observeEvent(input$chip_select_all, {
    feats <- isoform_features_to_plot()
    if (length(feats) > 0) rv$selected_isoform_names <- feats
  }, ignoreInit = TRUE)
  observeEvent(input$chip_deselect_all, {
    rv$selected_isoform_names <- character(0)
  }, ignoreInit = TRUE)

  # ── Recently analysed genes (auto-populated on each GO) ─────────────────
  recent_genes <- reactiveVal(character(0))

  # Add gene to recent list whenever a successful analysis completes
  observeEvent(filtered_data(), {
    gene <- tryCatch(filtered_data()$feature, error = function(e) NULL)
    if (!is.null(gene) && nchar(gene) > 0) {
      cur <- recent_genes()
      cur <- c(gene, setdiff(cur, gene))        # prepend (most recent first), deduplicate
      if (length(cur) > 12) cur <- cur[1:12]   # keep last 12
      recent_genes(cur)
    }
  }, ignoreInit = TRUE)

  # Click a recent gene → fill input + auto-run GO
  observeEvent(input$recent_gene_click, {
    gene <- input$recent_gene_click
    updateSelectizeInput(session, "feature", selected = gene)
    shinyjs::delay(350, { shinyjs::click("GO") })
  }, ignoreInit = TRUE)

  output$recent_genes_ui <- renderUI({
    genes <- recent_genes()
    if (length(genes) == 0) return(NULL)
    chips <- lapply(genes, function(g) {
      tags$span(
        style = paste0(
          "display:inline-block; background:#f0f4f8; color:#2c5282;",
          "border:1px solid #bee3f8; border-radius:10px; padding:1px 8px;",
          "font-size:11px; cursor:pointer; margin:2px; font-weight:500;",
          "transition:background 0.15s;"
        ),
        onclick = sprintf("Shiny.setInputValue('recent_gene_click','%s',{priority:'event'});", g),
        title  = paste0("Re-run analysis for ", g),
        g
      )
    })
    tags$div(
      style = "margin-bottom:5px;",
      tags$div(
        style = "font-size:10px; font-weight:700; color:#888; text-transform:uppercase; letter-spacing:0.05em; margin-bottom:3px;",
        HTML('<i class="fa fa-history" style="margin-right:3px;"></i>Recent')
      ),
      tags$div(style = "display:flex; flex-wrap:wrap;", chips)
    )
  })

  # Validate GO press — show friendly notification if data/gene missing
  observeEvent(input$GO, {
    if (is.null(seurat_obj())) {
      showNotification(
        HTML('<i class="fa fa-upload" style="margin-right:5px;"></i>Please upload a Seurat object first.'),
        type = "warning", duration = 5
      )
    } else if (is.null(input$feature) || nchar(trimws(input$feature)) == 0) {
      showNotification(
        HTML('<i class="fa fa-search" style="margin-right:5px;"></i>Please type and select a gene name first.'),
        type = "warning", duration = 5
      )
    }
  }, ignoreInit = TRUE)

  # Show report button only after a valid GO; hide on reset
  observeEvent(input$GO, {
    if (!is.null(seurat_obj()) && !is.null(input$feature) && nchar(trimws(input$feature)) > 0) {
      shinyjs::show("report_btn")
    }
  }, ignoreInit = TRUE)
  observeEvent(lastReset(), {
    shinyjs::hide("report_btn")
  }, ignoreInit = TRUE)

  # Report format modal
  observeEvent(input$report_btn, {
    showModal(modalDialog(
      title = tags$span(icon("file-alt"), " Download Report"),
      tags$p("Choose a format for your gene report:", style = "margin-bottom: 12px;"),
      shinyWidgets::prettyRadioButtons(
        inputId  = "report_format",
        label    = NULL,
        choices  = c("HTML (interactive)" = "html", "PDF" = "pdf"),
        selected = "html",
        inline   = FALSE,
        status   = "primary",
        shape    = "round"
      ),
      footer = tagList(
        modalButton("Cancel"),
        downloadButton("download_report", "Download",
                       icon  = icon("download"),
                       class = "btn btn-primary")
      ),
      easyClose = TRUE
    ))
  })

  #### Main logic ###
  # Hide everything except landing page initially
  shinyjs::hide("instructionsPage")
  shinyjs::hide("mainUI")
  shinyjs::hide("report_btn")   # hidden until GO is pressed
  
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
      recent_genes(character(0))
      lastReset(input$GO)
      shinyjs::reset("analysisForm")
      shinyjs::show("file_inputs_panel")
      updateSelectizeInput(session, "scatter_gene_x", choices = character(0), selected = NULL)
      updateSelectizeInput(session, "scatter_gene_y", choices = character(0), selected = NULL)
      updateSelectInput(session,    "scatter_iso_x",  choices = character(0), selected = NULL)
      updateSelectInput(session,    "scatter_iso_y",  choices = character(0), selected = NULL)
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
    output$isoforms_warning <- renderUI({ NULL })
    shinyjs::reset("seurat_file")
    shinyjs::reset("gtf")
    shinyjs::show("file_inputs_panel")
    updateSelectizeInput(session, "scatter_gene_x", choices = character(0), selected = NULL)
    updateSelectizeInput(session, "scatter_gene_y", choices = character(0), selected = NULL)
    updateSelectInput(session,    "scatter_iso_x",  choices = character(0), selected = NULL)
    updateSelectInput(session,    "scatter_iso_y",  choices = character(0), selected = NULL)
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
    
    join_and_store(seurat_obj_demo)
    gtf(gtf_obj_demo)
    update_seurat_ui()
    # server=TRUE registers choices lazily — selected= is ignored at call time.
    # shinyjs::delay waits 300ms for the selectize widget to initialise client-side,
    # then sends the selection. This is more reliable than onFlushed for server-mode.
    # Delay 300ms for selectize to init, set VIM, then after another 400ms
    # auto-click GO so the top 4 isoforms are pre-selected for the user.
    # server=TRUE selectize needs time to initialise before selected= takes effect.
    # 700ms for widget init, then 900ms for the selection to round-trip before GO fires.
    shinyjs::delay(700, {
      updateSelectizeInput(session, "feature", selected = "VIM")
      shinyjs::delay(900, {
        shinyjs::click("GO")
        # Fallback: if chips are still empty 1.2s after GO (e.g. slow round-trip),
        # populate from whatever filtered_data returned.
        shinyjs::delay(1200, {
          feats <- isolate(filtered_data()$isoform_features)
          if (length(feats) > 0 && length(isolate(rv$selected_isoform_names)) == 0) {
            rv$selected_isoform_names <- feats[seq_len(min(4, length(feats)))]
          }
          # Pre-populate the Compare Two Isoforms tab with VIM isoform 1 vs 2
          vim_feats <- iso_feats_for_gene("VIM")
          if (length(vim_feats) >= 2) {
            updateSelectizeInput(session, "scatter_gene_x", selected = "VIM")
            updateSelectizeInput(session, "scatter_gene_y", selected = "VIM")
            shinyjs::delay(400, {
              updateSelectInput(session, "scatter_iso_x", selected = vim_feats[1])
              updateSelectInput(session, "scatter_iso_y", selected = vim_feats[2])
            })
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
  
  # Clark Lab Data — password-gated loader
  observeEvent(input$clark_pw_attempt, {
    req(input$clark_pw_attempt)
    attempt <- input$clark_pw_attempt
    nm      <- attempt$dataset
    pw      <- attempt$pw

    # Wrong password — signal JS to show error, do nothing else
    if (!identical(pw, clark_password)) {
      session$sendCustomMessage("clark_pw_wrong", list())
      return()
    }

    # Correct — dismiss modal and load
    session$sendCustomMessage("clark_pw_ok", list())

    dataset <- clark_datasets[[nm]]
    if (is.null(dataset)) {
      showNotification(paste0("Dataset \u201c", nm, "\u201d not found."),
                       type = "error", duration = 8)
      return()
    }

    app_state$demo_mode <- TRUE
    rv$selected_isoform_names <- character(0)
    lastReset(input$GO)

    shinyjs::hide("file_inputs_panel")
    shinyjs::hide("resetBtn")
    shinyjs::show("spinner")

    notif_id <- showNotification(
      paste0("Loading \u201c", nm, "\u201d \u2014 this may take a moment\u2026"),
      type = "message", duration = NULL
    )

    tryCatch({
      obj <- load_seurat_object(dataset$seurat_path, basename(dataset$seurat_path))
      join_and_store(obj)

      if (!is.null(dataset$gtf_path)) {
        gtf(rtracklayer::import(dataset$gtf_path) %>% as_tibble())
      } else {
        gtf(NULL)
      }

      update_seurat_ui()

      removeNotification(notif_id)
      showNotification(paste0("\u201c", nm, "\u201d loaded successfully! Select a gene and press GO."),
                       type = "message", duration = 6)
    }, error = function(e) {
      removeNotification(notif_id)
      showNotification(paste0("Error loading \u201c", nm, "\u201d: ", conditionMessage(e)),
                       type = "error", duration = 10)
      shinyjs::hide("spinner")
      shinyjs::show("file_inputs_panel")
      return()
    })

    shinyjs::hide("landingPage")
    shinyjs::hide("spinner")
    shinyjs::show("mainUI")
  }, ignoreInit = TRUE, ignoreNULL = TRUE)

  #set the state of the app to control the demo mode resetting the input objects
  app_state <- reactiveValues(
    demo_mode = FALSE
  )
  
  #### Main logic ####
  seurat_obj    <- reactiveVal(NULL)
  seurat_joined <- reactiveVal(NULL)   # set once at load, never recomputed reactively
  gtf           <- reactiveVal(NULL)

  # Clear gene cache and recent genes whenever a new Seurat object is loaded
  observeEvent(seurat_joined(), {
    gene_cache(list())
    recent_genes(character(0))
  })
  
  # Call this everywhere seurat_obj is set — JoinLayers once and store.
  # Only join when there are genuinely split layers (e.g. "counts.sample1",
  # "counts.sample2"). A normal already-joined object ("counts", "data") is
  # passed through unchanged to avoid corrupting the assay structure.
  join_and_store <- function(obj) {
    if (is.null(obj)) { seurat_joined(NULL); return() }
    for (assay_name in names(obj@assays)) {
      tryCatch({
        lyrs  <- tryCatch(SeuratObject::Layers(obj@assays[[assay_name]]), error = function(e) character(0))
        # Split layers share a base name, e.g. "counts.s1" + "counts.s2" → base "counts" appears twice
        bases      <- sub("\\.[^.]+$", "", lyrs)
        has_splits <- length(lyrs) > 1 && any(duplicated(bases))
        if (has_splits) obj <- JoinLayers(obj, assay = assay_name)
      }, error = function(e) {
        message("Could not join layers for assay '", assay_name, "': ", conditionMessage(e))
      })
    }
    # Set RNA (or first gene-like assay) as default so rownames() returns gene features
    rna_assay <- names(obj@assays)[grepl("^(RNA|GEX|Spatial)$", names(obj@assays), ignore.case = TRUE)][1]
    if (!is.na(rna_assay)) Seurat::DefaultAssay(obj) <- rna_assay
    seurat_obj(obj)
    seurat_joined(obj)
  }
  
  # ── Sidebar collapse via shinyjs ─────────────────────────────────────────────
  observeEvent(input$sidebarToggle_click, {
    shinyjs::toggle(id = "sidebar_panel_wrap", anim = TRUE, animType = "slide", time = 0.3)
  }, ignoreInit = TRUE)
  
  # ── Helper: detect gene assay — prefers RNA/GEX by name, then any assay whose
  #    features don't look like isoform IDs (ENST…). Falls back to first non-iso assay.
  .detect_gene_assay <- function(obj, isoform_assay_name) {
    all_assays <- names(obj@assays)
    non_iso    <- setdiff(all_assays, isoform_assay_name)
    if (length(non_iso) == 0) return(all_assays[1])
    # Prefer an assay explicitly named RNA/GEX/Spatial
    by_name <- non_iso[grepl("^(RNA|GEX|Spatial)$", non_iso, ignore.case = TRUE)][1]
    if (!is.na(by_name)) return(by_name)
    # Otherwise pick the non-iso assay whose features least resemble ENST IDs
    for (a in non_iso) {
      feats <- tryCatch(rownames(obj@assays[[a]]), error = function(e) character(0))
      if (length(feats) > 0 && mean(grepl("^ENST", feats)) < 0.5) return(a)
    }
    non_iso[1]
  }

  # ── Helper: update all dropdowns after any Seurat object is loaded ──────────
  update_seurat_ui <- function() {
    req(seurat_obj())
    reductions    <- names(seurat_obj()@reductions)
    all_assays    <- names(seurat_obj()@assays)
    group_by      <- colnames(seurat_obj()@meta.data)[!vapply(seurat_obj()@meta.data, function(x) is.numeric(x) || is.list(x), logical(1))]

    default_isoform_assay <- all_assays[grepl("^iso(form)?$", all_assays, ignore.case = TRUE)][1]
    default_isoform_assay <- default_isoform_assay %||% all_assays[length(all_assays)]
    default_gene_assay    <- .detect_gene_assay(seurat_obj(), default_isoform_assay)

    # Prefer cell type over cluster for default grouping
    cell_type_col <- group_by[grepl("cell.*type|celltype", group_by, ignore.case = TRUE)][1]
    cluster_col   <- group_by[grepl("seurat_clusters|annotation", group_by, ignore.case = TRUE)][1]
    default_group_by  <- cell_type_col %||% cluster_col %||% group_by[1]
    default_reduction <- reductions[grepl("umap|tsne|iso", reductions, ignore.case = TRUE)][1]

    updateSelectInput(session, "reduction",     choices = reductions,                        selected = default_reduction  %||% reductions[1])
    updateSelectInput(session, "isoform_assay", choices = setdiff(all_assays, default_gene_assay),   selected = default_isoform_assay)
    updateSelectInput(session, "gene_assay",    choices = setdiff(all_assays, default_isoform_assay), selected = default_gene_assay)
    updateSelectInput(session, "group_by",      choices = group_by,                          selected = default_group_by)

    gene_features <- rownames(seurat_obj()@assays[[default_gene_assay]])
    updateSelectizeInput(session, "feature",        choices = gene_features, server = TRUE)
    updateSelectizeInput(session, "scatter_gene_x", choices = gene_features, server = TRUE)
    updateSelectizeInput(session, "scatter_gene_y", choices = gene_features, server = TRUE)
    shinyjs::hide("spinner")

    showNotification(paste("Auto-selected isoform assay:", default_isoform_assay), type = "message", duration = 4)
    showNotification(paste("Auto-selected gene assay:", default_gene_assay),    type = "message", duration = 4)
    if (!is.null(default_group_by))
      showNotification(paste("Auto-selected metadata column:", default_group_by), type = "message", duration = 4)
  }

  # ── Cross-filter: when isoform assay changes, remove it from gene assay choices
  observeEvent(input$isoform_assay, {
    req(seurat_obj(), input$isoform_assay)
    all_assays <- names(seurat_obj()@assays)
    gene_choices <- setdiff(all_assays, input$isoform_assay)
    cur_gene <- if (input$gene_assay %in% gene_choices) input$gene_assay else gene_choices[1]
    updateSelectInput(session, "gene_assay", choices = gene_choices, selected = cur_gene)
  }, ignoreInit = TRUE)

  # ── Cross-filter: when gene assay changes, remove it from isoform assay choices
  observeEvent(input$gene_assay, {
    req(seurat_obj(), input$gene_assay)
    all_assays <- names(seurat_obj()@assays)
    iso_choices <- setdiff(all_assays, input$gene_assay)
    cur_iso <- if (input$isoform_assay %in% iso_choices) input$isoform_assay else iso_choices[1]
    updateSelectInput(session, "isoform_assay", choices = iso_choices, selected = cur_iso)
  }, ignoreInit = TRUE)
  
  # Handle Seurat object file upload (local)
  observeEvent(input$seurat_file, {
    req(input$seurat_file)

    shinyjs::show("spinner")
    withProgress(message = "Loading Seurat object\u2026", value = 0, {
      setProgress(0.1, detail = paste0("Reading \u201c", input$seurat_file$name, "\u201d\u2026"))
      obj <- tryCatch(
        load_seurat_object(input$seurat_file$datapath, input$seurat_file$name),
        error = function(e) {
          showNotification(conditionMessage(e), type = "error", duration = 8)
          NULL
        }
      )
      if (is.null(obj)) return()
      setProgress(0.65, detail = "Joining assay layers\u2026")
      app_state$demo_mode <- FALSE
      rv$selected_isoform_names <- character(0)
      lastReset(input$GO)
      join_and_store(obj)
      setProgress(0.9, detail = "Updating interface\u2026")
      update_seurat_ui()
    })
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
  
  # Reactive expression triggered by "GO" button (updates everything)
  # Results are cached in gene_cache() by a key combining gene + settings,
  # so switching back to a previously-viewed gene is near-instant.
  # Isoform stats (rowSums / pct expressing / median) are pre-computed here
  # so the table render only needs to add the chip-selection column.
  filtered_data <- eventReactive(input$GO, {
    req(seurat_obj(), input$feature)

    cache_key <- paste(input$feature, input$isoform_assay,
                       input$reduction, input$group_by, sep = "||")
    cache <- gene_cache()
    if (!is.null(cache[[cache_key]])) return(cache[[cache_key]])

    assay_name    <- input$isoform_assay
    isoform_feats <- tryCatch(
      rownames(seurat_obj()@assays[[assay_name]]@features),
      error = function(e) rownames(seurat_obj()@assays[[assay_name]])
    )
    # Prefer the normalised "data" layer; fall back to "counts" if absent
    expr_matrix <- tryCatch(
      GetAssayData(seurat_joined(), assay = assay_name, layer = "data"),
      error = function(e)
        GetAssayData(seurat_joined(), assay = assay_name, layer = "counts")
    )

    matching_feats <- grep(
      paste0("(^|-|\\b)", input$feature, "($|\\b)"),
      isoform_feats, value = TRUE
    )
    subset_expr    <- expr_matrix[matching_feats, , drop = FALSE]
    total_expr_ord <- Matrix::rowSums(subset_expr)
    matching_feats <- names(sort(total_expr_ord, decreasing = TRUE))

    # ── Pre-compute isoform stats (sparse-safe — no dense matrix conversion) ──
    expr_sub    <- expr_matrix[matching_feats, , drop = FALSE]
    n_cells     <- ncol(expr_matrix)
    total_expr  <- as.numeric(Matrix::rowSums(expr_sub))
    num_cells   <- as.integer(Matrix::rowSums(expr_sub > 0))
    pct_cells   <- round(num_cells / n_cells * 100, 1)
    short_label <- sub("-[^-]+$", "", matching_feats)
    median_expr <- round(vapply(seq_len(nrow(expr_sub)), function(i) {
      nz <- as.numeric(expr_sub[i, ])   # force sparse row vector to plain numeric
      nz <- nz[nz > 0]
      if (length(nz) == 0) 0 else median(nz)
    }, numeric(1)), 2)
    # Raw counts per isoform (from counts layer, robust v4/v5)
    counts_mat <- tryCatch(
      GetAssayData(seurat_joined(), assay = assay_name, layer = "counts"),
      error = function(e) tryCatch(
        GetAssayData(seurat_joined(), assay = assay_name, slot  = "counts"),
        error = function(e2) NULL
      )
    )
    raw_counts <- if (!is.null(counts_mat) && all(matching_feats %in% rownames(counts_mat))) {
      as.integer(Matrix::rowSums(counts_mat[matching_feats, , drop = FALSE]))
    } else {
      rep(NA_integer_, length(matching_feats))
    }
    isoform_stats_df <- data.frame(
      feat        = matching_feats,
      short_label = short_label,
      total_expr  = total_expr,
      raw_counts  = raw_counts,
      pct_cells   = pct_cells,
      num_cells   = num_cells,
      median_expr = median_expr,
      stringsAsFactors = FALSE
    )

    # Ensure gene-level assay is default so FeaturePlot/VlnPlot find the gene immediately
    seu_g <- seurat_joined()
    gene_assay_name <- if (!is.null(input$gene_assay) && nzchar(input$gene_assay)) {
      input$gene_assay
    } else {
      setdiff(names(seu_g@assays), assay_name)[1]
    }
    if (!is.null(gene_assay_name) && !is.na(gene_assay_name)) Seurat::DefaultAssay(seu_g) <- gene_assay_name

    # Friendly placeholder when gene has no isoforms in the selected assay
    if (length(matching_feats) == 0) {
      showNotification(
        paste0("\u26a0\ufe0f \"", input$feature, "\" was not found in the \"",
               assay_name, "\" assay. Check the gene name and Isoform Assay setting."),
        type = "warning", duration = 8
      )
    }

    .gene_plot_or_msg <- function(expr, gene) {
      tryCatch(expr, error = function(e) {
        ggplot2::ggplot() +
          ggplot2::annotate("text", x=0.5, y=0.6, label=paste0("\"", gene, "\" not found in gene assay"),
                            size=5, colour="#c0392b", hjust=0.5) +
          ggplot2::annotate("text", x=0.5, y=0.4, label="Check gene name or assay selection",
                            size=3.5, colour="#aaa", hjust=0.5) +
          ggplot2::theme_void() + ggplot2::xlim(0,1) + ggplot2::ylim(0,1)
      })
    }

    .plot_or_err <- function(expr, label = "Plot") {
      tryCatch(expr, error = function(e) {
        ggplot2::ggplot() +
          ggplot2::annotate("text", x=0.5, y=0.6, label=paste0(label, " could not be rendered"),
                            size=5, colour="#c0392b", hjust=0.5) +
          ggplot2::annotate("text", x=0.5, y=0.4, label=conditionMessage(e),
                            size=3, colour="#aaa", hjust=0.5) +
          ggplot2::theme_void() + ggplot2::xlim(0,1) + ggplot2::ylim(0,1)
      })
    }

    result <- list(
      celltype_plot     = .plot_or_err(
        DimPlot(seu_g, reduction = input$reduction, group.by = input$group_by),
        "Cell type UMAP"),
      feature_plot_gene = .gene_plot_or_msg(
        FeaturePlot(seu_g, features = input$feature, reduction = input$reduction), input$feature),
      vln_plot          = .gene_plot_or_msg(
        VlnPlot(seu_g, features = input$feature, group.by = input$group_by,
                pt.size = 0.3, alpha = 0.4), input$feature),
      isoform_features  = matching_feats,
      feature           = input$feature,
      isoform_stats_df  = isoform_stats_df
    )

    # Free memory before caching — do not hold the full expression matrix in cache
    gc()

    new_cache <- cache
    new_cache[[cache_key]] <- result
    # LRU cap: keep the most recent 10 entries
    if (length(new_cache) > 10) new_cache <- new_cache[tail(seq_along(new_cache), 10)]
    gene_cache(new_cache)
    result
  }, ignoreNULL = FALSE)
  
  # Derived directly from filtered_data — no need for a second eventReactive
  isoform_features_to_plot <- reactive({
    req(filtered_data())
    filtered_data()$isoform_features
  })
  

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
  
  # Returns a plain list of per-isoform ggplot objects (no shared-scale applied,
  # no combining).  Shared-scale logic is applied downstream so that toggling
  # input$shared_scale triggers a fresh render without rerunning FeaturePlot.
  isoform_plot <- reactive({
    req(input$GO > lastReset(), iso_settings(), seurat_joined())
    req(length(selected_isoforms_debounced()) > 0)
    # Require a valid assay — iso_settings() can fire before input$isoform_assay
    # is populated (ignoreNULL=FALSE + demo auto-GO timing), and caching a plot
    # built with a NULL assay poisons subsequent shared-scale lookups for VIM.
    assay_name <- iso_settings()$assay
    req(!is.null(assay_name) && nzchar(assay_name))
    seu <- seurat_joined()
    Seurat::DefaultAssay(seu) <- assay_name

    FeaturePlot(
      seu,
      features  = selected_isoforms_debounced(),
      reduction = iso_settings()$reduction,
      order     = TRUE,
      combine   = FALSE
    )
  })

  # Apply shared colour scale (if requested) and combine into a patchwork.
  # Called from both renderPlot and the download handler so both always reflect
  # the current toggle state.
  .combine_iso_feature_plot <- function(plots_list, shared_scale, seu, assay, sel_isoforms, use_viridis = FALSE) {
    pal <- if (isTRUE(use_viridis)) viridis::viridis(20) else c("lightgrey", "blue")
    if (isTRUE(shared_scale)) {
      # Fall back to the object's current default assay if assay is NULL/empty
      assay <- if (!is.null(assay) && nzchar(assay)) assay else Seurat::DefaultAssay(seu)
      expr_mat <- tryCatch(
        GetAssayData(seu, assay = assay, layer = "data"),
        error = function(e) tryCatch(
          GetAssayData(seu, assay = assay, slot = "data"),
          error = function(e2) NULL
        )
      )
      message("[shared_scale] assay=", assay,
              " | sel_isoforms=", length(sel_isoforms),
              " | expr_mat rows=", if (is.null(expr_mat)) "NULL" else nrow(expr_mat))
      if (!is.null(expr_mat)) {
        sel  <- sel_isoforms[sel_isoforms %in% rownames(expr_mat)]
        message("[shared_scale] matched=", length(sel))
        vals <- as.numeric(expr_mat[sel, , drop = FALSE])
        vals <- vals[is.finite(vals)]
        if (length(vals) > 0) {
          g_min <- min(vals); g_max <- max(vals)
          message("[shared_scale] g_min=", round(g_min,3), " g_max=", round(g_max,3))
          plots_list <- lapply(plots_list, function(p) {
            p + scale_color_gradientn(colours = pal, limits = c(g_min, g_max))
          })
        } else {
          message("[shared_scale] vals empty — scale NOT applied")
        }
      }
    } else if (isTRUE(use_viridis)) {
      plots_list <- lapply(plots_list, function(p) {
        p + scale_color_gradientn(colours = pal)
      })
    }
    lapply(plots_list, function(p) p + theme(
      plot.title   = element_text(size = 14, face = "bold", family = "sans"),
      axis.text    = element_text(size = 12, family = "sans"),
      legend.text  = element_text(size = 12, family = "sans"),
      legend.title = element_text(size = 13, face = "bold", family = "sans")
    )) %>%
      wrap_plots(ncol = 4)
  }

  dotplot_isoform <- reactive({
    req(input$GO > lastReset(), iso_settings(), seurat_joined())
    req(length(selected_isoforms_debounced()) > 0)

    # Strip gene symbol suffix — keep only the ENST ID
    feats        <- selected_isoforms_debounced()
    feats_labels <- setNames(sub("-.*", "", feats), feats)

    p <- DotPlot(
      seurat_joined(),
      features = feats,
      assay    = iso_settings()$assay,
      group.by = iso_settings()$group_by,
      col.min  = 0,
      cols     = c("white", "#08306b")
    )

    # Relabel x-axis ticks and drop "Features" axis title
    p + scale_x_discrete(labels = feats_labels) +
      theme(
        axis.text.x  = element_text(angle = 80, hjust = 1, size = 13, family = "sans"),
        axis.text.y  = element_text(size = 13, family = "sans"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14, face = "bold", family = "sans"),
        legend.text  = element_text(size = 13, family = "sans"),
        legend.title = element_text(size = 14, face = "bold", family = "sans"),
        legend.key.size  = unit(1.1, "lines"),
        legend.spacing.y = unit(0.4, "cm"),
        plot.margin  = margin(10, 20, 10, 10)
      )
  })
  
  # Reactive expression to trigger the heatmap plot when the button is clicked
  reactive_heatmap <- reactive({
    req(input$GO > lastReset(), seurat_joined(), iso_settings())
    req(length(selected_isoforms_debounced()) > 0)
    validate(need(length(selected_isoforms_debounced()) >= 2,
                  "Please select at least 2 isoforms to display the heatmap."))
    plot_pseudobulk_heatmap(
      seurat_obj       = seurat_joined(),
      group.by         = iso_settings()$group_by,
      isoforms_to_plot = selected_isoforms_debounced(),
      isoform_assay    = iso_settings()$assay
    )
  })
  
  
  # ============================
  # Render the DataTable itself
  # ── Isoform stats table (read-only — selection via chips in sidebar) ─────
  # Heavy computation (dense matrix, row medians) is pre-computed inside
  # filtered_data(); here we only add the chip-selection column, so chip
  # toggles re-render instantly without repeating the expensive work.
  output$isoform_table <- DT::renderDataTable({
    if (!isTRUE(input$GO > lastReset())) {
      return(DT::datatable(
        data.frame(Message = "Select a gene and press \u25b6 Run Analysis to populate this table."),
        rownames = FALSE, options = list(dom = 't'), selection = "none"
      ))
    }
    req(filtered_data(), isoform_features_to_plot())
    stats_df    <- filtered_data()$isoform_stats_df
    is_selected <- stats_df$feat %in% rv$selected_isoform_names

    df_iso <- data.frame(
      `Selected`                = ifelse(is_selected, "\u2714", ""),
      `Transcript ID`           = stats_df$short_label,
      `Total Norm. Counts`      = round(stats_df$total_expr, 1),
      `Raw Counts`              = ifelse(is.na(stats_df$raw_counts), NA_integer_, stats_df$raw_counts),
      `% Cells Expressing`      = stats_df$pct_cells,
      `Cells (n)`               = stats_df$num_cells,
      `Median Expr. (non-zero)` = stats_df$median_expr,
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
        autoWidth      = FALSE,
        order          = list(list(2, 'desc')),
        initComplete   = JS("function(settings,json){this.api().columns.adjust().draw();}")
      )
    ) %>%
      DT::formatRound(columns = c("Total Norm. Counts", "Median Expr. (non-zero)"), digits = 1) %>%
      DT::formatRound(columns = "Raw Counts", digits = 0) %>%
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

  # Debounced version for heavy plots — prevents re-render on every chip toggle
  selected_isoforms_debounced <- debounce(selected_isoforms, 250)

  # Shared colour map: rank 1 = hue 1 — used by chips, transcript plot, pie
  isoform_colour_map <- reactive({
    feats <- isoform_features_to_plot()
    req(length(feats) > 0)
    setNames(scales::hue_pal()(length(feats)), feats)
  })
  
  
  # ── Dynamic plot box titles (use actual group_by column name) ───────────────
  .grp_label <- function() {
    grp <- tryCatch(iso_settings()$group_by, error = function(e) NULL)
    if (!is.null(grp) && nchar(grp) > 0) grp else "Group"
  }
  output$title_celltype_plot    <- renderText({ .grp_label() })
  output$title_vln_plot         <- renderText({
    gene <- tryCatch(filtered_data()$feature, error = function(e) "Gene")
    paste0(gene, " by ", .grp_label())
  })
  output$title_trajectory_gene  <- renderText({
    paste0("Gene Expression Across ", .grp_label())
  })
  output$title_isoform_trajectory <- renderText({
    paste0("Isoform Expression Across ", .grp_label())
  })
  output$title_pie_plot <- renderText({
    paste0("Isoform Proportions by ", .grp_label())
  })
  # Use input$group_by directly (pre-GO) so trajectory section headers update immediately
  output$traj_grp_label  <- renderText({ if (nchar(input$group_by) > 0) input$group_by else "Group" })
  output$traj_grp_label2 <- renderText({ if (nchar(input$group_by) > 0) input$group_by else "Group" })

  #### Render the main plots (with the new req) ####
  
  .empty_plot <- function(msg = "Select a gene and press \u25b6 Run Analysis to view this plot.") {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5, label = msg,
                        size = 5, colour = "#999", hjust = 0.5, vjust = 0.5) +
      ggplot2::theme_void() + ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1)
  }

  .safe_feature <- function() tryCatch(filtered_data()$feature, error = function(e) NULL)
  .safe_reduction <- function() tryCatch(iso_settings()$reduction, error = function(e) NULL)
  .safe_group_by  <- function() tryCatch(iso_settings()$group_by,  error = function(e) NULL)

  output$feature_plot_gene <- renderPlot({
    if (!isTRUE(input$GO > lastReset())) return(.empty_plot())
    req(filtered_data())
    p <- filtered_data()$feature_plot_gene
    if (isTRUE(input$gene_viridis)) {
      p + scale_color_gradientn(colours = viridis::viridis(20))
    } else {
      p
    }
  })
  output$vln_plot <- renderPlot({
    if (!isTRUE(input$GO > lastReset())) return(.empty_plot())
    req(filtered_data())
    filtered_data()$vln_plot
  })
  output$celltype_plot <- renderPlot({
    if (!isTRUE(input$GO > lastReset())) return(.empty_plot())
    req(filtered_data())
    filtered_data()$celltype_plot
  })

  .no_isoforms_plot <- function() .empty_plot(
    "No isoforms selected \u2014 click a chip in the sidebar to select at least one."
  )

  # 4) Isoform Feature Plot — no bindCache here so that toggling input$shared_scale
  # always triggers a fresh render.  The expensive FeaturePlot computation is
  # cached inside the isoform_plot() reactive, so performance is unchanged.
  output$feature_plot_iso <- renderPlot({
    if (!isTRUE(input$GO > lastReset())) return(.empty_plot())
    if (length(selected_isoforms()) == 0) return(.no_isoforms_plot())
    .combine_iso_feature_plot(
      plots_list   = isoform_plot(),
      shared_scale = input$shared_scale,
      seu          = seurat_joined(),
      assay        = iso_settings()$assay,
      sel_isoforms = selected_isoforms_debounced(),
      use_viridis  = isTRUE(input$feature_viridis)
    )
  })

  output$dot_plot_iso <- renderPlot({
    if (!isTRUE(input$GO > lastReset())) return(.empty_plot())
    if (length(selected_isoforms()) == 0) return(.no_isoforms_plot())
    dotplot_isoform()
  }) |> bindCache(selected_isoforms_debounced(), iso_settings()$assay, iso_settings()$group_by, lastReset())
  
  # Transcript Structure (from GTF)
  output$transcript_plot <- renderPlot({
    if (!isTRUE(input$GO > lastReset())) return(.empty_plot())
    if (length(selected_isoforms()) == 0) return(.no_isoforms_plot())
    req(seurat_joined())
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
    if (!isTRUE(input$GO > lastReset())) return(NULL)
    if (length(selected_isoforms()) == 0) return(NULL)
    req(reactive_heatmap())
    reactive_heatmap()
  }) |> bindCache(selected_isoforms_debounced(), iso_settings()$group_by, iso_settings()$assay, lastReset())
  
  
  
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
      req(input$GO > lastReset())
      .combine_iso_feature_plot(
        plots_list   = isoform_plot(),
        shared_scale = input$shared_scale,
        seu          = seurat_joined(),
        assay        = iso_settings()$assay,
        sel_isoforms = selected_isoforms_debounced(),
        use_viridis  = isTRUE(input$feature_viridis)
      )
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
  
  # ── Pseudobulk heatmap data download ──────────────────────────────────────
  output$heatmap_download_csv <- downloadHandler(
    filename = function() {
      gene <- tryCatch(isolate(filtered_data()$feature), error = function(e) "heatmap")
      paste0("pseudobulk_heatmap_", gene, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(seurat_joined(), iso_settings(), selected_isoforms_debounced())
      agg <- AggregateExpression(seurat_joined(),
                                  group.by   = iso_settings()$group_by,
                                  assays     = iso_settings()$assay,
                                  slot       = "counts")[[iso_settings()$assay]]
      mat     <- as.matrix(agg)
      cpm_mat <- t(t(mat) / colSums(mat) * 1e6)
      log_cpm <- log1p(cpm_mat)
      sub_mat <- log_cpm[rownames(log_cpm) %in% selected_isoforms_debounced(), , drop = FALSE]
      df_out  <- as.data.frame(sub_mat)
      # Strip -GENENAME suffix, keep ENST ID + version number
      transcript_ids <- sub("-[^-]+$", "", rownames(df_out))
      df_out  <- cbind(transcript_id = transcript_ids, df_out)
      write.csv(df_out, file, row.names = FALSE)
    }
  )

  # ── Isoform Proportions Pie Chart ──────────────────────────────────────────
  pie_plot <- reactive({
    req(input$GO > lastReset(), seurat_joined(), iso_settings(), filtered_data())
    req(length(selected_isoforms_debounced()) > 0)
    # Use the same colour map as chips — key by stripped transcript ID
    col_map        <- isoform_colour_map()
    stripped_names <- sub("-[^-]+$", "", names(col_map))  # strip -GENENAME
    stripped_names <- sub("[.][^.]*$", "", stripped_names) # strip version
    names(col_map) <- stripped_names
    seu_pie <- seurat_joined()
    sel_pie <- input$pie_cells_filter
    grp_pie <- iso_settings()$group_by
    if (!is.null(sel_pie) && length(sel_pie) > 0 && !is.null(grp_pie) &&
        nzchar(grp_pie) && grp_pie %in% colnames(seu_pie@meta.data)) {
      keep_pie <- colnames(seu_pie)[as.character(seu_pie@meta.data[[grp_pie]]) %in% sel_pie]
      if (length(keep_pie) >= 1) seu_pie <- subset(seu_pie, cells = keep_pie)
    }
    plotIsoformPieFromSeurat(
      seurat_obj        = seu_pie,
      gene              = filtered_data()$feature,
      selected_isoforms = selected_isoforms_debounced(),
      ncol              = input$pie_ncol %||% 4L,
      assay             = iso_settings()$assay,
      layer             = "counts",
      cell_type_col     = iso_settings()$group_by,
      min_counts        = input$pie_min_counts,
      colour_map        = col_map,
      group_by_label    = iso_settings()$group_by,
      display_mode      = input$pie_display_mode %||% "proportion"
    )
  })

  output$pie_plot <- renderPlot({
    if (!isTRUE(input$GO > lastReset())) return(.empty_plot())
    if (length(selected_isoforms()) == 0) return(.no_isoforms_plot())
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
    req(length(selected_isoforms_debounced()) > 0)
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
      isoforms       = selected_isoforms_debounced(),
      cell_type_col  = iso_settings()$group_by,
      ordered_levels = ordered_lvls,
      assay          = iso_settings()$assay,
      layer          = "data",
      shared_y       = isTRUE(input$traj_shared_y)
    )
  })
  
  output$trajectory_plot <- renderPlot({
    if (!isTRUE(input$GO > lastReset())) return(.empty_plot())
    if (length(selected_isoforms()) == 0) return(.no_isoforms_plot())
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
    if (!isTRUE(input$GO > lastReset())) return(.empty_plot())
    req(trajectory_gene_plot())
    trajectory_gene_plot()
  })
  
  # ── Isoform chip panel ────────────────────────────────────────────────────
  output$isoform_chips_panel <- renderUI({
    # Return nothing if no data loaded or no GO has run since last reset
    if (is.null(seurat_obj()) || input$GO <= lastReset()) return(NULL)
    feats <- isoform_features_to_plot()
    if (is.null(feats) || length(feats) == 0) return(NULL)
    
    # Expression totals for bar widths — use pre-computed stats (no full matrix in cache)
    stats_df   <- filtered_data()$isoform_stats_df
    totals     <- stats_df$total_expr[match(feats, stats_df$feat)]
    totals     <- ifelse(is.na(totals), 0, totals)
    max_total  <- max(totals, 1)
    pct_widths <- round(totals / max_total * 100)
    
    # Colour palette — reuse cached isoform_colour_map() (same as pie/transcript)
    n_feats  <- length(feats)
    col_map  <- isoform_colour_map()
    selected <- rv$selected_isoform_names
    
    # Short label: strip trailing -GENENAME
    short_labels <- sub("-[^-]+$", "", feats)
    
    # Build one chip div per isoform
    chips <- lapply(seq_along(feats), function(i) {
      feat      <- feats[i]
      is_sel    <- feat %in% selected
      bg_colour <- if (is_sel) col_map[[feat]] else "#f0f0f0"
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

    help_card <- if (!chip_help_dismissed()) {
      tags$div(
        style = paste0(
          "background:#eaf6f2; border:1px solid #b2dfdb; border-radius:6px;",
          "padding:7px 9px; margin-top:8px; font-size:11.5px; color:#2c5444;",
          "display:inline-flex; align-items:flex-start; gap:7px;"
        ),
        tags$span(style = "font-size:14px; flex-shrink:0;", "\U0001f4a1"),
        tags$div(
          style = "flex:1; line-height:1.45;",
          tags$strong("Tip: "),
          "Click chips to select/deselect isoforms for plotting."
        ),
        tags$button(
          id = "dismiss_chip_help",
          style = paste0(
            "background:none; border:none; color:#2c5444; cursor:pointer;",
            "font-size:13px; padding:0; flex-shrink:0; opacity:0.6; line-height:1;"
          ),
          onclick = "Shiny.setInputValue('dismiss_chip_help', Math.random());",
          "\u00d7"
        )
      )
    } else NULL

    tagList(
      tags$div(style = "display: flex; align-items: center; gap: 6px; margin-bottom: 4px; flex-wrap: wrap;",
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
               ),
               tags$button(
                 style = paste0("background:none; border:1px solid #ccc; border-radius:8px;",
                                "font-size:10px; color:#555; padding:1px 6px; cursor:pointer;",
                                "margin-left:auto;"),
                 onclick = "Shiny.setInputValue('chip_select_all', Math.random(), {priority:'event'});",
                 "All"
               ),
               tags$button(
                 style = paste0("background:none; border:1px solid #ccc; border-radius:8px;",
                                "font-size:10px; color:#555; padding:1px 6px; cursor:pointer;"),
                 onclick = "Shiny.setInputValue('chip_deselect_all', Math.random(), {priority:'event'});",
                 "None"
               )
      ),
      tags$div(style = "display: flex; flex-wrap: wrap;", chips),
      help_card
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

  # ── HTML Report ────────────────────────────────────────────────────────────
  output$download_report <- downloadHandler(
    filename = function() {
      gene <- tryCatch(isolate(filtered_data()$feature), error = function(e) "report")
      fmt  <- tryCatch(isolate(input$report_format), error = function(e) "html")
      ext  <- if (identical(fmt, "pdf")) "pdf" else "html"
      paste0("LongViewSC_", gene, "_", Sys.Date(), ".", ext)
    },
    content = function(file) {
      req(filtered_data(), seurat_joined())
      fmt <- tryCatch(isolate(input$report_format), error = function(e) "html")

      # Transcript plot — re-compute with current GTF & selected isoforms
      transcript_plot_obj <- tryCatch({
        sel     <- isolate(selected_isoforms())
        gtf_val <- isolate(gtf())
        if (!is.null(gtf_val) && is.data.frame(gtf_val) && nrow(gtf_val) > 0 && length(sel) > 0) {
          col_map  <- isolate(isoform_colour_map())
          gtf_keys <- sub("-[^-]+$", "", sel)
          fill_map <- stats::setNames(col_map[sel], gtf_keys)
          p <- plot_gene_transcripts(isoforms_to_plot = sel, gtf = gtf_val)
          p + ggplot2::scale_fill_manual(values = fill_map, na.value = "grey80")
        } else NULL
      }, error = function(e) NULL)

      # isoform_plot() now returns a raw list; apply shared scale + combine for report
      isoform_plot_combined <- tryCatch({
        pl <- isolate(isoform_plot())
        .combine_iso_feature_plot(
          plots_list   = pl,
          shared_scale = isolate(input$shared_scale),
          seu          = isolate(seurat_joined()),
          assay        = isolate(iso_settings()$assay),
          sel_isoforms = isolate(selected_isoforms_debounced()),
          use_viridis  = isTRUE(isolate(input$feature_viridis))
        )
      }, error = function(e) NULL)

      params_list <- list(
        gene                 = isolate(filtered_data()$feature),
        date                 = as.character(Sys.Date()),
        group_by             = isolate(iso_settings()$group_by),
        report_format        = fmt,
        selected_isoforms    = isolate(selected_isoforms()),
        celltype_plot        = isolate(filtered_data()$celltype_plot),
        feature_plot_gene    = isolate(filtered_data()$feature_plot_gene),
        vln_plot             = isolate(filtered_data()$vln_plot),
        isoform_plot         = isoform_plot_combined,
        dotplot_isoform      = tryCatch(isolate(dotplot_isoform()),          error = function(e) NULL),
        transcript_plot      = transcript_plot_obj,
        pie_plot             = tryCatch(isolate(pie_plot()),                  error = function(e) NULL),
        heatmap_plot         = tryCatch(isolate(reactive_heatmap()),          error = function(e) NULL),
        trajectory_gene_plot = tryCatch(isolate(trajectory_gene_plot()),      error = function(e) NULL),
        trajectory_plot      = tryCatch(isolate(trajectory_plot()),           error = function(e) NULL),
        isoform_stats_df     = isolate(filtered_data()$isoform_stats_df),
        scatter_plot         = tryCatch(isolate(scatter_plot_obj()),          error = function(e) NULL),
        scatter_iso_x        = tryCatch(isolate(input$scatter_iso_x),         error = function(e) ""),
        scatter_iso_y        = tryCatch(isolate(input$scatter_iso_y),         error = function(e) ""),
        compare_feature_plot = tryCatch(isolate(compare_feature_plot_obj()),  error = function(e) NULL),
        compare_gene_x       = tryCatch(isolate(input$scatter_gene_x),        error = function(e) ""),
        compare_gene_y       = tryCatch(isolate(input$scatter_gene_y),        error = function(e) "")
      )

      # Use string identifiers so rmarkdown reads settings from the Rmd YAML
      # and writes directly to the path Shiny has allocated (no extension
      # appending).  Format objects can cause rmarkdown to rename the output
      # file, leaving the download handler pointing at a non-existent path.
      out_format <- if (identical(fmt, "pdf")) "pdf_document" else "html_document"

      # Copy Rmd into a dedicated temp dir so figures, .tex, and output all
      # live in the same directory — prevents xelatex "file not found" errors.
      render_dir <- tempfile("lvsc_report_")
      dir.create(render_dir, recursive = TRUE)
      rmd_copy <- file.path(render_dir, "report.Rmd")
      file.copy(normalizePath("app/report_template.Rmd"), rmd_copy)

      session$sendCustomMessage("lvsc_toast", "Report ready — downloading…")
      withProgress(message = "Generating report…", value = 0.5, {
        tryCatch(
          rmarkdown::render(
            input         = rmd_copy,
            output_file   = file,
            output_format = out_format,
            params        = params_list,
            envir         = new.env(parent = globalenv()),
            quiet         = TRUE
          ),
          error = function(e) {
            # Try to grab the LaTeX log for a more useful error message
            msg <- conditionMessage(e)
            tex_log <- tryCatch({
              log_path <- sub("\\.tex$", ".log", regmatches(msg, regexpr("/tmp/[^ ]+\\.tex", msg)))
              if (length(log_path) && nzchar(log_path) && file.exists(log_path)) {
                log_lines <- readLines(log_path, warn = FALSE)
                errors <- grep("^!", log_lines, value = TRUE)
                if (length(errors)) paste(errors, collapse = "\n") else NULL
              } else NULL
            }, error = function(e2) NULL)
            full_msg <- if (!is.null(tex_log)) paste0(msg, "\n\nLaTeX errors:\n", tex_log) else msg
            showNotification(
              paste0("Report generation failed: ", full_msg),
              type = "error", duration = 20
            )
            writeLines(paste0(
              "<html><body><h2>Report generation failed</h2><pre>",
              htmltools::htmlEscape(full_msg),
              "</pre></body></html>"
            ), file)
          }
        )
      })
    }
  )






  # ── Isoform Scatter Tab ───────────────────────────────────────────────────

  # Robust expression-matrix getter: tries layer= (Seurat v5) then slot= (v4),
  # prefers normalised "data" and falls back to "counts".
  .get_expr_matrix <- function(seu, assay_name) {
    for (lyr in c("data", "counts")) {
      mat <- tryCatch(
        GetAssayData(seu, assay = assay_name, layer = lyr),
        error = function(e) tryCatch(
          GetAssayData(seu, assay = assay_name, slot  = lyr),
          error = function(e2) NULL
        )
      )
      if (!is.null(mat)) return(mat)
    }
    stop(paste("Cannot retrieve expression matrix from assay:", assay_name))
  }

  # Gene dropdowns are populated by update_seurat_ui() (same source as main sidebar).
  # Isoform lists update when gene selection changes — exact same logic as filtered_data().
  iso_feats_for_gene <- function(gene) {
    req(seurat_joined(), input$isoform_assay)
    assay_name    <- input$isoform_assay
    isoform_feats <- tryCatch(
      rownames(seurat_joined()@assays[[assay_name]]@features),
      error = function(e) rownames(seurat_joined()@assays[[assay_name]])
    )
    feats <- grep(paste0("(^|-|\\b)", gene, "($|\\b)"), isoform_feats, value = TRUE)
    # Sort by total expression descending (same order as chip selector)
    if (length(feats) > 1) {
      mat   <- tryCatch(.get_expr_matrix(seurat_joined(), assay_name), error = function(e) NULL)
      if (!is.null(mat) && all(feats %in% rownames(mat))) {
        totals <- Matrix::rowSums(mat[feats, , drop = FALSE])
        feats  <- names(sort(totals, decreasing = TRUE))
      }
    }
    feats
  }

  observeEvent(input$scatter_gene_x, {
    req(input$scatter_gene_x)
    feats <- iso_feats_for_gene(input$scatter_gene_x)
    updateSelectInput(session, "scatter_iso_x", choices = feats,
                      selected = if (length(feats) > 0) feats[1] else NULL)
    # Default gene Y to same gene so users compare isoforms within a gene by default
    if (is.null(input$scatter_gene_y) || input$scatter_gene_y == "") {
      updateSelectizeInput(session, "scatter_gene_y", selected = input$scatter_gene_x)
    }
  })

  observeEvent(input$scatter_gene_y, {
    req(input$scatter_gene_y)
    feats <- iso_feats_for_gene(input$scatter_gene_y)
    updateSelectInput(session, "scatter_iso_y", choices = feats,
                      selected = if (length(feats) > 0) feats[1] else NULL)
  })

  # Populate cell filter picker from the currently selected metadata column
  observe({
    req(seurat_joined(), input$group_by)
    grp_vals <- sort(unique(as.character(seurat_joined()@meta.data[[input$group_by]])))
    shinyWidgets::updatePickerInput(session, "scatter_cells_filter",
      choices  = grp_vals,
      selected = grp_vals
    )
  })

  # Populate pie cell filter from the currently selected metadata column
  observe({
    req(seurat_joined(), input$group_by)
    grp_vals <- sort(unique(as.character(seurat_joined()@meta.data[[input$group_by]])))
    shinyWidgets::updatePickerInput(session, "pie_cells_filter",
      choices  = grp_vals,
      selected = grp_vals
    )
  })

  # ── Fetch raw expression — auto-reactive (no button needed) ──────────────
  scatter_raw <- reactive({
    req(seurat_joined(), input$isoform_assay)
    scatter_type <- if (is.null(input$scatter_compare_type)) "isoform" else input$scatter_compare_type
    if (scatter_type == "gene") {
      req(input$scatter_gene_x, input$scatter_gene_y)
      gene_assay <- setdiff(names(seurat_joined()@assays), input$isoform_assay)[1]
      req(!is.na(gene_assay))
      expr_mat <- .get_expr_matrix(seurat_joined(), gene_assay)
      ix <- input$scatter_gene_x
      iy <- input$scatter_gene_y
      validate(
        need(ix %in% rownames(expr_mat), paste("Gene not found in assay:", ix)),
        need(iy %in% rownames(expr_mat), paste("Gene not found in assay:", iy))
      )
    } else {
      req(input$scatter_iso_x, input$scatter_iso_y)
      expr_mat <- .get_expr_matrix(seurat_joined(), input$isoform_assay)
      ix <- input$scatter_iso_x
      iy <- input$scatter_iso_y
      validate(
        need(ix %in% rownames(expr_mat), paste("Isoform not found in assay:", ix)),
        need(iy %in% rownames(expr_mat), paste("Isoform not found in assay:", iy))
      )
    }
    data.frame(
      feat_x = ix, feat_y = iy,
      expr_x = as.numeric(expr_mat[ix, ]),
      expr_y = as.numeric(expr_mat[iy, ]),
      stringsAsFactors = FALSE
    )
  })

  # ── Apply colour + cell filter + zero filter ─────────────────────────────
  scatter_data <- reactive({
    req(scatter_raw())
    df         <- scatter_raw()
    colour_col <- input$group_by
    if (!is.null(colour_col) && nzchar(colour_col) &&
        colour_col %in% colnames(seurat_joined()@meta.data)) {
      df$colour_by <- as.character(seurat_joined()@meta.data[[colour_col]])
    } else {
      df$colour_by <- NA_character_
    }
    validate(
      need(!identical(df$feat_x[1], df$feat_y[1]),
           "Please select two different features to compare.")
    )
    # Filter to selected cell groups
    sel <- input$scatter_cells_filter
    if (!is.null(sel) && length(sel) > 0 && !all(is.na(df$colour_by)))
      df <- df[df$colour_by %in% sel, ]
    if (!isTRUE(input$scatter_zero_cells))
      df <- df[df$expr_x > 0 | df$expr_y > 0, ]
    df
  })

  # ── Build plot ──────────────────────────────────────────────────────────
  scatter_plot_obj <- reactive({
    req(scatter_data())
    df  <- scatter_data()
    ix  <- df$feat_x[1]
    iy  <- df$feat_y[1]
    validate(need(nrow(df) >= 3, "Not enough cells with expression data to plot."))

    lm_fit <- lm(expr_y ~ expr_x, data = df)
    r2     <- round(summary(lm_fit)$r.squared, 3)
    r2_lab <- paste0("R\u00b2 = ", r2, "   n = ",
                     format(nrow(df), big.mark = ","), " cells")
    has_colour <- !all(is.na(df$colour_by))

    p_scatter <- ggplot(df, aes(x = expr_x, y = expr_y)) +
      { if (has_colour)
          geom_point(aes(colour = colour_by), alpha = 0.55, size = 2.0)
        else
          geom_point(colour = "#3d8b6e", alpha = 0.50, size = 2.0) } +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                  colour = "#e05c2a", fill = "#e05c2a", alpha = 0.12, linewidth = 1) +
      annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.8,
               label = r2_lab, size = 4, fontface = "bold", colour = "#333333") +
      labs(x = paste0(ix, "  (log-normalised)"),
           y = paste0(iy, "  (log-normalised)"),
           colour = if (has_colour) input$group_by else NULL) +
      theme_minimal(base_size = 13) +
      theme(legend.position  = if (has_colour) "right" else "none",
            panel.grid.minor = element_blank())

    if (isTRUE(input$scatter_show_marginal)) {
      p_top <- ggplot(df, aes(x = expr_x)) +
        geom_histogram(aes(y = after_stat(density)),
                       bins = 40, fill = "#3d8b6e", alpha = 0.35, colour = NA) +
        geom_density(fill = "#3d8b6e", alpha = 0.25, colour = "#3d8b6e", linewidth = 0.5) +
        theme_void() + theme(plot.margin = margin(2, 0, 2, 0))

      p_right <- ggplot(df, aes(x = expr_y)) +
        geom_histogram(aes(y = after_stat(density)),
                       bins = 40, fill = "#e05c2a", alpha = 0.35, colour = NA) +
        geom_density(fill = "#e05c2a", alpha = 0.25, colour = "#e05c2a", linewidth = 0.5) +
        coord_flip() +
        theme_void() + theme(plot.margin = margin(0, 2, 0, 2))

      (p_top + plot_spacer() + p_scatter + p_right) +
        plot_layout(ncol = 2, nrow = 2, widths = c(5, 1), heights = c(1, 5))
    } else {
      p_scatter
    }
  })

  output$scatter_plot <- renderPlot({
    scatter_type <- if (is.null(input$scatter_compare_type)) "isoform" else input$scatter_compare_type
    if (scatter_type == "gene") {
      if (is.null(input$scatter_gene_x) || !nzchar(input$scatter_gene_x))
        return(.empty_plot("Select Gene X and Gene Y above to compare."))
    } else {
      if (is.null(input$scatter_iso_x) || !nzchar(input$scatter_iso_x))
        return(.empty_plot("Select a gene and isoforms above to compare."))
    }
    scatter_plot_obj()
  }, res = 110)

  downloadModalServer(
    "scatter_plot",
    reactive({ req(scatter_plot_obj()); scatter_plot_obj() }),
    "Isoform scatter"
  )

  # ── Compare Feature Map (custom RGB co-expression) ────────────────────────
  # Per-cell colour: red channel = feature 1, blue channel = feature 2.
  # Purple = co-expression. Grey = neither. No Seurat blend artefacts.
  .rgb_coexpr_panel <- function(seu, feat1, feat2, assay, reduction, label1, label2, title = NULL) {
    emb      <- Seurat::Embeddings(seu, reduction = reduction)
    expr_mat <- tryCatch(
      Seurat::GetAssayData(seu, assay = assay, layer = "data"),
      error = function(e) Seurat::GetAssayData(seu, assay = assay, slot  = "data")
    )
    req(feat1 %in% rownames(expr_mat), feat2 %in% rownames(expr_mat))

    e1 <- as.numeric(expr_mat[feat1, ])
    e2 <- as.numeric(expr_mat[feat2, ])

    norm01 <- function(x) {
      hi <- quantile(x[x > 0], 0.99, na.rm = TRUE)
      if (is.na(hi) || hi == 0) return(rep(0, length(x)))
      pmin(x / hi, 1)
    }
    n1 <- norm01(e1)
    n2 <- norm01(e2)

    # 4-corner bilinear colour blend:
    #   (0,0) = grey (neither), (1,0) = rose-red (feat1 only),
    #   (0,1) = slate-blue (feat2 only), (1,1) = vibrant purple (co-expressed)
    col_00 <- c(0.88, 0.88, 0.88)   # grey
    col_10 <- c(0.80, 0.28, 0.25)   # rose-red
    col_01 <- c(0.25, 0.42, 0.78)   # slate-blue
    col_11 <- c(0.58, 0.08, 0.78)   # vibrant violet-purple

    bilinear <- function(ch) {
      (1-n1)*(1-n2)*col_00[ch] + n1*(1-n2)*col_10[ch] +
      (1-n1)*n2   *col_01[ch] + n1*n2    *col_11[ch]
    }
    cell_col <- rgb(
      red   = pmin(bilinear(1), 1),
      green = pmin(bilinear(2), 1),
      blue  = pmin(bilinear(3), 1),
      maxColorValue = 1
    )

    df <- data.frame(x = emb[, 1], y = emb[, 2],
                     score = n1 + n2, col = cell_col)
    df <- df[order(df$score), ]   # plot low-expr cells first

    ax_labels <- colnames(emb)

    # Legend: arc c1 → col_11 → c2 so the vibrant purple shows at the centre
    legend_n <- 200
    t_seq    <- seq(0, 1, length.out = legend_n)
    n1_leg   <- ifelse(t_seq <= 0.5, 1, 2*(1 - t_seq))
    n2_leg   <- ifelse(t_seq <= 0.5, 2*t_seq, 1)
    leg_col  <- function(ch) {
      (1-n1_leg)*(1-n2_leg)*col_00[ch] + n1_leg*(1-n2_leg)*col_10[ch] +
      (1-n1_leg)*n2_leg    *col_01[ch] + n1_leg*n2_leg    *col_11[ch]
    }
    legend_df <- data.frame(
      x     = t_seq,
      y     = 0.5,
      color = rgb(pmin(leg_col(1),1), pmin(leg_col(2),1), pmin(leg_col(3),1))
    )

    p_legend <- ggplot(legend_df, aes(x, y)) +
      geom_tile(aes(fill = color), height = 1, width = 1 / legend_n) +
      scale_fill_identity() +
      annotate("text", x = 0.08, y = -0.7, label = label1,
               colour = rgb(col_10[1], col_10[2], col_10[3]), size = 3.8, hjust = 0.5,
               fontface = "bold", family = "sans") +
      annotate("text", x = 0.5, y = -0.7, label = "Co-expressed",
               colour = rgb(col_11[1], col_11[2], col_11[3]), size = 3.5, hjust = 0.5,
               fontface = "bold", family = "sans") +
      annotate("text", x = 0.92, y = -0.7, label = label2,
               colour = rgb(col_01[1], col_01[2], col_01[3]), size = 3.8, hjust = 0.5,
               fontface = "bold", family = "sans") +
      scale_x_continuous(expand = c(0.02, 0.02)) +
      scale_y_continuous(limits = c(-1.2, 1)) +
      theme_void() +
      theme(plot.margin = margin(2, 20, 2, 20))

    p_umap <- ggplot(df, aes(x, y)) +
      geom_point(colour = df$col, size = 1.4, alpha = 0.85, stroke = 0) +
      labs(x = ax_labels[1], y = ax_labels[2], title = title) +
      theme_classic(base_size = 13) +
      theme(
        axis.text    = element_blank(),
        axis.ticks   = element_blank(),
        axis.line    = element_blank(),
        panel.border = element_rect(colour = "#cccccc", fill = NA, linewidth = 0.5),
        plot.title   = element_text(size = 14, face = "bold", hjust = 0.5, family = "sans")
      )

    p_umap / p_legend + plot_layout(heights = c(10, 1))
  }

  compare_feature_plot_obj <- reactive({
    req(seurat_joined(), iso_settings())
    req(input$scatter_iso_x, input$scatter_iso_y,
        input$scatter_gene_x, input$scatter_gene_y)

    seu        <- seurat_joined()
    # Filter to selected cell groups (same picker as scatter tab)
    sel <- input$scatter_cells_filter
    grp <- input$group_by
    if (!is.null(sel) && length(sel) > 0 && !is.null(grp) && nzchar(grp) &&
        grp %in% colnames(seu@meta.data)) {
      keep <- colnames(seu)[as.character(seu@meta.data[[grp]]) %in% sel]
      if (length(keep) >= 1) seu <- subset(seu, cells = keep)
    }
    assay_iso  <- iso_settings()$assay
    reduction  <- iso_settings()$reduction
    iso_x      <- input$scatter_iso_x
    iso_y      <- input$scatter_iso_y
    gene_x     <- input$scatter_gene_x
    gene_y     <- input$scatter_gene_y
    diff_genes <- !identical(gene_x, gene_y)

    use_blend       <- isTRUE(input$compare_view == "blend")
    use_categorical <- isTRUE(input$compare_view == "feature")

    .categorical_coexpr <- function(seu, feat1, feat2, assay, reduction, label1, label2, title) {
      emb      <- Seurat::Embeddings(seu, reduction = reduction)
      expr_mat <- tryCatch(
        Seurat::GetAssayData(seu, assay = assay, layer = "data"),
        error = function(e) Seurat::GetAssayData(seu, assay = assay, slot = "data")
      )
      req(feat1 %in% rownames(expr_mat), feat2 %in% rownames(expr_mat))
      e1 <- as.numeric(expr_mat[feat1, ]) > 0
      e2 <- as.numeric(expr_mat[feat2, ]) > 0
      category <- dplyr::case_when(
        e1 & e2   ~ "Both",
        e1 & !e2  ~ paste0("Only ", label1),
        !e1 & e2  ~ paste0("Only ", label2),
        TRUE      ~ "Neither"
      )
      lvls <- c("Both", paste0("Only ", label1), paste0("Only ", label2), "Neither")
      category <- factor(category, levels = lvls)
      pal <- c("Both" = "#e6b800",
               setNames("#2e5ba8", paste0("Only ", label1)),
               setNames("#b83232", paste0("Only ", label2)),
               "Neither" = "#d0d0d0")
      df <- data.frame(x = emb[,1], y = emb[,2], category = category)
      df <- df[order(df$category == "Neither"), ]
      ax <- colnames(emb)
      ggplot2::ggplot(df, ggplot2::aes(x, y, colour = category)) +
        ggplot2::geom_point(size = 1.2, alpha = 0.8, stroke = 0) +
        ggplot2::scale_colour_manual(values = pal, name = NULL,
          guide = ggplot2::guide_legend(override.aes = list(size = 4))) +
        ggplot2::labs(x = ax[1], y = ax[2], title = title) +
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::theme(
          axis.text = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(colour = "#cccccc", fill = NA, linewidth = 0.5),
          plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
          legend.position = "bottom",
          legend.text = ggplot2::element_text(size = 11)
        )
    }

    if (use_blend) {
      blend_type <- if (is.null(input$blend_compare_type)) "isoform" else input$blend_compare_type
      if (blend_type == "gene") {
        validate(
          need(!identical(gene_x, gene_y),
               "Blend mode requires two different genes. Please select a different Gene X and Gene Y.")
        )
        gene_assay <- setdiff(names(seu@assays), assay_iso)[1]
        req(!is.na(gene_assay))
        Seurat::DefaultAssay(seu) <- gene_assay
        p_blend <- Seurat::FeaturePlot(
          seu,
          features  = c(gene_x, gene_y),
          blend     = TRUE,
          cols      = c("grey90", "red", "blue"),
          reduction = reduction,
          blend.threshold = 0,
          pt.size   = 0.6
        )
        p_blend[[1]] <- p_blend[[1]] + ggplot2::ggtitle(gene_x)
        p_blend[[2]] <- p_blend[[2]] + ggplot2::ggtitle(gene_y)
        p_blend[[3]] <- p_blend[[3]] + ggplot2::ggtitle("Blend")
        p_blend
      } else {
        validate(
          need(!identical(iso_x, iso_y),
               "Blend mode requires two different isoforms. Please select a different Isoform X and Isoform Y.")
        )
        Seurat::DefaultAssay(seu) <- assay_iso
        p_blend <- Seurat::FeaturePlot(
          seu,
          features  = c(iso_x, iso_y),
          blend     = TRUE,
          cols      = c("grey90", "red", "blue"),
          reduction = reduction,
          blend.threshold = 0,
          pt.size   = 0.6
        )
        short_x <- sub("^(ENST[0-9]+\\.[0-9]+).*", "\\1", iso_x)
        short_y <- sub("^(ENST[0-9]+\\.[0-9]+).*", "\\1", iso_y)
        p_blend[[1]] <- p_blend[[1]] + ggplot2::ggtitle(short_x)
        p_blend[[2]] <- p_blend[[2]] + ggplot2::ggtitle(short_y)
        p_blend[[3]] <- p_blend[[3]] + ggplot2::ggtitle("Blend")
        p_blend
      }

    } else if (use_categorical) {
      iso_panel <- .categorical_coexpr(seu, iso_x, iso_y, assay_iso, reduction,
                                        sub("-[^-]+$", "", iso_x),
                                        sub("-[^-]+$", "", iso_y),
                                        title = "Isoform co-expression")
      if (diff_genes) {
        gene_assay <- setdiff(names(seu@assays), assay_iso)[1]
        req(!is.na(gene_assay))
        gene_panel <- .categorical_coexpr(seu, gene_x, gene_y, gene_assay, reduction,
                                           gene_x, gene_y, title = "Gene co-expression")
        patchwork::wrap_plots(gene_panel, iso_panel, ncol = 2)
      } else {
        iso_panel
      }
    } else {
      iso_panel <- .rgb_coexpr_panel(seu, iso_x, iso_y, assay_iso, reduction,
                                      sub("-[^-]+$", "", iso_x),
                                      sub("-[^-]+$", "", iso_y),
                                      title = "Isoform co-expression")
      if (diff_genes) {
        gene_assay <- setdiff(names(seu@assays), assay_iso)[1]
        req(!is.na(gene_assay))
        gene_panel <- .rgb_coexpr_panel(seu, gene_x, gene_y, gene_assay, reduction,
                                         gene_x, gene_y, title = "Gene co-expression")
        patchwork::wrap_plots(gene_panel, iso_panel, ncol = 2)
      } else {
        iso_panel
      }
    }
  })

  output$compare_feature_plot <- renderPlot({
    if (is.null(input$scatter_iso_x) || !nzchar(input$scatter_iso_x)) {
      return(.empty_plot("Select a gene and isoforms above to compare."))
    }
    compare_feature_plot_obj()
  }, res = 110)

  downloadModalServer(
    "compare_feature_plot",
    reactive({ req(compare_feature_plot_obj()); compare_feature_plot_obj() }),
    "Feature map"
  )


}
