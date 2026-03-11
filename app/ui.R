ui <- fluidPage(
  useShinyjs(),
  
  # Add custom styles for a header banner and plot titles
  tags$head(
    tags$style(HTML("
    /* Hide the blue background on a selected row—only the checkbox is visible */
    table.dataTable tr.selected td,
    table.dataTable td.selected {
      background-color: transparent !important;
    }
    .select-checkbox {
      text-align: center;
      color: #444;
    }
    /* ── Isoform chip selector ───────────────────────────────────────── */
    #isoform_chips_panel {
      margin-top: 12px;
    }
    /* Sidebar sections */
    .sidebar-section {
      background: #f8f9fa;
      border-left: 3px solid #3d8b6e;
      border-radius: 4px;
      padding: 10px 12px 6px 12px;
      margin-bottom: 10px;
    }
    .sidebar-section .form-group { margin-bottom: 6px; }
    .sidebar-section-label {
      font-size: 11px;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.06em;
      color: #3d8b6e;
      margin: 0 0 7px 0;
    }
    .sidebar-section-label .fa, .sidebar-section-label .glyphicon { margin-right: 4px; }
    .chip-panel-title {
      font-size: 12px;
      font-weight: bold;
      color: #555;
      text-transform: uppercase;
      letter-spacing: 0.04em;
      margin-bottom: 8px;
    }
    .input-label-row {
      display: flex;
      align-items: center;
      gap: 5px;
      margin-bottom: 2px;
    }
    .input-label-row label {
      margin: 0;
      font-size: 13px;
      font-weight: 500;
      color: #333;
    }
    .chip-info-icon {
      display: inline-flex;
      align-items: center;
      justify-content: center;
      width: 16px;
      height: 16px;
      border-radius: 50%;
      background: #aaa;
      color: #fff;
      font-size: 10px;
      font-weight: 700;
      font-style: normal;
      font-family: sans-serif;
      cursor: help;
      line-height: 1;
      opacity: 0.7;
      flex-shrink: 0;
    }
    .chip-info-icon:hover { opacity: 1; background: #3d8b6e; }
    .iso-chip {
      display: inline-flex;
      flex-direction: column;
      align-items: flex-start;
      margin: 3px 3px;
      padding: 5px 8px 4px 8px;
      border-radius: 16px;
      border: 1.5px solid #ccc;
      background: #f0f0f0;
      color: #888;
      font-size: 11px;
      font-family: monospace;
      cursor: pointer;
      transition: all 0.15s ease;
      min-width: 80px;
      max-width: 100%;
      user-select: none;
    }
    .iso-chip:hover {
      border-color: #999;
      background: #e4e4e4;
    }
    .iso-chip.chip-selected {
      color: white;
      border-color: transparent;
      box-shadow: 0 1px 4px rgba(0,0,0,0.18);
    }
    .iso-chip.chip-selected:hover {
      filter: brightness(0.92);
    }
    .chip-label {
      white-space: nowrap;
      overflow: hidden;
      text-overflow: ellipsis;
      max-width: 150px;
      display: block;
    }
    .chip-bar-wrap {
      width: 100%;
      height: 3px;
      background: rgba(255,255,255,0.35);
      border-radius: 2px;
      margin-top: 4px;
    }
    .chip-bar {
      height: 3px;
      border-radius: 2px;
      background: rgba(255,255,255,0.85);
      min-width: 2px;
    }
    .iso-chip:not(.chip-selected) .chip-bar {
      background: #bbb;
    }
    /* ── Bucket list styling for trajectory tab ────────────────────── */
    /* Bucket list styling for trajectory tab */
    .rank-list-container {
      min-height: 60px;
      border: 2px dashed #ccc;
      border-radius: 6px;
      padding: 6px;
      background: #fafafa;
    }
    .rank-list-item {
      background: white;
      border: 1px solid #ddd;
      border-radius: 4px;
      padding: 5px 10px;
      margin: 3px 0;
      cursor: grab;
      font-size: 13px;
    }
    .rank-list-item:hover { border-color: #2c3e50; background: #eaf0fb; }
    .rank-list-title { font-size: 13px; font-weight: bold; margin-bottom: 4px; }
    /* Horizontal drag list for trajectory ordering */
    #traj_cell_order.rank-list {
      display: flex;
      flex-direction: row;
      flex-wrap: wrap;
      gap: 4px;
      min-height: 36px;
      padding: 4px;
    }
    #traj_cell_order .rank-list-item {
      margin: 0;
      padding: 4px 10px;
      font-size: 12px;
      cursor: grab;
    }
      .banner {
        background-color: #2c3e50; 
        color: white; 
        padding: 15px; 
        text-align: center; 
        font-size: 24px; 
        font-weight: bold;
        border-radius: 10px;
        margin-bottom: 20px;
      }
      .sidebar-panel {
        background-color: #f8f9fa;
        padding: 15px;
        border-radius: 10px;
      }
      .plot-box {
        border: 2px solid #ddd;
        padding: 10px;
        border-radius: 10px;
        box-shadow: 2px 2px 10px rgba(0, 0, 0, 0.1);
        background-color: white;
        margin-bottom: 20px;
      }
      .plot-title {
        text-align: center;
        font-weight: bold;
        font-size: 18px;
        margin-bottom: 10px;
      }
      .mailbox-button {
        position: absolute;
        top: 15px;
        right: 20px;
        background: #2c3e50;
        border: none;
        color: white;
        padding: 10px;
        border-radius: 50%;
        font-size: 18px;
        cursor: pointer;
      }
          #landingPage {
      position: fixed; 
      top: 0; left: 0; width: 100%; height: 100%;
      background-color: #2c3e50; 
      display: flex; align-items: center; justify-content: center;
      flex-direction: column;
      color: white; text-align: center;
      z-index: 9999; /* Ensures it's on top */
    }

    /* Tab panel custom colour */
    .nav-tabs > li > a {
      color: #337ab7;
      font-weight: 500;
    }
    .nav-tabs > li.active > a,
    .nav-tabs > li.active > a:hover,
    .nav-tabs > li.active > a:focus {
      color: #fff;
      background-color: #337ab7;
      border-color: #337ab7;
    }
    .nav-tabs > li > a:hover {
      background-color: #e8f0fb;
      border-color: #337ab7;
    }
    .nav-tabs { border-bottom: 2px solid #337ab7; }
    /* Sidebar collapse */
    #sidebarToggle { transition: left 0.25s ease; }
    /* Collapsible plot boxes */
    .plot-title-row .collapse-chevron {
      margin-left: auto;
      font-size: 14px;
      color: #aaa;
      transition: transform 0.2s ease;
      flex-shrink: 0;
    }
    .plot-box.plot-collapsed .plot-body { display: none !important; }
    .plot-box.plot-collapsed .collapse-chevron { transform: rotate(-90deg); }
    .plot-title-row:hover .collapse-chevron { color: #337ab7; }
    #sidebarToggle {
      position: fixed;
      top: 80px;
      left: 12px;
      z-index: 9999;
      background: #2c3e50;
      color: white;
      border: none;
      border-radius: 50%;
      width: 32px;
      height: 32px;
      font-size: 14px;
      cursor: pointer;
      box-shadow: 1px 1px 6px rgba(0,0,0,0.3);
      display: flex;
      align-items: center;
      justify-content: center;
      padding: 0;
    }
    #sidebarToggle:hover { background: #337ab7; }
    "))
  ),
  
  
  # Landing Page
  div(id = "landingPage",
      h1(HTML("Welcome to LongViewSC: <br> The Long-read Single-cell Visualization App")),
      br(),
      # Wrap buttons in a flexbox container
      div(style = "display: flex; gap: 10px; align-items: center;",
          
          # Start Analysis Button
          actionButton("startBtn", "Start Analysis", class = "btn btn-primary btn-lg"),
          
          # Instructions Button
          actionButton("instructionsBtn", "Instructions", class = "btn btn-secondary btn-lg"),
          
          # Demo Data Button
          actionButton("DemoBtn", "Demo Data", class = "btn btn-success btn-lg"),
          
          # Mailbox/Support Button
          div(class = "dropdown",
              style = "margin-left: 5px;",
              tags$button(class = "btn btn-info btn-lg dropdown-toggle", 
                          type = "button", 
                          `data-toggle` = "dropdown",
                          HTML('<i class="fa fa-envelope"></i> Support')
              ),
              tags$ul(class = "dropdown-menu",
                      tags$li(HTML('<li><a href="https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/clark-lab" target="_blank"><i class="fa fa-user"></i> Clark Laboratory</a></li>')),
                      tags$li(HTML('<li><a href="mailto:sefi.prawer@unimelb.edu.au" target="_blank"><i class="fa fa-question-circle"></i> sefi.prawer@unimelb.edu.au</a></li>')),
                      tags$li(HTML('<li><a href="https://github.com/Sefi196/LongViewSC" target="_blank"><i class="fa fa-github"></i> GitHub repo </a></li>'))
                      
              )
          )
      ),
      br(),
      br(),
      h5("Developed by Sefi Prawer at the University of Melbourne")
  ),
  
  # Banner Section with Back Arrow and mail drop down menu  
  div(class = "banner", 
      div(style = "display: flex; align-items: center; justify-content: space-between; width: 100%;",
          
          # Back Button
          actionButton("backToLanding", label = "", icon = icon("arrow-left"), class = "btn btn-dark btn-lg"),
          
          # Title in Banner
          span("View Single Cell Isoform Expression", 
               style = "flex-grow: 1; text-align: center; font-size: 24px; font-weight: bold;"),
          
          # Mailbox Button (Dropdown Menu)
          div(class = "dropdown",
              style = "margin-left: auto;",  # Pushes it to the right
              tags$button(class = "btn btn-dark btn-lg", 
                          type = "button", 
                          `data-toggle` = "dropdown",
                          icon("envelope", style = "color: black;")  # Add style to make icon black
              ),
              tags$ul(class = "dropdown-menu dropdown-menu-right",
                      tags$li(HTML('<li><a href="https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/clark-lab" target="_blank"><i class="fa fa-user"></i> Clark Laboratory</a></li>')),
                      tags$li(HTML('<li><a href="mailto:sefi.prawer@unimelb.edu.au" target="_blank"><i class="fa fa-question-circle"></i> sefi.prawer@unimelb.edu.au</a></li>')),
                      tags$li(HTML('<li><a href="https://github.com/Sefi196/LongViewSC" target="_blank"><i class="fa fa-github"></i> GitHub repo </a></li>'))
              )
          )
      )
  ),
  
  # Spinner element (hidden by default)
  div(id = "spinner", 
      img(src = "https://i.gifer.com/ZZ5H.gif", height = "50px", width = "50px"),
      style = "display: none; 
               position: fixed; 
               bottom: 20px; 
               left: 50%; 
               transform: translateX(-50%); 
               text-align: center; 
               background-color: white; 
               padding: 20px; 
               border-radius: 8px; 
               box-shadow: 0 4px 8px rgba(0, 0, 0, 0.3); 
               z-index: 9999;"),
  
  
  div(id = "instructionsPage", style = "display: none;", 
      fluidPage(
        column(12,
               # Full-width image display
               tags$img(src = "Readme_home_app.png", style = "width: 100%; height: auto; display: block;"),
               includeMarkdown("README.md")
        )
      )
  ),
  
  
  # Main UI - Initially Hidden
  div(id = "mainUI",
      
      tags$script(HTML('
    $(document).ready(function() {
      // Bootstrap popovers for chip info icons
      $(document).on("mouseenter focus", "[data-toggle=popover]", function() {
        $(this).popover("show");
      }).on("mouseleave blur", "[data-toggle=popover]", function() {
        $(this).popover("hide");
      });
      // Collapsible plot boxes
      function addChevrons() {
        $(".plot-box .plot-title-row").each(function() {
          if ($(this).find(".collapse-chevron").length === 0) {
            $(this).append($("<span>").addClass("collapse-chevron").html("&#9662;"));
          }
        });
      }
      // Run after Shiny finishes each render cycle
      $(document).on("shiny:idle", addChevrons);
      // Collapse on title click
      $(document).on("click", ".plot-box .plot-title-row", function(e) {
        if ($(e.target).closest(".chip-info-icon").length) return;
        var $box = $(this).closest(".plot-box");
        var collapsing = !$box.hasClass("plot-collapsed");
        $box.toggleClass("plot-collapsed");
        $(this).nextAll().toggle(!collapsing);
      });
      Shiny.addCustomMessageHandler("reinit_popovers", function(x) {
        $("[data-toggle=popover]").popover();
      });
      // Sidebar collapse
      $(document).on("click", "#sidebarToggle", function() {
        Shiny.setInputValue("sidebarToggle_click", Math.random());
        // Use sidebar sibling — mainPanel id lands on inner div not the col
        var $sidebarCol = $(".sidebar-panel").parent();
        var $mainCol    = $sidebarCol.siblings().last();
        var icon = $(this).find("i");
        if (icon.hasClass("fa-chevron-left")) {
          icon.removeClass("fa-chevron-left").addClass("fa-chevron-right");
          $mainCol.removeClass("col-sm-8").addClass("col-sm-12");
        } else {
          icon.removeClass("fa-chevron-right").addClass("fa-chevron-left");
          $mainCol.removeClass("col-sm-12").addClass("col-sm-8");
        }
      });
    });
  ')),
      sidebarLayout(
        sidebarPanel(
          width = 4,
          id = "sidebar_panel_wrap",
          class = "sidebar-panel",
          
          # (1) File inputs — wrapped in one div so demo mode can hide everything cleanly
          div(id = "file_inputs_panel",
              # --- Seurat object ---
              tags$label("Seurat Object", style = "font-weight: bold;"),
              tags$small("Choose one method:", style = "color: #666;"),
              tags$br(),
              fileInput("seurat_file", "Upload from local machine (.rds / .qs / .qs2)",
                        accept = c(".rds", ".qs", ".qs2")),
              div(style = "margin-top: -10px; margin-bottom: 6px;",
                  tags$em("— or —", style = "color: #999; font-size: 13px;")
              ),
              div(id = "hpc_seurat_inputs",
                  div(style = "display: flex; gap: 6px; align-items: flex-end; margin-bottom: 15px;",
                      div(style = "flex: 1;",
                          textInput("seurat_hpc_path", "Load from HPC path:",
                                    placeholder = "/path/to/object.qs2")
                      ),
                      actionButton("load_hpc_seurat", "Load", class = "btn btn-primary btn-sm",
                                   style = "margin-bottom: 1px;")
                  )
              ),
              tags$hr(),
              # --- GTF ---
              tags$label("GTF File", style = "font-weight: bold;"),
              tags$small("Choose one method:", style = "color: #666;"),
              tags$br(),
              fileInput("gtf", "Upload from local machine (.gtf)", accept = ".gtf"),
              div(style = "margin-top: -10px; margin-bottom: 6px;",
                  tags$em("— or —", style = "color: #999; font-size: 13px;")
              ),
              div(id = "hpc_gtf_inputs",
                  div(style = "display: flex; gap: 6px; align-items: flex-end; margin-bottom: 15px;",
                      div(style = "flex: 1;",
                          textInput("gtf_hpc_path", "Load from HPC path:",
                                    placeholder = "/path/to/annotation.gtf")
                      ),
                      actionButton("load_hpc_gtf", "Load", class = "btn btn-primary btn-sm",
                                   style = "margin-bottom: 1px;")
                  )
              ),
              tags$hr()
          ),
          
          # Prompt shown before Seurat is loaded
          div(id = "analysisForm",
              
              # ── Section 1: Gene ────────────────────────────────────
              tags$div(class = "sidebar-section",
                       tags$p(class = "sidebar-section-label", "Gene"),
                       selectizeInput("feature", NULL, choices = NULL,
                                      options = list(placeholder = "Start typing gene name…",
                                                     maxOptions  = 3000))
              ),
              
              # ── Section 2: Plotting options ──────────────────────────
              tags$div(class = "sidebar-section",
                       tags$p(class = "sidebar-section-label", "Plotting Options"),
                       
                       # Reduction
                       tags$div(class = "input-label-row",
                                tags$label("Reduction", `for` = "reduction"),
                                tags$span(class = "chip-info-icon",
                                          `data-toggle` = "popover", `data-trigger` = "hover focus",
                                          `data-placement` = "right",
                                          `data-content` = "Controls which embedding is shown in the feature and cell-type plots (e.g. UMAP, tSNE).",
                                          tabindex = "0", "i")
                       ),
                       selectInput("reduction", NULL, choices = NULL),
                       
                       # Isoform Assay
                       tags$div(class = "input-label-row",
                                tags$label("Isoform Assay", `for` = "isoform_assay"),
                                tags$span(class = "chip-info-icon",
                                          `data-toggle` = "popover", `data-trigger` = "hover focus",
                                          `data-placement` = "right",
                                          `data-content` = "The Seurat assay containing isoform expression data (e.g. iso, tx_counts).",
                                          tabindex = "0", "i")
                       ),
                       selectInput("isoform_assay", NULL, choices = NULL),
                       
                       # Metadata Column
                       tags$div(class = "input-label-row",
                                tags$label("Metadata Column", `for` = "group_by"),
                                tags$span(class = "chip-info-icon",
                                          `data-toggle` = "popover", `data-trigger` = "hover focus",
                                          `data-placement` = "right",
                                          `data-content` = "A categorical metadata column to group cells. Used for colouring plots (e.g. cell type, cluster).",
                                          tabindex = "0", "i")
                       ),
                       selectInput("group_by", NULL, choices = NULL)
              ),
              
              # ── Run / Reset ───────────────────────────────────────────
              tags$div(style = "display: flex; gap: 8px; margin-bottom: 10px;",
                       actionButton("GO",       "Run Analysis", icon = icon("play"),
                                    class = "btn btn-success",
                                    style = "flex: 2; font-size: 14px; font-weight: 600;"),
                       actionButton("resetBtn", "Reset", icon = icon("undo"),
                                    class = "btn btn-secondary",
                                    style = "flex: 1; font-size: 13px; background: #bbb; border-color: #bbb; color: #fff;")
              ),
              
              # ── Section 3: Isoform chip selector ──────────────────────────
              uiOutput("isoform_chips_panel")
              
          )  # end of analysisForm
        ),
        
        mainPanel(
          id = "main_panel_wrap",
          tags$button(id = "sidebarToggle", title = "Collapse sidebar",
                      tags$i(class = "fa fa-chevron-left")),
          tabsetPanel(id = "main_tabs",
                      
                      # ── Tab 1: Overview ───────────────────────────────────────
                      tabPanel("Overview",
                               fluidRow(
                                 column(6,
                                        div(class = "plot-box",
                                            tags$div(class = "plot-title-row", style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:10px; cursor:pointer; user-select:none;",
                                                     h3(class = "plot-title", style="margin:0;", "Cell Types"),
                                                     tags$span(class = "chip-info-icon",
                                                               `data-toggle` = "popover",
                                                               `data-trigger` = "hover focus",
                                                               `data-placement` = "right",
                                                               `data-content` = "Cells plotted on the selected embedding (e.g. UMAP), coloured by the chosen Metadata Column.",
                                                               tabindex = "0", "i")
                                            ),
                                            shinycssloaders::withSpinner(plotOutput("celltype_plot"), type=4, color="#337ab7", size=0.7),
                                            div(style = "text-align:right; padding-right:10px; margin-top:4px;",
                                                actionButton("celltype_plot-download_celltype_plot",
                                                             "Download", icon = icon("download"),
                                                             class = "btn btn-sm btn-link",
                                                             style = "font-size:12px; color:#337ab7;"))
                                        )
                                 ),
                                 column(6,
                                        div(class = "plot-box",
                                            tags$div(class = "plot-title-row", style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:10px; cursor:pointer; user-select:none;",
                                                     h3(class = "plot-title", style="margin:0;", "Gene Expression"),
                                                     tags$span(class = "chip-info-icon",
                                                               `data-toggle` = "popover",
                                                               `data-trigger` = "hover focus",
                                                               `data-placement` = "right",
                                                               `data-content` = "Expression of the selected gene overlaid on the embedding. Brighter colour = higher expression.",
                                                               tabindex = "0", "i")
                                            ),
                                            shinycssloaders::withSpinner(plotOutput("feature_plot_gene"), type=4, color="#337ab7", size=0.7),
                                            div(style = "text-align:right; padding-right:10px; margin-top:4px;",
                                                actionButton("feature_plot_gene-download_feature_plot_gene",
                                                             "Download", icon = icon("download"),
                                                             class = "btn btn-sm btn-link",
                                                             style = "font-size:12px; color:#337ab7;"))
                                        )
                                 )
                               ),
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            tags$div(class = "plot-title-row", style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:10px; cursor:pointer; user-select:none;",
                                                     h3(class = "plot-title", style="margin:0;", "Gene Expression by Cell Type"),
                                                     tags$span(class = "chip-info-icon",
                                                               `data-toggle` = "popover",
                                                               `data-trigger` = "hover focus",
                                                               `data-placement` = "right",
                                                               `data-content` = "Distribution of gene expression across cells, grouped by the Metadata Column. Useful for comparing expression levels between cell types.",
                                                               tabindex = "0", "i")
                                            ),
                                            shinycssloaders::withSpinner(plotOutput("vln_plot"), type=4, color="#337ab7", size=0.7),
                                            div(style = "text-align:right; padding-right:10px; margin-top:4px;",
                                                actionButton("vln_plot-download_vln_plot",
                                                             "Download", icon = icon("download"),
                                                             class = "btn btn-sm btn-link",
                                                             style = "font-size:12px; color:#337ab7;"))
                                        )
                                 )
                               )
                      ),
                      
                      # ── Tab 2: Isoform Statistics ─────────────────────────────
                      tabPanel("Isoform Statistics",
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            tags$div(class = "plot-title-row", style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:10px; cursor:pointer; user-select:none;",
                                                     h3(class = "plot-title", style="margin:0;", "Isoform Statistics"),
                                                     tags$span(class = "chip-info-icon",
                                                               `data-toggle` = "popover",
                                                               `data-trigger` = "hover focus",
                                                               `data-placement` = "right",
                                                               `data-content` = "Shows all isoforms for the selected gene. Isoforms currently selected via the sidebar chips are marked with a tick.",
                                                               tabindex = "0", "i")
                                            ),
                                            DT::dataTableOutput("isoform_table")
                                        )
                                 )
                               )
                      ),
                      
                      # ── Tab 3: Isoform Maps ───────────────────────────────────
                      tabPanel("Isoform Expression",
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            tags$div(class = "plot-title-row", style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:10px; cursor:pointer; user-select:none;",
                                                     h3(class = "plot-title", style="margin:0;", "Isoform Feature Plots"),
                                                     tags$span(class = "chip-info-icon",
                                                               `data-toggle` = "popover",
                                                               `data-trigger` = "hover focus",
                                                               `data-placement` = "right",
                                                               `data-content` = "One plot per selected isoform, showing its expression on the embedding. Toggle chips in the sidebar to add or remove isoforms.",
                                                               tabindex = "0", "i")
                                            ),
                                            shinycssloaders::withSpinner(plotOutput("feature_plot_iso"), type=4, color="#337ab7", size=0.7),
                                            div(style = "text-align:right; padding-right:10px; margin-top:4px;",
                                                actionButton("feature_plot_isoform-download_feature_plot_isoform",
                                                             "Download", icon = icon("download"),
                                                             class = "btn btn-sm btn-link",
                                                             style = "font-size:12px; color:#337ab7;"))
                                        )
                                 )
                               ),
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            tags$div(class = "plot-title-row", style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:10px; cursor:pointer; user-select:none;",
                                                     h3(class = "plot-title", style="margin:0;", "Dot Plot"),
                                                     tags$span(class = "chip-info-icon",
                                                               `data-toggle` = "popover",
                                                               `data-trigger` = "hover focus",
                                                               `data-placement` = "right",
                                                               `data-content` = "Dot size = fraction of cells expressing the isoform. Dot colour = average normalised expression. Each row is a cell type; each column is a selected isoform.",
                                                               tabindex = "0", "i")
                                            ),
                                            shinycssloaders::withSpinner(plotOutput("dot_plot_iso"), type=4, color="#337ab7", size=0.7),
                                            div(style = "text-align:right; padding-right:10px; margin-top:4px;",
                                                actionButton("dot_plot_isoform-download_dot_plot_isoform",
                                                             "Download", icon = icon("download"),
                                                             class = "btn btn-sm btn-link",
                                                             style = "font-size:12px; color:#337ab7;"))
                                        )
                                 )
                               )
                      ),
                      
                      # ── Tab 4: Transcript Structure ───────────────────────────
                      tabPanel("Transcript Structure",
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            tags$div(class = "plot-title-row", style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:10px; cursor:pointer; user-select:none;",
                                                     h3(class = "plot-title", style="margin:0;", "Isoform Transcript Structure"),
                                                     tags$span(class = "chip-info-icon",
                                                               `data-toggle` = "popover",
                                                               `data-trigger` = "hover focus",
                                                               `data-placement` = "right",
                                                               `data-content` = "Exon-intron structure of selected isoforms drawn from the GTF annotation. Colours match the sidebar chips. Requires a GTF file to be loaded.",
                                                               tabindex = "0", "i")
                                            ),
                                            shinycssloaders::withSpinner(plotOutput("transcript_plot"), type=4, color="#337ab7", size=0.7),
                                            div(style = "text-align:right; padding-right:10px; margin-top:4px;",
                                                actionButton("Isoform_TranscriptStructure-download_Isoform_TranscriptStructure",
                                                             "Download", icon = icon("download"),
                                                             class = "btn btn-sm btn-link",
                                                             style = "font-size:12px; color:#337ab7;"))
                                        )
                                 )
                               ),
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            tags$div(class = "plot-title-row", style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:10px; cursor:pointer; user-select:none;",
                                                     h3(class = "plot-title", style="margin:0;", "Pseudobulk Expression Heatmap"),
                                                     tags$span(class = "chip-info-icon",
                                                               `data-toggle` = "popover",
                                                               `data-trigger` = "hover focus",
                                                               `data-placement` = "right",
                                                               `data-content` = "Mean normalised expression of each selected isoform, aggregated (pseudobulked) per cell type. Useful for comparing isoform usage across cell populations.",
                                                               tabindex = "0", "i")
                                            ),
                                            shinycssloaders::withSpinner(plotlyOutput("heatmap_plot"), type=4, color="#337ab7", size=0.7),
                                            div(style = "text-align:right; padding-right:10px; margin-top:4px;",
                                                actionButton("pseudobulk_heatmap-download_pseudobulk_heatmap",
                                                             "Download", icon = icon("download"),
                                                             class = "btn btn-sm btn-link",
                                                             style = "font-size:12px; color:#337ab7;"))
                                        )
                                 )
                               )
                      ),
                      
                      # ── Tab 5: Proportions ────────────────────────────────────
                      tabPanel("Proportions",
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            div(style = "display:flex; align-items:flex-end; gap:24px; flex-wrap:wrap; margin-bottom:4px;",
                                                div(
                                                  tags$label("Columns per row",
                                                             style = "font-size:13px; font-weight:500; display:block; margin-bottom:2px;"),
                                                  numericInput("pie_ncol", NULL, value = 4, min = 1, max = 10, step = 1, width = "90px")
                                                ),
                                                div(
                                                  tags$label("Min. counts threshold",
                                                             style = "font-size:13px; font-weight:500; display:block; margin-bottom:2px;"),
                                                  numericInput("pie_min_counts", NULL, value = 10, min = 0, step = 1, width = "90px")
                                                ),
                                                tags$p(style = "font-size:12px; color:#888; padding-bottom:6px; margin:0;",
                                                       "Cell types below the threshold are hidden.",
                                                       tags$br(),
                                                       "Grouping uses the Metadata Column.")
                                            )
                                        )
                                 )
                               ),
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            tags$div(class = "plot-title-row", style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:10px; cursor:pointer; user-select:none;",
                                                     h3(class = "plot-title", style="margin:0;", "Isoform Proportions by Cell Type"),
                                                     tags$span(class = "chip-info-icon",
                                                               `data-toggle` = "popover",
                                                               `data-trigger` = "hover focus",
                                                               `data-placement` = "right",
                                                               `data-content` = "Each pie shows the relative proportion of selected isoforms within a cell type. Only cell types above the minimum counts threshold are shown.",
                                                               tabindex = "0", "i")
                                            ),
                                            shinycssloaders::withSpinner(plotOutput("pie_plot", height = "600px"), type=4, color="#337ab7", size=0.7),
                                            div(style = "text-align:right; padding-right:10px; margin-top:4px;",
                                                actionButton("isoform_pie-download_isoform_pie",
                                                             "Download", icon = icon("download"),
                                                             class = "btn btn-sm btn-link",
                                                             style = "font-size:12px; color:#337ab7;"))
                                        )
                                 )
                               )
                      ),
                      
                      # ── Tab 6: Trajectory ─────────────────────────────────────
                      tabPanel("Trajectory",
                               
                               # ── Row 1: settings (left) + gene plot (right) ───
                               fluidRow(
                                 # Settings — stacked vertically, compact
                                 column(4,
                                        div(class = "plot-box",
                                            tags$div(class = "plot-title-row", style = "display:flex; align-items:center; gap:6px; margin-bottom:10px; cursor:pointer; user-select:none;",
                                                     h4("Trajectory Settings", style = "font-weight:600; margin:0;"),
                                                     tags$span(class = "chip-info-icon",
                                                               `data-toggle` = "popover", `data-trigger` = "hover focus",
                                                               `data-placement` = "right",
                                                               `data-content` = "Configure which cell types to include and their ordering. Uses the Metadata Column from the sidebar. Takes effect on next GO press.",
                                                               tabindex = "0", "i")
                                            ),
                                            tags$div(class = "input-label-row",
                                                     tags$label("1. Select cell types", style = "font-size:13px; font-weight:500;"),
                                                     tags$span(class = "chip-info-icon",
                                                               `data-toggle` = "popover", `data-trigger` = "hover focus",
                                                               `data-placement` = "right",
                                                               `data-content` = "Choose which cell types to include. Remove any you don't want to show.",
                                                               tabindex = "0", "i")
                                            ),
                                            uiOutput("traj_include_ui"),
                                            tags$hr(style = "margin:8px 0;"),
                                            tags$div(class = "input-label-row",
                                                     tags$label("2. Drag to set order", style = "font-size:13px; font-weight:500;"),
                                                     tags$span(class = "chip-info-icon",
                                                               `data-toggle` = "popover", `data-trigger` = "hover focus",
                                                               `data-placement` = "right",
                                                               `data-content` = "Drag cell types left-to-right to set their order on the plot x-axis.",
                                                               tabindex = "0", "i")
                                            ),
                                            uiOutput("traj_order_ui")
                                        )
                                 ),
                                 # Gene plot — compact, right side
                                 column(8,
                                        div(class = "plot-box",
                                            tags$div(class = "plot-title-row", style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:10px; cursor:pointer; user-select:none;",
                                                     h3(class = "plot-title", style="margin:0;", "Gene Expression Across Cell States"),
                                                     tags$span(class = "chip-info-icon",
                                                               `data-toggle` = "popover",
                                                               `data-trigger` = "hover focus",
                                                               `data-placement` = "right",
                                                               `data-content` = "Gene-level expression trajectory. Each point is a single cell, coloured by cell type. Uses the RNA assay and the cell ordering set in the settings panel.",
                                                               tabindex = "0", "i")
                                            ),
                                            shinycssloaders::withSpinner(plotOutput("trajectory_gene_plot", height = "300px"), type=4, color="#337ab7", size=0.7),
                                            div(style = "text-align:right; padding-right:10px; margin-top:4px;",
                                                actionButton("trajectory_gene_plot-download_trajectory_gene_plot",
                                                             "Download", icon = icon("download"),
                                                             class = "btn btn-sm btn-link",
                                                             style = "font-size:12px; color:#337ab7;"))
                                        )
                                 )
                               ),
                               
                               # ── Row 2: isoform plot (full width) ─────────────
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            tags$div(class = "plot-title-row", style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:10px; cursor:pointer; user-select:none;",
                                                     h3(class = "plot-title", style="margin:0;", "Isoform Expression Across Cell States"),
                                                     tags$span(class = "chip-info-icon",
                                                               `data-toggle` = "popover",
                                                               `data-trigger` = "hover focus",
                                                               `data-placement` = "right",
                                                               `data-content` = "Each point is a single cell, coloured by cell type and ordered left-to-right by the cell state order set above. The curved line shows the smoothed expression trend.",
                                                               tabindex = "0", "i")
                                            ),
                                            shinycssloaders::withSpinner(plotOutput("trajectory_plot", height = "500px"), type=4, color="#337ab7", size=0.7),
                                            div(style = "text-align:right; padding-right:10px; margin-top:4px;",
                                                actionButton("trajectory_plot-download_trajectory_plot",
                                                             "Download", icon = icon("download"),
                                                             class = "btn btn-sm btn-link",
                                                             style = "font-size:12px; color:#337ab7;"))
                                        )
                                 )
                               )
                      )
                      
                      
          ) # end of tabsetPanel
        ) # end of mainPanel
      ) # end of sidebarLayout
  ) # end of mainUI div
)