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
        position: relative;
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
    #sidebarToggle {
      position: fixed;
      bottom: 20px;
      left: 12px;
      z-index: 9999;
      background: #f8f9fa;
      color: #555;
      border: 1px solid #ddd;
      border-radius: 50%;
      width: 32px;
      height: 32px;
      font-size: 14px;
      cursor: pointer;
      box-shadow: 1px 1px 6px rgba(0,0,0,0.15);
      display: flex;
      align-items: center;
      justify-content: center;
      padding: 0;
    }
    #sidebarToggle:hover { background: #e2e6ea; color: #337ab7; }
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
    /* ── Tutorial / Instructions page ──────────────────────────────────── */
    #instructionsPage { background: #fff; padding: 0 !important; }
    html { scroll-behavior: smooth; }
    .wf-header {
      background: #f8f9fa; padding: 28px 30px 24px; border-bottom: 1px solid #e9ecef;
      text-align: center;
    }
    .wf-header h1 { font-size: 2em; font-weight: 700; color: #2c3e50; margin: 0 0 6px; }
    .wf-header > p { font-size: 1em; color: #666; margin: 0 0 18px; }
    .wf-header-btns { display: flex; gap: 8px; justify-content: center; flex-wrap: wrap; }
    .btn-wf-outline {
      background: transparent; color: #3d8b6e !important; padding: 9px 20px;
      border-radius: 6px; border: 2px solid #3d8b6e;
      text-decoration: none; font-weight: 600; font-size: 0.88em;
      display: inline-block; margin: 3px;
    }
    .btn-wf-outline:hover { background: #3d8b6e; color: white !important; }
    .btn-wf-outline.primary { background: #3d8b6e; color: white !important; }
    .btn-wf-outline.primary:hover { background: #2e6d56; }
    .workflow-steps {
      display: flex; align-items: stretch; padding: 30px 30px;
      background: #fff; justify-content: center; flex-wrap: wrap;
      border-bottom: 1px solid #e9ecef; gap: 0;
    }
    .wf-step {
      background: #fff; border: 1px solid #e9ecef; border-radius: 10px;
      padding: 22px 16px; flex: 1; min-width: 160px; max-width: 215px;
      text-align: center; border-top: 4px solid #3d8b6e;
      box-shadow: 0 2px 8px rgba(0,0,0,0.06);
    }
    .wf-num {
      display: inline-flex; align-items: center; justify-content: center;
      width: 40px; height: 40px; border-radius: 50%; background: #3d8b6e;
      color: white; font-size: 18px; font-weight: 700; margin-bottom: 10px;
    }
    .wf-step h4 { font-weight: 700; color: #2c3e50; font-size: 14px; margin: 0 0 7px; }
    .wf-step p  { font-size: 12px; color: #666; line-height: 1.55; margin: 0 0 10px; }
    .wf-btns { display: flex; flex-direction: column; gap: 6px; align-items: center; }
    .btn-wf-primary {
      background: #3d8b6e; color: white !important; padding: 7px 14px;
      border-radius: 6px; text-decoration: none; font-weight: 600; font-size: 11.5px;
      display: inline-block;
    }
    .btn-wf-primary:hover { background: #2e6d56; }
    .btn-wf-link {
      color: #3d8b6e !important; font-size: 11.5px; font-weight: 500;
      text-decoration: none;
    }
    .btn-wf-link:hover { text-decoration: underline; }
    .wf-connector {
      display: flex; align-items: center; padding: 0 10px; color: #adb5bd;
      font-size: 24px; font-weight: 300; align-self: center; flex-shrink: 0;
    }
    .tut-divider {
      text-align: center; padding: 24px 0 8px; font-size: 11.5px; font-weight: 700;
      letter-spacing: 0.12em; text-transform: uppercase; color: #3d8b6e; background: white;
      border-bottom: 1px solid #e9ecef;
    }
    .tut-section { padding: 38px 30px; border-bottom: 1px solid #e9ecef; background: white; }
    .tut-section:nth-child(even) { background: #f8f9fa; }
    .tut-section-inner {
      max-width: 1080px; margin: 0 auto; display: flex; align-items: flex-start; gap: 40px;
    }
    .tut-section-inner.reverse { flex-direction: row-reverse; }
    .tut-text { flex: 1; min-width: 0; }
    .tut-img  { flex: 1.25; min-width: 0; }
    .tut-step-badge {
      display: inline-flex; align-items: center; justify-content: center;
      width: 30px; height: 30px; border-radius: 50%; background: #2c3e50;
      color: white; font-size: 13px; font-weight: 700; margin-right: 10px; flex-shrink: 0;
    }
    .tut-section h2 {
      display: flex; align-items: center; font-size: 1.3em;
      font-weight: 700; color: #2c3e50; margin-bottom: 12px;
    }
    .tut-tab-badge {
      display: inline-block; background: #337ab7; color: white; font-size: 10.5px;
      font-weight: 600; padding: 2px 8px; border-radius: 12px;
      vertical-align: middle; margin-left: 8px; letter-spacing: 0.03em;
    }
    .tut-section p  { color: #444; line-height: 1.7; font-size: 13.5px; margin-bottom: 8px; }
    .tut-section ul { color: #444; line-height: 1.85; font-size: 13.5px; padding-left: 18px; margin-bottom: 8px; }
    .tut-section li { margin-bottom: 3px; }
    .tut-note {
      background: #eaf6f2; border-left: 4px solid #3d8b6e; border-radius: 0 6px 6px 0;
      padding: 9px 13px; margin: 12px 0; font-size: 13px; color: #2c5444;
    }
    .tut-warn {
      background: #fff8e1; border-left: 4px solid #f0ad4e; border-radius: 0 6px 6px 0;
      padding: 9px 13px; margin: 12px 0; font-size: 13px; color: #7a5c00;
    }
    .tut-screenshot {
      border-radius: 8px; box-shadow: 0 4px 20px rgba(0,0,0,0.12);
      border: 1px solid #dee2e6; max-width: 100%; height: auto; display: block;
    }
    .tut-screenshot-cap {
      font-size: 12px; color: #888; text-align: center; margin-top: 7px; font-style: italic;
    }
    .tut-screenshot-pair { display: flex; gap: 10px; margin-bottom: 10px; }
    .tut-screenshot-pair > div { flex: 1; }
    .tut-install { background: #2c3e50; color: #e0e0e0; padding: 40px 30px; }
    .tut-install > h2 { color: white; font-size: 1.4em; margin-bottom: 22px; text-align: center; }
    .tut-install-cards { display: flex; gap: 20px; max-width: 980px; margin: 0 auto; flex-wrap: wrap; }
    .tut-install-card {
      background: rgba(255,255,255,0.07); border-radius: 10px; padding: 22px; flex: 1;
      min-width: 270px; border: 1px solid rgba(255,255,255,0.1);
    }
    .tut-install-card h3 { color: #7ecba1; font-size: 1.05em; margin: 0 0 10px; }
    .tut-install-card p { font-size: 13px; line-height: 1.65; color: #bbb; margin-bottom: 6px; }
    .tut-code {
      background: #1a252f; color: #a8d8a8; padding: 10px 13px; border-radius: 6px;
      font-family: monospace; font-size: 12px; margin: 8px 0; white-space: pre;
      overflow-x: auto; display: block;
    }
    .tut-data { background: #fff; padding: 40px 30px; }
    .tut-data > h2 { text-align: center; color: #2c3e50; margin-bottom: 26px; font-weight: 700; }
    .tut-data-cards { display: flex; gap: 18px; max-width: 980px; margin: 0 auto; flex-wrap: wrap; }
    .tut-data-card {
      flex: 1; min-width: 240px; border: 1px solid #dee2e6; border-radius: 8px;
      padding: 20px; border-top: 3px solid #3d8b6e;
    }
    .tut-data-card h4 { color: #2c3e50; font-weight: 700; margin: 0 0 10px; font-size: 14px; }
    .tut-data-card p, .tut-data-card ul { font-size: 13px; color: #555; line-height: 1.7; margin: 0; padding-left: 16px; }
    .tut-footer {
      background: #f4f7f6; padding: 18px 20px; text-align: center;
      font-size: 13px; color: #888; border-top: 1px solid #e9ecef;
    }
    .tut-footer a { color: #3d8b6e; text-decoration: none; }
    @media (max-width: 768px) {
      .tut-section-inner, .tut-section-inner.reverse { flex-direction: column; }
      .tut-install-cards, .tut-data-cards, .workflow-steps { flex-direction: column; }
      .wf-connector { display: none; }
      .tut-screenshot-pair { flex-direction: column; }
      .wf-header h1 { font-size: 1.6em; }
    }

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

    # ── Page header ────────────────────────────────────────────────────────
    div(class = "wf-header",
      tags$h1("LongViewSC"),
      tags$p("Interactive visualisation of gene and isoform expression in long-read single-cell data"),
      div(class = "wf-header-btns",
        HTML('<a class="btn-wf-outline primary" href="https://longviewsc.researchsoftware.unimelb.edu.au/LongViewSC/" target="_blank">&#127758; Open App Online</a>
             <a class="btn-wf-outline" href="https://github.com/Sefi196/LongViewSC" target="_blank">&#128196; GitHub</a>
             <a class="btn-wf-outline" href="https://sefi196.github.io/FLAMESv2_LR_sc_tutorial/" target="_blank">&#128218; Data Prep Tutorial</a>
             <a class="btn-wf-outline" href="#installation">&#128295; Local Install</a>')
      )
    ),

    # ── 4-step workflow ────────────────────────────────────────────────────
    div(class = "workflow-steps",
      div(class = "wf-step",
        div(class = "wf-num", "1"),
        tags$h4("Prepare your data"),
        tags$p("Generate a compatible Seurat object using the FLAMES long-read single-cell pipeline."),
        div(class = "wf-btns",
          HTML('<a class="btn-wf-primary" href="https://sefi196.github.io/FLAMESv2_LR_sc_tutorial/" target="_blank">&#128218; FLAMESv2 Tutorial</a>'),
          HTML('<a class="btn-wf-link" href="#data-format">Data format details &#8595;</a>')
        )
      ),
      div(class = "wf-connector", HTML("&#8594;")),
      div(class = "wf-step",
        div(class = "wf-num", "2"),
        tags$h4("Configure"),
        tags$p("Upload your Seurat object, select gene, isoform assay and grouping variable, then pick isoforms from the chip panel.")
      ),
      div(class = "wf-connector", HTML("&#8594;")),
      div(class = "wf-step",
        div(class = "wf-num", "3"),
        tags$h4("Click GO"),
        tags$p("Press the green GO button to generate all visualisations across every tab at once.")
      ),
      div(class = "wf-connector", HTML("&#8594;")),
      div(class = "wf-step",
        div(class = "wf-num", "4"),
        tags$h4("Visualise"),
        tags$p("Explore 7 interactive tabs: UMAP, isoform stats, expression, transcript structure, heatmap, proportions, and trajectory."),
        div(class = "wf-btns",
          HTML('<a class="btn-wf-link" href="#visualisations">See examples &#8595;</a>')
        )
      )
    ),

    # ── Section label ─────────────────────────────────────────────────────
    div(id = "visualisations", class = "tut-divider", "Visualisation Tabs — Example Outputs (VIM, Demo Dataset)"),

    # ── Tab 1: Overview ───────────────────────────────────────────────────
    div(class = "tut-section",
      div(class = "tut-section-inner",
        div(class = "tut-text",
          tags$h2(span(class="tut-step-badge","1"), "Overview",
                  HTML('<span class="tut-tab-badge">Overview tab</span>')),
          tags$p("The Overview tab gives a high-level view of your dataset with a UMAP coloured by your chosen grouping variable, a gene-level feature plot, and a violin plot comparing expression across groups."),
          tags$ul(
            tags$li("Select a ", tags$strong("dimensional reduction"), " (e.g. UMAP, PCA) in the sidebar."),
            tags$li("Choose a ", tags$strong("metadata column"), " to colour cells by (cluster, cell type, condition)."),
            tags$li("Hit ", tags$strong("GO"), " to generate all plots — the same button updates every tab.")
          ),
          div(class = "tut-note",
            tags$strong("Quick start:"), " Click ", tags$strong("Demo Data"), " on the landing page to load the built-in dataset, then press GO."
          )
        ),
        div(class = "tut-img",
          div(class = "tut-screenshot-pair",
            div(
              tags$img(src = "example_images/01_umap_clusters.png", class = "tut-screenshot"),
              tags$p(class = "tut-screenshot-cap", "UMAP — Seurat clusters")
            ),
            div(
              tags$img(src = "example_images/02_vim_feature_plot.png", class = "tut-screenshot"),
              tags$p(class = "tut-screenshot-cap", "VIM gene expression on UMAP")
            )
          ),
          tags$img(src = "example_images/03_vim_violin.png", class = "tut-screenshot", style = "width:100%;"),
          tags$p(class = "tut-screenshot-cap", "Violin plot — VIM gene expression by cluster")
        )
      )
    ),

    # ── Tab 2: Isoform Statistics ─────────────────────────────────────────
    div(class = "tut-section",
      div(class = "tut-section-inner reverse",
        div(class = "tut-text",
          tags$h2(span(class="tut-step-badge","2"), "Isoform Statistics",
                  HTML('<span class="tut-tab-badge">Isoform Statistics tab</span>')),
          tags$p("A searchable, sortable table listing every isoform detected for your selected gene, with expression and detection statistics. Use this to identify the most relevant isoforms before diving into other plots."),
          tags$ul(
            tags$li("Isoforms are named ", tags$code("ENST00000XXXXX.X-GENE"), " (e.g. ", tags$code("ENST00000224237.9-VIM"), ")."),
            tags$li("A tick (", HTML("&#10003;"), ") marks isoforms currently selected in the sidebar chip panel."),
            tags$li("Clicking a row also selects/deselects that isoform.")
          ),
          div(class = "tut-note",
            tags$strong("Chip selector:"), " Coloured chips in the sidebar toggle which isoforms appear across the Expression, Heatmap, Proportions and Trajectory tabs. Click any chip to select or deselect it."
          )
        ),
        div(class = "tut-img",
          tags$img(src = "example_images/10_chip_selector.png", class = "tut-screenshot",
                   style = "width:100%; margin-bottom: 14px;"),
          tags$p(class = "tut-screenshot-cap", "Isoform chip selector — click chips to toggle isoforms on/off"),
          tags$img(src = "example_images/09_isoform_table.png", class = "tut-screenshot",
                   style = "width:100%;"),
          tags$p(class = "tut-screenshot-cap", "Isoform Statistics table — VIM isoforms with expression metrics")
        )
      )
    ),

    # ── Tab 3: Isoform Expression ─────────────────────────────────────────
    div(class = "tut-section",
      div(class = "tut-section-inner",
        div(class = "tut-text",
          tags$h2(span(class="tut-step-badge","3"), "Isoform Expression",
                  HTML('<span class="tut-tab-badge">Isoform Expression tab</span>')),
          tags$p("Single-cell resolution plots for each selected isoform: a feature plot projected onto the isoform-level UMAP embedding and a violin plot across groups."),
          tags$ul(
            tags$li(tags$strong("Feature plots"), " use the isoform-level UMAP (", tags$code("umap_iso"), ") if present, otherwise the gene-level UMAP."),
            tags$li(tags$strong("Violin plots"), " show expression per isoform, grouped by your chosen metadata column."),
            tags$li("Only isoforms selected via the chip panel are shown.")
          ),
          div(class = "tut-warn",
            tags$strong("Important:"), " Set the ", tags$strong("Isoform Assay"), " field in the sidebar to the correct assay name (e.g. ", tags$code("iso"), ") before clicking GO."
          )
        ),
        div(class = "tut-img",
          tags$img(src = "example_images/04_vim_isoform_feature.png", class = "tut-screenshot"),
          tags$p(class = "tut-screenshot-cap", "VIM isoforms ENST00000224237 and ENST00000487938 on isoform UMAP")
        )
      )
    ),

    # ── Tab 4: Transcript Structure ───────────────────────────────────────
    div(class = "tut-section",
      div(class = "tut-section-inner reverse",
        div(class = "tut-text",
          tags$h2(span(class="tut-step-badge","4"), "Transcript Structure",
                  HTML('<span class="tut-tab-badge">Transcript Structure tab</span>')),
          tags$p("Visualises the exon-intron architecture of each selected isoform from your uploaded GTF file, enabling direct structural comparison."),
          tags$ul(
            tags$li("Exons are drawn as filled blocks; introns as directional arrows showing strand."),
            tags$li("Upload a ", tags$strong(".gtf"), " file in the sidebar — if omitted this tab is empty."),
            tags$li("GTF transcript IDs must match the ENST IDs in your isoform assay.")
          ),
          div(class = "tut-note",
            tags$strong("Naming:"), " Isoforms must follow the format ", tags$code("ENST00000XXXXX.X-GENENAME"),
            " for automatic GTF matching."
          )
        ),
        div(class = "tut-img",
          tags$img(src = "example_images/05_vim_transcript_structure.png", class = "tut-screenshot"),
          tags$p(class = "tut-screenshot-cap", "Exon-intron structure of all four VIM isoforms")
        )
      )
    ),

    # ── Tab 5: Pseudobulk Heatmap ─────────────────────────────────────────
    div(class = "tut-section",
      div(class = "tut-section-inner",
        div(class = "tut-text",
          tags$h2(span(class="tut-step-badge","5"), "Pseudobulk Heatmap",
                  HTML('<span class="tut-tab-badge">Pseudobulk Heatmap tab</span>')),
          tags$p("Aggregates isoform counts per group, normalises to log(CPM+1), and renders an interactive clustered heatmap — ideal for comparing relative isoform usage across conditions or cell types at a glance."),
          tags$ul(
            tags$li("Rows are isoforms; columns are groups; colour intensity = log(CPM+1)."),
            tags$li("Rows are clustered by expression similarity."),
            tags$li("The heatmap is fully interactive — hover for exact values, zoom, pan, and download via the camera icon (", HTML("&#128247;"), ").")
          )
        ),
        div(class = "tut-img",
          tags$img(src = "example_images/06_vim_pseudobulk_heatmap.png", class = "tut-screenshot"),
          tags$p(class = "tut-screenshot-cap", "VIM isoform pseudobulk expression (log CPM+1) across 8 clusters")
        )
      )
    ),

    # ── Tab 6: Proportions ────────────────────────────────────────────────
    div(class = "tut-section",
      div(class = "tut-section-inner reverse",
        div(class = "tut-text",
          tags$h2(span(class="tut-step-badge","6"), "Isoform Proportions",
                  HTML('<span class="tut-tab-badge">Proportions tab</span>')),
          tags$p("Faceted pie charts — one per group — showing each selected isoform as a proportion of total gene expression. Ideal for detecting isoform switching between cell types or conditions."),
          tags$ul(
            tags$li("Selected isoforms are shown as distinct colours; all others are grouped as ", tags$strong("Other isoforms"), " (grey)."),
            tags$li("Percentages are displayed for slices > 5%."),
            tags$li("Set ", tags$strong("Min counts"), " in the sidebar to hide groups with low coverage.")
          ),
          div(class = "tut-note",
            tags$strong("Isoform switching:"), " Compare the dominant (largest) slice across groups — a shift in colour indicates a switch in the preferred isoform."
          )
        ),
        div(class = "tut-img",
          tags$img(src = "example_images/07_vim_isoform_pie.png", class = "tut-screenshot"),
          tags$p(class = "tut-screenshot-cap", "VIM isoform proportions across all 8 clusters")
        )
      )
    ),

    # ── Tab 7: Trajectory ─────────────────────────────────────────────────
    div(class = "tut-section",
      div(class = "tut-section-inner",
        div(class = "tut-text",
          tags$h2(span(class="tut-step-badge","7"), "Expression Trajectory",
                  HTML('<span class="tut-tab-badge">Trajectory tab</span>')),
          tags$p("Plots single-cell normalised isoform expression across an ordered sequence of groups (e.g. developmental stages), with a loess smooth overlaid per isoform panel. Reveals gradual changes in isoform usage along a trajectory."),
          tags$ul(
            tags$li(tags$strong("Drag and drop"), " the group labels in the sidebar to define the left-to-right ordering."),
            tags$li("Each isoform gets its own panel with an independent y-axis."),
            tags$li("The black loess curve highlights overall expression trends.")
          ),
          div(class = "tut-note",
            tags$strong("Pairing tip:"), " Use Proportions to see ", tags$em("which"), " isoform dominates, and Trajectory to see ", tags$em("when"), " it changes along the axis."
          )
        ),
        div(class = "tut-img",
          tags$img(src = "example_images/08_vim_trajectory.png", class = "tut-screenshot"),
          tags$p(class = "tut-screenshot-cap", "VIM isoform expression trajectory across clusters 0 \u2192 7")
        )
      )
    ),

    # ── Data Format ───────────────────────────────────────────────────────
    div(id = "data-format", class = "tut-data",
      tags$h2("Input Data Format"),
      div(class = "tut-data-cards",
        div(class = "tut-data-card",
          tags$h4("Seurat Object (.rds / .qs / .qs2)"),
          tags$ul(
            tags$li("Must contain a gene assay (e.g. ", tags$code("RNA"), ") and an isoform assay (e.g. ", tags$code("iso"), ")."),
            tags$li("Isoform features: ", tags$code("ENST00000XXXXX.X-GENENAME"), " (e.g. ", tags$code("ENST00000544301.7-VIM"), ")."),
            tags$li("Requires at least one dimensional reduction (UMAP, PCA)."),
            tags$li("Metadata must include a column for cell grouping.")
          ),
          div(class = "tut-note",
            HTML('See the <a href="https://sefi196.github.io/FLAMESv2_LR_sc_tutorial/" target="_blank">FLAMESv2 tutorial</a> for generating compatible objects.')
          )
        ),
        div(class = "tut-data-card",
          tags$h4("GTF Annotation (.gtf) — Optional"),
          tags$ul(
            tags$li("Standard GTF format with a ", tags$code("transcript_id"), " attribute."),
            tags$li("Transcript IDs must match the ENST IDs in your isoform assay."),
            tags$li("Only required for the ", tags$strong("Transcript Structure"), " tab."),
            tags$li("Files up to 5 GB supported locally.")
          )
        ),
        div(class = "tut-data-card",
          tags$h4("Performance Tips"),
          tags$ul(
            tags$li(tags$strong("Online:"), " Keep uploads under 200 MB for best responsiveness."),
            tags$li(tags$strong("Local:"), " Use ", tags$code("qs::qsave()"), " format for large objects — up to 10x faster loading than .rds."),
            tags$li(tags$strong("qs / qs2:"), " Both formats are supported; the app auto-detects which library to use.")
          )
        )
      )
    ),

    # ── Installation ──────────────────────────────────────────────────────
    div(id = "installation", class = "tut-install",
      tags$h2("Installation"),
      div(class = "tut-install-cards",
        div(class = "tut-install-card",
          tags$h3("Option 1 \u2014 Use Online"),
          tags$p("Access LongViewSC directly in your browser. No setup required."),
          HTML('<div class="tut-code">https://longviewsc.researchsoftware.unimelb.edu.au/LongViewSC/</div>'),
          tags$p(style="font-size:12px;color:#999;margin-top:6px;", "Recommended for files < 200 MB.")
        ),
        div(class = "tut-install-card",
          tags$h3("Option 2 \u2014 Local Installation"),
          tags$p("For large files (> 200 MB) or offline use. Requires conda or miniconda."),
          HTML('<div class="tut-code">git clone https://github.com/Sefi196/LongViewSC.git
cd LongViewSC
conda env create -f environment.yml
conda activate LongViewSC_env</div>'),
          tags$p(style="font-size:12px;color:#aaa;", "Install ggtranscript (not on conda):"),
          HTML('<div class="tut-code">Rscript -e "devtools::install_github(\'dzhang32/ggtranscript\')"</div>'),
          tags$p(style="font-size:12px;color:#aaa;", "Launch the app:"),
          HTML('<div class="tut-code">Rscript -e "shiny::runApp(\'app.R\', launch.browser=TRUE)"</div>'),
          tags$p(style="font-size:12px;color:#999;margin-top:6px;", "Runs entirely locally. No data is sent online.")
        )
      )
    ),

    # ── Footer ────────────────────────────────────────────────────────────
    div(class = "tut-footer",
      HTML('Developed by <strong>Sefi Prawer</strong> \u2014 Clark Laboratory, University of Melbourne. &nbsp;
           <a href="mailto:sefi.prawer@unimelb.edu.au">sefi.prawer@unimelb.edu.au</a> &nbsp;|&nbsp;
           <a href="https://github.com/Sefi196/LongViewSC" target="_blank">GitHub</a> &nbsp;|&nbsp;
           If you use LongViewSC please cite our work.')
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
        setTimeout(function() { $(window).trigger("resize"); }, 310);
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
              uiOutput("isoform_chips_panel"),
              
          )  # end of analysisForm
        ),
        
        mainPanel(
          id = "main_panel_wrap",
          tags$button(id = "sidebarToggle", title = "Toggle sidebar",
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