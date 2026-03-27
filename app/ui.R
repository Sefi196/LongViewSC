ui <- fluidPage(
  useShinyjs(),
  
  # Add custom styles for a header banner and plot titles
  tags$head(
    tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
    tags$link(rel = "preconnect", href = "https://fonts.gstatic.com", crossorigin = NA),
    tags$link(rel = "stylesheet",
              href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap"),
    tags$style(HTML("
    /* ── Global font ─────────────────────────────────────────────────── */
    body, input, select, textarea, button, .selectize-input, .shiny-input-container label {
      font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif !important;
    }
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
      position: relative;
    }
    .sidebar-section .form-group { margin-bottom: 6px; }
    .sidebar-section-label {
      font-size: 12.5px;
      font-weight: 700;
      text-transform: uppercase;
      letter-spacing: 0.06em;
      color: #3d8b6e;
      margin: 0 0 7px 0;
    }
    .sidebar-section-label .fa, .sidebar-section-label .glyphicon { margin-right: 4px; }
    .chip-panel-title {
      font-size: 11px;
      font-weight: 700;
      color: #3d8b6e;
      text-transform: uppercase;
      letter-spacing: 0.06em;
      margin-bottom: 8px;
    }
    .input-label-row {
      display: flex;
      align-items: center;
      gap: 6px;
      margin-bottom: 2px;
    }
    .input-label-row label {
      margin: 0;
      font-size: 13px;
      font-weight: 500;
      color: #444;
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
      font-style: italic;
      font-family: serif;
      cursor: help;
      line-height: 1;
      opacity: 0.7;
      flex-shrink: 0;
      transition: opacity 0.15s ease, background 0.15s ease;
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
        background: linear-gradient(180deg, #ffffff 0%, #f3f6fa 100%);
        padding: 15px;
        border-radius: 0;
        position: relative;
        border: none;
        border-right: 1.5px solid rgba(0,0,0,0.1);
        min-height: 100vh;
        box-shadow: none;
      }
      .plot-box {
        border: 1.5px solid rgba(100,100,100,0.3);
        padding: 14px;
        border-radius: 12px;
        box-shadow: 0 2px 12px rgba(0,0,0,0.07);
        background-color: white;
        margin-bottom: 20px;
        min-width: 220px;
        min-height: 120px;
        box-sizing: border-box;
        overflow: visible;
        position: relative;
      }
      /* Resize grip — three stacked lines in bottom-right */
      .plot-box .ui-resizable-se {
        width: 16px; height: 16px;
        right: 2px; bottom: 2px;
        cursor: se-resize;
        background: none;
        display: flex; align-items: flex-end; justify-content: flex-end;
      }
      .plot-box .ui-resizable-se::after {
        content: \"\";
        display: block;
        width: 10px; height: 10px;
        border-right: 2px solid rgba(120,120,120,0.5);
        border-bottom: 2px solid rgba(120,120,120,0.5);
        border-radius: 0 0 3px 0;
      }
      .plot-box:hover .ui-resizable-se::after {
        border-color: rgba(80,80,80,0.85);
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
      background: #2c3e50;
      color: #fff;
      border: none;
      border-radius: 8px;
      width: 42px;
      height: 42px;
      font-size: 18px;
      cursor: pointer;
      box-shadow: 0 2px 10px rgba(0,0,0,0.25);
      display: flex;
      align-items: center;
      justify-content: center;
      padding: 0;
      transition: background 0.15s ease;
    }
    #sidebarToggle:hover { background: #337ab7; }
    /* Collapsible plot boxes */
    .plot-title-row .collapse-chevron {
      margin-left: auto;
      font-size: 20px;
      color: #bbb;
      transition: transform 0.2s ease, color 0.15s ease;
      flex-shrink: 0;
      line-height: 1;
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
      text-align: center; border-top: 4px solid #337ab7;
      box-shadow: 0 2px 8px rgba(0,0,0,0.06);
    }
    .wf-num {
      display: inline-flex; align-items: center; justify-content: center;
      width: 40px; height: 40px; border-radius: 50%; background: #337ab7;
      color: white; font-size: 18px; font-weight: 700; margin-bottom: 10px;
    }
    .wf-step h4 { font-weight: 700; color: #2c3e50; font-size: 14px; margin: 0 0 7px; }
    .wf-step p  { font-size: 12px; color: #666; line-height: 1.55; margin: 0 0 10px; }
    .wf-btns { display: flex; flex-direction: column; gap: 6px; align-items: center; }
    .btn-wf-primary {
      background: #337ab7; color: white !important; padding: 7px 14px;
      border-radius: 6px; text-decoration: none; font-weight: 600; font-size: 11.5px;
      display: inline-block;
    }
    .btn-wf-primary:hover { background: #2a6496; }
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
    .tut-section p  { color: #444; line-height: 1.75; font-size: 14.5px; margin-bottom: 10px; }
    .tut-section ul { color: #444; line-height: 1.9; font-size: 14.5px; padding-left: 20px; margin-bottom: 10px; }
    .tut-section li { margin-bottom: 5px; }
    .tut-note {
      background: #eaf6f2; border-left: 4px solid #3d8b6e; border-radius: 0 6px 6px 0;
      padding: 10px 14px; margin: 14px 0; font-size: 13.5px; color: #2c5444; line-height: 1.65;
    }
    .tut-warn {
      background: #fff8e1; border-left: 4px solid #f0ad4e; border-radius: 0 6px 6px 0;
      padding: 10px 14px; margin: 14px 0; font-size: 13.5px; color: #7a5c00; line-height: 1.65;
    }
    .tut-screenshot {
      border-radius: 8px; box-shadow: 0 4px 20px rgba(0,0,0,0.12);
      border: 1px solid #dee2e6; max-width: 100%; height: auto; display: block;
    }
    .tut-screenshot-cap {
      font-size: 12.5px; color: #666; text-align: center; margin-top: 8px; font-style: italic; line-height: 1.4;
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
    /* ── Inline code styling — file formats etc ─────────────────────── */
    code {
      background: #e8f4fd;
      color: #1a6ea8;
      padding: 1px 5px;
      border-radius: 3px;
      font-size: 0.88em;
      font-family: monospace;
    }
    @media (max-width: 768px) {
      .tut-section-inner, .tut-section-inner.reverse { flex-direction: column; }
      .tut-install-cards, .tut-data-cards, .workflow-steps { flex-direction: column; }
      .wf-connector { display: none; }
      .tut-screenshot-pair { flex-direction: column; }
      .wf-header h1 { font-size: 1.6em; }
    }

    /* ── CSS spinner (no external URL) ─────────────────────────────────── */
    @keyframes lvsc-spin { to { transform: rotate(360deg); } }
    .lvsc-spinner-ring {
      width: 44px; height: 44px;
      border: 5px solid #e9ecef;
      border-top-color: #337ab7;
      border-radius: 50%;
      animation: lvsc-spin 0.8s linear infinite;
    }

    /* ── NProgress-style top loading bar ───────────────────────────────── */
    #lvsc-progress-bar {
      position: fixed; top: 0; left: 0; height: 3px;
      background: linear-gradient(90deg, #337ab7, #3d8b6e);
      z-index: 99999; width: 0%; pointer-events: none;
      transition: width 0.3s ease, opacity 0.4s ease;
      opacity: 0;
    }

    /* ── Responsive: stack Trajectory columns on medium screens ────────── */
    @media (max-width: 992px) {
      .traj-settings-col, .traj-main-col { width: 100% !important; float: none !important; }
    }

    /* ── Proportions settings gear panel ───────────────────────────────── */
    #pie-settings-panel {
      background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 6px;
      padding: 10px 14px; margin-bottom: 8px;
    }
    #pie-settings-toggle {
      font-size: 12px; padding: 3px 10px; color: #555;
      border: 1px solid #ccc; border-radius: 5px; background: #fff;
      cursor: pointer; display: inline-flex; align-items: center; gap: 5px;
    }
    #pie-settings-toggle:hover { background: #e9ecef; color: #337ab7; border-color: #337ab7; }
    #coexpr-settings-panel, #traj-settings-panel {
      background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 6px;
      padding: 8px 12px; margin-bottom: 8px;
    }
    .settings-gear-btn {
      font-size: 12px; padding: 3px 10px; color: #555;
      border: 1px solid #ccc; border-radius: 5px; background: #fff;
      cursor: pointer; display: inline-flex; align-items: center; gap: 5px;
    }
    .settings-gear-btn:hover { background: #e9ecef; color: #337ab7; border-color: #337ab7; }

    /* ── Input helper text ──────────────────────────────────────────────── */
    .input-helper-text {
      font-size: 11px; color: #999; margin-top: -4px; margin-bottom: 6px;
      line-height: 1.3;
    }

    /* ── Toast notification ─────────────────────────────────────────────── */
    #lvsc-toast {
      position: fixed; bottom: 28px; right: 28px;
      background: #2c3e50; color: #fff;
      padding: 11px 18px; border-radius: 8px;
      font-size: 13px; font-weight: 500;
      box-shadow: 0 4px 18px rgba(0,0,0,0.22);
      z-index: 99999; opacity: 0; pointer-events: none;
      transition: opacity 0.28s ease, transform 0.28s ease;
      transform: translateY(8px);
      display: flex; align-items: center; gap: 9px;
    }
    #lvsc-toast.lvsc-toast-show {
      opacity: 1; transform: translateY(0);
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

          # Clark Lab Data Dropdown
          div(class = "dropdown",
              style = "margin-left: 5px;",
              tags$button(
                class = "btn btn-warning btn-lg dropdown-toggle",
                type = "button",
                `data-toggle` = "dropdown",
                HTML('<i class="fa fa-database"></i> Clark Lab Data')
              ),
              tags$ul(
                class = "dropdown-menu",
                if (length(clark_datasets) == 0) {
                  tags$li(tags$span(style = "padding: 8px 16px; color: grey; display: block;",
                                    "No datasets available"))
                } else {
                  lapply(names(clark_datasets), function(nm) {
                    tags$li(
                      tags$a(href = "#",
                             onclick = sprintf(
                               "document.getElementById('clark_pw_dataset_name').innerText = '%s'; $('#clark_pw_modal').modal('show'); return false;",
                               nm),
                             nm)
                    )
                  })
                }
              )
          ),

          # ── Clark Lab password modal ───────────────────────────────────────
          tags$div(id = "clark_pw_modal", class = "modal fade", tabindex = "-1",
                   role = "dialog", `aria-labelledby` = "clark_pw_modal_label",
            tags$div(class = "modal-dialog modal-sm", role = "document",
              tags$div(class = "modal-content",
                tags$div(class = "modal-header",
                  tags$button(type = "button", class = "close", `data-dismiss` = "modal",
                              tags$span(`aria-hidden` = "true", HTML("&times;"))),
                  tags$h4(class = "modal-title", id = "clark_pw_modal_label",
                          HTML('<i class="fa fa-lock" style="margin-right:6px;"></i>Clark Lab Access'))
                ),
                tags$div(class = "modal-body",
                  tags$p(style = "font-size:13px; color:#555; margin-bottom:10px;",
                         "Enter the password to load ",
                         tags$strong(id = "clark_pw_dataset_name", "this dataset"), "."),
                  tags$div(class = "form-group", style = "margin-bottom:8px;",
                    tags$input(id = "clark_pw_input", type = "password",
                               class = "form-control", placeholder = "Password",
                               onkeydown = "if(event.key==='Enter') document.getElementById('clark_pw_confirm').click();")
                  ),
                  tags$div(id = "clark_pw_error",
                           style = "color:#c0392b; font-size:12px; display:none; margin-top:4px;",
                           HTML('<i class="fa fa-exclamation-circle"></i> Incorrect password.'))
                ),
                tags$div(class = "modal-footer",
                  tags$button(type = "button", class = "btn btn-default", `data-dismiss` = "modal", "Cancel"),
                  tags$button(id = "clark_pw_confirm", type = "button",
                              class = "btn btn-warning",
                              onclick = paste0(
                                "var pw = document.getElementById('clark_pw_input').value;",
                                "var nm = document.getElementById('clark_pw_dataset_name').innerText;",
                                "Shiny.setInputValue('clark_pw_attempt', {pw: pw, dataset: nm, nonce: Math.random()});"),
                              HTML('<i class="fa fa-unlock" style="margin-right:4px;"></i>Load Dataset'))
                )
              )
            )
          ),

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
          actionButton("backToLanding", label = "", icon = icon("arrow-left"), class = "btn btn-dark btn-lg",
                       `aria-label` = "Return to home screen"),
          
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
  # Top progress bar (shown on GO press, fades on idle)
  tags$div(id = "lvsc-progress-bar"),

  div(id = "spinner",
      tags$div(class = "lvsc-spinner-ring"),
      tags$p(style = "margin: 8px 0 0; font-size: 12px; color: #666;", "Loading…"),
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
      tags$div(style = "display:flex; justify-content:center; margin-bottom:6px;",
        tags$div(style = "display:inline-flex; align-items:stretch; padding: 8px 0;",
          # Vertical bar of the L — green spine
          tags$div(style = "width:12px; background:#3d8b6e; border-radius:3px; flex-shrink:0;"),
          # Right column: title text above, exon track below (= bottom of the L)
          tags$div(style = "display:flex; flex-direction:column;",
            tags$h1(style = "margin:0 0 4px 0; padding-left:10px; font-size:2em; font-weight:700; color:#2c3e50; line-height:1; white-space:nowrap;",
                    "ongViewSC"),
            HTML('<div style="position:relative; display:flex; align-items:center; height:11px; width:100%; padding-left:10px; box-sizing:border-box;">
              <div style="position:absolute; left:0; right:0; height:2px; background:#3d8b6e; top:50%; transform:translateY(-50%);"></div>
              <div style="position:relative; height:11px; flex:3; background:#337ab7; border-radius:2px; margin-right:5px;"></div>
              <div style="position:relative; height:11px; flex:2.2; background:#e67e22; border-radius:2px; margin-right:5px;"></div>
              <div style="position:relative; height:11px; flex:1.8; background:#8e44ad; border-radius:2px;"></div>
            </div>')
          )
        )
      ),
      tags$p(style = "color:#666; font-size:1em; margin:0 0 6px;",
             "Interactive visualisation of gene and isoform expression in long-read single-cell data"),
      tags$p(style = "color:#888; font-size:0.85em; margin:0;",
             HTML('Developed by the <a href="https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/clark-lab" target="_blank" style="color:#3d8b6e;">Clark Laboratory</a>, University of Melbourne \u2014
                  <a href="https://github.com/Sefi196/LongViewSC" target="_blank" style="color:#3d8b6e;">GitHub</a> \u2014
                  <a href="https://sefi196.github.io/FLAMESv2_LR_sc_tutorial/" target="_blank" style="color:#3d8b6e;">FLAMESv2 Tutorial</a>'))
    ),

    # ── 4-step workflow ────────────────────────────────────────────────────
    div(class = "workflow-steps",

      div(class = "wf-step", style = "border-top-color: #3d8b6e;",
        div(class = "wf-num", style = "background: #3d8b6e;", "1"),
        tags$h4("Load your data"),
        tags$p("Upload a Seurat object (", tags$code(".rds"), ", ", tags$code(".qs"), ", or ", tags$code(".qs2"), ") with gene and isoform expression data, or try the built-in demo dataset."),
        div(class = "wf-btns",
          HTML('<a class="btn-wf-primary" style="background:#3d8b6e;" href="https://sefi196.github.io/FLAMESv2_LR_sc_tutorial/" target="_blank">&#128218; FLAMESv2 Tutorial</a>'),
          HTML('<a class="btn-wf-link" href="#data-format">Data format &#8595;</a>')
        )
      ),

      div(class = "wf-connector", HTML("&#8594;")),

      div(class = "wf-step", style = "border-top-color: #337ab7;",
        div(class = "wf-num", style = "background: #337ab7;", "2"),
        tags$h4("Configure"),
        tags$p("Type a gene name in the sidebar. Select your isoform assay, dimensional reduction, and metadata grouping column. Optionally load a GTF for transcript structure.")
      ),

      div(class = "wf-connector", HTML("&#8594;")),

      div(class = "wf-step", style = "border-top-color: #e67e22;",
        div(class = "wf-num", style = "background: #e67e22;", "3"),
        tags$h4("Run Analysis"),
        tags$p("Press Run Analysis to generate all plots. The top 4 isoforms are pre-selected — use the coloured chips in the sidebar to toggle isoforms across every tab.")
      ),

      div(class = "wf-connector", HTML("&#8594;")),

      div(class = "wf-step", style = "border-top-color: #8e44ad;",
        div(class = "wf-num", style = "background: #8e44ad;", "4"),
        tags$h4("Explore & export"),
        tags$p("Navigate the tabs to explore cell-type maps, isoform statistics, expression plots, transcript structures, heatmaps, and proportions. Download individual plots or generate a full HTML/PDF report."),
        div(class = "wf-btns",
          HTML('<a class="btn-wf-link" href="#visualisations">See examples &#8595;</a>')
        )
      )
    ),

    # ── Section label ─────────────────────────────────────────────────────
    div(id = "visualisations", class = "tut-divider", "Visualisation Tabs — Example Outputs from the Demo Dataset"),

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
          tags$img(src = "example_images/01_overview_stack.png?v=202603270900", class = "tut-screenshot",
            style = "width:100%;"),
          tags$p(class = "tut-screenshot-cap", "Cell type UMAP (top) and VIM gene expression feature plot (bottom)")
        )
      )
    ),

    # ── Tab 2: Isoform Statistics ─────────────────────────────────────────
    div(class = "tut-section",
      div(class = "tut-section-inner reverse",
        div(class = "tut-text",
          tags$h2(span(class="tut-step-badge","2"), "Isoform Statistics",
                  HTML('<span class="tut-tab-badge">Isoform Summary tab</span>')),
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
          tags$img(src = "example_images/10_chip_selector.png?v=202603270442", class = "tut-screenshot",
                   style = "width:100%; margin-bottom: 14px;"),
          tags$p(class = "tut-screenshot-cap", "Isoform chip selector — click chips to toggle isoforms on/off"),
          tags$img(src = "example_images/09_isoform_table.png?v=202603270442", class = "tut-screenshot",
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
            tags$li("Only isoforms selected via the chip panel are shown.")
          ),
          div(class = "tut-warn",
            tags$strong("Important:"), " Set the ", tags$strong("Isoform Assay"), " field in the sidebar to the correct assay name (e.g. ", tags$code("iso"), ") before clicking GO."
          )
        ),
        div(class = "tut-img",
          tags$img(src = "example_images/04_vim_isoform_feature.png?v=202603270522", class = "tut-screenshot", style = "max-height:270px; width:auto; display:block; margin:0 auto;"),
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
            tags$strong("Note:"), " Any standard GTF file should work. The app matches transcript IDs automatically."
          )
        ),
        div(class = "tut-img",
          tags$img(src = "example_images/05_vim_transcript_structure.png?v=202603270442", class = "tut-screenshot"),
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
          tags$img(src = "example_images/06_vim_pseudobulk_heatmap.png?v=202603270522", class = "tut-screenshot"),
          tags$p(class = "tut-screenshot-cap", "VIM isoform pseudobulk expression (log CPM+1) across clusters")
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
          tags$img(src = "example_images/07_isoform_pie_rbfox1.png?v=202603270900", class = "tut-screenshot", style = "width:100%;"),
          tags$p(class = "tut-screenshot-cap", "RBFOX1 isoform proportions in Glutamatergic neurons — each slice shows the fraction of total gene expression contributed by each isoform")
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
            tags$li("A loess smooth with a shaded confidence ribbon highlights overall expression trends per isoform.")
          ),
          div(class = "tut-note",
            tags$strong("Pairing tip:"), " Use Proportions to see ", tags$em("which"), " isoform dominates, and Trajectory to see ", tags$em("when"), " it changes along the axis."
          )
        ),
        div(class = "tut-img",
          tags$img(src = "example_images/08_trajectory_wide.png?v=202603270900", class = "tut-screenshot", style = "width:100%;"),
          tags$p(class = "tut-screenshot-cap", "VIM isoform expression trajectory across cell types — loess smooth with 95% confidence ribbon per isoform")
        )
      )
    ),

    # ── Tab 8: Compare Two Isoforms ───────────────────────────────────────
    div(class = "tut-section",
      div(class = "tut-section-inner reverse",
        div(class = "tut-text",
          tags$h2(span(class="tut-step-badge","8"), "Compare Two Isoforms",
                  HTML('<span class="tut-tab-badge">Compare Two Isoforms tab</span>')),
          tags$p("Directly compare the expression of any two isoforms (or genes) at single-cell resolution across three complementary views:"),
          tags$ul(
            tags$li(tags$strong("Blend:"), " a Seurat-style UMAP blend — one isoform in red, the other in blue, co-expressing cells in purple. Works across genes (e.g. TBR1 vs MAPT) or within a gene. Toggle between ", tags$strong("Isoforms"), " and ", tags$strong("Genes"), " using the Compare selector."),
            tags$li(tags$strong("Scatter:"), " each cell is a dot — X axis = isoform 1 expression, Y axis = isoform 2 expression. Includes a linear fit, R², and optional marginal distributions. Switch between isoform or gene-level comparison."),
            tags$li(tags$strong("Feature Map:"), " a categorical UMAP classifying every cell as ",
              tags$strong("Both"), " (co-expressed), ",
              tags$strong("Only X"), ", ",
              tags$strong("Only Y"), ", or ",
              tags$strong("Neither"), " — ideal for mutual-exclusivity analysis.")
          ),
          div(class = "tut-note",
            tags$strong("Tip:"), " Use the ", tags$strong("Show cells"), " filter to restrict comparisons to specific cell groups. All views support the full download modal (PNG, PDF, JPEG)."
          )
        ),
        div(class = "tut-img",
          tags$img(src = "example_images/12_isoform_blend_tbr1_mapt.png?v=202603270900", class = "tut-screenshot",
                   style = "width:100%; margin-bottom:10px;"),
          tags$p(class = "tut-screenshot-cap", "Blend view — dominant isoforms of TBR1 (red) vs MAPT (blue) on the UMAP; purple cells co-express both"),
          tags$img(src = "example_images/13_isoform_coexpr_categorical.png?v=202603270800", class = "tut-screenshot",
                   style = "width:100%; margin-top:14px;"),
          tags$p(class = "tut-screenshot-cap", "Feature Map view — categorical co-expression: cells classified by which isoform(s) they express")
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
            tags$li("Isoform IDs must follow ", tags$code("ENST00000XXXXX.X-GENENAME"), " — the gene suffix must match gene assay names (e.g. ", tags$code("ENST00000544301.7-VIM"), " ↔ ", tags$code("VIM"), ")."),
            tags$li("Requires at least one dimensional reduction (UMAP, PCA)."),
            tags$li("Metadata must include a column for cell grouping.")
          ),
          div(class = "tut-note",
            HTML('See the <a class="btn-wf-primary" style="background:#3d8b6e; font-size:12px; padding:5px 12px;" href="https://sefi196.github.io/FLAMESv2_LR_sc_tutorial/" target="_blank">&#128218; FLAMESv2 Tutorial</a> for generating compatible objects.')
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

# Toast element (always present, shown via JS)
  tags$div(id = "lvsc-toast"),

# Main UI - Initially Hidden
  div(id = "mainUI",

      tags$script(HTML('
    $(document).ready(function() {
      // ── Bootstrap popovers for chip info icons ───────────────────────────
      $(document).on("mouseenter focus", "[data-toggle=popover]", function() {
        $(this).popover("show");
      }).on("mouseleave blur", "[data-toggle=popover]", function() {
        $(this).popover("hide");
      });

      // ── Collapsible plot boxes ────────────────────────────────────────────
      function addChevrons() {
        $(".plot-box .plot-title-row").each(function() {
          if ($(this).find(".collapse-chevron").length === 0) {
            $(this).append($("<span>").addClass("collapse-chevron").html("&#9662;"));
          }
        });
      }
      $(document).on("shiny:idle", addChevrons);
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

      // ── Sidebar collapse ──────────────────────────────────────────────────
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

      // ── NProgress-style top loading bar ──────────────────────────────────
      var $bar = $("#lvsc-progress-bar");
      var barTimer = null;
      $(document).on("click", "#GO", function() {
        clearTimeout(barTimer);
        $bar.stop(true).css({ opacity: 1, width: "20%" })
            .animate({ width: "75%" }, 1500);
      });
      $(document).on("shiny:idle", function() {
        $bar.stop(true).animate({ width: "100%" }, 250, function() {
          barTimer = setTimeout(function() {
            $bar.animate({ opacity: 0 }, 300, function() {
              $bar.css({ width: "0%", opacity: 1 });
            });
          }, 200);
        });
      });

      // ── Clark Lab password modal handlers ────────────────────────────────
      Shiny.addCustomMessageHandler("clark_pw_wrong", function(x) {
        document.getElementById("clark_pw_error").style.display = "block";
        document.getElementById("clark_pw_input").focus();
      });
      Shiny.addCustomMessageHandler("clark_pw_ok", function(x) {
        document.getElementById("clark_pw_error").style.display = "none";
        document.getElementById("clark_pw_input").value = "";
        $("#clark_pw_modal").modal("hide");
      });

      // ── Enter key fires GO button ─────────────────────────────────────────
      $(document).on("keydown", "#feature + .selectize-control input, #feature-selectized", function(e) {
        if (e.key === "Enter") { e.preventDefault(); $("#GO").trigger("click"); }
      });
      // Also handle direct keydown on the feature container (selectize wraps the input)
      $(document).on("keydown", ".selectize-input input", function(e) {
        if (e.key === "Enter" && $(this).closest(".form-group").find("select").attr("id") === "feature") {
          e.preventDefault(); $("#GO").trigger("click");
        }
      });

      // ── Resizable plot boxes ──────────────────────────────────────────────
      function initResizable() {
        $(".plot-box").each(function() {
          if (!$(this).hasClass("ui-resizable")) {
            $(this).resizable({
              handles: "se",
              minWidth: 220,
              minHeight: 120,
              stop: function() {
                // Tell Shiny plots to re-render at new size
                $(window).trigger("resize");
              }
            });
          }
        });
      }
      $(document).on("shiny:idle", initResizable);

      // ── Proportions settings gear toggle ─────────────────────────────────
      $(document).on("click", "#pie-settings-toggle", function() {
        $("#pie-settings-panel").slideToggle(150);
        var icon = $(this).find("i");
        icon.toggleClass("fa-cog fa-times");
      });

      // ── Co-expression settings gear toggle ────────────────────────────────
      $(document).on("click", "#coexpr-settings-toggle", function() {
        $("#coexpr-settings-panel").slideToggle(150);
        var icon = $(this).find("i");
        icon.toggleClass("fa-cog fa-times");
      });

      // ── Trajectory settings gear toggle ───────────────────────────────────
      $(document).on("click", "#traj-settings-toggle", function() {
        $("#traj-settings-panel").slideToggle(150);
        var icon = $(this).find("i");
        icon.toggleClass("fa-cog fa-times");
      });

      // ── Toast notification handler ────────────────────────────────────────
      Shiny.addCustomMessageHandler("lvsc_toast", function(msg) {
        var $t = $("#lvsc-toast");
        $t.html("<i class=\'fa fa-check-circle\' style=\'color:#3d8b6e\'></i> " + msg);
        $t.addClass("lvsc-toast-show");
        setTimeout(function() { $t.removeClass("lvsc-toast-show"); }, 3200);
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
              fileInput("seurat_file", "Upload from local machine (.rds / .qs / .qs2)",
                        accept = c(".rds", ".qs", ".qs2")),
              tags$hr(),
              # --- GTF ---
              tags$label("GTF File", style = "font-weight: bold;"),
              fileInput("gtf", "Upload from local machine (.gtf)", accept = ".gtf"),
              tags$hr()
          ),
          
          div(id = "analysisForm",
              
              # ── Section 1: Gene ────────────────────────────────────
              tags$div(class = "sidebar-section",
                       tags$p(class = "sidebar-section-label", "Gene"),
                       uiOutput("recent_genes_ui"),
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

                       # Gene Assay
                       tags$div(class = "input-label-row",
                                tags$label("Gene Assay", `for` = "gene_assay"),
                                tags$span(class = "chip-info-icon",
                                          `data-toggle` = "popover", `data-trigger` = "hover focus",
                                          `data-placement` = "right",
                                          `data-content` = "The Seurat assay containing gene-level expression data (e.g. RNA). Must differ from the Isoform Assay.",
                                          tabindex = "0", "i")
                       ),
                       selectInput("gene_assay", NULL, choices = NULL),

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
              tags$div(style = "display: flex; gap: 5px; margin-bottom: 4px; align-items: stretch;",
                       actionButton("GO", "Run Analysis", icon = icon("play"),
                                    class = "btn btn-success",
                                    style = "flex: 3; font-size: 13px; font-weight: 600;",
                                    `aria-label` = "Run analysis for selected gene"),
                       actionButton("report_btn", "Report",
                                    icon  = icon("file-alt"),
                                    class = "btn btn-default",
                                    style = "flex: 1.5; font-size: 12px; background:#6c757d; color:white; border-color:#6c757d; padding: 6px 8px;",
                                    `aria-label` = "Download report"),
                       actionButton("resetBtn", "Reset", icon = icon("undo"),
                                    class = "btn btn-warning",
                                    style = "flex: 1; font-size: 12px; padding: 6px 8px;",
                                    `aria-label` = "Reset analysis and clear data")
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
                               fluidRow(style = "display:flex; align-items:stretch; flex-wrap:wrap;",
                                 column(6, style = "display:flex; flex-direction:column;",
                                        div(class = "plot-box", style = "flex:1;",
                                            tags$div(class = "plot-title-row", style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:6px; cursor:pointer; user-select:none;",
                                                     h3(class = "plot-title", style="margin:0;", textOutput("title_celltype_plot", inline = TRUE)),
                                                     tags$span(class = "chip-info-icon",
                                                               `data-toggle` = "popover",
                                                               `data-trigger` = "hover focus",
                                                               `data-placement` = "right",
                                                               `data-content` = "Cells plotted on the selected embedding (e.g. UMAP), coloured by the chosen Metadata Column.",
                                                               tabindex = "0", "i")
                                            ),
                                            div(style = "height:37px;"),
                                            shinycssloaders::withSpinner(plotOutput("celltype_plot"), type=4, color="#337ab7", size=0.7),
                                            div(style = "text-align:right; padding-right:10px; margin-top:4px;",
                                                actionButton("celltype_plot-download_celltype_plot",
                                                             "Download", icon = icon("download"),
                                                             class = "btn btn-sm btn-link",
                                                             style = "font-size:12px; color:#337ab7;"))
                                        )
                                 ),
                                 column(6, style = "display:flex; flex-direction:column;",
                                        div(class = "plot-box", style = "flex:1;",
                                            tags$div(class = "plot-title-row", style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:6px; cursor:pointer; user-select:none;",
                                                     h3(class = "plot-title", style="margin:0;", "Gene Expression"),
                                                     tags$span(class = "chip-info-icon",
                                                               `data-toggle` = "popover",
                                                               `data-trigger` = "hover focus",
                                                               `data-placement` = "right",
                                                               `data-content` = "Expression of the selected gene overlaid on the embedding. Brighter colour = higher expression.",
                                                               tabindex = "0", "i")
                                            ),
                                            div(style = "display:flex; align-items:center; justify-content:flex-end; gap:6px; margin-bottom:6px;",
                                              shinyWidgets::materialSwitch(
                                                inputId = "gene_viridis",
                                                label   = "Viridis palette",
                                                value   = FALSE,
                                                status  = "primary",
                                                inline  = TRUE
                                              )
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
                                                     h3(class = "plot-title", style="margin:0;", textOutput("title_vln_plot", inline = TRUE)),
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
                      tabPanel("Isoform Summary",
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
                      tabPanel("Expression Maps",
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            tags$div(class = "plot-title-row", style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:6px; cursor:pointer; user-select:none;",
                                                     h3(class = "plot-title", style="margin:0;", "Isoform Expression Plots"),
                                                     tags$span(class = "chip-info-icon",
                                                               `data-toggle` = "popover",
                                                               `data-trigger` = "hover focus",
                                                               `data-placement` = "right",
                                                               `data-content` = "One plot per selected isoform, showing its expression on the embedding. Toggle chips in the sidebar to add or remove isoforms.",
                                                               tabindex = "0", "i")
                                            ),
                                            div(style = "display:flex; align-items:center; justify-content:flex-end; gap:6px; margin-bottom:6px;",
                                              shinyWidgets::materialSwitch(
                                                inputId = "shared_scale",
                                                label   = "Shared colour scale",
                                                value   = FALSE,
                                                status  = "success",
                                                inline  = TRUE
                                              ),
                                              shinyWidgets::materialSwitch(
                                                inputId = "feature_viridis",
                                                label   = "Viridis palette",
                                                value   = FALSE,
                                                status  = "primary",
                                                inline  = TRUE
                                              )
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
                                            tags$div(class = "plot-title-row", style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:6px; cursor:pointer; user-select:none;",
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
                               tabsetPanel(id = "transcript_subtabs", type = "pills",
                                 tabPanel("Structure",
                                   fluidRow(style = "margin-top:8px;",
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
                                                div(style = "text-align:right; padding-right:10px; margin-top:4px; display:flex; gap:8px; justify-content:flex-end;",
                                                    actionButton("pseudobulk_heatmap-download_pseudobulk_heatmap",
                                                                 "Download Plot", icon = icon("download"),
                                                                 class = "btn btn-sm btn-link",
                                                                 style = "font-size:12px; color:#337ab7;"),
                                                    downloadButton("heatmap_download_csv", "Download Values",
                                                                   class = "btn btn-sm btn-link",
                                                                   style = "font-size:12px; color:#337ab7;"))
                                            )
                                     )
                                   )
                                 ),
                                 tabPanel("IsoViz",
                                   fluidRow(style = "margin-top:8px;",
                                     column(12,
                                            div(class = "plot-box", style = "padding:0; overflow:hidden;",
                                                tags$iframe(
                                                  src    = "https://isomix-test.stemformatics.org/isovis/",
                                                  width  = "100%",
                                                  height = "800px",
                                                  style  = "border:none; display:block;",
                                                  title  = "IsoViz isoform visualisation"
                                                )
                                            )
                                     )
                                   )
                                 )
                               )
                      ),
                      
                      # ── Tab 5: Proportions ────────────────────────────────────
                      tabPanel("Proportions",
                               fluidRow(
                                 column(12,
                                        div(style = "display:flex; justify-content:flex-end; align-items:center; padding: 4px 4px 0; margin-bottom:4px;",
                                            tags$button(
                                              id = "pie-settings-toggle",
                                              `aria-label` = "Toggle plot settings",
                                              tags$i(class = "fa fa-cog", style = "margin-right:4px;"), "Plot Settings"
                                            )
                                        ),
                                        div(id = "pie-settings-panel", style = "display:block;",
                                            div(style = "display:flex; align-items:flex-start; gap:24px; flex-wrap:wrap; margin-bottom:4px;",
                                                conditionalPanel("input.pie_display_mode !== 'bar'",
                                                  tags$label("Columns per row",
                                                             style = "font-size:13px; font-weight:500; display:block; margin-bottom:2px;"),
                                                  numericInput("pie_ncol", NULL, value = 4, min = 1, max = 10, step = 1, width = "90px")
                                                ),
                                                div(
                                                  tags$label("Min. counts threshold",
                                                             style = "font-size:13px; font-weight:500; display:block; margin-bottom:2px;"),
                                                  numericInput("pie_min_counts", NULL, value = 10, min = 0, step = 1, width = "90px")
                                                ),
                                                div(
                                                  tags$label("Display mode",
                                                             style = "font-size:13px; font-weight:500; display:block; margin-bottom:2px;"),
                                                  shinyWidgets::radioGroupButtons(
                                                    inputId  = "pie_display_mode",
                                                    label    = NULL,
                                                    choices  = c("Pie" = "proportion", "Bar" = "bar"),
                                                    selected = "proportion",
                                                    size     = "sm",
                                                    status   = "default"
                                                  )
                                                ),
                                                div(style = "min-width:200px; flex:1;",
                                                  div(class = "input-label-row",
                                                    tags$label("Show groups", style = "font-size:13px; font-weight:500; display:block; margin-bottom:2px; margin-right:4px;"),
                                                    tags$span(class = "chip-info-icon",
                                                      `data-toggle` = "popover", `data-trigger` = "hover focus",
                                                      `data-placement` = "right",
                                                      `data-content` = "Choose which cell groups to include in the proportions plot. Groups come from the Metadata Column selected in the sidebar.",
                                                      tabindex = "0", "i")
                                                  ),
                                                  shinyWidgets::pickerInput(
                                                    inputId  = "pie_cells_filter",
                                                    label    = NULL,
                                                    choices  = NULL,
                                                    multiple = TRUE,
                                                    options  = shinyWidgets::pickerOptions(
                                                      actionsBox         = TRUE,
                                                      selectedTextFormat = "count > 2",
                                                      countSelectedText  = "{0} groups selected",
                                                      liveSearch         = FALSE,
                                                      size               = 8
                                                    ),
                                                    width = "100%"
                                                  )
                                                )
                                            )
                                        )
                                 )
                               ),
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box", style = "margin-top:8px;",
                                            tags$div(class = "plot-title-row", style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:10px; cursor:pointer; user-select:none;",
                                                     h3(class = "plot-title", style="margin:0;", textOutput("title_pie_plot", inline = TRUE)),
                                                     tags$span(class = "chip-info-icon",
                                                               `data-toggle` = "popover",
                                                               `data-trigger` = "hover focus",
                                                               `data-placement` = "right",
                                                               `data-content` = "Each pie (or bar) shows isoform proportions or absolute counts within a cell type. Switch modes in Plot Settings.",
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

                               # ── Settings toggle ───────────────────────────────
                               fluidRow(
                                 column(12,
                                   div(style = "display:flex; justify-content:flex-end; padding: 4px 4px 0; margin-bottom:4px;",
                                     tags$button(id = "traj-settings-toggle", class = "settings-gear-btn",
                                       `aria-label` = "Toggle trajectory settings",
                                       tags$i(class = "fa fa-cog", style = "margin-right:4px;"), "Plot Settings")
                                   ),
                                   div(id = "traj-settings-panel", style = "display:block;",
                                     fluidRow(
                                       column(4,
                                   div(class = "sidebar-section",
                                     div(class = "input-label-row",
                                       p(class = "sidebar-section-label", style = "margin:0;",
                                         icon("filter"), " 1. Select ", textOutput("traj_grp_label", inline = TRUE)),
                                       tags$span(class = "chip-info-icon",
                                         `data-toggle` = "popover", `data-trigger` = "hover focus",
                                         `data-placement` = "right",
                                         `data-content` = "Choose which cell types from the Metadata Column to include. Remove any you don't want to show on the trajectory.",
                                         tabindex = "0", "i")
                                     ),
                                     uiOutput("traj_include_ui")
                                   )
                                 ),
                                 column(4,
                                   div(class = "sidebar-section",
                                     div(class = "input-label-row",
                                       p(class = "sidebar-section-label", style = "margin:0;",
                                         icon("arrows-alt-h"), " 2. Set ", textOutput("traj_grp_label2", inline = TRUE), " order"),
                                       tags$span(class = "chip-info-icon",
                                         `data-toggle` = "popover", `data-trigger` = "hover focus",
                                         `data-placement` = "right",
                                         `data-content` = "Drag cell types to set their left-to-right order on the x-axis. This defines the trajectory direction.",
                                         tabindex = "0", "i")
                                     ),
                                     uiOutput("traj_order_ui")
                                   )
                                 ),
                                 column(4,
                                   div(class = "sidebar-section",
                                     p(class = "sidebar-section-label",
                                       icon("info-circle"), " How it works"),
                                     tags$ul(style = "font-size:13px; color:#444; padding-left:16px; margin:0;",
                                       tags$li("Type a gene and press ", tags$strong("Run Analysis"), " in the sidebar."),
                                       tags$li("Select which cell types to include, then drag them into order left-to-right."),
                                       tags$li("Each dot is a single cell; the curved line is a loess smooth showing the expression trend across the ordered groups."),
                                       tags$li("Each isoform gets its own panel — compare the shapes to see when isoform usage shifts along the trajectory.")
                                     )
                                   )
                                 )
                                     )
                                   ) # end traj-settings-panel div
                                 )
                               ),

                               # ── Row 2: gene trajectory plot (full width) ──────
                               fluidRow(
                                 column(12,
                                   div(class = "plot-box",
                                     tags$div(class = "plot-title-row", style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:10px;",
                                       h3(class = "plot-title", style = "margin:0;", textOutput("title_trajectory_gene", inline = TRUE)),
                                       tags$span(class = "chip-info-icon",
                                         `data-toggle` = "popover", `data-trigger` = "hover focus",
                                         `data-placement` = "right",
                                         `data-content` = "Gene-level expression trajectory across the ordered cell states.",
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

                               # ── Row 3: isoform trajectory plot (full width) ───
                               fluidRow(
                                 column(12,
                                   div(class = "plot-box",
                                     tags$div(class = "plot-title-row", style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:6px;",
                                       h3(class = "plot-title", style = "margin:0;", textOutput("title_isoform_trajectory", inline = TRUE)),
                                       tags$span(class = "chip-info-icon",
                                         `data-toggle` = "popover", `data-trigger` = "hover focus",
                                         `data-placement` = "right",
                                         `data-content` = "One panel per selected isoform. Each dot is a cell, ordered left-to-right by the cell state order set above.",
                                         tabindex = "0", "i")
                                     ),
                                     div(style = "display:flex; align-items:center; justify-content:flex-end; gap:6px; margin-bottom:6px;",
                                       shinyWidgets::materialSwitch(
                                         inputId = "traj_shared_y",
                                         label   = "Shared Y axis",
                                         value   = FALSE,
                                         status  = "success",
                                         inline  = TRUE
                                       ),
                                       tags$span(class = "chip-info-icon",
                                         `data-toggle` = "popover", `data-trigger` = "hover focus",
                                         `data-placement` = "left",
                                         `data-content` = "When on, all isoform panels share the same Y axis range, making expression levels directly comparable across isoforms.",
                                         tabindex = "0", "i")
                                     ),
                                     shinycssloaders::withSpinner(plotOutput("trajectory_plot", height = "650px"), type=4, color="#337ab7", size=0.7),
                                     div(style = "text-align:right; padding-right:10px; margin-top:4px;",
                                       actionButton("trajectory_plot-download_trajectory_plot",
                                                    "Download", icon = icon("download"),
                                                    class = "btn btn-sm btn-link",
                                                    style = "font-size:12px; color:#337ab7;"))
                                   )
                                 )
                               )
                      ),



                      # ── Tab: Compare Two Isoforms ────────────────────────────
                      tabPanel("Compare Two Isoforms",
                        fluidRow(
                          column(12,
                            div(style = "display:flex; justify-content:flex-end; padding: 4px 4px 0; margin-bottom:4px;",
                              tags$button(id = "coexpr-settings-toggle", class = "settings-gear-btn",
                                `aria-label` = "Toggle plot settings",
                                tags$i(class = "fa fa-cog", style = "margin-right:4px;"), "Plot Settings")
                            ),
                            div(id = "coexpr-settings-panel", style = "display:block;",
                              fluidRow(
                                # ── Isoform X selector ──
                                column(4,
                            div(class = "sidebar-section",
                              div(class = "input-label-row",
                                p(class = "sidebar-section-label", style = "margin:0;", "Gene X"),
                                tags$span(class = "chip-info-icon",
                                  `data-toggle` = "popover", `data-trigger` = "hover focus",
                                  `data-placement` = "right",
                                  `data-content` = "Select a gene, then choose which isoform to plot on the X axis. Isoforms are filtered to those present in the isoform assay. Format: TranscriptID-GeneName.",
                                  tabindex = "0", "i")
                              ),
                              selectizeInput("scatter_gene_x", NULL,
                                choices = NULL, options = list(placeholder = "Type gene name…")),
                              selectInput("scatter_iso_x", "Isoform X",
                                choices = NULL, width = "100%")
                            )
                          ),
                          # ── Isoform Y selector ──
                          column(4,
                            div(class = "sidebar-section",
                              div(class = "input-label-row",
                                p(class = "sidebar-section-label", style = "margin:0;", "Gene Y"),
                                tags$span(class = "chip-info-icon",
                                  `data-toggle` = "popover", `data-trigger` = "hover focus",
                                  `data-placement` = "right",
                                  `data-content` = "Select a gene, then choose which isoform to plot on the Y axis. The gene can be the same as Gene X to compare isoforms within a gene. Format: TranscriptID-GeneName.",
                                  tabindex = "0", "i")
                              ),
                              selectizeInput("scatter_gene_y", NULL,
                                choices = NULL, options = list(placeholder = "Type gene name…")),
                              selectInput("scatter_iso_y", "Isoform Y",
                                choices = NULL, width = "100%")
                            )
                          ),
                          # ── Controls ──
                          column(4,
                            div(class = "sidebar-section",
                              p(class = "sidebar-section-label", icon("sliders-h"), " Options"),
                              # View toggle
                              div(style = "margin-bottom:10px;",
                                shinyWidgets::radioGroupButtons(
                                  inputId  = "compare_view",
                                  label    = NULL,
                                  choices  = c("Blend" = "blend", "Scatter" = "scatter", "Feature Map" = "feature"),
                                  selected = "blend",
                                  size     = "sm",
                                  status   = "default",
                                  width    = "100%"
                                )
                              ),
                              # Blend / Scatter: compare genes or isoforms
                              conditionalPanel("input.compare_view === 'blend' || input.compare_view === 'scatter'",
                                div(style = "margin-bottom:10px;",
                                  tags$label("Compare", style = "font-size:13px; font-weight:500; color:#444; display:block; margin-bottom:4px;"),
                                  conditionalPanel("input.compare_view === 'blend'",
                                    shinyWidgets::radioGroupButtons(
                                      inputId  = "blend_compare_type",
                                      label    = NULL,
                                      choices  = c("Isoforms" = "isoform", "Genes" = "gene"),
                                      selected = "isoform",
                                      size     = "sm",
                                      status   = "default",
                                      width    = "100%"
                                    )
                                  ),
                                  conditionalPanel("input.compare_view === 'scatter'",
                                    shinyWidgets::radioGroupButtons(
                                      inputId  = "scatter_compare_type",
                                      label    = NULL,
                                      choices  = c("Isoforms" = "isoform", "Genes" = "gene"),
                                      selected = "isoform",
                                      size     = "sm",
                                      status   = "default",
                                      width    = "100%"
                                    )
                                  )
                                )
                              ),
                              # Cell filter — shown for both scatter and feature views
                              div(class = "input-label-row",
                                tags$label("Show cells", style = "font-size:13px; font-weight:500; color:#444; margin:0;"),
                                tags$span(class = "chip-info-icon",
                                  `data-toggle` = "popover", `data-trigger` = "hover focus",
                                  `data-placement` = "right",
                                  `data-content` = "Choose which cell groups to include. Groups come from the Metadata Column selected in the sidebar. Cells are coloured by that column automatically.",
                                  tabindex = "0", "i")
                              ),
                              shinyWidgets::pickerInput(
                                inputId  = "scatter_cells_filter",
                                label    = NULL,
                                choices  = NULL,
                                multiple = TRUE,
                                options  = shinyWidgets::pickerOptions(
                                  actionsBox         = TRUE,
                                  selectedTextFormat = "count > 2",
                                  countSelectedText  = "{0} groups selected",
                                  liveSearch         = FALSE,
                                  size               = 8
                                ),
                                width = "100%"
                              )
                            )
                          )
                              ) # end inner fluidRow
                            ) # end coexpr-settings-panel
                          )
                        ),
                        fluidRow(
                          column(12,
                            div(class = "plot-box",
                              conditionalPanel("input.compare_view === 'scatter'",
                                tags$div(class = "plot-title-row",
                                  style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:10px;",
                                  h3(class = "plot-title", style = "margin:0;", "Isoform Expression Scatter"),
                                  tags$span(class = "chip-info-icon",
                                    `data-toggle` = "popover", `data-trigger` = "hover focus",
                                    `data-placement` = "right",
                                    `data-content` = "Each dot is a single cell. The orange line is a linear fit with 95% confidence band. R\u00b2 measures correlation strength. Marginal histograms + density curves show the expression distribution along each axis. Select both isoforms and press Plot to run.",
                                    tabindex = "0", "i")
                                ),
                                div(style = "display:flex; align-items:center; gap:18px; justify-content:flex-end; margin-bottom:8px;",
                                  div(style = "display:flex; align-items:center; gap:5px;",
                                    shinyWidgets::materialSwitch("scatter_zero_cells",
                                      "Include zeros", value = FALSE, status = "success", inline = TRUE),
                                    tags$span(class = "chip-info-icon",
                                      `data-toggle` = "popover", `data-trigger` = "hover focus",
                                      `data-placement` = "left",
                                      `data-content` = "When off, cells where both isoforms have zero expression are excluded. This reduces overplotting at the origin.",
                                      tabindex = "0", "i")
                                  ),
                                  div(style = "display:flex; align-items:center; gap:5px;",
                                    shinyWidgets::materialSwitch("scatter_show_marginal",
                                      "Show distributions", value = TRUE, status = "success", inline = TRUE),
                                    tags$span(class = "chip-info-icon",
                                      `data-toggle` = "popover", `data-trigger` = "hover focus",
                                      `data-placement` = "left",
                                      `data-content` = "Toggle the histogram and density plots on the top and right margins of the scatter plot.",
                                      tabindex = "0", "i")
                                  )
                                ),
                                plotOutput("scatter_plot", height = "600px"),
                                div(style = "text-align:right; padding-right:10px; margin-top:4px;",
                                  actionButton("scatter_plot-download_scatter_plot", "Download",
                                    icon = icon("download"),
                                    class = "btn btn-sm btn-link",
                                    style = "font-size:12px; color:#337ab7;"))
                              ),
                              conditionalPanel("input.compare_view === 'feature'",
                                tags$div(class = "plot-title-row",
                                  style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:10px;",
                                  h3(class = "plot-title", style = "margin:0;", "Isoform Feature Maps")
                                ),
                                div(style = "display:flex; align-items:center; gap:18px; justify-content:flex-end; margin-bottom:8px;",
                                  tags$p(style = "font-size:11px; margin:0;",
                                    HTML('<span style="color:#e6b800;">&#9632;</span> Both &nbsp; <span style="color:#2e5ba8;">&#9632;</span> Only X &nbsp; <span style="color:#b83232;">&#9632;</span> Only Y &nbsp; <span style="color:#d0d0d0;">&#9632;</span> Neither'))
                                )
                              ),
                              conditionalPanel("input.compare_view === 'blend'",
                                tags$div(class = "plot-title-row",
                                  style = "display:flex; align-items:center; justify-content:center; gap:6px; margin-bottom:10px;",
                                  h3(class = "plot-title", style = "margin:0;", "Isoform Blend"),
                                  tags$span(class = "chip-info-icon",
                                    `data-toggle` = "popover", `data-trigger` = "hover focus",
                                    `data-placement` = "right",
                                    `data-content` = "Seurat blend mode: red = isoform X only, blue = isoform Y only, purple = co-expressed. Requires two different isoforms.",
                                    tabindex = "0", "i")
                                )
                              ),
                              conditionalPanel("input.compare_view === 'feature' || input.compare_view === 'blend'",
                                shinycssloaders::withSpinner(plotOutput("compare_feature_plot", height = "500px"), type=4, color="#337ab7", size=0.7),
                                div(style = "text-align:right; padding-right:10px; margin-top:4px;",
                                  actionButton("compare_feature_plot-download_compare_feature_plot", "Download",
                                    icon = icon("download"),
                                    class = "btn btn-sm btn-link",
                                    style = "font-size:12px; color:#337ab7;"))
                              )
                            )
                          )
                        )
                      ),

          ) # end of tabsetPanel
        ) # end of mainPanel
      ) # end of sidebarLayout
  ) # end of mainUI div
)