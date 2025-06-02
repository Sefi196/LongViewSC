


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
      sidebarLayout(
        sidebarPanel(
          width = 4,
          class = "sidebar-panel",
          
          # (1) File uploads stay OUTSIDE the resettable form:
          fileInput("seurat_file", "Upload Seurat Object (.rds)", accept = ".rds"),
          fileInput("gtf",         "Upload GTF (.gtf)",            accept = ".gtf"),
          
          # (2) EVERYTHING below goes inside the resettable “analysisForm”:
          div(id = "analysisForm",
              selectizeInput("feature", "Select a Gene:", choices = NULL, options  = list(placeholder = "Start typing…", maxOptions = 3000)),
              selectInput("reduction",     "Select Reduction Type",  choices = NULL),
              selectInput("isoform_assay", "Select Isoform Assay",   choices = NULL),
              selectInput("group_by",      "Select Metadata Column", choices = NULL),
              numericInput("number_of_isoforms", "Number of Isoforms to Plot", value = 4, min = 1, step = 1),
              #uiOutput("isoforms_warning"),
              actionButton("GO", "GO", class = "btn-block btn-lg btn-success"),
              actionButton("resetBtn", "Clear All", icon = icon("refresh"), class = "btn btn-warning btn-block")
          )  # end of analysisForm
        ),
    
        mainPanel(
          tabsetPanel(id = "main_tabs",
                      
                      # Gene Expression Tab
                      tabPanel("Gene Expression",
                               fluidRow(
                                 column(5,
                                        div(class = "plot-box",
                                            h3(class = "plot-title", "Gene Feature Plot"),
                                            plotOutput("feature_plot_gene"),
                                            div(
                                              actionButton(
                                                inputId = "feature_plot_gene-download_feature_plot_gene",
                                                label = "Download",
                                                icon = icon("download"),
                                                class = "btn btn-sm btn-link",
                                                title = "Download Plot",
                                                style = "font-size: 12px;"
                                              ),
                                              style = "text-align: right; padding-right: 10px; margin-top: 5px;"
                                            )
                                        )
                                 ),
                                 column(7,
                                        div(class = "plot-box",
                                            h3(class = "plot-title", "Cell Types"),
                                            plotOutput("celltype_plot"),
                                            div(
                                              actionButton(
                                                inputId = "celltype_plot-download_celltype_plot",
                                                label = "Download",
                                                icon = icon("download"),
                                                class = "btn btn-sm btn-link",
                                                title = "Download Plot",
                                                style = "font-size: 12px;"
                                              ),
                                              style = "text-align: right; padding-right: 10px; margin-top: 5px;"
                                            )
                                        )
                                 ),
                                 column(12,
                                        div(class = "plot-box",
                                            h3(class = "plot-title", "Violin Plot"),
                                            plotOutput("vln_plot"),
                                            div(
                                              actionButton(
                                                inputId = "vln_plot-download_vln_plot",
                                                label = "Download",
                                                icon = icon("download"),
                                                class = "btn btn-sm btn-link",
                                                title = "Download Plot",
                                                style = "font-size: 12px;"
                                              ),
                                              style = "text-align: right; padding-right: 10px; margin-top: 5px;"
                                            )
                                        )
                                 )
                               )
                      ),
                      
                      # Isoform Expression Tab
                      tabPanel("Isoform Expression",
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            h3(class = "plot-title", "Pick Isoform(s) to Plot"),
                                            DT::dataTableOutput("isoform_table")
                                        )
                                 )
                               ),
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            h3(class = "plot-title", "Isoform Feature Plot"),
                                            plotOutput("feature_plot_iso"),
                                            div(
                                              actionButton(
                                                inputId = "feature_plot_isoform-download_feature_plot_isoform",
                                                label = "Download",
                                                icon = icon("download"),
                                                class = "btn btn-sm btn-link",
                                                title = "Download Plot",
                                                style = "font-size: 12px;"
                                              ),
                                              style = "text-align: right; padding-right: 10px; margin-top: 5px;"
                                            )
                                        )
                                 )
                               ),
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            h3(class = "plot-title", "Dot Plot"),
                                            plotOutput("dot_plot_iso"),
                                            div(
                                              actionButton(
                                                inputId = "dot_plot_isoform-download_dot_plot_isoform",
                                                label = "Download",
                                                icon = icon("download"),
                                                class = "btn btn-sm btn-link",
                                                title = "Download Plot",
                                                style = "font-size: 12px;"
                                              ),
                                              style = "text-align: right; padding-right: 10px; margin-top: 5px;"
                                            ))
                                        )
                                 )
                               ),
                      
                      # Isoform Transcript Structure Tab
                      tabPanel("Isoform Transcript Structure",
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            h3(class = "plot-title", "Isoform Transcript Structure"),
                                            plotOutput("transcript_plot"),
                                            div(
                                              actionButton(
                                                inputId = "Isoform_TranscriptStructure-download_Isoform_TranscriptStructure",
                                                label = "Download",
                                                icon = icon("download"),
                                                class = "btn btn-sm btn-link",
                                                title = "Download Plot",
                                                style = "font-size: 12px;"
                                              ),
                                              style = "text-align: right; padding-right: 10px; margin-top: 5px;"
                                            )
                                        )
                                 )
                               ),
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            h3(class = "plot-title", "Pseudobulk Heatmap"),
                                            #downloadButton("download_heatmap_plot", "Download"),
                                            plotlyOutput("heatmap_plot"),
                                            div(
                                              actionButton(
                                                inputId = "pseudobulk_heatmap-download_pseudobulk_heatmap",
                                                label = "Download",
                                                icon = icon("download"),
                                                class = "btn btn-sm btn-link",
                                                title = "Download Plot",
                                                style = "font-size: 12px;"
                                              ),
                                              style = "text-align: right; padding-right: 10px; margin-top: 5px;"
                                            )
                                        )
                                 )
                               )
                               )
                      
          ) # end of tabsetPanel
        ) # end of mainPanel
      ) # end of sidebarLayout
  ) # end of mainUI div
)