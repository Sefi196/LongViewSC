


ui <- fluidPage(
  useShinyjs(),
  
  # Add custom styles for a header banner and plot titles
  tags$head(
    tags$style(HTML("
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
        sidebarPanel(width = 4,
                     class = "sidebar-panel",
                     div(class = "sidebar-panel",
                         fileInput("seurat_file", "Upload Seurat Object (.rds)", accept = ".rds"),
                         fileInput("gtf", "Upload gtf (.gtf)", accept = ".gtf"),
                         selectizeInput("feature", "Select a Gene:", choices = NULL, options = list(placeholder = "Start typing...", maxOptions = 1000)),  
                         selectInput("reduction", "Select Reduction Type", choices = NULL),
                         selectInput("isoform_assay", "Select Isoform Assay", choices = NULL, selected = "RNA"),
                         selectInput("group_by", "Select Metadata Column", choices = NULL, selected = "seurat_clusters"),
                         numericInput("number_of_isoforms", "Number of Isoforms to Plot", value = 3, min = 1, step = 1),
                         uiOutput("isoforms_warning"),
                         numericInput("plot_width", "Plot width (inches)", value = 8, min = 4, max = 20),
                         numericInput("plot_height", "Plot height (inches)", value = 6, min = 4, max = 20),
                         actionButton("GO", "GO", class = "btn-block btn-lg btn-success")
                     )
        ),
        
        mainPanel(
          tabsetPanel(id = "main_tabs",
                      
                      # Gene Expression Tab
                      tabPanel("Gene Expression",
                               fluidRow(
                                 column(5,
                                        div(class = "plot-box",
                                            h3(class = "plot-title", "Gene Feature Plot"),
                                            downloadButton("download_feature_plot_gene", "Download"),
                                            plotOutput("feature_plot_gene")
                                        )
                                 ),
                                 column(7,
                                        div(class = "plot-box",
                                            h3(class = "plot-title", "Cell Types"),
                                            downloadButton("download_celltype_plot", "Download"),
                                            plotOutput("celltype_plot")
                                        )
                                 ),
                                 column(12,
                                        div(class = "plot-box",
                                            h3(class = "plot-title", "Violin Plot"),
                                            downloadButton("download_vln_plot", "Download"),
                                            plotOutput("vln_plot")
                                        )
                                 )
                               )
                      ),
                      
                      # Isoform Expression Tab
                      tabPanel("Isoform Expression",
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            h3(class = "plot-title", "Isoform Feature Plot"),
                                            downloadButton("download_feature_plot_isoform", "Download"),
                                            plotOutput("feature_plot_iso")
                                        )
                                 )
                               ),
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            h3(class = "plot-title", "Dot Plot"),
                                            downloadButton("download_dot_plot_isoform", "Download"),
                                            plotOutput("dot_plot_iso")
                                        )
                                 )
                               )
                      ),
                      
                      # Isoform Transcript Structure Tab
                      tabPanel("Isoform Transcript Structure",
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            h3(class = "plot-title", "Isoform Transcript Structure"),
                                            downloadButton("download_Isoform_TranscriptStructure", "Download"),
                                            plotOutput("transcript_plot")
                                        )
                                 )
                               ),
                               fluidRow(
                                 column(12,
                                        div(class = "plot-box",
                                            h3(class = "plot-title", "Pseudobulk Heatmap"),
                                            downloadButton("download_heatmap_plot", "Download"),
                                            plotlyOutput("heatmap_plot")
                                        )
                                 )
                               )
                      )
                      
          ) # end of tabsetPanel
        ) # end of mainPanel
      ) # end of sidebarLayout
  ) # end of mainUI div
)