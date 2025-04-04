library(shiny)
library(shinydashboard)

ui <- fluidPage(
  # Include the necessary styles for the dropdown menu and page
  tags$head(
    tags$style(HTML("
      .dropdown-menu {
        min-width: 300px;
      }
    "))
  ),
  
  # Create a dashboardHeader using dropdownMenu inside fluidPage
  tags$header(
    div(class = "navbar navbar-default", style = "background-color: #2c3e50; color: white; padding: 10px;",
        div(style = "display: flex; justify-content: space-between; align-items: center;",
            h2("Gene Expression Analysis App", style = "color: white;"),
            div(class = "dropdown", style = "position: relative; display: inline-block;",
                actionButton("mailbox", label = "", icon = icon("envelope"), class = "mailbox-button"),
                tags$ul(class = "dropdown-menu", style = "position: absolute; right: 0; top: 35px; list-style: none; background-color: white; box-shadow: 0px 8px 16px rgba(0,0,0,0.2); padding: 10px;",
                        tags$li(HTML('<a href="https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/physiology/Parker-laboratory-Metabolic-Proteomics" target="_blank"><i class="fa fa-user"></i><h4>About us</h4><p>Parker Laboratory</p></a>')),
                        tags$li(HTML('<a href="https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/clark-lab" target="_blank"><i class="fa fa-user"></i><h4>About us</h4><p>Clark Laboratory</p></a>')),
                        tags$li(HTML('<a href="mailto:sefi.prawer@unimelb.edu.au" target="_blank"><i class="fa fa-question"></i><h4>Support</h4><p>sefi.prawer@unimelb.edu.au</p></a>'))
                )
            )
        )
    )
  ),
  
  # Main Content Area Using fluidPage
  fluidPage(
    # Landing Page content
    div(id = "landingPage",
        h1("Welcome to the Gene Expression Analysis App"),
        div(style = "display: flex; gap: 10px;",
            actionButton("startBtn", "Start Analysis", class = "btn btn-primary btn-lg"),
            actionButton("instructionsBtn", "Instructions", class = "btn btn-secondary btn-lg")
        )
    ),
  
    # Instructions Page
    div(id = "instructionsPage", style = "display: none;", 
        fluidPage(
            h1("App Instructions")
        )
    ),
  
    # Main UI - Initially Hidden
    div(id = "mainUI", style = "display: none;",
        sidebarLayout(
            sidebarPanel(width = 4,
                         class = "sidebar-panel",
                         fileInput("seurat_file", "Upload Seurat Object (.rds)", accept = ".rds"),
                         fileInput("gtf", "Upload gtf (.gtf)", accept = ".gtf"),
                         selectizeInput("feature", "Select a Gene:", choices = NULL, options = list(placeholder = "Start typing...", maxOptions = 1000)),  
                         selectInput("reduction", "Select Reduction Type", choices = NULL),
                         selectInput("isoform_assay", "Select Isoform Assay", choices = NULL, selected = "RNA"),
                         selectInput("group_by", "Select Metadata Column", choices = NULL, selected = "seurat_clusters"),
                         numericInput("number_of_isoforms", "Number of Isoforms to Plot", value = 3, min = 1, step = 1),
                         uiOutput("isoforms_warning"),
                         actionButton("GO", "GO", class = "btn-block btn-lg btn-success")
            ),
    
            mainPanel(
                tabsetPanel(id = "main_tabs",  # Add an ID to control tab switching
                            tabPanel("Gene Expression",
                                     fluidRow(
                                         column(5,
                                                div(class = "plot-box",
                                                    h3(class = "plot-title", "Gene Feature Plot"),
                                                    plotOutput("feature_plot_gene")
                                                )
                                         ), 
                                         column(7,
                                                div(class = "plot-box",
                                                    h3(class = "plot-title", "Cell Types"),
                                                    plotOutput("celltype_plot")
                                                )
                                         ),
                                         column(12,
                                                div(class = "plot-box",
                                                    h3(class = "plot-title", "Violin Plot"),
                                                    plotOutput("vln_plot")))
                                     )
                            ),
    
                            tabPanel("Isoform Expression",
                                     fluidRow(
                                         column(12,
                                                div(class = "plot-box",
                                                    h3(class = "plot-title", "Isoform Feature Plot"),
                                                    plotOutput("feature_plot_iso")
                                                )
                                         )
                                     ),
                                     fluidRow(
                                         column(12,
                                                div(class = "plot-box",
                                                    h3(class = "plot-title", "Dot Plot"),
                                                    plotOutput("dot_plot_iso")
                                                )
                                         )
                                     )
                            ),
    
                            tabPanel("Isoform Transcript Structure",  
                                     fluidRow(
                                         column(12,
                                                div(class = "plot-box",
                                                    h3(class = "plot-title", "Isoform Transcript Structure"),
                                                    plotOutput("transcript_plot")
                                                )
                                         ),
                                         column(12,
                                                div(class = "plot-box",
                                                    h3(class = "plot-title", "Pseudobulk Heatmap"),
                                                    plotlyOutput("heatmap_plot")))))))))
)

server <- function(input, output, session) {
  # Hide everything except landing page initially
  shinyjs::hide("instructionsPage")
  shinyjs::hide("mainUI")
  
  # When the mailbox button is clicked, show a modal with email info
  observeEvent(input$mailbox, {
    showModal(modalDialog(
      title = "Contact Information",
      tags$ul(
        tags$li(HTML('<a href="https://biomedicalsciences.unimelb.edu.au/sbs-research-groups/anatomy-and-physiology-research/stem-cell-and-developmental-biology/clark-lab" target="_blank"><i class="fa fa-user"></i> About us: Clark Laboratory</a>')),
        tags$li(HTML('<a href="mailto:sefi.prawer@unimelb.edu.au" target="_blank"><i class="fa fa-question"></i> Support: sefi.prawer@unimelb.edu.au</a>'))
      ),
      easyClose = TRUE,   # Allows closing the modal
      footer = NULL       # No footer button in the modal
    ))
  })
  
  # When 'Start Analysis' is clicked, hide landing page and show main UI
  observeEvent(input$startBtn, {
    shinyjs::hide("landingPage")
    shinyjs::hide("instructionsPage")  # Ensure instructions page is hidden
    shinyjs::show("mainUI")
  })
  
  # When 'Instructions' is clicked, hide landing page and show instructions
  observeEvent(input$instructionsBtn, {
    shinyjs::hide("landingPage")
    shinyjs::hide("mainUI")  # Ensure main UI is hidden
    shinyjs::show("instructionsPage")
  })
  
  # Add functionality for back button from instructions to landing page
  observeEvent(input$backToLanding, {
    shinyjs::hide("instructionsPage")
    shinyjs::hide("mainUI")  # Ensure main UI is hidden
    shinyjs::show("landingPage")
  })
}

shinyApp(ui = ui, server = server)


#tabPanel("Welcome",
#         fluidRow(
#           column(12,
#                  includeMarkdown("welcome.md"),  # Embed Markdown file
#                  actionButton("start_analysis", "Start Analysis", 
#                               class = "btn btn-primary btn-lg")  # Button to switch tabs
#           )
#         )
#)