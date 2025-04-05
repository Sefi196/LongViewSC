#server side

server <- function(input, output, session) {
  # Hide everything except landing page initially
  shinyjs::hide("instructionsPage")
  shinyjs::hide("mainUI")
  
  # When 'Start Analysis' is clicked, hide landing page and show main UI
  observeEvent(input$startBtn, {
    shinyjs::hide("landingPage")
    shinyjs::hide("instructionsPage")  # Ensure instructions page is hidden
    shinyjs::show("mainUI")
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
  observeEvent(input$backToLanding, {
    shinyjs::hide("instructionsPage")
    shinyjs::hide("mainUI")  # Ensure main UI is hidden
    shinyjs::show("landingPage")
    
    # If in demo mode, reset the Seurat object and UI
    if (app_state$demo_mode) {
      seurat_obj(NULL)  # Reset the Seurat object
      gtf(NULL)  # Reset the Seurat object
      
      # Reset all plots when the demo data is loaded
      output$heatmap_plot <- renderPlotly({ NULL })
      output$feature_plot_gene <- renderPlot({ NULL })
      output$vln_plot <- renderPlot({ NULL })
      output$celltype_plot <- renderPlot({ NULL })
      output$feature_plot_iso <- renderPlot({ NULL })
      output$dot_plot_iso <- renderPlot({ NULL })
      output$transcript_plot <- renderPlot({ NULL })
      
      app_state$demo_mode <- FALSE  # Reset the flag
      
      # clear UI inputs
      updateSelectInput(session, "reduction", choices = character(0))
      updateSelectInput(session, "isoform_assay", choices = character(0))
      updateSelectInput(session, "group_by", choices = character(0))
      updateSelectizeInput(session, "feature", choices = NULL, options = list(placeholder = "Start typing...", maxOptions = 1000))
    }
  })
  
  # Demo Button functionality
  observeEvent(input$DemoBtn, {
    app_state$demo_mode <- TRUE  # <--- Track that demo data was loaded
    
    # Load the Seurat object from a predefined demo file
    demo_seurat_path <- "demo_data/Day55_tutorial_gene_and_isoform_seurat.rds"  # Update this path as necessary
    
    # Show spinner
    shinyjs::show("spinner")  # Assuming you have a div with id="spinner"
    demo_seurat_obj <- readRDS(demo_seurat_path)
    # Hide spinner after loading data
    
    # Update the reactive value with the demo Seurat object
    seurat_obj(demo_seurat_obj)
    
    # Dynamically update dropdown menus for reductions, assays, and features after demo Seurat object is loaded
    observe({
      req(seurat_obj())  # Ensure the Seurat object is available before updating UI
      print("Updating UI based on Demo Seurat object.")  # Debugging print
      
      # Get available reductions, assays, and group_by
      reductions <- names(seurat_obj()@reductions)
      isoform_assay <- names(seurat_obj()@assays)
      group_by <- colnames(seurat_obj()@meta.data)[!sapply(seurat_obj()@meta.data, is.numeric)]
      
      # Update UI elements based on the Seurat object
      updateSelectInput(session, "reduction", choices = reductions, selected = reductions[1])
      updateSelectInput(session, "isoform_assay", choices = isoform_assay, selected = isoform_assay[1])
      updateSelectInput(session, "group_by", choices = group_by, selected = group_by[1])
      updateSelectizeInput(session, "feature", choices = rownames(seurat_obj()), server = TRUE)
    })
    
    gtf_path_demo <- "demo_data/demo_isoform_annotated.gtf"  # Get the uploaded file path
    # Import GTF file using rtracklayer
    gtf_obj_demo <- rtracklayer::import(gtf_path_demo) %>% as_tibble()
    print("reading in demo GTF.")  # Debugging print
    # Update the reactive value with the new GTF object
    gtf(gtf_obj_demo)
    
    # Show notification when demo data is loaded
    showNotification("Demo data loaded successfully! Ready for analysis.", type = "message", duration = 5)
    
    # Show the main UI after demo data is loaded
    shinyjs::hide("landingPage")
    shinyjs::hide("spinner")
    shinyjs::show("mainUI")
  })
  
  #set the state of the app to control the demo mode reseting the input objects
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
    obj <- readRDS(input$seurat_file$datapath)  # Load the Seurat object
    
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
      
      # Update UI elements based on the Seurat object
      updateSelectInput(session, "reduction", choices = reductions, selected = reductions[1])
      updateSelectInput(session, "isoform_assay", choices = isoform_assay, selected = isoform_assay[1])
      updateSelectInput(session, "group_by", choices = group_by, selected = group_by[1])
      updateSelectizeInput(session, "feature", choices = rownames(seurat_obj()), server = TRUE)
    })
  })
  
  # Handle GTF file upload
  observeEvent(input$gtf, {
    req(input$gtf)  # Ensure a file is uploaded
    gtf_path <- input$gtf$datapath  # Get the uploaded file path
    
    # Import GTF file using rtracklayer
    gtf_obj <- rtracklayer::import(gtf_path) %>% as_tibble()
    
    # Update the reactive value with the new GTF object
    gtf(gtf_obj)
  })
  
  
  # Reactive expression triggered by "GO" button (updates everything)
  # Extract features related to isoforms
  filtered_data <- eventReactive(input$GO, {
    req(seurat_obj(), input$feature)
    
    assay_name <- input$isoform_assay
    isoform_features <-rownames(seurat_obj()@assays[[assay_name]]@features)
    
    matching_features <- grep(paste0("(^|-|\\b)", input$feature, "($|\\b)"), isoform_features, value = TRUE)
    
    list(
      celltype_plot = DimPlot(seurat_obj(), reduction = input$reduction, group.by = input$group_by),
      feature_plot_gene = FeaturePlot(seurat_obj(), features = input$feature, reduction = input$reduction), # Gene Feature Plot
      vln_plot = VlnPlot(seurat_obj(), features = input$feature, group.by = input$group_by),  # Gene Violin Plot
      isoform_features = matching_features  # Store matching isoform features
    )
  })
  
  # Reactive function to get the top N isoform features
  # provide warning if too many isoforms selected 
  isoform_features_to_plot <- eventReactive(input$GO,{
    req(filtered_data()$isoform_features, input$number_of_isoforms)
    
    # Get the matching isoform features
    isoform_features <- filtered_data()$isoform_features
    matching_features <- isoform_features  # Assuming you have already filtered by gene or feature
    
    # Define gene_name dynamically, assuming it's input by user
    gene_name <- input$feature  # Adjust this to your actual input for the gene name
    
    # Calculate the number of isoforms
    num_isoforms <- length(matching_features)
    
    # Check if the requested number of isoforms exceeds the available matching isoforms
    if (input$number_of_isoforms > num_isoforms) {
      # Create the dynamic message
      message <- paste(gene_name, "has", num_isoforms, "isoforms.")
      
      # Display the warning message in the sidebar
      output$isoforms_warning <- renderUI({
        HTML(paste("<p style='color:red;font-size:20px; text-align:center;'>", message, "</p>"))
      })
    } else {
      # Clear the warning message if the condition is not met
      output$isoforms_warning <- renderUI({
        NULL
      })
    }
    
    # Return the top N isoform features to plot
    head(matching_features, input$number_of_isoforms)
  })
  
  # Reactive function for the isoform Feature Plot
  isoform_plot <- eventReactive(input$GO,{
    req(isoform_features_to_plot())
    FeaturePlot(seurat_obj(), features = isoform_features_to_plot(), reduction = input$reduction)
  })
  
  # Reactive function for the isoform DotPlot
  dotplot_isoform <- eventReactive(input$GO,{
    req(isoform_features_to_plot())
    DotPlot(seurat_obj(), features = isoform_features_to_plot(), assay = input$isoform_assay, group.by = input$group_by) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels to 45 degrees
  })
  
  
  # Reactive expression to trigger the heatmap plot when the button is clicked
  reactive_heatmap <- eventReactive(input$GO, {
    # Ensure necessary data is available before proceeding
    req(isoform_features_to_plot())
    
    # Call the plotting function with the necessary inputs
    plot_pseudobulk_heatmap(
      seurat_obj = seurat_obj(), 
      group.by = input$group_by, 
      isoforms_to_plot = isoform_features_to_plot(),
      isoform_assay = input$isoform_assay
    )
  })
  
  # Render the main plots
  #### Genes ####
  output$feature_plot_gene <- renderPlot({
    filtered_data()$feature_plot })
  
  output$vln_plot <- renderPlot({
    filtered_data()$vln_plot})
  
  output$celltype_plot <- renderPlot({
    filtered_data()$celltype_plot})
  
  #### Isoforms ####
  output$feature_plot_iso <- renderPlot({
    isoform_plot()})
  
  output$dot_plot_iso <- renderPlot({
    dotplot_isoform()})
  
  #### Stack 
  output$transcript_plot <- renderPlot({
    # Ensure the necessary data is available before proceeding
    req(isoform_features_to_plot)  
    # Check if GTF file is empty or missing
    if (is.null(gtf()) || nrow(gtf()) == 0) {
      # Show a modal with a custom message if the GTF is blank
      showModal(modalDialog(
        title = "Missing GTF File",
        "Please provide a valid GTF to build hte isoform stack.",
        easyClose = TRUE,
        footer = NULL
      ))
    } else {
      # If GTF file is not blank, proceed with the plotting function
      plot_gene_transcripts(
        isoforms_to_plot = isoform_features_to_plot(),  # Pass the isoform features
        gtf = gtf()  # Pass the GTF data for plotting
      )
    }
  })
  
  # Render the heatmap plot when the button is clicked
  output$heatmap_plot <- renderPlotly({
    reactive_heatmap()  # Render the heatmap when the button is pressed
  })
}