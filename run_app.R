# app.R (or run_app.R)
source("app/ui.R")      # Path to your UI file
source("app/server.R")  # Path to your server file

# Start the Shiny app
print("Starting ...")
print("please cite ....")
shinyApp(ui = ui, server = server)