# app.R
# Start the Shiny app
print("LongViewSC Starting ....")
print("please cite ...")

source("app/global.R")
source("app/ui.R")      # Path to your UI file
source("app/server.R")  # Path to your server file


shinyApp(ui = ui, server = server)