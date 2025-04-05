# A helper to register a downloadHandler for a given reactive plot and ID
add_download_pdf <- function(id, plot_reactive, width = 8, height = 6) {
  downloadHandler(
    filename = function() {
      paste0(id, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      pdf(file, width = width, height = height)
      print(plot_reactive())
      dev.off()
    }
  )
}