# A helper to register a downloadHandler for a given reactive plot and ID
# Helper function
create_plot_download_handler <- function(plot_expr, prefix, width_input, height_input) {
  downloadHandler(
    filename = function() {
      paste0(prefix, "_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      ggsave(
        file,
        plot = plot_expr(),
        width = width_input(),
        height = height_input()
      )
    }
  )
}