
# Define plot for isoform visualizations 
plot_gene_transcripts <- function(isoforms_to_plot, gtf) {
  
  # remove gene symbol to get just ENST ID to each element in the vector
  isoforms_to_plot <- sub("-.*", "", isoforms_to_plot)
  
  # Filter for the specified isoforms
  filtered_gtf <- gtf %>% 
    filter(transcript_id %in% isoforms_to_plot)
  
  # Filter for exons only
  exon_data <- filtered_gtf %>% filter(type == "exon")
  
  # Create transcript structure plot
  plot <- ggplot(exon_data, aes(xstart = start, xend = end, y = transcript_id)) +
    geom_intron(data = to_intron(exon_data, "transcript_id"),
                aes(strand = strand), arrow.min.intron.length = 0) +
    geom_range(aes(fill = transcript_id), show.legend = FALSE, height = 0.3) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text.y = element_blank())
  
  return(plot)
}
