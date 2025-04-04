source("plot_gene_transcripts.R")
# Define plot for psedubulk visualizations 

#### testing function ###
#lets look at PKM as an example 
#ENSG00000067225 ENST00000335181.10 ENST00000319622.10
#seurat_obj = readRDS("/data/scratch/users/yairp/FLAMES_Day55/analysis/FLAMES_Single_sample_tutorial/bookdown_project/Day55_Long_read_tutorial/output_files/seu_objects/Day55_tutorial_gene_and_isoform_seurat.rds")
#group.by = "seurat_clusters"
#isoform_assay <- "iso"

########
# Define fucntion for psedubulk visualizations 
plot_pseudobulk_heatmap <- function(seurat_obj, group.by, isoforms_to_plot, isoform_assay) {
  # Aggregate expression (pseudobulk)

  #isoforms_to_plot = c("ENST00000335181.10-PKM", "ENST00000319622.10-PKM")
  
  aggregated_expr <- AggregateExpression(seurat_obj, 
                                         group.by = group.by, 
                                         assays = isoform_assay, 
                                         slot = "counts")$iso
  # Convert to matrix
  expr_mat <- as.matrix(aggregated_expr)
  
  # Normalize using CPM (Counts Per Million)
  cpm_mat <- t(t(expr_mat) / colSums(expr_mat) * 1e6)
  
  # Log-transform for better visualization
  log_cpm_mat <- log1p(cpm_mat)
  
  # Subset selected isoforms
  log_cpm_subset <- log_cpm_mat[rownames(log_cpm_mat) %in% isoforms_to_plot, , drop = FALSE]
  
  # Scale per isoform (row-wise z-score normalization)
  #scaled_mat <- t(scale(t(log_cpm_subset)))
  
  # Generate interactive heatmap
  p <- heatmaply(
    log_cpm_subset,
    scale = "none",  # No per-row normalization
    Showticklabels = c(TRUE, TRUE),
    dendrogram = "row")
  
  heatmaply(
    log_cpm_subset,
    scale = "none",  # No per-row normalization
    Showticklabels = c(TRUE, TRUE),
    dendrogram = "row")
}

# Call the plotting function with the necessary inputs
#plot_pseudobulk_heatmap(seurat_obj, 
#                        group.by = group.by, 
#                        isoforms_to_plot = c("ENST00000335181.10-PKM", "ENST00000319622.10-PKM"),
#                        isoform_assay = "iso")
