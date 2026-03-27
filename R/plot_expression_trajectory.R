# ─────────────────────────────────────────────────────────────────────────────
# plotExpressionTrajectory
#
# Plots isoform expression across an ordered set of cell types.
# Each dot = one cell. A loess smooth overlays each panel.
# Faceted: one panel per isoform.
#
# Arguments:
#   seurat_obj     : Seurat object (layers already joined)
#   isoforms       : character vector of isoform feature names to plot
#   cell_type_col  : metadata column whose values define the x-axis groups
#   ordered_levels : character vector giving the desired left→right order
#   assay          : assay name (e.g. "iso")
#   layer          : layer to pull (default "data" = normalised)
# ─────────────────────────────────────────────────────────────────────────────
plotExpressionTrajectory <- function(
    seurat_obj,
    isoforms,
    cell_type_col,
    ordered_levels,
    assay = "iso",
    layer = "data",
    shared_y = FALSE
) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # Ensure ordered_levels is a plain character vector (bucket_list returns a list)
  ordered_levels <- as.character(unlist(ordered_levels))
  
  # ── 1. Pull expression matrix ──────────────────────────────────────────────
  mat     <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
  present <- intersect(isoforms, rownames(mat))
  if (length(present) == 0) stop("None of the requested isoforms found in assay.")
  
  mat_sub <- mat[present, , drop = FALSE]
  
  # ── 2. Long-format expression table ───────────────────────────────────────
  # check.names=FALSE preserves dashes in feature names (e.g. ENST00000380394.9-RPS6)
  df_wide <- as.data.frame(as.matrix(t(mat_sub)), check.names = FALSE)
  df_wide$cell <- rownames(df_wide)
  
  meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(cell, cell_type = dplyr::all_of(cell_type_col))
  
  df_long <- df_wide %>%
    tidyr::pivot_longer(-cell, names_to = "isoform", values_to = "expression") %>%
    dplyr::left_join(meta, by = "cell") %>%
    dplyr::filter(cell_type %in% ordered_levels) %>%
    dplyr::mutate(
      cell_type = factor(cell_type, levels = ordered_levels),
      # Label: strip the trailing -GENENAME to show just the transcript ID
      isoform_label = sub("-[^-]+$", "", isoform)
    )
  
  # ── 3. Colour palette (one hue per ordered level) ─────────────────────────
  pal <- setNames(
    scales::hue_pal()(length(ordered_levels)),
    ordered_levels
  )
  
  # ── 4. Build numeric x for the smooth so geom_smooth works across groups ──
  df_long <- df_long %>%
    dplyr::mutate(x_num = as.integer(cell_type))
  
  # ── 5. Plot ────────────────────────────────────────────────────────────────
  ggplot2::ggplot(df_long,
                  ggplot2::aes(x = cell_type, y = expression, colour = cell_type)) +
    
    # jittered cells
    ggplot2::geom_jitter(alpha = 0.35, size = 0.9, width = 0.22, height = 0) +
    
    # smooth across the numeric x so it spans all groups
    ggplot2::geom_smooth(
      ggplot2::aes(x = x_num, y = expression, group = 1),
      method  = "loess",
      formula = y ~ x,
      colour  = "black",
      fill    = "grey70",
      linewidth = 1,
      se      = TRUE,
      alpha   = 0.25
    ) +
    
    ggplot2::scale_colour_manual(values = pal, guide = "none") +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.5)) +
    
    ggplot2::facet_wrap(~ isoform_label, scales = if (shared_y) "fixed" else "free_y", ncol = 2) +
    
    ggplot2::labs(
      x = cell_type_col,
      y = "Normalised Expression"
    ) +
    
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "#2c3e50"),
      strip.text       = ggplot2::element_text(colour = "white", face = "bold", size = 13, family = "sans"),
      axis.text.x      = ggplot2::element_text(angle = 25, hjust = 1, size = 13, family = "sans"),
      axis.text.y      = ggplot2::element_text(size = 13, family = "sans"),
      axis.title       = ggplot2::element_text(face = "bold", size = 14, family = "sans"),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin      = ggplot2::margin(10, 15, 10, 10)
    )
}
