# ─────────────────────────────────────────────────────────────────────────────
# plotIsoformPie
#
# Faceted pie chart showing each isoform's proportion of gene expression,
# one pie per cell type / grouping level.
#
# Arguments:
#   raw_data  : long-format df with columns: gene_id, cell_type,
#               transcript_id, expression
#   gene      : character. Gene symbol to plot (e.g. "NSD1")
#   top_n     : integer. Number of highest-proportion isoforms to highlight
#   ncol      : integer. Pies per row
# ─────────────────────────────────────────────────────────────────────────────
plotIsoformPie <- function(raw_data, gene, selected_isoforms = NULL, ncol = 4, min_counts = 0, colour_map = NULL) {
  
  # 1) Filter to this gene; sum expression per (cell_type, transcript_id)
  df_gene <- raw_data %>%
    dplyr::filter(gene_id == gene) %>%
    dplyr::group_by(gene_id, cell_type, transcript_id) %>%
    dplyr::summarise(isoform_expression = sum(expression, na.rm = TRUE),
                     .groups = "drop")
  
  # 2) Compute proportions within each cell_type
  df_gene <- df_gene %>%
    dplyr::group_by(gene_id, cell_type) %>%
    dplyr::mutate(
      total_expression = sum(isoform_expression, na.rm = TRUE),
      proportion       = ifelse(total_expression > 0,
                                isoform_expression / total_expression, 0)
    ) %>%
    dplyr::ungroup()
  
  totals <- dplyr::distinct(df_gene, cell_type, total_expression)
  
  # Filter out cell types below min_counts threshold
  if (min_counts > 0) {
    keep_types <- totals$cell_type[totals$total_expression >= min_counts]
    df_gene    <- dplyr::filter(df_gene, cell_type %in% keep_types)
    totals     <- dplyr::filter(totals, cell_type %in% keep_types)
    if (nrow(df_gene) == 0)
      stop("No cell types remaining after min_counts filter (threshold = ", min_counts, ")")
  }
  
  # 3) Normalise transcript_id to bare ID — strip dot-version (ENST00000x.10)
  #    AND -GENENAME suffix (BambuTx13041-PKM) so both conventions match
  #    the same format used in top_isoforms below.
  df_gene <- dplyr::mutate(df_gene,
                           stripped_id = sub("-[^-]+$", "",   # strip -GENE
                                             sub("[.].*$",  "", transcript_id)))  # strip .version
  
  # 4) Determine which isoforms to highlight (from chip selection, stripped to match stripped_id)
  if (!is.null(selected_isoforms) && length(selected_isoforms) > 0) {
    # Strip -GENENAME suffix from selected feature names to match stripped_id
    top_isoforms <- sub("-[^-]+$", "", selected_isoforms)
    # Also handle version suffixes (e.g. .9) — strip after last dot
    top_isoforms <- sub("[.][^.]*$", "", top_isoforms)
    top_isoforms <- unique(top_isoforms)
  } else {
    # Fallback: show top 5 by peak proportion
    peak_prop <- df_gene %>%
      dplyr::group_by(stripped_id) %>%
      dplyr::summarise(peak_prop = max(proportion, na.rm = TRUE), .groups = "drop")
    top_isoforms <- peak_prop %>%
      dplyr::arrange(dplyr::desc(peak_prop)) %>%
      dplyr::slice_head(n = 5) %>%
      dplyr::pull(stripped_id)
  }
  
  # 5) Group non-selected isoforms into "Other isoforms"
  df5 <- dplyr::mutate(df_gene,
                       transcript_group = ifelse(stripped_id %in% top_isoforms,
                                                 stripped_id, "Other isoforms"))
  
  df5_summed <- df5 %>%
    dplyr::group_by(gene_id, cell_type, transcript_group) %>%
    dplyr::summarise(proportion = sum(proportion, na.rm = TRUE), .groups = "drop")
  
  # Clamp proportions to [0, 1] — floating point can push sums just above 1,
  # which causes ggplot to silently clip the last slice off the pie
  df5_summed <- df5_summed %>%
    dplyr::mutate(proportion = pmin(pmax(proportion, 0), 1))
  
  # 6) Build colour map — use externally-supplied map if provided (keeps colours
  #    consistent with the sidebar chips), otherwise compute from hue_pal
  non_grey_levels <- sort(top_isoforms)
  if (!is.null(colour_map) && all(non_grey_levels %in% names(colour_map))) {
    hues <- colour_map[non_grey_levels]
  } else {
    hues <- stats::setNames(scales::hue_pal()(length(non_grey_levels)), non_grey_levels)
  }
  color_map <- c(hues, "Other isoforms" = "grey50")
  
  # 7) Add cell label with total counts
  totals <- dplyr::mutate(totals,
                          cell_label = paste0(cell_type, "\n(counts = ",
                                              scales::comma(total_expression), ")"))
  
  df5_summed <- dplyr::left_join(
    df5_summed,
    dplyr::select(totals, cell_type, cell_label),
    by = "cell_type"
  )
  
  # 8) Plot
  ggplot2::ggplot(df5_summed,
                  ggplot2::aes(x = "", y = proportion, fill = transcript_group)) +
    ggplot2::geom_col(color = "black", size = 0.2, width = 1) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
                                oob = scales::squish) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::facet_wrap(~ cell_label, ncol = ncol, strip.position = "bottom") +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(proportion > 0.05,
                                  paste0(round(proportion * 100, 1), "%"), "")),
      position = ggplot2::position_stack(vjust = 0.5),
      color = "white", size = 3
    ) +
    ggplot2::scale_fill_manual(values = color_map,
                               breaks = c(non_grey_levels, "Other isoforms")) +
    ggplot2::labs(
      title = paste0('Isoform Proportions for \u201c', gene,
                     '\u201d (', length(top_isoforms), ' selected isoforms shown)'),
      fill  = "Isoform"
    ) +
    ggplot2::theme_void(base_size = 14) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
      strip.text      = ggplot2::element_text(size = 10, face = "bold", vjust = 0.5),
      strip.placement = "outside",
      legend.position = "right",
      legend.title    = ggplot2::element_text(size = 12, face = "bold"),
      legend.text     = ggplot2::element_text(size = 9),
      panel.spacing.x = ggplot2::unit(1, "lines"),
      panel.spacing.y = ggplot2::unit(1, "lines"),
      plot.margin     = ggplot2::margin(10, 10, 10, 10)
    )
}


# ─────────────────────────────────────────────────────────────────────────────
# plotIsoformPieFromSeurat
#
# Wrapper around plotIsoformPie that pulls data directly from a Seurat object.
#
# Arguments:
#   seurat_obj    : Seurat object (layers already joined)
#   gene          : character. Gene symbol
#   top_n         : integer. Top isoforms to highlight
#   ncol          : integer. Pies per row
#   assay         : character. Assay name (default "iso")
#   layer         : character. Layer to use (default "counts")
#   cell_type_col : character. metadata column to group by
# ─────────────────────────────────────────────────────────────────────────────
plotIsoformPieFromSeurat <- function(
    seurat_obj,
    gene,
    selected_isoforms = NULL,
    ncol          = 4,
    assay         = "iso",
    layer         = "counts",
    cell_type_col = "BroadType",
    min_counts    = 0,
    colour_map    = NULL
) {
  Seurat::DefaultAssay(seurat_obj) <- assay
  
  mat  <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
  meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(cell, cell_type = dplyr::all_of(cell_type_col))
  
  # Match isoforms belonging to this gene. Two naming conventions:
  #   1. Standard: ENST00000123.1-VIM  (ends with -GENENAME)
  #   2. Bambu novel: BambuTx12041      (no suffix — gene IS the transcript id)
  gene_regex    <- paste0("(^|-)", gene, "$")           # ends with -GENE
  matched_feats <- grep(gene_regex, rownames(mat), value = TRUE)
  if (length(matched_feats) == 0) {
    # Fallback: exact match or prefix match (Bambu-style IDs)
    matched_feats <- grep(paste0("^", gene, "($|\\.)"), rownames(mat), value = TRUE)
  }
  if (length(matched_feats) == 0)
    stop("No isoforms found for gene: ", gene)
  
  # Build long-format table — only non-zero cells for speed
  mat_sub  <- mat[matched_feats, , drop = FALSE]
  df_long  <- as.data.frame(as.table(as.matrix(mat_sub)))
  colnames(df_long) <- c("transcript_id", "cell", "expression")
  
  df_long <- df_long %>%
    dplyr::left_join(meta, by = "cell") %>%
    dplyr::filter(!is.na(cell_type), expression > 0) %>%
    dplyr::mutate(gene_id = gene)
  
  plotIsoformPie(raw_data = df_long, gene = gene, selected_isoforms = selected_isoforms, ncol = ncol, min_counts = min_counts, colour_map = colour_map)
}