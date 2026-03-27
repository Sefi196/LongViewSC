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
plotIsoformPie <- function(raw_data, gene, selected_isoforms = NULL, ncol = 4, min_counts = 0, colour_map = NULL, group_by_label = NULL, display_mode = "proportion") {

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

  # 3) Normalise transcript_id
  df_gene <- dplyr::mutate(df_gene,
                           stripped_id = sub("-[^-]+$", "",
                                             sub("[.].*$",  "", transcript_id)))

  # 4) Determine which isoforms to highlight
  if (!is.null(selected_isoforms) && length(selected_isoforms) > 0) {
    top_isoforms <- sub("-[^-]+$", "", selected_isoforms)
    top_isoforms <- sub("[.][^.]*$", "", top_isoforms)
    top_isoforms <- unique(top_isoforms)
  } else {
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

  df5_summed <- df5_summed %>%
    dplyr::mutate(proportion = pmin(pmax(proportion, 0), 1))

  # 6) Build colour map
  non_grey_levels <- sort(top_isoforms)
  if (!is.null(colour_map) && all(non_grey_levels %in% names(colour_map))) {
    hues <- colour_map[non_grey_levels]
  } else {
    hues <- stats::setNames(scales::hue_pal()(length(non_grey_levels)), non_grey_levels)
  }
  color_map <- c(hues, "Other isoforms" = "#cccccc")

  # 7) Add cell label with total counts
  totals <- dplyr::mutate(totals,
                          cell_label = paste0(cell_type,
                                              "\n(n = ", scales::comma(total_expression), ")"))

  df5_summed <- dplyr::left_join(
    df5_summed,
    dplyr::select(totals, cell_type, cell_label),
    by = "cell_type"
  )

  plot_title <- if (!is.null(group_by_label) && nchar(group_by_label) > 0)
    paste0("Isoform proportions \u2014 ", gene, " by ", group_by_label)
  else
    paste0("Isoform proportions \u2014 ", gene)

  base_theme <- ggplot2::theme(
    plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold", size = 17,
                                             family = "sans", color = "#2c3e50",
                                             margin = ggplot2::margin(b = 16)),
    strip.text       = ggplot2::element_text(size = 13, color = "#333333", face = "bold",
                                             family = "sans",
                                             margin = ggplot2::margin(t = 8, b = 2)),
    strip.background = ggplot2::element_rect(fill = "#f4f6f8", color = NA),
    legend.position  = "right",
    legend.title     = ggplot2::element_text(size = 13, face = "bold", family = "sans", color = "#333333"),
    legend.text      = ggplot2::element_text(size = 12, family = "sans", color = "#333333"),
    legend.key.size  = ggplot2::unit(1.2, "lines"),
    panel.spacing.x  = ggplot2::unit(1.6, "lines"),
    panel.spacing.y  = ggplot2::unit(1.2, "lines"),
    plot.margin      = ggplot2::margin(16, 16, 16, 16),
    plot.background  = ggplot2::element_rect(fill = "white", color = NA)
  )

  if (isTRUE(display_mode == "bar")) {
    # ── 8b) Bar: stacked proportions with raw total count labels on top ───────
    df_bar <- df5_summed %>%
      dplyr::left_join(dplyr::select(totals, cell_type, total_expression), by = "cell_type")

    top_labels <- dplyr::distinct(df_bar, cell_label, total_expression)

    ggplot2::ggplot(df_bar,
                    ggplot2::aes(x = cell_label, y = proportion, fill = transcript_group)) +
      ggplot2::geom_col(color = "black", linewidth = 0.3, width = 0.7) +
      ggplot2::geom_text(
        data = top_labels,
        ggplot2::aes(x = cell_label, y = 1.01,
                     label = scales::comma(total_expression, accuracy = 1)),
        inherit.aes = FALSE, size = 3.2, hjust = 0.5, vjust = 0, color = "#333333"
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::percent,
        expand = ggplot2::expansion(mult = c(0, 0.12))
      ) +
      ggplot2::scale_fill_manual(values = color_map,
                                 breaks = c(non_grey_levels, "Other isoforms")) +
      ggplot2::labs(title = plot_title, x = NULL, y = "Proportion", fill = "Isoform") +
      ggplot2::theme_minimal(base_size = 14) +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_text(angle = 35, hjust = 1, size = 13, family = "sans"),
        axis.text.y  = ggplot2::element_text(size = 13, family = "sans"),
        axis.title.y = ggplot2::element_text(size = 14, face = "bold", family = "sans")
      ) +
      base_theme
  } else {
    # ── 8a) Proportions: pie chart ────────────────────────────────────────────
    ggplot2::ggplot(df5_summed,
                    ggplot2::aes(x = "", y = proportion, fill = transcript_group)) +
      ggplot2::geom_col(color = "black", linewidth = 0.4, width = 1) +
      ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
                                  oob = scales::squish) +
      ggplot2::coord_polar(theta = "y") +
      ggplot2::facet_wrap(~ cell_label, ncol = ncol, strip.position = "bottom") +
      ggplot2::geom_text(
        ggplot2::aes(label = ifelse(proportion > 0.07,
                                    paste0(round(proportion * 100, 0), "%"), "")),
        position = ggplot2::position_stack(vjust = 0.5),
        color = "black", size = 4.5, fontface = "bold"
      ) +
      ggplot2::scale_fill_manual(values = color_map,
                                 breaks = c(non_grey_levels, "Other isoforms")) +
      ggplot2::labs(title = plot_title, fill = "Isoform") +
      ggplot2::theme_void(base_size = 14) +
      base_theme
  }
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
    colour_map    = NULL,
    group_by_label = NULL,
    display_mode  = "proportion"
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
  
  plotIsoformPie(raw_data = df_long, gene = gene, selected_isoforms = selected_isoforms,
                 ncol = ncol, min_counts = min_counts, colour_map = colour_map,
                 group_by_label = group_by_label, display_mode = display_mode)
}


# ─────────────────────────────────────────────────────────────────────────────
# plotIsoformSankey
#
# Interactive Sankey diagram: isoforms (left) → cell types (right).
# Flow width = proportion or absolute counts.
# Uses plotly for an interactive output.
# ─────────────────────────────────────────────────────────────────────────────
plotIsoformSankeyFromSeurat <- function(
    seurat_obj,
    gene,
    selected_isoforms = NULL,
    assay         = "iso",
    layer         = "counts",
    cell_type_col = "BroadType",
    min_counts    = 0,
    colour_map    = NULL,
    display_mode  = "proportion"
) {
  mat  <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer)
  meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(cell, cell_type = dplyr::all_of(cell_type_col))

  gene_regex    <- paste0("(^|-)", gene, "$")
  matched_feats <- grep(gene_regex, rownames(mat), value = TRUE)
  if (length(matched_feats) == 0)
    matched_feats <- grep(paste0("^", gene, "($|\\.)"), rownames(mat), value = TRUE)
  if (length(matched_feats) == 0)
    stop("No isoforms found for gene: ", gene)

  mat_sub <- mat[matched_feats, , drop = FALSE]
  df_long <- as.data.frame(as.table(as.matrix(mat_sub)))
  colnames(df_long) <- c("transcript_id", "cell", "expression")
  df_long <- df_long %>%
    dplyr::left_join(meta, by = "cell") %>%
    dplyr::filter(!is.na(cell_type), expression > 0) %>%
    dplyr::mutate(
      stripped_id = sub("-[^-]+$", "", sub("[.].*$", "", transcript_id))
    )

  # Determine isoforms to highlight (same logic as pie)
  if (!is.null(selected_isoforms) && length(selected_isoforms) > 0) {
    top_isoforms <- unique(sub("-[^-]+$", "", sub("[.][^.]*$", "", selected_isoforms)))
  } else {
    top_isoforms <- df_long %>%
      dplyr::group_by(stripped_id) %>%
      dplyr::summarise(tot = sum(expression), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(tot)) %>%
      dplyr::slice_head(n = 5) %>%
      dplyr::pull(stripped_id)
  }

  df_long <- dplyr::mutate(df_long,
    isoform_group = ifelse(stripped_id %in% top_isoforms, stripped_id, "Other isoforms"))

  # Aggregate per (isoform_group, cell_type)
  df_agg <- df_long %>%
    dplyr::group_by(isoform_group, cell_type) %>%
    dplyr::summarise(total = sum(expression, na.rm = TRUE), .groups = "drop")

  # Filter by min_counts (per cell_type total)
  if (min_counts > 0) {
    ct_totals <- df_agg %>%
      dplyr::group_by(cell_type) %>%
      dplyr::summarise(ct_total = sum(total), .groups = "drop")
    keep_ct <- ct_totals$cell_type[ct_totals$ct_total >= min_counts]
    df_agg  <- dplyr::filter(df_agg, cell_type %in% keep_ct)
  }
  if (nrow(df_agg) == 0)
    stop("No data remaining after filtering.")

  if (display_mode == "proportion") {
    df_agg <- df_agg %>%
      dplyr::group_by(cell_type) %>%
      dplyr::mutate(value = total / sum(total)) %>%
      dplyr::ungroup()
    value_fmt <- "proportion"
  } else {
    df_agg <- dplyr::mutate(df_agg, value = total)
    value_fmt <- "counts"
  }

  # Build node lists
  iso_nodes  <- sort(unique(df_agg$isoform_group))
  ct_nodes   <- sort(unique(df_agg$cell_type))
  all_nodes  <- c(iso_nodes, ct_nodes)
  node_idx   <- stats::setNames(seq_along(all_nodes) - 1L, all_nodes)

  # Node colours
  non_grey <- setdiff(iso_nodes, "Other isoforms")
  if (!is.null(colour_map)) {
    stripped_map <- colour_map
    names(stripped_map) <- sub("-[^-]+$", "", sub("[.][^.]*$", "", names(colour_map)))
    iso_colours <- sapply(iso_nodes, function(n) {
      if (n == "Other isoforms") "#cccccc"
      else if (n %in% names(stripped_map)) stripped_map[[n]]
      else scales::hue_pal()(length(non_grey))[match(n, non_grey)]
    })
  } else {
    pal <- stats::setNames(scales::hue_pal()(length(non_grey)), non_grey)
    iso_colours <- sapply(iso_nodes, function(n) {
      if (n == "Other isoforms") "#cccccc" else pal[[n]]
    })
  }
  ct_colours <- rep("#a8c5da", length(ct_nodes))
  node_colours <- c(iso_colours, ct_colours)

  sources <- node_idx[df_agg$isoform_group]
  targets <- node_idx[df_agg$cell_type]

  hover_label <- if (display_mode == "proportion")
    paste0(round(df_agg$value * 100, 1), "% of ", df_agg$cell_type)
  else
    paste0(scales::comma(round(df_agg$value)), " counts in ", df_agg$cell_type)

  plotly::plot_ly(
    type = "sankey",
    orientation = "h",
    node = list(
      label = all_nodes,
      color = node_colours,
      pad   = 15,
      thickness = 20,
      line  = list(color = "black", width = 0.5)
    ),
    link = list(
      source      = sources,
      target      = targets,
      value       = df_agg$value,
      customdata  = hover_label,
      hovertemplate = "%{customdata}<extra></extra>"
    )
  ) %>%
    plotly::layout(
      title = list(
        text = paste0("Isoform flow \u2014 ", gene,
                      " (", if (display_mode == "proportion") "proportions" else "counts", ")"),
        font = list(size = 15, color = "#2c3e50")
      ),
      font = list(size = 11),
      paper_bgcolor = "white"
    )
}