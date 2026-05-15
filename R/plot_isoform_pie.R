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
  color_map <- c(hues, "Other isoforms" = "#d6d6d6")

  # 7) Keep totals as-is; cell_type used directly for labels

  df5_summed <- dplyr::left_join(
    df5_summed,
    dplyr::select(totals, cell_type, total_expression),
    by = "cell_type"
  )

  plot_title <- if (!is.null(group_by_label) && nchar(group_by_label) > 0)
    paste0("Isoform proportions \u2014 ", gene, " by ", group_by_label)
  else
    paste0("Isoform proportions \u2014 ", gene)

  # Shared legend theme
  legend_theme <- ggplot2::theme(
    legend.position      = "right",
    legend.title         = ggplot2::element_text(size = 12, face = "bold",
                                                  family = "sans", color = "#2c3e50"),
    legend.text          = ggplot2::element_text(size = 11, family = "sans",
                                                  color = "#444444"),
    legend.key.size      = ggplot2::unit(1.2, "lines"),
    legend.background    = ggplot2::element_blank(),
    legend.key           = ggplot2::element_blank(),
    legend.spacing.y     = ggplot2::unit(0.4, "lines")
  )

  if (isTRUE(display_mode == "bar")) {
    # ── Bar: horizontal stacked proportions ──────────────────────────────────
    lvls <- c(non_grey_levels, "Other isoforms")
    df5_summed$transcript_group <- factor(df5_summed$transcript_group, levels = lvls)

    df_bar <- df5_summed %>%
      dplyr::group_by(cell_type) %>%
      dplyr::mutate(cum_prop = cumsum(proportion),
                    label_y  = cum_prop - proportion / 2) %>%
      dplyr::ungroup()

    ggplot2::ggplot(df_bar,
                    ggplot2::aes(y = cell_type, x = proportion,
                                 fill = transcript_group)) +
      ggplot2::geom_col(color = "#555555", linewidth = 0.4, width = 0.68) +
      ggplot2::geom_text(
        ggplot2::aes(x     = label_y,
                     label = ifelse(proportion >= 0.07,
                                    paste0(round(proportion * 100), "%"), "")),
        color = "#111111", size = 4, fontface = "bold", family = "sans"
      ) +
      ggplot2::geom_text(
        data = dplyr::distinct(df_bar, cell_type, total_expression),
        ggplot2::aes(y = cell_type, x = 1.01,
                     label = paste0("n = ", scales::comma(total_expression, accuracy = 1))),
        inherit.aes = FALSE, size = 3.4, hjust = 0,
        color = "#000000", fontface = "italic", family = "sans"
      ) +
      ggplot2::scale_x_continuous(
        labels = scales::percent,
        expand = ggplot2::expansion(mult = c(0, 0.18)),
        limits = c(0, NA)
      ) +
      ggplot2::scale_fill_manual(values = color_map,
                                 breaks = c(non_grey_levels, "Other isoforms")) +
      ggplot2::labs(title = plot_title, x = "Proportion", y = NULL, fill = "Isoform") +
      ggplot2::theme_minimal(base_size = 13, base_family = "sans") +
      ggplot2::theme(
        plot.title         = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16,
                                                    color = "#2c3e50",
                                                    margin = ggplot2::margin(b = 16)),
        axis.text.y        = ggplot2::element_text(size = 13, color = "#111111",
                                                    face = "bold"),
        axis.text.x        = ggplot2::element_text(size = 13, color = "#111111"),
        axis.title.x       = ggplot2::element_text(size = 13, face = "bold",
                                                     color = "#111111",
                                                     margin = ggplot2::margin(t = 8)),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_line(color = "#e0e0e0", linewidth = 0.5),
        plot.background    = ggplot2::element_rect(fill = "white", color = NA),
        plot.margin        = ggplot2::margin(16, 30, 16, 16)
      ) +
      legend_theme

  } else {
    # ── Pie chart ─────────────────────────────────────────────────────────────
    lvls <- c(non_grey_levels, "Other isoforms")
    df5_summed$transcript_group <- factor(df5_summed$transcript_group, levels = lvls)

    # Include n count in the facet strip label
    n_lookup <- dplyr::distinct(df5_summed, cell_type, total_expression) %>%
      dplyr::mutate(cell_type_label = paste0(cell_type, "\nn = ",
                                              scales::comma(total_expression, accuracy = 1)))
    df5_summed <- dplyr::left_join(df5_summed,
                                    dplyr::select(n_lookup, cell_type, cell_type_label),
                                    by = "cell_type")
    ordered_labels <- n_lookup$cell_type_label[match(unique(df5_summed$cell_type),
                                                      n_lookup$cell_type)]
    df5_summed$cell_type_label <- factor(df5_summed$cell_type_label,
                                          levels = ordered_labels)

    ggplot2::ggplot(df5_summed,
                    ggplot2::aes(x = 0.5, y = proportion, fill = transcript_group)) +
      ggplot2::geom_col(color = "grey65", linewidth = 0.35, width = 1) +
      ggplot2::geom_text(
        ggplot2::aes(label = ifelse(proportion >= 0.06,
                                    paste0(round(proportion * 100), "%"), "")),
        position  = ggplot2::position_stack(vjust = 0.5),
        color     = "#111111", size = 4.2, fontface = "bold", family = "sans"
      ) +
      ggplot2::coord_polar(theta = "y", start = 0) +
      ggplot2::xlim(c(0, 1)) +
      ggplot2::facet_wrap(~ cell_type_label, ncol = ncol) +
      ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0),
                                   oob = scales::squish) +
      ggplot2::scale_fill_manual(values = color_map,
                                  breaks = c(non_grey_levels, "Other isoforms")) +
      ggplot2::labs(title = plot_title, fill = "Isoform") +
      ggplot2::theme_void(base_size = 13, base_family = "sans") +
      ggplot2::theme(
        plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16,
                                                  color = "#2c3e50",
                                                  margin = ggplot2::margin(b = 16)),
        strip.text       = ggplot2::element_text(size = 11, color = "#444444",
                                                   face = "bold", family = "sans",
                                                   margin = ggplot2::margin(b = 4, t = 10)),
        strip.background = ggplot2::element_blank(),
        panel.spacing.x  = ggplot2::unit(2.2, "lines"),
        panel.spacing.y  = ggplot2::unit(1.6, "lines"),
        plot.margin      = ggplot2::margin(16, 16, 16, 16),
        plot.background  = ggplot2::element_rect(fill = "white", color = NA)
      ) +
      legend_theme
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
