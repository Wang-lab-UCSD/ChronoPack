# ChronoPack dimensionality-reduction plotting helpers.
# Extracted from the consolidated puromycin DimensionalityReduction&Clustering notebook.

feature_plot_pc_simple <- function(
  pc_df,
  feature,
  x = "PC1",
  y = "PC2",
  palette = "viridis",
  title = NULL,
  colorbar_title = NULL,
  point_size = 1.2,      # scalar, or named vector mapping categorical feature levels -> sizes
  viridis_direction = 1,
  cap = NULL,
  floor = NULL,
  alpha_col = NULL,
  midpoint = NULL,
  mid_gamma = 1,
  plot_order = NULL,     # NULL, "plot_highest_last", "plot_nonzero_last", or "plot_highest_alpha_last"
  cat_colors = NULL,     # named vector for categorical feature colors and gene highlight colors
  cat_plot_order = NULL,
  base_size = 12,
  title_size = 12,
  axis_title_size = 11,
  axis_text_size = 10,
  legend_title_size = 10,
  legend_text_size = 9,
  show_legend = TRUE,
  show_ticks = TRUE,
  show_genes = NULL,         # character vector of genes to label
  gene_color_col = "broad_cluster",
  gene_label_style = c("box", "shadow"),
  gene_label_size = 3,
  gene_point_size = NULL,
  gene_label_color = "black",
  gene_label_nudge_x = NULL,
  gene_label_nudge_y = NULL,
  gene_label_fill = grDevices::adjustcolor("white", alpha.f = 0.88),
  gene_label_padding = 0.10,     # lines, for box labels
  gene_box_padding = 0.30,
  gene_point_padding = 0.15,
  gene_label_bg_col = "white",   # halo color for shadow mode
  gene_label_bg_r = 0.12,
  gene_label_segment_color = "grey45",
  gene_label_segment_size = 0.25,
  gene_label_seed = 1,
  check_overlap = FALSE          # only relevant for shadow mode
) {
  stopifnot(is.data.frame(pc_df))
  stopifnot(all(c(x, y, feature) %in% colnames(pc_df)))
  if (!is.null(alpha_col)) stopifnot(alpha_col %in% colnames(pc_df))

  if (!is.null(plot_order)) {
    stopifnot(plot_order %in% c(
      "plot_highest_last",
      "plot_nonzero_last",
      "plot_highest_alpha_last"
    ))
  }

    gene_label_style <- match.arg(gene_label_style)
    
    palette_input <- palette
    palette_key <- tolower(palette_input)
    
    # keep your existing shorthand aliases working
    brewer_aliases <- c(
      greys   = "Greys",
      blues   = "Blues",
      reds    = "Reds",
      greens  = "Greens",
      purples = "Purples",
      ylorrd  = "YlOrRd",
      ylgnbu  = "YlGnBu",
      gnbu    = "GnBu",
      pubu    = "PuBu",
      orrd    = "OrRd",
      bugn    = "BuGn"
    )
    
    brewer_name <- NULL
    
    if (palette_key %in% c("viridis", "magma", "plasma", "inferno", "cividis", "rwb")) {
      palette <- palette_key
    } else if (palette_key %in% names(brewer_aliases)) {
      palette <- palette_key
      brewer_name <- unname(brewer_aliases[palette_key])
    } else if (requireNamespace("RColorBrewer", quietly = TRUE) &&
               palette_input %in% rownames(RColorBrewer::brewer.pal.info)) {
      palette <- palette_input
      brewer_name <- palette_input
    } else if (requireNamespace("RColorBrewer", quietly = TRUE) &&
               palette_key %in% tolower(rownames(RColorBrewer::brewer.pal.info))) {
      brewer_name <- rownames(RColorBrewer::brewer.pal.info)[
        match(palette_key, tolower(rownames(RColorBrewer::brewer.pal.info)))
      ]
      palette <- brewer_name
    } else {
      stop(
        "Unsupported palette: ", palette_input,
        ". Use a viridis option, 'rwb', one of the existing aliases, ",
        "or any valid RColorBrewer palette name."
      )
    }

  if (!is.numeric(mid_gamma) || length(mid_gamma) != 1 || mid_gamma <= 0) {
    stop("mid_gamma must be a single positive number.")
  }

  if (is.null(title)) {
    title <- paste0("PCA (", feature, ")")
  }
  if (is.null(colorbar_title)) {
    colorbar_title <- feature
  }

  df_plot <- pc_df
  df_plot$.row_id <- if (!is.null(rownames(df_plot))) rownames(df_plot) else as.character(seq_len(nrow(df_plot)))

  # Parse rownames of form {gene}_{cell_type} using the last underscore
  df_plot$.parsed_gene <- sub("_(?=[^_]+$).*", "", df_plot$.row_id, perl = TRUE)
  df_plot$.parsed_gene[df_plot$.parsed_gene == ""] <- df_plot$.row_id

  v <- df_plot[[feature]]
  is_num <- is.numeric(v)

  if (is_num) {
    if (!is.null(cap))   v[v > cap]   <- cap
    if (!is.null(floor)) v[v < floor] <- floor
    df_plot[[feature]] <- v

    lo <- min(df_plot[[feature]], na.rm = TRUE)
    hi <- max(df_plot[[feature]], na.rm = TRUE)
  } else {
    df_plot[[feature]] <- as.factor(df_plot[[feature]])
    levs <- levels(df_plot[[feature]])

    if (!is.null(cat_plot_order)) {
      if (!is.character(cat_plot_order)) {
        stop("cat_plot_order must be a character vector of factor levels.")
      }

      unknown_levs <- setdiff(cat_plot_order, levs)
      if (length(unknown_levs) > 0) {
        stop(
          "cat_plot_order contains unknown levels: ",
          paste(unknown_levs, collapse = ", ")
        )
      }

      levs <- c(cat_plot_order, setdiff(levs, cat_plot_order))
      df_plot[[feature]] <- factor(df_plot[[feature]], levels = levs)
    }

    levs <- levels(df_plot[[feature]])

    if (!is.null(cat_colors)) {
      if (!is.character(cat_colors)) {
        stop("cat_colors must be a character vector of color values.")
      }

      if (!is.null(names(cat_colors)) && all(nzchar(names(cat_colors)))) {
        missing_levs <- setdiff(levs, names(cat_colors))
        if (length(missing_levs) > 0) {
          stop(
            "cat_colors is missing colors for these levels: ",
            paste(missing_levs, collapse = ", ")
          )
        }
        cat_colors <- cat_colors[levs]
      } else {
        if (length(cat_colors) < length(levs)) {
          stop(
            "cat_colors has fewer colors than the number of categorical levels (",
            length(levs), ")."
          )
        }
        cat_colors <- cat_colors[seq_along(levs)]
        names(cat_colors) <- levs
      }
    }
  }

  # Resolve point sizes
  scalar_point_size <- is.numeric(point_size) && length(point_size) == 1 && !is.na(point_size)

  if (scalar_point_size) {
    point_size_scalar <- point_size
    use_size_mapping <- FALSE
  } else {
    if (is_num) {
      stop("If feature is numeric, point_size must be a single numeric value.")
    }
    if (!is.numeric(point_size)) {
      stop("point_size must be either a single numeric value or a numeric vector.")
    }

    levs <- levels(df_plot[[feature]])

    if (!is.null(names(point_size)) && all(nzchar(names(point_size)))) {
      missing_levs <- setdiff(levs, names(point_size))
      if (length(missing_levs) > 0) {
        stop(
          "point_size is missing sizes for these levels: ",
          paste(missing_levs, collapse = ", ")
        )
      }
      point_size_map <- point_size[levs]
    } else {
      if (length(point_size) < length(levs)) {
        stop(
          "point_size has fewer values than the number of categorical levels (",
          length(levs), ")."
        )
      }
      point_size_map <- point_size[seq_along(levs)]
      names(point_size_map) <- levs
    }

    df_plot$.point_size <- unname(point_size_map[as.character(df_plot[[feature]])])
    use_size_mapping <- TRUE
  }

  # Plotting order
  if (!is.null(plot_order)) {
    if (plot_order == "plot_highest_alpha_last") {
      if (is.null(alpha_col)) {
        stop("plot_order = 'plot_highest_alpha_last' requires alpha_col to be provided.")
      }
      ord <- order(df_plot[[alpha_col]], na.last = TRUE)
      df_plot <- df_plot[ord, , drop = FALSE]

    } else {
      if (!is_num) {
        warning("plot_order based on feature magnitude is only applied for numeric features.")
      } else {
        if (plot_order == "plot_highest_last") {
          ord <- order(df_plot[[feature]], na.last = TRUE)
        } else if (plot_order == "plot_nonzero_last") {
          ord <- order(abs(df_plot[[feature]]), na.last = TRUE)
        }
        df_plot <- df_plot[ord, , drop = FALSE]
      }
    }
  }

  if (!is_num && !is.null(cat_plot_order)) {
    ord <- order(match(as.character(df_plot[[feature]]), levels(df_plot[[feature]])), na.last = TRUE)
    df_plot <- df_plot[ord, , drop = FALSE]
  }

  if (is.null(alpha_col) && !use_size_mapping) {
    aes_map <- ggplot2::aes_string(x = x, y = y, color = feature)
  } else if (!is.null(alpha_col) && !use_size_mapping) {
    aes_map <- ggplot2::aes_string(x = x, y = y, color = feature, alpha = alpha_col)
  } else if (is.null(alpha_col) && use_size_mapping) {
    aes_map <- ggplot2::aes_string(x = x, y = y, color = feature, size = ".point_size")
  } else {
    aes_map <- ggplot2::aes_string(x = x, y = y, color = feature, alpha = alpha_col, size = ".point_size")
  }

  p <- ggplot2::ggplot(df_plot, aes_map)

  if (use_size_mapping) {
    p <- p + ggplot2::geom_point(shape = 16)
  } else {
    p <- p + ggplot2::geom_point(size = point_size_scalar, shape = 16)
  }

  p <- p +
    ggplot2::labs(
      title = title,
      x = x,
      y = y,
      color = colorbar_title
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(-10, 0, 10),
      expand = ggplot2::expansion(mult = c(0.06, 0.06))
    ) +
    ggplot2::scale_y_continuous(
      breaks = c(-10, 0, 10),
      expand = ggplot2::expansion(mult = c(0.06, 0.06))
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = if (show_legend) "right" else "none",
      plot.title   = ggplot2::element_text(size = title_size),
      axis.title.x = ggplot2::element_text(size = axis_title_size),
      axis.title.y = ggplot2::element_text(size = axis_title_size),
      axis.text.x  = if (show_ticks) ggplot2::element_text(size = axis_text_size) else ggplot2::element_blank(),
      axis.text.y  = if (show_ticks) ggplot2::element_text(size = axis_text_size) else ggplot2::element_blank(),
      axis.ticks.x = if (show_ticks) ggplot2::element_line() else ggplot2::element_blank(),
      axis.ticks.y = if (show_ticks) ggplot2::element_line() else ggplot2::element_blank(),
      legend.title = ggplot2::element_text(size = legend_title_size),
      legend.text  = ggplot2::element_text(size = legend_text_size),
      plot.margin  = ggplot2::margin(8, 14, 8, 14)
    ) +
    ggplot2::coord_cartesian(clip = "off")

  if (!is.null(alpha_col)) {
    p <- p + ggplot2::scale_alpha_identity()
  }

  if (use_size_mapping) {
    p <- p + ggplot2::scale_size_identity()
  }

  if (!is_num) {
    if (is.null(cat_colors)) {
      p <- p + ggplot2::scale_color_discrete(drop = FALSE)
    } else {
      p <- p + ggplot2::scale_color_manual(values = cat_colors, drop = FALSE)
    }
  } else if (palette == "rwb") {
    if (is.null(midpoint)) midpoint <- 0

    vals_raw <- seq(lo, hi, length.out = 256)
    s <- numeric(length(vals_raw))

    left <- vals_raw <= midpoint
    right <- vals_raw >= midpoint

    if (any(left)) {
      denom_left <- max(midpoint - lo, .Machine$double.eps)
      s[left] <- -((midpoint - vals_raw[left]) / denom_left)
    }
    if (any(right)) {
      denom_right <- max(hi - midpoint, .Machine$double.eps)
      s[right] <- ((vals_raw[right] - midpoint) / denom_right)
    }

    s2 <- sign(s) * abs(s)^mid_gamma
    vals_scaled <- (s2 + 1) / 2

    p <- p + ggplot2::scale_color_gradientn(
      colours = grDevices::colorRampPalette(c("blue", "white", "red"))(256),
      values = vals_scaled,
      limits = c(lo, hi),
      guide = ggplot2::guide_colorbar(title = colorbar_title)
    )
    } else if (!is.null(brewer_name)) {
      p <- p + ggplot2::scale_color_distiller(
        palette = brewer_name,
        direction = 1,
        limits = c(lo, hi),
        guide = ggplot2::guide_colorbar(title = colorbar_title))
      } else {
    p <- p + ggplot2::scale_color_viridis_c(
      option = palette,
      direction = viridis_direction,
      limits = c(lo, hi),
      guide = ggplot2::guide_colorbar(title = colorbar_title)
    )
  }

  if (!show_legend) {
    p <- p + ggplot2::guides(color = "none")
  }

  # Gene highlighting and labels
  if (!is.null(show_genes)) {
    if (!is.character(show_genes)) {
      stop("show_genes must be NULL or a character vector.")
    }

    label_match <- df_plot$.row_id %in% show_genes | df_plot$.parsed_gene %in% show_genes
    label_text <- df_plot$.parsed_gene

    if ("gene_sym" %in% colnames(df_plot)) {
      gene_sym_chr <- as.character(df_plot$gene_sym)
      label_match <- label_match | gene_sym_chr %in% show_genes
      label_text <- ifelse(!is.na(gene_sym_chr) & nzchar(gene_sym_chr), gene_sym_chr, df_plot$.parsed_gene)
    }

    label_df <- df_plot[label_match, , drop = FALSE]
    label_df$.label <- label_text[label_match]

    if (nrow(label_df) > 0) {
      # highlighted dot colors come from category colors if available
      label_df$.dot_col <- rep("black", nrow(label_df))

      if (!is.null(cat_colors) &&
          !is.null(gene_color_col) &&
          gene_color_col %in% colnames(label_df) &&
          !is.null(names(cat_colors))) {
        gene_cat <- as.character(label_df[[gene_color_col]])
        ok <- !is.na(gene_cat) & gene_cat %in% names(cat_colors)
        label_df$.dot_col[ok] <- unname(cat_colors[gene_cat[ok]])
      }

      if (is.null(gene_point_size)) {
        if (scalar_point_size) {
          gene_point_size_use <- point_size_scalar * 1.8
        } else {
          gene_point_size_use <- max(df_plot$.point_size, na.rm = TRUE) * 1.35
        }
      } else {
        gene_point_size_use <- gene_point_size
      }

      # larger colored highlight dots
      p <- p + ggplot2::geom_point(
        data = label_df,
        mapping = ggplot2::aes_string(x = x, y = y),
        inherit.aes = FALSE,
        size = gene_point_size_use,
        shape = 16,
        colour = label_df$.dot_col
      )

      if (gene_label_style == "box") {
        if (!requireNamespace("ggrepel", quietly = TRUE)) {
          stop("Please install the ggrepel package to use gene_label_style = 'box'.")
        }

        p <- p + ggrepel::geom_label_repel(
          data = label_df,
          mapping = ggplot2::aes_string(x = x, y = y, label = ".label"),
          inherit.aes = FALSE,
          colour = gene_label_color,
          fill = gene_label_fill,
          label.size = 0,
          label.padding = grid::unit(gene_label_padding, "lines"),
          box.padding = gene_box_padding,
          point.padding = gene_point_padding,
          force = 1,
          force_pull = 0.2,
          max.overlaps = Inf,
          min.segment.length = 0,
          seed = gene_label_seed,
          segment.size = gene_label_segment_size,
          segment.color = gene_label_segment_color,
          size = gene_label_size / ggplot2::.pt
        )

      } else if (gene_label_style == "shadow") {
        if (!requireNamespace("shadowtext", quietly = TRUE)) {
          stop("Please install the shadowtext package to use gene_label_style = 'shadow'.")
        }

        xr <- range(df_plot[[x]], na.rm = TRUE)
        yr <- range(df_plot[[y]], na.rm = TRUE)
        xmid <- mean(xr)
        ymid <- mean(yr)

        dx_auto <- 0.04 * diff(xr)
        dy_auto <- 0.04 * diff(yr)
        pad_x <- 0.05 * diff(xr)
        pad_y <- 0.05 * diff(yr)

        dx_use <- if (is.null(gene_label_nudge_x)) dx_auto else gene_label_nudge_x
        dy_use <- if (is.null(gene_label_nudge_y)) dy_auto else gene_label_nudge_y

        label_df$.x_lab <- label_df[[x]] + ifelse(label_df[[x]] <= xmid, dx_use, -dx_use)
        label_df$.y_lab <- label_df[[y]] + ifelse(label_df[[y]] <= ymid, dy_use, -dy_use)

        # clamp inward so labels do not collide with axes / panel edge
        label_df$.x_lab <- pmin(pmax(label_df$.x_lab, xr[1] + pad_x), xr[2] - pad_x)
        label_df$.y_lab <- pmin(pmax(label_df$.y_lab, yr[1] + pad_y), yr[2] - pad_y)

        p <- p +
          ggplot2::geom_segment(
            data = label_df,
            mapping = ggplot2::aes_string(x = x, y = y, xend = ".x_lab", yend = ".y_lab"),
            inherit.aes = FALSE,
            colour = gene_label_segment_color,
            linewidth = gene_label_segment_size
          ) +
        shadowtext::geom_shadowtext(
          data = label_df,
          mapping = ggplot2::aes(x = .x_lab, y = .y_lab, label = .label),
          inherit.aes = FALSE,
          colour = gene_label_color,
          bg.colour = gene_label_bg_col,
          bg.r = gene_label_bg_r,
          size = gene_label_size / ggplot2::.pt,
          fontface = "plain",
          check_overlap = check_overlap
        )
      }
    }
  }

  p
}

compute_tradeSeq_reference_offset <- function(sizeFactors, off_ref = NULL) {
  if (!is.null(off_ref)) return(as.numeric(off_ref)[1])
  sf_vec <- as.numeric(sizeFactors)
  if (!length(sf_vec) || all(!is.finite(sf_vec))) {
    stop("sizeFactors must contain at least one finite value.")
  }
  log(exp(mean(log(sf_vec), na.rm = TRUE)))
}

get_tradeSeq_sce_prediction_components <- function(sce, gene) {
  rd_ts <- S4Vectors::metadata(SummarizedExperiment::rowData(sce))$tradeSeq
  if (is.null(rd_ts)) {
    rd_ts <- SummarizedExperiment::rowData(sce)$tradeSeq
  }
  cd_ts <- SummarizedExperiment::colData(sce)$tradeSeq
  if (is.null(rd_ts) || is.null(cd_ts)) {
    stop("SCE does not contain tradeSeq rowData/colData metadata.")
  }

  gene <- as.character(gene)[1]
  gene_idx <- match(gene, rownames(sce))
  if (is.na(gene_idx)) {
    stop("Gene not present in SCE rownames: ", gene)
  }

  beta_tbl <- rd_ts$beta[[1]]
  if (is.null(beta_tbl)) {
    stop("Could not locate tradeSeq beta matrix in rowData(sce)$tradeSeq$beta[[1]].")
  }
  beta_mat <- as.matrix(beta_tbl)
  beta_idx <- gene_idx
  if (!is.null(rownames(beta_mat)) && gene %in% rownames(beta_mat)) {
    beta_idx <- match(gene, rownames(beta_mat))
  }
  beta <- as.numeric(beta_mat[beta_idx, ])
  names(beta) <- colnames(beta_mat)

  Sigma <- rd_ts$Sigma[[gene_idx]]
  if (is.null(Sigma)) {
    stop("Could not locate tradeSeq covariance matrix for gene: ", gene)
  }
  Sigma <- as.matrix(Sigma)

  list(beta = beta, Sigma = Sigma, dm = cd_ts$dm, X = cd_ts$X)
}

predict_tradeSeq_sce_with_se <- function(sce, gene, lineageId = 1, nPoints = 200,
                                         off_ref = 0, se_factor = 1.96,
                                         return_log1p = TRUE) {
  comp <- get_tradeSeq_sce_prediction_components(sce, gene)
  dm <- comp$dm
  X <- comp$X
  if (is.null(dm) || is.null(X)) {
    stop("SCE is missing tradeSeq design matrix metadata.")
  }

  t_cols <- grep("^t[0-9]+$", colnames(dm), value = TRUE)
  avail_lineages <- as.integer(sub("^t", "", t_cols))
  if (!length(avail_lineages)) {
    stop("Could not find t1/t2/... columns in tradeSeq design metadata.")
  }
  ell <- as.integer(lineageId)
  if (!(ell %in% avail_lineages)) ell <- avail_lineages[1]

  df_pred <- tradeSeq:::.getPredictRangeDf(dm, lineageId = ell, nPoints = nPoints)
  slingshot_coldata <- SummarizedExperiment::colData(sce)$crv
  if (!is.null(slingshot_coldata)) {
    pseudotime <- slingshot_coldata[, grep(x = colnames(slingshot_coldata), pattern = "pseudotime"), drop = FALSE]
  } else {
    pseudotime <- as.matrix(dm[, t_cols, drop = FALSE])
  }
  X_pred <- tradeSeq:::predictGAM(lpmatrix = X, df = df_pred, pseudotime = pseudotime)

  beta <- comp$beta
  Sigma <- comp$Sigma
  if (!is.null(names(beta)) && !is.null(colnames(X_pred)) && all(colnames(X_pred) %in% names(beta))) {
    beta <- beta[colnames(X_pred)]
  }
  if (!is.null(rownames(Sigma)) && !is.null(colnames(X_pred)) && all(colnames(X_pred) %in% rownames(Sigma))) {
    Sigma <- Sigma[colnames(X_pred), colnames(X_pred), drop = FALSE]
  }
  if (length(beta) != ncol(X_pred)) {
    stop("Beta length does not match prediction design columns for gene: ", gene)
  }
  if (!all(dim(Sigma) == c(ncol(X_pred), ncol(X_pred)))) {
    stop("Covariance dimensions do not match prediction design columns for gene: ", gene)
  }

  eta <- as.numeric(X_pred %*% beta + as.numeric(off_ref)[1])
  se_eta <- sqrt(pmax(0, rowSums((X_pred %*% Sigma) * X_pred)))
  out <- data.frame(
    gene = as.character(gene)[1],
    lineage = ell,
    pt = df_pred[[paste0("t", ell)]],
    fit = eta,
    se.fit = se_eta,
    stringsAsFactors = FALSE
  )
  if (isTRUE(return_log1p)) {
    out$mu <- log1p(exp(eta))
    out$lo <- log1p(exp(eta - se_factor * se_eta))
    out$hi <- log1p(exp(eta + se_factor * se_eta))
  } else {
    out$mu <- exp(eta)
    out$lo <- exp(eta - se_factor * se_eta)
    out$hi <- exp(eta + se_factor * se_eta)
  }
  out
}

plot_tradeSeq_gene <- function(
    gene, sce_ls, q_by_gene = NULL,
    keys = NULL, key_labels = NULL,
    offset_map = NULL, time_max_map = NULL,
    plot_normalized_pt = FALSE, trim_before_zero = FALSE,
    point_tmax_map = NULL,
    time_unit = "min", x_tick_step = 20,
    x_tick_labels = NULL, xlim = NULL,
    omit_points = FALSE, nPoints = 1000,
    ny_ticks = 5, ylim = NULL,
    color_points_by_run = TRUE,
    points_color = "#274C77", points_alpha = 0.30, points_size = 0.75,
    se_factor = 1.96,
    line_width = 1.7, ribbon_alpha = 0.35,
    show_title = TRUE, show_legend = TRUE,
    axis_text_size = 14, axis_title_size = 16,
    colours = NULL, ribbon_alpha_mult = 0.30
) {
    if (is.null(keys)) keys <- names(sce_ls)
    keys <- intersect(keys, names(sce_ls))
    keys <- keys[vapply(keys, function(k) gene %in% rownames(sce_ls[[k]]), logical(1))]
    if (!length(keys)) stop("Gene not found in any requested SCE.")

    if (is.null(key_labels)) key_labels <- setNames(keys, keys)
    labels <- unname(unlist(key_labels[keys]))

    if (is.null(colours)) {
        colours <- rep(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a",
                         "#66a61e", "#e6ab02", "#a6761d", "#666666"),
                       length.out = length(keys))
        names(colours) <- keys
    } else if (is.null(names(colours))) {
        colours <- setNames(rep(colours, length.out = length(keys)), keys)
    } else {
        colours <- colours[keys]
        if (any(is.na(colours))) stop("Named `colours` must include all requested keys.")
    }
    fill_cols <- setNames(grDevices::adjustcolor(colours, alpha.f = ribbon_alpha_mult), keys)

    map_num <- function(x, k, default = NA_real_) {
        if (is.null(x) || is.null(names(x)) || !(k %in% names(x))) return(default)
        v <- x[[k]]
        if (!length(v) || is.na(v)) default else as.numeric(v)
    }
    get_time_max <- function(k) {
        if (plot_normalized_pt) return(1)
        v <- map_num(time_max_map, k)
        if (is.na(v)) stop("`time_max_map` must provide a value for key: ", k)
        v
    }
    get_offset <- function(k) map_num(offset_map, k, 0)
    pt_to_x <- function(pt, tm, off) (if (plot_normalized_pt) pt else pt * tm) + off
    wmean <- function(x, w) {
        s <- sum(w, na.rm = TRUE)
        if (is.finite(s) && s > 0) sum(x * w, na.rm = TRUE) / s else NA_real_
    }

    fit_list <- obs_list <- setNames(vector("list", length(keys)), keys)

    for (k in keys) {
        sce <- sce_ls[[k]]
        cd <- SummarizedExperiment::colData(sce)
        rd <- SummarizedExperiment::rowData(sce)
        crv <- cd$crv
        ts_cd <- cd$tradeSeq
        ts_rd <- S4Vectors::metadata(rd)$tradeSeq
        if (is.null(ts_rd)) ts_rd <- rd$tradeSeq
        if (is.null(crv)) stop("[", k, "] Missing slingshot curve metadata in colData(sce)$crv.")
        if (is.null(ts_cd) || is.null(ts_rd)) stop("[", k, "] Missing tradeSeq metadata in SCE.")

        dm <- ts_cd$dm
        if (is.null(dm)) stop("[", k, "] Missing tradeSeq design matrix metadata.")

        t_cols <- grep("^t[0-9]+$", colnames(dm), value = TRUE)
        if (!length(t_cols)) stop("[", k, "] No lineage columns found in tradeSeq design matrix.")
        ell <- if ("t1" %in% t_cols) 1L else as.integer(sub("^t", "", t_cols[1]))

        pt_col <- paste0("pseudotime.lineage", ell)
        cw_col <- paste0("cellWeights.lineage", ell)
        if (!(pt_col %in% colnames(crv))) stop("[", k, "] Missing ", pt_col)

        pt <- as.numeric(crv[[pt_col]])
        w <- if (cw_col %in% colnames(crv)) as.numeric(crv[[cw_col]]) else rep(1, ncol(sce))
        w[!is.finite(w)] <- 0

        offset_col <- grep("offset", colnames(dm), value = TRUE)[1]
        if (is.na(offset_col)) stop("[", k, "] No offset column found in tradeSeq design matrix.")

        offset_vec <- as.numeric(dm[[offset_col]])
        off_ref <- mean(offset_vec, na.rm = TRUE)

        df_pred <- tradeSeq:::.getPredictRangeDf(dm, lineageId = ell, nPoints = nPoints)
        t_pt <- df_pred[[paste0("t", ell)]]

        pt_cols <- grep("^pseudotime\\.lineage", colnames(crv), value = TRUE)
        pt_mat <- as.matrix(crv[, pt_cols, drop = FALSE])

        X_pred <- tradeSeq:::predictGAM(
            lpmatrix = ts_cd$X,
            df = df_pred,
            pseudotime = pt_mat
        )

        gene_i <- match(gene, rownames(sce))
        beta_mat <- as.matrix(ts_rd$beta[[1]])
        beta_i <- gene_i
        if (!is.null(rownames(beta_mat)) && gene %in% rownames(beta_mat)) {
            beta_i <- match(gene, rownames(beta_mat))
        }
        beta <- as.numeric(beta_mat[beta_i, , drop = TRUE])
        names(beta) <- colnames(beta_mat)
        Sigma <- as.matrix(ts_rd$Sigma[[gene_i]])

        if (!is.null(names(beta)) && !is.null(colnames(X_pred)) && all(colnames(X_pred) %in% names(beta))) {
            beta <- beta[colnames(X_pred)]
        }
        if (!is.null(rownames(Sigma)) && !is.null(colnames(X_pred)) && all(colnames(X_pred) %in% rownames(Sigma))) {
            Sigma <- Sigma[colnames(X_pred), colnames(X_pred), drop = FALSE]
        }
        if (length(beta) != ncol(X_pred)) {
            stop("[", k, "] Beta length does not match prediction design columns for gene: ", gene)
        }
        if (!all(dim(Sigma) == c(ncol(X_pred), ncol(X_pred)))) {
            stop("[", k, "] Covariance dimensions do not match prediction design columns for gene: ", gene)
        }

        eta <- drop(X_pred %*% beta) + off_ref
        se <- sqrt(pmax(0, rowSums((X_pred %*% Sigma) * X_pred)))

        tm <- get_time_max(k)
        off <- get_offset(k)
        run <- factor(unname(key_labels[[k]]), levels = labels)

        fit_list[[k]] <- data.frame(
            key = k, run = run, pt = t_pt, x = pt_to_x(t_pt, tm, off),
            mu = log1p(exp(eta)),
            lo = log1p(exp(eta - se_factor * se)),
            hi = log1p(exp(eta + se_factor * se))
        )

        counts <- SummarizedExperiment::assay(sce, "counts")
        y <- as.numeric(counts[gene, ]) * exp(off_ref - offset_vec)

        obs_df <- data.frame(
            key = k, run = run, pt = pt,
            x = pt_to_x(pt, tm, off), y = log1p(y)
        )

        nb <- min(nPoints, 200)
        bin <- cut(pt, breaks = seq(0, 1, length.out = nb + 1),
                   include.lowest = TRUE, right = TRUE)

        pt_mid <- tapply(pt, bin, mean, na.rm = TRUE)
        y_mean <- tapply(seq_along(y), bin, function(ii) wmean(y[ii], w[ii]))

        mean_df <- data.frame(
            key = k, run = run,
            pt_mid = as.numeric(pt_mid),
            x_mid = pt_to_x(as.numeric(pt_mid), tm, off),
            mean_norm_log1p = log1p(as.numeric(y_mean))
        )
        mean_df <- mean_df[
            is.finite(mean_df$pt_mid) & is.finite(mean_df$mean_norm_log1p), ,
            drop = FALSE
        ]

        obs_list[[k]] <- list(points = obs_df, means = mean_df)
    }

    fit_all <- do.call(rbind, fit_list)
    obs_pts <- do.call(rbind, lapply(obs_list, `[[`, "points"))
    obs_mea <- do.call(rbind, lapply(obs_list, `[[`, "means"))

    if (trim_before_zero) {
        fit_all <- fit_all[fit_all$x >= 0, , drop = FALSE]
        obs_pts <- obs_pts[obs_pts$x >= 0, , drop = FALSE]
        obs_mea <- obs_mea[obs_mea$x_mid >= 0, , drop = FALSE]
    }

    if (!is.null(point_tmax_map)) {
        cutoff <- vapply(as.character(obs_pts$key),
                         function(k) map_num(point_tmax_map, k), numeric(1))
        obs_pts <- obs_pts[is.na(cutoff) | obs_pts$x <= cutoff, , drop = FALSE]

        cutoff <- vapply(as.character(obs_mea$key),
                         function(k) map_num(point_tmax_map, k), numeric(1))
        obs_mea <- obs_mea[is.na(cutoff) | obs_mea$x_mid <= cutoff, , drop = FALSE]
    }

    sub_txt <- NULL
    if (!is.null(q_by_gene) && gene %in% names(q_by_gene)) {
        qv <- unname(q_by_gene[gene])
        if (is.finite(qv)) sub_txt <- paste0("FDR = ", formatC(qv, format = "fg", digits = 3))
    }

    if (is.null(ylim)) {
        ylo <- min(c(fit_all$lo, obs_mea$mean_norm_log1p,
                     if (!omit_points) obs_pts$y), na.rm = TRUE)
        yhi <- max(c(fit_all$hi, obs_mea$mean_norm_log1p,
                     if (!omit_points) obs_pts$y), na.rm = TRUE)

        if (is.finite(ylo) && is.finite(yhi)) {
            pad <- 0.04 * (yhi - ylo)
            if (!is.finite(pad) || pad == 0) pad <- 0.2
            ylim <- c(max(0, ylo - pad), 1.3 * (yhi + pad))
        } else {
            ylim <- c(0, 1)
        }
    }

    ylim[2] <- 1.4 * max(fit_all$mu, na.rm = TRUE)

    if (plot_normalized_pt) {
        x_title <- "Pseudotime"
        x_breaks <- x_labels <- seq(0, 1, by = 0.2)
    } else {
        x_title <- paste0("Time", if (!is.null(time_unit)) paste0(" (", time_unit, ")") else "")
        if (is.null(x_tick_labels)) {
            tmax <- max(vapply(keys, get_time_max, numeric(1)), na.rm = TRUE)
            x_tick_labels <- seq(0, tmax, by = x_tick_step)
        }
        x_breaks <- x_labels <- x_tick_labels
    }

    if (is.null(xlim)) {
        xlim <- range(c(fit_all$x, obs_mea$x_mid,
                        if (!omit_points) obs_pts$x), na.rm = TRUE)
    }

    oob_fun <- if (requireNamespace("scales", quietly = TRUE))
        scales::squish else function(x, ...) x

    p <- ggplot2::ggplot()

    if (!omit_points) {
        if (color_points_by_run) {
            p <- p + ggplot2::geom_point(
                data = obs_pts, ggplot2::aes(x, y, colour = run),
                size = points_size, alpha = points_alpha
            )
        } else {
            p <- p + ggplot2::geom_point(
                data = obs_pts, ggplot2::aes(x, y),
                size = points_size, alpha = points_alpha, colour = points_color
            )
        }
    }

    p +
        ggplot2::geom_ribbon(
            data = fit_all, ggplot2::aes(x, ymin = lo, ymax = hi, fill = run),
            alpha = ribbon_alpha, colour = NA
        ) +
        ggplot2::geom_line(
            data = fit_all, ggplot2::aes(x, mu, colour = run),
            linewidth = line_width
        ) +
        ggplot2::geom_point(
            data = obs_mea, ggplot2::aes(x_mid, mean_norm_log1p, colour = run),
            shape = 21, fill = "white", size = 3.2, stroke = 0.9
        ) +
        ggplot2::scale_colour_manual(
            values = setNames(colours[keys], labels), name = "Run"
        ) +
        ggplot2::scale_fill_manual(
            values = setNames(fill_cols[keys], labels), name = "Run"
        ) +
        ggplot2::scale_y_continuous(
            breaks = if (requireNamespace("scales", quietly = TRUE))
                scales::pretty_breaks(n = ny_ticks) else ggplot2::waiver(),
            name = "Expression"
        ) +
        ggplot2::scale_x_continuous(
            breaks = x_breaks, labels = x_labels, limits = xlim, oob = oob_fun
        ) +
        ggplot2::coord_cartesian(expand = FALSE, ylim = ylim) +
        ggplot2::labs(
            title = if (show_title) gene else NULL,
            subtitle = if (show_title) sub_txt else NULL,
            x = x_title, y = "Expression"
        ) +
        ggplot2::theme_minimal(base_size = axis_text_size) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_line(colour = "grey88"),
            panel.grid.major.y = ggplot2::element_line(colour = "grey88"),
            axis.text = ggplot2::element_text(colour = "black", size = axis_text_size),
            axis.title = ggplot2::element_text(colour = "black", size = axis_title_size),
            plot.title = ggplot2::element_text(face = "bold", hjust = 0, colour = "black"),
            plot.subtitle = ggplot2::element_text(size = axis_text_size, colour = "black"),
            legend.position = if (show_legend) "right" else "none"
        )
}

# TESTING.ipynb's SCE-only multi-run gene plot (its plot_tradeSeq_gene, cell 3), verbatim
# under a new name; plot_tradeSeq_gene above is the same plot with extra metadata checks.
plot_tradeSeq_gene_sce_pub <- function(
    gene, sce_ls, q_by_gene = NULL,
    keys = NULL, key_labels = NULL,
    offset_map = NULL, time_max_map = NULL,
    plot_normalized_pt = FALSE, trim_before_zero = FALSE,
    point_tmax_map = NULL,
    time_unit = "min", x_tick_step = 20,
    x_tick_labels = NULL, xlim = NULL,
    omit_points = FALSE, nPoints = 1000,
    ny_ticks = 5, ylim = NULL,
    color_points_by_run = TRUE,
    points_color = "#274C77", points_alpha = 0.30, points_size = 0.75,
    se_factor = 1.96,
    line_width = 1.7, ribbon_alpha = 0.35,
    show_title = TRUE, show_legend = TRUE,
    axis_text_size = 14, axis_title_size = 16,
    colours = NULL, ribbon_alpha_mult = 0.30
) {
    if (is.null(keys)) keys <- names(sce_ls)
    keys <- intersect(keys, names(sce_ls))
    keys <- keys[vapply(keys, function(k) gene %in% rownames(sce_ls[[k]]), logical(1))]
    if (!length(keys)) stop("Gene not found in any requested SCE.")

    if (is.null(key_labels)) key_labels <- setNames(keys, keys)
    labels <- unname(unlist(key_labels[keys]))

    if (is.null(colours)) {
        colours <- rep(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a",
                         "#66a61e", "#e6ab02", "#a6761d", "#666666"),
                       length.out = length(keys))
        names(colours) <- keys
    } else if (is.null(names(colours))) {
        colours <- setNames(rep(colours, length.out = length(keys)), keys)
    } else {
        colours <- colours[keys]
        if (any(is.na(colours))) stop("Named `colours` must include all requested keys.")
    }
    fill_cols <- setNames(grDevices::adjustcolor(colours, alpha.f = ribbon_alpha_mult), keys)

    map_num <- function(x, k, default = NA_real_) {
        if (is.null(x) || is.null(names(x)) || !(k %in% names(x))) return(default)
        v <- x[[k]]
        if (!length(v) || is.na(v)) default else as.numeric(v)
    }
    get_time_max <- function(k) {
        if (plot_normalized_pt) return(1)
        v <- map_num(time_max_map, k)
        if (is.na(v)) stop("`time_max_map` must provide a value for key: ", k)
        v
    }
    get_offset <- function(k) map_num(offset_map, k, 0)
    pt_to_x <- function(pt, tm, off) (if (plot_normalized_pt) pt else pt * tm) + off
    wmean <- function(x, w) {
        s <- sum(w, na.rm = TRUE)
        if (is.finite(s) && s > 0) sum(x * w, na.rm = TRUE) / s else NA_real_
    }

    fit_list <- obs_list <- setNames(vector("list", length(keys)), keys)

    for (k in keys) {
        sce <- sce_ls[[k]]
        cd <- SummarizedExperiment::colData(sce)
        rd <- SummarizedExperiment::rowData(sce)
        crv <- cd$crv
        ts_cd <- cd$tradeSeq
        ts_rd <- rd$tradeSeq
        dm <- ts_cd$dm

        t_cols <- grep("^t[0-9]+$", colnames(dm), value = TRUE)
        if (!length(t_cols)) stop("[", k, "] No lineage columns found in tradeSeq design matrix.")
        ell <- if ("t1" %in% t_cols) 1L else as.integer(sub("^t", "", t_cols[1]))

        pt_col <- paste0("pseudotime.lineage", ell)
        cw_col <- paste0("cellWeights.lineage", ell)
        if (!(pt_col %in% colnames(crv))) stop("[", k, "] Missing ", pt_col)

        pt <- as.numeric(crv[[pt_col]])
        w <- if (cw_col %in% colnames(crv)) as.numeric(crv[[cw_col]]) else rep(1, ncol(sce))
        w[!is.finite(w)] <- 0

        offset_col <- grep("offset", colnames(dm), value = TRUE)[1]
        if (is.na(offset_col)) stop("[", k, "] No offset column found in tradeSeq design matrix.")

        offset_vec <- as.numeric(dm[[offset_col]])
        off_ref <- mean(offset_vec, na.rm = TRUE)

        df_pred <- tradeSeq:::.getPredictRangeDf(dm, lineageId = ell, nPoints = nPoints)
        t_pt <- df_pred[[paste0("t", ell)]]

        pt_cols <- grep("^pseudotime\\.lineage", colnames(crv), value = TRUE)
        pt_mat <- as.matrix(crv[, pt_cols, drop = FALSE])

        X_pred <- tradeSeq:::predictGAM(
            lpmatrix = ts_cd$X,
            df = df_pred,
            pseudotime = pt_mat
        )

        gene_i <- match(gene, rownames(sce))
        beta_mat <- ts_rd$beta[[1]]
        beta <- as.numeric(beta_mat[gene_i, , drop = TRUE])
        Sigma <- ts_rd$Sigma[[gene_i]]

        eta <- drop(X_pred %*% beta) + off_ref
        se <- sqrt(pmax(0, rowSums((X_pred %*% Sigma) * X_pred)))

        tm <- get_time_max(k)
        off <- get_offset(k)
        run <- factor(unname(key_labels[[k]]), levels = labels)

        fit_list[[k]] <- data.frame(
            key = k, run = run, pt = t_pt, x = pt_to_x(t_pt, tm, off),
            mu = log1p(exp(eta)),
            lo = log1p(exp(eta - se_factor * se)),
            hi = log1p(exp(eta + se_factor * se))
        )

        counts <- SummarizedExperiment::assay(sce, "counts")
        y <- as.numeric(counts[gene, ]) * exp(off_ref - offset_vec)

        obs_df <- data.frame(
            key = k, run = run, pt = pt,
            x = pt_to_x(pt, tm, off), y = log1p(y)
        )

        nb <- min(nPoints, 200)
        bin <- cut(pt, breaks = seq(0, 1, length.out = nb + 1),
                   include.lowest = TRUE, right = TRUE)

        pt_mid <- tapply(pt, bin, mean, na.rm = TRUE)
        y_mean <- tapply(seq_along(y), bin, function(ii) wmean(y[ii], w[ii]))

        mean_df <- data.frame(
            key = k, run = run,
            pt_mid = as.numeric(pt_mid),
            x_mid = pt_to_x(as.numeric(pt_mid), tm, off),
            mean_norm_log1p = log1p(as.numeric(y_mean))
        )
        mean_df <- mean_df[
            is.finite(mean_df$pt_mid) & is.finite(mean_df$mean_norm_log1p), ,
            drop = FALSE
        ]

        obs_list[[k]] <- list(points = obs_df, means = mean_df)
    }

    fit_all <- do.call(rbind, fit_list)
    obs_pts <- do.call(rbind, lapply(obs_list, `[[`, "points"))
    obs_mea <- do.call(rbind, lapply(obs_list, `[[`, "means"))

    if (trim_before_zero) {
        fit_all <- fit_all[fit_all$x >= 0, , drop = FALSE]
        obs_pts <- obs_pts[obs_pts$x >= 0, , drop = FALSE]
        obs_mea <- obs_mea[obs_mea$x_mid >= 0, , drop = FALSE]
    }

    if (!is.null(point_tmax_map)) {
        cutoff <- vapply(as.character(obs_pts$key),
                         function(k) map_num(point_tmax_map, k), numeric(1))
        obs_pts <- obs_pts[is.na(cutoff) | obs_pts$x <= cutoff, , drop = FALSE]

        cutoff <- vapply(as.character(obs_mea$key),
                         function(k) map_num(point_tmax_map, k), numeric(1))
        obs_mea <- obs_mea[is.na(cutoff) | obs_mea$x_mid <= cutoff, , drop = FALSE]
    }

    sub_txt <- NULL
    if (!is.null(q_by_gene) && gene %in% names(q_by_gene)) {
        qv <- unname(q_by_gene[gene])
        if (is.finite(qv)) sub_txt <- paste0("FDR = ", formatC(qv, format = "fg", digits = 3))
    }

    if (is.null(ylim)) {
        ylo <- min(c(fit_all$lo, obs_mea$mean_norm_log1p,
                     if (!omit_points) obs_pts$y), na.rm = TRUE)
        yhi <- max(c(fit_all$hi, obs_mea$mean_norm_log1p,
                     if (!omit_points) obs_pts$y), na.rm = TRUE)

        if (is.finite(ylo) && is.finite(yhi)) {
            pad <- 0.04 * (yhi - ylo)
            if (!is.finite(pad) || pad == 0) pad <- 0.2
            ylim <- c(max(0, ylo - pad), 1.3 * (yhi + pad))
        } else {
            ylim <- c(0, 1)
        }
    }

    # Preserve behavior of the original function.
    ylim[2] <- 1.4 * max(fit_all$mu, na.rm = TRUE)

    if (plot_normalized_pt) {
        x_title <- "Pseudotime"
        x_breaks <- x_labels <- seq(0, 1, by = 0.2)
    } else {
        x_title <- paste0("Time", if (!is.null(time_unit)) paste0(" (", time_unit, ")") else "")
        if (is.null(x_tick_labels)) {
            tmax <- max(vapply(keys, get_time_max, numeric(1)), na.rm = TRUE)
            x_tick_labels <- seq(0, tmax, by = x_tick_step)
        }
        x_breaks <- x_labels <- x_tick_labels
    }

    if (is.null(xlim)) {
        xlim <- range(c(fit_all$x, obs_mea$x_mid,
                        if (!omit_points) obs_pts$x), na.rm = TRUE)
    }

    oob_fun <- if (requireNamespace("scales", quietly = TRUE))
        scales::squish else function(x, ...) x

    p <- ggplot2::ggplot()

    if (!omit_points) {
        if (color_points_by_run) {
            p <- p + ggplot2::geom_point(
                data = obs_pts, ggplot2::aes(x, y, colour = run),
                size = points_size, alpha = points_alpha
            )
        } else {
            p <- p + ggplot2::geom_point(
                data = obs_pts, ggplot2::aes(x, y),
                size = points_size, alpha = points_alpha, colour = points_color
            )
        }
    }

    p +
        ggplot2::geom_ribbon(
            data = fit_all, ggplot2::aes(x, ymin = lo, ymax = hi, fill = run),
            alpha = ribbon_alpha, colour = NA
        ) +
        ggplot2::geom_line(
            data = fit_all, ggplot2::aes(x, mu, colour = run),
            linewidth = line_width
        ) +
        ggplot2::geom_point(
            data = obs_mea, ggplot2::aes(x_mid, mean_norm_log1p, colour = run),
            shape = 21, fill = "white", size = 3.2, stroke = 0.9
        ) +
        ggplot2::scale_colour_manual(
            values = setNames(colours[keys], labels), name = "Run"
        ) +
        ggplot2::scale_fill_manual(
            values = setNames(fill_cols[keys], labels), name = "Run"
        ) +
        ggplot2::scale_y_continuous(
            breaks = if (requireNamespace("scales", quietly = TRUE))
                scales::pretty_breaks(n = ny_ticks) else ggplot2::waiver(),
            name = "Expression"
        ) +
        ggplot2::scale_x_continuous(
            breaks = x_breaks, labels = x_labels, limits = xlim, oob = oob_fun
        ) +
        ggplot2::coord_cartesian(expand = FALSE, ylim = ylim) +
        ggplot2::labs(
            title = if (show_title) gene else NULL,
            subtitle = if (show_title) sub_txt else NULL,
            x = x_title, y = "Expression"
        ) +
        ggplot2::theme_minimal(base_size = axis_text_size) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_line(colour = "grey88"),
            panel.grid.major.y = ggplot2::element_line(colour = "grey88"),
            axis.text = ggplot2::element_text(colour = "black", size = axis_text_size),
            axis.title = ggplot2::element_text(colour = "black", size = axis_title_size),
            plot.title = ggplot2::element_text(face = "bold", hjust = 0, colour = "black"),
            plot.subtitle = ggplot2::element_text(size = axis_text_size, colour = "black"),
            legend.position = if (show_legend) "right" else "none"
        )
}

plot_tradeSeq_gene_multirun_pub <- function(
  gene,
  sce_ls, counts_sub_ls, PT_ls, sf_ls,
  CW_ls = NULL,                    # optional named list of cellWeights per run (nCells x 1 ok)
  q_by_gene = NULL,

  keys = NULL,                     # which runs to include; default = names(sce_ls)
  key_labels = NULL,               # optional named vector/list mapping key -> legend label

  # time handling
  offset_map = NULL,
  time_max_map = NULL,             # named numeric vector/list: e.g. c(THP1=110, K562_1=90, K562_2=120)
  plot_normalized_pt = FALSE,      # FALSE => real time axis (pt * time_max); TRUE => 0..1
  trim_before_zero = FALSE,        # remove anything plotted before t=0
  point_tmax_map = NULL,           # NEW: per-run max plotted time for observed/mean points only

  time_unit = "min",
  x_tick_step = 20,                # in real-time units if plot_normalized_pt=FALSE
  x_tick_labels = NULL,            # optional: custom ticks in real-time units
  xlim = NULL,                     # optional numeric c(lo, hi) in plotted x-units

  # plot elements
  omit_points = FALSE,
  nPoints = 1000,
  ny_ticks = 5,
  ylim = NULL,

  # raw points style
  color_points_by_run = TRUE,
  points_color = "#274C77",        # used only if color_points_by_run=FALSE
  points_alpha = 0.30,
  points_size  = 0.75,

  se_factor = 1.96,

  # fitted curve style
  line_width = 1.7,
  ribbon_alpha = 0.35,

  # labeling
  show_title = TRUE,
  show_legend = TRUE,

  axis_text_size = 14,
  axis_title_size = 16,

  # colors per run
  colours = NULL,                  # named vector by key or recycled
  ribbon_alpha_mult = 0.30         # fill alpha multiplier via adjustcolor
) {

  # -----------------------------
  # 0) keys + basic checks
  # -----------------------------
  if (is.null(keys)) keys <- names(sce_ls)

  keys <- intersect(
    keys,
    Reduce(
      intersect,
      list(
        names(sce_ls),
        names(counts_sub_ls),
        names(PT_ls),
        names(sf_ls)
      )
    )
  )

  if (!length(keys)) {
    stop("No overlapping keys across sce_ls/counts_sub_ls/PT_ls/sf_ls.")
  }

  # keep only runs where gene exists
  keys <- keys[
    vapply(
      keys,
      function(k) {
        gene %in% rownames(counts_sub_ls[[k]]) &&
          gene %in% rownames(sce_ls[[k]])
      },
      logical(1)
    )
  ]

  if (!length(keys)) {
    stop("Gene not found in counts_sub and SCE for any provided key.")
  }

  # labels
  if (is.null(key_labels)) {
    key_labels <- setNames(keys, keys)
    print(key_labels)
  }

  # colours
  if (is.null(colours)) {

    base_cols <- c(
      "#1b9e77", "#d95f02", "#7570b3", "#e7298a",
      "#66a61e", "#e6ab02", "#a6761d", "#666666"
    )

    colours <- rep(base_cols, length.out = length(keys))
    names(colours) <- keys

  } else {

    if (is.null(names(colours))) {

      colours <- rep(colours, length.out = length(keys))
      names(colours) <- keys

    } else {

      colours <- colours[keys]

      if (any(is.na(colours))) {
        stop("Named `colours` must include all requested keys.")
      }
    }
  }

  # fill colors for ribbons
  fill_cols <- setNames(
    grDevices::adjustcolor(
      colours,
      alpha.f = ribbon_alpha_mult
    ),
    keys
  )

  # time_max_map parsing
  get_time_max <- function(k) {

    if (isTRUE(plot_normalized_pt)) return(1)

    if (
      is.null(time_max_map) ||
      is.na(time_max_map[[k]])
    ) {

      stop(
        paste0(
          "time_max_map must provide a value for key: ",
          k,
          " (or set plot_normalized_pt=TRUE)."
        )
      )
    }

    as.numeric(time_max_map[[k]])
  }

  # per-run x offset (in plotted x units)
  get_offset <- function(k) {

    if (is.null(offset_map)) return(0)

    v <- offset_map[[k]]

    if (is.null(v) || is.na(v)) {
      0
    } else {
      as.numeric(v)
    }
  }

  # helper: scale pt -> x
  pt_to_x <- function(pt, tm, off) {

    x0 <- if (isTRUE(plot_normalized_pt)) {
      pt
    } else {
      pt * tm
    }

    x0 + off
  }

  # optional scales::squish for oob
  oob_fun <- if (
    requireNamespace("scales", quietly = TRUE)
  ) {
    scales::squish
  } else {
    function(x, ...) x
  }

  # -----------------------------
  # 1) build fitted + observed data per run
  # -----------------------------
  fit_list <- vector("list", length(keys))
  names(fit_list) <- keys

  obs_list <- vector("list", length(keys))
  names(obs_list) <- keys

  for (k in keys) {

    sce <- sce_ls[[k]]
    counts_sub <- counts_sub_ls[[k]]
    PT <- PT_ls[[k]]
    sf <- sf_ls[[k]]

    # tradeSeq dm for this run
    dm <- SummarizedExperiment::colData(sce)$tradeSeq$dm

    if (!("t1" %in% colnames(dm))) {

      # still try to detect t columns; but you said single-lineage, so expect t1
      t_cols <- grep(
        "^t[0-9]+$",
        colnames(dm),
        value = TRUE
      )

      if (!length(t_cols)) {
        stop(
          paste0(
            "[",
            k,
            "] Could not find t1/t2/... columns in dm."
          )
        )
      }

      # pick the first available lineage as "1" for plotting
      ell <- as.integer(sub("^t", "", t_cols[1]))

    } else {

      ell <- 1L
    }

    # size-factor ref (normalize to geometric mean library size)
    sf_vec <- as.numeric(sf)
    s_ref <- exp(mean(log(sf_vec), na.rm = TRUE))
    off_ref <- log(s_ref)

    pr <- predict_tradeSeq_sce_with_se(
      sce = sce,
      gene = gene,
      lineageId = ell,
      nPoints = nPoints,
      off_ref = off_ref,
      se_factor = se_factor,
      return_log1p = TRUE
    )
    t_pt <- pr$pt

    tm <- get_time_max(k)
    off <- get_offset(k)

    fit_df <- data.frame(
      key = k,
      run = factor(
        unname(key_labels[[k]]),
        levels = unname(key_labels[keys])
      ),
      pt = t_pt,
      x = pt_to_x(t_pt, tm, off),
      mu = pr$mu,
      lo = pr$lo,
      hi = pr$hi
    )

    fit_list[[k]] <- fit_df

    # observed points (normalized counts)
    y_cnt <- as.numeric(counts_sub[gene, ])
    y_cnt_norm <- y_cnt * (s_ref / sf_vec)
    pt_vec <- as.numeric(PT)  # single-lineage PT vector

    obs_df <- data.frame(
      key = k,
      run = factor(
        unname(key_labels[[k]]),
        levels = unname(key_labels[keys])
      ),
      pt = pt_vec,
      x = pt_to_x(pt_vec, tm, off),
      y = log1p(y_cnt_norm)
    )

    # optional weights => binned weighted mean curve
    # (helps when PT is continuous)
    CW <- if (
      !is.null(CW_ls) &&
      !is.null(CW_ls[[k]])
    ) {
      CW_ls[[k]]
    } else {
      NULL
    }

    if (is.null(CW)) {

      w <- rep(1, nrow(obs_df))

    } else {

      CWm <- as.matrix(CW)

      # accept nCells x 1 OR nCells x >=1
      # (take column 1 for single-lineage)
      w <- as.numeric(CWm[, 1])
      w[!is.finite(w)] <- 0
    }

    # bin PT into nPoints bins
    # (same count as fitted resolution)
    nb <- min(nPoints, 200)

    br <- seq(
      0,
      1,
      length.out = nb + 1
    )

    bin <- cut(
      obs_df$pt,
      breaks = br,
      include.lowest = TRUE,
      right = TRUE
    )

    wmean <- function(x, w) {

      s <- sum(w, na.rm = TRUE)

      if (is.finite(s) && s > 0) {
        sum(x * w, na.rm = TRUE) / s
      } else {
        NA_real_
      }
    }

    pt_mid <- tapply(
      obs_df$pt,
      bin,
      function(v) mean(v, na.rm = TRUE)
    )

    y_mean <- tapply(
      seq_len(nrow(obs_df)),
      bin,
      function(ii) {
        wmean(
          exp(obs_df$y[ii]) - 1,
          w[ii]
        )
      }
    )

    mean_df <- data.frame(
      key = k,
      run = factor(
        unname(key_labels[[k]]),
        levels = unname(key_labels[keys])
      ),
      pt_mid = as.numeric(pt_mid),
      x_mid = pt_to_x(
        as.numeric(pt_mid),
        tm,
        off
      ),
      mean_norm_log1p = log1p(
        as.numeric(y_mean)
      )
    )

    mean_df <- mean_df[
      is.finite(mean_df$pt_mid) &
        is.finite(mean_df$mean_norm_log1p),
      ,
      drop = FALSE
    ]

    obs_list[[k]] <- list(
      points = obs_df,
      means = mean_df
    )
  }

  fit_all <- do.call(
    rbind,
    fit_list
  )

  obs_pts <- do.call(
    rbind,
    lapply(obs_list, `[[`, "points")
  )

  obs_mea <- do.call(
    rbind,
    lapply(obs_list, `[[`, "means")
  )

  # -----------------------------
  # optionally remove anything before t = 0
  # -----------------------------
  if (isTRUE(trim_before_zero)) {

    fit_all <- fit_all[
      fit_all$x >= 0,
      ,
      drop = FALSE
    ]

    obs_pts <- obs_pts[
      obs_pts$x >= 0,
      ,
      drop = FALSE
    ]

    obs_mea <- obs_mea[
      obs_mea$x_mid >= 0,
      ,
      drop = FALSE
    ]
  }

  # -----------------------------
  # optionally truncate observed/mean points
  # after a condition-specific plotted t_max
  # -----------------------------
  if (!is.null(point_tmax_map)) {

    keep_pts <- vapply(
      seq_len(nrow(obs_pts)),
      function(i) {
        k <- as.character(obs_pts$key[i])
        cutoff <- unname(point_tmax_map[k])

        is.na(cutoff) ||
          obs_pts$x[i] <= as.numeric(cutoff)
      },
      logical(1)
    )

    obs_pts <- obs_pts[
      keep_pts,
      ,
      drop = FALSE
    ]

    keep_mea <- vapply(
      seq_len(nrow(obs_mea)),
      function(i) {
        k <- as.character(obs_mea$key[i])
        cutoff <- unname(point_tmax_map[k])

        is.na(cutoff) ||
          obs_mea$x_mid[i] <= as.numeric(cutoff)
      },
      logical(1)
    )

    obs_mea <- obs_mea[
      keep_mea,
      ,
      drop = FALSE
    ]
  }

  # -----------------------------
  # 2) subtitle with FDR (optional)
  # -----------------------------
  sub_txt <- NULL

  if (
    !is.null(q_by_gene) &&
    gene %in% names(q_by_gene)
  ) {

    qv <- unname(q_by_gene[gene])

    if (is.finite(qv)) {
      sub_txt <- paste0(
        "FDR = ",
        formatC(
          qv,
          format = "fg",
          digits = 3
        )
      )
    }
  }

  # -----------------------------
  # 3) auto ylim if NULL
  # (include ribbon + mean dots + points if present)
  # -----------------------------
  if (is.null(ylim)) {

    ylo <- suppressWarnings(
      min(fit_all$lo, na.rm = TRUE)
    )

    yhi <- suppressWarnings(
      max(fit_all$hi, na.rm = TRUE)
    )

    ylo2 <- suppressWarnings(
      min(obs_mea$mean_norm_log1p, na.rm = TRUE)
    )

    yhi2 <- suppressWarnings(
      max(obs_mea$mean_norm_log1p, na.rm = TRUE)
    )

    ylo <- min(
      ylo,
      ylo2,
      na.rm = TRUE
    )

    yhi <- max(
      yhi,
      yhi2,
      na.rm = TRUE
    )

    if (!omit_points) {

      ylo3 <- suppressWarnings(
        min(obs_pts$y, na.rm = TRUE)
      )

      yhi3 <- suppressWarnings(
        max(obs_pts$y, na.rm = TRUE)
      )

      ylo <- min(
        ylo,
        ylo3,
        na.rm = TRUE
      )

      yhi <- max(
        yhi,
        yhi3,
        na.rm = TRUE
      )
    }

    if (
      is.finite(ylo) &&
      is.finite(yhi)
    ) {

      pad <- 0.04 * (yhi - ylo)

      if (
        !is.finite(pad) ||
        pad == 0
      ) {
        pad <- 0.2
      }

      ylim <- c(
        max(0, ylo - pad),
        1.3 * (yhi + pad)
      )

    } else {

      ylim <- c(0, 1)
    }
  }

  ylim[2] <- 1.4 * max(
    fit_all$mu,
    na.rm = TRUE
  )

  # -----------------------------
  # 4) x axis: real time vs normalized
  # -----------------------------
  if (isTRUE(plot_normalized_pt)) {

    x_title <- "Pseudotime"
    x_breaks <- seq(0, 1, by = 0.2)
    x_labels <- x_breaks
    x_limits <- c(0, 1)

  } else {

    x_title <- paste0(
      "Time",
      if (!is.null(time_unit)) {
        paste0(" (", time_unit, ")")
      } else {
        ""
      }
    )

    # choose global max time across requested keys for axis defaults
    tmax_all <- max(
      vapply(
        keys,
        get_time_max,
        numeric(1)
      ),
      na.rm = TRUE
    )

    if (is.null(x_tick_labels)) {
      x_tick_labels <- seq(
        0,
        tmax_all,
        by = x_tick_step
      )
    }

    x_breaks <- x_tick_labels
    x_labels <- x_tick_labels
    x_limits <- c(
      min(x_breaks),
      max(x_breaks)
    )
  }

  # if (!is.null(xlim)) x_limits <- xlim

  if (is.null(xlim)) {

    x_min <- min(
      fit_all$x,
      obs_mea$x_mid,
      if (!omit_points) obs_pts$x else Inf,
      na.rm = TRUE
    )

    x_max <- max(
      fit_all$x,
      obs_mea$x_mid,
      if (!omit_points) obs_pts$x else -Inf,
      na.rm = TRUE
    )

    x_limits <- c(
      x_min,
      x_max
    )

  } else {

    x_limits <- xlim
  }

  # y scale
  y_scale <- scale_y_continuous(
    breaks = if (
      requireNamespace("scales", quietly = TRUE)
    ) {
      scales::pretty_breaks(n = ny_ticks)
    } else {
      waiver()
    },
    name = "Expression"
  )

  legend_pos <- if (show_legend) {
    "right"
  } else {
    "none"
  }

  title_val <- if (show_title) {
    gene
  } else {
    NULL
  }

  subtitle_val <- if (show_title) {
    sub_txt
  } else {
    NULL
  }

  # -----------------------------
  # 5) plot
  # -----------------------------
  p <- ggplot()

  if (!omit_points) {

    if (isTRUE(color_points_by_run)) {

      p <- p +
        geom_point(
          data = obs_pts,
          aes(x, y, colour = run),
          size = points_size,
          alpha = points_alpha
        )

    } else {

      p <- p +
        geom_point(
          data = obs_pts,
          aes(x, y),
          size = points_size,
          alpha = points_alpha,
          colour = points_color
        )
    }
  }

  p <- p +
    geom_ribbon(
      data = fit_all,
      aes(
        x,
        ymin = lo,
        ymax = hi,
        fill = run
      ),
      alpha = ribbon_alpha,
      colour = NA
    ) +
    geom_line(
      data = fit_all,
      aes(
        x,
        mu,
        colour = run
      ),
      linewidth = line_width
    ) +
    geom_point(
      data = obs_mea,
      aes(
        x_mid,
        mean_norm_log1p,
        colour = run
      ),
      shape = 21,
      fill = "white",
      size = 3.2,
      stroke = 0.9
    ) +
    scale_colour_manual(
      values = setNames(
        colours[keys],
        unname(key_labels[keys])
      ),
      name = "Run"
    ) +
    scale_fill_manual(
      values = setNames(
        fill_cols[keys],
        unname(key_labels[keys])
      ),
      name = "Run"
    ) +
    y_scale +
    coord_cartesian(
      expand = FALSE,
      ylim = ylim
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      limits = x_limits,
      oob = oob_fun
    ) +
    labs(
      title = title_val,
      subtitle = subtitle_val,
      x = x_title,
      y = "Expression"
    ) +
    theme_minimal(
      base_size = axis_text_size
    ) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(
        colour = "grey88"
      ),
      panel.grid.major.y = element_line(
        colour = "grey88"
      ),
      axis.text = element_text(
        colour = "black",
        size = axis_text_size
      ),
      axis.title = element_text(
        colour = "black",
        size = axis_title_size
      ),
      plot.title = element_text(
        face = "bold",
        hjust = 0,
        colour = "black"
      ),
      plot.subtitle = element_text(
        size = axis_text_size,
        colour = "black"
      ),
      legend.position = legend_pos
    )

  return(p)
}

# Verbatim from TESTING.ipynb (cell 4), except that cluster::silhouette() is called through
# chronopack_silhouette(), which gives the same widths when the 'cluster' package is missing.
plot_cluster_medoids_and_top_genes <- function(
  traj_mat,
  X,
  tt,
  cl,
  D,
  top_n = 9,
  gene_selection = c("random", "wald"),
  gene_score_map = NULL,
  show_median_iqr = FALSE,
  show_sil_text = TRUE,
  use_silhouette_alpha = TRUE,
  show_gene_titles = FALSE,
  gene_title_mode = c("rownames", "gene_sym"),
  gene_title_size = 2.6,
  gene_title_y = NULL,
  seed = 13,
  cat_colors = NULL,
  cluster_order = NULL,
  show_cluster_labels = TRUE,
  title = NULL,
  base_size = 11,
  sil_text_y = 2.8,
  return_data = FALSE
) {
  gene_selection <- match.arg(gene_selection)
  gene_title_mode <- match.arg(gene_title_mode)

  stopifnot(is.matrix(traj_mat), is.matrix(X))
  stopifnot(!is.null(rownames(traj_mat)), !is.null(rownames(X)))
  stopifnot(length(tt) == ncol(X))
  stopifnot(length(cl) == nrow(traj_mat))

  if (!is.factor(cl)) cl <- factor(cl)
  stopifnot(nlevels(cl) >= 1)

  if (is.null(cluster_order)) {
    cl <- factor(cl)
  } else {
    cl <- factor(cl, levels = cluster_order)
    if (any(is.na(cl))) {
      stop("cluster_order does not include all cluster labels in cl.")
    }
  }

  nm <- rownames(traj_mat)

  # ------------------------------------------------------------
  # 0) Parse "{gene}_{key}" and attach score
  # ------------------------------------------------------------
  m <- regexec("^(.*)_([^_]+)$", nm, perl = TRUE)
  mm <- regmatches(nm, m)
  stopifnot(all(lengths(mm) == 3))

  gene_sym <- vapply(mm, `[[`, "", 2)
  key <- vapply(mm, `[[`, "", 3)

  waldStat_1 <- rep(NA_real_, length(nm))

  if (!is.null(gene_score_map)) {
    if (is.list(gene_score_map)) {
      gene_score_map <- unlist(gene_score_map, use.names = TRUE)
    }

    nm_map <- names(gene_score_map)

    if (is.null(nm_map) || any(!nzchar(nm_map))) {
      stop(
        "gene_score_map must be a named numeric vector or named list ",
        "with names = rownames(traj_mat)."
      )
    }

    gene_score_map <- as.numeric(gene_score_map)
    names(gene_score_map) <- nm_map

    waldStat_1 <- unname(gene_score_map[nm])
  }

  if (gene_selection == "wald" && all(is.na(waldStat_1))) {
    stop("gene_selection = 'wald' requires gene_score_map.")
  }

  wald <- data.frame(
    gene = nm,
    gene_sym = gene_sym,
    key = key,
    waldStat_1 = waldStat_1,
    stringsAsFactors = FALSE
  )

  # ------------------------------------------------------------
  # 1) Silhouette width per gene
  # ------------------------------------------------------------
  sil_obj <- chronopack_silhouette(as.integer(cl), D)

  sil_df <- data.frame(
    gene = nm,
    sil_width = sil_obj[, "sil_width"],
    stringsAsFactors = FALSE
  )

  sil_to_alpha <- function(s) {
    s <- pmax(-0.2, pmin(0.8, s))
    0.15 + (s + 0.2) / 1.0 * (1.0 - 0.15)
  }

  # ------------------------------------------------------------
  # 2) gene_info + medoids + selected genes
  # ------------------------------------------------------------
  gene_info <- data.frame(
    gene = nm,
    cluster = cl,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(
      wald %>%
        dplyr::select(
          gene,
          gene_sym,
          key,
          waldStat_1
        ),
      by = "gene"
    ) %>%
    dplyr::left_join(
      sil_df,
      by = "gene"
    )

  Dmat <- as.matrix(D)

  medoids <- gene_info %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
      medoid = {
        g <- gene
        subD <- Dmat[g, g, drop = FALSE]
        g[which.min(rowSums(subD))]
      },
      .groups = "drop"
    )

  set.seed(seed)

  select_cluster_genes <- function(
    df_cluster,
    medoid_gene,
    top_n,
    mode
  ) {
    pool <- df_cluster %>%
      dplyr::filter(gene != medoid_gene)

    if (nrow(pool) == 0) {
      return(character(0))
    }

    if (mode == "random") {
      return(
        sample(
          pool$gene,
          size = min(top_n, nrow(pool)),
          replace = FALSE
        )
      )
    }

    pool <- pool %>%
      dplyr::mutate(
        .ord_na = is.na(waldStat_1)
      ) %>%
      dplyr::arrange(
        .ord_na,
        dplyr::desc(waldStat_1),
        gene
      )

    head(
      pool$gene,
      n = min(top_n, nrow(pool))
    )
  }

  top_genes <- gene_info %>%
    dplyr::group_by(cluster) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(function(df_cl) {
      cl_val <- df_cl$cluster[1]

      medoid_gene <- medoids$medoid[
        match(cl_val, medoids$cluster)
      ]

      chosen <- select_cluster_genes(
        df_cluster = df_cl,
        medoid_gene = medoid_gene,
        top_n = top_n,
        mode = gene_selection
      )

      tibble::tibble(
        cluster = cl_val,
        medoid = medoid_gene,
        top = list(chosen),
        panel_genes = list(chosen)
      )
    })

  # If median/IQR is shown, reserve panel 1 and shift genes right
  panel_offset <- as.integer(show_median_iqr)
  n_panels <- top_n + panel_offset

  panel_tbl <- top_genes %>%
    dplyr::select(
      cluster,
      panel_genes
    ) %>%
    tidyr::unnest(panel_genes) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(
      panel_pos = dplyr::row_number() + panel_offset
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      gene = panel_genes
    )

  # ------------------------------------------------------------
  # 3) Individual-gene plot dataframe
  # ------------------------------------------------------------
  panel_gene_set <- unique(panel_tbl$gene)

  stopifnot(
    all(panel_gene_set %in% rownames(X))
  )

  Xsub <- X[
    panel_gene_set,
    ,
    drop = FALSE
  ]

  tp_levels <- colnames(Xsub)

  tt_lookup <- data.frame(
    tp = tp_levels,
    t = suppressWarnings(as.numeric(tp_levels)),
    stringsAsFactors = FALSE
  )

  if (any(!is.finite(tt_lookup$t))) {
    stopifnot(length(tt) == length(tp_levels))
    tt_lookup$t <- tt
  }

  plot_df <- as.data.frame(Xsub) %>%
    dplyr::mutate(
      gene = rownames(Xsub)
    ) %>%
    tidyr::pivot_longer(
      cols = -gene,
      names_to = "tp",
      values_to = "z"
    ) %>%
    dplyr::left_join(
      tt_lookup,
      by = "tp"
    ) %>%
    dplyr::inner_join(
      panel_tbl %>%
        dplyr::select(
          cluster,
          gene,
          panel_pos
        ),
      by = "gene"
    ) %>%
    dplyr::left_join(
      gene_info %>%
        dplyr::select(
          gene,
          sil_width,
          waldStat_1,
          gene_sym,
          key
        ),
      by = "gene"
    ) %>%
    dplyr::mutate(
      panel = factor(
        panel_pos,
        levels = seq_len(n_panels)
      ),
      alpha = if (use_silhouette_alpha) {
        sil_to_alpha(sil_width)
      } else {
        1
      },
      sil_lab = sprintf(
        "sil=%.2f",
        sil_width
      )
    )

  stopifnot(
    length(tt) == ncol(Xsub)
  )

  plot_df2 <- plot_df %>%
    dplyr::mutate(
      t_idx = match(
        tp,
        colnames(Xsub)
      ),
      t = tt[t_idx],
      z = as.numeric(z)
    )

  # ------------------------------------------------------------
  # 4) Median + IQR cluster summary for panel 1
  # ------------------------------------------------------------
  summary_df <- NULL

  if (show_median_iqr) {
    summary_genes <- intersect(
      gene_info$gene,
      rownames(X)
    )

    Xall <- X[
      summary_genes,
      ,
      drop = FALSE
    ]

    summary_df <- as.data.frame(Xall) %>%
      dplyr::mutate(
        gene = rownames(Xall),
        cluster = factor(
          gene_info$cluster[
            match(
              rownames(Xall),
              gene_info$gene
            )
          ],
          levels = levels(cl)
        )
      ) %>%
      tidyr::pivot_longer(
        cols = -c(gene, cluster),
        names_to = "tp",
        values_to = "z"
      ) %>%
      dplyr::mutate(
        t_idx = match(
          tp,
          colnames(Xall)
        ),
        t = tt[t_idx],
        z = as.numeric(z)
      ) %>%
      dplyr::group_by(
        cluster,
        t
      ) %>%
      dplyr::summarise(
        median = median(
          z,
          na.rm = TRUE
        ),
        q25 = quantile(
          z,
          0.25,
          na.rm = TRUE
        ),
        q75 = quantile(
          z,
          0.75,
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        panel = factor(
          1,
          levels = seq_len(n_panels)
        )
      )
  }

  # ------------------------------------------------------------
  # 5) Colors
  # ------------------------------------------------------------
  clus_levels <- levels(cl)
  pal <- NULL

  if (!is.null(cat_colors)) {
    if (
      !is.null(names(cat_colors)) &&
        all(nzchar(names(cat_colors)))
    ) {
      if (!all(clus_levels %in% names(cat_colors))) {
        missing_cols <- setdiff(
          clus_levels,
          names(cat_colors)
        )

        stop(
          "cat_colors is named, but missing colors for cluster(s): ",
          paste(
            missing_cols,
            collapse = ", "
          )
        )
      }

      pal <- unname(
        cat_colors[clus_levels]
      )

      names(pal) <- clus_levels
    } else {
      if (length(cat_colors) < length(clus_levels)) {
        stop(
          "cat_colors must have at least ",
          length(clus_levels),
          " entries, or be a named vector covering all clusters."
        )
      }

      pal <- setNames(
        cat_colors[
          seq_along(clus_levels)
        ],
        clus_levels
      )
    }
  } else if (length(clus_levels) == 3) {
    pal <- setNames(
      c(
        "#1f77b4",
        "#d62728",
        "#2ca02c"
      ),
      clus_levels
    )
  }

  # ------------------------------------------------------------
  # 6) Figure size
  # ------------------------------------------------------------
  n_clusters <- nlevels(cl)

  options(
    repr.plot.width = max(
      14,
      1.55 * n_panels
    ),
    repr.plot.height = max(
      4,
      1.55 * n_clusters
    )
  )

  if (is.null(title)) {
    pick_txt <- if (gene_selection == "random") {
      paste0(
        "random ",
        top_n
      )
    } else {
      paste0(
        "top ",
        top_n,
        " by score"
      )
    }

    if (show_median_iqr) {
      pick_txt <- paste0(
        "median + IQR; ",
        pick_txt
      )
    }

    title <- paste0(
      "Approach C clusters: ",
      pick_txt,
      " genes (color=cluster; opacity=membership confidence via silhouette)"
    )
  }

  strip_y_el <- if (show_cluster_labels) {
    element_text(
      angle = 0,
      size = base_size,
      face = "bold"
    )
  } else {
    element_blank()
  }

  p <- ggplot(
    plot_df2,
    aes(
      x = t,
      y = z,
      group = gene
    )
  ) +
    geom_line(
      aes(
        color = cluster,
        alpha = alpha
      ),
      linewidth = 0.95
    ) +
    facet_grid(
      cluster ~ panel,
      switch = "y"
    ) +
    scale_alpha_identity() +
    theme_classic(
      base_size = base_size
    ) +
    theme(
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = strip_y_el,
      strip.text.x = element_blank(),
      strip.text.x.top = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(),
      axis.line.y = element_line(),
      panel.spacing.x = grid::unit(
        0.25,
        "lines"
      ),
      panel.spacing.y = grid::unit(
        0.35,
        "lines"
      ),
      legend.position = "none",
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 14
      )
    ) +
    labs(
      title = title,
      x = NULL,
      y = NULL
    )

  # Median + IQR in first panel
  if (show_median_iqr) {
    p <- p +
      geom_ribbon(
        data = summary_df,
        aes(
          x = t,
          ymin = q25,
          ymax = q75,
          fill = cluster
        ),
        inherit.aes = FALSE,
        alpha = 0.22
      ) +
      geom_line(
        data = summary_df,
        aes(
          x = t,
          y = median,
          color = cluster,
          group = cluster
        ),
        inherit.aes = FALSE,
        linewidth = 1.25
      )
  }

  if (!is.null(pal)) {
    p <- p +
      scale_color_manual(
        values = pal
      )

    if (show_median_iqr) {
      p <- p +
        scale_fill_manual(
          values = pal
        )
    }
  }

  if (show_sil_text) {
    ann <- plot_df2 %>%
      dplyr::group_by(
        cluster,
        panel,
        gene
      ) %>%
      dplyr::summarise(
        alpha = dplyr::first(alpha),
        sil_lab = dplyr::first(sil_lab),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        t = min(tt) +
          0.02 * diff(range(tt)),
        z = sil_text_y
      )

    p <- p +
      geom_text(
        data = ann,
        aes(
          x = t,
          y = z,
          label = sil_lab,
          alpha = alpha,
          color = cluster
        ),
        size = 2.2,
        hjust = 0,
        vjust = 1
      )
  }

  if (show_gene_titles) {
    y_rng <- range(
      plot_df2$z,
      na.rm = TRUE
    )

    if (is.null(gene_title_y)) {
      gene_title_y <- y_rng[2] -
        0.04 * diff(y_rng)
    }

    ann_titles <- plot_df2 %>%
      dplyr::group_by(
        cluster,
        panel,
        gene,
        gene_sym
      ) %>%
      dplyr::summarise(
        alpha = dplyr::first(alpha),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        t = mean(range(tt)),
        z = gene_title_y,
        label = ifelse(
          gene_title_mode == "gene_sym" &
            !is.na(gene_sym) &
            nzchar(gene_sym),
          gene_sym,
          gene
        )
      )

    p <- p +
      geom_text(
        data = ann_titles,
        aes(
          x = t,
          y = z,
          label = label,
          alpha = alpha
        ),
        inherit.aes = FALSE,
        size = gene_title_size,
        fontface = "bold",
        hjust = 0.5,
        vjust = 1
      )
  }

  if (return_data) {
    return(
      list(
        plot = p,
        gene_info = gene_info,
        medoids = medoids,
        top_genes = top_genes,
        panel_tbl = panel_tbl,
        plot_df = plot_df2,
        summary_df = summary_df,
        pal = pal
      )
    )
  }

  p
}

plot_traj_trunc_gene <- function(traj_trunc, gene, keys = NULL,
                                 colours = NULL, lw = 1.7, show_points = FALSE,
                                 point_size = 1.2, point_alpha = 0.8,
                                 mode = c("expr", "deriv", "deriv2", "curvature")) {
  stopifnot(is.list(traj_trunc), length(traj_trunc) >= 1)
  mode <- match.arg(mode)

  if (is.null(keys)) keys <- names(traj_trunc)
  keys <- intersect(keys, names(traj_trunc))
  stopifnot(length(keys) >= 1)

  df_list <- lapply(keys, function(k) {
    mat <- traj_trunc[[k]]
    if (!(gene %in% rownames(mat))) return(NULL)

    t <- suppressWarnings(as.numeric(colnames(mat)))
    if (!all(is.finite(t))) t <- seq(0, ncol(mat) - 1)

    y <- as.numeric(mat[gene, , drop = TRUE])

    if (length(t) != length(y)) {
      stop(sprintf("Time grid and expression length mismatch for key '%s'.", k))
    }
    if (any(!is.finite(t)) || any(!is.finite(y))) {
      stop(sprintf("Non-finite values detected for key '%s'.", k))
    }

    dt <- diff(t)
    if (any(dt <= 0)) {
      stop(sprintf("Non-increasing or invalid time grid for key '%s'.", k))
    }

    if (mode == "expr") {
      return(data.frame(
        key  = k,
        time = t,
        y    = y,
        stringsAsFactors = FALSE
      ))
    }

    if (mode == "deriv") {
      # first derivative dy/dt on interval midpoints
      dy <- diff(y) / dt
      t_mid <- 0.5 * (t[-1] + t[-length(t)])

      return(data.frame(
        key  = k,
        time = t_mid,
        y    = dy,
        stringsAsFactors = FALSE
      ))
    }

    if (length(t) < 3) {
      stop(sprintf("At least 3 time points are required for mode '%s' (key '%s').",
                   mode, k))
    }

    # interior-point quantities
    t_in  <- t[2:(length(t) - 1)]
    dt_l  <- t[2:(length(t) - 1)] - t[1:(length(t) - 2)]
    dt_r  <- t[3:length(t)]       - t[2:(length(t) - 1)]

    slope_l <- (y[2:(length(y) - 1)] - y[1:(length(y) - 2)]) / dt_l
    slope_r <- (y[3:length(y)]       - y[2:(length(y) - 1)]) / dt_r

    # second derivative at interior points; works for nonuniform spacing
    d2y <- 2 * (slope_r - slope_l) / (dt_l + dt_r)

    if (mode == "deriv2") {
      return(data.frame(
        key  = k,
        time = t_in,
        y    = d2y,
        stringsAsFactors = FALSE
      ))
    }

    if (mode == "curvature") {
      # first derivative at interior points (simple central secant)
      dy_c <- (y[3:length(y)] - y[1:(length(y) - 2)]) /
              (t[3:length(t)] - t[1:(length(t) - 2)])

      kappa <- abs(d2y) / (1 + dy_c^2)^(3/2)

      return(data.frame(
        key  = k,
        time = t_in,
        y    = kappa,
        stringsAsFactors = FALSE
      ))
    }

    stop("Unsupported mode.")
  })

  df <- do.call(rbind, df_list)
  stopifnot(!is.null(df) && nrow(df) > 0)

  df$key <- factor(df$key, levels = keys)

  # colors
  if (is.null(colours)) {
    base_cols <- c("#222021", "#386cb0", "#DC143C", "#1b9e77", "#d95f02", "#7570b3")
    colours <- setNames(rep(base_cols, length.out = length(keys)), keys)
  } else {
    if (is.null(names(colours))) {
      colours <- setNames(rep(colours, length.out = length(keys)), keys)
    } else {
      colours <- colours[keys]
    }
  }

  title_txt <- switch(
    mode,
    expr = gene,
    deriv = paste0(gene, " (first derivative)"),
    deriv2 = paste0(gene, " (second derivative)"),
    curvature = paste0(gene, " (curvature)")
  )

  y_lab <- switch(
    mode,
    expr = "Expression",
    deriv = "d(Expression)/dt",
    deriv2 = "d2(Expression)/dt2",
    curvature = "Curvature"
  )

  p <- ggplot(df, aes(time, y, colour = key, group = key)) +
    geom_line(linewidth = lw)

  if (isTRUE(show_points)) {
    p <- p + geom_point(size = point_size, alpha = point_alpha)
  }

  if (mode %in% c("deriv", "deriv2")) {
    p <- p + geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50")
  }

  p +
    scale_colour_manual(values = colours, name = "Run") +
    theme_classic(base_size = 13) +
    labs(
      title = title_txt,
      x = "Time (min)",
      y = y_lab
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(face = "bold", hjust = 0)
    )
}

plot_cluster_violin <- function(
  pc_df,
  feature = "t_peak",
  y_title = NULL,
  plot_title = NULL,
  compare_mode = c("across_keys_for_cluster", "across_clusters_for_key"),
  cluster_col = "cluster",
  keys_use = c("K562a", "K562b", "THP1"),
  cluster_use = 1,
  key_use = "K562b",
  clusters_use = NULL,
  use_top_k = TRUE,
  top_k = 50,
  topk_mode = c("per_group", "shared_from_ref", "union", "per_key"),
  ref_group = NULL,
  ref_key = "K562b",   # kept for backward compatibility
  trim_percentiles = FALSE,
  trim_lower = 0.05,
  trim_upper = 0.95,
  print_summary = TRUE,
  print_top_genes = TRUE,
  plot_width = 7.2,
  plot_height = 4.8
) {
  compare_mode <- match.arg(compare_mode)
  topk_mode <- match.arg(topk_mode)

  # backward compatibility
  if (identical(topk_mode, "per_key")) {
    topk_mode <- "per_group"
  }

  df0 <- pc_df

  if (!(feature %in% colnames(df0))) {
    stop("feature column not found in pc_df: ", feature)
  }
  if (!(cluster_col %in% colnames(df0))) {
    stop("cluster_col not found in pc_df: ", cluster_col)
  }

  # ---- ensure key exists ----
  if (!"key" %in% colnames(df0)) {
    df0$key <- df0$cell_type
  }

  # ---- parse base gene symbol (remove trailing _celltype) ----
  if (!"gene_symbol" %in% colnames(df0)) {
    df0$gene_symbol <- df0$gene
  }

  # ---- establish comparison mode ----
  if (compare_mode == "across_keys_for_cluster") {
    group_col <- "key"
    group_levels <- keys_use

    df_cl <- df0 %>%
      dplyr::filter(
        .data[[cluster_col]] == cluster_use,
        key %in% group_levels,
        !is.na(.data[[feature]]), is.finite(.data[[feature]]),
        !is.na(waldStat_1), is.finite(waldStat_1)
      ) %>%
      dplyr::mutate(
        .group = factor(as.character(key), levels = as.character(group_levels))
      )

    compare_desc <- paste0("compare keys within ", cluster_col, " = ", cluster_use)

    if (is.null(ref_group)) {
      ref_group <- ref_key
    }

  } else if (compare_mode == "across_clusters_for_key") {
    group_col <- cluster_col

    if (is.null(clusters_use)) {
      clusters_use <- df0 %>%
        dplyr::filter(key == key_use) %>%
        dplyr::pull(.data[[cluster_col]]) %>%
        unique()

      if (is.numeric(clusters_use)) {
        clusters_use <- sort(clusters_use)
      }
    }

    group_levels <- clusters_use

    df_cl <- df0 %>%
      dplyr::filter(
        key == key_use,
        .data[[cluster_col]] %in% group_levels,
        !is.na(.data[[feature]]), is.finite(.data[[feature]]),
        !is.na(waldStat_1), is.finite(waldStat_1)
      ) %>%
      dplyr::mutate(
        .group = factor(as.character(.data[[cluster_col]]), levels = as.character(group_levels))
      )

    compare_desc <- paste0("compare ", cluster_col, " within key = ", key_use)

    if (is.null(ref_group)) {
      ref_group <- as.character(group_levels[[1]])
    }
  }

  if (nrow(df_cl) == 0) {
    stop("No rows remain after filtering. Check compare_mode / cluster_col / cluster_use / key_use / keys_use.")
  }

  selected_shared_genes <- NULL
  selected_union_genes <- NULL

  # ---- choose plotting genes ----
  if (!isTRUE(use_top_k)) {
    df_plot <- df_cl
    selection_info <- "all genes"

  } else if (topk_mode == "per_group") {
    df_plot <- df_cl %>%
      dplyr::group_by(.group) %>%
      dplyr::arrange(dplyr::desc(waldStat_1), .by_group = TRUE) %>%
      dplyr::slice_head(n = top_k) %>%
      dplyr::ungroup()

    selection_info <- paste0("top ", top_k, " genes/group by waldStat_1")

  } else if (topk_mode == "shared_from_ref") {
    if (!(as.character(ref_group) %in% as.character(group_levels))) {
      stop("ref_group must be one of the groups being compared when topk_mode = 'shared_from_ref'.")
    }

    ref_genes <- df_cl %>%
      dplyr::filter(as.character(.group) == as.character(ref_group)) %>%
      dplyr::arrange(dplyr::desc(waldStat_1)) %>%
      dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
      dplyr::slice_head(n = top_k) %>%
      dplyr::pull(gene_symbol)

    if (length(ref_genes) == 0) {
      stop("No reference genes found for ref_group after filtering.")
    }

    df_plot <- df_cl %>%
      dplyr::filter(gene_symbol %in% ref_genes)

    selection_info <- paste0(
      "shared top ", top_k, " genes from group ", ref_group, " by waldStat_1"
    )
    selected_shared_genes <- ref_genes

  } else if (topk_mode == "union") {
    union_genes <- df_cl %>%
      dplyr::group_by(.group) %>%
      dplyr::arrange(dplyr::desc(waldStat_1), .by_group = TRUE) %>%
      dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
      dplyr::slice_head(n = top_k) %>%
      dplyr::ungroup() %>%
      dplyr::pull(gene_symbol) %>%
      unique()

    if (length(union_genes) == 0) {
      stop("No union genes found after filtering.")
    }

    df_plot <- df_cl %>%
      dplyr::filter(gene_symbol %in% union_genes)

    selection_info <- paste0(
      "union of top ", top_k, " genes/group by waldStat_1"
    )
    selected_union_genes <- union_genes
  }

  # ---- optionally trim plotting range by percentiles ----
  n_removed_trim <- 0
  if (isTRUE(trim_percentiles)) {
    if (!is.numeric(trim_lower) || !is.numeric(trim_upper) ||
        length(trim_lower) != 1 || length(trim_upper) != 1 ||
        is.na(trim_lower) || is.na(trim_upper) ||
        trim_lower < 0 || trim_upper > 1 || trim_lower >= trim_upper) {
      stop("trim_lower and trim_upper must satisfy 0 <= trim_lower < trim_upper <= 1.")
    }

    q_lo <- stats::quantile(df_plot[[feature]], trim_lower, na.rm = TRUE)
    q_hi <- stats::quantile(df_plot[[feature]], trim_upper, na.rm = TRUE)

    n_before <- nrow(df_plot)
    df_plot <- df_plot %>%
      dplyr::filter(.data[[feature]] >= q_lo, .data[[feature]] <= q_hi)
    n_removed_trim <- n_before - nrow(df_plot)

    selection_info <- paste0(
      selection_info,
      "; keeping ",
      100 * trim_lower, "th-",
      100 * trim_upper, "th percentile of ", feature
    )
  }

  if (nrow(df_plot) == 0) {
    stop("No rows remain in df_plot after applying the selection filters.")
  }

  # ---- summary table ----
  summary_tbl <- df_plot %>%
    dplyr::group_by(.group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      n_unique_gene_symbols = dplyr::n_distinct(gene_symbol),
      mean_feature = mean(.data[[feature]], na.rm = TRUE),
      median_feature = stats::median(.data[[feature]], na.rm = TRUE),
      q25_feature = stats::quantile(.data[[feature]], 0.25, na.rm = TRUE),
      q75_feature = stats::quantile(.data[[feature]], 0.75, na.rm = TRUE),
      mean_wald = mean(waldStat_1, na.rm = TRUE),
      min_wald = min(waldStat_1, na.rm = TRUE),
      max_wald = max(waldStat_1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(group = .group)

  # ---- Wilcoxon tests ----
  wilcox_tbl <- NULL
  if (length(unique(df_plot$.group)) >= 2) {
    group_pairs <- combn(as.character(unique(df_plot$.group)), 2, simplify = FALSE)

    wilcox_tbl <- dplyr::bind_rows(lapply(group_pairs, function(gp) {
      x <- df_plot[[feature]][as.character(df_plot$.group) == gp[1]]
      y <- df_plot[[feature]][as.character(df_plot$.group) == gp[2]]

      wt <- suppressWarnings(stats::wilcox.test(x, y, exact = FALSE))

      data.frame(
        group1 = gp[1],
        group2 = gp[2],
        n1 = length(x),
        n2 = length(y),
        median1 = stats::median(x, na.rm = TRUE),
        median2 = stats::median(y, na.rm = TRUE),
        p_value = wt$p.value,
        stringsAsFactors = FALSE
      )
    })) %>%
      dplyr::mutate(
        p_adj_BH = stats::p.adjust(p_value, method = "BH")
      ) %>%
      dplyr::arrange(p_value)
  }

  if (isTRUE(print_summary)) {
    cat("\n===", compare_desc, ":", feature, "summary by", group_col, "===\n")
    cat("Selection mode:", selection_info, "\n")
    if (isTRUE(trim_percentiles)) {
      cat("Rows removed by percentile trimming:", n_removed_trim, "\n")
    }
    print(summary_tbl)

    if (!is.null(wilcox_tbl)) {
      cat("\n=== Pairwise Wilcoxon rank-sum tests for", feature, "===\n")
      print(wilcox_tbl)
    }
  }

  # ---- print selected genes ----
  if (isTRUE(print_top_genes) && isTRUE(use_top_k)) {
    if (topk_mode == "per_group") {
      cat("\n=== Top selected genes per group (first 20 shown) ===\n")
      print(
        df_plot %>%
          dplyr::group_by(.group) %>%
          dplyr::arrange(dplyr::desc(waldStat_1), .by_group = TRUE) %>%
          dplyr::summarise(
            top_genes = paste(head(gene, 20), collapse = ", "),
            .groups = "drop"
          ) %>%
          dplyr::rename(group = .group)
      )

    } else if (topk_mode == "shared_from_ref") {
      cat("\n=== Shared reference gene set (first 50 shown) ===\n")
      cat(paste(head(selected_shared_genes, 50), collapse = ", "), "\n")

      cat("\n=== Genes from shared set present in each group after filters (first 20 shown) ===\n")
      print(
        df_plot %>%
          dplyr::group_by(.group) %>%
          dplyr::summarise(
            genes_present = paste(head(unique(gene_symbol), 20), collapse = ", "),
            .groups = "drop"
          ) %>%
          dplyr::rename(group = .group)
      )

    } else if (topk_mode == "union") {
      cat("\n=== Union gene set (first 50 shown) ===\n")
      cat(paste(head(selected_union_genes, 50), collapse = ", "), "\n")

      cat("\n=== Genes from union set present in each group after filters (first 20 shown) ===\n")
      print(
        df_plot %>%
          dplyr::group_by(.group) %>%
          dplyr::summarise(
            genes_present = paste(head(unique(gene_symbol), 20), collapse = ", "),
            .groups = "drop"
          ) %>%
          dplyr::rename(group = .group)
      )

      cat("\n=== Top-k contributors from each group used to build the union (first 20 shown) ===\n")
      print(
        df_cl %>%
          dplyr::group_by(.group) %>%
          dplyr::arrange(dplyr::desc(waldStat_1), .by_group = TRUE) %>%
          dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
          dplyr::slice_head(n = top_k) %>%
          dplyr::summarise(
            top_genes = paste(head(gene_symbol, 20), collapse = ", "),
            .groups = "drop"
          ) %>%
          dplyr::rename(group = .group)
      )
    }
  }

  # ---- title ----
  if (is.null(plot_title)) {
    if (compare_mode == "across_keys_for_cluster") {
      plot_title <- paste0(
        cluster_col, " ", cluster_use, ": ", feature,
        " across keys (", selection_info, ")"
      )
    } else {
      plot_title <- paste0(
        "key ", key_use, ": ", feature,
        " across ", cluster_col, " (", selection_info, ")"
      )
    }
  }

  if (is.null(y_title)) {
    y_title <- feature
  }

  # ---- plot ----
  options(repr.plot.width = plot_width, repr.plot.height = plot_height)

  p <- ggplot(df_plot, aes(x = .group, y = .data[[feature]], fill = .group)) +
    geom_violin(trim = FALSE, alpha = 0.35, color = NA) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.75) +
    geom_jitter(width = 0.12, height = 0, size = 0.5, alpha = 0.35) +
    theme_classic(base_size = 12) +
    labs(
      title = plot_title,
      x = NULL,
      y = y_title
    ) +
    guides(fill = "none")

  if (identical(feature, "t_peak") && !isTRUE(trim_percentiles)) {
    p <- p + scale_y_continuous(limits = c(0, 1))
  }

  print(p)

  invisible(list(
    plot = p,
    df_plot = df_plot,
    summary = summary_tbl,
    wilcox = wilcox_tbl,
    shared_genes = selected_shared_genes,
    union_genes = selected_union_genes,
    n_removed_trim = n_removed_trim,
    compare_mode = compare_mode,
    cluster_col = cluster_col,
    group_col = group_col,
    group_levels = group_levels
  ))
}

PlotGeneKeyHighlight <- function(
  query_gene,
  pc_df,
  pt_size_bg = 0.9,
  pt_size_hi = 2.2,
  lab_size   = 10,
  base_size  = 11,
  width      = 14,
  height     = 5.5,
  key_colors = NULL,   # named vector c(K562a="#..", K562b="#..", THP1="#..") to override
  legend_key_width = 0.7,
  quiet      = FALSE
) {
  stopifnot(is.character(query_gene), length(query_gene) == 1)
  stopifnot(all(c("gene","PC1","PC2") %in% names(pc_df)))

get_key <- function(x) sub("^[^_]+_", "", as.character(x))

    
  pc_df <- pc_df %>% mutate(key = get_key(rownames(.)))

  keys_present <- sort(unique(pc_df$key))
  if (length(keys_present) < 2) stop("Need >=2 shared keys between pc_df")

  target_ids <- paste0(query_gene, "_", keys_present)

  print(target_ids)

  pc_hi <- pc_df %>% filter(rownames(.) %in% target_ids)

  if (!nrow(pc_hi)) {
    stop("Query gene not found in pc_df for the expected keys.")
  }

  # ---- colors ----
  if (is.null(key_colors)) {
    base_cols <- c("#1f77b4", "#d62728", "#2ca02c", "#FFA500")
    key_colors <- setNames(
      base_cols[seq_len(min(3, length(keys_present)))],
      keys_present[seq_len(min(3, length(keys_present)))]
    )
  } else {
    stopifnot(is.character(key_colors), !is.null(names(key_colors)))
  }

  col_scale <- scale_color_manual(
    values = key_colors,
    breaks = intersect(names(key_colors), keys_present),
    drop = FALSE
  )

  add_labels <- function(p, df, x, y) {
    if (!nrow(df)) return(p)

    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p + ggrepel::geom_text_repel(
        data = df,
        aes_string(x = x, y = y, label = "gene", color = "key"),
        size = lab_size / 3.0,
        max.overlaps = Inf,
        box.padding = 0.25,
        point.padding = 0.15,
        min.segment.length = 0,
        show.legend = FALSE
      )
    } else {
      p + geom_text(
        data = df,
        aes_string(x = x, y = y, label = "gene", color = "key"),
        size = lab_size / 3.0,
        nudge_x = 0.02,
        nudge_y = 0.02,
        check_overlap = TRUE,
        show.legend = FALSE
      )
    }
  }

  p_pca <- ggplot(pc_df, aes(PC1, PC2)) +
    geom_point(color = "grey80", size = pt_size_bg, shape = 16, stroke = 0) +
    geom_point(
      data = pc_hi,
      aes(color = key),
      size = pt_size_hi,
      shape = 16,
      stroke = 0
    ) +
    col_scale +
    scale_x_continuous(breaks = c(-10, 0, 10)) +
    scale_y_continuous(breaks = c(-10, 0, 10)) +
    guides(color = guide_legend(override.aes = list(shape = 16, linetype = 0))) +
    theme_classic(base_size = base_size) +
    theme(
      plot.title        = element_text(face = "bold"),
      axis.text         = element_text(color = "black"),
      axis.ticks        = element_line(color = "black"),
      axis.title.y      = element_text(margin = margin(r = -6)),
      legend.key.width  = grid::unit(legend_key_width, "lines")
    ) +
    labs(title = query_gene, color = "key", x = "PC1", y = "PC2")

  p_pca <- add_labels(p_pca, pc_hi, "PC1", "PC2")

  options(repr.plot.width = width, repr.plot.height = height)

  if (!quiet) {
    print(p_pca)
  }

  invisible(list(p_pca = p_pca, pc_hi = pc_hi, keys = keys_present))
}

plot_selected_genes_same_plot <- function(
  traj_mat,
  tt,
  genes,
  normalize = FALSE,
  normalize_method = c("range01", "zscore"),
  match_by = c("gene_sym", "rownames"),
  keys = NULL,
  colours = NULL,
  pc_df = NULL,
  cat_to_color = NULL,
  linewidth = 1.1,
  alpha = 0.95,
  label_size = 4,
  label_nudge_x = NULL,
  base_size = 11,
  title = NULL,
  enforce_nonnegative_raw = TRUE,
  highlight_pos = NULL,
  return_data = FALSE
) {
  normalize_method <- match.arg(normalize_method)
  match_by <- match.arg(match_by)

  stopifnot(is.matrix(traj_mat) || is.data.frame(traj_mat))
  traj_mat <- as.matrix(traj_mat)
  stopifnot(!is.null(rownames(traj_mat)))
  stopifnot(length(tt) == ncol(traj_mat))
  stopifnot(length(genes) >= 1)

  if (is.null(label_nudge_x)) {
    label_nudge_x <- 0.03 * diff(range(tt))
  }

  if (!is.null(highlight_pos)) {
    if (length(highlight_pos) != 1 || is.na(highlight_pos) ||
        highlight_pos < 1 || highlight_pos > length(genes)) {
      stop("highlight_pos must be a single integer between 1 and length(genes).")
    }
    highlight_pos <- as.integer(highlight_pos)
  }

  nm <- rownames(traj_mat)

  # Parse rownames like "{gene}_{key}" using the last underscore
  m  <- regexec("^(.*?)(?:_([^_]+))?$", nm, perl = TRUE)
  mm <- regmatches(nm, m)

  gene_sym <- vapply(mm, function(x) x[2], character(1))
  key <- vapply(
    mm,
    function(x) if (length(x) >= 3 && nzchar(x[3])) x[3] else NA_character_,
    character(1)
  )

  gene_map <- data.frame(
    row_id = nm,
    gene_sym = gene_sym,
    key = key,
    stringsAsFactors = FALSE
  )

  keep <- if (match_by == "gene_sym") {
    gene_map$gene_sym %in% genes
  } else {
    gene_map$row_id %in% genes
  }

  if (!is.null(keys)) {
    keep <- keep & gene_map$key %in% keys
  }

  sel <- gene_map[keep, , drop = FALSE]

  if (nrow(sel) == 0) {
    stop("No requested genes were found in traj_mat.")
  }

  # Preserve requested order
  if (match_by == "gene_sym") {
    sel$requested_idx <- match(sel$gene_sym, genes)
  } else {
    sel$requested_idx <- match(sel$row_id, genes)
  }
  sel <- sel[order(sel$requested_idx, sel$key, sel$row_id), , drop = FALSE]

  # If a requested gene maps to multiple rows, label with full rowname
  dup_gene_sym <- duplicated(sel$gene_sym) | duplicated(sel$gene_sym, fromLast = TRUE)
  sel$label <- ifelse(dup_gene_sym | match_by == "rownames", sel$row_id, sel$gene_sym)

  found_ids <- if (match_by == "gene_sym") unique(sel$gene_sym) else unique(sel$row_id)
  missing_ids <- setdiff(genes, found_ids)
  if (length(missing_ids) > 0) {
    warning(
      "These requested entries were not found: ",
      paste(missing_ids, collapse = ", ")
    )
  }

  # Optional broad_cluster lookup from pc_df
  sel$broad_cluster <- NA_character_
  if (!is.null(pc_df)) {
    if (!is.data.frame(pc_df)) stop("pc_df must be a data.frame.")
    if (!"broad_cluster" %in% colnames(pc_df)) {
      stop("pc_df must contain a column named 'broad_cluster'.")
    }
    if (is.null(rownames(pc_df))) {
      stop("pc_df must have rownames matching rownames(traj_mat).")
    }

    idx <- match(sel$row_id, rownames(pc_df))
    sel$broad_cluster <- as.character(pc_df$broad_cluster[idx])
  }

  # Color grouping:
  # 1) if pc_df + cat_to_color provided -> color by broad_cluster
  # 2) else color by gene label
  use_broad_cluster_colors <- !is.null(pc_df) && !is.null(cat_to_color)

  if (use_broad_cluster_colors) {
    if (is.null(names(cat_to_color)) || any(!nzchar(names(cat_to_color)))) {
      stop("cat_to_color must be a named vector: names are broad_cluster categories.")
    }

    missing_bc <- is.na(sel$broad_cluster) | !nzchar(sel$broad_cluster)
    if (any(missing_bc)) {
      if ("other" %in% names(cat_to_color)) {
        sel$broad_cluster[missing_bc] <- "other"
      } else {
        stop(
          "Some selected genes are missing broad_cluster in pc_df, and cat_to_color has no 'other' entry."
        )
      }
    }

    missing_cols <- setdiff(unique(sel$broad_cluster), names(cat_to_color))
    if (length(missing_cols) > 0) {
      stop(
        "cat_to_color is missing colors for: ",
        paste(missing_cols, collapse = ", ")
      )
    }

    sel$colour_group <- sel$broad_cluster
  } else {
    sel$colour_group <- sel$label
  }

  # Highlight logic
  if (is.null(highlight_pos)) {
    sel$is_highlight <- TRUE
  } else {
    target_val <- genes[highlight_pos]
    if (match_by == "gene_sym") {
      sel$is_highlight <- sel$gene_sym == target_val
    } else {
      sel$is_highlight <- sel$row_id == target_val
    }
  }

  sel$line_alpha  <- ifelse(sel$is_highlight, alpha, 0.35)
  sel$label_alpha <- ifelse(sel$is_highlight, 1.0, 0.75)

  sub_mat <- traj_mat[sel$row_id, , drop = FALSE]

  plot_df <- as.data.frame(sub_mat) %>%
    mutate(row_id = rownames(sub_mat)) %>%
    pivot_longer(cols = -row_id, names_to = "tp", values_to = "value") %>%
    left_join(sel, by = "row_id") %>%
    mutate(
      tp_idx = match(tp, colnames(sub_mat)),
      t = tt[tp_idx],
      value = as.numeric(value)
    )

  norm_fun <- function(x, method = c("range01", "zscore")) {
    method <- match.arg(method)

    if (method == "range01") {
      rng <- range(x, na.rm = TRUE)
      if (!all(is.finite(rng)) || diff(rng) == 0) {
        return(rep(0, length(x)))
      }
      return((x - rng[1]) / diff(rng))
    }

    mu <- mean(x, na.rm = TRUE)
    s  <- sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) {
      return(rep(0, length(x)))
    }
    (x - mu) / s
  }

  if (normalize) {
    plot_df <- plot_df %>%
      group_by(row_id) %>%
      mutate(y = norm_fun(value, normalize_method)) %>%
      ungroup()
  } else {
    plot_df$y <- plot_df$value

    if (enforce_nonnegative_raw && any(plot_df$y < 0, na.rm = TRUE)) {
      stop(
        "normalize = FALSE expects non-negative expression trajectories.\n",
        "You likely passed a standardized/scaled matrix instead of raw expression.\n",
        "Please supply the raw trajectory matrix to `traj_mat`."
      )
    }
  }

  if (is.null(title)) {
    title <- if (normalize) {
      paste0("Selected genes (normalized: ", normalize_method, ")")
    } else {
      "Selected genes (raw expression)"
    }
  }

  # Label positions: place labels at the right end of each curve
  label_df <- plot_df %>%
    group_by(row_id, label, colour_group, label_alpha) %>%
    filter(t == max(t, na.rm = TRUE)) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    mutate(x_lab = t + label_nudge_x)

  p <- ggplot(plot_df, aes(x = t, y = y, group = row_id, color = colour_group)) +
    geom_line(aes(alpha = line_alpha), linewidth = linewidth) +
    geom_text(
      data = label_df,
      aes(
        x = x_lab,
        y = y,
        label = label,
        color = colour_group,
        alpha = label_alpha
      ),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 0.5,
      size = label_size
    ) +
    scale_alpha_identity() +
    theme_classic(base_size = base_size) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.margin = margin(5.5, 60, 5.5, 5.5)
    ) +
    labs(
      title = title,
      x = "Time",
      y = if (normalize) "Expression z-score" else "Expression"
    ) +
    coord_cartesian(clip = "off") +
    expand_limits(x = max(label_df$x_lab, na.rm = TRUE))

  if (use_broad_cluster_colors) {
    p <- p + scale_color_manual(values = cat_to_color)
  } else if (!is.null(colours)) {
    lab_levels <- unique(sel$label)

    if (!is.null(names(colours)) && all(nzchar(names(colours)))) {
      if (!all(lab_levels %in% names(colours))) {
        missing_cols <- setdiff(lab_levels, names(colours))
        stop(
          "Named colours vector is missing: ",
          paste(missing_cols, collapse = ", ")
        )
      }
      p <- p + scale_color_manual(values = colours[lab_levels])
    } else {
      if (length(colours) < length(lab_levels)) {
        stop("Need at least ", length(lab_levels), " colours.")
      }
      p <- p + scale_color_manual(values = setNames(colours[seq_along(lab_levels)], lab_levels))
    }
  }

  if (return_data) {
    return(list(
      plot = p,
      selected_rows = sel,
      plot_df = plot_df,
      label_df = label_df
    ))
  } else {
    return(p)
  }
}

# Verbatim from TESTING.ipynb (cell 15), except one change: TESTING called
# plot_tradeSeq_gene_multirun_pub(sce_ls = ...), which it never defined (a kernel
# leftover; the counts-based package function of that name needs counts_sub_ls etc.).
# The call goes to plot_tradeSeq_gene_sce_pub, TESTING's SCE-only gene plot.
plot_cluster_pair_representative_genes <- function(
  pc_df,
  cluster_list = c(1, 2),
  row_cell_type = "K562a",
  col_cell_type = "THP1",
  cluster_col = "cluster",
  gamList_ls, sce_ls, counts_sub_ls, PT_ls, sf_ls,
  CW_ls = NULL,
  colours = c(THP1 = "#222021", K562a = "#386cb0"),
  cluster_palette = NULL,
  rank_by = c("row", "col", "sum", "min"),
  qval_cut = 0.05,
  exclude_genes = NULL,
  show_legend = FALSE,
  print_gene_table = TRUE,
  overall_title = NULL,
  border_colour = "grey65",
  border_linewidth = 0.7,
  annotation_box_linewidth = 2.2,
  annotation_box_rel_size = 0.16,
  top_annotation_box_rel_size = 0.24,
  axis_label_rel_size = 0.08,
  gene_title_size = 9,
  label_size = 5.2,
  hide_subplot_axes = TRUE,
  hide_mean_fit_points = FALSE,
  ...
) {
  rank_by <- match.arg(rank_by)

  stopifnot(is.data.frame(pc_df))
  stopifnot(
    all(
      c(
        "gene",
        "cell_type",
        cluster_col,
        "waldStat_1",
        "qval_1"
      ) %in% colnames(pc_df)
    )
  )
  stopifnot(length(cluster_list) >= 1)

  if (
    is.na(colours[[row_cell_type]]) ||
    is.na(colours[[col_cell_type]])
  ) {
    stop(
      "colours must be a named vector containing both ",
      "row_cell_type and col_cell_type."
    )
  }

  row_clusters <- as.character(cluster_list)
  col_clusters <- as.character(cluster_list)
  n_cl <- length(cluster_list)
  cluster_ids <- as.character(unique(cluster_list))

  if (is.null(cluster_palette)) {
    pal_vals <- hcl.colors(
      length(cluster_ids),
      palette = "Dark 3"
    )

    cluster_palette <- stats::setNames(
      pal_vals,
      cluster_ids
    )
  } else {
    if (!all(cluster_ids %in% names(cluster_palette))) {
      stop(
        "cluster_palette must be a named vector with names ",
        "matching as.character(cluster_list)."
      )
    }

    cluster_palette <- cluster_palette[cluster_ids]
  }

  # -----------------------------------
  # prepare shared-gene assignment table
  # -----------------------------------
  df_row <- pc_df %>%
    dplyr::filter(
      .data$cell_type == row_cell_type
    ) %>%
    dplyr::transmute(
      gene = .data$gene,
      cluster_row = as.character(.data[[cluster_col]]),
      wald_row = .data$waldStat_1,
      qval_row = .data$qval_1
    )

  df_col <- pc_df %>%
    dplyr::filter(
      .data$cell_type == col_cell_type
    ) %>%
    dplyr::transmute(
      gene = .data$gene,
      cluster_col = as.character(.data[[cluster_col]]),
      wald_col = .data$waldStat_1,
      qval_col = .data$qval_1
    )

  shared_df <- dplyr::inner_join(
    df_row,
    df_col,
    by = "gene"
  ) %>%
    dplyr::mutate(
      pass_qval_both =
        !is.na(.data$qval_row) &
        !is.na(.data$qval_col) &
        (.data$qval_row < qval_cut) &
        (.data$qval_col < qval_cut),

      score = dplyr::case_when(
        rank_by == "row" ~ .data$wald_row,
        rank_by == "col" ~ .data$wald_col,
        rank_by == "sum" ~
          .data$wald_row + .data$wald_col,
        rank_by == "min" ~
          pmin(.data$wald_row, .data$wald_col),
        TRUE ~ .data$wald_row
      )
    )

  # -----------------------------------
  # helpers
  # -----------------------------------
  add_box_theme <- function(p) {
    p +
      theme(
        panel.border = element_rect(
          colour = border_colour,
          fill = NA,
          linewidth = border_linewidth
        ),
        plot.margin = margin(2, 2, 2, 2)
      )
  }

  hide_axes_theme <- theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank()
  )

  remove_point_layers <- function(p) {
    keep <- vapply(
      p$layers,
      function(lyr) {
        !inherits(lyr$geom, "GeomPoint")
      },
      logical(1)
    )

    p$layers <- p$layers[keep]
    p
  }

  make_empty_panel <- function() {
    add_box_theme(
      ggplot() +
        xlim(0, 1) +
        ylim(0, 1) +
        theme_void()
    )
  }

  make_annotation_panel <- function(cluster_id) {
    fill_col <- cluster_palette[[as.character(cluster_id)]]

    ggplot() +
      annotate(
        "rect",
        xmin = 0,
        xmax = 1,
        ymin = 0,
        ymax = 1,
        fill = fill_col,
        colour = fill_col,
        linewidth = annotation_box_linewidth
      ) +
      xlim(0, 1) +
      ylim(0, 1) +
      theme_void() +
      theme(
        plot.margin = margin(2, 2, 2, 2)
      )
  }

  make_axis_label_panel <- function(
    label,
    angle = 0,
    face = "bold"
  ) {
    ggplot() +
      annotate(
        "text",
        x = 0.5,
        y = 0.5,
        label = label,
        angle = angle,
        hjust = 0.5,
        vjust = 0.5,
        size = label_size,
        fontface = face
      ) +
      xlim(0, 1) +
      ylim(0, 1) +
      theme_void() +
      theme(
        plot.margin = margin(0, 0, 0, 0)
      )
  }

  # -----------------------------------
  # construct ranked candidate lists
  # for every row/column cluster pair
  #
  # Preferred order:
  #   1. qval-passing, non-excluded genes
  #   2. qval-passing excluded genes as fallback
  #
  # Within each group, genes are ranked using
  # the same score as before.
  # -----------------------------------
  candidate_lists <- list()
  chosen_rows <- list()

  pair_idx <- 1L

  for (rcl in row_clusters) {
    for (ccl in col_clusters) {

      cand_all <- shared_df %>%
        dplyr::filter(
          .data$cluster_row == as.character(rcl),
          .data$cluster_col == as.character(ccl)
        )

      cand_qpass <- cand_all %>%
        dplyr::filter(
          .data$pass_qval_both
        )

      # Preferred candidates: non-excluded
      cand_primary <- cand_qpass

      if (!is.null(exclude_genes)) {
        cand_primary <- cand_primary %>%
          dplyr::filter(
            !.data$gene %in% exclude_genes
          )
      }

      cand_primary <- cand_primary %>%
        dplyr::arrange(
          dplyr::desc(.data$score),
          dplyr::desc(.data$wald_row)
        )

      # Fallback candidates: excluded genes
      if (!is.null(exclude_genes)) {
        cand_fallback <- cand_qpass %>%
          dplyr::filter(
            .data$gene %in% exclude_genes
          ) %>%
          dplyr::arrange(
            dplyr::desc(.data$score),
            dplyr::desc(.data$wald_row)
          )
      } else {
        cand_fallback <- cand_qpass[0, , drop = FALSE]
      }

      # Keep the complete ordered list.
      cand_ranked <- dplyr::bind_rows(
        cand_primary,
        cand_fallback
      ) %>%
        dplyr::distinct(
          .data$gene,
          .keep_all = TRUE
        )

      candidate_lists[[pair_idx]] <- cand_ranked

      # Placeholder row; populated with the gene that
      # actually plots successfully below.
      chosen_rows[[pair_idx]] <- data.frame(
        row_cluster = rcl,
        col_cluster = ccl,
        gene = NA_character_,
        wald_row = NA_real_,
        wald_col = NA_real_,
        qval_row = NA_real_,
        qval_col = NA_real_,
        score = NA_real_,
        n_candidates_all = nrow(cand_all),
        n_candidates_qpass = nrow(cand_ranked),
        candidate_rank = NA_integer_,
        stringsAsFactors = FALSE
      )

      pair_idx <- pair_idx + 1L
    }
  }

  chosen_df <- dplyr::bind_rows(chosen_rows)

  # -----------------------------------
  # interior plot panels
  #
  # If a candidate gene is absent from all
  # requested SCEs, try the next-ranked gene.
  # Other plotting errors are NOT swallowed.
  # -----------------------------------
  keys_use <- c(
    row_cell_type,
    col_cell_type
  )

  colours_use <- colours[keys_use]

  panel_lookup <- vector(
    "list",
    length(candidate_lists)
  )

  skipped_missing_genes <- vector(
    "list",
    length(candidate_lists)
  )

  for (i in seq_along(candidate_lists)) {

    candidates <- candidate_lists[[i]]

    if (nrow(candidates) == 0) {
      panel_lookup[[i]] <- make_empty_panel()
      skipped_missing_genes[[i]] <- character(0)
      next
    }

    successful_plot <- NULL
    successful_row <- NULL
    missing_genes <- character(0)

    for (j in seq_len(nrow(candidates))) {

      g <- candidates$gene[j]

      p_try <- tryCatch(
        {
          p_tmp <- NULL

          utils::capture.output({
            p_tmp <- suppressMessages(
              suppressWarnings(
                plot_tradeSeq_gene_sce_pub(
                  sce_ls = sce_ls,
                  keys = keys_use,
                  colours = colours_use,
                  key_labels = NULL,
                  gene = g,
                  show_legend = show_legend,
                  ...
                )
              )
            )
          })

          p_tmp
        },
        error = function(e) {

          msg <- conditionMessage(e)

          # This is the specific recoverable error.
          if (
            grepl(
              "Gene not found in any requested SCE",
              msg,
              fixed = TRUE
            )
          ) {
            missing_genes <<- c(
              missing_genes,
              g
            )

            return(NULL)
          }

          # Any other error should still stop execution,
          # since it likely represents a real bug.
          stop(e)
        }
      )

      # Missing gene -> try next-ranked candidate.
      if (is.null(p_try)) {
        next
      }

      # Successful candidate.
      successful_plot <- p_try
      successful_row <- candidates[j, , drop = FALSE]

      chosen_df$candidate_rank[i] <- j

      break
    }

    skipped_missing_genes[[i]] <- missing_genes

    # If every candidate failed because it was absent
    # from the requested SCEs, leave an empty panel.
    if (is.null(successful_plot)) {
      panel_lookup[[i]] <- make_empty_panel()
      next
    }

    # Record the gene that was actually plotted.
    chosen_df$gene[i] <- successful_row$gene
    chosen_df$wald_row[i] <- successful_row$wald_row
    chosen_df$wald_col[i] <- successful_row$wald_col
    chosen_df$qval_row[i] <- successful_row$qval_row
    chosen_df$qval_col[i] <- successful_row$qval_col
    chosen_df$score[i] <- successful_row$score

    p <- successful_plot
    g <- successful_row$gene

    if (isTRUE(hide_mean_fit_points)) {
      p <- remove_point_layers(p)
    }

    p <- p +
      ggtitle(g) +
      theme(
        plot.title = element_text(
          size = gene_title_size,
          face = "bold",
          hjust = 0.5
        ),
        plot.margin = margin(2, 2, 2, 2)
      )

    if (isTRUE(hide_subplot_axes)) {
      p <- p + hide_axes_theme
    }

    panel_lookup[[i]] <- add_box_theme(p)
  }

  # Print after plotting so this reflects
  # the gene that actually rendered.
  if (print_gene_table) {
    print(chosen_df)
  }

  # -----------------------------------
  # assemble central grid
  # -----------------------------------
  full_panels <- list()

  full_panels[[1]] <- patchwork::plot_spacer()

  for (ccl in col_clusters) {
    full_panels[[length(full_panels) + 1]] <-
      make_annotation_panel(ccl)
  }

  idx <- 1L

  for (rcl in row_clusters) {

    full_panels[[length(full_panels) + 1]] <-
      make_annotation_panel(rcl)

    for (ccl in col_clusters) {

      full_panels[[length(full_panels) + 1]] <-
        panel_lookup[[idx]]

      idx <- idx + 1L
    }
  }

  main_grid <- wrap_plots(
    full_panels,
    ncol = n_cl + 1,
    byrow = TRUE
  ) +
    plot_layout(
      widths = c(
        annotation_box_rel_size,
        rep(1, n_cl)
      ),
      heights = c(
        top_annotation_box_rel_size,
        rep(1, n_cl)
      )
    )

  y_axis_panel <- make_axis_label_panel(
    label = row_cell_type,
    angle = 90
  )

  x_axis_panel <- make_axis_label_panel(
    label = col_cell_type,
    angle = 0
  )

  x_label_row <- wrap_plots(
    list(
      patchwork::plot_spacer(),
      x_axis_panel
    ),
    ncol = 2,
    widths = c(
      axis_label_rel_size,
      1
    )
  )

  grid_row <- wrap_plots(
    list(
      y_axis_panel,
      main_grid
    ),
    ncol = 2,
    widths = c(
      axis_label_rel_size,
      1
    )
  )

  p_out <- x_label_row / grid_row +
    plot_layout(
      heights = c(
        axis_label_rel_size,
        1
      )
    )

  if (!is.null(overall_title)) {
    p_out <- p_out +
      plot_annotation(
        title = overall_title
      )
  }

  print(p_out)

  invisible(
    list(
      plot = p_out,
      chosen_genes = chosen_df,
      candidate_genes = candidate_lists,
      skipped_missing_genes = skipped_missing_genes,
      shared_df = shared_df,
      cluster_palette = cluster_palette
    )
  )
}

plot_cluster_reassignment_sankey <- function(
  pc_df,
  source_cell_type = "K562a",
  target_cell_type = "THP1",
  cluster_col = "cluster",
  gene_col = "gene",
  rank_col = "waldStat_1",
  top_k = 300,
  rank_desc = TRUE,
  clusters_use = c(1, 2, 3, 4, 9),
  min_flow = 1,
  stratum_gap = NULL,
  stratum_width = 0.28,
  show_legend = FALSE,
  cluster_palette = NULL,
  low_significance_label = "low_significance",
  low_significance_col = "grey80",
  treat_missing_as_low_significance = FALSE,
  low_significance_action = c("drop", "hide", "show"),
  title = NULL,
  subtitle = NULL,
  base_size = 13
) {

  low_significance_action <- match.arg(low_significance_action)

  stopifnot(
    is.data.frame(pc_df),
    all(c("cell_type", cluster_col, gene_col, rank_col) %in% names(pc_df))
  )

  # ------------------------------------------------------------
  # Helpers
  # ------------------------------------------------------------

  arrange_rank <- function(x) {
    if (rank_desc) dplyr::arrange(x, dplyr::desc(.data$rank_val))
    else dplyr::arrange(x, .data$rank_val)
  }

  get_ann <- function(ct, cl_name, rank_name) {
    pc_df %>%
      dplyr::filter(.data$cell_type == ct) %>%
      dplyr::transmute(
        gene = as.character(.data[[gene_col]]),
        cluster = as.character(.data[[cluster_col]]),
        rank_val = .data[[rank_col]]
      ) %>%
      dplyr::filter(!is.na(.data$rank_val), is.finite(.data$rank_val)) %>%
      arrange_rank() %>%
      dplyr::distinct(.data$gene, .keep_all = TRUE) %>%
      dplyr::rename(
        !!cl_name := .data$cluster,
        !!rank_name := .data$rank_val
      )
  }

  # ------------------------------------------------------------
  # Gene annotations + cluster universe
  # ------------------------------------------------------------

  df_src_all <- get_ann(source_cell_type, "source_cluster", "source_rank")
  df_tgt_all <- get_ann(target_cell_type, "target_cluster", "target_rank")

  base_clusters <- if (is.null(clusters_use)) {
    sort(unique(c(df_src_all$source_cluster, df_tgt_all$target_cluster)))
  } else {
    as.character(clusters_use)
  }

  base_clusters <- setdiff(unique(base_clusters), low_significance_label)

  plot_clusters <- if (low_significance_action == "show") {
    c(base_clusters, low_significance_label)
  } else {
    base_clusters
  }

  if (!length(plot_clusters)) stop("No clusters remain for plotting.")

  # ------------------------------------------------------------
  # Top-k union
  # ------------------------------------------------------------

  df_src_topk <- dplyr::slice_head(df_src_all, n = top_k)
  df_tgt_topk <- dplyr::slice_head(df_tgt_all, n = top_k)

  union_genes <- union(df_src_topk$gene, df_tgt_topk$gene)
  if (!length(union_genes)) stop("No genes remain in the top-k union.")

  flow_genes_pre <- data.frame(gene = union_genes) %>%
    dplyr::left_join(
      dplyr::filter(df_src_all, .data$gene %in% union_genes),
      by = "gene"
    ) %>%
    dplyr::left_join(
      dplyr::filter(df_tgt_all, .data$gene %in% union_genes),
      by = "gene"
    ) %>%
    dplyr::mutate(
      missing_source_assignment = is.na(.data$source_cluster),
      missing_target_assignment = is.na(.data$target_cluster)
    )

  genes_with_missing_assignments <- flow_genes_pre %>%
    dplyr::filter(
      .data$missing_source_assignment |
        .data$missing_target_assignment
    ) %>%
    dplyr::pull(.data$gene) %>%
    unique()

  imputed_missing_assignments_df <- dplyr::bind_rows(
    flow_genes_pre %>%
      dplyr::filter(.data$missing_source_assignment) %>%
      dplyr::transmute(
        gene, cell_type = source_cell_type,
        side = "source",
        imputed_cluster = low_significance_label
      ),
    flow_genes_pre %>%
      dplyr::filter(.data$missing_target_assignment) %>%
      dplyr::transmute(
        gene, cell_type = target_cell_type,
        side = "target",
        imputed_cluster = low_significance_label
      )
  )

  if (treat_missing_as_low_significance) {
    flow_genes_pre$source_cluster[
      flow_genes_pre$missing_source_assignment
    ] <- low_significance_label

    flow_genes_pre$target_cluster[
      flow_genes_pre$missing_target_assignment
    ] <- low_significance_label
  }

  # ------------------------------------------------------------
  # Classify assignments
  # ------------------------------------------------------------

  flow_genes_pre <- flow_genes_pre %>%
    dplyr::mutate(
      source_low = !is.na(.data$source_cluster) &
        .data$source_cluster == low_significance_label,
      target_low = !is.na(.data$target_cluster) &
        .data$target_cluster == low_significance_label,

      source_base = !is.na(.data$source_cluster) &
        .data$source_cluster %in% base_clusters,
      target_base = !is.na(.data$target_cluster) &
        .data$target_cluster %in% base_clusters,

      source_any = .data$source_base | .data$source_low,
      target_any = .data$target_base | .data$target_low,

      base_pair = .data$source_base & .data$target_base,
      show_pair = .data$source_any & .data$target_any,

      one_sided_low =
        (.data$source_base & .data$target_low) |
        (.data$source_low & .data$target_base)
    )

  # Genes represented by actual ribbons
  flow_genes_visible <- if (low_significance_action == "show") {
    dplyr::filter(flow_genes_pre, .data$show_pair)
  } else {
    dplyr::filter(flow_genes_pre, .data$base_pair)
  }

  hidden_low_significance_genes <- if (low_significance_action == "hide") {
    dplyr::filter(flow_genes_pre, .data$one_sided_low)
  } else {
    flow_genes_pre[0, ]
  }

  flow_genes <- switch(
    low_significance_action,
    drop = dplyr::filter(flow_genes_pre, .data$base_pair),
    hide = dplyr::filter(
      flow_genes_pre,
      .data$base_pair | .data$one_sided_low
    ),
    show = dplyr::filter(flow_genes_pre, .data$show_pair)
  )

  # ------------------------------------------------------------
  # Visible ribbons
  # ------------------------------------------------------------

  flow_df_all <- flow_genes_visible %>%
    dplyr::count(
      .data$source_cluster,
      .data$target_cluster,
      name = "n"
    )

  flow_df <- flow_df_all %>%
    dplyr::filter(.data$n >= min_flow)

  if (!nrow(flow_df_all))
    stop("No source-target cluster pairs remain.")

  # ------------------------------------------------------------
  # Bar membership
  #
  # hide: valid blocks retain genes whose other side is low
  # show: low-significance is an ordinary displayed block
  # drop: only base-to-base genes remain
  # ------------------------------------------------------------

  if (low_significance_action == "hide") {

    source_bar_genes <- flow_genes_pre %>%
      dplyr::filter(
        .data$source_base &
          (.data$target_base | .data$target_low)
      )

    target_bar_genes <- flow_genes_pre %>%
      dplyr::filter(
        .data$target_base &
          (.data$source_base | .data$source_low)
      )

  } else {

    source_bar_genes <- flow_genes_visible
    target_bar_genes <- flow_genes_visible
  }

  cluster_template <- data.frame(
    cluster = plot_clusters,
    stringsAsFactors = FALSE
  )

  source_totals <- source_bar_genes %>%
    dplyr::count(.data$source_cluster, name = "n") %>%
    dplyr::rename(cluster = .data$source_cluster) %>%
    dplyr::right_join(cluster_template, by = "cluster") %>%
    dplyr::mutate(n = tidyr::replace_na(.data$n, 0))

  target_totals <- target_bar_genes %>%
    dplyr::count(.data$target_cluster, name = "n") %>%
    dplyr::rename(cluster = .data$target_cluster) %>%
    dplyr::right_join(cluster_template, by = "cluster") %>%
    dplyr::mutate(n = tidyr::replace_na(.data$n, 0))

  # Restore explicit plotting order after joins
  source_totals <- cluster_template %>%
    dplyr::left_join(source_totals, by = "cluster")
  target_totals <- cluster_template %>%
    dplyr::left_join(target_totals, by = "cluster")

  # ------------------------------------------------------------
  # Palette
  # ------------------------------------------------------------

  if (is.null(cluster_palette)) {
    vals <- grDevices::hcl.colors(
      max(length(base_clusters), 1),
      "Dark 3"
    )[seq_along(base_clusters)]

    cluster_palette <- stats::setNames(vals, base_clusters)
  } else {
    cluster_palette <- unlist(cluster_palette)
  }

  missing_base <- setdiff(base_clusters, names(cluster_palette))
  if (length(missing_base)) {
    stop(
      "cluster_palette is missing: ",
      paste(missing_base, collapse = ", ")
    )
  }

  if (low_significance_action == "show") {
    cluster_palette[low_significance_label] <- low_significance_col
  }

  cluster_palette <- cluster_palette[plot_clusters]

  # ------------------------------------------------------------
  # Bar heights / vertical layout
  # ------------------------------------------------------------

  max_n <- max(c(source_totals$n, target_totals$n, 1))

  if (is.null(stratum_gap))
    stratum_gap <- max(0.015 * max_n, 0.6)

  empty_box_size <- max(0.15, 0.25 * stratum_gap)

  source_totals <- source_totals %>%
    dplyr::mutate(
      n_display = dplyr::if_else(
        .data$n > 0,
        as.numeric(.data$n),
        empty_box_size
      )
    )

  target_totals <- target_totals %>%
    dplyr::mutate(
      n_display = dplyr::if_else(
        .data$n > 0,
        as.numeric(.data$n),
        empty_box_size
      )
    )

  src_sum <- sum(source_totals$n_display)
  tgt_sum <- sum(target_totals$n_display)
  common_height <- max(src_sum, tgt_sum, 1)

  source_scale <- common_height / src_sum
  target_scale <- common_height / tgt_sum

  source_totals$height <- source_totals$n_display * source_scale
  target_totals$height <- target_totals$n_display * target_scale

  total_plot_height <-
    common_height +
    stratum_gap * (length(plot_clusters) - 1)

  position_blocks <- function(df, side, x) {
    top <- total_plot_height
    out <- vector("list", nrow(df))

    for (i in seq_len(nrow(df))) {
      ymax <- top
      ymin <- ymax - df$height[i]

      out[[i]] <- data.frame(
        side = side,
        cluster = df$cluster[i],
        n = df$n[i],
        height = df$height[i],
        xmin = x - stratum_width / 2,
        xmax = x + stratum_width / 2,
        ymin = ymin,
        ymax = ymax,
        x = x,
        y = (ymin + ymax) / 2
      )

      top <- ymin - stratum_gap
    }

    dplyr::bind_rows(out)
  }

  source_blocks <- position_blocks(source_totals, "source", 1)
  target_blocks <- position_blocks(target_totals, "target", 2)

  stratum_df <- dplyr::bind_rows(source_blocks, target_blocks) %>%
    dplyr::mutate(
      ymin = pmax(.data$ymin, 0),
      ymax = pmin(.data$ymax, total_plot_height),
      y = (.data$ymin + .data$ymax) / 2
    )

  # ------------------------------------------------------------
  # Ribbon coordinates
  # ------------------------------------------------------------

  cluster_order <- stats::setNames(
    seq_along(plot_clusters),
    plot_clusters
  )

  flow_layout <- flow_df %>%
    dplyr::mutate(
      source_order = unname(cluster_order[.data$source_cluster]),
      target_order = unname(cluster_order[.data$target_cluster]),
      source_height = .data$n * source_scale,
      target_height = .data$n * target_scale
    )

  src_coords <- flow_layout %>%
    dplyr::left_join(
      dplyr::select(
        source_blocks,
        source_cluster = cluster,
        source_bar_ymax = ymax
      ),
      by = "source_cluster"
    ) %>%
    dplyr::group_by(.data$source_cluster) %>%
    dplyr::arrange(.data$target_order, .by_group = TRUE) %>%
    dplyr::mutate(
      cum = cumsum(.data$source_height),
      source_ymax =
        .data$source_bar_ymax - .data$cum + .data$source_height,
      source_ymin =
        .data$source_ymax - .data$source_height
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      source_cluster, target_cluster,
      source_ymin, source_ymax
    )

  tgt_coords <- flow_layout %>%
    dplyr::left_join(
      dplyr::select(
        target_blocks,
        target_cluster = cluster,
        target_bar_ymax = ymax
      ),
      by = "target_cluster"
    ) %>%
    dplyr::group_by(.data$target_cluster) %>%
    dplyr::arrange(.data$source_order, .by_group = TRUE) %>%
    dplyr::mutate(
      cum = cumsum(.data$target_height),
      target_ymax =
        .data$target_bar_ymax - .data$cum + .data$target_height,
      target_ymin =
        .data$target_ymax - .data$target_height
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(
      source_cluster, target_cluster,
      target_ymin, target_ymax
    )

  flow_layout <- flow_layout %>%
    dplyr::left_join(
      src_coords,
      by = c("source_cluster", "target_cluster")
    ) %>%
    dplyr::left_join(
      tgt_coords,
      by = c("source_cluster", "target_cluster")
    ) %>%
    dplyr::mutate(flow_id = dplyr::row_number())

  make_ribbon <- function(
    x0, x1, y0lo, y0hi, y1lo, y1hi, id, n = 60
  ) {
    t <- seq(0, 1, length.out = n)
    s <- 3 * t^2 - 2 * t^3
    x <- x0 + (x1 - x0) * t

    lo <- y0lo + (y1lo - y0lo) * s
    hi <- y0hi + (y1hi - y0hi) * s

    data.frame(
      flow_id = id,
      x = c(x, rev(x)),
      y = c(lo, rev(hi))
    )
  }

  ribbon_df <- if (nrow(flow_layout)) {
    dplyr::bind_rows(lapply(seq_len(nrow(flow_layout)), function(i) {
      make_ribbon(
        1 + stratum_width / 2,
        2 - stratum_width / 2,
        flow_layout$source_ymin[i],
        flow_layout$source_ymax[i],
        flow_layout$target_ymin[i],
        flow_layout$target_ymax[i],
        flow_layout$flow_id[i]
      )
    }))
  } else {
    data.frame(
      flow_id = integer(),
      x = numeric(),
      y = numeric()
    )
  }

  # ------------------------------------------------------------
  # Bookkeeping
  # ------------------------------------------------------------

  genes_touching_low_significance <- flow_genes_pre %>%
    dplyr::filter(.data$source_low | .data$target_low) %>%
    dplyr::pull(.data$gene) %>%
    unique()

  dropped_low_significance_genes <- switch(
    low_significance_action,
    drop = genes_touching_low_significance,
    hide = setdiff(
      genes_touching_low_significance,
      hidden_low_significance_genes$gene
    ),
    show = character(0)
  )

  dropped_incomplete_genes <- if (treat_missing_as_low_significance) {
    character(0)
  } else {
    genes_with_missing_assignments
  }

  dropped_invalid_cluster_genes <- setdiff(
    flow_genes_pre$gene,
    unique(c(
      flow_genes$gene,
      dropped_incomplete_genes,
      dropped_low_significance_genes
    ))
  )

  hidden_flow_df <- if (low_significance_action == "hide") {
    hidden_low_significance_genes %>%
      dplyr::count(
        .data$source_cluster,
        .data$target_cluster,
        name = "n"
      )
  } else {
    flow_df_all[0, ]
  }

  # ------------------------------------------------------------
  # Plot
  # ------------------------------------------------------------

  if (is.null(title)) {
    title <- paste0(
      source_cell_type, " \u2192 ",
      target_cell_type, " mode reassignments"
    )
  }

  if (is.null(subtitle)) {
    subtitle <- paste0(
      "Union of top ", top_k,
      " genes = ", scales::comma(length(union_genes))
    )
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = ribbon_df,
      ggplot2::aes(x, y, group = flow_id),
      fill = "grey70",
      alpha = 0.65,
      color = NA
    ) +
    ggplot2::geom_rect(
      data = stratum_df,
      ggplot2::aes(
        xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax,
        fill = cluster
      ),
      color = "grey25",
      linewidth = 0.6
    ) +
    ggplot2::scale_fill_manual(
      values = cluster_palette,
      breaks = plot_clusters,
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(1, 2),
      labels = c(source_cell_type, target_cell_type),
      limits = c(0.82, 2.18),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(
      ylim = c(0, total_plot_height),
      clip = "off"
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = NULL, y = NULL,
      fill = "Cluster"
    ) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = if (show_legend) "right" else "none",
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      plot.title.position = "plot",
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    )

  invisible(list(
    plot = p,
    flow_df = flow_df,
    flow_df_all = flow_df_all,
    flow_layout = flow_layout,
    ribbon_df = ribbon_df,
    flow_genes = flow_genes,
    flow_genes_visible = flow_genes_visible,
    flow_genes_pre = flow_genes_pre,
    hidden_low_significance_genes = hidden_low_significance_genes,
    hidden_flow_df = hidden_flow_df,
    union_genes = union_genes,
    df_src_topk = df_src_topk,
    df_tgt_topk = df_tgt_topk,
    df_src_all = df_src_all,
    df_tgt_all = df_tgt_all,
    stratum_df = stratum_df,
    source_totals = source_totals,
    target_totals = target_totals,
    source_scale = source_scale,
    target_scale = target_scale,
    stratum_gap = stratum_gap,
    empty_box_size = empty_box_size,
    genes_with_missing_assignments = genes_with_missing_assignments,
    genes_touching_low_significance = genes_touching_low_significance,
    imputed_missing_assignments_df = imputed_missing_assignments_df,
    dropped_incomplete_genes = dropped_incomplete_genes,
    dropped_low_significance_genes = dropped_low_significance_genes,
    dropped_invalid_cluster_genes = dropped_invalid_cluster_genes,
    low_significance_label = low_significance_label,
    treat_missing_as_low_significance = treat_missing_as_low_significance,
    low_significance_action = low_significance_action
  ))
}

plot_cluster_reassignment_chord_pairwise <- function(
  pc_df,
  compare_cell_types,
  cluster_col = "broad_cluster_ls",
  gene_col = "gene",
  rank_col = "waldStat_1",
  top_k = 300,
  topk_set_mode = c("union", "intersection"),
  rank_desc = TRUE,
  clusters_use = NULL,
  keep_empty_clusters = TRUE,
  empty_sector_size = 0.5,
  cluster_palette = NULL,
  cell_type_palette = NULL,
  chord_color = TRUE,   # FALSE => all grey; TRUE => same-category links get category color
  chord_alpha = 0.55,
  link_lwd = 0.5,
  min_flow = 1,

  # Low-significance / missing-assignment handling
  show_low_significance = TRUE,
  low_significance_label = "low_significance",
  treat_missing_as_low_significance = FALSE,
  hidden_low_significance_mode = c("none", "bar"),
  hidden_bar_col = "grey70",
  hidden_bar_alpha = 0.9,
  hidden_bar_track_height = 0.03,

  big_gap = 14,
  small_gap = 2,
  start_degree = 90,
  cluster_track_height = 0.10,
  celltype_track_height = 0.06,
  ring_gap = 0.01,
  show_legend = TRUE,
  legend_title = "Cluster",
  legend_cex = 0.9,
  legend_ncol = 1,
  legend_inset = c(-0.50, 0),
  legend_mar_right = 16,
  title = NULL,
  return_data = FALSE
) {
  stopifnot(is.data.frame(pc_df))
  stopifnot(all(c("cell_type", cluster_col, gene_col, rank_col) %in% colnames(pc_df)))

  topk_set_mode <- match.arg(topk_set_mode)
  hidden_low_significance_mode <- match.arg(hidden_low_significance_mode)
  compare_cell_types <- as.character(compare_cell_types)

  if (length(compare_cell_types) < 2) {
    stop("compare_cell_types must contain at least two cell types.")
  }

  if (!all(compare_cell_types %in% unique(as.character(pc_df$cell_type)))) {
    stop("All compare_cell_types must exist in pc_df$cell_type.")
  }

  arrange_by_rank <- function(df) {
    ord <- order(df$rank_val, decreasing = isTRUE(rank_desc), na.last = NA)
    df[ord, , drop = FALSE]
  }

  get_gene_annotations <- function(df, cell_type_value) {
    out <- df %>%
      dplyr::filter(.data$cell_type == cell_type_value) %>%
      dplyr::transmute(
        gene = as.character(.data[[gene_col]]),
        cluster_value = as.character(.data[[cluster_col]]),
        rank_val = .data[[rank_col]]
      ) %>%
      dplyr::filter(!is.na(.data$rank_val), is.finite(.data$rank_val))

    out <- arrange_by_rank(out)

    out %>%
      dplyr::distinct(.data$gene, .keep_all = TRUE)
  }

  # ----------------------------
  # Ranked annotations per cell type
  # ----------------------------
  all_ann <- lapply(compare_cell_types, function(ct) get_gene_annotations(pc_df, ct))
  names(all_ann) <- compare_cell_types

  # ----------------------------
  # Full cluster universe used for assignment validity
  # ----------------------------
  observed_clusters_chr <- sort(unique(unlist(lapply(all_ann, function(df) {
    as.character(df$cluster_value)
  }))))

  if (is.null(clusters_use)) {
    all_clusters_chr <- observed_clusters_chr
  } else {
    all_clusters_chr <- as.character(clusters_use)
  }

  # If missing genes are treated as low_significance, make sure that label
  # is considered a valid cluster even if it was not supplied in clusters_use.
  if (isTRUE(treat_missing_as_low_significance) &&
      !(low_significance_label %in% all_clusters_chr)) {
    all_clusters_chr <- c(all_clusters_chr, low_significance_label)
  }

  all_clusters_chr <- unique(all_clusters_chr)

  # Displayed clusters
  if (isTRUE(show_low_significance)) {
    display_clusters_chr <- all_clusters_chr
  } else {
    display_clusters_chr <- setdiff(all_clusters_chr, low_significance_label)
  }

  if (length(display_clusters_chr) == 0) {
    stop("No displayed clusters remain.")
  }

  # ----------------------------
  # Top-k gene sets define global universe
  # ----------------------------
  topk_ann <- lapply(all_ann, function(df) dplyr::slice_head(df, n = top_k))
  topk_gene_sets <- lapply(topk_ann, function(df) as.character(df$gene))

  gene_universe <- switch(
    topk_set_mode,
    union = Reduce(union, topk_gene_sets),
    intersection = Reduce(intersect, topk_gene_sets)
  )

  if (length(gene_universe) == 0) {
    stop("No genes remain in the selected global gene universe.")
  }

  # ----------------------------
  # Use all available annotations for selected universe
  # ----------------------------
  gene_assign_df <- data.frame(
    gene = gene_universe,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  missing_assignment_mat <- matrix(
    FALSE,
    nrow = length(gene_universe),
    ncol = length(compare_cell_types),
    dimnames = list(NULL, compare_cell_types)
  )

  for (ct in compare_cell_types) {
    vec <- stats::setNames(
      as.character(all_ann[[ct]]$cluster_value),
      all_ann[[ct]]$gene
    )

    vals <- unname(vec[gene_assign_df$gene])

    missing_assignment_mat[, ct] <- is.na(vals)
    gene_assign_df[[ct]] <- vals
  }

  genes_with_missing_assignments <- unique(
    gene_assign_df$gene[rowSums(missing_assignment_mat) > 0]
  )

  imputed_missing_assignments_df <- dplyr::bind_rows(
    lapply(compare_cell_types, function(ct) {
      idx <- missing_assignment_mat[, ct]

      if (!any(idx)) {
        return(data.frame(
          gene = character(),
          cell_type = character(),
          imputed_cluster = character(),
          stringsAsFactors = FALSE
        ))
      }

      data.frame(
        gene = gene_assign_df$gene[idx],
        cell_type = ct,
        imputed_cluster = low_significance_label,
        stringsAsFactors = FALSE
      )
    })
  )

  if (isTRUE(treat_missing_as_low_significance)) {
    for (ct in compare_cell_types) {
      idx <- missing_assignment_mat[, ct]
      gene_assign_df[[ct]][idx] <- low_significance_label
    }
  }

  # Keep genes with valid assignments in every compared cell type.
  # If treat_missing_as_low_significance = TRUE, missing assignments have
  # already been replaced by low_significance_label.
  keep_valid <- apply(
    gene_assign_df[, compare_cell_types, drop = FALSE],
    1,
    function(x) all(!is.na(x) & x %in% all_clusters_chr)
  )

  if (isTRUE(treat_missing_as_low_significance)) {
    dropped_incomplete_genes <- character(0)
  } else {
    dropped_incomplete_genes <- genes_with_missing_assignments
  }

  dropped_invalid_cluster_genes <- setdiff(
    gene_assign_df$gene[!keep_valid],
    dropped_incomplete_genes
  )

  gene_assign_df <- gene_assign_df[keep_valid, , drop = FALSE]

  if (nrow(gene_assign_df) == 0) {
    if (isTRUE(treat_missing_as_low_significance)) {
      stop("No genes have valid assignments across compare_cell_types after imputing missing assignments as low_significance.")
    } else {
      stop("No genes have complete valid assignments across all compare_cell_types.")
    }
  }

  n_ct <- length(compare_cell_types)

  # ----------------------------
  # Sector sizes for displayed clusters only
  # ----------------------------
  sector_df <- expand.grid(
    cell_type = compare_cell_types,
    cluster   = display_clusters_chr,
    stringsAsFactors = FALSE,
    KEEP.OUT.ATTRS = FALSE
  )

  sector_df$sector <- paste0(sector_df$cell_type, "::", sector_df$cluster)

  sector_df$unique_n <- vapply(seq_len(nrow(sector_df)), function(i) {
    ct <- sector_df$cell_type[i]
    cl <- sector_df$cluster[i]
    sum(gene_assign_df[[ct]] == cl)
  }, numeric(1))

  if (isTRUE(keep_empty_clusters)) {
    sector_df$block_n <- ifelse(
      sector_df$unique_n > 0,
      sector_df$unique_n,
      empty_sector_size
    )
    sector_df$plot_n <- sector_df$block_n * (n_ct - 1)
  } else {
    sector_df <- sector_df[sector_df$unique_n > 0, , drop = FALSE]
    sector_df$block_n <- sector_df$unique_n
    sector_df$plot_n <- sector_df$block_n * (n_ct - 1)
  }

  if (nrow(sector_df) == 0) {
    stop("No sectors remain after filtering.")
  }

  # Sector order: grouped by cell type, then cluster
  sector_order_all <- unlist(lapply(compare_cell_types, function(ct) {
    paste0(ct, "::", display_clusters_chr)
  }), use.names = FALSE)

  sector_order <- sector_order_all[sector_order_all %in% sector_df$sector]
  sector_df <- sector_df[match(sector_order, sector_df$sector), , drop = FALSE]

  # ----------------------------
  # Palettes
  # ----------------------------
  if (is.null(cluster_palette)) {
    pal_vals <- grDevices::hcl.colors(length(all_clusters_chr), palette = "Dark 3")
    cluster_palette <- stats::setNames(pal_vals, all_clusters_chr)
  } else {
    if (is.list(cluster_palette)) {
      cluster_palette <- unlist(cluster_palette)
    }

    missing_palette_clusters <- setdiff(all_clusters_chr, names(cluster_palette))

    # Be permissive if the only missing color is the low-significance label.
    # This matters when missing assignments are imputed as low_significance
    # but that label is hidden and represented by hidden_bar_col.
    if (length(missing_palette_clusters) > 0) {
      if (all(missing_palette_clusters %in% low_significance_label)) {
        cluster_palette <- c(
          cluster_palette,
          stats::setNames(hidden_bar_col, low_significance_label)
        )
      } else {
        stop("cluster_palette must be a named vector containing all clusters in clusters_use/data.")
      }
    }

    cluster_palette <- cluster_palette[all_clusters_chr]
  }

  if (is.null(cell_type_palette)) {
    ct_cols <- grDevices::hcl.colors(length(compare_cell_types), palette = "Set 2")
    cell_type_palette <- stats::setNames(ct_cols, compare_cell_types)
  } else {
    if (is.list(cell_type_palette)) {
      cell_type_palette <- unlist(cell_type_palette)
    }

    if (!all(compare_cell_types %in% names(cell_type_palette))) {
      stop("cell_type_palette must be a named vector containing all compare_cell_types.")
    }

    cell_type_palette <- cell_type_palette[compare_cell_types]
  }

  grid.col <- stats::setNames(
    unname(cluster_palette[sub("^.*::", "", sector_df$sector)]),
    sector_df$sector
  )

  # ----------------------------
  # Split each displayed sector into one block per OTHER cell type
  # ----------------------------
  block_rows <- list()
  block_idx  <- 1L

  for (i in seq_len(nrow(sector_df))) {
    ct      <- sector_df$cell_type[i]
    sec     <- sector_df$sector[i]
    cl      <- sector_df$cluster[i]
    block_n <- sector_df$block_n[i]
    partners <- compare_cell_types[compare_cell_types != ct]

    x_left <- 0

    for (partner in partners) {
      block_rows[[block_idx]] <- data.frame(
        sector    = sec,
        cell_type = ct,
        cluster   = cl,
        partner   = partner,
        x1        = x_left,
        x2        = x_left + block_n,
        block_n   = block_n,
        stringsAsFactors = FALSE
      )

      x_left <- x_left + block_n
      block_idx <- block_idx + 1L
    }
  }

  block_df <- dplyr::bind_rows(block_rows)

  # Running offsets inside each sector-partner block
  block_offset <- stats::setNames(
    block_df$x1,
    paste(block_df$sector, block_df$partner, sep = "||")
  )

  # ----------------------------
  # Build links
  # ----------------------------
  ct_pairs <- utils::combn(compare_cell_types, 2, simplify = FALSE)

  link_rows <- list()
  stub_rows <- list()
  link_idx  <- 1L
  stub_idx  <- 1L

  for (pair in ct_pairs) {
    ct1 <- pair[1]
    ct2 <- pair[2]

    for (cl1 in all_clusters_chr) {
      for (cl2 in all_clusters_chr) {
        n_ij <- sum(gene_assign_df[[ct1]] == cl1 & gene_assign_df[[ct2]] == cl2)
        if (n_ij <= 0) next

        vis1 <- cl1 %in% display_clusters_chr
        vis2 <- cl2 %in% display_clusters_chr

        # Case 1: both endpoints displayed -> normal chord
        if (vis1 && vis2) {
          sec1 <- paste0(ct1, "::", cl1)
          sec2 <- paste0(ct2, "::", cl2)

          key1 <- paste(sec1, ct2, sep = "||")
          key2 <- paste(sec2, ct1, sep = "||")

          x1 <- block_offset[[key1]]
          x2 <- x1 + n_ij
          block_offset[[key1]] <- x2

          y1 <- block_offset[[key2]]
          y2 <- y1 + n_ij
          block_offset[[key2]] <- y2

          if (isTRUE(chord_color) &&
              identical(cl1, cl2) &&
              !identical(cl1, low_significance_label)) {
            col_ij <- unname(cluster_palette[cl1])
          } else {
            col_ij <- "grey70"
          }

          link_rows[[link_idx]] <- data.frame(
            cell_type_1 = ct1,
            cluster_1   = cl1,
            sector1     = sec1,
            x1          = x1,
            x2          = x2,
            cell_type_2 = ct2,
            cluster_2   = cl2,
            sector2     = sec2,
            y1          = y1,
            y2          = y2,
            n           = n_ij,
            visible     = n_ij >= min_flow,
            col         = col_ij,
            stringsAsFactors = FALSE
          )

          link_idx <- link_idx + 1L
        }

        # Case 2: visible -> hidden low_significance
        else if (!show_low_significance &&
                 vis1 &&
                 identical(cl2, low_significance_label)) {
          sec1 <- paste0(ct1, "::", cl1)
          key1 <- paste(sec1, ct2, sep = "||")

          x1 <- block_offset[[key1]]
          x2 <- x1 + n_ij
          block_offset[[key1]] <- x2

          if (hidden_low_significance_mode == "bar") {
            stub_rows[[stub_idx]] <- data.frame(
              sector = sec1,
              x1 = x1,
              x2 = x2,
              n = n_ij,
              cell_type = ct1,
              cluster = cl1,
              partner_cell_type = ct2,
              hidden_cluster = cl2,
              stringsAsFactors = FALSE
            )

            stub_idx <- stub_idx + 1L
          }
        }

        # Case 3: hidden low_significance -> visible
        else if (!show_low_significance &&
                 identical(cl1, low_significance_label) &&
                 vis2) {
          sec2 <- paste0(ct2, "::", cl2)
          key2 <- paste(sec2, ct1, sep = "||")

          y1 <- block_offset[[key2]]
          y2 <- y1 + n_ij
          block_offset[[key2]] <- y2

          if (hidden_low_significance_mode == "bar") {
            stub_rows[[stub_idx]] <- data.frame(
              sector = sec2,
              x1 = y1,
              x2 = y2,
              n = n_ij,
              cell_type = ct2,
              cluster = cl2,
              partner_cell_type = ct1,
              hidden_cluster = cl1,
              stringsAsFactors = FALSE
            )

            stub_idx <- stub_idx + 1L
          }
        }

        # Case 4: hidden -> hidden, or any non-displayed category pair
        else {
          next
        }
      }
    }
  }

  link_df <- if (length(link_rows) > 0) {
    dplyr::bind_rows(link_rows)
  } else {
    data.frame(
      cell_type_1 = character(),
      cluster_1   = character(),
      sector1     = character(),
      x1          = numeric(),
      x2          = numeric(),
      cell_type_2 = character(),
      cluster_2   = character(),
      sector2     = character(),
      y1          = numeric(),
      y2          = numeric(),
      n           = numeric(),
      visible     = logical(),
      col         = character(),
      stringsAsFactors = FALSE
    )
  }

  hidden_bar_df <- if (length(stub_rows) > 0) {
    dplyr::bind_rows(stub_rows)
  } else {
    data.frame(
      sector = character(),
      x1 = numeric(),
      x2 = numeric(),
      n = numeric(),
      cell_type = character(),
      cluster = character(),
      partner_cell_type = character(),
      hidden_cluster = character(),
      stringsAsFactors = FALSE
    )
  }

  if (nrow(link_df) > 0) {
    link_df <- link_df[order(link_df$visible, link_df$n), , drop = FALSE]
  }

  # ----------------------------
  # Circos layout
  # ----------------------------
  gap.after <- unlist(lapply(compare_cell_types, function(ct) {
    n_ct_sectors <- sum(sector_df$cell_type == ct)
    c(rep(small_gap, max(n_ct_sectors - 1, 0)), big_gap)
  }), use.names = FALSE)

  if (is.null(title)) {
    title <- ""
  }

  link_lwd_scalar <- as.numeric(link_lwd)[1]

  if (!is.finite(link_lwd_scalar) || link_lwd_scalar < 0) {
    link_lwd_scalar <- 0.5
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  if (isTRUE(show_legend)) {
    par(mar = c(1, 1, 3, legend_mar_right), xpd = NA)
  } else {
    par(mar = c(1, 1, 3, 1), xpd = NA)
  }

  circlize::circos.clear()
  on.exit(circlize::circos.clear(), add = TRUE)

  circlize::circos.par(
    start.degree = start_degree,
    gap.after = gap.after,
    cell.padding = c(0, 0, 0, 0),
    track.margin = c(ring_gap, ring_gap),
    points.overflow.warning = FALSE
  )

  circlize::circos.initialize(
    factors = factor(sector_df$sector, levels = sector_df$sector),
    xlim = cbind(rep(0, nrow(sector_df)), sector_df$plot_n)
  )

  # ----------------------------
  # Outer cell-type band
  # ----------------------------
  circlize::circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = celltype_track_height,
    bg.border = NA,
    bg.col = NA,
    panel.fun = function(x, y) {}
  )

  for (ct in compare_cell_types) {
    ct_sectors <- sector_df$sector[sector_df$cell_type == ct]
    if (length(ct_sectors) == 0) next

    circlize::highlight.sector(
      sector.index = ct_sectors,
      track.index  = 1,
      col          = grDevices::adjustcolor(cell_type_palette[ct], alpha.f = 0.95),
      border       = NA,
      text         = ct,
      cex          = 1.0,
      text.col     = "black",
      niceFacing   = TRUE,
      font         = 2
    )
  }

  # ----------------------------
  # Cluster band
  # ----------------------------
  circlize::circos.trackPlotRegion(
    ylim = c(0, 1),
    track.height = cluster_track_height,
    bg.border = NA,
    bg.col = NA,
    panel.fun = function(x, y) {
      sec <- circlize::get.cell.meta.data("sector.index")
      xl  <- circlize::get.cell.meta.data("xlim")

      circlize::circos.rect(
        xleft  = xl[1],
        ybottom = 0,
        xright = xl[2],
        ytop    = 1,
        col    = grid.col[sec],
        border = "grey30",
        lwd    = 0.6
      )
    }
  )

  # ----------------------------
  # Optional inner gray bars for hidden low_significance flows
  # ----------------------------
  if (!show_low_significance &&
      hidden_low_significance_mode == "bar" &&
      nrow(hidden_bar_df) > 0 &&
      hidden_bar_track_height > 0) {

    circlize::circos.trackPlotRegion(
      ylim = c(0, 1),
      track.height = hidden_bar_track_height,
      bg.border = NA,
      bg.col = NA,
      panel.fun = function(x, y) {
        sec <- circlize::get.cell.meta.data("sector.index")
        rows <- hidden_bar_df[hidden_bar_df$sector == sec, , drop = FALSE]

        if (nrow(rows) == 0) return(NULL)

        for (j in seq_len(nrow(rows))) {
          circlize::circos.rect(
            xleft  = rows$x1[j],
            ybottom = 0,
            xright = rows$x2[j],
            ytop   = 1,
            col    = grDevices::adjustcolor(hidden_bar_col, alpha.f = hidden_bar_alpha),
            border = NA
          )
        }
      }
    )
  }

  # ----------------------------
  # Draw links
  # ----------------------------
  if (nrow(link_df) > 0) {
    link_df_draw <- link_df[link_df$visible, , drop = FALSE]

    for (i in seq_len(nrow(link_df_draw))) {
      circlize::circos.link(
        sector.index1 = link_df_draw$sector1[i],
        point1        = c(link_df_draw$x1[i], link_df_draw$x2[i]),
        sector.index2 = link_df_draw$sector2[i],
        point2        = c(link_df_draw$y1[i], link_df_draw$y2[i]),
        col           = grDevices::adjustcolor(link_df_draw$col[i], alpha.f = chord_alpha),
        border        = NA,
        lwd           = link_lwd_scalar
      )
    }
  }

  title(main = title)

  # ----------------------------
  # Legend
  # ----------------------------
  if (isTRUE(show_legend)) {
    if (!show_low_significance && hidden_low_significance_mode == "bar") {
      legend_labels <- c(display_clusters_chr, low_significance_label)
      legend_fills <- c(
        unname(cluster_palette[display_clusters_chr]),
        hidden_bar_col
      )
    } else {
      legend_labels <- display_clusters_chr
      legend_fills <- unname(cluster_palette[legend_labels])
    }

    legend(
      "right",
      inset  = legend_inset,
      legend = legend_labels,
      fill   = legend_fills,
      border = "grey30",
      bty    = "n",
      title  = legend_title,
      cex    = legend_cex,
      ncol   = legend_ncol,
      xpd    = NA
    )
  }

  # ----------------------------
  # Directed summary table
  # ----------------------------
  directed_link_df <- dplyr::bind_rows(
    link_df %>%
      dplyr::transmute(
        src_cell_type = cell_type_1,
        src_cluster   = cluster_1,
        dst_cell_type = cell_type_2,
        dst_cluster   = cluster_2,
        n
      ),
    link_df %>%
      dplyr::transmute(
        src_cell_type = cell_type_2,
        src_cluster   = cluster_2,
        dst_cell_type = cell_type_1,
        dst_cluster   = cluster_1,
        n
      )
  )

  out <- list(
    gene_universe = gene_universe,
    dropped_incomplete_genes = dropped_incomplete_genes,
    genes_with_missing_assignments = genes_with_missing_assignments,
    imputed_missing_assignments_df = imputed_missing_assignments_df,
    dropped_invalid_cluster_genes = dropped_invalid_cluster_genes,
    gene_assign_df = gene_assign_df,
    sector_df = sector_df,
    block_df = block_df,
    link_df = link_df,
    hidden_bar_df = hidden_bar_df,
    directed_link_df = directed_link_df,
    compare_cell_types = compare_cell_types,
    topk_set_mode = topk_set_mode,
    cluster_palette = cluster_palette,
    cell_type_palette = cell_type_palette,
    chord_color = chord_color,
    show_low_significance = show_low_significance,
    low_significance_label = low_significance_label,
    treat_missing_as_low_significance = treat_missing_as_low_significance,
    hidden_low_significance_mode = hidden_low_significance_mode
  )

  if (isTRUE(return_data)) {
    return(out)
  } else {
    invisible(out)
  }
}

plot_cluster_confusion <- function(
  df,
  cluster_x,
  cluster_y,
  normalize = c("none", "row", "column")
) {
  normalize <- match.arg(normalize)

  tab <- table(df[[cluster_y]], df[[cluster_x]])

  if (normalize == "row") {
    tab <- prop.table(tab, margin = 1)
  } else if (normalize == "column") {
    tab <- prop.table(tab, margin = 2)
  }

  plot_df <- as.data.frame(tab)
  colnames(plot_df) <- c("cluster_y", "cluster_x", "value")

  ggplot(plot_df, aes(x = cluster_x, y = cluster_y, fill = value)) +
    geom_tile() +
    geom_text(
      aes(label = if (normalize == "none") value else sprintf("%.2f", value)),
      size = 3
    ) +
    scale_fill_viridis_c() +
    labs(
      x = cluster_x,
      y = cluster_y,
      fill = if (normalize == "none") "Count" else "Proportion"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      axis.ticks = element_blank()
    )
}

plot_pca_reconstruction <- function(gene, X_mat, pca_X, time = tt, n_pcs = c(3, 4)) {
  stopifnot(gene %in% rownames(X_mat))

  gene_idx <- match(gene, rownames(X_mat))
  original <- X_mat[gene_idx, ]

  reconstruct_k <- function(K) {
    scores   <- pca_X$x[gene_idx, seq_len(K), drop = FALSE]
    loadings <- pca_X$rotation[, seq_len(K), drop = FALSE]

    reconstructed <- as.numeric(scores %*% t(loadings))

    if (!is.null(pca_X$scale) && !identical(pca_X$scale, FALSE)) {
      reconstructed <- reconstructed * pca_X$scale
    }

    if (!is.null(pca_X$center) && !identical(pca_X$center, FALSE)) {
      reconstructed <- reconstructed + pca_X$center
    }

    reconstructed
  }

  plot_df <- data.frame(time = time, Original = original)

  for (K in n_pcs) {
    plot_df[[paste0("PC1-", K)]] <- reconstruct_k(K)
  }

  plot_df_long <- plot_df %>%
    tidyr::pivot_longer(
      cols = -time,
      names_to = "trajectory",
      values_to = "value"
    )

  ggplot(plot_df_long, aes(x = time, y = value, linetype = trajectory)) +
    geom_line(linewidth = 1.1) +
    labs(
      title = paste0(gene, ": PCA trajectory reconstruction"),
      x = "Time",
      y = "Normalized trajectory",
      linetype = NULL
    ) +
    theme_classic()
}


# ---------------------------------------------------------------------------
# Python (cp.pl) entry points. Each reads its inputs from CSV files written by
# the wrapper, draws with the helpers above and saves to every path in `output`
# (the requested file plus a PNG preview for inline display).
# ---------------------------------------------------------------------------

# cluster::silhouette() widths; same definition when the 'cluster' package is not installed
# (a = mean distance to the rest of the own cluster, b = smallest mean distance to another
# cluster, s = (b - a) / max(a, b), 0 for singletons).
chronopack_silhouette <- function(cl, D) {
  if (requireNamespace("cluster", quietly = TRUE)) {
    return(cluster::silhouette(cl, D))
  }
  Dm <- as.matrix(D)
  ids <- sort(unique(cl))
  idx <- match(cl, ids)
  n <- length(idx)
  sizes <- tabulate(idx, nbins = length(ids))
  sums <- t(rowsum(t(Dm), group = idx))
  own <- cbind(seq_len(n), idx)
  a <- sums[own] / pmax(sizes[idx] - 1, 1)
  means <- sweep(sums, 2, sizes, "/")
  means[own] <- Inf
  b <- if (length(ids) > 1) apply(means, 1, min) else rep(NA_real_, n)
  neighbor <- if (length(ids) > 1) ids[apply(means, 1, which.min)] else rep(NA, n)
  denom <- pmax(a, b)
  s <- ifelse(is.finite(denom) & denom > 0, (b - a) / denom, 0)
  s[sizes[idx] == 1] <- 0
  cbind(cluster = cl, neighbor = neighbor, sil_width = s)
}

# The TESTING-derived helpers call ggplot2, dplyr and patchwork functions unqualified.
chronopack_attach_plot_packages <- function() {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(patchwork)
  })
}

chronopack_read_plot_table <- function(path) {
  df <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = "")
  if ("row_id" %in% colnames(df)) rownames(df) <- df$row_id
  df
}

chronopack_open_device <- function(output, width, height, dpi = 300) {
  output <- path.expand(output)
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  ext <- tolower(tools::file_ext(output))
  if (ext == "pdf") {
    grDevices::pdf(output, width = width, height = height)
  } else if (ext == "svg") {
    grDevices::svg(output, width = width, height = height)
  } else if (ext == "png") {
    grDevices::png(output, width = width, height = height, units = "in", res = dpi)
  } else {
    stop("Unsupported plot output extension: ", ext, ". Use pdf, svg, or png.")
  }
}

# `plot` is a ggplot / patchwork object, or a function drawing with base graphics
# (its value from the last draw is returned).
chronopack_save_plot <- function(plot, output, width = 7, height = 5, dpi = 300) {
  value <- NULL
  for (path in output) {
    chronopack_open_device(path, width, height, dpi)
    tryCatch(
      if (is.function(plot)) value <- plot() else print(plot),
      finally = grDevices::dev.off()
    )
  }
  invisible(value)
}

# Evaluate `expr` with a null device open, for helpers that print their plot themselves.
chronopack_without_display <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expr
}

chronopack_mode_palette <- function(levels) {
  base <- c(
    "impulse-up_1" = "#F8766D",
    "impulse-up_2" = "#D84A42",
    "trend-up" = "#C29A33",
    "transient-down" = "#8B78C6",
    "trend-down_1" = "#619CFF",
    "trend-down_2" = "#3C6FB6"
  )
  levels <- as.character(levels)
  missing <- setdiff(levels, names(base))
  if (length(missing) > 0) {
    extra <- grDevices::hcl.colors(length(missing), palette = "Dark 3")
    names(extra) <- missing
    base <- c(base, extra)
  }
  base[levels]
}

chronopack_plot_gene_trajectory_csv <- function(
  data_path,
  output,
  width = 9,
  height = 6,
  title = "",
  normalize = FALSE
) {
  df <- utils::read.csv(data_path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("row_id", "gene", "cell_type", "time", "value")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) stop("Missing required columns: ", paste(missing, collapse = ", "))

  df$label <- ifelse(
    is.na(df$cell_type) | df$cell_type == "",
    df$gene,
    paste0(df$gene, " (", df$cell_type, ")")
  )
  if (!"broad_cluster" %in% colnames(df)) df$broad_cluster <- "trajectory"
  df$broad_cluster <- as.character(df$broad_cluster)
  pal <- chronopack_mode_palette(sort(unique(df$broad_cluster)))
  title_value <- if (nzchar(title)) title else NULL

  p <- ggplot2::ggplot(df, ggplot2::aes(time, value, color = broad_cluster, group = row_id)) +
    ggplot2::geom_line(linewidth = 0.9, alpha = 0.9) +
    ggplot2::facet_wrap(stats::as.formula("~ label"), scales = "free_y", ncol = 3) +
    ggplot2::scale_color_manual(values = pal) +
    ggplot2::labs(
      title = title_value,
      x = "Aligned time (min)",
      y = ifelse(isTRUE(normalize), "z trajectory", "Trajectory"),
      color = NULL
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 8),
      legend.position = "bottom"
    )

  chronopack_save_plot(p, output, width = width, height = height)
}

# cp.pl.trajectory_features: feature_plot_pc_simple() on a pc_df-style table.
chronopack_plot_trajectory_features_csv <- function(data_path, output, feature, width = 6.5, height = 5, dpi = 300, ...) {
  chronopack_attach_plot_packages()
  pc_df <- chronopack_read_plot_table(data_path)
  p <- feature_plot_pc_simple(pc_df, feature = feature, ...)
  chronopack_save_plot(p, output, width = width, height = height, dpi = dpi)
}

# cp.pl.cluster_medoids_and_top_genes. As in TESTING, traj_mat and X are the row z-scored
# trajectories and D the Euclidean distances of the joint clustering features.
chronopack_plot_cluster_medoids_csv <- function(
  traj_path, features_path, obs_path, time_path, output,
  width = 13, height = 7, dpi = 300, ...
) {
  chronopack_attach_plot_packages()
  traj_mat <- as.matrix(utils::read.csv(traj_path, row.names = 1, check.names = FALSE))
  X_use <- as.matrix(utils::read.csv(features_path, row.names = 1, check.names = FALSE))[rownames(traj_mat), , drop = FALSE]
  obs <- utils::read.csv(obs_path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = "")
  obs <- obs[match(rownames(traj_mat), obs$row_id), , drop = FALSE]
  tt <- utils::read.csv(time_path)$time

  p <- plot_cluster_medoids_and_top_genes(
    traj_mat = traj_mat,
    X = traj_mat,
    tt = tt,
    cl = obs$cluster,
    D = stats::dist(X_use, method = "euclidean"),
    gene_score_map = stats::setNames(obs$score, obs$row_id),
    ...
  )
  chronopack_save_plot(p, output, width = width, height = height, dpi = dpi)
}

# cp.pl.mode_proportions: the pie charts of TESTING.ipynb cell 82 (pie_df is computed in Python).
chronopack_plot_mode_proportions_csv <- function(data_path, output, cluster_order, cluster_colors,
                                                 width = 12, height = 5, dpi = 300) {
  chronopack_attach_plot_packages()
  pie_df <- chronopack_read_plot_table(data_path)
  pie_df$broad_cluster <- factor(pie_df$broad_cluster, levels = cluster_order)

  p_pie <- ggplot(pie_df, aes(x = "", y = prop, fill = broad_cluster)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y", direction = -1) +
    facet_wrap(~cell_type, nrow = 1) +
    theme_void(base_size = 16) +
    theme(
      strip.text = element_text(face = "bold", size = 18),
      plot.title = element_text(hjust = 0.5, size = 20),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    ) +
    labs(
      title = "",
      fill = "broad_cluster"
    ) +
    scale_fill_manual(
      breaks = cluster_order,
      drop = FALSE,
      values = cluster_colors
    )

  chronopack_save_plot(p_pie, output, width = width, height = height, dpi = dpi)
}

# cp.pl.mode_sankey: TESTING.ipynb cell 90, one Sankey per context pair (legend on the last),
# combined with patchwork. Returns the flow tables.
chronopack_plot_mode_sankey_csv <- function(data_path, output, source_cell_types, target_cell_types,
                                            width = 11, height = 5, dpi = 300, ...) {
  chronopack_attach_plot_packages()
  pc_df_subset <- chronopack_read_plot_table(data_path)
  n <- length(source_cell_types)
  res <- lapply(seq_len(n), function(i) {
    plot_cluster_reassignment_sankey(
      pc_df = pc_df_subset,
      source_cell_type = source_cell_types[[i]],
      target_cell_type = target_cell_types[[i]],
      show_legend = i == n,
      ...
    )
  })
  p_all <- Reduce(`+`, lapply(res, function(x) x$plot + ggtitle("")))
  chronopack_save_plot(p_all, output, width = width, height = height, dpi = dpi)
  dplyr::bind_rows(lapply(seq_len(n), function(i) {
    data.frame(
      source_cell_type = source_cell_types[[i]],
      target_cell_type = target_cell_types[[i]],
      as.data.frame(res[[i]]$flow_df),
      check.names = FALSE
    )
  }))
}

# cp.pl.mode_chord: TESTING.ipynb cell 102 (circlize, base graphics). Returns the pairwise links.
chronopack_plot_mode_chord_csv <- function(data_path, output, width = 8, height = 6, dpi = 300, ...) {
  chronopack_attach_plot_packages()
  pc_df_subset <- chronopack_read_plot_table(data_path)
  res <- chronopack_save_plot(
    function() plot_cluster_reassignment_chord_pairwise(pc_df = pc_df_subset, return_data = TRUE, ...),
    output, width = width, height = height, dpi = dpi
  )
  as.data.frame(res$link_df)
}

# cp.pl.cluster_pair_representative_genes: TESTING.ipynb cell 107; `sce_paths` names the
# tradeSeq SCE (.rds) of each context. Returns the chosen gene per cluster pair.
chronopack_plot_cluster_pair_genes_csv <- function(data_path, output, sce_paths,
                                                   width = 10, height = 10, dpi = 300, ...) {
  chronopack_attach_plot_packages()
  suppressPackageStartupMessages(requireNamespace("SingleCellExperiment"))
  pc_df <- chronopack_read_plot_table(data_path)
  sce.ls <- lapply(as.list(sce_paths), readRDS)
  res <- chronopack_without_display(
    plot_cluster_pair_representative_genes(pc_df = pc_df, sce_ls = sce.ls, ...)
  )
  chronopack_save_plot(res$plot, output, width = width, height = height, dpi = dpi)
  as.data.frame(res$chosen_genes)
}

# cp.pl.cluster_violin: TESTING.ipynb cell 111, one plot_cluster_violin() panel per feature,
# combined with patchwork. Returns the per-group summaries.
chronopack_plot_cluster_violin_csv <- function(data_path, output, features, plot_titles,
                                               width = 7, height = 3.5, dpi = 300, ...) {
  chronopack_attach_plot_packages()
  pc_df_subset <- chronopack_read_plot_table(data_path)
  res <- chronopack_without_display(lapply(seq_along(features), function(i) {
    plot_cluster_violin(pc_df_subset, feature = features[[i]], plot_title = plot_titles[[i]], ...)
  }))
  fig <- Reduce(`+`, lapply(res, `[[`, "plot"))
  chronopack_save_plot(fig, output, width = width, height = height, dpi = dpi)
  dplyr::bind_rows(lapply(seq_along(features), function(i) {
    data.frame(feature = features[[i]], as.data.frame(res[[i]]$summary), check.names = FALSE)
  }))
}

# Volcano of a single-lineage tradeSeq association table (TradeSeqFits.ipynb): waldStat_1 on a
# two-slope log10 axis (compressed below y_break) against signedMaxLog2FC_start, coloured by qval_1.
plot_association_volcano <- function(
  asso_df,
  fdr_thresh        = 0.01,
  l2fc_thresh       = 0.5,       # magnitude cutoff for coloring & guide lines
  top_labels        = 12,        # number of labels
  point_alpha       = 0.85,
  point_size        = 2.2,
  axis_title_size   = 26,
  axis_text_size    = 18,
  label_size        = 4.6,
  label_only_pass   = FALSE,
  # limits & flooring
  y_log_floor       = 1e-6,      # floor on Wald to keep log defined
  x_lim             = NULL,      # c(xmin, xmax) in log2FC units
  y_lim             = NULL,      # c(ymin, ymax) in Wald units (pre-log)
  # --- two-slope scaling controls (in Wald/log10 space) ---
  y_break           = NULL,      # break on Wald (pre-log); NULL = auto from non-sig
  y_break_quantile  = 0.95,      # quantile of non-sig Wald to set break if NULL
  y_lower_scale     = 0.35,      # slope below break (<=1; 1 = no compression)
  # --- label placement: "stack" (labels stacked vertically near their points, then
  # overlaps nudged apart with ggrepel) or "repel" (ggrepel layout from scratch) ---
  label_method        = c("stack", "repel"),
  stack_resolve_overlaps = TRUE,   # "stack": ggrepel pass starting from the stacked positions
  stack_box_padding   = 0.12,
  stack_force         = 1,
  repel_box_padding   = 0.6,
  repel_point_padding = 0.15,
  repel_force         = 4,
  repel_nudge_x       = 0.25,
  repel_seed          = 1,
  # --- "stack" placement ---
  label_x_nudge       = 0.05,
  label_y_offset_frac = 0.05,
  label_min_vgap_frac = 0.012,
  label_segment_alpha = 0.6,
  label_segment_size  = 0.4,
  # --- how to choose top-k labels ---
  label_rank_by     = c("wald","logfc"),
  label_logfc_abs   = TRUE,
  save_prefix       = NULL
) {
  label_rank_by <- match.arg(label_rank_by)
  label_method <- match.arg(label_method)

  stopifnot(all(c("signedMaxLog2FC_start","qval_1","waldStat_1") %in% names(asso_df)))
  df <- asso_df

  # ---- raw axes ----
  x_raw <- df$signedMaxLog2FC_start
  y_raw <- df$waldStat_1
  q_raw <- df$qval_1

  # ---- winsorize / clamp to limits so points are kept ----
  y_plot <- ifelse(is.finite(y_raw), pmax(y_raw, y_log_floor), y_log_floor)
  if (!is.null(y_lim) && length(y_lim) == 2) {
    if (is.finite(y_lim[1])) y_plot <- pmax(y_plot, max(y_log_floor, y_lim[1]))
    if (is.finite(y_lim[2])) y_plot <- pmin(y_plot, y_lim[2])
  }

  x_plot <- x_raw
  if (!is.null(x_lim) && length(x_lim) == 2 && all(is.finite(x_lim))) {
    xmin <- min(x_lim); xmax <- max(x_lim)
    msk <- is.finite(x_plot)
    x_plot[msk] <- pmin(pmax(x_plot[msk], xmin), xmax)
  }

  df$x_plot <- x_plot
  df$y_plot <- y_plot

  # ---- color categories: other / sig_mid / neg / pos ----
  pass_fdr <- is.finite(q_raw) & (q_raw < fdr_thresh)
  n <- nrow(df)
  side <- rep("other", n)

  if (l2fc_thresh > 0) {
    neg_msk <- pass_fdr & is.finite(x_raw) & (x_raw <= -l2fc_thresh)
    pos_msk <- pass_fdr & is.finite(x_raw) & (x_raw >=  l2fc_thresh)
    mid_msk <- pass_fdr & is.finite(x_raw) & (abs(x_raw) < l2fc_thresh)  # NEW darker gray
    side[neg_msk] <- "neg"
    side[pos_msk] <- "pos"
    side[mid_msk] <- "sig_mid"
  } else {
    side[ pass_fdr & is.finite(x_raw) & (x_raw < 0) ] <- "neg"
    side[ pass_fdr & is.finite(x_raw) & (x_raw > 0) ] <- "pos"
    # no mid band when no magnitude cutoff
  }
  df$side <- factor(side, levels = c("other","sig_mid","neg","pos"))

  # ---- labels: top-k by chosen metric ----
  df$label <- ""
  if (top_labels > 0) {
    rank_val <- if (label_rank_by == "wald") y_raw else if (label_logfc_abs) abs(x_raw) else x_raw
    cand <- which(is.finite(rank_val))
    if (label_only_pass) cand <- cand[df$side[cand] %in% c("neg","pos")]
    if (length(cand) > 0) {
      tiebreak <- if (label_rank_by == "wald") abs(x_raw) else y_raw
      ord <- cand[order(rank_val[cand], tiebreak[cand], decreasing = TRUE)]
      labs_idx <- head(ord, top_labels)
      rn <- rownames(df)
      if (!is.null(rn))               df$label[labs_idx] <- rn[labs_idx]
      else if ("gene" %in% names(df)) df$label[labs_idx] <- as.character(df$gene[labs_idx])
    }
  }
  label_df <- df[df$label != "" & is.finite(df$x_plot) & is.finite(df$y_plot), , drop = FALSE]

  # ---- piecewise scaling in log10 space ----
  logy <- log10(df$y_plot)

  if (is.null(y_break) || !is.finite(y_break)) {
    non_sig <- df$y_plot[!pass_fdr & is.finite(df$y_plot)]
    if (length(non_sig)) {
      y_break <- as.numeric(stats::quantile(non_sig, probs = y_break_quantile, na.rm = TRUE))
    } else {
      y_break <- as.numeric(stats::quantile(df$y_plot, probs = 0.60, na.rm = TRUE))
    }
  }
  y_break <- max(y_break, y_log_floor)
  t_break <- log10(y_break)
  y_lower_scale <- max(1e-6, min(1, y_lower_scale))

  ps <- function(t) ifelse(t <= t_break, y_lower_scale * t, y_lower_scale * t_break + (t - t_break))
  df$y_scaled <- ps(logy)

  # ---- axis limits & ticks ----
  t_lo <- if (!is.null(y_lim) && is.finite(y_lim[1])) log10(max(y_log_floor, y_lim[1])) else min(logy[is.finite(logy)])
  t_hi <- if (!is.null(y_lim) && is.finite(y_lim[2])) log10(y_lim[2])                 else max(logy[is.finite(logy)])
  y_limits_scaled <- ps(c(t_lo, t_hi))
  tick_t   <- seq(floor(t_lo), ceiling(t_hi), by = 1)
  tick_pos <- ps(tick_t)
  tick_lab <- as.character(tick_t)

  # ---- compute the Wald value closest to the FDR cutoff & map to scaled y ----
  fdr_line_scaled <- NA_real_
  m_fdr <- is.finite(q_raw) & is.finite(y_raw)
  if (any(m_fdr)) {
    idx <- which.min(abs(q_raw[m_fdr] - fdr_thresh))
    y_fdr <- y_raw[m_fdr][idx]
    # clamp to plotting domain, then map to scaled coordinates
    y_fdr <- max(y_fdr, y_log_floor)
    if (!is.null(y_lim) && length(y_lim) == 2) {
      if (is.finite(y_lim[1])) y_fdr <- max(y_fdr, y_lim[1])
      if (is.finite(y_lim[2])) y_fdr <- min(y_fdr, y_lim[2])
    }
    fdr_line_scaled <- ps(log10(y_fdr))
    # clip inside final scale limits
    if (is.finite(y_limits_scaled[1])) fdr_line_scaled <- max(fdr_line_scaled, y_limits_scaled[1])
    if (is.finite(y_limits_scaled[2])) fdr_line_scaled <- min(fdr_line_scaled, y_limits_scaled[2])
  }

  # ---- colors ----
  col_other   <- "#B8B8B8"
  col_sigmid  <- "#6E6E6E"   # darker gray for significant but not |log2FC| strong
  col_neg     <- "#1f4e79"
  col_pos     <- "#b80d2d"
  col_line    <- "#6E6E6E"

  # ---- plot base ----
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x_plot, y = y_scaled)) +
    ggplot2::geom_point(
      ggplot2::aes(color = side),
      shape = 16, alpha = point_alpha, size = point_size, show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = c(other = col_other, sig_mid = col_sigmid, neg = col_neg, pos = col_pos),
                                guide = "none") +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.6, color = col_line) +
    { if (l2fc_thresh > 0)
        ggplot2::geom_vline(xintercept = c(-l2fc_thresh, l2fc_thresh),
                            linetype = "dotted", linewidth = 0.7, color = col_line)
      else NULL } +
    # dashed horizontal line at Wald closest to FDR cutoff (on scaled y)
    { if (is.finite(fdr_line_scaled))
        ggplot2::geom_hline(yintercept = fdr_line_scaled, linetype = "dashed",
                            linewidth = 0.6, color = col_line)
      else NULL } +
    ggplot2::scale_x_continuous(limits = if (!is.null(x_lim) && length(x_lim)==2) c(min(x_lim), max(x_lim)) else NULL) +
    ggplot2::scale_y_continuous(
      limits = y_limits_scaled,
      breaks = tick_pos,
      labels = tick_lab
    ) +
    ggplot2::labs(
      x = "Signed max log2FC (vs. start)",
      y = expression(log[10]("Wald"))
    ) +
    ggplot2::theme_classic(base_size = 20) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "grey85", linewidth = 0.4),
      panel.grid.minor.y = ggplot2::element_blank(),
      axis.title.x       = ggplot2::element_text(size = axis_title_size,
                                                 margin = ggplot2::margin(t = 12)),
      axis.title.y       = ggplot2::element_text(size = axis_title_size, face = "bold",
                                                 margin = ggplot2::margin(r = 12)),
      axis.text          = ggplot2::element_text(size = axis_text_size),
      plot.title         = ggplot2::element_blank(),
      plot.margin        = ggplot2::margin(14, 18, 14, 16)
    )

  # ---- ggrepel labels: kept off each other and off every point (unlabelled points
  # are passed with empty labels so the labels avoid them too) ----
  if (label_method == "repel" && nrow(label_df) > 0) {
    if (!requireNamespace("ggrepel", quietly = TRUE)) stop("label_method = 'repel' needs the ggrepel package.")
    rep_df <- df[is.finite(df$x_plot) & is.finite(df$y_scaled), , drop = FALSE]
    sgn <- ifelse(rep_df$x_plot < 0, -1, 1)
    x_rng <- if (!is.null(x_lim) && length(x_lim) == 2 && all(is.finite(x_lim))) sort(x_lim) else c(NA, NA)
    # Nudge labels outward, but not past the x limits (a label starting outside the
    # scale range would be dropped).
    x_to <- rep_df$x_plot + ifelse(rep_df$label != "", sgn * repel_nudge_x, 0)
    if (all(is.finite(x_rng))) x_to <- pmin(pmax(x_to, x_rng[1]), x_rng[2])
    p <- p + ggrepel::geom_label_repel(
      data = rep_df,
      ggplot2::aes(x = x_plot, y = y_scaled, label = label),
      inherit.aes = FALSE,
      size = label_size,
      label.size = 0.15,
      label.padding = ggplot2::unit(0.10, "lines"),
      fill = scales::alpha("white", 0.75), color = "black",
      box.padding = ggplot2::unit(repel_box_padding, "lines"),
      point.padding = ggplot2::unit(repel_point_padding, "lines"),
      force = repel_force,
      nudge_x = x_to - rep_df$x_plot,
      min.segment.length = 0,
      segment.colour = "#555555", segment.alpha = label_segment_alpha, segment.size = label_segment_size,
      max.overlaps = Inf, max.time = 5, max.iter = 1e5,
      seed = repel_seed
    )
  }

  # ---- dependency-free label "repel": above + spaced vertically ----
  if (label_method == "stack" && nrow(label_df) > 0) {
    label_df$y_scaled <- ps(log10(label_df$y_plot))

    y_rng   <- diff(range(y_limits_scaled[is.finite(y_limits_scaled)]))
    if (!is.finite(y_rng) || y_rng == 0) y_rng <- diff(range(df$y_scaled, na.rm = TRUE))
    y_off   <- label_y_offset_frac * y_rng
    y_gap   <- label_min_vgap_frac * y_rng

    label_df$y_lab <- label_df$y_scaled + y_off
    ord <- order(label_df$y_lab, decreasing = FALSE)
    label_df <- label_df[ord, , drop = FALSE]
    if (nrow(label_df) > 1) {
      for (i in 2:nrow(label_df)) {
        min_allowed <- label_df$y_lab[i-1] + y_gap
        if (label_df$y_lab[i] < min_allowed) label_df$y_lab[i] <- min_allowed
      }
    }
    y_top <- if (is.finite(y_limits_scaled[2])) y_limits_scaled[2] else max(df$y_scaled, na.rm = TRUE)
    if (nrow(label_df) > 0 && is.finite(y_top)) {
      overflow <- max(0, label_df$y_lab[nrow(label_df)] - (y_top - y_gap*0.5))
      if (overflow > 0) label_df$y_lab <- label_df$y_lab - overflow
    }

    # sign-based horizontal nudge and clamp
    sgn <- ifelse(is.finite(label_df$x_plot) & label_df$x_plot < 0, -1, 1)
    label_df$x_lab <- label_df$x_plot + sgn * label_x_nudge
    if (!is.null(x_lim) && length(x_lim) == 2 && all(is.finite(x_lim))) {
      label_df$x_lab <- pmin(pmax(label_df$x_lab, min(x_lim)), max(x_lim))
    }

    if (isTRUE(stack_resolve_overlaps) && requireNamespace("ggrepel", quietly = TRUE)) {
      # Start each label at its stacked position (nudge from its point) and let ggrepel
      # push apart only the labels that overlap; springs pull labels back to that start.
      p <- p + ggrepel::geom_label_repel(
        data = label_df,
        ggplot2::aes(x = x_plot, y = y_scaled, label = label),
        inherit.aes = FALSE,
        nudge_x = label_df$x_lab - label_df$x_plot,
        nudge_y = label_df$y_lab - label_df$y_scaled,
        size = label_size, label.size = 0.15,
        label.padding = ggplot2::unit(0.10, "lines"),
        fill = scales::alpha("white", 0.75), color = "black",
        box.padding = ggplot2::unit(stack_box_padding, "lines"),
        point.padding = ggplot2::unit(0, "lines"),
        force = stack_force, force_pull = 1,
        min.segment.length = 0,
        segment.colour = "#555555", segment.alpha = label_segment_alpha, segment.size = label_segment_size,
        max.overlaps = Inf, seed = repel_seed
      )
    } else {
    # leader lines + labels
    p <- p + ggplot2::geom_segment(
      data = label_df,
      ggplot2::aes(x = x_plot, xend = x_lab, y = y_scaled, yend = y_lab),
      inherit.aes = FALSE,
      linewidth = label_segment_size, alpha = label_segment_alpha, color = "#555555"
    ) +
    ggplot2::geom_label(
      data = label_df,
      ggplot2::aes(x = x_lab, y = y_lab, label = label),
      size = label_size, label.size = 0.15,
      label.padding = ggplot2::unit(0.10, "lines"),
      fill = scales::alpha("white", 0.75), color = "black"
    )
    }
  }

  if (!is.null(save_prefix)) {
    ggplot2::ggsave(paste0(save_prefix, ".png"), p, width = 7.5, height = 6.5, dpi = 500)
    ggplot2::ggsave(paste0(save_prefix, ".pdf"), p, width = 7.5, height = 6.5)
  }
  p
}

# cp.pl.association_volcano: plot_association_volcano() on an association table with
# signedMaxLog2FC_start (rownames from row_id). Returns the labelled genes.
chronopack_plot_association_volcano_csv <- function(data_path, output, width = 7.5, height = 6.5, dpi = 500, ...) {
  asso_df <- chronopack_read_plot_table(data_path)
  asso_df$row_id <- NULL
  p <- plot_association_volcano(asso_df, ...)
  chronopack_save_plot(p, output, width = width, height = height, dpi = dpi)
  label_layers <- Filter(function(l) inherits(l$geom, c("GeomLabel", "GeomLabelRepel")), p$layers)
  cols <- c("signedMaxLog2FC_start", "waldStat_1", "qval_1")
  if (!length(label_layers)) {
    return(data.frame(gene = character(), signedMaxLog2FC_start = numeric(), waldStat_1 = numeric(),
                      qval_1 = numeric(), side = character()))
  }
  d <- label_layers[[1]]$data
  d <- d[d$label != "", , drop = FALSE]
  d <- d[order(d$waldStat_1, decreasing = TRUE), , drop = FALSE]
  data.frame(gene = d$label, d[, cols, drop = FALSE], side = as.character(d$side), row.names = NULL)
}

# Two-context response-specificity volcano (VolcanoPlot.ipynb plot_association_volcano):
# |log2FC(condA)| - |log2FC(condB)| against the larger (or difference of the) Wald statistics,
# coloured by condition-specific induction / repression. Input: cp.tl.combine_context_associations().
plot_response_specificity_volcano <- function(
  asso_df,

  # ---------- condition display names ----------
  condition_key = c(
    condA = "Condition A",
    condB = "Condition B"
  ),

  # ---------- thresholds ----------
  fdr_thresh        = 0.01,
  require_fdr       = TRUE,
  dyn_thresh_condA  = 0.75,
  flat_thresh_condA = 0.25,
  dyn_thresh_condB  = 0.75,
  flat_thresh_condB = 0.25,

  # ---------- plotting ----------
  wald_mode         = c("max", "difference"),
  point_alpha       = 0.95,
  point_size        = 2.2,
  axis_title_size   = 26,
  axis_text_size    = 18,
  legend_text_size  = 16,
  axis_title_pad_x  = 12,
  axis_title_pad_y  = 12,
  y_log_floor       = 1e-6,
  x_lim             = NULL,
  y_lim             = NULL,
  y_break           = NULL,
  y_break_quantile  = 0.95,
  y_lower_scale     = 0.35,

  # ---------- labels ----------
  top_labels        = 12,
  label_size        = 4.6,
  label_only_pass   = FALSE,
  label_rank_by     = c("wald", "x"),
  label_top_wald    = NULL,
  label_top_condA   = NULL,
  label_top_condB   = NULL,

  # ---------- ggrepel controls ----------
  repel_box_padding        = 0.25,
  repel_point_padding      = 0.20,
  repel_min_segment_length = 0,
  repel_force              = 1.2,
  repel_max_overlaps       = Inf,
  repel_seed               = 1,

  # ---------- optional legend-label overrides ----------
  lab_condA_ind = NULL,
  lab_condA_rep = NULL,
  lab_condB_ind = NULL,
  lab_condB_rep = NULL,

  # ---------- colors ----------
  col_condA_ind = "#E35335",
  col_condA_rep = "#FFE5B4",
  col_condB_ind = "#008000",
  col_condB_rep = "#90EE90",

  col_sigmid = "#6E6E6E",
  col_other  = "#B8B8B8",
  col_line   = "#6E6E6E",

  # ---------- legend ----------
  legend_pos      = c(0.03, 0.05),
  legend_just     = c(0, 0),
  legend_bg_alpha = 0.75,
  legend_bg_col   = "white",
  legend_key_size = 0.8,

  save_prefix = NULL
) {
  label_rank_by <- match.arg(label_rank_by)
  wald_mode <- match.arg(wald_mode)

  required_cols <- c(
    "waldStat_1", "waldStat_2", "qval",
    "log2FC_t1", "log2FC_t2"
  )

  stopifnot(all(required_cols %in% names(asso_df)))

  if (!all(c("condA", "condB") %in% names(condition_key))) {
    stop(
      "condition_key must be a named vector containing ",
      "'condA' and 'condB'."
    )
  }

  condA_name <- as.character(condition_key[["condA"]])
  condB_name <- as.character(condition_key[["condB"]])

  # Construct legend labels unless explicitly overridden
  if (is.null(lab_condA_ind)) {
    lab_condA_ind <- paste0(condA_name, "-specific induction")
  }

  if (is.null(lab_condA_rep)) {
    lab_condA_rep <- paste0(condA_name, "-specific repression")
  }

  if (is.null(lab_condB_ind)) {
    lab_condB_ind <- paste0(condB_name, "-specific induction")
  }

  if (is.null(lab_condB_rep)) {
    lab_condB_rep <- paste0(condB_name, "-specific repression")
  }

  x_axis_label <- sprintf(
    "Response specificity (|log2FC(%s)| - |log2FC(%s)|)",
    condA_name,
    condB_name
  )

  y_axis_label <- switch(
    wald_mode,
    max = expression(log[10]("Maximum Wald")),
    difference = expression(log[10]("Wald difference"))
  )

  df <- asso_df

  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop(
      "Package 'ggrepel' is required. ",
      "Install via install.packages('ggrepel')."
    )
  }

  # ---------- condition mapping ----------
  # t1 = condA; t2 = condB
  condA_eff <- as.numeric(df$log2FC_t1)
  condB_eff <- as.numeric(df$log2FC_t2)

  # Right side means condA is more dynamic
  x_raw <- abs(condA_eff) - abs(condB_eff)

  # ---------- Wald metric ----------
  wald_1 <- as.numeric(df$waldStat_1)
  wald_2 <- as.numeric(df$waldStat_2)

  if (wald_mode == "max") {
    y_raw <- pmax(wald_1, wald_2, na.rm = TRUE)
    y_raw[!is.finite(y_raw)] <- NA_real_
  } else {
    # Higher Wald minus lower Wald
    y_raw <- ifelse(
      is.finite(wald_1) & is.finite(wald_2),
      abs(wald_1 - wald_2),
      NA_real_
    )
  }

  q_raw <- as.numeric(df$qval)

  pass_fdr <- is.finite(q_raw) & q_raw < fdr_thresh

  pass_gate <- if (isTRUE(require_fdr)) {
    pass_fdr
  } else {
    rep(TRUE, nrow(df))
  }

  # ---------- clamp y ----------
  y_plot <- ifelse(
    is.finite(y_raw),
    pmax(y_raw, y_log_floor),
    y_log_floor
  )

  if (!is.null(y_lim) && length(y_lim) == 2) {
    if (is.finite(y_lim[1])) {
      y_plot <- pmax(
        y_plot,
        max(y_log_floor, y_lim[1])
      )
    }

    if (is.finite(y_lim[2])) {
      y_plot <- pmin(y_plot, y_lim[2])
    }
  }

  # ---------- clamp x ----------
  x_plot <- x_raw

  if (
    !is.null(x_lim) &&
    length(x_lim) == 2 &&
    all(is.finite(x_lim))
  ) {
    xmin <- min(x_lim)
    xmax <- max(x_lim)
    msk <- is.finite(x_plot)

    x_plot[msk] <- pmin(
      pmax(x_plot[msk], xmin),
      xmax
    )
  }

  df$x_plot <- x_plot
  df$y_plot <- y_plot

  # ---------- specificity and direction ----------
  a_dyn <- is.finite(condA_eff) &
    abs(condA_eff) >= dyn_thresh_condA

  a_flat <- is.finite(condA_eff) &
    abs(condA_eff) <= flat_thresh_condA

  b_dyn <- is.finite(condB_eff) &
    abs(condB_eff) >= dyn_thresh_condB

  b_flat <- is.finite(condB_eff) &
    abs(condB_eff) <= flat_thresh_condB

  condA_spec <- pass_gate & a_dyn & b_flat
  condB_spec <- pass_gate & b_dyn & a_flat

  a_pos <- is.finite(condA_eff) & condA_eff > 0
  a_neg <- is.finite(condA_eff) & condA_eff < 0
  b_pos <- is.finite(condB_eff) & condB_eff > 0
  b_neg <- is.finite(condB_eff) & condB_eff < 0

  side <- rep("other", nrow(df))
  side[pass_fdr] <- "sig_mid"

  side[condA_spec & a_pos] <- "condA_ind"
  side[condA_spec & a_neg] <- "condA_rep"
  side[condB_spec & b_pos] <- "condB_ind"
  side[condB_spec & b_neg] <- "condB_rep"

  df$side <- factor(
    side,
    levels = c(
      "other",
      "sig_mid",
      "condA_ind",
      "condA_rep",
      "condB_ind",
      "condB_rep"
    )
  )

  # ---------- label selection ----------
  df$label <- ""

  rn <- rownames(df)

  gene_names <- if (!is.null(rn)) {
    rn
  } else if ("gene" %in% names(df)) {
    as.character(df$gene)
  } else {
    as.character(seq_len(nrow(df)))
  }

  pass_side <- if (label_only_pass) {
    df$side %in% c(
      "condA_ind",
      "condA_rep",
      "condB_ind",
      "condB_rep"
    )
  } else {
    rep(TRUE, nrow(df))
  }

  labs_idx <- integer(0)

  if (!is.null(label_top_wald) && as.integer(label_top_wald) > 0) {
    cand <- which(is.finite(y_raw) & pass_side)

    if (length(cand)) {
      ord <- cand[
        order(
          y_raw[cand],
          abs(x_raw[cand]),
          decreasing = TRUE
        )
      ]

      labs_idx <- c(
        labs_idx,
        head(ord, as.integer(label_top_wald))
      )
    }
  }

  if (!is.null(label_top_condA) && as.integer(label_top_condA) > 0) {
    cand <- which(
      df$side %in% c("condA_ind", "condA_rep") &
        pass_side &
        is.finite(condA_eff) &
        is.finite(y_raw)
    )

    if (length(cand)) {
      ord <- cand[
        order(
          abs(condA_eff[cand]),
          y_raw[cand],
          decreasing = TRUE
        )
      ]

      labs_idx <- c(
        labs_idx,
        head(ord, as.integer(label_top_condA))
      )
    }
  }

  if (!is.null(label_top_condB) && as.integer(label_top_condB) > 0) {
    cand <- which(
      df$side %in% c("condB_ind", "condB_rep") &
        pass_side &
        is.finite(condB_eff) &
        is.finite(y_raw)
    )

    if (length(cand)) {
      ord <- cand[
        order(
          abs(condB_eff[cand]),
          y_raw[cand],
          decreasing = TRUE
        )
      ]

      labs_idx <- c(
        labs_idx,
        head(ord, as.integer(label_top_condB))
      )
    }
  }

  if (length(labs_idx) == 0 && top_labels > 0) {
    rank_val <- if (label_rank_by == "wald") {
      y_raw
    } else {
      x_raw
    }

    cand <- which(
      is.finite(rank_val) & pass_side
    )

    if (length(cand)) {
      ord <- cand[
        order(
          rank_val[cand],
          y_raw[cand],
          decreasing = TRUE
        )
      ]

      labs_idx <- head(
        ord,
        as.integer(top_labels)
      )
    }
  }

  labs_idx <- unique(labs_idx)
  df$label[labs_idx] <- gene_names[labs_idx]

  label_df <- df[
    df$label != "" &
      is.finite(df$x_plot) &
      is.finite(df$y_plot),
    ,
    drop = FALSE
  ]

  # ---------- piecewise y scaling ----------
  logy <- log10(df$y_plot)

  if (is.null(y_break) || !is.finite(y_break)) {
    non_sig <- df$y_plot[
      !pass_fdr & is.finite(df$y_plot)
    ]

    y_break <- if (length(non_sig)) {
      as.numeric(
        stats::quantile(
          non_sig,
          probs = y_break_quantile,
          na.rm = TRUE
        )
      )
    } else {
      as.numeric(
        stats::quantile(
          df$y_plot,
          probs = 0.60,
          na.rm = TRUE
        )
      )
    }
  }

  y_break <- max(y_break, y_log_floor)
  t_break <- log10(y_break)

  y_lower_scale <- max(
    1e-6,
    min(1, y_lower_scale)
  )

  ps <- function(t) {
    ifelse(
      t <= t_break,
      y_lower_scale * t,
      y_lower_scale * t_break + (t - t_break)
    )
  }

  df$y_scaled <- ps(log10(df$y_plot))

  t_lo <- if (!is.null(y_lim) && is.finite(y_lim[1])) {
    log10(max(y_log_floor, y_lim[1]))
  } else {
    min(logy[is.finite(logy)])
  }

  t_hi <- if (!is.null(y_lim) && is.finite(y_lim[2])) {
    log10(y_lim[2])
  } else {
    max(logy[is.finite(logy)])
  }

  y_limits_scaled <- ps(c(t_lo, t_hi))

  tick_t <- seq(
    floor(t_lo),
    ceiling(t_hi),
    by = 1
  )

  tick_pos <- ps(tick_t)
  tick_lab <- as.character(tick_t)

  label_df$y_scaled <- ps(
    log10(label_df$y_plot)
  )

  # ---------- colors ----------
  col_map <- c(
    other     = col_other,
    sig_mid   = col_sigmid,
    condA_ind = col_condA_ind,
    condA_rep = col_condA_rep,
    condB_ind = col_condB_ind,
    condB_rep = col_condB_rep
  )

  # ---------- plot ----------
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = x_plot,
      y = y_scaled
    )
  ) +
    ggplot2::geom_point(
      ggplot2::aes(colour = side),
      shape = 16,
      alpha = point_alpha,
      size = point_size,
      show.legend = TRUE
    ) +
    ggplot2::scale_colour_manual(
      values = col_map,
      breaks = c(
        "condA_ind",
        "condA_rep",
        "condB_ind",
        "condB_rep"
      ),
      labels = c(
        lab_condA_ind,
        lab_condA_rep,
        lab_condB_ind,
        lab_condB_rep
      ),
      drop = FALSE,
      name = NULL
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "solid",
      linewidth = 0.6,
      color = col_line
    ) +
    ggplot2::scale_x_continuous(
      limits = if (
        !is.null(x_lim) &&
        length(x_lim) == 2
      ) {
        c(min(x_lim), max(x_lim))
      } else {
        NULL
      }
    ) +
    ggplot2::scale_y_continuous(
      limits = y_limits_scaled,
      breaks = tick_pos,
      labels = tick_lab
    ) +
    ggplot2::labs(
      x = x_axis_label,
      y = y_axis_label
    ) +
    ggplot2::theme_classic(
      base_size = 20
    ) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(
        color = "grey85",
        linewidth = 0.4
      ),
      panel.grid.minor.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(
        size = axis_title_size,
        margin = ggplot2::margin(
          t = axis_title_pad_x
        )
      ),
      axis.title.y = ggplot2::element_text(
        size = axis_title_size,
        face = "bold",
        margin = ggplot2::margin(
          r = axis_title_pad_y
        )
      ),
      axis.text = ggplot2::element_text(
        size = axis_text_size
      ),
      plot.margin = ggplot2::margin(
        14, 18, 14, 16
      ),
      legend.position = legend_pos,
      legend.justification = legend_just,
      legend.background = ggplot2::element_rect(
        fill = scales::alpha(
          legend_bg_col,
          legend_bg_alpha
        ),
        color = "grey60",
        linewidth = 0.3
      ),
      legend.key = ggplot2::element_rect(
        fill = NA,
        color = NA
      ),
      legend.key.size = ggplot2::unit(
        legend_key_size,
        "lines"
      ),
      legend.text = ggplot2::element_text(
        size = legend_text_size
      )
    )

  # ---------- ggrepel labels ----------
  if (nrow(label_df) > 0) {
    p <- p +
      ggrepel::geom_label_repel(
        data = label_df,
        ggplot2::aes(
          x = x_plot,
          y = y_scaled,
          label = label
        ),
        inherit.aes = FALSE,
        size = label_size,
        label.size = 0.15,
        label.padding = ggplot2::unit(
          0.10,
          "lines"
        ),
        fill = scales::alpha(
          "white",
          0.75
        ),
        color = "black",
        box.padding = ggplot2::unit(
          repel_box_padding,
          "lines"
        ),
        point.padding = ggplot2::unit(
          repel_point_padding,
          "lines"
        ),
        min.segment.length = repel_min_segment_length,
        force = repel_force,
        max.overlaps = repel_max_overlaps,
        seed = repel_seed
      )
  }

  if (!is.null(save_prefix)) {
    ggplot2::ggsave(
      paste0(save_prefix, ".png"),
      p,
      width = 10,
      height = 8.75,
      dpi = 300
    )

    ggplot2::ggsave(
      paste0(save_prefix, ".pdf"),
      p,
      width = 10,
      height = 8.75,
      dpi = 300
    )
  }

  p
}


# cp.pl.response_specificity_volcano: plot_response_specificity_volcano() on a combined table
# (rownames from row_id). Returns every gene's side and whether it is labelled.
chronopack_plot_response_specificity_volcano_csv <- function(data_path, output, width = 10, height = 8.75,
                                                              dpi = 300, ...) {
  asso_df <- chronopack_read_plot_table(data_path)
  asso_df$row_id <- NULL
  p <- plot_response_specificity_volcano(asso_df, ...)
  chronopack_save_plot(p, output, width = width, height = height, dpi = dpi)
  d <- p$data
  data.frame(gene = rownames(d), x = d$x_plot, y = d$y_plot, side = as.character(d$side),
             labelled = d$label != "", row.names = NULL)
}
