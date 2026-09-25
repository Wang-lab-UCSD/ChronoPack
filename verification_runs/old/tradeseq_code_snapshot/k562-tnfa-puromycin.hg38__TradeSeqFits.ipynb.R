### CELL 0
library(tradeSeq)
library(RColorBrewer)
library(SingleCellExperiment)
library(slingshot)
library(Matrix)
library(SingleCellExperiment)
library(tradeSeq)
library(BiocParallel)
library(ggplot2)
library(dplyr)

RNGversion("3.5.0") # For reproducibility
palette(brewer.pal(8, "Dark2"))

setwd('../intmd/sc-K562-1e6-to-tradeSeq/')

### CELL 1
# ============== single-gene plot (flexible point color, ticks, time scaling) ==============
plot_tradeSeq_gene_pub_v2 <- function(
  gene, gamList, sce, counts_sub, PT, sizeFactors, cellWeights = NULL,
  q_by_gene = NULL,
  lineages   = NULL,
  nPoints    = 1000,
  colours    = NULL,
  lineage_labels = NULL,
  base_size = 18,
  title_size = 10,
  legend_size = 1.25,
  # NEW: point-color behavior
  color_points_by_lineage = FALSE,
  points_color = "#BBBBBB",
  points_alpha = 0.22,
  points_size  = 1.8,   # NEW: allow setting raw point size
  # NEW: spacing + axes
  inner_pad_pt = 6,
  time_max = NULL,          # e.g., 120 (minutes). If NULL, show pseudotime.
  nx_ticks = 5,
  ny_ticks = 5,
  ylim = NULL
) {
  stopifnot(gene %in% rownames(counts_sub))
  gm <- gamList[[gene]]

  dm <- SummarizedExperiment::colData(sce)$tradeSeq$dm
  t_cols <- grep("^t[0-9]+$", colnames(dm), value = TRUE)
  avail_lineages <- as.integer(sub("^t", "", t_cols))
  if (!length(avail_lineages)) stop("Could not find t1/t2/... columns in tradeSeq design (dm).")

  if (is.null(lineages)) lineages <- avail_lineages
  lineages <- intersect(lineages, avail_lineages)
  if (!length(lineages)) stop("Requested 'lineages' not present in dm.")

  ln_ids   <- lineages
  ln_codes <- paste0("L", ln_ids)
  if (is.null(lineage_labels)) {
    ln_labels <- ln_codes
  } else {
    stopifnot(length(lineage_labels) == length(lineages))
    ln_labels <- as.character(lineage_labels)
  }
  lab_map <- setNames(ln_labels, ln_ids)

  mk_lineage_vec <- function(ell, n) {
    lab <- unname(lab_map[as.character(ell)])
    factor(rep(lab, n), levels = ln_labels)
  }

  # reference library size
  sf_vec <- as.numeric(sizeFactors)
  s_ref  <- exp(mean(log(sf_vec), na.rm = TRUE))
  off_ref <- log(s_ref)

  # fitted curves
  fits <- lapply(ln_ids, function(ell) {
    df <- tradeSeq:::.getPredictRangeDf(dm, lineageId = ell, nPoints = nPoints)
    pr <- predict(gm, newdata = df, type = "link", se.fit = TRUE,
                  offset = rep(off_ref, nrow(df)))
    data.frame(
      lineage = mk_lineage_vec(ell, nrow(df)),
      t  = df[[paste0("t", ell)]],
      mu = log1p(exp(pr$fit)),
      lo = log1p(exp(pr$fit - 1.96 * pr$se.fit)),
      hi = log1p(exp(pr$fit + 1.96 * pr$se.fit))
    )
  })
  fit_all <- do.call(rbind, fits)

  # observed points normalized to s_ref
  y_cnt <- as.numeric(counts_sub[gene, ])
  y_cnt_norm <- y_cnt * (s_ref / sf_vec)
  get_pt_vec <- function(ell) if (is.matrix(PT)) as.numeric(PT[, ell]) else as.numeric(PT)

  # default cellWeights if missing
  if (is.null(cellWeights)) {
    cellWeights <- matrix(1, nrow = length(y_cnt_norm), ncol = max(avail_lineages))
    colnames(cellWeights) <- paste0("L", seq_len(ncol(cellWeights)))
  }

  # assign each cell to a lineage (argmax over selected lineages)
  cw_sub <- as.matrix(cellWeights[, ln_ids, drop = FALSE])
  cw_sub[!is.finite(cw_sub)] <- -Inf
  idx <- max.col(cw_sub, ties.method = "first")
  bad <- rowSums(cw_sub, na.rm = TRUE) <= 0
  idx[bad] <- NA_integer_
  lineage_obs <- factor(ifelse(is.na(idx), NA,
                               lab_map[as.character(ln_ids[idx])]),
                        levels = ln_labels)

  obs_all <- data.frame(
    pt = get_pt_vec(ln_ids[1]),
    y  = log1p(y_cnt_norm),
    lineage = lineage_obs
  )

  # weighted observed means per lineage
  wmean <- function(x, w) { s <- sum(w, na.rm = TRUE); if (is.finite(s) && s > 0) sum(x*w, na.rm=TRUE)/s else NA_real_ }
  obs_u_all <- lapply(ln_ids, function(ell) {
    pt_vec <- get_pt_vec(ell); w_col <- as.numeric(cellWeights[, ell])
    u_pt <- sort(unique(pt_vec))
    m_w <- vapply(u_pt, function(v) wmean(y_cnt_norm[pt_vec==v], w_col[pt_vec==v]), numeric(1))
    data.frame(
      lineage = mk_lineage_vec(ell, length(u_pt)),
      pt_mid  = u_pt,
      mean_norm_log1p = log1p(as.numeric(m_w))
    )
  })
  obs_u_all <- do.call(rbind, obs_u_all)

  # subtitle with FDR
  sub_txt <- NULL
  if (!is.null(q_by_gene) && gene %in% names(q_by_gene)) {
    qv <- unname(q_by_gene[gene]); if (is.finite(qv)) sub_txt <- paste0("FDR = ", formatC(qv, format="fg", digits=3))
  }

  # colors
  if (is.null(colours)) {
    base_cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a",
                   "#66a61e", "#e6ab02", "#a6761d", "#666666")
    colours <- rep(base_cols, length.out = length(ln_labels))
  } else {
    colours <- rep(colours, length.out = length(ln_labels))
  }
  names(colours) <- ln_labels
  fill_cols <- grDevices::adjustcolor(colours, alpha.f = 0.30)

  # legend size scalars
  lt  <- 12 * legend_size
  lx  <- 11 * legend_size
  lks <- grid::unit(12 * legend_size, "pt")
  lbm <- margin(l = 16 * legend_size, r = 16 * legend_size,
                t =  6 * legend_size, b =  6 * legend_size)

  # x-scale: map to minutes if time_max provided
  pt_min <- min(c(obs_all$pt, fit_all$t), na.rm = TRUE)
  pt_max <- max(c(obs_all$pt, fit_all$t), na.rm = TRUE)
  to_minutes <- function(x) (x - pt_min) / max(pt_max - pt_min, .Machine$double.eps) * time_max

  x_scale <- if (!is.null(time_max)) {
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = nx_ticks),
      labels = function(x) round(to_minutes(x)),
      name   = "Time (min)"
    )
  } else {
    scale_x_continuous(
      breaks = scales::pretty_breaks(n = nx_ticks),
      name   = "Time"
    )
  }

  y_scale <- scale_y_continuous(breaks = scales::pretty_breaks(n = ny_ticks), name = "")

  # NEW: auto ylim if NULL (captures both curves and lo/hi ribbons)
  if (is.null(ylim)) {
    ylo <- suppressWarnings(min(fit_all$lo, na.rm = TRUE))
    yhi <- suppressWarnings(max(fit_all$hi, na.rm = TRUE))
    if (is.finite(ylo) && is.finite(yhi)) {
      pad <- 0.04 * (yhi - ylo)
      if (!is.finite(pad) || pad == 0) pad <- 0.2
      ylim <- c(max(0, ylo - pad), (yhi + pad)*1.25)
    }
  }

  # ------------------ build plot ------------------
  p <- ggplot()

  if (isTRUE(color_points_by_lineage)) {
    p <- p +
      geom_point(data = subset(obs_all, !is.na(lineage)),
                 aes(pt, y, colour = lineage),
                 size = points_size, alpha = points_alpha) +
      geom_point(data = subset(obs_all, is.na(lineage)),
                 aes(pt, y),
                 size = points_size, alpha = 0.20, colour = "#BBBBBB")
  } else {
    p <- p +
      geom_point(data = obs_all, aes(pt, y),
                 size = points_size, alpha = points_alpha, colour = points_color)
  }

  p <- p +
    geom_ribbon(data = fit_all, aes(t, ymin = lo, ymax = hi, fill = lineage),
                alpha = 0.35, colour = NA) +
    geom_line(data = fit_all, aes(t, mu, colour = lineage),
              linewidth = 1.7) +
    geom_point(data = obs_u_all, aes(pt_mid, mean_norm_log1p, colour = lineage),
               shape = 21, fill = "white", size = 3.8, stroke = 1.0) +
    scale_colour_manual(values = colours, name = "Lineage") +
    scale_fill_manual(values = fill_cols,  name = "Lineage") +
    x_scale + y_scale +
    labs(title = gene, subtitle = sub_txt) +
    coord_cartesian(expand = FALSE, ylim = ylim) +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "grey88"),
      panel.grid.major.y = element_line(colour = "grey88"),
      axis.text  = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(size = title_size, face = "bold", hjust = 0, colour = "black"),
      plot.subtitle = element_text(size = 10, colour = "black"),
      legend.position = "right",
      legend.title = element_text(size = lt),
      legend.text  = element_text(size = lx),
      legend.key.size = lks,
      legend.box.margin = lbm,
      plot.margin = margin(inner_pad_pt, inner_pad_pt, inner_pad_pt, inner_pad_pt)
    )

  p
}


library(ggplot2)

plot_tradeSeq_gene_pub <- function(
  gene, gamList, sce, counts_sub, PT, sizeFactors, cellWeights = NULL,
  q_by_gene = NULL,
  lineages   = NULL,
  omit_points = FALSE,
  nPoints    = 1000,
  colours    = NULL,
  lineage_labels = NULL,

  time_max = NULL,
  time_unit = "min",

  axis_text_size  = 14,
  axis_title_size = 16,

  # NEW: if NULL and time_max is set, defaults to seq(0, time_max, by = x_tick_step)
  x_tick_labels   = NULL,
  x_tick_step     = 20,

  # NEW: y tick density (like v2)
  ny_ticks = 5,

  show_title  = TRUE,
  show_legend = TRUE,

  # NEW: raw-point controls (match v2)
  color_points_by_lineage = FALSE,
  points_color = "#274C77",
  points_alpha = 0.35,
  points_size  = 1.8,

  # NEW: ylim controls (match v2 style; NULL => auto)
  ylim = NULL,

  save_pdf    = FALSE,
  pdf_filename = NULL,
  pdf_width   = 4,
  pdf_height  = 4
) {
  stopifnot(gene %in% rownames(counts_sub))
  gm <- gamList[[gene]]

  dm <- SummarizedExperiment::colData(sce)$tradeSeq$dm
  t_cols <- grep("^t[0-9]+$", colnames(dm), value = TRUE)
  avail_lineages <- as.integer(sub("^t", "", t_cols))
  if (length(avail_lineages) == 0L)
    stop("Could not find t1/t2/... columns in tradeSeq design (dm).")

  ## ---- choose which lineages to plot -----------------------------------
  if (is.null(lineages)) {
    lineages <- avail_lineages
  } else {
    lineages <- as.integer(lineages)
    lineages <- intersect(lineages, avail_lineages)
  }
  if (!length(lineages))
    stop("Requested 'lineages' not present in dm.")

  # legend labels
  lineage_codes <- paste0("L", lineages)
  if (is.null(lineage_labels)) {
    lineage_labels <- lineage_codes
  } else {
    if (length(lineage_labels) != length(lineages))
      stop("'lineage_labels' must have same length as 'lineages'.")
    lineage_labels <- as.character(lineage_labels)
  }
  lab_map <- setNames(lineage_labels, lineage_codes)

  # library size / offset
  sf_vec <- as.numeric(sizeFactors)
  s_ref  <- exp(mean(log(sf_vec), na.rm = TRUE))
  off_ref <- log(s_ref)

  # fitted curves
  fits <- lapply(lineages, function(ell) {
    df <- tradeSeq:::.getPredictRangeDf(dm, lineageId = ell, nPoints = nPoints)
    pr <- predict(gm, newdata = df, type = "link", se.fit = TRUE,
                  offset = rep(off_ref, nrow(df)))
    lab <- unname(lab_map[paste0("L", ell)])
    data.frame(
      lineage = factor(rep(lab, nrow(df)), levels = lineage_labels),
      t  = df[[paste0("t", ell)]],
      mu = log1p(exp(pr$fit)),
      lo = log1p(exp(pr$fit - 1.96 * pr$se.fit)),
      hi = log1p(exp(pr$fit + 1.96 * pr$se.fit))
    )
  })
  fit_all <- do.call(rbind, fits)

  # observed points normalized to s_ref
  y_cnt <- as.numeric(counts_sub[gene, ])
  y_cnt_norm <- y_cnt * (s_ref / sf_vec)
  get_pt_vec <- function(ell) if (is.matrix(PT)) as.numeric(PT[, ell]) else as.numeric(PT)

  # default cellWeights if missing
  if (is.null(cellWeights)) {
    cellWeights <- matrix(1, nrow = length(y_cnt_norm), ncol = max(avail_lineages))
    colnames(cellWeights) <- paste0("L", seq_len(ncol(cellWeights)))
  }

  # scatter points use pseudotime of the first selected lineage
  obs_all <- data.frame(
    pt = get_pt_vec(lineages[1]),
    y  = log1p(y_cnt_norm)
  )

  # if coloring points by lineage, assign each cell to a lineage (argmax over selected lineages)
  if (isTRUE(color_points_by_lineage) && !omit_points) {
    cw_sub <- as.matrix(cellWeights[, lineages, drop = FALSE])
    cw_sub[!is.finite(cw_sub)] <- -Inf
    idx <- max.col(cw_sub, ties.method = "first")
    bad <- rowSums(cw_sub, na.rm = TRUE) <= 0
    idx[bad] <- NA_integer_

    lab_by_id <- setNames(unname(lab_map[paste0("L", lineages)]), lineages)
    obs_all$lineage <- factor(
      ifelse(is.na(idx), NA, lab_by_id[as.character(lineages[idx])]),
      levels = lineage_labels
    )
  }

  # weighted observed means per lineage
  wmean <- function(x, w) {
    s <- sum(w, na.rm = TRUE)
    if (is.finite(s) && s > 0) sum(x*w, na.rm=TRUE)/s else NA_real_
  }
  obs_u_all <- lapply(lineages, function(ell) {
    pt_vec <- get_pt_vec(ell)
    w_col  <- as.numeric(cellWeights[, ell])
    u_pt   <- sort(unique(pt_vec))
    m_w    <- vapply(u_pt, function(v) wmean(y_cnt_norm[pt_vec==v], w_col[pt_vec==v]), numeric(1))
    lab <- unname(lab_map[paste0("L", ell)])
    data.frame(
      lineage = factor(rep(lab, length(u_pt)), levels = lineage_labels),
      pt_mid  = u_pt,
      mean_norm_log1p = log1p(as.numeric(m_w))
    )
  })
  obs_u_all <- do.call(rbind, obs_u_all)

  # subtitle with FDR
  sub_txt <- NULL
  if (!is.null(q_by_gene) && gene %in% names(q_by_gene)) {
    qv <- unname(q_by_gene[gene])
    if (is.finite(qv)) sub_txt <- paste0("FDR = ", formatC(qv, format = "fg", digits = 3))
  }

  # colors
  if (is.null(colours)) {
    base_cols <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a",
                   "#66a61e", "#e6ab02", "#a6761d", "#666666")
    colours <- rep(base_cols, length.out = length(lineage_labels))
  } else {
    colours <- rep(colours, length.out = length(lineage_labels))
  }
  names(colours) <- lineage_labels
  fill_cols <- grDevices::adjustcolor(colours, alpha.f = 0.30)

  legend_pos <- if (show_legend) "right" else "none"

  # ---- auto ylim if NULL: include points + lineage means + fitted lo/hi (matches your desired look)
    if (is.null(ylim)) {
      ylo <- suppressWarnings(min(fit_all$lo, na.rm = TRUE))
      yhi <- suppressWarnings(max(fit_all$hi, na.rm = TRUE))
    
      # (optional) also ensure the lineage mean dots fit
      ylo2 <- suppressWarnings(min(obs_u_all$mean_norm_log1p, na.rm = TRUE))
      yhi2 <- suppressWarnings(max(obs_u_all$mean_norm_log1p, na.rm = TRUE))
      ylo <- min(ylo, ylo2, na.rm = TRUE)
      yhi <- max(yhi, yhi2, na.rm = TRUE)
    
      if (is.finite(ylo) && is.finite(yhi)) {
        pad <- 0.04 * (yhi - ylo)
        if (!is.finite(pad) || pad == 0) pad <- 0.2
        ylim <- c(max(0, ylo - pad), 1.3 * (yhi + pad))  # <-- your old style
      } else {
        ylim <- c(0, 1)
      }
    }


  # y scale: show expression ticks (no limits here to avoid dropping points)
  y_scale <- scale_y_continuous(
    breaks = scales::pretty_breaks(n = ny_ticks),
    name = "Expression"
  )

  # build plot
  p <- ggplot()

  if (!omit_points) {
    if (isTRUE(color_points_by_lineage)) {
      p <- p +
        geom_point(data = subset(obs_all, !is.na(lineage)),
                   aes(pt, y, colour = lineage),
                   size = points_size, alpha = points_alpha) +
        geom_point(data = subset(obs_all, is.na(lineage)),
                   aes(pt, y),
                   size = points_size, alpha = points_alpha, colour = "#BBBBBB")
    } else {
      p <- p +
        geom_point(data = obs_all, aes(pt, y),
                   size = points_size, alpha = points_alpha, colour = points_color)
    }
  }

  p <- p +
    geom_ribbon(data = fit_all,
                aes(t, ymin = lo, ymax = hi, fill = lineage),
                alpha = 0.35, colour = NA) +
    geom_line(data = fit_all,
              aes(t, mu, colour = lineage),
              linewidth = 1.7) +
    geom_point(data = obs_u_all,
               aes(pt_mid, mean_norm_log1p, colour = lineage),
               shape = 21, fill = "white", size = 3.8, stroke = 1.0) +
    scale_colour_manual(values = colours, name = "Lineage") +
    scale_fill_manual(values = fill_cols,  name = "Lineage") +
    y_scale +
    coord_cartesian(expand = FALSE, ylim = ylim) +
    theme_minimal(base_size = axis_text_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "grey88"),
      panel.grid.major.y = element_line(colour = "grey88"),
      axis.text  = element_text(colour = "black", size = axis_text_size),
      axis.title = element_text(colour = "black", size = axis_title_size),
      plot.title = element_text(face = "bold", hjust = 0, colour = "black"),
      plot.subtitle = element_text(size = axis_text_size, colour = "black"),
      legend.position = legend_pos
    )

  # titles
  title_val    <- if (show_title) gene    else NULL
  subtitle_val <- if (show_title) sub_txt else NULL

  # x-axis labeling: convert 0..1 to 0..time_max if provided
  if (!is.null(time_max)) {
    stopifnot(is.numeric(time_max), length(time_max) == 1, is.finite(time_max), time_max > 0)

    if (is.null(x_tick_labels)) {
      x_tick_labels <- seq(0, time_max, by = x_tick_step)
    }

    br_time <- x_tick_labels
    br_pt   <- br_time / time_max
    keep    <- br_pt >= 0 & br_pt <= 1
    br_time <- br_time[keep]; br_pt <- br_pt[keep]

    x_title <- paste0("Time", if (!is.null(time_unit)) paste0(" (", time_unit, ")") else "")
    p <- p +
      scale_x_continuous(breaks = br_pt, labels = br_time, limits = c(0, 1)) +
      labs(title = title_val, subtitle = subtitle_val, x = x_title, y = "Expression")
  } else {
    p <- p +
      labs(
        title    = title_val,
        subtitle = subtitle_val,
        x = "Pseudotime",
        y = "log(normalized count + 1) @ reference library size"
      )
  }

  # optional PDF export
  if (save_pdf) {
    if (is.null(pdf_filename)) pdf_filename <- paste0(gene, "_tradeSeq.pdf")
    ggplot2::ggsave(filename = pdf_filename, plot = p,
                    width = pdf_width, height = pdf_height, units = "in")
  }

  return(p)
}

# ============== panel wrapper (size, spacing, ticks passed through) ==============
plot_tradeSeq_panel <- function(
  genes,
  gamList, sce, counts_sub, PT,
  sizeFactors, sizefactors = NULL,
  cellWeights = NULL, q_by_gene = NULL,
  lineages = NULL, nPoints = 1000, colours = NULL,
  lineage_labels = NULL,
  base_size = 18,
  title_size = 10,
  legend_size = 1.25,
  # NEW controls
  color_points_by_lineage = FALSE,
  points_color = "#BBBBBB",
  points_alpha = 0.22,
  points_size  = 1.8,   # NEW: pass-through
  inner_pad_pt = 6,
  panel_outer_pad_pt = 2,   # outer margin around the whole patchwork
  time_max = NULL,
  nx_ticks = 5,
  ny_ticks = 5,
  ylim = NULL,          # NEW: pass-through; NULL -> auto per gene
  ncol = 3,
  width = NULL, height = NULL, dpi = 300,
  save_path = NULL
) {
  if (missing(sizeFactors) && !is.null(sizefactors)) sizeFactors <- sizefactors
  if (is.function(sizeFactors))
    stop("You passed a function (likely DESeq2::sizeFactors). Provide the numeric vector used in fitGAM.")
  if (!is.numeric(sizeFactors) || length(sizeFactors) != ncol(counts_sub))
    stop("sizeFactors must be numeric, length ncol(counts_sub).")

  if (!is.null(width) && !is.null(height)) {
    local_options(list(repr.plot.width = width, repr.plot.height = height))
  }

  pls <- lapply(genes, function(g) {
    if (!g %in% rownames(counts_sub)) { warning(sprintf("Skipping '%s': not found in counts_sub", g)); return(NULL) }
    plot_tradeSeq_gene_pub_v2(
      gene = g,
      gamList = gamList, sce = sce, counts_sub = counts_sub,
      PT = PT, sizeFactors = sizeFactors, cellWeights = cellWeights,
      q_by_gene = q_by_gene, lineages = lineages, nPoints = nPoints,
      colours = colours, lineage_labels = lineage_labels,
      base_size = base_size, title_size = title_size, legend_size = legend_size,
      color_points_by_lineage = color_points_by_lineage,
      points_color = points_color, points_alpha = points_alpha,
      points_size = points_size,
      inner_pad_pt = inner_pad_pt,
      time_max = time_max, nx_ticks = nx_ticks, ny_ticks = ny_ticks,
      ylim = ylim
    )
  })
  pls <- Filter(Negate(is.null), pls)
  if (!length(pls)) stop("No valid genes to plot.")

  # panel legend sizing (same scale)
  lt  <- 12 * legend_size
  lx  <- 11 * legend_size
  lks <- grid::unit(12 * legend_size, "pt")
  lbm <- margin(l = 18 * legend_size, r = 18 * legend_size,
                t =  8 * legend_size, b =  8 * legend_size)

  panel <- wrap_plots(pls, ncol = ncol, guides = "collect") &
    theme(
      legend.position   = "right",
      legend.title      = element_text(size = lt),
      legend.text       = element_text(size = lx),
      legend.key.size   = lks,
      legend.box.margin = lbm,
      plot.margin       = margin(panel_outer_pad_pt, panel_outer_pad_pt,
                                 panel_outer_pad_pt, panel_outer_pad_pt)
    )

  if (!is.null(save_path)) {
    if (is.null(width) || is.null(height))
      stop("Provide 'width' and 'height' (inches) with save_path.")
    ggsave(save_path, panel, width = width, height = height, dpi = dpi)
  }

  panel
}

### CELL 3
# load counts (genes x cells), times, subset genes, and full libsizes
counts <- as.matrix(read.csv(gzfile("counts_genes_x_cells.csv.gz"),
                             row.names = 1, check.names = FALSE))
storage.mode(counts) <- "integer"

ct = read.csv(gzfile("cell_times.csv.gz"), stringsAsFactors = FALSE)
libs_full  = read.csv("cell_libsize_full.csv", stringsAsFactors = FALSE)

# align vectors to the count matrix cell order
stopifnot(all(colnames(counts) %in% ct$cell_id),
          all(colnames(counts) %in% libs_full$cell_id))

time_vec = setNames(ct$time, ct$cell_id)[ colnames(counts) ]
lib_full = setNames(libs_full$libsize_full, libs_full$cell_id)[ colnames(counts) ]

genes_sub = read.csv("subset_genes.csv", stringsAsFactors = FALSE)[,1]
genes_keep = intersect(genes_sub, rownames(counts))
counts_sub = counts[genes_keep, , drop = FALSE]

sf = as.numeric(lib_full)
sf = sf/mean(sf)

hist(sf)

### CELL 4
dim(counts_sub)

### CELL 5
# pseudotime (single trajectory) rescaled to [0,1]
pt01 = (time_vec - min(time_vec, na.rm = TRUE)) / diff(range(time_vec, na.rm = TRUE))
PT   = matrix(pt01, ncol = 1); colnames(PT) <- "lineage1"
CW = matrix(1, nrow = ncol(counts_sub), ncol = 1); colnames(CW) <- "lineage1" # weights: all ones (single lineage)

# ---- tradeSeq fit ----
suppressPackageStartupMessages({
  library(tradeSeq); library(SingleCellExperiment); library(BiocParallel)
})

dim(CW)

offset_mat = matrix(log(sf), nrow = nrow(counts_sub), ncol = ncol(counts_sub), byrow = TRUE)

### CELL 7
sce = fitGAM(counts = counts_sub, pseudotime  = PT, cellWeights = CW, offset = offset_mat, nknots = 4, 
             BPPARAM  = BiocParallel::MulticoreParam(workers = 4), verbose = TRUE)

### CELL 8
assoRes = associationTest(sce, global = FALSE, lineages = TRUE, inverse='Chol', nPoints=9)
assoRes = assoRes[order(assoRes[, 'waldStat_1'], decreasing = T), ]
asso_df = as.data.frame(assoRes)
if (!"meanLogFC" %in% names(asso_df) && ncol(asso_df) == 4)
  names(asso_df) = c("waldStat_1","df_1","pvalue_1","meanLogFC")
asso_df$qval_1 = p.adjust(asso_df$pvalue_1, method = "BH") # Benjamini–Hochberg FDR for lineage 1
sig_df = subset(asso_df, !is.na(qval_1) & qval_1 < 0.05) # keep significant genes (ignore NAs)
dim(sig_df)
asso_df['NFKBIA',]

### CELL 9
asso_df = asso_df[!is.na(asso_df$waldStat_1), ]
genes_keep = rownames(asso_df[asso_df$pvalue_1 < 0.05,])
length(genes_keep)

### CELL 10
write.csv(asso_df, 'asso_df.1trajectory.nknots4.csv', row.names = TRUE)
write.csv(sig_df, 'sig_df.1trajectory.nknots4.csv', row.names = TRUE)

### CELL 11
counts_sub2 = counts_sub[genes_keep, , drop = FALSE]

stopifnot(all(colnames(counts_sub2) %in% ct$cell_id),
          all(colnames(counts_sub2) %in% libs_full$cell_id))

dim(counts_sub2)

### CELL 12
sce = fitGAM(counts = counts_sub2, pseudotime  = PT, cellWeights = CW, offset = offset_mat, nknots = 8, 
             BPPARAM  = BiocParallel::MulticoreParam(workers = 4), verbose = TRUE)

saveRDS(sce, 'sce.1trajectory.nknots8.rds')

### CELL 13
gamList = fitGAM(counts = counts_sub2, pseudotime = PT, cellWeights = CW, offset = offset_mat, nknots = 8,
            sce = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 4), verbose = TRUE)

saveRDS(gamList, 'gamList.1trajectory.nknots8.rds')

### CELL 14
options(repr.plot.width = 5, repr.plot.height = 3, repr.plot.res = 300)

gene_symbol = 'NFKBIA'

plot_tradeSeq_gene_pub_v2(
  gene = gene_symbol, gamList = gamList, sce = sce,
  counts_sub = counts_sub, PT = PT, sizeFactors = sf,
  lineages = c(1, 2), cellWeights = CW,
  colours  = c("#222021", "#222021"),
  color_points_by_lineage = TRUE,
  points_size = 0.3, points_alpha = 0.15,
)

### CELL 15
options(repr.plot.width = 5.25, repr.plot.height = 3.5)

gene_symbol = 'CXCL3'

plot_tradeSeq_gene_pub_v2(
  gene = gene_symbol, gamList = gamList, sce = sce,
  counts_sub = counts_sub, PT = PT, sizeFactors = sf,
  lineages = c(1, 2), cellWeights = CW,
  colours  = c("#222021", "#222021"),
  color_points_by_lineage = TRUE,
  points_size = 0.3, points_alpha = 0.15,
)

### CELL 16
options(repr.plot.width = 5.25, repr.plot.height = 3.5)

gene_symbol = 'JUN'

plot_tradeSeq_gene_pub_v2(
  gene = gene_symbol, gamList = gamList, sce = sce,
  counts_sub = counts_sub, PT = PT, sizeFactors = sf,
  lineages = c(1, 2), cellWeights = CW,
  colours  = c("#222021", "#222021"),
  color_points_by_lineage = TRUE,
  points_size = 0.3, points_alpha = 0.15,
)

### CELL 17
options(repr.plot.width = 5.25, repr.plot.height = 3.5)

gene_symbol = 'IRF1'

plot_tradeSeq_gene_pub_v2(
  gene = gene_symbol, gamList = gamList, sce = sce,
  counts_sub = counts_sub, PT = PT, sizeFactors = sf,
  lineages = c(1, 2), cellWeights = CW,
  colours  = c("#222021", "#222021"),
  color_points_by_lineage = TRUE,
  points_size = 0.3, points_alpha = 0.15,
)

### CELL 18
options(repr.plot.width = 5.25, repr.plot.height = 3.5)

gene_symbol = 'TNFAIP3'

plot_tradeSeq_gene_pub_v2(
  gene = gene_symbol, gamList = gamList, sce = sce,
  counts_sub = counts_sub, PT = PT, sizeFactors = sf,
  lineages = c(1, 2), cellWeights = CW,
  colours  = c("#222021", "#222021"),
  color_points_by_lineage = TRUE,
  points_size = 0.3, points_alpha = 0.15,
)
