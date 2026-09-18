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

RNGversion("3.5.0") # For reproducibility
palette(brewer.pal(8, "Dark2"))

setwd('../intmd/sc-K562-1e4_to_tradeSeq//')

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
# pseudotime (single trajectory) rescaled to [0,1]
pt01 = (time_vec - min(time_vec, na.rm = TRUE)) / diff(range(time_vec, na.rm = TRUE))
PT   = matrix(pt01, ncol = 1); colnames(PT) <- "lineage1"
CW = matrix(1, nrow = ncol(counts_sub), ncol = 1); colnames(CW) <- "lineage1" # weights: all ones (single lineage)
CW_1traj = CW

# ---- tradeSeq fit ----
suppressPackageStartupMessages({
  library(tradeSeq); library(SingleCellExperiment); library(BiocParallel)
})

### CELL 5
library(dplyr)

### CELL 6
# tmp = readRDS('../../intmd/global_analysis/asso_df_ls.rds')

# genes1 = rownames(tmp$K562a %>% filter(pvalue_1 < 0.05))
# genes2 = rownames(tmp$THP1 %>% filter(pvalue_1 < 0.05))

# write.csv(union(genes1, genes2), 'ora_bg.csv')

### CELL 8
write.csv(genes1, 'ora_bg.K562a.csv')
write.csv(genes2, 'ora_bg.THP1.csv')

### CELL 9
getwd()

### CELL 13
offset_mat = matrix(log(sf), nrow = nrow(counts_sub), ncol = ncol(counts_sub), byrow = TRUE)

sce = fitGAM(counts = counts_sub, pseudotime  = PT, cellWeights = CW, offset = offset_mat, nknots = 4, 
             BPPARAM  = BiocParallel::MulticoreParam(workers = 4), verbose = TRUE)

### CELL 14
assoRes = associationTest(sce, global = FALSE, lineages = TRUE, inverse='Chol', nPoints=9)
assoRes = assoRes[order(assoRes[, 'waldStat_1'], decreasing = T), ]
asso_df = as.data.frame(assoRes)
if (!"meanLogFC" %in% names(asso_df) && ncol(asso_df) == 4)
  names(asso_df) = c("waldStat_1","df_1","pvalue_1","meanLogFC")
asso_df$qval_1 = p.adjust(asso_df$pvalue_1, method = "BH") # Benjamini–Hochberg FDR for lineage 1
sig_df = subset(asso_df, !is.na(qval_1) & qval_1 < 0.05) # keep significant genes (ignore NAs)
dim(sig_df)
asso_df['NFKBIA',]

### CELL 15
# order by q-value (then optionally by effect size)
sig_df <- sig_df[order(-sig_df$waldStat_1, -sig_df$meanLogFC), ]
nrow(sig_df)
head(sig_df, 10)

# vector of significant gene names
sig_genes <- rownames(sig_df)

### CELL 16
write.csv(asso_df, 'asso_df.1trajectory.csv', row.names = TRUE)
write.csv(sig_df, 'sig_df.1trajectory.csv', row.names = TRUE)

### CELL 18
asso_df = asso_df[!is.na(asso_df$waldStat_1), ]
genes_keep = rownames(asso_df[asso_df$pvalue_1 < 0.05,])
counts_sub2 = counts_sub[genes_keep, , drop = FALSE]

stopifnot(all(colnames(counts_sub2) %in% ct$cell_id),
          all(colnames(counts_sub2) %in% libs_full$cell_id))

dim(counts_sub2)

### CELL 19
sce = fitGAM(counts = counts_sub2, pseudotime  = PT, cellWeights = CW, offset = offset_mat, nknots = 8, 
             BPPARAM  = BiocParallel::MulticoreParam(workers = 4), verbose = TRUE)

saveRDS(sce, 'sce.1trajectory.nknots8.csv')

### CELL 20
gamList = fitGAM(counts = counts_sub2, pseudotime = PT, cellWeights = CW, offset = offset_mat, nknots = 8,
            sce = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 4), verbose = TRUE)

saveRDS(gamList, 'gamList.1trajectory.csv')

### CELL 22
# p-value column can be "pvalue_1" (multi-lineage) or "pvalue" (single-lineage)
pv_col <- grep("^pvalue", names(asso_df), value = TRUE)[1]

q_by_gene <- setNames(
  p.adjust(asso_df[[pv_col]], method = "BH"),
  rownames(asso_df)
)

### CELL 23
options(repr.plot.width = 5.25, repr.plot.height = 3.5)

gene_symbol = 'NFKBIA'

plot_tradeSeq_gene_pub_v2(
  gene = gene_symbol, gamList = gamList, sce = sce,
  counts_sub = counts_sub, PT = PT, sizeFactors = sf,
  lineages = c(1, 2), cellWeights = CW_1traj,
  colours  = c("#222021", "#222021"),
  color_points_by_lineage = TRUE,
  points_size = 0.3, points_alpha = 0.15,
)

### CELL 24
# save.image()
# load('.RData')

### CELL 27
retentions = read.csv("cluster_memberships.retentions.csv", row.names = 1, check.names = FALSE)

CW = as.matrix(retentions[, c("softmem_0", "softmem_1")])
colnames(CW) = c("L1", "L2")
CW = CW[colnames(counts_sub), , drop = FALSE]
stopifnot(all(rownames(CW) == colnames(counts_sub)))
rs = rowSums(CW)

if (any(abs(rs - 1) > 1e-8)) {
  CW <- sweep(CW, 1, rowSums(CW), "/")  # renormalize just in case
}

PT_mat = cbind(L1 = PT, L2 = PT)
rownames(PT_mat) = colnames(counts_sub)

# tradeSeq expects NA pseudotime when a cell isn't in a lineage
eps = 1e-6
CWz = CW
CWz[CWz < eps] = 0
# PT_mat[CWz == 0] <- 0

CWz_ret = CWz

offset_mat = matrix(log(sf), nrow = nrow(counts_sub), ncol = ncol(counts_sub), byrow = TRUE)

### CELL 28
sce2 = fitGAM(counts = counts_sub, pseudotime  = PT_mat, cellWeights = CWz, offset = offset_mat,
              nknots = 4, BPPARAM = MulticoreParam(workers = 4), verbose = TRUE)

### CELL 29
assoRes2 = associationTest(sce2, global = TRUE, lineages = TRUE, inverse='Chol', nPoints=9)
assoRes2 = assoRes2[order(assoRes2[, 'waldStat_1'], decreasing = T), ]

sum(is.na(assoRes2['waldStat_1']))
sum(is.na(assoRes2['waldStat_2']))
sum(is.na(assoRes2['waldStat_1']) & is.na(assoRes2['waldStat_2']))

asso_df2 = as.data.frame(assoRes2)
asso_df2$qval = p.adjust(asso_df2$pvalue, method = "BH")
asso_df2$qval_1 = p.adjust(asso_df2$pvalue_1, method = "BH")
asso_df2$qval_2 = p.adjust(asso_df2$pvalue_2, method = "BH")

### CELL 30
mask = (asso_df2$qval_1 < 0.01) | (asso_df2$qval_2 < 0.01)
sig_genes2 = rownames(asso_df2)[mask]

sig_df1 = sig_df[sig_df$qval_1 < 0.01, ]
sig_genes1 = rownames(sig_df1)
sig_genes = union(sig_genes2, sig_genes1)
sig_df2 = asso_df2[sig_genes, ]

### CELL 31
# write.csv(asso_df2, 'asso_df.2trajectories.retention.v2.csv', row.names = TRUE)
# write.csv(sig_df2, 'sig_df.2trajectories.retention.v2.csv', row.names = TRUE)

### CELL 32
asso_df2 = asso_df2[!(is.na(asso_df2$waldStat) & is.na(asso_df2$waldStat_1) & is.na(asso_df2$waldStat_2)), ]

### CELL 33
mask <- (asso_df2$pvalue   < 0.05 & !is.na(asso_df2$pvalue)) |
        (asso_df2$pvalue_1 < 0.05 & !is.na(asso_df2$pvalue_1)) |
        (asso_df2$pvalue_2 < 0.05 & !is.na(asso_df2$pvalue_2))

genes_keep2 = rownames(asso_df2[mask, ])
genes_keep2 = union(genes_keep, genes_keep2)

### CELL 34
counts_sub2 = counts_sub[genes_keep2, , drop = FALSE]

stopifnot(all(colnames(counts_sub2) %in% ct$cell_id),
          all(colnames(counts_sub2) %in% libs_full$cell_id))

dim(counts_sub2)

### CELL 35
gamList2 = fitGAM(counts = counts_sub2, pseudotime = PT_mat, cellWeights = CWz, offset = offset_mat, nknots = 8,
            sce = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 4), verbose = TRUE)

### CELL 36
counts_sub3 = counts_sub[c('HBE1', 'HBA2'), , drop = FALSE]
dim(counts_sub3)

sce3 = fitGAM(counts = counts_sub3, pseudotime  = PT_mat, cellWeights = CWz, offset = offset_mat,
              nknots = 4, BPPARAM = MulticoreParam(workers = 4), verbose = TRUE)

gamList3 = fitGAM(counts = counts_sub3, pseudotime = PT_mat, cellWeights = CWz, offset = offset_mat, nknots = 8,
            sce = FALSE, BPPARAM = BiocParallel::MulticoreParam(workers = 4), verbose = TRUE)

### CELL 37
saveRDS(gamList2, 'gamList.2trajectory.csv')

### CELL 38
signed_max_log2fc_start_many <- function(
  sce, genes, nPoints = 100,
  fit_scale = c("auto","log","response"),
  pseudocount = 0.1   # <-- your shrink knob, in mean-count units
) {
  fit_scale <- match.arg(fit_scale)
  pc <- max(0, as.numeric(pseudocount))

  pr <- tryCatch(
    as.data.frame(predictSmooth(sce, gene = genes, nPoints = nPoints, tidy = TRUE),
                  stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(pr) || nrow(pr) == 0) return(setNames(rep(NA_real_, length(genes)), genes))

  gene_col <- intersect(c("gene","Gene","feature","Feature"), names(pr))
  fit_col  <- intersect(c("yhat","mu","fit","y","estimate"), names(pr))
  time_col <- intersect(c("pseudotime","time","t","x","lambda","arclength"), names(pr))
  lin_col  <- intersect(c("lineage","Lineage","lin","curve"), names(pr))

  if (length(gene_col) == 0) stop("No gene column found in predictSmooth() output.")
  if (length(fit_col)  == 0) stop("No fitted-value column found in predictSmooth() output.")
  if (length(lin_col)  == 0) { pr$.__tmp_lineage__ <- "L1"; lin_col <- ".__tmp_lineage__" }

  gene_col <- gene_col[1]
  fit_col  <- fit_col[1]
  time_col <- if (length(time_col)) time_col[1] else NA_character_
  lin_col  <- lin_col[1]

  infer_scale <- function(colname, v) {
    if (fit_scale != "auto") return(fit_scale)
    if (colname %in% c("mu","yhat","y","estimate")) return("response")
    if (any(v < 0, na.rm = TRUE)) return("log")
    "log"
  }

  v_all <- as.numeric(pr[[fit_col]])
  sc <- infer_scale(fit_col, v_all)

  # Convert to eta (natural log mean) and mu (mean counts)
  if (sc == "response") {
    mu_all  <- pmax(v_all, 0)
    eta_all <- log(mu_all + 1e-12)
  } else {
    eta_all <- v_all
    mu_all  <- exp(eta_all)
  }

  # order by gene, lineage, time so "start" is first row within each (gene,lineage)
  if (!is.na(time_col)) {
    ord <- order(pr[[gene_col]], pr[[lin_col]], pr[[time_col]])
  } else {
    ord <- order(pr[[gene_col]], pr[[lin_col]])
  }
  pr <- pr[ord, , drop = FALSE]
  mu_all <- mu_all[ord]

  # Split by (gene,lineage)
  gl <- paste(pr[[gene_col]], pr[[lin_col]], sep = "||")
  idx <- split(seq_along(gl), gl)

  # For each (gene,lineage): compute signed max | log2((mu+pc)/(mu0+pc)) |
  per_gl <- vapply(idx, function(ii) {
    mu <- mu_all[ii]
    if (!length(mu) || all(!is.finite(mu))) return(NA_real_)

    mu0 <- mu[1L]
    log2fc_series <- log2((mu + pc) / (mu0 + pc))
    log2fc_series[ which.max(abs(log2fc_series)) ]
  }, numeric(1))

  # Pick the lineage with largest |effect| for each gene
  gene_names <- sub("\\|\\|.*$", "", names(per_gl))
  idx2 <- split(seq_along(per_gl), gene_names)
  out <- vapply(idx2, function(jj) {
    vv <- per_gl[jj]
    if (all(is.na(vv))) return(NA_real_)
    vv[ which.max(abs(vv)) ]
  }, numeric(1))

  # align to input genes
  out2 <- setNames(rep(NA_real_, length(genes)), genes)
  out2[names(out)] <- out
  out2
}

res <- signed_max_log2fc_start_many(sce, rownames(asso_df), nPoints=100, pseudocount=0.01)
asso_df$signedMaxLog2FC_start<- res[rownames(asso_df)]

### CELL 39
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
  # --- label placement without ggrepel ---
  label_x_nudge       = 0.8,
  label_y_offset_frac = 0.015,
  label_min_vgap_frac = 0.018,
  label_segment_alpha = 0.6,
  label_segment_size  = 0.4,
  # --- how to choose top-k labels ---
  label_rank_by     = c("wald","logfc"),
  label_logfc_abs   = TRUE,
  save_prefix       = NULL
) {
  label_rank_by <- match.arg(label_rank_by)

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
    ggplot2::geom_vline(xintercept = 0, linetype = "solid", size = 0.6, color = col_line) +
    { if (l2fc_thresh > 0)
        ggplot2::geom_vline(xintercept = c(-l2fc_thresh, l2fc_thresh),
                            linetype = "dotted", size = 0.7, color = col_line)
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

  # ---- dependency-free label "repel": above + spaced vertically ----
  if (nrow(label_df) > 0) {
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

  if (!is.null(save_prefix)) {
    ggplot2::ggsave(paste0(save_prefix, ".png"), p, width = 7.5, height = 6.5, dpi = 500)
    ggplot2::ggsave(paste0(save_prefix, ".pdf"), p, width = 7.5, height = 6.5)
  }
  p
}

### CELL 40
# closest gene to the qval_1 = 0.05 cutoff (by absolute distance)
i <- which.min(abs(asso_df$qval_1 - 0.01))
asso_df[i, , drop = FALSE]

# (optional) show a few nearest on both sides of the cutoff
asso_df[order(abs(asso_df$qval_1 - 0.01)), , drop = FALSE][1:10, ]

### CELL 41
options(repr.plot.width = 7.5, repr.plot.height = 6.5, repr.plot.res = 300)

# Compress the non-significant cloud strongly
p <- plot_association_volcano(
  asso_df,
  x_lim = c(-3.5, 5.5),
  l2fc_thresh = 0.5,
  y_lim = c(1, 1000),
  y_lower_scale = 0.1,   # <-- more compression below the break
  y_break = 18.31693,         # auto: 95th pct of non-sig Wald
  top_labels = 30,
  label_rank_by = "wald",      # use x-axis
  label_logfc_abs = TRUE,       # rank by magnitude
  label_only_pass = FALSE,
  label_x_nudge       = 0.05,
  label_y_offset_frac = 0.05,
  label_min_vgap_frac = 0.012,
  label_size = 4,
  axis_text_size = 18,
  axis_title_size   = 20,
  fdr_thresh = 0.01,
  save_prefix = 'l2fc05.q01'
)
print(p)

### CELL 42
options(repr.plot.width = 5.25, repr.plot.height = 3.5)

gene_symbol = 'NFKBIA'

plot_tradeSeq_gene_pub_v2(
  gene = gene_symbol, gamList = gamList, sce = sce,
  counts_sub = counts_sub, PT = PT_mat, sizeFactors = sf,
  lineages = c(1, 2), cellWeights = CW_1traj,
  colours  = c("#222021", "#222021"),
  color_points_by_lineage = TRUE,
  points_size = 0.3, points_alpha = 0.15,
)

### CELL 43
# options(repr.plot.width = 12, repr.plot.height = 8)

# library(ggplot2)
# library(patchwork)

# # ---- pick top 100 genes ----
# genes_top <- if ("gene" %in% colnames(asso_df3)) as.character(asso_df3$gene) else rownames(asso_df3)
# genes_top <- genes_top[!is.na(genes_top) & nzchar(genes_top)]
# genes_top <- unique(genes_top)[1:100]

# # ---- helper to compact a ggplot for panels ----
# .compact_panel_theme <- function(p, title_size = 8, base_size = 7, show_legend = FALSE) {
#   p <- p +
#     theme_classic(base_size = base_size) +
#     theme(
#       plot.title = element_text(size = title_size, hjust = 0.5, margin = margin(b = 1)),
#       axis.title = element_text(size = base_size),
#       axis.text  = element_text(size = base_size),
#       legend.title = element_text(size = base_size),
#       legend.text  = element_text(size = base_size),
#       plot.margin  = margin(1, 1, 1, 1)
#     )
#   if (!show_legend) p <- p + theme(legend.position = "none")
#   p
# }

# # ---- build plots ----
# plots <- lapply(genes_top, function(g) {
#   p <- plot_tradeSeq_gene_pub_v2(
#     gene = g, gamList = gamList3, sce = sce3,
#     counts_sub = counts_sub, PT = PT_mat, sizeFactors = sf,
#     lineages = c(1, 2), cellWeights = CWz_ret,
#     colours  = c("#386cb0", "#DC143C"),
#     color_points_by_lineage = TRUE,
#     points_size = 0.15, points_alpha = 0.10
#   )
#   p <- p + ggtitle(g)
#   .compact_panel_theme(p, title_size = 8, base_size = 7, show_legend = FALSE)
# })

# # ---- function to make a page grid that actually fits ----
# make_page <- function(plots, i0, n_per_page = 20, ncol = 5, collect_legend = TRUE) {
#   i1 <- min(i0 + n_per_page - 1, length(plots))
#   page_plots <- plots[i0:i1]

#   # Collect one legend for the whole page (turn legend ON only for first plot)
#   if (collect_legend) {
#     page_plots[[1]] <- .compact_panel_theme(
#       page_plots[[1]] + theme(legend.position = "right"),
#       title_size = 8, base_size = 7, show_legend = TRUE
#     )
#     wrap_plots(page_plots, ncol = ncol, guides = "collect") &
#       theme(legend.position = "right")
#   } else {
#     wrap_plots(page_plots, ncol = ncol)
#   }
# }

# # ---- show page 1 (20 genes) ----
# n_per_page <- 20
# ncol <- 5
# make_page(plots, i0 = 1, n_per_page = n_per_page, ncol = ncol, collect_legend = TRUE)

# # ---- multi-page PDF with sane sizing ----
# pdf("asso_df3_top100_tradeSeq_panels.pdf", width = 11, height = 8.5)
# for (i in seq(1, length(plots), by = n_per_page)) {
#   print(make_page(plots, i0 = i, n_per_page = n_per_page, ncol = ncol, collect_legend = TRUE))
# }
# dev.off()
# message("Wrote: asso_df3_top100_tradeSeq_panels.pdf")

### CELL 45
options(repr.plot.width = 3.75, repr.plot.height = 3.5)

plot_tradeSeq_gene_pub(
  "NFKBIA", gamList2, sce2, counts_sub, PT_mat, sf,
  cellWeights = CWz_ret,
  lineages = c(1, 2),
  colours = c("#386cb0", "#DC143C"),
  time_max = 110,
  x_tick_step = 20,        # <-- every 20 min
  color_points_by_lineage = TRUE,
  points_size = 0.75,
  points_alpha = 0.3,
  axis_text_size = 14,
  axis_title_size = 16,
  ylim = NULL,
  lineage_labels = c("Inflammatory-hi", "Erythroid-like"),
  show_legend = FALSE,
  save_pdf = FALSE
)

### CELL 46
options(repr.plot.width = 3.75, repr.plot.height = 3.5)

query_gene = 'HBA2'

plot_tradeSeq_gene_pub(
  query_gene, gamList3, sce3, counts_sub3, PT_mat, sf,
  cellWeights = CWz_ret,
  lineages = c(1,2),
  colours = c("#386cb0", "#DC143C"),
  time_max = 110,
  x_tick_step = 30,        # <-- every 20 min
  color_points_by_lineage = TRUE,
  points_size = 0,
  points_alpha = 0,
  axis_text_size = 14,
  axis_title_size = 16,
  ylim = NULL,
  lineage_labels = c("Inflammatory-hi", "Erythroid-like"),
  show_legend = FALSE,
  save_pdf = TRUE,
  pdf_filename = paste0(query_gene, '_Fig3c.pdf')
)

### CELL 49
compute_pattern_only_metrics_simple <- function(
  asso         = NULL,
  gamList      = NULL,
  sce          = NULL,
  lineages     = c(1, 2),
  nPoints      = 300,
  eps          = 1e-8,
  t_range      = NULL,        # NULL => full overlap
  pseudo_count = 0.1,
  flat_tol     = 1e-10,
  return_diag  = TRUE
) {
    
  genes = rownames(asso)
  if (is.null(genes) || !length(genes)) stop("No genes supplied / detected.")

  dm <- SummarizedExperiment::colData(sce)$tradeSeq$dm
  t_cols <- grep("^t[0-9]+$", colnames(dm), value = TRUE)
  avail  <- as.integer(sub("^t","", t_cols))
  if (!all(lineages %in% avail)) {
    stop("Requested 'lineages' not found in tradeSeq design (dm).")
  }

  # Predict curves for a lineage; return (t, mu_log1p, mu_mean)
  predict_curves <- function(g, ell) {
    gm <- gamList[[g]]
    if (is.null(gm)) {
      return(list(t = rep(NA_real_, nPoints),
                  mu_log1p = rep(NA_real_, nPoints),
                  mu_mean  = rep(NA_real_, nPoints)))
    }
    newdf <- tradeSeq:::.getPredictRangeDf(dm, lineageId = ell, nPoints = nPoints)
    tvec  <- as.numeric(newdf[[paste0("t", ell)]])
    eta   <- as.numeric(predict(gm, newdata = newdf, type = "link", se.fit = FALSE))
    mu    <- exp(eta)
    list(t = tvec, mu_log1p = log1p(mu), mu_mean = mu)
  }

  # SAME core log2FC logic as before
  signed_max_log2fc_from_start <- function(mu, pc) {
    if (!length(mu) || any(!is.finite(mu))) return(NA_real_)
    mu0 <- mu[1L]
    log2fc_series <- log2((mu + pc) / (mu0 + pc))
    log2fc_series[ which.max(abs(log2fc_series)) ]
  }

  pc <- max(0, as.numeric(pseudo_count))

  # Pre-allocate diagnostics so the column always exists
  diag_df <- if (isTRUE(return_diag)) {
    data.frame(
      gene = as.character(genes),
      m1_range = NA_real_,
      m2_range = NA_real_,
      m2_flat  = NA,
      stringsAsFactors = FALSE
    )
  } else NULL

  res_list <- vector("list", length(genes))

  for (i in seq_along(genes)) {
    g <- genes[i]

    c1 <- predict_curves(g, lineages[1])
    c2 <- predict_curves(g, lineages[2])

    # If missing/invalid predictions, keep NA row but still fill diag gene
    if (any(!is.finite(c1$t)) || any(!is.finite(c2$t))) {
      res_list[[i]] <- c(shape_sd = NA, shape_l2 = NA, shape_amp = NA,
                         mean_diff = NA, end_diff = NA, shape_ratio = NA,
                         log2FC_t1 = NA, log2FC_t2 = NA, diff_log2FC = NA)
      next
    }

    # Overlap window (no lag)
    lo <- max(min(c1$t, na.rm = TRUE), min(c2$t, na.rm = TRUE))
    hi <- min(max(c1$t, na.rm = TRUE), max(c2$t, na.rm = TRUE))

    # Optional absolute restriction
    if (!is.null(t_range) && length(t_range) == 2L && all(is.finite(t_range))) {
      lo <- max(lo, min(t_range))
      hi <- min(hi, max(t_range))
    }

    if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
      res_list[[i]] <- c(shape_sd = NA, shape_l2 = NA, shape_amp = NA,
                         mean_diff = NA, end_diff = NA, shape_ratio = NA,
                         log2FC_t1 = NA, log2FC_t2 = NA, diff_log2FC = NA)
      next
    }

    t_common <- seq(lo, hi, length.out = nPoints)

    mu1_log1p <- approx(c1$t, c1$mu_log1p, xout = t_common, ties = "ordered")$y
    mu2_log1p <- approx(c2$t, c2$mu_log1p, xout = t_common, ties = "ordered")$y
    m1        <- approx(c1$t, c1$mu_mean,  xout = t_common, ties = "ordered")$y
    m2        <- approx(c2$t, c2$mu_mean,  xout = t_common, ties = "ordered")$y

    if (any(!is.finite(mu1_log1p)) || any(!is.finite(mu2_log1p)) ||
        any(!is.finite(m1))        || any(!is.finite(m2))) {
      res_list[[i]] <- c(shape_sd = NA, shape_l2 = NA, shape_amp = NA,
                         mean_diff = NA, end_diff = NA, shape_ratio = NA,
                         log2FC_t1 = NA, log2FC_t2 = NA, diff_log2FC = NA)
      next
    }

    # pattern-only
    d         <- mu1_log1p - mu2_log1p
    mean_diff <- mean(d)
    end_diff  <- tail(d, 1)
    d_ctr     <- d - mean_diff

    shape_sd    <- stats::sd(d_ctr)
    shape_l2    <- sqrt(mean(d_ctr^2))
    shape_amp   <- max(abs(d_ctr))
    shape_ratio <- shape_sd / (abs(mean_diff) + eps)

    # SAME log2FC
    log2FC_t1   <- signed_max_log2fc_from_start(m1, pc)
    log2FC_t2   <- signed_max_log2fc_from_start(m2, pc)
    diff_log2FC <- log2FC_t1 - log2FC_t2

    res_list[[i]] <- c(shape_sd = shape_sd,
                       shape_l2 = shape_l2,
                       shape_amp = shape_amp,
                       mean_diff = mean_diff,
                       end_diff  = end_diff,
                       shape_ratio = shape_ratio,
                       log2FC_t1 = log2FC_t1,
                       log2FC_t2 = log2FC_t2,
                       diff_log2FC = diff_log2FC)

    if (isTRUE(return_diag)) {
      r1 <- diff(range(m1))
      r2 <- diff(range(m2))
      diag_df$m1_range[i] <- r1
      diag_df$m2_range[i] <- r2
      diag_df$m2_flat[i]  <- is.finite(r2) && (r2 <= flat_tol)
    }
  }

  out <- as.data.frame(do.call(rbind, res_list))
  rownames(out) <- as.character(genes)

  out$rank_shape_l2 <- rank(-out$shape_l2,    na.last = "keep", ties.method = "average")
  out$rank_shape_sd <- rank(-out$shape_sd,    na.last = "keep", ties.method = "average")
  out$rank_ratio    <- rank(-out$shape_ratio, na.last = "keep", ties.method = "average")

  if (isTRUE(return_diag)) attr(out, "diag") <- diag_df
  out
}

# ---- Run + show top pattern-only candidates -----------------------------------
pattern_only_tab2 <- compute_pattern_only_metrics_simple(
  asso          = asso_df2,
  gamList       = gamList2,
  sce           = sce2,
  lineages      = c(1, 2),
  nPoints       = 300,
  t_range       = NULL,
  pseudo_count  = 0.1,
  return_diag   = FALSE
)

combined2 <- cbind(asso_df2, pattern_only_tab2[rownames(asso_df2), ])
head(combined2)
saveRDS(combined2, 'combined2.rds')

### CELL 50
compute_single_lineage_metrics <- function(
  asso         = NULL,
  gamList      = NULL,
  sce          = NULL,
  lineage      = 1,
  nPoints      = 300,
  t_range      = NULL,      # NULL => full available range
  pseudo_count = 0.1,
  flat_tol     = 1e-10,
  return_diag  = TRUE
) {
  genes <- rownames(asso)
  if (is.null(genes) || !length(genes)) stop("No genes supplied / detected.")

  dm <- SummarizedExperiment::colData(sce)$tradeSeq$dm
  t_cols <- grep("^t[0-9]+$", colnames(dm), value = TRUE)
  avail  <- as.integer(sub("^t","", t_cols))
  if (!(lineage %in% avail)) stop("Requested 'lineage' not found in tradeSeq design (dm).")

  predict_curves <- function(g, ell) {
    gm <- gamList[[g]]
    if (is.null(gm)) {
      return(list(t = rep(NA_real_, nPoints),
                  mu_log1p = rep(NA_real_, nPoints),
                  mu_mean  = rep(NA_real_, nPoints)))
    }
    newdf <- tradeSeq:::.getPredictRangeDf(dm, lineageId = ell, nPoints = nPoints)
    tvec  <- as.numeric(newdf[[paste0("t", ell)]])
    eta   <- as.numeric(predict(gm, newdata = newdf, type = "link", se.fit = FALSE))
    mu    <- exp(eta)
    list(t = tvec, mu_log1p = log1p(mu), mu_mean = mu)
  }

  # SAME core log2FC logic as your function
  signed_max_log2fc_from_start <- function(mu, pc) {
    if (!length(mu) || any(!is.finite(mu))) return(NA_real_)
    mu0 <- mu[1L]
    log2fc_series <- log2((mu + pc) / (mu0 + pc))
    log2fc_series[ which.max(abs(log2fc_series)) ]
  }

  pc <- max(0, as.numeric(pseudo_count))

  diag_df <- if (isTRUE(return_diag)) {
    data.frame(gene = as.character(genes),
               m_range = NA_real_,
               m_flat  = NA,
               stringsAsFactors = FALSE)
  } else NULL

  res_list <- vector("list", length(genes))

  for (i in seq_along(genes)) {
    g <- genes[i]
    c1 <- predict_curves(g, lineage)

    if (any(!is.finite(c1$t))) {
      res_list[[i]] <- c(log2FC = NA, amp_log1p = NA, mean_mu = NA, end_mu = NA)
      next
    }

    lo <- min(c1$t, na.rm = TRUE); hi <- max(c1$t, na.rm = TRUE)
    if (!is.null(t_range) && length(t_range) == 2L && all(is.finite(t_range))) {
      lo <- max(lo, min(t_range)); hi <- min(hi, max(t_range))
    }
    if (!is.finite(lo) || !is.finite(hi) || hi <= lo) {
      res_list[[i]] <- c(log2FC = NA, amp_log1p = NA, mean_mu = NA, end_mu = NA)
      next
    }

    t_common <- seq(lo, hi, length.out = nPoints)
    mu_log1p <- approx(c1$t, c1$mu_log1p, xout = t_common, ties = "ordered")$y
    m        <- approx(c1$t, c1$mu_mean,  xout = t_common, ties = "ordered")$y
    if (any(!is.finite(mu_log1p)) || any(!is.finite(m))) {
      res_list[[i]] <- c(log2FC = NA, amp_log1p = NA, mean_mu = NA, end_mu = NA)
      next
    }

    # single-lineage metrics
    log2FC  <- signed_max_log2fc_from_start(m, pc)
    amp_l1p <- diff(range(mu_log1p))               # amplitude on log1p scale
    mean_mu <- mean(m)
    end_mu  <- tail(m, 1)

    res_list[[i]] <- c(log2FC = log2FC,
                       amp_log1p = amp_l1p,
                       mean_mu = mean_mu,
                       end_mu = end_mu)

    if (isTRUE(return_diag)) {
      r <- diff(range(m))
      diag_df$m_range[i] <- r
      diag_df$m_flat[i]  <- is.finite(r) && (r <= flat_tol)
    }
  }

  out <- as.data.frame(do.call(rbind, res_list))
  rownames(out) <- as.character(genes)

  # preserve asso order
  out <- out[rownames(asso), , drop = FALSE]

  # ranks (bigger = more dynamic)
  out$rank_amp   <- rank(-out$amp_log1p, na.last = "keep")
  out$rank_log2FC <- rank(-abs(out$log2FC), na.last = "keep")

  if (isTRUE(return_diag)) attr(out, "diag") <- diag_df
  out
}

# Example:
single_tab <- compute_single_lineage_metrics(
  asso = asso_df,
  gamList = gamList,
  sce = sce,
  lineage = 1,
  nPoints = 300,
  t_range = NULL,
  pseudo_count = 0.1,
  return_diag = TRUE
)


combined1 <- cbind(asso_df, single_tab[rownames(asso_df), ])
head(combined1)
saveRDS(combined1, 'combined1.rds')

### CELL 56
# load('.RData')

### CELL 57
length(gamList)
dim(sce)

### CELL 58
suppressPackageStartupMessages({
  library(tradeSeq)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(matrixStats)
  library(cluster)
})

sum(asso_df$pvalue_1 < 0.05)

### CELL 59
# -------------------------------
# 1) pick top genes by waldStat_1
# -------------------------------
topN <- 1820
top_genes_raw <- rownames(asso_df)[order(asso_df$waldStat_1, decreasing = TRUE)][seq_len(min(topN, nrow(asso_df)))]
message("Top genes (n=", length(top_genes_raw), ") e.g.: ", paste(head(top_genes_raw, 8), collapse = ", "))

### CELL 60
# -------------------------------
# 2) Map asso_df gene IDs -> features used by sce/gamList
# -------------------------------
.map_to_sce_features <- function(genes, sce) {
  feats <- rownames(sce)
  if (is.null(feats)) stop("sce has NULL rownames; cannot map genes to features.")

  mapped <- setNames(character(0), character(0))

  # direct match
  direct <- genes[genes %in% feats]
  if (length(direct) > 0) mapped <- c(mapped, setNames(direct, direct))

  remaining <- setdiff(genes, names(mapped))
  if (length(remaining) == 0) return(mapped)

  rd <- tryCatch(as.data.frame(SummarizedExperiment::rowData(sce)), error = function(e) NULL)
  if (!is.null(rd) && nrow(rd) == nrow(sce)) {
    cand_cols <- unique(c(
      intersect(c("symbol","gene","gene_name","geneSymbol","feature_name","Feature","Gene","GeneName"), colnames(rd)),
      grep("symbol|gene.*name|feature", colnames(rd), ignore.case = TRUE, value = TRUE)
    ))
    for (cc in cand_cols) {
      v <- as.character(rd[[cc]])
      idx <- match(remaining, v)
      ok <- which(!is.na(idx))
      if (length(ok) > 0) {
        newly <- feats[idx[ok]]
        names(newly) <- remaining[ok] # original -> feature
        mapped <- c(mapped, newly)
        remaining <- setdiff(genes, names(mapped))
        if (length(remaining) == 0) break
      }
    }
  }
  mapped
}

gene_map <- .map_to_sce_features(top_genes_raw, sce)
if (length(gene_map) == 0) stop("No genes mapped to sce features. Check ID mismatch.")
top_features <- unname(gene_map)
message("Mapped ", length(top_features), "/", length(top_genes_raw), " genes to sce features.")

# feature -> original label (for nicer tables)
feat_to_label <- setNames(names(gene_map), unname(gene_map))

### CELL 61
# # -------------------------------
# # 3) get fitted trajectories robustly across tradeSeq versions
# # -------------------------------
# nPoints <- 200

# .get_traj_from_predictSmooth <- function(gamList, features, nPoints = 200, lineage_pick = 1) {

#   # try tidy interface
#   pr <- tryCatch({
#     as.data.frame(tradeSeq::predictSmooth(gamList, gene = features, nPoints = nPoints, tidy = TRUE),
#                   stringsAsFactors = FALSE)
#   }, error = function(e) NULL)

#   if (!is.null(pr) && nrow(pr) > 0) {
#     gene_col <- intersect(c("gene","Gene","feature","Feature"), names(pr))[1]
#     time_col <- intersect(c("t","time","pseudotime","T"), names(pr))[1]
#     y_col    <- intersect(c("yhat","y","fit","fitted"), names(pr))[1]
#     lin_col  <- intersect(c("lineage","Lineage","l","curve"), names(pr))[1]
#     if (any(is.na(c(gene_col, time_col, y_col))))
#       stop("Tidy predictSmooth output missing required columns. Columns: ", paste(names(pr), collapse = ", "))

#     if (!is.na(lin_col)) {
#       lin_levels <- sort(unique(pr[[lin_col]]))
#       lp <- if (any(lin_levels %in% c(lineage_pick, as.character(lineage_pick))))
#         lin_levels[lin_levels %in% c(lineage_pick, as.character(lineage_pick))][1]
#       else lin_levels[1]
#       pr <- pr[pr[[lin_col]] == lp, , drop = FALSE]
#       message("Using lineage = ", lp, " (tidy)")
#     }

#     tt <- sort(unique(pr[[time_col]]))

#     traj_df <- pr %>%
#       select(feature = all_of(gene_col), t = all_of(time_col), y = all_of(y_col)) %>%
#       mutate(feature = as.character(feature)) %>%
#       arrange(feature, t)

#     keep <- traj_df %>% count(feature) %>% filter(n == length(tt)) %>% pull(feature)
#     traj_df <- traj_df %>% filter(feature %in% keep)

#     traj_mat <- traj_df %>%
#       tidyr::pivot_wider(names_from = t, values_from = y) %>%
#       tibble::column_to_rownames("feature") %>%
#       as.matrix()

#     tt_num <- suppressWarnings(as.numeric(colnames(traj_mat)))
#     if (any(!is.finite(tt_num))) tt_num <- seq(0, 1, length.out = ncol(traj_mat))

#     return(list(traj_mat = traj_mat, tt = tt_num))
#   }

#   # fallback interface
#   ps <- tradeSeq::predictSmooth(gamList, gene = features, nPoints = nPoints)
#   if (is.list(ps)) {
#     if (!is.null(ps$yhat)) ps <- ps$yhat
#     else if (length(ps) == 1) ps <- ps[[1]]
#   }

#   # If array with singleton 3rd dim, drop it cleanly
#   if (is.array(ps) && length(dim(ps)) == 3 && dim(ps)[3] == 1) {
#     ps <- ps[, , 1, drop = TRUE]
#   }

#   # Handle 3D array (time x gene x lineage) or (gene x time x lineage)
#   if (is.array(ps) && length(dim(ps)) == 3) {
#     d <- dim(ps)
#     lin <- min(lineage_pick, d[3])
#     ps2 <- ps[, , lin, drop = TRUE]
#     message("Using lineage slice = ", lin, " (array fallback)")

#     if (d[1] == nPoints) {
#       mat <- ps2              # time x gene
#     } else if (d[2] == nPoints) {
#       mat <- t(ps2)           # time x gene
#     } else {
#       mat <- ps2
#     }
#   } else if (is.matrix(ps)) {
#     d <- dim(ps)
#     if (d[1] == nPoints) {
#       mat <- ps               # time x gene
#     } else if (d[2] == nPoints) {
#       mat <- t(ps)            # time x gene
#     } else {
#       mat <- ps
#     }
#   } else {
#     stop("Unsupported predictSmooth() output type. Class: ", paste(class(ps), collapse = ", "),
#          " dim: ", paste(dim(ps), collapse = "x"))
#   }

#   # Ensure colnames are gene/feature IDs
#   if (is.null(colnames(mat)) || all(colnames(mat) == "")) {
#     if (length(features) == ncol(mat)) colnames(mat) <- features
#     else stop("No colnames in fallback matrix and length(features) != ncol(mat).")
#   }

#   # Convert to feature x time
#   traj_mat <- t(mat)
#   rownames(traj_mat) <- colnames(mat)

#   keep <- intersect(features, rownames(traj_mat))
#   if (length(keep) == 0) {
#     message("Diagnostics output feature names: ", paste(head(rownames(traj_mat), 10), collapse = ", "))
#     stop("Fallback predictSmooth produced 0 matching genes. Likely ID mismatch.")
#   }
#   traj_mat <- traj_mat[keep, , drop = FALSE]

#   # Time grid
#   tt <- attr(ps, "t")
#   if (is.null(tt)) tt <- seq(0, 1, length.out = ncol(traj_mat))

#   return(list(traj_mat = traj_mat, tt = tt))
# }

# out <- .get_traj_from_predictSmooth(gamList, top_features, nPoints = nPoints, lineage_pick = 1)
# traj_mat <- out$traj_mat
# tt <- out$tt

# # *** CRITICAL FIX: force to plain numeric matrix (no 3D arrays, no weird classes) ***
# if (is.array(traj_mat) && length(dim(traj_mat)) == 3 && dim(traj_mat)[3] == 1) {
#   traj_mat <- traj_mat[, , 1, drop = TRUE]
# }
# traj_mat <- as.matrix(traj_mat)
# storage.mode(traj_mat) <- "double"

# # sanity
# message("Valid trajectories: ", nrow(traj_mat), " genes x ", ncol(traj_mat), " timepoints.")

# orig_label <- feat_to_label[rownames(traj_mat)]
# orig_label[is.na(orig_label)] <- rownames(traj_mat)

### CELL 62
# -------------------------------
# get fitted trajectories robustly across tradeSeq versions + multi-lineage shapes
# -------------------------------
.get_traj_from_predictSmooth <- function(gamList, features, nPoints = 200, lineage_pick = 1, verbose = TRUE) {
  features <- as.character(features)

  # -------------------------------
  # helper: parse col/row names like "lineage2_17" or "lineage2.17" or "t2_17"
  # -------------------------------
  .parse_lin_idx <- function(nms) {
    if (is.null(nms)) return(NULL)
    nms <- as.character(nms)
    m <- regexec("^(?:lineage|lin|l|t)([0-9]+)[_.]([0-9]+)$", nms, perl = TRUE)
    mm <- regmatches(nms, m)
    ok <- lengths(mm) == 3
    if (!any(ok)) return(NULL)
    data.frame(
      name = nms[ok],
      lin  = as.integer(vapply(mm[ok], `[[`, "", 2)),
      idx  = as.integer(vapply(mm[ok], `[[`, "", 3)),
      stringsAsFactors = FALSE
    )
  }

  # -------------------------------
  # 1) try tidy interface
  # -------------------------------
  pr <- tryCatch({
    as.data.frame(
      tradeSeq::predictSmooth(gamList, gene = features, nPoints = nPoints, tidy = TRUE),
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)

  if (!is.null(pr) && nrow(pr) > 0) {
    gene_col <- intersect(c("gene","Gene","feature","Feature"), names(pr))[1]
    time_col <- intersect(c("t","time","pseudotime","T"), names(pr))[1]
    y_col    <- intersect(c("yhat","y","fit","fitted"), names(pr))[1]
    lin_col  <- intersect(c("lineage","Lineage","l","curve"), names(pr))[1]

    if (any(is.na(c(gene_col, time_col, y_col))))
      stop("Tidy predictSmooth output missing required columns. Columns: ", paste(names(pr), collapse = ", "))

    if (!is.na(lin_col)) {
      lin_levels <- sort(unique(pr[[lin_col]]))
      # choose requested lineage if present, else first
      lp <- if (any(lin_levels %in% c(lineage_pick, as.character(lineage_pick))))
        lin_levels[lin_levels %in% c(lineage_pick, as.character(lineage_pick))][1]
      else lin_levels[1]
      pr <- pr[pr[[lin_col]] == lp, , drop = FALSE]
      if (verbose) message("Using lineage = ", lp, " (tidy)")
    }

    tt <- sort(unique(pr[[time_col]]))

    traj_df <- data.frame(
      feature = as.character(pr[[gene_col]]),
      t       = pr[[time_col]],
      y       = pr[[y_col]],
      stringsAsFactors = FALSE
    )
    traj_df <- traj_df[order(traj_df$feature, traj_df$t), , drop = FALSE]

    # keep only genes that have all timepoints
    tab <- table(traj_df$feature)
    keep <- names(tab)[tab == length(tt)]
    traj_df <- traj_df[traj_df$feature %in% keep, , drop = FALSE]

    # pivot (base R)
    # rows: feature, cols: t
    wide <- reshape(traj_df, idvar = "feature", timevar = "t", direction = "wide")
    rn <- wide$feature
    wide$feature <- NULL
    # columns are y.<t>; strip "y."
    colnames(wide) <- sub("^y\\.", "", colnames(wide))
    traj_mat <- as.matrix(wide)
    rownames(traj_mat) <- rn

    tt_num <- suppressWarnings(as.numeric(colnames(traj_mat)))
    if (any(!is.finite(tt_num))) tt_num <- seq(0, 1, length.out = ncol(traj_mat))
    return(list(traj_mat = traj_mat, tt = tt_num))
  }

  # -------------------------------
  # 2) fallback interface (many versions / shapes)
  # -------------------------------
  ps_raw <- tryCatch(
    tradeSeq::predictSmooth(gamList, gene = features, nPoints = nPoints),
    error = function(e) e
  )
  if (inherits(ps_raw, "error")) stop("predictSmooth() failed: ", ps_raw$message)

  ps <- ps_raw
  # common: list(yhat=...)
  if (is.list(ps) && !is.null(ps$yhat)) ps <- ps$yhat
  # sometimes: list with single element
  if (is.list(ps) && length(ps) == 1) ps <- ps[[1]]

  # -------------------------------
  # case A: 3D array: time x gene x lineage  OR gene x time x lineage
  # -------------------------------
  if (is.array(ps) && length(dim(ps)) == 3) {
    d <- dim(ps)
    lin <- min(lineage_pick, d[3])
    ps2 <- ps[, , lin, drop = TRUE]
    if (verbose) message("Using lineage slice = ", lin, " (3D array)")

    # determine orientation without trusting nPoints too much
    # pick the axis that matches features more strongly
    rn <- rownames(ps2); cn <- colnames(ps2)
    ov_r <- if (!is.null(rn)) sum(rn %in% features) else 0
    ov_c <- if (!is.null(cn)) sum(cn %in% features) else 0

    # want genes x time in mat_gene_time
    mat_gene_time <- if (ov_r >= ov_c) ps2 else t(ps2)

    keep <- intersect(features, rownames(mat_gene_time))
    if (!length(keep)) {
      if (verbose) message("Diagnostics rownames: ", paste(head(rownames(mat_gene_time), 10), collapse = ", "))
      stop("No matching genes after slicing 3D array output. Likely ID mismatch.")
    }
    mat_gene_time <- mat_gene_time[keep, , drop = FALSE]

    # time grid
    tt <- attr(ps_raw, "t")
    if (is.null(tt) || length(tt) != ncol(mat_gene_time)) {
      tt <- if (ncol(mat_gene_time) > 1) (0:(ncol(mat_gene_time)-1)) / (ncol(mat_gene_time)-1) else 0
    }
    return(list(traj_mat = mat_gene_time, tt = tt))
  }

  # -------------------------------
  # case B: list of lineage-specific matrices (common in some versions)
  # -------------------------------
  if (is.list(ps) && length(ps) >= 1 && all(vapply(ps, function(x) is.matrix(x) || (is.array(x) && length(dim(x))==2), TRUE))) {
    lin <- min(lineage_pick, length(ps))
    if (verbose) message("Using element ", lin, " from list predictSmooth output (list-of-matrices).")
    ps <- ps[[lin]]
  }

  # now expect a 2D matrix/array
  if (!(is.matrix(ps) || (is.array(ps) && length(dim(ps)) == 2))) {
    stop("Unsupported predictSmooth() output type. Class: ", paste(class(ps), collapse = ", "),
         " dim: ", paste(dim(ps), collapse = "x"))
  }
  ps <- as.matrix(ps)

  # -------------------------------
  # Normalize to genes x something by choosing gene axis via overlap with 'features'
  # -------------------------------
  rn <- rownames(ps); cn <- colnames(ps)
  ov_r <- if (!is.null(rn)) sum(rn %in% features) else 0
  ov_c <- if (!is.null(cn)) sum(cn %in% features) else 0

  mat_gene_timeish <- if (ov_r >= ov_c) ps else t(ps)

  # -------------------------------
  # NEW: detect concatenated multi-lineage columns like "lineage1_1"
  # -------------------------------
  info <- .parse_lin_idx(colnames(mat_gene_timeish))
  if (!is.null(info)) {
    lin_max <- suppressWarnings(max(info$lin, na.rm = TRUE))
    lp <- min(as.integer(lineage_pick), as.integer(lin_max))
    sel <- info$name[info$lin == lp]
    ord <- sel[order(info$idx[info$lin == lp])]

    if (!length(ord)) {
      stop("Detected lineage-coded columns but found 0 columns for lineage_pick=", lineage_pick)
    }
    if (verbose) message("Using lineage = ", lp, " (concat columns like lineage", lp, "_k).")

    mat_gene_timeish <- mat_gene_timeish[, ord, drop = FALSE]

    # time grid from idx if reasonable, else 0..1
    idx <- info$idx[info$lin == lp]
    idx <- idx[order(idx)]
    if (length(idx) == ncol(mat_gene_timeish) && all(is.finite(idx)) && length(unique(idx)) > 1) {
      tt <- (idx - min(idx)) / (max(idx) - min(idx))
    } else {
      tt <- if (ncol(mat_gene_timeish) > 1) (0:(ncol(mat_gene_timeish)-1)) / (ncol(mat_gene_timeish)-1) else 0
    }

    keep <- intersect(features, rownames(mat_gene_timeish))
    if (!length(keep)) {
      if (verbose) message("Diagnostics rownames (after lineage slicing): ",
                           paste(head(rownames(mat_gene_timeish), 10), collapse = ", "))
      stop("After lineage slicing, 0 matching genes. Likely ID mismatch (symbols vs Ensembl, etc.).")
    }
    mat_gene_timeish <- mat_gene_timeish[keep, , drop = FALSE]
    return(list(traj_mat = mat_gene_timeish, tt = tt))
  }

  # -------------------------------
  # Otherwise: assume mat_gene_timeish is already genes x time
  # -------------------------------
  keep <- intersect(features, rownames(mat_gene_timeish))
  if (!length(keep)) {
    if (verbose) message("Diagnostics rownames: ", paste(head(rownames(mat_gene_timeish), 10), collapse = ", "))
    stop("Fallback predictSmooth produced 0 matching genes. Likely ID mismatch.")
  }
  mat_gene_timeish <- mat_gene_timeish[keep, , drop = FALSE]

  tt <- attr(ps_raw, "t")
  if (is.null(tt) || length(tt) != ncol(mat_gene_timeish)) {
    tt <- if (ncol(mat_gene_timeish) > 1) (0:(ncol(mat_gene_timeish)-1)) / (ncol(mat_gene_timeish)-1) else 0
  }
  return(list(traj_mat = mat_gene_timeish, tt = tt))
}

### CELL 63
nPoints

### CELL 64
out1 <- .get_traj_from_predictSmooth(gamList, top_features, nPoints = nPoints, lineage_pick = 1)
out2a <- .get_traj_from_predictSmooth(gamList2, top_features, nPoints = nPoints, lineage_pick = 1)
out2b <- .get_traj_from_predictSmooth(gamList2, top_features, nPoints = nPoints, lineage_pick = 2)

traj_mat1 <- as.matrix(out1$traj_mat)
traj_mat2a <- as.matrix(out2a$traj_mat)
traj_mat2b <- as.matrix(out2b$traj_mat)

storage.mode(traj_mat1) <- "double"
storage.mode(traj_mat2a) <- "double"
storage.mode(traj_mat2b) <- "double"

message("Valid trajectories (gamList2): ", nrow(traj_mat1), " genes x ", ncol(traj_mat1), " timepoints.")
message("Valid trajectories (gamList2): ", nrow(traj_mat2a), " genes x ", ncol(traj_mat2a), " timepoints.")
message("Valid trajectories (gamList2): ", nrow(traj_mat2b), " genes x ", ncol(traj_mat2b), " timepoints.")

orig_label <- feat_to_label[rownames(traj_mat1)]
orig_label[is.na(orig_label)] <- rownames(traj_mat1)

# L1 --> immune
# L2 --> erythroid

### CELL 65
# ---- simplest: compare one gene across traj_mat1 / traj_mat2a / traj_mat2b ----
query_gene <- "H2AC12"   # <-- change me (must match rownames(traj_mat*))

stopifnot(query_gene %in% rownames(traj_mat1),
          query_gene %in% rownames(traj_mat2a),
          query_gene %in% rownames(traj_mat2b))

tt1 <- if (exists("out1")  && !is.null(out1$tt))  out1$tt  else seq(0, 1, length.out = ncol(traj_mat1))
tt2 <- if (exists("out2a") && !is.null(out2a$tt)) out2a$tt else seq(0, 1, length.out = ncol(traj_mat2a))
tt3 <- if (exists("out2b") && !is.null(out2b$tt)) out2b$tt else seq(0, 1, length.out = ncol(traj_mat2b))

df <- rbind(
  data.frame(t = tt1, y = as.numeric(traj_mat1 [query_gene, ]), curve = "gamList (1 lineage)"),
  data.frame(t = tt2, y = as.numeric(traj_mat2a[query_gene, ]), curve = "gamList2 (L1)"),
  data.frame(t = tt3, y = as.numeric(traj_mat2b[query_gene, ]), curve = "gamList2 (L2)")
)

library(ggplot2)
ggplot(df, aes(t, y, colour = curve)) +
  geom_line(linewidth = 1.2) +
  theme_classic(base_size = 14) +
  labs(x = "Pseudotime", y = "Fitted expression", colour = NULL, title = query_gene)

### CELL 66
get_descriptors <- function(traj_mat) {
    scale_rows <- function(M) {
      mu <- rowMeans(M, na.rm = TRUE)
      sd <- matrixStats::rowSds(M, na.rm = TRUE)
      sd[sd == 0 | !is.finite(sd)] <- 1
      sweep(sweep(M, 1, mu, "-"), 1, sd, "/")
    }
    
    baseline <- traj_mat[, 1]
    traj0 <- sweep(traj_mat, 1, baseline, "-")   # SAFE baseline subtraction
    
    dt <- diff(tt)
    if (any(dt <= 0)) stop("Non-increasing time grid in tt.")
    
    dtraj <- (traj_mat[, -1, drop = FALSE] - traj_mat[, -ncol(traj_mat), drop = FALSE])
    dtraj <- sweep(dtraj, 2, dt, "/")
    
    ddt <- diff(tt[-1])
    ddtraj <- (dtraj[, -1, drop = FALSE] - dtraj[, -ncol(dtraj), drop = FALSE])
    ddtraj <- sweep(ddtraj, 2, ddt, "/")
    
    onset_delta <- 0.25
    
    peak_idx <- max.col(traj_mat, ties.method = "first")
    t_peak <- tt[peak_idx]
    y_peak <- traj_mat[cbind(seq_len(nrow(traj_mat)), peak_idx)]
    peak_amp <- y_peak - baseline
    
    onset_idx <- rep(NA_integer_, nrow(traj_mat))
    for (i in seq_len(nrow(traj_mat))) {
      thr <- baseline[i] + onset_delta
      j <- which(traj_mat[i, ] >= thr)[1]
      onset_idx[i] <- if (!is.na(j)) j else NA_integer_
    }
    t_onset <- ifelse(is.na(onset_idx), NA_real_, tt[onset_idx])
    
    k_early <- max(3, floor(0.30 * length(tt)))
    k_late  <- max(3, floor(0.30 * length(tt)))
    
    early_slope <- vapply(seq_len(nrow(traj_mat)), function(i) {
      coef(lm(traj_mat[i, 1:k_early] ~ tt[1:k_early]))[2]
    }, numeric(1))
    
    late_slope <- vapply(seq_len(nrow(traj_mat)), function(i) {
      idx <- (length(tt)-k_late+1):length(tt)
      coef(lm(traj_mat[i, idx] ~ tt[idx]))[2]
    }, numeric(1))
    
    decay_slope <- vapply(seq_len(nrow(traj_mat)), function(i) {
      pk <- peak_idx[i]
      if (pk >= length(tt)-1) return(0)
      coef(lm(traj_mat[i, pk:length(tt)] ~ tt[pk:length(tt)]))[2]
    }, numeric(1))
    
    mean_abs_curv <- rowMeans(abs(ddtraj), na.rm = TRUE)
    mid <- max(1, floor(ncol(ddtraj)/2))
    curv_early <- rowMeans(ddtraj[, 1:mid, drop = FALSE], na.rm = TRUE)
    curv_late  <- rowMeans(ddtraj[, (mid+1):ncol(ddtraj), drop = FALSE], na.rm = TRUE)
    curv_delta <- curv_late - curv_early
    
    auc_above0 <- vapply(seq_len(nrow(traj0)), function(i) {
      y <- traj0[i, ]
      sum(0.5 * (y[-1] + y[-length(y)]) * diff(tt))
    }, numeric(1))
    
    desc <- data.frame(
      feature = rownames(traj_mat),
      gene = orig_label,
      t_onset = t_onset,
      t_peak = t_peak,
      peak_amp = peak_amp,
      early_slope = early_slope,
      late_slope = late_slope,
      decay_slope = decay_slope,
      auc_above0 = auc_above0,
      mean_abs_curv = mean_abs_curv,
      curv_delta = curv_delta,
      stringsAsFactors = FALSE
    )

    return(desc)
}

desc1 = get_descriptors(traj_mat1)
desc2a = get_descriptors(traj_mat2a)
desc2b = get_descriptors(traj_mat2b)

### CELL 67
# =========================
# CONFIG
# =========================
alpha  <- 0.70      # weight on expression shape vs derivative shape (0.6–0.85 is typical)
K_grid <- 3:30
set.seed(1)

stopifnot(exists("traj_mat1"), exists("tt"), exists("asso_df"))

# =========================
# helpers
# =========================
z_row <- function(M) {
  mu <- rowMeans(M, na.rm = TRUE)
  sd <- matrixStats::rowSds(M, na.rm = TRUE)
  sd[!is.finite(sd) | sd == 0] <- 1
  sweep(sweep(M, 1, mu, "-"), 1, sd, "/")
}

# =========================
# 1) Build standardized curve + standardized derivative
# =========================
X  <- z_row(traj_mat1)

dt <- diff(tt)
Xd <- traj_mat1[, -1, drop = FALSE] - traj_mat1[, -ncol(traj_mat1), drop = FALSE]
Xd <- sweep(Xd, 2, dt, "/")
Xd <- z_row(Xd)

# sanity
stopifnot(nrow(X) == nrow(Xd))

# =========================
# 2) Correlation distances
# =========================
D_x  <- as.dist(1 - cor(t(X),  use = "pairwise.complete.obs"))
D_dx <- as.dist(1 - cor(t(Xd), use = "pairwise.complete.obs"))

# combined distance
D <- alpha * as.matrix(D_x) + (1 - alpha) * as.matrix(D_dx)
D <- as.dist(D)

# =========================
# 3) Hierarchical clustering + choose K by silhouette
# =========================
hc <- hclust(D, method = "average")

sil <- sapply(K_grid, function(K) {
  cl <- cutree(hc, k = K)
  mean(silhouette(cl, D)[, 3])
})

### CELL 68
K_best <- K_grid[which.max(sil)]
cat("Approach C: alpha =", alpha, " | silhouette-selected K =", K_best, "\n")

plot(K_grid, sil, type = "b", xlab = "K", ylab = "Mean silhouette",
     main = paste0("Approach C: K selection (alpha=", alpha, ")"))

### CELL 69
K_best = 19
clC <- factor(cutree(hc, k = K_best))

### CELL 70
# =========================
# 5) Summaries + top genes per cluster
# =========================
peak_idx <- max.col(traj_mat1, ties.method = "first")
t_peak <- tt[peak_idx]

post_slope <- vapply(seq_len(nrow(traj_mat1)), function(i) {
  pk <- peak_idx[i]
  if (pk >= ncol(traj_mat1) - 2) return(0)
  coef(lm(traj_mat1[i, pk:ncol(traj_mat1)] ~ tt[pk:ncol(traj_mat1)]))[2]
}, numeric(1))

wald <- data.frame(gene = rownames(asso_df), waldStat_1 = asso_df$waldStat_1, stringsAsFactors = FALSE)

sumC <- data.frame(
  gene = rownames(traj_mat1),
  cluster = clC,
  t_peak = t_peak,
  post_slope = post_slope
) %>%
  left_join(wald, by = "gene") %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    med_t_peak = median(t_peak, na.rm = TRUE),
    med_post_slope = median(post_slope, na.rm = TRUE),
    .groups = "drop"
  )
print(sumC)

topC <- data.frame(
  gene = rownames(traj_mat1),
  cluster = clC,
  t_peak = t_peak,
  post_slope = post_slope
) %>%
  left_join(wald, by = "gene") %>%
  group_by(cluster) %>%
  arrange(desc(waldStat_1), t_peak) %>%
  slice_head(n = 12) %>%
  ungroup()

print(topC, n = Inf)

# Save labels for downstream
cluster_C <- clC

### CELL 71
stopifnot(exists("traj_mat1"), exists("tt"), exists("asso_df"), exists("cluster_C"))

# -------------------------
# CONFIG
# -------------------------
alpha_mix <- 0.70     # must match Approach C distance definition
top_n <- 19
show_sil_text <- TRUE

# -------------------------
# Helpers
# -------------------------
z_row <- function(M) {
  mu <- rowMeans(M, na.rm = TRUE)
  sd <- matrixStats::rowSds(M, na.rm = TRUE)
  sd[!is.finite(sd) | sd == 0] <- 1
  sweep(sweep(M, 1, mu, "-"), 1, sd, "/")
}

# -------------------------
# Standardized curves
# -------------------------
X <- z_row(traj_mat1)

# -------------------------
# Distance matrix D (recompute if needed)
# -------------------------
if (!exists("D")) {
  dt <- diff(tt)
  Xd <- traj_mat1[, -1, drop = FALSE] - traj_mat1[, -ncol(traj_mat1), drop = FALSE]
  Xd <- sweep(Xd, 2, dt, "/")
  Xd <- z_row(Xd)

  D_x  <- as.dist(1 - cor(t(X),  use = "pairwise.complete.obs"))
  D_dx <- as.dist(1 - cor(t(Xd), use = "pairwise.complete.obs"))
  D <- alpha_mix * as.matrix(D_x) + (1 - alpha_mix) * as.matrix(D_dx)
  D <- as.dist(D)
  cat("Recomputed D with alpha =", alpha_mix, "\n")
}

# -------------------------
# Silhouette width per gene = membership confidence
# -------------------------
cl <- factor(cluster_C)
sil_obj <- silhouette(as.integer(cl), D)
sil_df <- data.frame(
  gene = rownames(traj_mat1),
  sil_width = sil_obj[, "sil_width"],
  stringsAsFactors = FALSE
)

# Map silhouette to opacity:
# clamp [-0.2, 0.8] to [0.15, 1.0] for readability
sil_to_alpha <- function(s) {
  s <- pmax(-0.2, pmin(0.8, s))
  0.15 + (s + 0.2) / (0.8 + 0.2) * (1.0 - 0.15)
}

# -------------------------
# Top genes per cluster (by waldStat_1), + medoid
# -------------------------
wald <- data.frame(gene = rownames(asso_df), waldStat_1 = asso_df$waldStat_1, stringsAsFactors = FALSE)
gene_info <- data.frame(gene = rownames(traj_mat1), cluster = cl, stringsAsFactors = FALSE) %>%
  left_join(wald, by = "gene") %>%
  left_join(sil_df, by = "gene")

Dmat <- as.matrix(D)

medoids <- gene_info %>%
  group_by(cluster) %>%
  summarise(
    medoid = {
      g <- gene
      subD <- Dmat[g, g, drop = FALSE]
      g[which.min(rowSums(subD))]
    },
    .groups = "drop"
  )

top_genes <- gene_info %>%
  group_by(cluster) %>%
  arrange(desc(waldStat_1), gene, .by_group = TRUE) %>%
  slice_head(n = top_n) %>%
  summarise(top = list(gene), .groups = "drop") %>%
  left_join(medoids, by = "cluster") %>%
  mutate(panel_genes = purrr::map2(medoid, top, ~ c(paste0("MEDOID:", .x), .y)))

panel_tbl <- top_genes %>%
  select(cluster, panel_genes) %>%
  tidyr::unnest(panel_genes) %>%
  group_by(cluster) %>%
  mutate(panel_pos = row_number()) %>%
  ungroup() %>%
  mutate(
    panel_label = ifelse(grepl("^MEDOID:", panel_genes),
                         sub("^MEDOID:", "", panel_genes),
                         panel_genes),
    is_medoid = grepl("^MEDOID:", panel_genes),
    gene = panel_label
  )

panel_gene_set <- unique(panel_tbl$gene)
Xsub <- X[panel_gene_set, , drop = FALSE]

plot_df <- as.data.frame(Xsub) %>%
  mutate(gene = rownames(Xsub)) %>%
  pivot_longer(cols = -gene, names_to = "tp", values_to = "z") %>%
  mutate(t_idx = as.integer(gsub(".*_", "", tp)),
         t = tt[t_idx]) %>%
  inner_join(panel_tbl %>% select(cluster, gene, panel_pos, is_medoid), by = "gene") %>%
  left_join(gene_info %>% select(gene, sil_width), by = "gene") %>%
  mutate(
    panel = factor(panel_pos, levels = seq_len(top_n + 1)),
    alpha = sil_to_alpha(sil_width),
    sil_lab = sprintf("sil=%.2f", sil_width)
  )

# -------------------------
# Colors: if exactly 3 clusters -> red/blue/green
# -------------------------
clus_levels <- levels(cl)
if (length(clus_levels) == 3) {
  pal <- setNames(c("#1f77b4", "#d62728", "#2ca02c"), clus_levels)  # blue, red, green
  # If you want specifically "red, blue, green" in that order, swap as desired:
  # pal <- setNames(c("#d62728", "#1f77b4", "#2ca02c"), clus_levels)
} else {
  pal <- NULL  # default ggplot hues
}

# figure sizing
n_clusters <- nlevels(cl)
options(repr.plot.width = max(14, 1.55 * (top_n + 1)),
        repr.plot.height = max(4, 1.55 * n_clusters))

p <- ggplot(plot_df, aes(x = t, y = z, group = interaction(gene, cluster))) +
  geom_hline(yintercept = 0, linewidth = 0.25, alpha = 0.35) +
  geom_line(aes(color = cluster, alpha = alpha), linewidth = 0.95) +
  facet_grid(cluster ~ panel, switch = "y") +
  scale_alpha_identity() +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, size = 11, face = "bold"),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 11),
    panel.spacing.x = unit(0.25, "lines"),
    panel.spacing.y = unit(0.35, "lines"),
    legend.position = "none"
  ) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
  labs(
    title = paste0("Approach C clusters: medoid + top ", top_n,
                   " genes (color=cluster; opacity=membership confidence via silhouette)"),
    x = "Time (0→1)",
    y = "z-scored fitted expression"
  )

if (!is.null(pal)) p <- p + scale_color_manual(values = pal)

if (show_sil_text) {
  # Put silhouette text in the upper-left of each panel
  ann <- plot_df %>%
    group_by(cluster, panel, gene) %>%
    summarise(alpha = first(alpha), sil_lab = first(sil_lab), .groups = "drop") %>%
    mutate(t = 0.02, z = 2.8)
  p <- p + geom_text(
    data = ann,
    aes(x = t, y = z, label = sil_lab, alpha = alpha, color = cluster),
    size = 2.2, hjust = 0, vjust = 1
  )
}

print(p)

### CELL 72
options(repr.plot.width = 1.5, repr.plot.height = 1.5)

# =========================
# REQUIREMENTS
# =========================
stopifnot(exists("traj_mat1"), exists("tt"), exists("cluster_C"))

# =========================
# CONFIG
# =========================
alpha_mix      <- 0.70     # must match Approach C
npc_umap       <- 10       # FPCA PCs to feed into UMAP
umap_neighbors <- 25
umap_min_dist  <- 0.15
set.seed(1)

# =========================
# Helpers
# =========================
z_row <- function(M) {
  mu <- rowMeans(M, na.rm = TRUE)
  sd <- matrixStats::rowSds(M, na.rm = TRUE)
  sd[!is.finite(sd) | sd == 0] <- 1
  sweep(sweep(M, 1, mu, "-"), 1, sd, "/")
}

sil_to_alpha <- function(s) {
  # map silhouette [-0.2, 0.8] -> alpha [0.15, 1.0]
  s <- pmax(-0.2, pmin(0.8, s))
  0.15 + (s + 0.2) / (1.0) * (1.0 - 0.15)
}

# =========================
# Standardized curves (genes x time)
# =========================
X <- z_row(traj_mat1)
genes <- rownames(traj_mat1)
cl <- factor(cluster_C)
stopifnot(length(cl) == nrow(X))

# =========================
# Distance D (recompute if missing) + silhouette
# =========================
if (!exists("D")) {
  dt <- diff(tt)
  Xd <- traj_mat1[, -1, drop = FALSE] - traj_mat1[, -ncol(traj_mat1), drop = FALSE]
  Xd <- sweep(Xd, 2, dt, "/")
  Xd <- z_row(Xd)

  D_x  <- as.dist(1 - cor(t(X),  use = "pairwise.complete.obs"))
  D_dx <- as.dist(1 - cor(t(Xd), use = "pairwise.complete.obs"))
  D <- alpha_mix * as.matrix(D_x) + (1 - alpha_mix) * as.matrix(D_dx)
  D <- as.dist(D)
  cat("Recomputed D with alpha =", alpha_mix, "\n")
}

sil_obj <- silhouette(as.integer(cl), D)
sil_df <- data.frame(
  gene = genes,
  cluster = cl,
  sil_width = sil_obj[, "sil_width"],
  alpha = sil_to_alpha(sil_obj[, "sil_width"]),
  stringsAsFactors = FALSE
)

# ============================================================
# (1) PCA on standardized curves
# ============================================================
pc <- prcomp(X, center = FALSE, scale. = FALSE)
pc_df <- data.frame(
  gene = genes,
  cluster = cl,
  PC1 = pc$x[, 1],
  PC2 = pc$x[, 2],
  PC3 = pc$x[, 3]
) %>% left_join(sil_df, by = c("gene","cluster"))

# Add an interpretable kinetic feature: t_peak on the standardized curve
peak_idx <- max.col(X, ties.method = "first")
pc_df$t_peak <- tt[peak_idx]

# ============================================================
# (2) MDS (cmdscale) on D (distance used by clustering)
# ============================================================
mds <- cmdscale(D, k = 2, eig = TRUE)
mds_df <- data.frame(
  gene = genes,
  cluster = cl,
  MDS1 = mds$points[, 1],
  MDS2 = mds$points[, 2]
) %>% left_join(sil_df, by = c("gene","cluster"))

### CELL 73
# ---- 2x2 panels, with optional top-k gene labels + colorbars for feature plots ----
suppressPackageStartupMessages({
  library(ggplot2)
})

# ========= label controls =========
SHOW_LABELS   <- TRUE
TOPK_LABELS   <- 100
LABEL_SIZE    <- 3.6
LABEL_NUDGE_X <- 0.02   # only used if ggrepel not installed
LABEL_NUDGE_Y <- 0.02

# ========= requirements =========
stopifnot(exists("pc_df"), exists("mds_df"))
stopifnot(all(c("gene","PC1","PC2") %in% names(pc_df)))
stopifnot(all(c("MDS1","MDS2") %in% names(mds_df)))
stopifnot(all(c("cluster","alpha","t_peak","log2FC") %in% names(pc_df)))

# ========= color controls for log2FC =========
FC_CAP      <- 0.5   # hard cap +/-2
FC_WHITE_IN <- 0     # values with |log2FC| <= 1 look ~white

# Clip log2FC just for coloring (keeps original pc_df$log2FC unchanged)
pc_df$log2FC_clip <- pmax(-FC_CAP, pmin(FC_CAP, pc_df$log2FC))

# choose which genes to label (top-k rows of pc_df by default)
genes_lab <- if (isTRUE(SHOW_LABELS)) {
  unique(as.character(pc_df$gene))[seq_len(min(TOPK_LABELS, nrow(pc_df)))]
} else character(0)

# helper to add labels if requested (uses ggrepel if available, else geom_text)
add_labels <- function(p, df, x, y, gene_col = "gene") {
  if (!isTRUE(SHOW_LABELS) || length(genes_lab) == 0) return(p)
  dlab <- df[df[[gene_col]] %in% genes_lab, , drop = FALSE]
  if (!nrow(dlab)) return(p)

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p + ggrepel::geom_text_repel(
      data = dlab,
      aes_string(x = x, y = y, label = gene_col),
      size = LABEL_SIZE,
      max.overlaps = Inf,
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0
    )
  } else {
    p + geom_text(
      data = dlab,
      aes_string(x = x, y = y, label = gene_col),
      size = LABEL_SIZE,
      nudge_x = LABEL_NUDGE_X,
      nudge_y = LABEL_NUDGE_Y,
      check_overlap = TRUE
    )
  }
}

# ========= plots =========
p1 <- ggplot(pc_df, aes(PC1, PC2, color = cluster, alpha = alpha)) +
  geom_point(size = 1.0) +
  scale_alpha_identity() +
  theme_classic(base_size = 10) +
  theme(legend.position = "none") +
  labs(title = "PCA (cluster)", x = "PC1", y = "PC2")
p1 <- add_labels(p1, pc_df, "PC1", "PC2")

p4 <- ggplot(mds_df, aes(MDS1, MDS2, color = cluster, alpha = alpha)) +
  geom_point(size = 1.0) +
  scale_alpha_identity() +
  theme_classic(base_size = 10) +
  theme(legend.position = "none") +
  labs(title = "MDS (cluster)", x = "MDS1", y = "MDS2")
p4 <- add_labels(p4, mds_df, "MDS1", "MDS2")

# Feature plot: t_peak (KEEP a colorbar)  ---- 3-color gradient ----
p2 <- ggplot(pc_df, aes(PC1, PC2, color = t_peak, alpha = alpha)) +
  geom_point(size = 1.0) +
  scale_alpha_identity() +
  theme_classic(base_size = 10) +
  theme(legend.position = "right",
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 8)) +
  labs(title = "PCA (t_peak)", x = "PC1", y = "PC2", color = "t_peak")

tp_rng <- range(pc_df$t_peak, na.rm = TRUE)
tp_lo  <- tp_rng[1]; tp_hi <- tp_rng[2]
tp_mid <- mean(tp_rng)

p2 <- p2 + scale_color_gradientn(
  colours = c("#2c7bb6", "#fdae61", "#d7191c"),   # early, mid, late
  values  = c(0, 0.5, 1),
  limits  = c(tp_lo, tp_hi),
  oob     = if (requireNamespace("scales", quietly = TRUE)) scales::squish else scales::squish,
  breaks  = c(tp_lo, tp_mid, tp_hi),
  labels  = c("early", "mid", "late")
)
p2 <- add_labels(p2, pc_df, "PC1", "PC2")

# Feature plot: log2FC (KEEP a colorbar; cap +/-2; mostly white for |log2FC| < 1)
p3 <- ggplot(pc_df, aes(PC1, PC2, color = log2FC_clip)) +
  geom_point(size = 1.0, alpha = 1.0) +   # fully opaque; no alpha mapping
  theme_classic(base_size = 10) +
  theme(legend.position = "right",
        legend.title = element_text(size = 9),
        legend.text  = element_text(size = 8)) +
  labs(title = "PCA (log2FC)", x = "PC1", y = "PC2", color = "log2FC")

if (requireNamespace("scales", quietly = TRUE)) {
  v_mid_low  <- (FC_WHITE_IN - (-FC_CAP)) / (2 * FC_CAP)
  v_mid_high <- (FC_WHITE_IN + FC_CAP) / (2 * FC_CAP)

  p3 <- p3 + scale_color_gradientn(
    colours = c("#2c7bb6", "grey95", "grey95", "#d7191c"),
    values  = c(0, v_mid_low, v_mid_high, 1),
    limits  = c(-FC_CAP, FC_CAP),
    oob     = scales::squish,
    breaks  = c(-FC_CAP, 0, FC_CAP)
  )
} else {
  p3 <- p3 + scale_color_gradient2(
    low = "#2c7bb6", mid = "grey95", high = "#d7191c",
    midpoint = 0,
    limits = c(-FC_CAP, FC_CAP)
  )
}
p3 <- add_labels(p3, pc_df, "PC1", "PC2")

# ========= layout (TOP: cluster plots, BOTTOM: feature plots) =========
options(repr.plot.width = 13.5, repr.plot.height = 9.2)

if (requireNamespace("patchwork", quietly = TRUE)) {
  print((p1 + p4) / (p2 + p3) + patchwork::plot_layout(guides = "collect"))
} else if (requireNamespace("gridExtra", quietly = TRUE)) {
  gridExtra::grid.arrange(p1, p4, p2, p3, ncol = 2)
} else {
  print(p1); print(p4); print(p2); print(p3)
}

### CELL 76
# =========================
# Compare ORIGINAL PCA vs OPTION B (Z-curves + log2FC feature)
# - ORIGINAL: PCA on z-scored curves only (your current default)
# - OPTION B: PCA on z-scored curves with log2FC appended as an extra feature
# Both visualized with the same 2x2 PCA-only panels (cluster/t_peak/log2FC), ignoring MDS
# =========================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# =========================
# REQUIREMENTS
# =========================
stopifnot(exists("traj_mat1"), exists("cluster_C"), exists("sil_df"), exists("tt"))
stopifnot(exists("z_row"))
stopifnot(exists("combined1"), "log2FC" %in% colnames(combined1))

# =========================
# CONFIG
# =========================
SHOW_LABELS   <- TRUE
TOPK_LABELS   <- 100
LABEL_SIZE    <- 3.2
LABEL_NUDGE_X <- 0.02
LABEL_NUDGE_Y <- 0.02

# t_peak: 3-color gradient anchors (early/mid/late)
TPEAK_COLS <- c("early" = "#2c7bb6", "mid" = "#fdae61", "late" = "#d7191c")

# log2FC coloring: cap +/-2 and mostly white for |log2FC| <= 1
FC_CAP      <- 2
FC_WHITE_IN <- 1

# PCA axes
PCX <- "PC1"
PCY <- "PC2"

# Figure size (per 2x2 set)
options(repr.plot.width = 13.5, repr.plot.height = 9.2)

# =========================
# Helpers
# =========================
genes_all <- rownames(traj_mat1)
stopifnot(!is.null(genes_all), length(genes_all) > 0)

genes_lab <- if (isTRUE(SHOW_LABELS)) {
  unique(as.character(genes_all))[seq_len(min(TOPK_LABELS, length(genes_all)))]
} else character(0)

add_labels <- function(p, df, x, y, gene_col = "gene") {
  if (!isTRUE(SHOW_LABELS) || length(genes_lab) == 0) return(p)
  dlab <- df[df[[gene_col]] %in% genes_lab, , drop = FALSE]
  if (!nrow(dlab)) return(p)

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p + ggrepel::geom_text_repel(
      data = dlab,
      aes_string(x = x, y = y, label = gene_col),
      size = LABEL_SIZE,
      max.overlaps = Inf,
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0
    )
  } else {
    p + geom_text(
      data = dlab,
      aes_string(x = x, y = y, label = gene_col),
      size = LABEL_SIZE,
      nudge_x = LABEL_NUDGE_X,
      nudge_y = LABEL_NUDGE_Y,
      check_overlap = TRUE
    )
  }
}

plot_pca_2x2 <- function(pc_df, title_prefix = "") {
  # t_peak 3-color gradient anchors
  tp_rng <- range(pc_df$t_peak, na.rm = TRUE)
  tp_mid <- mean(tp_rng)
  tp_vals <- scales::rescale(c(tp_rng[1], tp_mid, tp_rng[2]))

  # log2FC color prep
  pc_df$log2FC_clip <- pmax(-FC_CAP, pmin(FC_CAP, pc_df$log2FC))
  v_mid_low  <- (FC_WHITE_IN - (-FC_CAP)) / (2 * FC_CAP)
  v_mid_high <- (FC_WHITE_IN + FC_CAP) / (2 * FC_CAP)

  p_cluster <- ggplot(pc_df, aes_string(PCX, PCY, color = "cluster", alpha = "alpha")) +
    geom_point(size = 1.0) +
    scale_alpha_identity() +
    theme_classic(base_size = 10) +
    theme(legend.position = "none") +
    labs(title = paste0(title_prefix, "PCA (cluster)"), x = PCX, y = PCY)
  p_cluster <- add_labels(p_cluster, pc_df, PCX, PCY)

  p_tpeak <- ggplot(pc_df, aes_string(PCX, PCY, color = "t_peak", alpha = "alpha")) +
    geom_point(size = 1.0) +
    scale_alpha_identity() +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8)
    ) +
    labs(title = paste0(title_prefix, "PCA (t_peak)"), x = PCX, y = PCY, color = "t_peak") +
    scale_color_gradientn(colours = unname(TPEAK_COLS), values = tp_vals)
  p_tpeak <- add_labels(p_tpeak, pc_df, PCX, PCY)

  p_fc <- ggplot(pc_df, aes_string(PCX, PCY, color = "log2FC_clip")) +
    geom_point(size = 1.0, alpha = 1.0) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8)
    ) +
    labs(title = paste0(title_prefix, "PCA (log2FC)"), x = PCX, y = PCY, color = "log2FC") +
    scale_color_gradientn(
      colours = c("#2c7bb6", "grey95", "grey95", "#d7191c"),
      values  = c(0, v_mid_low, v_mid_high, 1),
      limits  = c(-FC_CAP, FC_CAP),
      oob     = scales::squish,
      breaks  = c(-FC_CAP, -1, 0, 1, FC_CAP)
    )
  p_fc <- add_labels(p_fc, pc_df, PCX, PCY)

  p_blank <- ggplot() + theme_void() + labs(title = "")

  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_cluster <- p_cluster + theme(legend.position = "none")
    print((p_cluster + p_tpeak) / (p_fc + p_blank) + patchwork::plot_layout(guides = "collect"))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(p_cluster, p_tpeak, p_fc, p_blank, ncol = 2)
  } else {
    print(p_cluster); print(p_tpeak); print(p_fc)
  }
}

# =========================
# Build shared metadata
# =========================
cl <- factor(cluster_C)[match(genes_all, rownames(traj_mat1))]
stopifnot(!anyNA(cl))

meta_df <- data.frame(
  gene    = genes_all,
  cluster = cl,
  stringsAsFactors = FALSE
) %>%
  left_join(sil_df, by = c("gene","cluster"))

meta_df$log2FC <- as.numeric(combined1[meta_df$gene, "log2FC"])

# =========================
# ORIGINAL PCA: z-scored curves only
# =========================
X_z <- z_row(traj_mat1)
storage.mode(X_z) <- "double"

pc_orig <- prcomp(X_z, center = FALSE, scale. = FALSE)

pc_df_orig <- meta_df
pc_df_orig$PC1 <- pc_orig$x[, 1]
pc_df_orig$PC2 <- pc_orig$x[, 2]
pc_df_orig$PC3 <- pc_orig$x[, 3]

peak_idx_z <- max.col(X_z, ties.method = "first")
pc_df_orig$t_peak <- tt[peak_idx_z]



# =========================
# OPTION B: z-scored curves + append log2FC as extra feature
#   NEW: allow weighting of log2FC_feature
# =========================

# ---- knob: target fraction of total variance attributed to log2FC_feature ----
# e.g. 0.10 (subtle), 0.25 (noticeable), 0.50 (strong)
FC_FRACTION <- 0.15

fc <- pc_df_orig$log2FC

# z-score log2FC across genes (same as before)
fc_z <- fc
if (any(is.finite(fc))) {
  s <- stats::sd(fc, na.rm = TRUE); if (!is.finite(s) || s == 0) s <- 1
  m <- mean(fc, na.rm = TRUE)
  fc_z <- (fc - m) / s
}
fc_z[!is.finite(fc_z)] <- 0

# ---- compute an automatic weight so log2FC contributes ~FC_FRACTION of variance ----
# PCA on an uncentered matrix still depends on second moments; this variance proxy is fine here.
v_time <- sum(apply(X_z, 2, stats::var, na.rm = TRUE))     # total variance across the 200 time features
v_fc   <- stats::var(fc_z, na.rm = TRUE)                   # variance of z-scored log2FC feature
if (!is.finite(v_fc) || v_fc <= 0) v_fc <- 1

# Solve: (w^2 v_fc) / (v_time + w^2 v_fc) = FC_FRACTION
FC_FRACTION <- min(max(FC_FRACTION, 1e-6), 1 - 1e-6)
FC_WEIGHT   <- sqrt((FC_FRACTION * v_time) / ((1 - FC_FRACTION) * v_fc))

fc_w <- FC_WEIGHT * fc_z

# build augmented matrix + PCA (same as before)
X_z_plus <- cbind(X_z, log2FC_feature = fc_w)
pc_b <- prcomp(X_z_plus, center = FALSE, scale. = FALSE)

pc_df_b <- meta_df
pc_df_b$PC1 <- pc_b$x[, 1]
pc_df_b$PC2 <- pc_b$x[, 2]
pc_df_b$PC3 <- pc_b$x[, 3]
pc_df_b$t_peak <- tt[peak_idx_z]  # keep comparable definition

cat(sprintf("Option B weighting: FC_FRACTION=%.3f, FC_WEIGHT=%.3f (v_time=%.3f, v_fc=%.3f)\n",
            FC_FRACTION, FC_WEIGHT, v_time, v_fc))

# =========================
# Plot comparison: ORIGINAL vs OPTION B (two 2x2 blocks)
# =========================
if (requireNamespace("patchwork", quietly = TRUE)) {
  pA <- patchwork::wrap_elements(full = grid::textGrob("")) # placeholder
}

# Render ORIGINAL then OPTION B
plot_pca_2x2(pc_df_orig, title_prefix = "ORIGINAL (z-curves only): ")
plot_pca_2x2(pc_df_b,    title_prefix = "OPTION B (z-curves + log2FC feat): ")

### CELL 77
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ---------- REQUIREMENTS ----------
stopifnot(exists("traj_mat1"), exists("tt"), exists("z_row"))
stopifnot(exists("pc_orig"), exists("pc_b"))

Tn <- ncol(traj_mat1)
stopifnot(length(tt) == Tn)

# ---------- CONFIG ----------
set.seed(1)

NOISE_SD_SET   <- c(0.01, 0.03, 0.06)
RAMP_MAGS      <- c(0.06, 0.12, 0.24)
PULSE_MAGS     <- c(0.10, 0.30, 0.50, 0.80, 1.00)
PULSE_CENTERS  <- c(0.20, 0.50, 0.80)
PULSE_WIDTHS   <- c(0.1)

NEARFLAT_FC_JITTER_FRAC <- 0.10

POINT_SIZE_BG   <- 0.8
ALPHA_BG        <- 0.18
POINT_SIZE_TEST <- 2.8
ALPHA_TEST      <- 0.95

# choose how to color synthetic points:
#   "abs"   -> unsigned magnitude
#   "signed"-> signed log2FC_feature (direction)
COLOR_MODE <- "signed"

options(repr.plot.width = 15.5, repr.plot.height = 6.2)

# ---------- helpers ----------
make_flat <- function(level = 1) rep(level, Tn)
make_nearflat <- function(level = 1, sd = 0.03) level + rnorm(Tn, 0, sd)

make_ramp <- function(level = 1, mag = 0.15, up = TRUE) {
  x <- seq(-0.5, 0.5, length.out = Tn)
  level + (if (up) +1 else -1) * mag * x
}

# signed pulse: sign=+1 bump, sign=-1 dip
make_pulse <- function(level = 1, mag = 0.20, center = 0.5, width = 0.12, sign = +1) {
  x <- tt
  level + (sign * mag) * exp(-0.5 * ((x - center) / width)^2)
}

apply_train_preproc <- function(X, pc_obj) {
  Xp <- X
  if (!is.null(pc_obj$center)) {
    if (length(pc_obj$center) == ncol(Xp)) Xp <- sweep(Xp, 2, pc_obj$center, "-")
  }
  if (!is.null(pc_obj$scale)) {
    if (length(pc_obj$scale) == ncol(Xp)) {
      sc <- pc_obj$scale
      bad <- (!is.finite(sc)) | (sc == 0)
      if (any(bad)) sc[bad] <- 1
      Xp <- sweep(Xp, 2, sc, "/")
    }
  }
  Xp
}

# Parse features like "lineage1_17" / "lineage2_003" / "t17" / "17"
parse_time_any <- function(nm) {
  m <- regmatches(nm, regexpr("[0-9]+$", nm))
  suppressWarnings(as.integer(m))
}
parse_lineage_any <- function(nm) {
  m <- regmatches(nm, regexpr("lineage[0-9]+", nm))
  if (length(m) == 0 || is.na(m) || m == "") return(1L)
  suppressWarnings(as.integer(gsub("^lineage", "", m)))
}

# ---------- build synthetic library ----------
tests <- list()
meta  <- list()

add_test <- function(x, family, mag, sign = 0, extra = "") {
  tests[[length(tests) + 1]] <<- x
  meta[[length(meta) + 1]]  <<- data.frame(
    family = family,
    mag = mag,
    sign = sign,
    extra = extra,
    stringsAsFactors = FALSE
  )
}

# flat
for (k in 1:8) add_test(make_flat(level = runif(1, 0.3, 2.5)),
                        family = "flat", mag = 0, sign = 0)

# nearflat
for (sd in NOISE_SD_SET) {
  for (k in 1:6) add_test(make_nearflat(level = runif(1, 0.3, 2.5), sd = sd),
                          family = "nearflat_noise", mag = sd, sign = 0, extra = paste0("sd=", sd))
}

# ramps
for (mag in RAMP_MAGS) {
  for (k in 1:4) add_test(make_ramp(level = runif(1, 0.3, 2.5), mag = mag, up = TRUE),
                          family = "ramp_up", mag = mag, sign = +1)
  for (k in 1:4) add_test(make_ramp(level = runif(1, 0.3, 2.5), mag = mag, up = FALSE),
                          family = "ramp_down", mag = mag, sign = -1)
}

# pulses: BOTH positive bumps and negative dips
for (mag in PULSE_MAGS) {
  for (cc in PULSE_CENTERS) {
    for (wd in PULSE_WIDTHS) {
      add_test(make_pulse(level = runif(1, 0.3, 2.5), mag = mag, center = cc, width = wd, sign = +1),
               family = "pulse_up", mag = mag, sign = +1, extra = paste0("c=", cc, ", w=", wd))
      add_test(make_pulse(level = runif(1, 0.3, 2.5), mag = mag, center = cc, width = wd, sign = -1),
               family = "pulse_down", mag = mag, sign = -1, extra = paste0("c=", cc, ", w=", wd))
    }
  }
}

X_test_raw <- do.call(rbind, tests)
meta_df <- do.call(rbind, meta)
rownames(X_test_raw) <- paste0("test_", seq_len(nrow(X_test_raw)))
meta_df$gene <- rownames(X_test_raw)

# z-score rows (ORIGINAL PCA convention)
X_test_z <- z_row(X_test_raw)
storage.mode(X_test_z) <- "double"

# ---------- ORIGINAL PCA projection ----------
S_test_orig <- X_test_z %*% pc_orig$rotation
test_df_orig <- cbind(meta_df, PC1 = S_test_orig[,1], PC2 = S_test_orig[,2])
real_df_orig <- if (exists("pc_df_orig")) pc_df_orig else data.frame(PC1 = pc_orig$x[,1], PC2 = pc_orig$x[,2])

# ---------- OPTION B: build feature matrix aligned to pc_b$rotation ----------
rot_names <- rownames(pc_b$rotation)
stopifnot(!is.null(rot_names))

# Identify log2FC feature row name
fc_cands <- c("log2FC_feature","log2FC","log2fc","logFC","fc")
fc_name <- intersect(fc_cands, rot_names)[1]

time_idx <- parse_time_any(rot_names)
is_time  <- is.finite(time_idx) & (time_idx >= 1) & (time_idx <= Tn)

time_cols <- rot_names[is_time]
time_map  <- time_idx[is_time]
time_lin  <- vapply(time_cols, parse_lineage_any, integer(1))  # not used, but kept for debugging

cat("Option B rotation dims:", paste(dim(pc_b$rotation), collapse=" x "), "\n")
cat("Parsed time-like features:", length(time_cols), "of", length(rot_names), "\n")
cat("Example time feature names:", paste(utils::head(time_cols, 8), collapse=", "), "\n")

if (length(time_cols) == 0) stop("Found 0 time-like rotation features; inspect head(rot_names).")

if (is.na(fc_name) || is.null(fc_name)) {
  non_time <- rot_names[!is_time]
  if (length(non_time) < 1) stop("No non-time features found to use as log2FC feature.")
  fc_name <- non_time[1]
  cat("WARNING: could not match log2FC feature name; using first non-time feature:", fc_name, "\n")
} else {
  cat("Option B feature name used for log2FC:", fc_name, "\n")
}

# Signed log2FC_feature for synthetic points
meta_df$log2FC_feature <- 0
is_nearflat <- meta_df$family == "nearflat_noise"
if (any(is_nearflat)) {
  meta_df$log2FC_feature[is_nearflat] <- runif(sum(is_nearflat), -1, 1) *
    (NEARFLAT_FC_JITTER_FRAC * meta_df$mag[is_nearflat])
}
idx_signed <- meta_df$family %in% c("ramp_up","ramp_down","pulse_up","pulse_down")
meta_df$log2FC_feature[idx_signed] <- meta_df$sign[idx_signed] * meta_df$mag[idx_signed]

# Build X_b with correct columns (MUST match rotation features exactly)
X_b <- matrix(0, nrow = nrow(X_test_z), ncol = length(rot_names))
colnames(X_b) <- rot_names
rownames(X_b) <- rownames(X_test_z)

# Fill time features one-by-one (handles lineage1_*, lineage2_*, etc.)
for (j in seq_along(time_cols)) {
  X_b[, time_cols[j]] <- X_test_z[, time_map[j]]
}

# Fill log2FC feature
X_b[, fc_name] <- meta_df$log2FC_feature

# Apply training preprocessing + project
X_bp <- apply_train_preproc(X_b, pc_b)
S_test_b <- X_bp %*% pc_b$rotation

test_df_b <- cbind(
  meta_df,
  PC1 = S_test_b[,1],
  PC2 = S_test_b[,2]
)

real_df_b <- if (exists("pc_df_b")) pc_df_b else data.frame(PC1 = pc_b$x[,1], PC2 = pc_b$x[,2])

# ---------- DIAGNOSTICS ----------
bad_b <- which(!is.finite(test_df_b$PC1) | !is.finite(test_df_b$PC2))
cat("Non-finite rows in Option B test projection:", length(bad_b), "\n")

test_df_b_finite <- test_df_b[is.finite(test_df_b$PC1) & is.finite(test_df_b$PC2), , drop = FALSE]

# Create a SINGLE numeric column to map to color (avoids aes list-combining)
if (COLOR_MODE == "signed") {
  test_df_b_finite$colval <- test_df_b_finite$log2FC_feature
  col_lab <- paste0(fc_name, " (signed)")
} else {
  test_df_b_finite$colval <- abs(test_df_b_finite$log2FC_feature)  # unsigned
  col_lab <- paste0("|", fc_name, "| (unsigned)")
}

# ---------- PLOTS ----------
p_orig <- ggplot() +
  geom_point(data = real_df_orig, aes(PC1, PC2), size = POINT_SIZE_BG, alpha = ALPHA_BG) +
  geom_point(data = test_df_orig, aes(PC1, PC2, color = mag, shape = family),
             size = POINT_SIZE_TEST, alpha = ALPHA_TEST) +
  theme_classic(base_size = 11) +
  labs(title = "Synthetic shapes projected into ORIGINAL PCA (z-curves only)",
       x = "PC1", y = "PC2", color = "magnitude (unsigned)", shape = "family")

p_b <- ggplot() +
  geom_point(data = real_df_b, aes(PC1, PC2), size = POINT_SIZE_BG, alpha = 0.10) +
  geom_point(data = test_df_b_finite,
             aes(PC1, PC2, shape = family, color = colval),
             size = POINT_SIZE_TEST, alpha = ALPHA_TEST) +
  theme_classic(base_size = 11) +
  labs(title = paste0("Synthetic shapes projected into OPTION B PCA (time features + ", fc_name, ")"),
       x = "PC1", y = "PC2", color = col_lab, shape = "family")

# Print side-by-side
if (requireNamespace("patchwork", quietly = TRUE)) {
  print(p_orig + p_b + patchwork::plot_layout(guides = "collect"))
} else {
  print(p_orig); print(p_b)
}

# ---------- NEW: visualize ALL synthetic trajectories (second figure) ----------
# Build long-form table for plotting (raw + z-scored)
traj_long_raw <- as.data.frame(X_test_raw) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(cols = -gene, names_to = "t_idx", values_to = "y") %>%
  mutate(t_idx = as.integer(gsub("[^0-9]+", "", t_idx)),
         t = tt[t_idx],
         space = "raw")

traj_long_z <- as.data.frame(X_test_z) %>%
  tibble::rownames_to_column("gene") %>%
  tidyr::pivot_longer(cols = -gene, names_to = "t_idx", values_to = "y") %>%
  mutate(t_idx = as.integer(gsub("[^0-9]+", "", t_idx)),
         t = tt[t_idx],
         space = "z")

traj_long <- bind_rows(traj_long_raw, traj_long_z) %>%
  left_join(meta_df %>% select(gene, family, mag, sign, extra, log2FC_feature), by = "gene") %>%
  mutate(family = factor(family, levels = c("flat","nearflat_noise","ramp_up","ramp_down","pulse_up","pulse_down")),
         space  = factor(space, levels = c("raw","z")))

# nicer label per synthetic curve
traj_long$label <- with(traj_long, paste0(gene, " | ", family, " | mag=", mag, " | sign=", sign))

options(repr.plot.width = 15.5, repr.plot.height = 8.8)

p_traj <- ggplot(traj_long, aes(x = t, y = y, group = gene)) +
  geom_line(aes(color = log2FC_feature), linewidth = 0.6, alpha = 0.85) +
  facet_grid(space ~ family, scales = "free_y") +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", color = NA),
    panel.spacing = unit(0.9, "lines")
  ) +
  labs(
    title = "All synthetic trajectories (raw and z-scored)",
    x = "time",
    y = "value",
    color = paste0(fc_name, " (signed)")
  )

print(p_traj)

### CELL 78
# ============================================================
# (3) Cluster median ± IQR ribbon panels (standardized space)
# ============================================================
df_long <- as.data.frame(X) %>%
  mutate(gene = genes, cluster = cl) %>%
  pivot_longer(cols = -c(gene, cluster), names_to = "tp", values_to = "z") %>%
  mutate(t_idx = as.integer(gsub(".*_", "", tp)),
         t = tt[t_idx])

summ_cl <- df_long %>%
  group_by(cluster, t) %>%
  summarise(
    q25 = quantile(z, 0.25, na.rm = TRUE),
    q50 = quantile(z, 0.50, na.rm = TRUE),
    q75 = quantile(z, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Plot
ncl <- nlevels(cl)
options(repr.plot.width = 9.5, repr.plot.height = max(4.5, 0.55 * ncl))
print(
  ggplot(summ_cl, aes(t, q50, color = cluster, fill = cluster)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.18, color = NA) +
    geom_line(linewidth = 1.05) +
    facet_wrap(~ cluster, ncol = 5, scales = "fixed") +
    theme_classic(base_size = 12) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(face = "bold")) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0, 1)) +
    labs(title = "3) Cluster trajectories: median ± IQR (standardized curves)",
         x = "Time (0→1)", y = "z-scored fitted expression")
)

### CELL 79
table(pc_df$cluster)

### CELL 80
# =========================
# Cluster composition vs WALD cutoff (stacked area)
#   Fixes "mostly white" by sweeping thresholds over the *observed* Wald range.
#   Convention: pass = waldStat_1 >= cutoff  (higher Wald = more significant)
# =========================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

stopifnot(exists("pc_df"), exists("asso_df"))
stopifnot(all(c("gene","cluster") %in% names(pc_df)))
stopifnot("waldStat_1" %in% names(asso_df))

# -------------------------
# CONFIG
# -------------------------
WALD_COL <- "waldStat_1"
N_T   <- 250                  # number of thresholds
CONDITIONAL_ON_PASS <- TRUE   # proportions among passing genes (recommended)
Q_LO <- 0.001                 # lower quantile for threshold grid
Q_HI <- 0.999                 # upper quantile for threshold grid

# Optional: add some "nice" cutoffs (Wald). These are safe; they just ensure inclusion.
ADD_CUTOFFS <- c(1, 3, 10, 30, 100)

# -------------------------
# Join pc_df + asso_df
# -------------------------
asso2 <- asso_df %>%
  tibble::rownames_to_column("gene") %>%
  select(gene, !!WALD_COL)

dat <- pc_df %>%
  select(gene, cluster, alpha) %>%
  left_join(asso2, by = "gene") %>%
  mutate(
    cluster = as.factor(cluster),
    wald = as.numeric(.data[[WALD_COL]]),
    wald = ifelse(is.finite(wald), wald, NA_real_)
  ) %>%
  filter(is.finite(wald))

if (nrow(dat) == 0) stop("No finite wald values after join. Check gene names + WALD_COL.")

clusters_all <- levels(dat$cluster)

# -------------------------
# Threshold grid from empirical Wald range (quantiles)
# Use a log-spaced grid to cover orders of magnitude if Wald is heavy-tailed.
# -------------------------
w_lo <- as.numeric(quantile(dat$wald, probs = Q_LO, na.rm = TRUE, names = FALSE))
w_hi <- as.numeric(quantile(dat$wald, probs = Q_HI, na.rm = TRUE, names = FALSE))

# safety: make strictly positive for log-grid; if Wald can be 0, bump a tiny epsilon
eps <- 1e-12
w_lo_pos <- max(w_lo, eps)
w_hi_pos <- max(w_hi, w_lo_pos * 1.01)

thr <- exp(seq(log(w_lo_pos), log(w_hi_pos), length.out = N_T))
thr <- sort(unique(c(thr, ADD_CUTOFFS)))
thr <- thr[thr > 0]

# -------------------------
# Sweep thresholds -> proportions by cluster
# pass = wald >= cutoff
# -------------------------
prop_df <- bind_rows(lapply(thr, function(cut) {
  sub <- dat %>% filter(wald >= cut)

  if (nrow(sub) == 0) {
    tibble(
      wald_cut = cut,
      n_pass   = 0L,
      cluster  = clusters_all,
      prop     = 0
    )
  } else {
    counts <- sub %>% count(cluster, name = "n")
    denom  <- if (CONDITIONAL_ON_PASS) sum(counts$n) else nrow(dat)

    counts %>%
      mutate(prop = n / denom) %>%
      right_join(tibble(cluster = clusters_all), by = "cluster") %>%
      mutate(
        n    = ifelse(is.na(n), 0L, n),
        prop = ifelse(is.na(prop), 0, prop),
        wald_cut = cut,
        n_pass   = nrow(sub)
      ) %>%
      select(wald_cut, n_pass, cluster, prop)
  }
})) %>%
  mutate(
    log10_wald_cut = log10(wald_cut),
    cluster = factor(cluster, levels = clusters_all)
  )

# boundaries (for the "flow" / contour feel)
bound <- prop_df %>%
  group_by(wald_cut) %>%
  arrange(cluster, .by_group = TRUE) %>%
  mutate(cum = cumsum(prop)) %>%
  ungroup()

# annotate n_pass at a few cutoffs (avoid clutter)
ann <- prop_df %>%
  group_by(wald_cut) %>%
  summarise(n_pass = max(n_pass), .groups = "drop") %>%
  arrange(wald_cut) %>%
  slice(round(seq(1, n(), length.out = 6))) %>%
  mutate(log10_wald_cut = log10(wald_cut))

# -------------------------
# Plot
# -------------------------
options(repr.plot.width = 8.2, repr.plot.height = 4.9)

p <- ggplot(prop_df, aes(x = log10_wald_cut, y = prop, fill = cluster)) +
  geom_area(alpha = 0.95, color = NA) +
  geom_line(
    data = bound,
    aes(x = log10_wald_cut, y = cum, group = cluster),
    inherit.aes = FALSE,
    linewidth = 0.25,
    alpha = 0.35
  ) +
  geom_text(
    data = ann,
    aes(x = log10_wald_cut, y = 1.02, label = paste0("n=", n_pass)),
    inherit.aes = FALSE,
    size = 3,
    vjust = 0
  ) +
  theme_classic(base_size = 12) +
  coord_cartesian(ylim = c(0, 1.08), clip = "off") +
  labs(
    title = "Cluster composition vs Wald cutoff",
    subtitle = paste0(
      "Pass = Wald >= cutoff. Cutoffs span observed Wald range (q", Q_LO, "–q", Q_HI, ").  ",
      "y = proportion (", ifelse(CONDITIONAL_ON_PASS, "among passing genes", "relative to all genes"), "); ",
      "x = log10(Wald cutoff)."
    ),
    x = "log10(Wald cutoff)",
    y = "Proportion",
    fill = "Cluster"
  )

print(p)

### CELL 81
# =========================
# PCA colored by log10(Wald) (cap=2) + label top-10 Wald genes
# (FIXED: geom_text now has x/y mapped)
# =========================
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

stopifnot(exists("pc3"), exists("asso_df"))
stopifnot(all(c("gene","PC1","PC2","alpha") %in% names(pc3)))
stopifnot("waldStat_1" %in% names(asso_df))

# ---- CONFIG ----
WALD_LOG10_CAP <- 2
PSEUDOCOUNT    <- 1e-9
POINT_SIZE     <- 1
LABEL_TOP_N    <- 30
LABEL_SIZE     <- 2.6
LABEL_ALPHA    <- 0.95
LABEL_NUDGE_Y  <- 0.01

# ---- attach Wald ----
wdf <- asso_df %>%
  rownames_to_column("gene") %>%
  select(gene, waldStat_1)

pc_wald <- pc_df %>%
  left_join(wdf, by = "gene") %>%
  mutate(
    wald_log10     = log10(pmax(waldStat_1, PSEUDOCOUNT)),
    wald_log10_cap = pmin(wald_log10, WALD_LOG10_CAP)
  )

# ---- pick top genes (by Wald) that exist in the embedding ----
top_genes <- wdf %>%
  filter(gene %in% pc_wald$gene, is.finite(waldStat_1)) %>%
  arrange(desc(waldStat_1)) %>%
  slice_head(n = LABEL_TOP_N) %>%
  pull(gene)

lab_df <- pc_wald %>% filter(gene %in% top_genes)

options(repr.plot.width = 5.8, repr.plot.height = 5)
print(
  ggplot(pc_wald, aes(PC1, PC2, color = wald_log10_cap, alpha = alpha)) +
    geom_point(size = POINT_SIZE) +
    scale_alpha_identity() +
    scale_color_continuous(
      limits = c(min(pc_wald$wald_log10_cap, na.rm = TRUE), WALD_LOG10_CAP)
    ) +
    geom_text(
      data = lab_df,
      aes(x = PC1, y = PC2, label = gene),  # <- FIX: provide x/y explicitly
      inherit.aes = FALSE,
      size = LABEL_SIZE,
      alpha = LABEL_ALPHA,
      nudge_y = LABEL_NUDGE_Y
    ) +
    theme_classic(base_size = 12) +
    labs(
      title = "PCA on standardized fitted curves",
      subtitle = paste0(
        "Color = log10(Wald) capped at ", WALD_LOG10_CAP,
        "; labeled top ", LABEL_TOP_N, " Wald genes"
      ),
      x = "PC1", y = "PC2",
      color = paste0("log10(Wald)\n(capped @ ", WALD_LOG10_CAP, ")")
    )
)

### CELL 82
options(repr.plot.width = 2.75, repr.plot.height = 2.5)

plot_tradeSeq_gene_pub(
  "ZFP36", gamList2, sce2, counts_sub, PT_mat, sf,
  cellWeights = CWz_ret,
  lineages = c(1,2),
  colours = c("#386cb0", "#DC143C"),
  time_max = 110,
  x_tick_step = 20,        # <-- every 20 min
  color_points_by_lineage = TRUE,
  points_size = 0.75,
  points_alpha = 0.3,
  axis_text_size = 14,
  axis_title_size = 16,
  # ylim = c(0, 0.1),
  lineage_labels = c("Inflammatory-hi", "Erythroid-like"),
  show_legend = FALSE,
  save_pdf = FALSE
)

### CELL 84
# =========================
# Dictionary/library-of-functions (UPDATED + post-hoc impulse gating)
#   Models fit per gene:
#     - linear (signed slope)
#     - logistic (signed amplitude; covers up/down)
#     - impulse (signed amplitude; bump/dip)
#
#   Selection:
#     1) Find best non-impulse model (min BIC among linear/logistic)
#     2) Allow "impulse" ONLY IF BOTH:
#          (a) ΔBIC = BIC_nonimp - BIC_imp >= BIC_GATE  (default 10)
#          (b) rel_SSE_gain = (SSE_nonimp - SSE_imp)/SSE_nonimp >= REL_IMPROVE_GATE
#              (default 0.20 = 20%, i.e. stricter than your ~17% example)
#     3) Otherwise: pick best non-impulse
#
#   Output:
#     - arch_df with archetype_final + archetype_raw (raw min-BIC winner)
#     - keeps BIC/SSE for all models so you can audit decisions
# =========================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

stopifnot(exists("traj_mat"), exists("tt"), exists("cluster_C"))
stopifnot(nrow(traj_mat) == length(cluster_C))

genes <- rownames(traj_mat)
cl <- factor(cluster_C)

# -------------------------
# CONFIG
# -------------------------
set.seed(1)
max_genes_fit <- Inf
use_top_wald  <- FALSE          # set TRUE if you have asso_df$waldStat_1
bic_noise_floor <- 1e-6
n_starts <- 6                   # restarts for nonlinear models

# Post-hoc gating for calling IMPULSE
BIC_GATE          <- 10         # require impulse to beat best non-impulse by at least this much
REL_IMPROVE_GATE  <- 0.20       # require >=20% SSE reduction vs best non-impulse (stricter than 17%)
SSE_EPS           <- 1e-12      # numerical guard

# -------------------------
# Optional: pick subset of genes to fit
# -------------------------
genes_fit <- genes
if (use_top_wald && exists("asso_df") && is.finite(max_genes_fit)) {
  ord <- order(asso_df$waldStat_1, decreasing = TRUE, na.last = TRUE)
  genes_fit <- intersect(rownames(asso_df)[ord], genes)
  genes_fit <- head(genes_fit, max_genes_fit)
} else if (is.finite(max_genes_fit)) {
  genes_fit <- head(genes_fit, max_genes_fit)
}
cat("Fitting models for", length(genes_fit), "genes\n")

t <- tt
n <- length(t)

# -------------------------
# Archetypes (SIGNED forms)
# -------------------------
sigmoid <- function(z) 1/(1+exp(-z))

# (A) Signed logistic (covers monotone up/down)
# y(t) = b + A * sigmoid(k*(t - t0)),  k>0, t0 in (0,1), A is signed
f_logistic_signed <- function(par, t) {
  b  <- par[1]
  A  <- par[2]        # signed
  k  <- exp(par[3])   # >0
  t0 <- sigmoid(par[4])
  b + A * sigmoid(k * (t - t0))
}

# (B) Linear trend (monotone-but-not-sigmoid; can be up or down)
# y(t) = b + m*(t - 0.5)
f_linear <- function(par, t) {
  b <- par[1]
  m <- par[2]  # signed slope
  b + m * (t - 0.5)
}

# (C) Signed impulse (difference of sigmoids; bump or dip), t2>t1
f_impulse_signed <- function(par, t) {
  b  <- par[1]
  A  <- par[2]        # signed
  k1 <- exp(par[3])
  k2 <- exp(par[4])
  t1 <- sigmoid(par[5])
  dt <- sigmoid(par[6]) * (1 - t1)  # dt in (0, 1-t1)
  t2 <- t1 + dt
  b + A * (sigmoid(k1 * (t - t1)) - sigmoid(k2 * (t - t2)))
}

# -------------------------
# Fit helper (SSE) + BIC
# -------------------------
fit_model <- function(y, t, f, p, par_init_list, maxit = 2000) {
  obj <- function(par) {
    mu <- f(par, t)
    sum((y - mu)^2, na.rm = TRUE)
  }
  best <- list(par = rep(NA_real_, p), sse = Inf, conv = 1)
  for (par0 in par_init_list) {
    out <- tryCatch(
      optim(par0, obj, method = "BFGS", control = list(maxit = maxit)),
      error = function(e) NULL
    )
    if (!is.null(out) && is.finite(out$value) && out$value < best$sse) {
      best$par <- out$par
      best$sse <- out$value
      best$conv <- out$convergence
    }
  }
  best
}

bic_from_sse <- function(sse, n, p) {
  sigma2 <- max(sse / n, bic_noise_floor)
  n * log(sigma2) + p * log(n)
}

# Initializations
make_inits_logistic <- function(y, t, n = 6) {
  y0 <- y[1]; y1 <- y[length(y)]
  amp <- (y1 - y0)
  sd_y <- sd(y, na.rm = TRUE) + 1e-3
  out <- vector("list", n)
  for (i in seq_len(n)) {
    b  <- y0 + rnorm(1, 0, 0.2 * sd_y)
    A  <- amp + rnorm(1, 0, 0.5 * sd_y)  # signed
    k  <- log(runif(1, 2, 60))
    t0 <- qlogis(runif(1, 0.1, 0.9))
    out[[i]] <- c(b, A, k, t0)
  }
  out
}

make_inits_impulse <- function(y, t, n = 6) {
  y0 <- y[1]
  yr <- range(y, finite = TRUE)
  amp <- max(yr[2] - yr[1], 1e-3)
  sd_y <- sd(y, na.rm = TRUE) + 1e-3
  out <- vector("list", n)
  for (i in seq_len(n)) {
    b  <- y0 + rnorm(1, 0, 0.2 * sd_y)
    A  <- (runif(1, -1, 1) * amp) + rnorm(1, 0, 0.2 * sd_y)  # signed
    k1 <- log(runif(1, 2, 80))
    k2 <- log(runif(1, 2, 80))
    t1 <- qlogis(runif(1, 0.05, 0.75))
    dt <- qlogis(runif(1, 0.05, 0.8))
    out[[i]] <- c(b, A, k1, k2, t1, dt)
  }
  out
}

# Linear: closed-form LS (fast + stable)
fit_linear_ls <- function(y, t) {
  X <- cbind(1, (t - 0.5))
  ok <- is.finite(y) & is.finite(X[,1]) & is.finite(X[,2])
  if (sum(ok) < 3) return(list(par = c(NA_real_, NA_real_), sse = Inf, conv = 1))
  beta <- solve(crossprod(X[ok, , drop = FALSE]), crossprod(X[ok, , drop = FALSE], y[ok]))
  yhat <- as.numeric(X %*% beta)
  sse <- sum((y - yhat)^2, na.rm = TRUE)
  list(par = as.numeric(beta), sse = sse, conv = 0)
}

# Decode interpretable params
decode_logistic <- function(par) {
  list(b = par[1], A = par[2], k = exp(par[3]), t0 = sigmoid(par[4]),
       direction = ifelse(par[2] >= 0, "up", "down"))
}
decode_linear <- function(par) {
  list(b = par[1], m = par[2], direction = ifelse(par[2] >= 0, "up", "down"))
}
decode_impulse <- function(par) {
  b  <- par[1]; A <- par[2]
  k1 <- exp(par[3]); k2 <- exp(par[4])
  t1 <- sigmoid(par[5]); dt <- sigmoid(par[6])*(1 - t1); t2 <- t1 + dt
  list(b=b, A=A, k1=k1, k2=k2, t1=t1, t2=t2, dt=(t2-t1),
       polarity = ifelse(A >= 0, "bump", "dip"))
}

# -------------------------
# Main loop (fit all 3, then post-hoc select)
# -------------------------
res <- vector("list", length(genes_fit))
names(res) <- genes_fit
pb <- txtProgressBar(min = 0, max = length(genes_fit), style = 3)

for (ii in seq_along(genes_fit)) {
  g <- genes_fit[ii]
  y <- as.numeric(traj_mat[g, ])

  # Fit all candidates
  fit_lin <- fit_linear_ls(y, t)
  fit_log <- fit_model(y, t, f_logistic_signed, p = 4,
                       par_init_list = make_inits_logistic(y, t, n_starts))
  fit_imp <- fit_model(y, t, f_impulse_signed, p = 6,
                       par_init_list = make_inits_impulse(y, t, n_starts))

  # Score
  bic_lin <- bic_from_sse(fit_lin$sse, n, 2)
  bic_log <- bic_from_sse(fit_log$sse, n, 4)
  bic_imp <- bic_from_sse(fit_imp$sse, n, 6)

  sse_lin <- fit_lin$sse
  sse_log <- fit_log$sse
  sse_imp <- fit_imp$sse

  # Raw min-BIC winner (for auditing)
  bics_all <- c(linear = bic_lin, logistic = bic_log, impulse = bic_imp)
  best_raw <- names(which.min(bics_all))

  # Best non-impulse
  bics_non <- c(linear = bic_lin, logistic = bic_log)
  sse_non  <- c(linear = sse_lin, logistic = sse_log)
  best_non <- names(which.min(bics_non))

  bic_best_non <- unname(bics_non[best_non])
  sse_best_non <- unname(sse_non[best_non])

  # Gating metrics for impulse
  dBIC_non_minus_imp <- bic_best_non - bic_imp
  rel_sse_gain <- (sse_best_non - sse_imp) / max(sse_best_non, SSE_EPS)

  impulse_ok <- is.finite(dBIC_non_minus_imp) && is.finite(rel_sse_gain) &&
    (dBIC_non_minus_imp >= BIC_GATE) && (rel_sse_gain >= REL_IMPROVE_GATE)

  best_final <- if (impulse_ok) "impulse" else best_non

  # Decode params for FINAL winner
  pars_best <- switch(best_final,
    linear   = decode_linear(fit_lin$par),
    logistic = decode_logistic(fit_log$par),
    impulse  = decode_impulse(fit_imp$par)
  )

  res[[g]] <- list(
    gene = g,
    cluster = as.character(cl[match(g, genes)]),

    archetype_raw   = best_raw,
    archetype_final = best_final,

    bic_linear = bic_lin, bic_logistic = bic_log, bic_impulse = bic_imp,
    sse_linear = sse_lin, sse_logistic = sse_log, sse_impulse = sse_imp,

    best_nonimp = best_non,
    dBIC_nonimp_minus_imp = dBIC_non_minus_imp,
    rel_sse_gain_imp_vs_nonimp = rel_sse_gain,
    impulse_pass = impulse_ok,

    par_best = pars_best
  )

  setTxtProgressBar(pb, ii)
}
close(pb)

# -------------------------
# Pack results into a data.frame
# -------------------------
arch_df <- bind_rows(lapply(res, function(x) {
  p <- x$par_best
  out <- data.frame(
    gene = x$gene,
    cluster = x$cluster,

    archetype_raw   = x$archetype_raw,
    archetype_final = x$archetype_final,

    best_nonimp = x$best_nonimp,
    dBIC_nonimp_minus_imp = x$dBIC_nonimp_minus_imp,
    rel_sse_gain_imp_vs_nonimp = x$rel_sse_gain_imp_vs_nonimp,
    impulse_pass = x$impulse_pass,

    bic_linear = x$bic_linear,
    bic_logistic = x$bic_logistic,
    bic_impulse = x$bic_impulse,

    sse_linear = x$sse_linear,
    sse_logistic = x$sse_logistic,
    sse_impulse = x$sse_impulse,

    stringsAsFactors = FALSE
  )

  # shared-ish fields
  out$b <- if (!is.null(p$b)) p$b else NA_real_
  out$A <- if (!is.null(p$A)) p$A else NA_real_
  out$direction <- if (!is.null(p$direction)) p$direction else NA_character_

  # logistic fields
  out$k  <- if (!is.null(p$k))  p$k  else NA_real_
  out$t0 <- if (!is.null(p$t0)) p$t0 else NA_real_

  # linear field
  out$m  <- if (!is.null(p$m))  p$m  else NA_real_

  # impulse fields
  out$k1 <- if (!is.null(p$k1)) p$k1 else NA_real_
  out$k2 <- if (!is.null(p$k2)) p$k2 else NA_real_
  out$t1 <- if (!is.null(p$t1)) p$t1 else NA_real_
  out$t2 <- if (!is.null(p$t2)) p$t2 else NA_real_
  out$dt <- if (!is.null(p$dt)) p$dt else NA_real_
  out$polarity <- if (!is.null(p$polarity)) p$polarity else NA_character_
  out
}))

arch_df$cluster <- factor(arch_df$cluster)
arch_df$archetype_raw   <- factor(arch_df$archetype_raw,   levels = c("linear","logistic","impulse"))
arch_df$archetype_final <- factor(arch_df$archetype_final, levels = c("linear","logistic","impulse"))

cat("\nArchetype counts (raw min-BIC):\n")
print(table(arch_df$archetype_raw))

cat("\nArchetype counts (final after impulse gating):\n")
print(table(arch_df$archetype_final))

cat("\nImpulse gate settings:\n")
cat("  BIC_GATE =", BIC_GATE, "\n")
cat("  REL_IMPROVE_GATE =", REL_IMPROVE_GATE, " (fractional SSE reduction)\n")

cat("\nSaved archetype assignments to: arch_df\n")

### CELL 85
length(tt)

### CELL 86
# ---- Composite figure: independent per-gene plots, stitched into rows by category ----
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

stopifnot(exists("arch_df"), exists("traj_mat"), exists("tt"))

# patchwork is the cleanest way to stitch ggplots
if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Please install patchwork: install.packages('patchwork')")
}

sigmoid <- function(z) 1/(1+exp(-z))

# ----------------- CONFIG -----------------
N_PER_CAT <- 10
SEED      <- 3
CATS      <- c("linear_up","linear_down","logistic_up","logistic_down","impulse")

PLOT_DELTA_FOR_IMPULSE <- TRUE  # plot y - y(t0) for impulse panels

SMOOTH_ALPHA <- 0.85
SMOOTH_LWD   <- 0.80
FIT_ALPHA    <- 1.00
FIT_LWD      <- 1.05
FIT_COLOR    <- "#D62728"

set.seed(SEED)

# ----------------- archetype column detection -----------------
detect_type_col <- function(df) {
  cand <- c("archetype_final","archetype_raw","archetype","best","model","type","arch")
  hit <- cand[cand %in% names(df)]
  if (!length(hit)) stop("Couldn't find an archetype column in arch_df. Tried: ", paste(cand, collapse=", "))
  hit[1]
}
TYPE_COL <- detect_type_col(arch_df)

get_type_scalar <- function(row) as.character(row[[TYPE_COL]][1])

infer_direction <- function(type, row) {
  if ("direction" %in% names(row)) {
    d <- as.character(row$direction[1])
    if (isTRUE(nzchar(d))) return(d)
  }
  if (type == "linear" && "m" %in% names(row) && is.finite(row$m[1])) {
    return(ifelse(row$m[1] >= 0, "up", "down"))
  }
  if (type == "logistic" && "A" %in% names(row) && is.finite(row$A[1])) {
    return(ifelse(row$A[1] >= 0, "up", "down"))
  }
  NA_character_
}

cat_from_type_dir <- function(type, direction) {
  if (type == "linear"   && direction == "up")   return("linear_up")
  if (type == "linear"   && direction == "down") return("linear_down")
  if (type == "logistic" && direction == "up")   return("logistic_up")
  if (type == "logistic" && direction == "down") return("logistic_down")
  if (type == "impulse")                         return("impulse")
  NA_character_
}

# ----------------- reconstruct best-fit curve -----------------
reconstruct_yhat <- function(row, t) {
  type <- get_type_scalar(row)
  b <- row$b[1]

  if (type == "linear") {
    if (!("m" %in% names(row))) return(rep(NA_real_, length(t)))
    m <- row$m[1]
    if (!is.finite(b) || !is.finite(m)) return(rep(NA_real_, length(t)))
    return(b + m * (t - 0.5))
  }

  if (type == "logistic") {
    need <- c("A","k","t0")
    if (!all(need %in% names(row))) return(rep(NA_real_, length(t)))
    A <- row$A[1]; k <- row$k[1]; t0 <- row$t0[1]
    if (!is.finite(b) || !is.finite(A) || !is.finite(k) || !is.finite(t0)) return(rep(NA_real_, length(t)))
    return(b + A * sigmoid(k * (t - t0)))
  }

  if (type == "impulse") {
    need <- c("A","k1","k2","t1","t2")
    if (!all(need %in% names(row))) return(rep(NA_real_, length(t)))
    A <- row$A[1]
    k1 <- row$k1[1]; k2 <- row$k2[1]; t1 <- row$t1[1]; t2 <- row$t2[1]
    if (!is.finite(b) || !is.finite(A) || !is.finite(k1) || !is.finite(k2) || !is.finite(t1) || !is.finite(t2)) {
      return(rep(NA_real_, length(t)))
    }
    return(b + A * (sigmoid(k1 * (t - t1)) - sigmoid(k2 * (t - t2))))
  }

  rep(NA_real_, length(t))
}

# ----------------- build sampling table -----------------
pick0 <- arch_df %>%
  mutate(type = as.character(.data[[TYPE_COL]]))

pick0$direction2 <- vapply(seq_len(nrow(pick0)), function(i) {
  infer_direction(pick0$type[i], pick0[i, , drop = FALSE])
}, character(1))

pick <- pick0 %>%
  mutate(cat = mapply(cat_from_type_dir, type, direction2)) %>%
  filter(cat %in% CATS, gene %in% rownames(traj_mat)) %>%
  group_by(cat) %>%
  slice_sample(prop = 1) %>%
  slice_head(n = N_PER_CAT) %>%
  ungroup()

cat("Using archetype column: ", TYPE_COL, "\n", sep="")
cat("Sampled genes per category:\n")
print(pick %>% count(cat))

# ----------------- per-gene ggplot -----------------
gene_panel_plot <- function(g, row, cat_i) {
  t <- tt
  y <- as.numeric(traj_mat[g, ])
  yhat <- reconstruct_yhat(row, t)

  if (isTRUE(PLOT_DELTA_FOR_IMPULSE) && cat_i == "impulse") {
    y0 <- y[1]
    y <- y - y0
    yhat <- yhat - y0
  }

  df <- tibble(t = t, y = y, yhat = yhat)

  ggplot(df, aes(t, y)) +
    geom_line(linewidth = SMOOTH_LWD, alpha = SMOOTH_ALPHA) +
    geom_line(aes(y = yhat),
              linewidth = FIT_LWD,
              linetype  = 2,
              alpha     = FIT_ALPHA,
              color     = FIT_COLOR,
              na.rm     = TRUE) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(size = 9, hjust = 0.5, margin = margin(b = 2)),
      axis.title = element_blank(),
      axis.text  = element_text(size = 7),
      plot.margin = margin(2, 2, 2, 2)
    ) +
    ggtitle(g)
}

# ----------------- stitch plots into rows -----------------
row_plots <- lapply(CATS, function(cat_i) {
  rows_cat <- pick %>% filter(cat == cat_i)
  if (nrow(rows_cat) == 0) return(NULL)

  plist <- lapply(seq_len(nrow(rows_cat)), function(i) {
    g <- rows_cat$gene[i]
    row <- rows_cat[i, , drop = FALSE]
    gene_panel_plot(g, row, cat_i)
  })

  patchwork::wrap_plots(plist, nrow = 1) +
    patchwork::plot_annotation(
      title = cat_i,
      theme = theme(plot.title = element_text(size = 12, face = "bold", hjust = 0))
    )
})

row_plots <- Filter(Negate(is.null), row_plots)

final_plot <- Reduce(`/`, row_plots) +
  patchwork::plot_annotation(
    title = paste0("Random genes per category (N=", N_PER_CAT, " per row; dashed=best-fit; seed=", SEED, ")"),
    subtitle = paste0(
      "Composite = independent per-gene plots (true per-panel scaling). archetype col = ", TYPE_COL,
      if (PLOT_DELTA_FOR_IMPULSE) " | impulse plotted as Δ from baseline" else ""
    )
  )

print(final_plot)

### CELL 87
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

stopifnot(exists("arch_df"), exists("traj_mat"), exists("tt"))

# ---------------------------
# Shared helpers
# ---------------------------
sigmoid <- function(z) 1/(1+exp(-z))

mask_finite <- function(y, yhat) is.finite(y) & is.finite(yhat)

sse_masked <- function(y, yhat) {
  ok <- mask_finite(y, yhat)
  if (sum(ok) < 3) return(NA_real_)
  sum((y[ok] - yhat[ok])^2)
}

tss_masked <- function(y) {
  ok <- is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  sum((y[ok] - mean(y[ok]))^2)
}

bic_from_sse <- function(sse_val, n, p, noise_floor = 1e-6) {
  if (!is.finite(sse_val)) return(NA_real_)
  sigma2 <- max(sse_val / n, noise_floor)
  n * log(sigma2) + p * log(n)
}

# ---------------------------
# Model definitions (MUST match your fitting code)
# ---------------------------
pred_linear <- function(par, t) {
  b <- par[1]; m <- par[2]
  b + m * (t - 0.5)
}

pred_logistic_signed <- function(par, t) {
  b  <- par[1]
  A  <- par[2]        # signed
  k  <- exp(par[3])   # >0
  t0 <- sigmoid(par[4])
  b + A * sigmoid(k * (t - t0))
}

pred_impulse_signed <- function(par, t) {
  b  <- par[1]
  A  <- par[2]        # signed
  k1 <- exp(par[3])
  k2 <- exp(par[4])
  t1 <- sigmoid(par[5])
  dt <- sigmoid(par[6]) * (1 - t1)
  t2 <- t1 + dt
  b + A * (sigmoid(k1 * (t - t1)) - sigmoid(k2 * (t - t2)))
}

# ---------------------------
# Fitting routines (refit per gene for fair comparison)
# ---------------------------

# Linear: closed-form LS (same as your code)
fit_linear_ls <- function(y, t) {
  X <- cbind(1, (t - 0.5))
  ok <- is.finite(y) & is.finite(X[,1]) & is.finite(X[,2])
  if (sum(ok) < 3) return(list(par = c(NA_real_, NA_real_), sse = NA_real_, conv = 1))
  beta <- solve(crossprod(X[ok, , drop = FALSE]), crossprod(X[ok, , drop = FALSE], y[ok]))
  yhat <- as.numeric(X %*% beta)
  list(par = as.numeric(beta), sse = sse_masked(y, yhat), conv = 0)
}

# Nonlinear multi-start optim (same objective style as your code)
fit_model <- function(y, t, pred_fun, p, par_init_list, maxit = 2000) {
  obj <- function(par) {
    mu <- pred_fun(par, t)
    sse_masked(y, mu)
  }
  best <- list(par = rep(NA_real_, p), sse = Inf, conv = 1)
  for (par0 in par_init_list) {
    out <- tryCatch(
      optim(par0, obj, method = "BFGS", control = list(maxit = maxit)),
      error = function(e) NULL
    )
    if (!is.null(out) && is.finite(out$value) && out$value < best$sse) {
      best$par  <- out$par
      best$sse  <- out$value
      best$conv <- out$convergence
    }
  }
  if (!is.finite(best$sse) || is.infinite(best$sse)) best$sse <- NA_real_
  best
}

# Initializations (match your updated code closely)
make_inits_logistic <- function(y, t, n = 6) {
  y0 <- y[1]; y1 <- y[length(y)]
  amp <- (y1 - y0)
  sd_y <- sd(y, na.rm = TRUE) + 1e-3
  out <- vector("list", n)
  for (i in seq_len(n)) {
    b  <- y0 + rnorm(1, 0, 0.2 * sd_y)
    A  <- amp + rnorm(1, 0, 0.5 * sd_y)  # signed
    k  <- log(runif(1, 2, 60))
    t0 <- qlogis(runif(1, 0.1, 0.9))
    out[[i]] <- c(b, A, k, t0)
  }
  out
}

make_inits_impulse <- function(y, t, n = 6) {
  y0 <- y[1]
  yr <- range(y, finite = TRUE)
  amp <- max(yr[2] - yr[1], 1e-3)
  sd_y <- sd(y, na.rm = TRUE) + 1e-3
  out <- vector("list", n)
  for (i in seq_len(n)) {
    b  <- y0 + rnorm(1, 0, 0.2 * sd_y)
    A  <- (runif(1, -1, 1) * amp) + rnorm(1, 0, 0.2 * sd_y)  # signed
    k1 <- log(runif(1, 2, 80))
    k2 <- log(runif(1, 2, 80))
    t1 <- qlogis(runif(1, 0.05, 0.75))
    dt <- qlogis(runif(1, 0.05, 0.8))
    out[[i]] <- c(b, A, k1, k2, t1, dt)
  }
  out
}

# Decode interpretable params (for display)
decode_logistic <- function(par) {
  list(b = par[1], A = par[2], k = exp(par[3]), t0 = sigmoid(par[4]),
       direction = ifelse(par[2] >= 0, "up", "down"))
}
decode_linear <- function(par) {
  list(b = par[1], m = par[2], direction = ifelse(par[2] >= 0, "up", "down"))
}
decode_impulse <- function(par) {
  b  <- par[1]; A <- par[2]
  k1 <- exp(par[3]); k2 <- exp(par[4])
  t1 <- sigmoid(par[5]); dt <- sigmoid(par[6])*(1 - t1); t2 <- t1 + dt
  list(b=b, A=A, k1=k1, k2=k2, t1=t1, t2=t2, dt=(t2-t1),
       polarity = ifelse(A >= 0, "bump", "dip"))
}

# ---------------------------
# Core: refit all candidate models for a queried gene
# ---------------------------
fit_candidates_for_gene <- function(gene, n_starts = 6, noise_floor = 1e-6) {
  if (!gene %in% rownames(traj_mat)) stop("Gene not in traj_mat: ", gene)

  y <- as.numeric(traj_mat[gene, ])
  t <- tt
  n <- length(t)
  tss <- tss_masked(y)

  # 1) linear
  fit_lin <- fit_linear_ls(y, t)
  bic_lin <- bic_from_sse(fit_lin$sse, n, p = 2, noise_floor = noise_floor)

  # 2) logistic signed
  fit_log <- fit_model(y, t, pred_logistic_signed, p = 4,
                       par_init_list = make_inits_logistic(y, t, n_starts))
  bic_log <- bic_from_sse(fit_log$sse, n, p = 4, noise_floor = noise_floor)

  # 3) impulse signed
  fit_imp <- fit_model(y, t, pred_impulse_signed, p = 6,
                       par_init_list = make_inits_impulse(y, t, n_starts))
  bic_imp <- bic_from_sse(fit_imp$sse, n, p = 6, noise_floor = noise_floor)

  # Build predictions (for plotting)
  yhat_lin <- if (all(is.finite(fit_lin$par))) pred_linear(fit_lin$par, t) else rep(NA_real_, n)
  yhat_log <- if (all(is.finite(fit_log$par))) pred_logistic_signed(fit_log$par, t) else rep(NA_real_, n)
  yhat_imp <- if (all(is.finite(fit_imp$par))) pred_impulse_signed(fit_imp$par, t) else rep(NA_real_, n)

  tab <- tibble(
    model = c("linear", "logistic", "impulse"),
    p     = c(2, 4, 6),
    sse   = c(fit_lin$sse, fit_log$sse, fit_imp$sse),
    bic   = c(bic_lin, bic_log, bic_imp),
    r2    = ifelse(is.finite(tss) & tss > 0, 1 - sse / tss, NA_real_),
    conv  = c(fit_lin$conv, fit_log$conv, fit_imp$conv)
  ) %>% arrange(bic)

  decoded <- list(
    linear  = if (all(is.finite(fit_lin$par))) decode_linear(fit_lin$par) else NULL,
    logistic= if (all(is.finite(fit_log$par))) decode_logistic(fit_log$par) else NULL,
    impulse = if (all(is.finite(fit_imp$par))) decode_impulse(fit_imp$par) else NULL
  )

  list(
    gene = gene,
    y = y, t = t,
    tab = tab,
    yhat = list(linear=yhat_lin, logistic=yhat_log, impulse=yhat_imp),
    decoded = decoded
  )
}

# ---------------------------
# User-facing function
# ---------------------------
compare_gene_models <- function(gene, show_plot = TRUE, n_starts = 6, noise_floor = 1e-6) {
  out <- fit_candidates_for_gene(gene, n_starts = n_starts, noise_floor = noise_floor)
  tab <- out$tab

  cat("\nGene:", gene, "\n")
  cat("Best by BIC:", tab$model[1], "\n\n")
  print(tab)

  best_non_imp <- tab %>% filter(model != "impulse") %>% slice(1)
  best_imp     <- tab %>% filter(model == "impulse") %>% slice(1)

  if (nrow(best_imp) == 1 && nrow(best_non_imp) == 1 &&
      is.finite(best_imp$bic) && is.finite(best_non_imp$bic)) {
    cat("\nImpulse vs best non-impulse:\n")
    cat("  ΔSSE  (nonimp - imp): ", best_non_imp$sse - best_imp$sse, "\n", sep="")
    cat("  ΔBIC  (nonimp - imp): ", best_non_imp$bic - best_imp$bic, "\n", sep="")
    cat("  ΔR²   (imp - nonimp): ", best_imp$r2 - best_non_imp$r2, "\n", sep="")
    cat("  Note: ΔBIC > 10 is often considered strong evidence for the lower-BIC model.\n")
  }

  # Show decoded params for the winning model
  win <- tab$model[1]
  if (!is.null(out$decoded[[win]])) {
    cat("\nDecoded params (winner = ", win, "):\n", sep = "")
    print(out$decoded[[win]])
  }

  if (show_plot) {
    top2 <- tab$model[1:min(2, nrow(tab))]
    dfp <- bind_rows(lapply(top2, function(mn) {
      data.frame(t = out$t, yhat = out$yhat[[mn]], model = mn)
    })) %>% filter(is.finite(yhat))

    dfy <- data.frame(t = out$t, y = out$y) %>% filter(is.finite(y))

    print(
      ggplot() +
        geom_line(data = dfy, aes(t, y), linewidth = 0.9, alpha = 0.55) +
        geom_line(data = dfp, aes(t, yhat, color = model, linetype = model), linewidth = 1.0) +
        theme_classic(base_size = 12) +
        labs(
          title = paste0("Model comparison for ", gene),
          subtitle = paste0("Top by BIC: ", paste(top2, collapse = ", "),
                            " (refit per model; conv=", paste(tab$conv[match(top2, tab$model)], collapse = ","), ")"),
          x = "Time (0→1)", y = "Expression (per-gene scale)"
        )
    )
  }

  invisible(out)
}


### CELL 88
compare_gene_models("BIRC3")

### CELL 89
# =========================
# Re-plot MDS three ways:
#   (i) color = original clusters
#   (ii) color = t_peak (computed on standardized curve X)
#   (iii) color = best-fit function (with up/down split)
# =========================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

stopifnot(exists("mds_df"), exists("pc_df"))   # created in your big cell above
stopifnot(exists("arch_df"))                  # your model-fit results
stopifnot(all(c("gene","cluster","MDS1","MDS2","alpha") %in% names(mds_df)))
stopifnot(all(c("gene","t_peak") %in% names(pc_df)))

# ---- attach t_peak to mds_df (from pc_df) ----
mds2 <- mds_df %>%
  left_join(pc_df %>% select(gene, t_peak), by = "gene")

# ---- attach archetype label (best-fit function + direction) ----
# This is robust to different column names you’ve used across versions.
get_arch_type <- function(df) {
  if ("archetype_final" %in% names(df)) return(as.character(df$archetype_final))
  if ("archetype_raw"   %in% names(df)) return(as.character(df$archetype_raw))
  if ("archetype"       %in% names(df)) return(as.character(df$archetype))
  stop("arch_df must contain one of: archetype_final, archetype_raw, archetype")
}

arch_type <- get_arch_type(arch_df)

arch_map <- arch_df %>%
  mutate(
    archetype_any = arch_type,
    direction_any = if ("direction" %in% names(.)) as.character(direction) else NA_character_,
    polarity_any  = if ("polarity"  %in% names(.)) as.character(polarity)  else NA_character_,
    # Build a display label that splits monotone models by up/down.
    # For impulse, we optionally label bump/dip if polarity exists.
    archetype_label = case_when(
      archetype_any == "linear"   & direction_any %in% c("up","down") ~ paste0("linear_", direction_any),
      archetype_any == "logistic" & direction_any %in% c("up","down") ~ paste0("logistic_", direction_any),
      archetype_any == "impulse"  & !is.na(polarity_any)              ~ paste0("impulse_", polarity_any),
      archetype_any == "impulse"                                   ~ "impulse",
      TRUE                                                         ~ as.character(archetype_any)
    )
  ) %>%
  select(gene, archetype_label)

mds2 <- mds2 %>% left_join(arch_map, by = "gene")

# ---- make ordering stable in the legend ----
arch_levels <- c("linear_up","linear_down","logistic_up","logistic_down",
                 "impulse_bump","impulse_dip","impulse")
mds2$archetype_label <- factor(mds2$archetype_label, levels = intersect(arch_levels, unique(mds2$archetype_label)))

# =========================
# Plot helpers
# =========================
base_theme <- theme_classic(base_size = 12)

p_cluster <- ggplot(mds2, aes(MDS1, MDS2, color = cluster, alpha = alpha)) +
  geom_point(size = 1.9) +
  scale_alpha_identity() +
  base_theme +
  labs(
    title = "MDS on Approach C distance D",
    subtitle = "Colored by original clusters; alpha = silhouette confidence",
    x = "MDS1", y = "MDS2", color = "Cluster"
  )

p_tpeak <- ggplot(mds2, aes(MDS1, MDS2, color = t_peak, alpha = alpha)) +
  geom_point(size = 1.9) +
  scale_alpha_identity() +
  base_theme +
  labs(
    title = "MDS on Approach C distance D",
    subtitle = "Colored by t_peak (on standardized curve); alpha = silhouette confidence",
    x = "MDS1", y = "MDS2", color = "t_peak"
  )

p_arch <- ggplot(mds2, aes(MDS1, MDS2, color = archetype_label, alpha = alpha)) +
  geom_point(size = 1.9) +
  scale_alpha_identity() +
  base_theme +
  labs(
    title = "MDS on Approach C distance D",
    subtitle = "Colored by best-fit function (up/down split); alpha = silhouette confidence",
    x = "MDS1", y = "MDS2", color = "Best-fit"
  )

# =========================
# Print three versions
# =========================
options(repr.plot.width = 4.8, repr.plot.height = 4.5)
print(p_cluster)

options(repr.plot.width = 5.8, repr.plot.height = 4.5)
print(p_tpeak)

options(repr.plot.width = 5.8, repr.plot.height = 4.5)
print(p_arch)

# =========================
# Extra MDS figure:
#   merge linear/logistic into monotone_up vs monotone_down
#   keep impulse separate (optionally bump/dip if available)
# =========================

stopifnot(exists("mds2"))
stopifnot("archetype_label" %in% names(mds2))

mds3 <- mds2 %>%
  mutate(
    archetype_coarse = case_when(
      archetype_label %in% c("linear_up","logistic_up")     ~ "monotone_up",
      archetype_label %in% c("linear_down","logistic_down") ~ "monotone_down",
      grepl("^impulse", as.character(archetype_label))      ~ as.character(archetype_label), # impulse / impulse_bump / impulse_dip
      TRUE                                                  ~ "other"
    ),
    archetype_coarse = factor(
      archetype_coarse,
      levels = intersect(c("monotone_up","monotone_down","impulse_bump","impulse_dip","impulse","other"),
                         unique(archetype_coarse))
    )
  )

options(repr.plot.width = 4.8, repr.plot.height = 4.5)
print(
  ggplot(mds3, aes(MDS1, MDS2, color = archetype_coarse, alpha = alpha)) +
    geom_point(size = 1.9) +
    scale_alpha_identity() +
    theme_classic(base_size = 12) +
    labs(
      title = "MDS on Approach C distance D",
      subtitle = "Colored by coarse function class (monotone_up/monotone_down/impulse); alpha = silhouette confidence",
      x = "MDS1", y = "MDS2", color = "Coarse class"
    )
)


### CELL 90
# =========================
# Re-plot PCA (PC1–PC2) four ways:
#   (i)  color = original clusters
#   (ii) color = t_peak (computed on standardized curve X)
#   (iii)color = best-fit function (with up/down split)
#   (iv) color = coarse class (monotone_up / monotone_down / impulse*)
# =========================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

stopifnot(exists("pc_df"))
stopifnot(exists("arch_df"))
stopifnot(all(c("gene","cluster","PC1","PC2","alpha","t_peak") %in% names(pc_df)))

# ---- attach archetype label (best-fit function + direction) ----
get_arch_type <- function(df) {
  if ("archetype_final" %in% names(df)) return(as.character(df$archetype_final))
  if ("archetype_raw"   %in% names(df)) return(as.character(df$archetype_raw))
  if ("archetype"       %in% names(df)) return(as.character(df$archetype))
  stop("arch_df must contain one of: archetype_final, archetype_raw, archetype")
}

arch_map <- arch_df %>%
  mutate(
    archetype_any = get_arch_type(.),
    direction_any = if ("direction" %in% names(.)) as.character(direction) else NA_character_,
    polarity_any  = if ("polarity"  %in% names(.)) as.character(polarity)  else NA_character_,
    # Build display label (split monotone by up/down; impulse optionally bump/dip)
    archetype_label = case_when(
      archetype_any == "linear"   & direction_any %in% c("up","down") ~ paste0("linear_", direction_any),
      archetype_any == "logistic" & direction_any %in% c("up","down") ~ paste0("logistic_", direction_any),
      archetype_any == "impulse"  & !is.na(polarity_any)              ~ paste0("impulse_", polarity_any),
      archetype_any == "impulse"                                   ~ "impulse",
      TRUE                                                         ~ as.character(archetype_any)
    )
  ) %>%
  select(gene, archetype_label)

pc2 <- pc_df %>%
  left_join(arch_map, by = "gene")

# ---- stable legend ordering ----
arch_levels <- c("linear_up","linear_down","logistic_up","logistic_down",
                 "impulse_bump","impulse_dip","impulse")
pc2$archetype_label <- factor(pc2$archetype_label, levels = intersect(arch_levels, unique(pc2$archetype_label)))

# ---- coarse categories ----
pc3 <- pc2 %>%
  mutate(
    archetype_coarse = case_when(
      archetype_label %in% c("linear_up","logistic_up")     ~ "monotone_up",
      archetype_label %in% c("linear_down","logistic_down") ~ "monotone_down",
      grepl("^impulse", as.character(archetype_label))      ~ as.character(archetype_label), # impulse / impulse_bump / impulse_dip
      TRUE                                                  ~ "other"
    ),
    archetype_coarse = factor(
      archetype_coarse,
      levels = intersect(c("monotone_up","monotone_down","impulse_bump","impulse_dip","impulse","other"),
                         unique(archetype_coarse))
    )
  )

base_theme <- theme_classic(base_size = 12)

p_cluster <- ggplot(pc3, aes(PC1, PC2, color = cluster, alpha = alpha)) +
  geom_point(size = 1.9) +
  scale_alpha_identity() +
  base_theme +
  labs(
    title = "PCA on standardized fitted curves",
    subtitle = "Colored by original clusters; alpha = silhouette confidence",
    x = "PC1", y = "PC2", color = "Cluster"
  )

p_tpeak <- ggplot(pc3, aes(PC1, PC2, color = t_peak, alpha = alpha)) +
  geom_point(size = 1.9) +
  scale_alpha_identity() +
  base_theme +
  labs(
    title = "PCA on standardized fitted curves",
    subtitle = "Colored by t_peak (on standardized curve); alpha = silhouette confidence",
    x = "PC1", y = "PC2", color = "t_peak"
  )

p_arch <- ggplot(pc3, aes(PC1, PC2, color = archetype_label, alpha = alpha)) +
  geom_point(size = 1.9) +
  scale_alpha_identity() +
  base_theme +
  labs(
    title = "PCA on standardized fitted curves",
    subtitle = "Colored by best-fit function (up/down split); alpha = silhouette confidence",
    x = "PC1", y = "PC2", color = "Best-fit"
  )

p_coarse <- ggplot(pc3, aes(PC1, PC2, color = archetype_coarse, alpha = alpha)) +
  geom_point(size = 1.9) +
  scale_alpha_identity() +
  base_theme +
  labs(
    title = "PCA on standardized fitted curves",
    subtitle = "Colored by coarse function class (monotone_up/monotone_down/impulse); alpha = silhouette confidence",
    x = "PC1", y = "PC2", color = "Coarse class"
  )

# =========================
# Print four versions
# =========================
options(repr.plot.width = 4.8, repr.plot.height = 4.5); print(p_cluster)
# options(repr.plot.width = 4.8, repr.plot.height = 4.5); print(p_tpeak)
# options(repr.plot.width = 5.8, repr.plot.height = 4.5); print(p_arch)
options(repr.plot.width = 5.8, repr.plot.height = 4.5); print(p_coarse)

### CELL 91
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

stopifnot(exists("arch_df"), exists("cluster_C"), exists("traj_mat"))
stopifnot(nrow(traj_mat) == length(cluster_C))

# ---------- robust column pickers ----------
get_arch_type <- function(df) {
  if ("archetype_final" %in% names(df)) return(as.character(df$archetype_final))
  if ("archetype_raw"   %in% names(df)) return(as.character(df$archetype_raw))
  if ("archetype"       %in% names(df)) return(as.character(df$archetype))
  stop("arch_df must contain one of: archetype_final, archetype_raw, archetype")
}

# ---------- build a single summary table keyed by gene ----------
genes <- rownames(traj_mat)
cl_df <- data.frame(gene = genes, cluster = factor(cluster_C), stringsAsFactors = FALSE)

summ_df <- cl_df %>%
  left_join(
    arch_df %>%
      mutate(
        archetype_any = get_arch_type(.),
        direction_any = if ("direction" %in% names(.)) as.character(direction) else NA_character_,
        polarity_any  = if ("polarity"  %in% names(.)) as.character(polarity)  else NA_character_,
        bestfit_label = case_when(
          archetype_any == "linear"   & direction_any %in% c("up","down") ~ paste0("linear_", direction_any),
          archetype_any == "logistic" & direction_any %in% c("up","down") ~ paste0("logistic_", direction_any),
          archetype_any == "impulse"  & !is.na(polarity_any)              ~ paste0("impulse_", polarity_any),
          archetype_any == "impulse"                                   ~ "impulse",
          TRUE                                                         ~ as.character(archetype_any)
        ),
        coarse_class = case_when(
          bestfit_label %in% c("linear_up","logistic_up")     ~ "monotone_up",
          bestfit_label %in% c("linear_down","logistic_down") ~ "monotone_down",
          grepl("^impulse", bestfit_label)                    ~ "impulse",
          TRUE                                                ~ "other"
        )
      ) %>%
      select(gene, bestfit_label, coarse_class),
    by = "gene"
  )

# sanity: keep only genes that have a fit label
summ_df2 <- summ_df %>% filter(!is.na(bestfit_label))

cat("Genes with fits:", nrow(summ_df2), " / ", nrow(summ_df), "\n\n")

# ---------- proportions ----------
prop_tbl <- function(x, name) {
  tibble(!!name := x) %>%
    count(!!sym(name), name = "n") %>%
    mutate(prop = n / sum(n)) %>%
    arrange(desc(n))
}

cat("=== Overall proportions: best-fit classes (with up/down + optional bump/dip) ===\n")
print(prop_tbl(summ_df2$bestfit_label, "bestfit_label"))

cat("\n=== Overall proportions: original cluster assignments ===\n")
print(prop_tbl(summ_df2$cluster, "cluster"))

cat("\n=== Overall proportions: coarse classes ===\n")
print(prop_tbl(summ_df2$coarse_class, "coarse_class"))

### CELL 95
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

stopifnot(exists("traj_mat1"), exists("traj_mat2a"), exists("traj_mat2b"),
          exists("combined1"), exists("combined2"),
          exists("z_row"))

# =========================
# CONFIG
# =========================
TOPK_LABELS  <- 100
LABEL_SIZE   <- 2.6
SHOW_LABELS  <- TRUE
POINT_SIZE   <- 1.0

# log2FC coloring
FC_CAP    <- 3     # clip log2FC to [-FC_CAP, FC_CAP]
FC_GAMMA  <- 1     # smaller = more aggressive, larger = more white

options(repr.plot.width = 13.0, repr.plot.height = 6.2)

# =========================
# 0) Align genes across all 3 matrices
# =========================
common_genes <- Reduce(intersect, list(rownames(traj_mat1), rownames(traj_mat2a), rownames(traj_mat2b)))
stopifnot(length(common_genes) > 0)

traj1  <- traj_mat1 [common_genes, , drop = FALSE]   # used only to TRAIN pc_orig
traj2a <- traj_mat2a[common_genes, , drop = FALSE]
traj2b <- traj_mat2b[common_genes, , drop = FALSE]

# =========================
# 1) Fit ORIGINAL PCA on traj1 z-curves only (if not already present)
# =========================
if (!exists("pc_orig")) {
  X1_z <- z_row(traj1)
  storage.mode(X1_z) <- "double"
  pc_orig <- prcomp(X1_z, center = FALSE, scale. = FALSE)
}

time_cols_train <- rownames(pc_orig$rotation)
stopifnot(length(time_cols_train) > 0)

# =========================
# 2) Project helper into ORIGINAL PCA space
# =========================
project_into_pc_orig <- function(traj_mat, pc_model, time_cols_train) {
  Xz <- z_row(traj_mat)
  storage.mode(Xz) <- "double"

  if (ncol(Xz) != length(time_cols_train)) {
    stop(sprintf("Time dimension mismatch: ncol(Xz)=%d but pc_orig expects %d time features",
                 ncol(Xz), length(time_cols_train)))
  }

  colnames(Xz) <- time_cols_train
  S <- Xz %*% pc_model$rotation
  S <- as.matrix(S)
  colnames(S) <- colnames(pc_model$x)
  S
}

S2a <- project_into_pc_orig(traj2a, pc_orig, time_cols_train)
S2b <- project_into_pc_orig(traj2b, pc_orig, time_cols_train)

# =========================
# 3) Dataset-specific log2FC to color by
# =========================
fc2a <- as.numeric(combined2[common_genes, "log2FC_t1"])
fc2b <- as.numeric(combined2[common_genes, "log2FC_t2"])

# =========================
# 4) Build plotting df (ONLY 2A + 2B)
# =========================
mk_df <- function(S, fc, tag) {
  data.frame(
    gene   = common_genes,
    PC1    = S[,1],
    PC2    = S[,2],
    log2FC = fc,
    panel  = tag,
    stringsAsFactors = FALSE
  )
}

pc_all <- bind_rows(
  mk_df(S2a, fc2a, "traj_mat2a"),
  mk_df(S2b, fc2b, "traj_mat2b")
)

# =========================
# Color transform (clip -> scale to [-1,1] -> gamma contrast)
# =========================
pc_all$log2FC_clip <- pmax(-FC_CAP, pmin(FC_CAP, pc_all$log2FC))
x <- pc_all$log2FC_clip / FC_CAP
pc_all$log2FC_col <- sign(x) * (abs(x) ^ FC_GAMMA)

# =========================
# 5) Robust axis limits
# =========================
robust_lim <- function(x, q = c(0.002, 0.998), pad = 0.08) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(-1, 1))
  lohi <- as.numeric(stats::quantile(x, probs = q, names = FALSE, na.rm = TRUE))
  r <- diff(lohi); if (!is.finite(r) || r == 0) r <- 1
  lohi + c(-pad, pad) * r
}
xlim <- robust_lim(pc_all$PC1)
ylim <- robust_lim(pc_all$PC2)

# =========================
# 6) Label SAME top K genes in both panels: choose from traj_mat2a only
#    (ranking by |log2FC| in traj_mat2a)
# =========================
top_genes <- pc_all %>%
  filter(panel == "traj_mat2a", is.finite(log2FC)) %>%
  slice_max(order_by = abs(log2FC), n = TOPK_LABELS, with_ties = FALSE) %>%
  pull(gene) %>%
  unique()

dlab <- pc_all %>%
  filter(gene %in% top_genes)

# =========================
# 7) Plot: ORIGINAL PCA space, color by log2FC (contrast-tuned), 1x2 facets
# =========================
p <- ggplot(pc_all, aes(PC1, PC2, color = log2FC_col)) +
  geom_point(size = POINT_SIZE, alpha = 1.0) +
  coord_cartesian(xlim = xlim, ylim = ylim) +
  facet_wrap(~panel, nrow = 1) +
  theme_classic(base_size = 11) +
  labs(
    title = sprintf("ORIGINAL PCA (shared pc_orig) | log2FC color (FC_CAP=%.2f, FC_GAMMA=%.2f) | top %d labels from traj_mat2a",
                    FC_CAP, FC_GAMMA, TOPK_LABELS),
    x = "PC1", y = "PC2", color = "log2FC"
  ) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    breaks = c(-1, -0.5, 0, 0.5, 1),
    labels = sprintf("%.1f", c(-1, -0.5, 0, 0.5, 1) * FC_CAP)
  )

if (isTRUE(SHOW_LABELS) && nrow(dlab) > 0) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = dlab,
      aes(label = gene),
      size = LABEL_SIZE,
      max.overlaps = Inf,
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0
    )
  } else {
    p <- p + geom_text(
      data = dlab,
      aes(label = gene),
      size = LABEL_SIZE,
      check_overlap = TRUE
    )
  }
}

print(p)

### CELL 96
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

stopifnot(exists("pc_all"))
stopifnot(all(c("gene","panel","PC1","PC2") %in% colnames(pc_all)))

TOPN <- 50

pc_delta <- pc_all %>%
  select(gene, panel, PC1, PC2) %>%
  filter(panel %in% c("traj_mat2a", "traj_mat2b")) %>%
  distinct(gene, panel, .keep_all = TRUE) %>%      # keep 1 row per gene×panel
  pivot_wider(
    names_from  = panel,
    values_from = c(PC1, PC2),
    names_sep   = "__"
  ) %>%
  # keep only genes present in both panels
  filter(!is.na(PC1__traj_mat2a), !is.na(PC1__traj_mat2b),
         !is.na(PC2__traj_mat2a), !is.na(PC2__traj_mat2b)) %>%
  mutate(
    dPC1 = PC1__traj_mat2b - PC1__traj_mat2a,
    dPC2 = PC2__traj_mat2b - PC2__traj_mat2a,
    abs_dPC1 = abs(dPC1),
    abs_dPC2 = abs(dPC2),
    # “largest change” = Euclidean shift in (PC1, PC2)
    delta_L2 = sqrt(dPC1^2 + dPC2^2)
  )

# 1) Largest overall (PC1,PC2) shift
top_overall <- pc_delta %>%
  arrange(desc(delta_L2)) %>%
  slice_head(n = TOPN)

# 2) Largest PC1-only shift
top_pc1 <- pc_delta %>%
  arrange(desc(abs_dPC1)) %>%
  slice_head(n = TOPN)

# 3) Largest PC2-only shift
top_pc2 <- pc_delta %>%
  arrange(desc(abs_dPC2)) %>%
  slice_head(n = TOPN)

top_overall
top_pc1
top_pc2

### CELL 97
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# =========================
# REQUIREMENTS
# =========================
stopifnot(exists("pc_all"), exists("asso_df2"))
stopifnot(all(c("gene","panel","PC1","PC2") %in% colnames(pc_all)))
stopifnot(all(c("waldStat","waldStat_1","waldStat_2") %in% colnames(asso_df2)))

# =========================
# CONFIG
# =========================
SIG_Q          <- 0.20
TOPK_WALD      <- 600    # top-k by max waldStat
TOPK_DELTA     <- 600    # top-k by delta_L2
POINT_SIZE     <- 1.8
POINT_ALPHA    <- 0.85
LABEL_SIZE     <- 3.0
WALD_EPS       <- 1e-12

# bigger plot
options(repr.plot.width = 10.5, repr.plot.height = 8.0)

# =========================
# 1) delta_L2 per gene between traj_mat2a vs traj_mat2b
# =========================
pc_delta <- pc_all %>%
  select(gene, panel, PC1, PC2) %>%
  filter(panel %in% c("traj_mat2a", "traj_mat2b")) %>%
  distinct(gene, panel, .keep_all = TRUE) %>%
  pivot_wider(
    names_from  = panel,
    values_from = c(PC1, PC2),
    names_sep   = "__"
  ) %>%
  filter(!is.na(PC1__traj_mat2a), !is.na(PC1__traj_mat2b),
         !is.na(PC2__traj_mat2a), !is.na(PC2__traj_mat2b)) %>%
  mutate(
    dPC1     = PC1__traj_mat2b - PC1__traj_mat2a,
    dPC2     = PC2__traj_mat2b - PC2__traj_mat2a,
    delta_L2 = sqrt(dPC1^2 + dPC2^2)
  ) %>%
  select(gene, delta_L2, dPC1, dPC2)

# =========================
# 2) max Wald stat per gene + significance
# =========================
asso_tbl <- asso_df2 %>%
  tibble::as_tibble(rownames = "gene") %>%
  mutate(
    wald_max = pmax(waldStat, waldStat_1, waldStat_2, na.rm = TRUE),
    sig      = if ("qval" %in% names(.)) (qval <= SIG_Q) else if ("pvalue" %in% names(.)) (pvalue <= SIG_Q) else FALSE,
    log10_wald = log10(wald_max + WALD_EPS)
  )

plot_df <- asso_tbl %>%
  inner_join(pc_delta, by = "gene")

# =========================
# 3) label genes = intersection( topK by wald_max, topK by delta_L2 ), within significant
# =========================
top_by_wald <- plot_df %>%
  filter(sig) %>%
  arrange(desc(wald_max)) %>%
  slice_head(n = TOPK_WALD) %>%
  pull(gene)

top_by_delta <- plot_df %>%
  filter(sig) %>%
  arrange(desc(delta_L2)) %>%
  slice_head(n = TOPK_DELTA) %>%
  pull(gene)

label_genes <- intersect(top_by_wald, top_by_delta)

# =========================
# 4) plot
# =========================
p <- ggplot(plot_df, aes(x = log10_wald, y = delta_L2)) +
  geom_point(aes(alpha = sig), size = POINT_SIZE) +
  scale_alpha_manual(values = c(`TRUE` = 1.0, `FALSE` = 0.18), guide = "none") +
  labs(
    x = "log10(max waldStat)",
    y = "PC shift (delta_L2): traj_mat2a → traj_mat2b",
    title = "Association strength vs. PC shift"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 26, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 14)
  )

if (length(label_genes) > 0) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = plot_df %>% filter(gene %in% label_genes),
      aes(label = gene),
      size = LABEL_SIZE,
      max.overlaps = Inf,
      box.padding = 0.30,
      point.padding = 0.18,
      min.segment.length = 0
    )
  } else {
    p <- p + geom_text(
      data = plot_df %>% filter(gene %in% label_genes),
      aes(label = gene),
      size = LABEL_SIZE,
      vjust = -0.6
    )
  }
} else {
  message("No genes in the intersection. Try increasing TOPK_WALD and/or TOPK_DELTA.")
}

print(p)

### CELL 98
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# =========================
# REQUIREMENTS
# =========================
stopifnot(exists("pc_all"), exists("asso_df2"))
stopifnot(all(c("gene","panel","PC1","PC2") %in% colnames(pc_all)))
stopifnot(all(c("waldStat","waldStat_1","waldStat_2") %in% colnames(asso_df2)))

# =========================
# CONFIG
# =========================
SIG_Q          <- 0.20
TOPK_WALD      <- 200    # now: top-k by |waldStat_1 - waldStat_2|
TOPK_DELTA     <- 200    # top-k by delta_L2
POINT_SIZE     <- 1.8
POINT_ALPHA    <- 0.85
LABEL_SIZE     <- 3.0
WALD_EPS       <- 1e-12

# bigger plot
options(repr.plot.width = 10.5, repr.plot.height = 8.0)

# =========================
# 1) delta_L2 per gene between traj_mat2a vs traj_mat2b
# =========================
pc_delta <- pc_all %>%
  select(gene, panel, PC1, PC2) %>%
  filter(panel %in% c("traj_mat2a", "traj_mat2b")) %>%
  distinct(gene, panel, .keep_all = TRUE) %>%
  pivot_wider(
    names_from  = panel,
    values_from = c(PC1, PC2),
    names_sep   = "__"
  ) %>%
  filter(!is.na(PC1__traj_mat2a), !is.na(PC1__traj_mat2b),
         !is.na(PC2__traj_mat2a), !is.na(PC2__traj_mat2b)) %>%
  mutate(
    dPC1     = PC1__traj_mat2b - PC1__traj_mat2a,
    dPC2     = PC2__traj_mat2b - PC2__traj_mat2a,
    delta_L2 = sqrt(dPC1^2 + dPC2^2)
  ) %>%
  select(gene, delta_L2, dPC1, dPC2)

# =========================
# 2) |waldStat_1 - waldStat_2| for x-axis + significance
# =========================
# NOTE on log space:
# We compute the difference in *linear* Wald space, then log10 it for plotting.
# (This answers “how big is the magnitude gap between the two Wald stats?”)
asso_tbl <- asso_df2 %>%
  tibble::as_tibble(rownames = "gene") %>%
  mutate(
    wald_diff_abs = abs(waldStat_1 - waldStat_2),
    sig      = if ("qval" %in% names(.)) (qval <= SIG_Q) else if ("pvalue" %in% names(.)) (pvalue <= SIG_Q) else FALSE,
    log10_wald_diff = log10(wald_diff_abs + WALD_EPS)
  )

plot_df <- asso_tbl %>%
  inner_join(pc_delta, by = "gene")

# =========================
# 3) label genes = intersection( topK by |wald_1-wald_2|, topK by delta_L2 ), within significant
# =========================
top_by_wald <- plot_df %>%
  filter(sig) %>%
  arrange(desc(wald_diff_abs)) %>%
  slice_head(n = TOPK_WALD) %>%
  pull(gene)

top_by_delta <- plot_df %>%
  filter(sig) %>%
  arrange(desc(delta_L2)) %>%
  slice_head(n = TOPK_DELTA) %>%
  pull(gene)

label_genes <- intersect(top_by_wald, top_by_delta)

# =========================
# 4) plot
# =========================
p <- ggplot(plot_df, aes(x = log10_wald_diff, y = delta_L2)) +
  geom_point(aes(alpha = sig), size = POINT_SIZE) +
  scale_alpha_manual(values = c(`TRUE` = 1.0, `FALSE` = 0.18), guide = "none") +
  labs(
    x = "log10(|waldStat_1 - waldStat_2|)",
    y = "PC shift (delta_L2): traj_mat2a → traj_mat2b",
    title = "Association-difference vs. PC shift"
  ) +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(size = 26, face = "bold"),
    axis.title = element_text(size = 20),
    axis.text  = element_text(size = 14)
  )

if (length(label_genes) > 0) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = plot_df %>% filter(gene %in% label_genes),
      aes(label = gene),
      size = LABEL_SIZE,
      max.overlaps = Inf,
      box.padding = 0.30,
      point.padding = 0.18,
      min.segment.length = 0
    )
  } else {
    p <- p + geom_text(
      data = plot_df %>% filter(gene %in% label_genes),
      aes(label = gene),
      size = LABEL_SIZE,
      vjust = -0.6
    )
  }
} else {
  message("No genes in the intersection. Try increasing TOPK_WALD and/or TOPK_DELTA.")
}

print(p)

### CELL 99
options(repr.plot.width = 2.75, repr.plot.height = 2.5)

plot_tradeSeq_gene_pub(
  "LTB", gamList2, sce2, counts_sub, PT_mat, sf,
  cellWeights = CWz_ret,
  lineages = c(1,2),
  colours = c("#386cb0", "#DC143C"),
  time_max = 110,
  x_tick_step = 20,        # <-- every 20 min
  color_points_by_lineage = TRUE,
  points_size = 0.75,
  points_alpha = 0.3,
  axis_text_size = 14,
  axis_title_size = 16,
  # ylim = c(0, 0.1),
  lineage_labels = c("Inflammatory-hi", "Erythroid-like"),
  show_legend = FALSE,
  save_pdf = FALSE
)

### CELL 102
X_z <- z_row(traj_mat1)
storage.mode(X_z) <- "double"

pc_orig <- prcomp(X_z, center = FALSE, scale. = FALSE)

pc_df_orig <- meta_df
pc_df_orig$PC1 <- pc_orig$x[, 1]
pc_df_orig$PC2 <- pc_orig$x[, 2]
pc_df_orig$PC3 <- pc_orig$x[, 3]

peak_idx_z <- max.col(X_z, ties.method = "first")
pc_df_orig$t_peak <- tt[peak_idx_z]

### CELL 103
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ---- build pc_orig if needed (matches your code) ----
if (!exists("pc_orig")) {
  stopifnot(exists("traj_mat1"), exists("z_row"))
  X_z <- z_row(traj_mat1)
  storage.mode(X_z) <- "double"
  pc_orig <- prcomp(X_z, center = FALSE, scale. = FALSE)
}

# ---- variance explained ----
ve <- (pc_orig$sdev^2) / sum(pc_orig$sdev^2)
df_ve <- tibble(
  PC = seq_along(ve),
  var_explained = ve,
  cum_explained = cumsum(ve)
)

# ---- plot 1: scree (individual variance explained) ----
p1 <- ggplot(df_ve, aes(x = PC, y = var_explained)) +
  geom_line() +
  geom_point(size = 1.6) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
  labs(x = "PC dimension", y = "Variance explained", title = "Scree plot") +
  theme_classic(base_size = 14)

# ---- plot 2: cumulative variance explained ----
p2 <- ggplot(df_ve, aes(x = PC, y = cum_explained)) +
  geom_line() +
  geom_point(size = 1.6) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1), limits = c(0, 1)) +
  labs(x = "PC dimension", y = "Cumulative variance explained", title = "Cumulative variance") +
  theme_classic(base_size = 14)

print(p1)
print(p2)

# (optional) quick table
df_ve %>% mutate(across(c(var_explained, cum_explained), ~ round(.x, 4))) %>% head(20)

### CELL 104
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

stopifnot(exists("asso_df2"))
stopifnot(all(c("waldStat_1","waldStat_2") %in% colnames(asso_df2)))

# =========================
# CONFIG
# =========================
TOPK_LABELS   <- 50
POINT_SIZE    <- 1.6
POINT_ALPHA   <- 0.55
HILITE_SIZE   <- 2.6
LABEL_SIZE    <- 3.0

# make plot bigger
options(repr.plot.width = 9.5, repr.plot.height = 8.0)

# =========================
# PREP
# =========================
df <- asso_df2 %>%
  tibble::as_tibble(rownames = "gene") %>%
  mutate(
    w1 = as.numeric(waldStat_1),
    w2 = as.numeric(waldStat_2),
    d  = w2 - w1,
    abs_d = abs(d)
  ) %>%
  filter(is.finite(w1), is.finite(w2))

# top-k by absolute difference
top_genes <- df %>%
  arrange(desc(abs_d)) %>%
  slice_head(n = TOPK_LABELS) %>%
  pull(gene)

df <- df %>% mutate(is_top = gene %in% top_genes)

# =========================
# PLOT
# =========================
p <- ggplot(df, aes(x = w1, y = w2)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.7) +
  geom_point(alpha = POINT_ALPHA, size = POINT_SIZE) +
  geom_point(
    data = df %>% filter(is_top),
    size = HILITE_SIZE
  ) +
  labs(
    x = "waldStat_1",
    y = "waldStat_2",
    title = sprintf("waldStat_2 vs waldStat_1 (top %d by |difference| labeled)", TOPK_LABELS)
  ) +
  theme_classic(base_size = 16)

# labels (ggrepel if available)
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p <- p + ggrepel::geom_text_repel(
    data = df %>% filter(is_top),
    aes(label = gene),
    size = LABEL_SIZE,
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.15,
    min.segment.length = 0
  )
} else {
  p <- p + geom_text(
    data = df %>% filter(is_top),
    aes(label = gene),
    size = LABEL_SIZE,
    vjust = -0.6
  )
}

print(p)

# optional: inspect the top table (robust to missing columns)
df %>%
  arrange(desc(abs_d)) %>%
  mutate(
    pvalue_1 = if ("pvalue_1" %in% names(.)) .data[["pvalue_1"]] else NA_real_,
    pvalue_2 = if ("pvalue_2" %in% names(.)) .data[["pvalue_2"]] else NA_real_,
    qval_1   = if ("qval_1"   %in% names(.)) .data[["qval_1"]]   else NA_real_,
    qval_2   = if ("qval_2"   %in% names(.)) .data[["qval_2"]]   else NA_real_
  ) %>%
  transmute(
    gene,
    waldStat_1 = w1,
    waldStat_2 = w2,
    diff       = d,
    abs_diff   = abs_d,
    pvalue_1, pvalue_2, qval_1, qval_2
  ) %>%
  head(10)

### CELL 105
save.image()
