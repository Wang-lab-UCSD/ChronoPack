# ChronoPack dimensionality reduction and clustering helpers.
# Extracted from the consolidated puromycin DimensionalityReduction&Clustering notebook.

.get_traj_from_predictSmooth <- function(
  gamList,
  features,
  nPoints = 200,
  lineage_pick = 1,
  verbose = TRUE,
  sce = NULL,
  sizeFactors = NULL,
  # optional override (rarely needed): directly provide reference offset
  off_ref = NULL,
  return_log1p = TRUE,
  # if TRUE and sce is provided, try to compute curves via predict(gam, ...) using dm + offset
  prefer_predict_gam = TRUE
) {
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
  # 0) Reference offset logic (matches)
  # -------------------------------
  .compute_off_ref <- function(sizeFactors, off_ref) {
    if (!is.null(off_ref)) return(as.numeric(off_ref)[1])
    if (is.null(sizeFactors)) return(NULL)
    sf_vec <- as.numeric(sizeFactors)
    if (!length(sf_vec) || all(!is.finite(sf_vec))) return(NULL)
    s_ref <- exp(mean(log(sf_vec), na.rm = TRUE))
    log(s_ref)
  }
  off_ref_use <- .compute_off_ref(sizeFactors = sizeFactors, off_ref = off_ref)

  # -------------------------------
  #     predict(gam, newdata=.getPredictRangeDf(dm,...), offset=off_ref)
  #     -> mu = log1p(exp(eta)) or exp(eta) depending on return_log1p
  # -------------------------------
  if (isTRUE(prefer_predict_gam) && !is.null(sce) && is.finite(off_ref_use)) {
    dm <- tryCatch(SummarizedExperiment::colData(sce)$tradeSeq$dm, error = function(e) NULL)
    if (!is.null(dm)) {
      t_cols <- grep("^t[0-9]+$", colnames(dm), value = TRUE)
      avail_lineages <- as.integer(sub("^t", "", t_cols))
      if (length(avail_lineages) > 0L) {
        ell <- as.integer(lineage_pick)
        if (!(ell %in% avail_lineages)) ell <- avail_lineages[1]
        if (verbose) message("Using lineage = ", ell, " (predict(gam,...) w/ reference offset)")

        # build t-grid
        df_grid <- tradeSeq:::.getPredictRangeDf(dm, lineageId = ell, nPoints = nPoints)
        tcol <- paste0("t", ell)
        tt <- df_grid[[tcol]]

        traj_mat <- matrix(NA_real_, nrow = 0, ncol = length(tt))
        rownames(traj_mat) <- character(0)
        colnames(traj_mat) <- as.character(tt)

        for (g in features) {
          if (!g %in% names(gamList)) next
          gm <- gamList[[g]]
          pr <- tryCatch(
            predict(gm, newdata = df_grid, type = "link", se.fit = FALSE,
                    offset = rep(off_ref_use, nrow(df_grid))),
            error = function(e) NULL
          )
          if (is.null(pr)) next

          # pr might be numeric vector or list; standardize to eta
          eta <- if (is.list(pr) && !is.null(pr$fit)) pr$fit else as.numeric(pr)

          mu <- exp(eta)                # mean on count scale at reference offset
          y  <- if (isTRUE(return_log1p)) log1p(mu) else mu

          traj_mat <- rbind(traj_mat, y)
          rownames(traj_mat)[nrow(traj_mat)] <- g
        }

        if (nrow(traj_mat) > 0) {
          return(list(traj_mat = traj_mat, tt = as.numeric(tt)))
        }
        # if this path fails to produce anything, fall through to predictSmooth with correction below
      }
    }
  }

  # -------------------------------
  # 1) try tidy interface (predictSmooth) and then re-anchor to reference offset
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
      lp <- if (any(lin_levels %in% c(lineage_pick, as.character(lineage_pick))))
        lin_levels[lin_levels %in% c(lineage_pick, as.character(lineage_pick))][1]
      else lin_levels[1]
      pr <- pr[pr[[lin_col]] == lp, , drop = FALSE]
      if (verbose) message("Using lineage = ", lp, " (tidy predictSmooth)")
    }

    tt <- sort(unique(pr[[time_col]]))

    traj_df <- data.frame(
      feature = as.character(pr[[gene_col]]),
      t       = pr[[time_col]],
      y       = as.numeric(pr[[y_col]]),
      stringsAsFactors = FALSE
    )
    traj_df <- traj_df[order(traj_df$feature, traj_df$t), , drop = FALSE]

    # keep only genes that have all timepoints
    tab <- table(traj_df$feature)
    keep <- names(tab)[tab == length(tt)]
    traj_df <- traj_df[traj_df$feature %in% keep, , drop = FALSE]

    # pivot (base R)
    wide <- reshape(traj_df, idvar = "feature", timevar = "t", direction = "wide")
    rn <- wide$feature
    wide$feature <- NULL
    colnames(wide) <- sub("^y\\.", "", colnames(wide))
    traj_mat <- as.matrix(wide)
    rownames(traj_mat) <- rn

    # --- re-anchor to reference offset ---
    # predictSmooth uses an implicit offset (often the first cell). If we can get the
    # original offset used by tradeSeq, we can adjust: mu_ref = mu_ps * exp(off_ref - off_used)
    if (is.finite(off_ref_use) && !is.null(sce)) {
      dm <- tryCatch(SummarizedExperiment::colData(sce)$tradeSeq$dm, error = function(e) NULL)
      off_used <- if (!is.null(dm) && "offset" %in% colnames(dm)) as.numeric(dm$offset[1]) else NA_real_

      if (is.finite(off_used)) {
        if (verbose) message(sprintf("Re-anchoring predictSmooth using exp(off_ref - off_used) = exp(%.4f - %.4f)",
                                     off_ref_use, off_used))
        # traj_mat is on count scale (mu); rescale to reference offset
        traj_mat <- traj_mat * exp(off_ref_use - off_used)
      } else if (verbose) {
        message("Could not find dm$offset[1]; returning predictSmooth scale unchanged.")
      }
    } else if (verbose) {
      message("No reference offset available (need sizeFactors/off_ref and sce); returning predictSmooth scale unchanged.")
    }

    if (isTRUE(return_log1p)) traj_mat <- log1p(traj_mat)

    tt_num <- suppressWarnings(as.numeric(colnames(traj_mat)))
    if (any(!is.finite(tt_num))) tt_num <- seq(0, 1, length.out = ncol(traj_mat))
    return(list(traj_mat = traj_mat, tt = tt_num))
  }

  # -------------------------------
  # 2) fallback interface (many versions / shapes) from predictSmooth, then re-anchor if possible
  # -------------------------------
  ps_raw <- tryCatch(
    tradeSeq::predictSmooth(gamList, gene = features, nPoints = nPoints),
    error = function(e) e
  )
  if (inherits(ps_raw, "error")) stop("predictSmooth() failed: ", ps_raw$message)

  ps <- ps_raw
  if (is.list(ps) && !is.null(ps$yhat)) ps <- ps$yhat
  if (is.list(ps) && length(ps) == 1) ps <- ps[[1]]

  # helper to apply offset correction + scale
  .postprocess_mat <- function(mat_gene_time, tt) {
    # re-anchor if possible
    if (is.finite(off_ref_use) && !is.null(sce)) {
      dm <- tryCatch(SummarizedExperiment::colData(sce)$tradeSeq$dm, error = function(e) NULL)
      off_used <- if (!is.null(dm) && "offset" %in% colnames(dm)) as.numeric(dm$offset[1]) else NA_real_
      if (is.finite(off_used)) {
        if (verbose) message(sprintf("Re-anchoring predictSmooth using exp(off_ref - off_used) = exp(%.4f - %.4f)",
                                     off_ref_use, off_used))
        mat_gene_time <- mat_gene_time * exp(off_ref_use - off_used)
      } else if (verbose) {
        message("Could not find dm$offset[1]; returning predictSmooth scale unchanged.")
      }
    } else if (verbose) {
      message("No reference offset available (need sizeFactors/off_ref and sce); returning predictSmooth scale unchanged.")
    }

    if (isTRUE(return_log1p)) mat_gene_time <- log1p(mat_gene_time)

    list(traj_mat = mat_gene_time, tt = tt)
  }

  # -------------------------------
  # case A: 3D array: time x gene x lineage  OR gene x time x lineage
  # -------------------------------
  if (is.array(ps) && length(dim(ps)) == 3) {
    d <- dim(ps)
    lin <- min(lineage_pick, d[3])
    ps2 <- ps[, , lin, drop = TRUE]
    if (verbose) message("Using lineage slice = ", lin, " (3D array)")

    rn <- rownames(ps2); cn <- colnames(ps2)
    ov_r <- if (!is.null(rn)) sum(rn %in% features) else 0
    ov_c <- if (!is.null(cn)) sum(cn %in% features) else 0
    mat_gene_time <- if (ov_r >= ov_c) ps2 else t(ps2)

    keep <- intersect(features, rownames(mat_gene_time))
    if (!length(keep)) {
      if (verbose) message("Diagnostics rownames: ", paste(head(rownames(mat_gene_time), 10), collapse = ", "))
      stop("No matching genes after slicing 3D array output. Likely ID mismatch.")
    }
    mat_gene_time <- mat_gene_time[keep, , drop = FALSE]

    tt <- attr(ps_raw, "t")
    if (is.null(tt) || length(tt) != ncol(mat_gene_time)) {
      tt <- if (ncol(mat_gene_time) > 1) (0:(ncol(mat_gene_time)-1)) / (ncol(mat_gene_time)-1) else 0
    }
    return(.postprocess_mat(mat_gene_time, tt))
  }

  # -------------------------------
  # case B: list of lineage-specific matrices
  # -------------------------------
  if (is.list(ps) && length(ps) >= 1 &&
      all(vapply(ps, function(x) is.matrix(x) || (is.array(x) && length(dim(x))==2), TRUE))) {
    lin <- min(lineage_pick, length(ps))
    if (verbose) message("Using element ", lin, " from list predictSmooth output (list-of-matrices).")
    ps <- ps[[lin]]
  }

  if (!(is.matrix(ps) || (is.array(ps) && length(dim(ps)) == 2))) {
    stop("Unsupported predictSmooth() output type. Class: ", paste(class(ps), collapse = ", "),
         " dim: ", paste(dim(ps), collapse = "x"))
  }
  ps <- as.matrix(ps)

  rn <- rownames(ps); cn <- colnames(ps)
  ov_r <- if (!is.null(rn)) sum(rn %in% features) else 0
  ov_c <- if (!is.null(cn)) sum(cn %in% features) else 0
  mat_gene_timeish <- if (ov_r >= ov_c) ps else t(ps)

  # concatenated multi-lineage columns like "lineage1_1"
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
    return(.postprocess_mat(mat_gene_timeish, tt))
  }

  # Otherwise: assume genes x time
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
  return(.postprocess_mat(mat_gene_timeish, tt))
}

get_descriptors <- function(
  traj_mat, tt,
  pseudocount = 1e-5,
  N_minutes = 10,
  z_diff_threshold = 0.10,
  require_consecutive = 3,
  onset_mode = c("early_lm", "horizontal"),
  owp_onset_fraction = 0.05
) {
  onset_mode <- match.arg(onset_mode)

  trapz_simple <- function(x, y) {
    if (length(x) < 2) return(0)
    sum(diff(x) * (y[-1] + y[-length(y)]) / 2)
  }

  # OWP = positive area / area under running-max envelope.
  # Expression: subtract initial value first.
  # Derivative: zero is the natural baseline.
  compute_owp <- function(y, t, onset_fraction = 0.05, baseline = c("first", "zero")) {
    baseline <- match.arg(baseline)
    keep <- is.finite(y) & is.finite(t)
    y <- y[keep]; t <- t[keep]
    if (length(y) < 2) return(NA_real_)

    if (baseline == "first") y <- y - y[1]
    yp <- pmax(y, 0)

    ymax <- max(yp, na.rm = TRUE)
    if (!is.finite(ymax) || ymax <= 0) return(NA_real_)

    onset <- which(yp >= onset_fraction * ymax)[1]
    if (is.na(onset) || onset >= length(yp)) return(NA_real_)

    yp <- yp[onset:length(yp)]
    t  <- t[onset:length(t)]
    env <- cummax(yp)

    denom <- trapz_simple(t, env)
    if (!is.finite(denom) || denom <= 0) return(NA_real_)

    trapz_simple(t, yp) / denom
  }

  calc_shape_metrics <- function(y, tt) {
    if (length(tt) < 3) return(list(mean_curv = NA_real_, max_curv = NA_real_, mean_turn = NA_real_, max_turn = NA_real_))

    tmin <- min(tt); tmax <- max(tt)
    if (!is.finite(tmin) || !is.finite(tmax) || tmax <= tmin)
      return(list(mean_curv = NA_real_, max_curv = NA_real_, mean_turn = NA_real_, max_turn = NA_real_))

    x <- (tt - tmin) / (tmax - tmin)
    ymin <- min(y, na.rm = TRUE); ymax <- max(y, na.rm = TRUE)
    if (!is.finite(ymin) || !is.finite(ymax) || ymax <= ymin)
      return(list(mean_curv = 0, max_curv = 0, mean_turn = 0, max_turn = 0))

    z <- (y - ymin) / (ymax - ymin)
    n <- length(x)

    dx_l <- x[2:(n - 1)] - x[1:(n - 2)]
    dx_r <- x[3:n] - x[2:(n - 1)]
    slope_l <- (z[2:(n - 1)] - z[1:(n - 2)]) / dx_l
    slope_r <- (z[3:n] - z[2:(n - 1)]) / dx_r
    dz  <- (z[3:n] - z[1:(n - 2)]) / (x[3:n] - x[1:(n - 2)])
    d2z <- 2 * (slope_r - slope_l) / (dx_l + dx_r)
    kappa <- abs(d2z) / (1 + dz^2)^(3/2)

    u <- cbind(x[2:(n - 1)] - x[1:(n - 2)], z[2:(n - 1)] - z[1:(n - 2)])
    v <- cbind(x[3:n] - x[2:(n - 1)], z[3:n] - z[2:(n - 1)])
    denom <- sqrt(rowSums(u^2)) * sqrt(rowSums(v^2))

    cosang <- rep(1, length(denom))
    ok <- is.finite(denom) & denom > 0
    cosang[ok] <- rowSums(u * v)[ok] / denom[ok]
    turn_angle <- acos(pmax(-1, pmin(1, cosang))) / pi

    list(
      mean_curv = mean(kappa, na.rm = TRUE),
      max_curv = max(kappa, na.rm = TRUE),
      mean_turn = mean(turn_angle, na.rm = TRUE),
      max_turn = max(turn_angle, na.rm = TRUE)
    )
  }

  calc_onset_zdiff <- function(y, tt, N_minutes = 10, z_diff_threshold = 0.25,
                               require_consecutive = 1, mode = c("early_lm", "horizontal")) {
    mode <- match.arg(mode)
    keep <- is.finite(y) & is.finite(tt)
    y <- y[keep]; tt <- tt[keep]
    if (length(y) < 2) return(NA_real_)

    y_sd <- sd(y, na.rm = TRUE)
    if (!is.finite(y_sd) || y_sd == 0) y_sd <- 1
    df <- data.frame(time = tt, expr_z = (y - mean(y, na.rm = TRUE)) / y_sd)
    early_df <- df[df$time <= N_minutes, , drop = FALSE]
    if (nrow(early_df) < 2) return(NA_real_)

    if (mode == "early_lm") {
      fit <- lm(expr_z ~ time, data = early_df)
      df$blue_fit <- predict(fit, newdata = df)
    } else {
      df$blue_fit <- mean(early_df$expr_z, na.rm = TRUE)
    }

    above <- abs(df$expr_z - df$blue_fit) > z_diff_threshold
    if (require_consecutive <= 1) {
      onset_idx <- which(above)[1]
    } else {
      r <- rle(above)
      ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
      good <- which(r$values & r$lengths >= require_consecutive)
      onset_idx <- if (length(good)) starts[good[1]] else NA_integer_
    }

    if (is.na(onset_idx)) NA_real_ else df$time[onset_idx]
  }

  # ----------------------------------------------------------------
  # Trajectory + derivatives
  # ----------------------------------------------------------------

  baseline <- traj_mat[, 1]
  traj0 <- sweep(traj_mat, 1, baseline, "-")

  dt <- diff(tt)
  if (any(dt <= 0)) stop("Non-increasing time grid in tt.")

  dtraj <- sweep(
    traj_mat[, -1, drop = FALSE] - traj_mat[, -ncol(traj_mat), drop = FALSE],
    2, dt, "/"
  )
  t_dtraj <- (tt[-1] + tt[-length(tt)]) / 2

  ddt <- diff(t_dtraj)
  ddtraj <- sweep(
    dtraj[, -1, drop = FALSE] - dtraj[, -ncol(dtraj), drop = FALSE],
    2, ddt, "/"
  )

  # ----------------------------------------------------------------
  # Peak descriptors
  # ----------------------------------------------------------------

    peak_idx <- max.col(abs(traj0), ties.method = "first")
    t_peak <- tt[peak_idx]
    y_peak <- traj_mat[cbind(seq_len(nrow(traj_mat)), peak_idx)]
    
    peak_dexp <- expm1(y_peak) - expm1(baseline)
    
    peak_log2fc <- log2(
      (expm1(y_peak) + pseudocount) /
      (expm1(baseline) + pseudocount)
    )

  # ----------------------------------------------------------------
  # Onset + AUC + OWP
  # ----------------------------------------------------------------

  t_onset <- vapply(seq_len(nrow(traj_mat)), function(i)
    calc_onset_zdiff(
      as.numeric(traj_mat[i, ]), tt,
      N_minutes, z_diff_threshold, require_consecutive, onset_mode
    ), numeric(1))

  auc_norm <- vapply(seq_len(nrow(traj_mat)), function(i) {
    y <- as.numeric(traj_mat[i, ])
    ymin <- min(y, na.rm = TRUE); ymax <- max(y, na.rm = TRUE)
    if (!is.finite(ymin) || !is.finite(ymax) || ymax <= ymin) return(0)
    trapz_simple(tt, (y - ymin) / (ymax - ymin)) / diff(range(tt))
  }, numeric(1))

  owp <- vapply(seq_len(nrow(traj_mat)), function(i)
    compute_owp(
      as.numeric(traj_mat[i, ]), tt,
      onset_fraction = owp_onset_fraction,
      baseline = "first"
    ), numeric(1))

  owp_ddt <- vapply(seq_len(nrow(dtraj)), function(i)
    compute_owp(
      as.numeric(dtraj[i, ]), t_dtraj,
      onset_fraction = owp_onset_fraction,
      baseline = "zero"
    ), numeric(1))

  # ----------------------------------------------------------------
  # Shape descriptors
  # ----------------------------------------------------------------

  shape_metrics <- lapply(seq_len(nrow(traj_mat)), function(i)
    calc_shape_metrics(as.numeric(traj_mat[i, ]), tt))

  mean_curv <- vapply(shape_metrics, `[[`, numeric(1), "mean_curv")
  max_curv  <- vapply(shape_metrics, `[[`, numeric(1), "max_curv")
  mean_turn <- vapply(shape_metrics, `[[`, numeric(1), "mean_turn")
  max_turn  <- vapply(shape_metrics, `[[`, numeric(1), "max_turn")

  data.frame(
    gene_celltype = rownames(traj_mat),
    t_onset = t_onset,
    t_peak = t_peak,
    baseline_exp = baseline,
    peak_exp = y_peak,
    peak_dlogexp = y_peak - baseline,
    peak_dexp = peak_dexp,              # NEW
    peak_log2fc = peak_log2fc,
    auc_norm = auc_norm,
    owp = owp,
    owp_ddt = owp_ddt,
    mean_curv = mean_curv,
    max_curv = max_curv,
    mean_turn = mean_turn,
    max_turn = max_turn,
    stringsAsFactors = FALSE
  )
}

# ============================================================
# 2. Plot clustering metrics
# ============================================================

plot_trajectory_clustering_metrics <- function(
  x,
  metrics_to_plot = c(
    "WSS",
    "p90_centroid_distance",
    "resolution",
    "mean_silhouette"
  ),
  overlay = FALSE,
  show_k = c("resolution", "elbow", "silhouette"),
  point_size = 1.8,
  line_width = 0.6
) {
  metrics <- x$metrics
  sel <- x$selection

  metric_labels <- c(
    WSS = "Within-cluster sum of squares",
    p90_centroid_distance = "90th percentile distance to centroid",
    resolution = "Fraction of achievable WSS improvement",
    mean_silhouette = "Mean silhouette"
  )

  df_long <- metrics %>%
    dplyr::select(K, dplyr::all_of(metrics_to_plot)) %>%
    tidyr::pivot_longer(-K, names_to = "metric", values_to = "value") %>%
    dplyr::mutate(
      metric = factor(
        metric,
        levels = metrics_to_plot,
        labels = metric_labels[metrics_to_plot]
      )
    )

  selection_df <- tibble::tibble(
    selection = c("resolution", "elbow", "silhouette"),
    K = c(sel$k_resolution, sel$k_elbow, sel$k_silhouette)
  ) %>%
    dplyr::filter(selection %in% show_k, !is.na(K)) %>%
    dplyr::mutate(
      selection = dplyr::recode(
        selection,
        resolution = paste0(round(100 * sel$wss_fraction), "% WSS resolution"),
        elbow = "Elbow",
        silhouette = "Maximum silhouette"
      )
    )

  if (!overlay) {
    return(
      ggplot2::ggplot(df_long, ggplot2::aes(K, value)) +
        ggplot2::geom_line(linewidth = line_width) +
        ggplot2::geom_point(size = point_size, shape = 1) +
        ggplot2::geom_vline(
          data = selection_df,
          ggplot2::aes(xintercept = K, linetype = selection),
          linewidth = 0.7
        ) +
        ggplot2::facet_wrap(~ metric, scales = "free_y", ncol = 1) +
        ggplot2::labs(x = "K", y = NULL, linetype = "K selection") +
        ggplot2::theme_bw() +
        ggplot2::theme(
          strip.background = ggplot2::element_blank(),
          strip.text = ggplot2::element_text(face = "bold"),
          legend.position = "top"
        )
    )
  }

  df_overlay <- df_long %>%
    dplyr::group_by(metric) %>%
    dplyr::mutate(
      value_scaled = {
        rng <- range(value, na.rm = TRUE)
        if (diff(rng) == 0) rep(0, dplyr::n()) else (value - rng[1]) / diff(rng)
      }
    ) %>%
    dplyr::ungroup()

  ggplot2::ggplot(
    df_overlay,
    ggplot2::aes(K, value_scaled, linetype = metric)
  ) +
    ggplot2::geom_line(linewidth = line_width) +
    ggplot2::geom_point(size = point_size, shape = 1) +
    ggplot2::geom_vline(
      data = selection_df,
      ggplot2::aes(xintercept = K),
      inherit.aes = FALSE,
      linewidth = 0.6,
      linetype = 2
    ) +
    ggplot2::labs(
      x = "K",
      y = "Metric scaled to [0,1]",
      linetype = "Metric"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "top")
}

# ============================================================
# 1. Calculate clustering metrics across K
# ============================================================

calculate_trajectory_clustering_metrics <- function(
  X_mat,
  D,
  hc = NULL,
  kmeans_fits = NULL,
  K_grid = 3:200,
  wss_fraction = 0.95,
  select_k_by = c("resolution", "elbow", "silhouette")
) {
  select_k_by <- match.arg(select_k_by)

  X_mat <- as.matrix(X_mat)
  n <- nrow(X_mat)

  if (!inherits(D, "dist")) stop("D must be a 'dist' object.")
  if (attr(D, "Size") != n) stop("X_mat and D contain different numbers of observations.")

  # Provide either hierarchical clustering OR k-means fits
  if (is.null(hc) == is.null(kmeans_fits)) {
    stop("Provide exactly one of 'hc' or 'kmeans_fits'.")
  }

  if (!is.null(hc)) {
    if (!inherits(hc, "hclust")) stop("hc must be an 'hclust' object.")
    if (length(hc$order) != n) stop("X_mat and hc contain different numbers of observations.")
  }

  if (!is.null(kmeans_fits)) {
    if (is.null(names(kmeans_fits))) {
      stop("kmeans_fits must be named by K.")
    }

    # If K_grid wasn't supplied, use the fitted K values
    if (missing(K_grid)) {
      K_grid <- as.integer(names(kmeans_fits))
    }

    if (!all(as.character(K_grid) %in% names(kmeans_fits))) {
      stop("Some values in K_grid are not present in kmeans_fits.")
    }
  }

  if (max(K_grid) >= n) stop("max(K_grid) must be smaller than the number of genes.")

  silhouette_widths_from_dist <- function(cl, D) {
    cl <- as.factor(cl)
    Dm <- as.matrix(D)
    n <- length(cl)
    out <- numeric(n)
    levels_cl <- levels(cl)
    for (i in seq_len(n)) {
      same <- which(cl == cl[i])
      same <- same[same != i]
      if (!length(same)) {
        out[i] <- 0
        next
      }
      a <- mean(Dm[i, same], na.rm = TRUE)
      other_levels <- setdiff(levels_cl, as.character(cl[i]))
      b_vals <- vapply(other_levels, function(k) {
        idx <- which(cl == k)
        if (!length(idx)) return(Inf)
        mean(Dm[i, idx], na.rm = TRUE)
      }, numeric(1))
      b <- min(b_vals, na.rm = TRUE)
      denom <- max(a, b)
      out[i] <- if (is.finite(denom) && denom > 0) (b - a) / denom else 0
    }
    out
  }

  calc_for_k <- function(K) {

    if (!is.null(hc)) {
      cl <- cutree(hc, k = K)

      centroids <- rowsum(X_mat, group = cl) / as.vector(table(cl))
      gene_centroids <- centroids[as.character(cl), , drop = FALSE]

    } else {
      fit <- kmeans_fits[[as.character(K)]]
      if (is.null(fit)) stop("kmeans_fits[['", K, "']] is NULL.")

      cl <- fit$cluster
      if (length(cl) != n) stop("K-means fit and X_mat contain different numbers of observations.")
      if (ncol(fit$centers) != ncol(X_mat)) stop("K-means fit and X_mat have different feature dimensions.")

      gene_centroids <- fit$centers[cl, , drop = FALSE]
    }

    d_to_centroid <- sqrt(rowSums((X_mat - gene_centroids)^2))
    sil_width <- silhouette_widths_from_dist(cl, D)

    data.frame(
      K = K,
      WSS = sum(d_to_centroid^2),
      p90_centroid_distance = unname(quantile(d_to_centroid, 0.90, na.rm = TRUE)),
      mean_silhouette = mean(sil_width, na.rm = TRUE)
    )
  }

  metrics <- dplyr::bind_rows(lapply(K_grid, calc_for_k))

  WSS_minK <- metrics$WSS[1]
  WSS_maxK <- metrics$WSS[nrow(metrics)]

  metrics <- dplyr::mutate(
    metrics,
      resolution = (WSS_minK - WSS) / (WSS_minK - WSS_maxK),
      k_norm = (K - min(K)) / (max(K) - min(K)),
      elbow_score = resolution - k_norm
  )

  eligible <- which(metrics$resolution >= wss_fraction)

  k_resolution <- if (length(eligible)) metrics$K[min(eligible)] else NA_integer_
  k_elbow <- metrics$K[which.max(metrics$elbow_score)]
  k_silhouette <- metrics$K[which.max(metrics$mean_silhouette)]

  k_selected <- switch(
    select_k_by,
    resolution = k_resolution,
    elbow = k_elbow,
    silhouette = k_silhouette
  )

  out <- list(
    metrics = metrics,
    selection = list(
      select_k_by = select_k_by,
      k_selected = k_selected,
      k_resolution = k_resolution,
      k_elbow = k_elbow,
      k_silhouette = k_silhouette,
      max_mean_silhouette = max(metrics$mean_silhouette, na.rm = TRUE),
      wss_fraction = wss_fraction
    )
  )

  class(out) <- c("trajectory_clustering_metrics", class(out))
  out
}


chronopack_dimred_config <- function(root_puro, root_k562a, root_lps) {
  list(
    K562_puro = list(
      dir = file.path(root_puro, "intmd", "sc-K562-1e6-to-tradeSeq.fixed"),
      asso = c("asso_df.1trajectory.nknots4.csv"),
      sig = c("sig_df.1trajectory.nknots4.csv"),
      gamList = c("gamList.1trajectory.nknots8.rds"),
      sce = c("sce.1trajectory.nknots8.rds")
    ),
    K562a = list(
      dir = file.path(root_k562a, "intmd", "sc-K562-1e4_to_tradeSeq.fixed"),
      asso = c("asso_df.1trajectory.nknots4.csv", "asso_df.1trajectory.csv"),
      sig = c("sig_df.1trajectory.nknots4.csv", "sig_df.1trajectory.csv"),
      gamList = c("gamList.1trajectory.nknots8.rds", "gamList.1trajectory.csv", "gamList.1trajectory.rds"),
      sce = c("sce.1trajectory.nknots8.rds", "sce.1trajectory.nknots8.csv")
    ),
    K562b = list(
      dir = file.path(root_lps, "intmd", "sc-K562_to_tradeSeq.fixed"),
      asso = c("asso_df.1trajectory.nknots4.csv", "asso_df.k562.1trajectory.nknots4.expanded.csv"),
      sig = c("sig_df.1trajectory.nknots4.csv", "sig_df.k562.1trajectory.nknots4.expanded.csv"),
      gamList = c("gamList.1trajectory.nknots8.rds", "gamList.knots8.rds"),
      sce = c("sce.1trajectory.nknots8.rds", "sce.knots8.rds")
    ),
    THP1 = list(
      dir = file.path(root_lps, "intmd", "sc-THP1_to_tradeSeq.fixed"),
      asso = c("asso_df.1trajectory.nknots4.csv", "asso_df.thp1.1trajectory.nknots4.expanded.csv"),
      sig = c("sig_df.1trajectory.nknots4.csv", "sig_df.thp1.1trajectory.nknots4.expanded.csv"),
      gamList = c("gamList.1trajectory.nknots8.rds", "gamList.knots8.rds"),
      sce = c("sce.1trajectory.nknots8.rds", "sce.knots8.rds")
    )
  )
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

  list(
    beta = beta,
    Sigma = Sigma,
    dm = cd_ts$dm,
    X = cd_ts$X
  )
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
  out$df_pred <- I(list(df_pred))
  out
}

get_traj_from_tradeSeq_sce <- function(sce, features, nPoints = 200, lineage_pick = 1,
                                       sizeFactors = NULL, off_ref = NULL,
                                       return_log1p = TRUE, verbose = FALSE) {
  features <- intersect(as.character(features), rownames(sce))
  if (!length(features)) {
    stop("No requested features are present in the supplied SCE.")
  }

  rd_ts <- S4Vectors::metadata(SummarizedExperiment::rowData(sce))$tradeSeq
  if (is.null(rd_ts)) {
    rd_ts <- SummarizedExperiment::rowData(sce)$tradeSeq
  }
  cd_ts <- SummarizedExperiment::colData(sce)$tradeSeq
  if (is.null(rd_ts) || is.null(cd_ts)) {
    stop("SCE does not contain tradeSeq rowData/colData metadata.")
  }
  beta_tbl <- rd_ts$beta[[1]]
  if (is.null(beta_tbl)) {
    stop("Could not locate tradeSeq beta matrix in rowData(sce)$tradeSeq$beta[[1]].")
  }
  beta_all <- as.matrix(beta_tbl)
  beta_idx <- match(features, rownames(sce))
  if (!is.null(rownames(beta_all)) && all(features %in% rownames(beta_all))) {
    beta_idx <- match(features, rownames(beta_all))
  }
  if (any(is.na(beta_idx))) {
    stop("Could not match all requested features to SCE beta rows.")
  }
  beta_mat <- beta_all[beta_idx, , drop = FALSE]
  rownames(beta_mat) <- features

  dm <- cd_ts$dm
  X <- cd_ts$X
  if (is.null(dm) || is.null(X)) {
    stop("SCE is missing tradeSeq design matrix metadata.")
  }
  t_cols <- grep("^t[0-9]+$", colnames(dm), value = TRUE)
  avail_lineages <- as.integer(sub("^t", "", t_cols))
  if (!length(avail_lineages)) {
    stop("Could not find t1/t2/... columns in tradeSeq design metadata.")
  }
  ell <- as.integer(lineage_pick)
  if (!(ell %in% avail_lineages)) ell <- avail_lineages[1]
  df_pred <- tradeSeq:::.getPredictRangeDf(dm, lineageId = ell, nPoints = nPoints)
  slingshot_coldata <- SummarizedExperiment::colData(sce)$crv
  if (!is.null(slingshot_coldata)) {
    pseudotime <- slingshot_coldata[, grep(x = colnames(slingshot_coldata), pattern = "pseudotime"), drop = FALSE]
  } else {
    pseudotime <- as.matrix(dm[, t_cols, drop = FALSE])
  }
  X_pred <- tradeSeq:::predictGAM(lpmatrix = X, df = df_pred, pseudotime = pseudotime)
  if (!is.null(colnames(X_pred)) && !is.null(colnames(beta_mat)) && all(colnames(X_pred) %in% colnames(beta_mat))) {
    beta_mat <- beta_mat[, colnames(X_pred), drop = FALSE]
  }
  if (ncol(beta_mat) != ncol(X_pred)) {
    stop("Beta matrix column count does not match prediction design columns.")
  }

  off_ref_use <- compute_tradeSeq_reference_offset(sizeFactors, off_ref = off_ref)
  eta_mat <- beta_mat %*% t(X_pred)
  eta_mat <- sweep(eta_mat, 2, rep(as.numeric(off_ref_use)[1], ncol(eta_mat)), "+")
  traj_mat <- if (isTRUE(return_log1p)) log1p(exp(eta_mat)) else exp(eta_mat)
  rownames(traj_mat) <- features
  tt <- as.numeric(df_pred[[paste0("t", ell)]])
  colnames(traj_mat) <- as.character(tt)
  if (isTRUE(verbose)) {
    message("Extracted ", nrow(traj_mat), " SCE-backed fitted trajectories.")
  }
  list(traj_mat = traj_mat, tt = tt)
}

first_existing <- function(...) {
  paths <- c(...)
  hit <- paths[file.exists(paths)]
  if (!length(hit)) {
    stop("None of these files exist: ", paste(paths, collapse = ", "))
  }
  hit[[1]]
}

resolve_output <- function(dir, candidates) {
  first_existing(file.path(dir, candidates))
}

read_counts_matrix <- function(path) {
  counts <- as.matrix(read.csv(gzfile(path), row.names = 1, check.names = FALSE))
  storage.mode(counts) <- "integer"
  counts
}

make_pt_cw_sf <- function(counts, dir) {
  cell_times_path <- first_existing(
    file.path(dir, "cell_times.csv.gz"),
    file.path(dir, "cell_times..csv.gz"),
    file.path(dir, "cell_times.full.csv.gz")
  )
  ct <- read.csv(gzfile(cell_times_path), stringsAsFactors = FALSE)
  libs_full <- read.csv(file.path(dir, "cell_libsize_full.csv"), stringsAsFactors = FALSE)

  stopifnot(all(colnames(counts) %in% ct$cell_id))
  stopifnot(all(colnames(counts) %in% libs_full$cell_id))

  time_vec <- setNames(ct$time, ct$cell_id)[colnames(counts)]
  lib_full <- setNames(libs_full$libsize_full, libs_full$cell_id)[colnames(counts)]

  sf <- as.numeric(lib_full)
  sf <- sf / mean(sf)

  pt01 <- (time_vec - min(time_vec, na.rm = TRUE)) / diff(range(time_vec, na.rm = TRUE))
  PT <- matrix(pt01, ncol = 1, dimnames = list(NULL, "lineage1"))
  CW <- matrix(1, nrow = ncol(counts), ncol = 1, dimnames = list(NULL, "lineage1"))

  list(PT = PT, CW = CW, sf = sf, time_vec = time_vec, lib_full = lib_full)
}

load_dimred_tradeSeq_dataset <- function(entry, load_gamList = FALSE) {
  dir <- normalizePath(entry$dir, mustWork = TRUE)

  asso_df <- read.csv(resolve_output(dir, entry$asso), row.names = 1, check.names = FALSE)
  asso_df <- asso_df[!is.na(asso_df$waldStat_1), , drop = FALSE]
  sig_df <- read.csv(resolve_output(dir, entry$sig), row.names = 1, check.names = FALSE)
  gamList <- NULL
  if (isTRUE(load_gamList)) {
    gamList_path <- tryCatch(resolve_output(dir, entry$gamList), error = function(e) NULL)
    if (!is.null(gamList_path)) {
      gamList <- readRDS(gamList_path)
    }
  }
  sce <- readRDS(resolve_output(dir, entry$sce))

  counts <- read_counts_matrix(first_existing(
    file.path(dir, "counts_genes_x_cells.csv.gz"),
    file.path(dir, "counts_genes_x_cells..csv.gz")
  ))
  if (is.null(rownames(sce)) && nrow(sce) == nrow(counts)) {
    rownames(sce) <- rownames(counts)
  }
  meta <- make_pt_cw_sf(counts, dir)

  list(
    asso_df = asso_df,
    sig_df = sig_df,
    sce = sce,
    gamList = gamList,
    counts_sub = counts,
    PT = meta$PT,
    CW = meta$CW,
    sf = meta$sf,
    time_vec = meta$time_vec,
    lib_full = meta$lib_full
  )
}

load_dimred_fixed_tradeSeq_outputs <- function(cfg, load_gamList = FALSE) {
  out <- lapply(names(cfg), function(key) {
    message("Loading fixed tradeSeq objects: ", key)
    load_dimred_tradeSeq_dataset(cfg[[key]], load_gamList = load_gamList)
  })
  names(out) <- names(cfg)

  list(
    asso_df.ls = lapply(out, `[[`, "asso_df"),
    sig_df.ls = lapply(out, `[[`, "sig_df"),
    sce.ls = lapply(out, `[[`, "sce"),
    gamList.ls = lapply(out, `[[`, "gamList"),
    counts_sub.ls = lapply(out, `[[`, "counts_sub"),
    PT.ls = lapply(out, `[[`, "PT"),
    CW.ls = lapply(out, `[[`, "CW"),
    sf.ls = lapply(out, `[[`, "sf")
  )
}

select_dimred_features <- function(asso_df.ls, sce.ls, counts_sub.ls = NULL, keys, pvalue_threshold = 0.05) {
  features.ls <- list()
  avail_genes.ls <- list()
  usable_features.ls <- list()

  for (key in keys) {
    tmp <- asso_df.ls[[key]]
    avail <- character(0)
    if (!is.null(counts_sub.ls) && !is.null(counts_sub.ls[[key]])) {
      avail <- rownames(counts_sub.ls[[key]])
    }
    sce_names <- rownames(sce.ls[[key]])
    if (!is.null(sce_names) && length(sce_names)) {
      avail <- if (length(avail)) intersect(avail, sce_names) else sce_names
    }
    avail_genes.ls[[key]] <- avail
    tmp <- tmp[tmp$pvalue_1 < pvalue_threshold, , drop = FALSE]
    features.ls[[key]] <- rownames(tmp)
    usable_features.ls[[key]] <- intersect(features.ls[[key]], avail_genes.ls[[key]])
  }
  names(features.ls) <- keys
  names(avail_genes.ls) <- keys
  names(usable_features.ls) <- keys

  features <- Reduce(union, features.ls[keys])
  avail_genes <- Reduce(union, avail_genes.ls[keys])
  usable_features <- Reduce(union, usable_features.ls[keys])

  list(
    features = features,
    avail_genes = avail_genes,
    usable_features = usable_features,
    features.ls = features.ls,
    avail_genes.ls = avail_genes.ls,
    usable_features.ls = usable_features.ls
  )
}

extract_dimred_smooth_trajectories <- function(sce.ls, sf.ls, PT.ls, features,
                                                keys, t_interval_map, lineage_pick = 1) {
  traj.ls <- list()

  for (key in keys) {
    features_key <- if (is.list(features)) features[[key]] else features
    if (is.null(features_key)) {
      stop("No feature vector supplied for key: ", key)
    }
    nPoints <- 1 + ((length(unique(PT.ls[[key]])) - 1) * t_interval_map[[key]])
    tmp <- get_traj_from_tradeSeq_sce(
      sce.ls[[key]], features_key, nPoints = nPoints, lineage_pick = lineage_pick,
      sizeFactors = sf.ls[[key]]
    )
    traj.ls[[key]] <- tmp
    rm(tmp)
  }

  traj.ls
}

.interpolate_dimred_rows <- function(mat, x, xout) {
  stopifnot(is.matrix(mat), length(x) == ncol(mat), length(xout) >= 1)

  ord <- order(x)
  x <- x[ord]
  mat <- mat[, ord, drop = FALSE]

  if (anyDuplicated(x)) {
    stop("Duplicated time points encountered after rescaling; choose a different scale_factor or collapse duplicates first.")
  }

  out <- t(vapply(
    seq_len(nrow(mat)),
    function(i) {
      stats::approx(
        x = x,
        y = mat[i, ],
        xout = xout,
        method = "linear",
        rule = 1,
        ties = "ordered"
      )$y
    },
    numeric(length(xout))
  ))

  rownames(out) <- rownames(mat)
  colnames(out) <- as.character(xout)
  out
}

align_rescale_dimred_trajectories <- function(traj.ls, time_map, offset_map, scale_map = NULL,
                                               default_scale_overrides = c(THP1 = 0.65)) {
  keys_use <- intersect(intersect(names(traj.ls), names(time_map)), names(offset_map))
  stopifnot(length(keys_use) >= 2)

  scale_map_use <- stats::setNames(rep(1, length(keys_use)), keys_use)
  for (nm in intersect(names(default_scale_overrides), names(scale_map_use))) {
    scale_map_use[[nm]] <- as.numeric(default_scale_overrides[[nm]])
  }
  if (!is.null(scale_map)) {
    nm <- intersect(names(scale_map_use), names(scale_map))
    scale_map_use[nm] <- as.numeric(unlist(scale_map[nm]))
  }
  stopifnot(all(is.finite(scale_map_use)), all(scale_map_use > 0))

  traj_rescaled <- list()
  t_aligned_by_key <- list()

  for (k in keys_use) {
    obj <- traj.ls[[k]]
    mat <- obj$traj_mat
    tt <- as.numeric(obj$tt)
    stopifnot(is.matrix(mat), length(tt) == ncol(mat))

    TMAX <- as.numeric(time_map[[k]])
    off <- as.numeric(offset_map[[k]])
    sf <- as.numeric(scale_map_use[[k]])
    t_real <- round(tt * TMAX, 6)
    t_aligned <- (t_real * sf) + off

    keep <- which(is.finite(t_aligned) & t_aligned >= 0)
    if (!length(keep)) {
      stop(sprintf("[%s] scale_factor=%g and offset=%g remove all timepoints.", k, sf, off))
    }

    mat_obs <- mat[, keep, drop = FALSE]
    t_obs <- t_aligned[keep]
    t_grid <- seq(from = max(1, ceiling(min(t_obs) - 1e-9)), to = floor(max(t_obs) + 1e-9), by = 1)
    if (!length(t_grid)) {
      stop(sprintf("[%s] no positive integer-minute grid remains after rescaling/offset.", k))
    }

    traj_rescaled[[k]] <- .interpolate_dimred_rows(mat_obs, x = t_obs, xout = t_grid)
    t_aligned_by_key[[k]] <- t_grid

    cat(sprintf(
      "[%s] TMAX=%g, scale=%g, offset=%g => aligned range %.2f..%.2f, resampled to %d minute points\n",
      k, TMAX, sf, off, min(t_obs), max(t_obs), length(t_grid)
    ))
  }

  t_end <- floor(min(vapply(t_aligned_by_key, max, numeric(1), na.rm = TRUE)) + 1e-9)
  cat("Shortest aligned end time across keys after rescaling:", t_end, "min\n")

  traj_trunc <- lapply(keys_use, function(k) {
    m <- traj_rescaled[[k]]
    t <- as.numeric(colnames(m))
    keep <- which(is.finite(t) & t >= 1 & t <= t_end)
    m2 <- m[, keep, drop = FALSE]
    colnames(m2) <- as.character(round(as.numeric(colnames(m2))))
    m2
  })
  names(traj_trunc) <- keys_use

  shared_minutes <- Reduce(intersect, lapply(traj_trunc, colnames))
  shared_minutes <- as.character(sort(as.numeric(shared_minutes)))
  if (length(shared_minutes) == 0) {
    stop("No shared aligned minute grid remains across trajectories after rescaling.")
  }

  traj_trunc_rescaled.ls <- lapply(traj_trunc, function(m) m[, shared_minutes, drop = FALSE])
  cat("Final shared aligned minute range:", min(as.numeric(shared_minutes)), "..", max(as.numeric(shared_minutes)),
      " (n =", length(shared_minutes), ")\n")

  list(
    traj_trunc_rescaled.ls = traj_trunc_rescaled.ls,
    scale_map_use = scale_map_use,
    shared_minutes = shared_minutes,
    keys_use = keys_use
  )
}

construct_dimred_traj_matrix <- function(traj) {
  t_max <- dim(traj[[1]])[2]
  tt <- seq_len(t_max) / t_max

  cn0 <- colnames(traj[[1]])
  stopifnot(all(vapply(traj, function(m) is.matrix(m) && identical(colnames(m), cn0), logical(1))))

  traj_all <- do.call(rbind, Map(function(mat, key) {
    rownames(mat) <- paste0(rownames(mat), "_", key)
    mat
  }, traj, names(traj)))

  traj_mat <- as.matrix(traj_all)
  storage.mode(traj_mat) <- "double"
  message("Valid trajectories: ", nrow(traj_mat), " genes x ", ncol(traj_mat), " timepoints.")

  list(tt = tt, traj_mat = traj_mat)
}

z_row <- function(M) {
  mu <- rowMeans(M, na.rm = TRUE)
  sd <- matrixStats::rowSds(M, na.rm = TRUE)
  sd[!is.finite(sd) | sd == 0] <- 1
  sweep(sweep(M, 1, mu, "-"), 1, sd, "/")
}

build_dimred_pc_df <- function(traj_mat, tt, asso_df.ls, scale_map_use,
                               descriptor_args = list(
                                 N_minutes = 10 / length(tt),
                                 z_diff_threshold = 0.2,
                                 require_consecutive = 30,
                                 onset_mode = "horizontal",
                                 pseudocount = 0
                               )) {
  desc <- do.call(get_descriptors, c(list(traj_mat = traj_mat, tt = tt), descriptor_args))
  gene_celltype <- rownames(traj_mat)
  descriptor_key <- sub("^[^_]+_", "", gene_celltype)
  descriptor_scale <- length(tt) / as.numeric(scale_map_use[descriptor_key])
  descriptor_scale[!is.finite(descriptor_scale)] <- length(tt)
  desc$t_onset <- desc$t_onset * descriptor_scale
  desc$t_peak <- desc$t_peak * descriptor_scale

  pc_df <- data.frame(
    gene_celltype = gene_celltype,
    gene = sapply(gene_celltype, function(x) strsplit(x, "_")[[1]][1], USE.NAMES = FALSE),
    cell_type = sub("^[^_]+_", "", gene_celltype)
  ) |>
    dplyr::left_join(desc, by = c("gene_celltype"))
  rownames(pc_df) <- pc_df$gene_celltype
  pc_df <- dplyr::select(pc_df, -gene_celltype)

  for (nm in c("waldStat_1", "meanLogFC", "pvalue_1", "qval_1")) {
    pc_df[[nm]] <- NA_real_
  }

  lookup_df <- do.call(rbind, lapply(names(asso_df.ls), function(ct) {
    df <- asso_df.ls[[ct]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    out <- data.frame(cell_type = ct, gene = rownames(df), stringsAsFactors = FALSE, row.names = NULL)
    for (nm in c("waldStat_1", "meanLogFC", "pvalue_1", "qval_1")) {
      out[[paste0(nm, "_lookup")]] <- if (nm %in% colnames(df)) df[[nm]] else NA_real_
    }
    out
  }))

  pc_df <- pc_df |>
    tibble::rownames_to_column("row_id") |>
    dplyr::left_join(lookup_df, by = c("cell_type", "gene")) |>
    dplyr::mutate(
      waldStat_1 = ifelse(is.na(waldStat_1), waldStat_1_lookup, waldStat_1),
      meanLogFC = ifelse(is.na(meanLogFC), meanLogFC_lookup, meanLogFC),
      pvalue_1 = ifelse(is.na(pvalue_1), pvalue_1_lookup, pvalue_1),
      qval_1 = ifelse(is.na(qval_1), qval_1_lookup, qval_1)
    ) |>
    dplyr::select(-waldStat_1_lookup, -meanLogFC_lookup, -pvalue_1_lookup, -qval_1_lookup) |>
    tibble::column_to_rownames("row_id")

  pc_df
}

prepare_dimred_matrices <- function(traj_mat, tt, pc_df, significance_col = "pvalue_1", significance_threshold = 0.05) {
  X <- z_row(traj_mat)
  dt <- diff(tt)
  Xd <- traj_mat[, -1, drop = FALSE] - traj_mat[, -ncol(traj_mat), drop = FALSE]
  Xd <- sweep(Xd, 2, dt, "/")
  Xd <- z_row(Xd)
  stopifnot(nrow(X) == nrow(Xd))

  significant_rows <- rownames(pc_df[!is.na(pc_df[[significance_col]]) & pc_df[[significance_col]] < significance_threshold, , drop = FALSE])
  X_mat <- X[significant_rows, , drop = FALSE]
  Xd_mat <- Xd[significant_rows, , drop = FALSE]
  stopifnot(identical(rownames(X_mat), rownames(Xd_mat)))

  list(X = X, Xd = Xd, X_mat = X_mat, Xd_mat = Xd_mat, significant_rows = significant_rows)
}

compute_dimred_pca <- function(X_mat, Xd_mat, variance_threshold = 0.99) {
  pca_X <- prcomp(X_mat, center = TRUE, scale. = FALSE)
  cumvar <- cumsum(pca_X$sdev^2) / sum(pca_X$sdev^2)
  n_pc_X <- which(cumvar >= variance_threshold)[1]
  X_pc <- pca_X$x[, seq_len(n_pc_X), drop = FALSE]

  pca_Xd <- prcomp(Xd_mat, center = TRUE, scale. = FALSE)
  cumvar_d <- cumsum(pca_Xd$sdev^2) / sum(pca_Xd$sdev^2)
  n_pc_Xd <- which(cumvar_d >= variance_threshold)[1]
  Xd_pc <- pca_Xd$x[, seq_len(n_pc_Xd), drop = FALSE]

  list(pca_X = pca_X, pca_Xd = pca_Xd, X_pc = X_pc, Xd_pc = Xd_pc, n_pc_X = n_pc_X, n_pc_Xd = n_pc_Xd)
}

add_dimred_pca_columns <- function(pc_df, pca_X, n_columns = 4) {
  pc_df <- pc_df[rownames(pca_X$x), , drop = FALSE]
  for (i in seq_len(n_columns)) {
    pc_df[[paste0("PC", i)]] <- pca_X$x[, i]
  }
  pc_df
}

scale_by_total_column_variance <- function(M) {
  M <- as.matrix(M)
  scale_factor <- sqrt(sum(apply(M, 2, var)))
  if (!is.finite(scale_factor) || scale_factor == 0) {
    stop("Cannot scale matrix: total column variance is zero or non-finite.")
  }
  M / scale_factor
}

prepare_dimred_clustering_inputs <- function(X_mat, Xd_mat) {
  list(
    X_n = scale_by_total_column_variance(X_mat),
    Xd_n = scale_by_total_column_variance(Xd_mat)
  )
}

cluster_centroid_distances <- function(X_use, cl) {
  gene_cluster_distance <- rep(NA_real_, nrow(X_use))
  names(gene_cluster_distance) <- rownames(X_use)
  for (k in sort(unique(cl))) {
    idx <- which(cl == k)
    Xi <- X_use[idx, , drop = FALSE]
    centroid <- colMeans(Xi)
    gene_cluster_distance[idx] <- sqrt(rowSums(sweep(Xi, 2, centroid, "-")^2))
  }
  data.frame(gene = rownames(X_use), cluster = cl, centroid_dist = gene_cluster_distance, stringsAsFactors = FALSE)
}

run_dimred_signed_clustering <- function(pc_df, X_n, Xd_n, sign = c("induced", "repressed"),
                                         alpha = 0.25, log2fc_threshold = 0,
                                         K_grid = 2:100, wss_fraction = 0.95,
                                         select_k_by = "elbow") {
  sign <- match.arg(sign)
  if (sign == "induced") {
    subset.rownames <- rownames(pc_df[pc_df$peak_log2fc > log2fc_threshold, , drop = FALSE])
  } else {
    subset.rownames <- rownames(pc_df[pc_df$peak_log2fc < -log2fc_threshold, , drop = FALSE])
  }

  X_joint <- cbind(sqrt(1 - alpha) * X_n, sqrt(alpha) * Xd_n)[subset.rownames, , drop = FALSE]
  D_joint <- dist(X_joint, method = "euclidean")
  hc <- hclust(D_joint, method = "ward.D2")

  clust_eval <- calculate_trajectory_clustering_metrics(
    X_mat = X_joint, D = D_joint, hc = hc, K_grid = K_grid,
    wss_fraction = wss_fraction, select_k_by = select_k_by
  )

  broad_cl <- cutree(hc, k = clust_eval$selection$k_silhouette)
  fine_cl <- cutree(hc, k = clust_eval$selection$k_elbow)

  list(
    X_use = X_joint,
    D_use = D_joint,
    hc = hc,
    clust_eval = clust_eval,
    broad = broad_cl,
    fine = fine_cl,
    broad_df = cluster_centroid_distances(X_joint, broad_cl),
    fine_df = cluster_centroid_distances(X_joint, fine_cl)
  )
}

add_dimred_cluster_assignment <- function(pc_df, cl, disp_df, name) {
  cluster_col <- paste0("cluster_", name)
  dist_col <- paste0("cluster_dist_", name)
  pc_df[[cluster_col]] <- as.factor(cl[rownames(pc_df)])
  pc_df[[dist_col]] <- disp_df$centroid_dist[match(rownames(pc_df), disp_df$gene)]
  pc_df[[cluster_col]] <- factor(pc_df[[cluster_col]], levels = c(levels(pc_df[[cluster_col]]), "-1"))
  pc_df[[cluster_col]][is.na(pc_df[[cluster_col]])] <- "-1"
  pc_df[[cluster_col]] <- factor(pc_df[[cluster_col]])
  pc_df
}

add_dimred_signed_clusters <- function(pc_df, induced_res, repressed_res) {
  pc_df <- add_dimred_cluster_assignment(pc_df, induced_res$broad, induced_res$broad_df, "broad_induced")
  pc_df <- add_dimred_cluster_assignment(pc_df, induced_res$fine, induced_res$fine_df, "fine_induced")
  pc_df <- add_dimred_cluster_assignment(pc_df, repressed_res$broad, repressed_res$broad_df, "broad_repressed")
  pc_df <- add_dimred_cluster_assignment(pc_df, repressed_res$fine, repressed_res$fine_df, "fine_repressed")
  pc_df
}

rename_dimred_broad_clusters <- function(pc_df) {
  rename_induced <- function(x) {
    x <- as.character(x)
    if (x == "1") return("impulse-up_1")
    if (x == "2") return("trend-up")
    if (x == "3") return("impulse-up_2")
    x
  }
  rename_repressed <- function(x) {
    x <- as.character(x)
    if (x == "1") return("trend-down_1")
    if (x == "2") return("trend-down_2")
    if (x == "3") return("transient-down")
    x
  }
  pc_df$cluster_broad_induced <- sapply(pc_df$cluster_broad_induced, rename_induced, USE.NAMES = FALSE)
  pc_df$cluster_broad_repressed <- sapply(pc_df$cluster_broad_repressed, rename_repressed, USE.NAMES = FALSE)
  pc_df
}

finalize_dimred_pc_df <- function(pc_df, cell_type_order = c("K562a", "K562b", "K562_puro", "THP1")) {
  pc_df <- rename_dimred_broad_clusters(pc_df)
  pc_df$broad_cluster <- factor(ifelse(
    pc_df$cluster_broad_induced != "-1",
    as.character(pc_df$cluster_broad_induced),
    as.character(pc_df$cluster_broad_repressed)
  ))
  pc_df$cluster_dist_broad <- dplyr::coalesce(pc_df$cluster_dist_broad_induced, pc_df$cluster_dist_broad_repressed)
  pc_df$cluster_dist_fine <- dplyr::coalesce(pc_df$cluster_dist_fine_induced, pc_df$cluster_dist_fine_repressed)
  pc_df[order(match(as.character(pc_df$cell_type), cell_type_order)), , drop = FALSE]
}

write_matrix_csv <- function(x, path) {
  write.csv(as.data.frame(x, check.names = FALSE), path, row.names = TRUE)
  invisible(normalizePath(path, mustWork = FALSE))
}

chronopack_characterize_gene_trajectories <- function(
  cfg,
  output_dir,
  keys = c("K562a", "K562_puro", "K562b", "THP1"),
  t_interval_map = c(K562a = 10, K562_puro = 10, K562b = 15, THP1 = 15),
  time_map = c(K562a = 110, K562_puro = 110, K562b = 165, THP1 = 165),
  offset_map = c(K562a = 0, K562_puro = 0, K562b = -68, THP1 = -10),
  scale_map = c(THP1 = 0.65),
  pvalue_threshold = 0.05,
  descriptor_args = NULL,
  lineage_pick = 1,
  write_outputs = TRUE
) {
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  loaded <- load_dimred_fixed_tradeSeq_outputs(cfg)
  feature_info <- select_dimred_features(
    asso_df.ls = loaded$asso_df.ls,
    sce.ls = loaded$sce.ls,
    counts_sub.ls = loaded$counts_sub.ls,
    keys = keys,
    pvalue_threshold = pvalue_threshold
  )

  traj.ls <- extract_dimred_smooth_trajectories(
    sce.ls = loaded$sce.ls,
    sf.ls = loaded$sf.ls,
    PT.ls = loaded$PT.ls,
    features = feature_info$usable_features.ls,
    keys = keys,
    t_interval_map = t_interval_map,
    lineage_pick = lineage_pick
  )

  align_res <- align_rescale_dimred_trajectories(
    traj.ls = traj.ls,
    time_map = time_map,
    offset_map = offset_map,
    scale_map = scale_map
  )

  traj_res <- construct_dimred_traj_matrix(align_res$traj_trunc_rescaled.ls)
  traj_mat <- traj_res$traj_mat
  tt <- traj_res$tt

  if (is.null(descriptor_args)) {
    descriptor_args <- list(
      N_minutes = 10 / length(tt),
      z_diff_threshold = 0.2,
      require_consecutive = 30,
      onset_mode = "horizontal",
      pseudocount = 0
    )
  }

  pc_df <- build_dimred_pc_df(
    traj_mat = traj_mat,
    tt = tt,
    asso_df.ls = loaded$asso_df.ls,
    scale_map_use = align_res$scale_map_use,
    descriptor_args = descriptor_args
  )

  matrix_res <- prepare_dimred_matrices(
    traj_mat = traj_mat,
    tt = tt,
    pc_df = pc_df,
    significance_col = "pvalue_1",
    significance_threshold = pvalue_threshold
  )

  paths <- list()
  if (isTRUE(write_outputs)) {
    paths$traj_mat <- file.path(output_dir, "gene_trajectory_fit.csv")
    paths$time_grid <- file.path(output_dir, "gene_trajectory_time_grid.csv")
    paths$pc_df_characterized <- file.path(output_dir, "pc_df.characterized.csv")
    paths$X <- file.path(output_dir, "gene_trajectory_fit_z.csv")
    paths$Xd <- file.path(output_dir, "gene_trajectory_fit_ddt_z.csv")
    paths$X_mat <- file.path(output_dir, "gene_trajectory_fit_significant_z.csv")
    paths$Xd_mat <- file.path(output_dir, "gene_trajectory_fit_ddt_significant_z.csv")

    write_matrix_csv(traj_mat, paths$traj_mat)
    write.csv(data.frame(time = tt), paths$time_grid, row.names = FALSE)
    write.csv(pc_df, paths$pc_df_characterized)
    write_matrix_csv(matrix_res$X, paths$X)
    write_matrix_csv(matrix_res$Xd, paths$Xd)
    write_matrix_csv(matrix_res$X_mat, paths$X_mat)
    write_matrix_csv(matrix_res$Xd_mat, paths$Xd_mat)
  }

  out <- list(
    loaded = loaded,
    feature_info = feature_info,
    traj = traj.ls,
    align = align_res,
    tt = tt,
    traj_mat = traj_mat,
    pc_df = pc_df,
    matrix_res = matrix_res,
    paths = paths,
    parameters = list(
      keys = keys,
      t_interval_map = t_interval_map,
      time_map = time_map,
      offset_map = offset_map,
      scale_map = scale_map,
      pvalue_threshold = pvalue_threshold,
      descriptor_args = descriptor_args,
      lineage_pick = lineage_pick
    )
  )
  class(out) <- c("chronopack_gene_trajectory_characterization", class(out))
  out
}

chronopack_cluster_gene_trajectories <- function(
  characterization,
  output_dir,
  alpha = 0.25,
  log2fc_threshold = 0,
  K_grid = 2:100,
  wss_fraction = 0.95,
  select_k_by = "elbow",
  variance_threshold = 0.99,
  pca_columns = 4,
  cell_type_order = c("K562a", "K562b", "K562_puro", "THP1"),
  write_outputs = TRUE
) {
  output_dir <- normalizePath(output_dir, mustWork = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  pc_df <- characterization$pc_df
  matrix_res <- characterization$matrix_res
  X_mat <- matrix_res$X_mat
  Xd_mat <- matrix_res$Xd_mat

  pca_res <- compute_dimred_pca(
    X_mat = X_mat,
    Xd_mat = Xd_mat,
    variance_threshold = variance_threshold
  )
  pc_df <- add_dimred_pca_columns(pc_df, pca_res$pca_X, n_columns = pca_columns)

  clustering_inputs <- prepare_dimred_clustering_inputs(X_mat, Xd_mat)
  induced_res <- run_dimred_signed_clustering(
    pc_df = pc_df,
    X_n = clustering_inputs$X_n,
    Xd_n = clustering_inputs$Xd_n,
    sign = "induced",
    alpha = alpha,
    log2fc_threshold = log2fc_threshold,
    K_grid = K_grid,
    wss_fraction = wss_fraction,
    select_k_by = select_k_by
  )
  repressed_res <- run_dimred_signed_clustering(
    pc_df = pc_df,
    X_n = clustering_inputs$X_n,
    Xd_n = clustering_inputs$Xd_n,
    sign = "repressed",
    alpha = alpha,
    log2fc_threshold = log2fc_threshold,
    K_grid = K_grid,
    wss_fraction = wss_fraction,
    select_k_by = select_k_by
  )

  pc_df <- add_dimred_signed_clusters(pc_df, induced_res, repressed_res)
  pc_df <- finalize_dimred_pc_df(pc_df, cell_type_order = cell_type_order)

  paths <- list()
  if (isTRUE(write_outputs)) {
    paths$pc_df <- file.path(output_dir, "pc_df.new.csv")
    paths$pc_df_rds <- file.path(output_dir, "pc_df.new.rds")
    paths$pca_X <- file.path(output_dir, "gene_trajectory_X_pca.csv")
    paths$pca_Xd <- file.path(output_dir, "gene_trajectory_Xd_pca.csv")
    paths$induced_metrics <- file.path(output_dir, "cluster_metrics.induced.csv")
    paths$repressed_metrics <- file.path(output_dir, "cluster_metrics.repressed.csv")

    write.csv(pc_df, paths$pc_df)
    saveRDS(pc_df, paths$pc_df_rds)
    write_matrix_csv(pca_res$pca_X$x, paths$pca_X)
    write_matrix_csv(pca_res$pca_Xd$x, paths$pca_Xd)
    write.csv(induced_res$clust_eval$metrics, paths$induced_metrics, row.names = FALSE)
    write.csv(repressed_res$clust_eval$metrics, paths$repressed_metrics, row.names = FALSE)
  }

  out <- list(
    pc_df = pc_df,
    pca_res = pca_res,
    clustering_inputs = clustering_inputs,
    induced = induced_res,
    repressed = repressed_res,
    paths = paths,
    parameters = list(
      alpha = alpha,
      log2fc_threshold = log2fc_threshold,
      K_grid = K_grid,
      wss_fraction = wss_fraction,
      select_k_by = select_k_by,
      variance_threshold = variance_threshold,
      pca_columns = pca_columns,
      cell_type_order = cell_type_order
    )
  )
  class(out) <- c("chronopack_gene_trajectory_clustering", class(out))
  out
}
