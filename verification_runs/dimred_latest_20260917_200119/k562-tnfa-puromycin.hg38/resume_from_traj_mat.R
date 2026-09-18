
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(matrixStats)
  library(cluster)
  library(SingleCellExperiment)
})
RNGversion("4.3.0")
set.seed(1)
SELF_CONTAINED_OUTDIR <- "/new-stg/home/alex/tmp/ChronoPack/verification_runs/dimred_latest_20260917_200119/output"
setwd(SELF_CONTAINED_OUTDIR)
root_puro <- "/new-stg/home/alex/projects/chrono_seq/analyses/k562-tnfa-puromycin.hg38"
traj_mat <- readRDS(file.path(SELF_CONTAINED_OUTDIR, "traj_mat.rds"))
tt <- readRDS(file.path(SELF_CONTAINED_OUTDIR, "tt.rds"))
scale_map_use <- c(K562a = 1, K562_puro = 1, K562b = 1, THP1 = 0.65)
read_asso <- function(dir, candidates) {
  hit <- file.path(dir, candidates)
  hit <- hit[file.exists(hit)]
  if (!length(hit)) stop("No association file in ", dir)
  df <- read.csv(hit[[1]], row.names = 1, check.names = FALSE)
  df[!is.na(df$waldStat_1), , drop = FALSE]
}
base <- "/new-stg/home/alex/tmp/ChronoPack/verification_runs/filename_migration_20260917_074750"
asso_df.ls <- list(
  K562_puro = read_asso(file.path(base, "k562-tnfa-puromycin.hg38/intmd/sc-K562-1e6-to-tradeSeq.fixed"), c("asso_df.1trajectory.nknots4.csv")),
  K562a = read_asso(file.path(base, "k562-tnfa.hg38.rep/intmd/sc-K562-1e4_to_tradeSeq.fixed"), c("asso_df.1trajectory.nknots4.csv", "asso_df.1trajectory.csv")),
  K562b = read_asso(file.path(base, "thp1-k562-lps.hg38.rep/intmd/sc-K562_to_tradeSeq.fixed"), c("asso_df.1trajectory.nknots4.csv", "asso_df.k562.1trajectory.nknots4.expanded.csv")),
  THP1 = read_asso(file.path(base, "thp1-k562-lps.hg38.rep/intmd/sc-THP1_to_tradeSeq.fixed"), c("asso_df.1trajectory.nknots4.csv", "asso_df.thp1.1trajectory.nknots4.expanded.csv"))
)


# ---- notebook cell 12 ----

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

# ---- notebook cell 29 ----

set.seed(1)

z_row <- function(M) {
  mu <- rowMeans(M, na.rm = TRUE)
  sd <- matrixStats::rowSds(M, na.rm = TRUE)
  sd[!is.finite(sd) | sd == 0] <- 1
  sweep(sweep(M, 1, mu, "-"), 1, sd, "/")
}

sil_to_alpha <- function(s) {
  s <- pmax(-0.2, pmin(0.8, s))
  0.15 + (s + 0.2) / (0.8 + 0.2) * (1.0 - 0.15)
}

desc = get_descriptors(traj_mat, tt, N_minutes = 10/length(tt), z_diff_threshold = 0.2, require_consecutive = 30, onset_mode = 'horizontal', pseudocount=0)
gene_celltype <- rownames(traj_mat)
descriptor_key <- sub("^[^_]+_", "", gene_celltype)
descriptor_scale <- length(tt) / as.numeric(scale_map_use[descriptor_key])
descriptor_scale[!is.finite(descriptor_scale)] <- length(tt)
desc$t_onset <- desc$t_onset * descriptor_scale
desc$t_peak <- desc$t_peak * descriptor_scale

# ---- notebook cell 30 ----

pc_df <- data.frame(gene_celltype = gene_celltype,
                    gene = sapply(gene_celltype, function(x) { strsplit(x, '_')[[1]][1] }, USE.NAMES = F),
                    cell_type = sub("^[^_]+_", "", gene_celltype)) %>% left_join(desc, by = c('gene_celltype'))
rownames(pc_df) = pc_df$gene_celltype
pc_df = pc_df %>% select(-gene_celltype)

for (nm in c("waldStat_1", "meanLogFC", "pvalue_1", "qval_1")) {
    pc_df[[nm]] <- NA_real_
}

# Build lookup table from asso_df.ls
lookup_df <- do.call(rbind, lapply(names(asso_df.ls), function(ct) {
                                df <- asso_df.ls[[ct]]
                                if (is.null(df) || nrow(df) == 0) return(NULL)
                            
                                out <- data.frame(cell_type = ct, gene = rownames(df), stringsAsFactors = FALSE, row.names = NULL)
                            
                                for (nm in c("waldStat_1", "meanLogFC", "pvalue_1", "qval_1")) {
                                  out[[paste0(nm, "_lookup")]] <- if (nm %in% colnames(df)) df[[nm]] else NA_real_
                                }
                        
                                out
                              })
                    )

pc_df <- pc_df %>% rownames_to_column("row_id") %>% left_join(lookup_df, by = c("cell_type", "gene")) %>%
                   mutate(
                     waldStat_1 = ifelse(is.na(waldStat_1), waldStat_1_lookup, waldStat_1),
                     meanLogFC  = ifelse(is.na(meanLogFC),  meanLogFC_lookup,  meanLogFC),
                     pvalue_1   = ifelse(is.na(pvalue_1),   pvalue_1_lookup,   pvalue_1),
                     qval_1     = ifelse(is.na(qval_1),     qval_1_lookup,     qval_1)) %>%
                   select(-waldStat_1_lookup, -meanLogFC_lookup, -pvalue_1_lookup, -qval_1_lookup) %>% column_to_rownames("row_id")

# Quick sanity check
cat("Remaining NA counts:\n")
print(colSums(is.na(pc_df[, c("waldStat_1", "meanLogFC", "pvalue_1", "qval_1")])))
print(head(pc_df[, c("gene", "waldStat_1", "meanLogFC", "pvalue_1", "qval_1")], 10))

# ---- notebook cell 31 ----

X <- z_row(traj_mat)

dt <- diff(tt)
Xd <- traj_mat[, -1, drop = FALSE] - traj_mat[, -ncol(traj_mat), drop = FALSE]
Xd <- sweep(Xd, 2, dt, "/")
Xd <- z_row(Xd)
stopifnot(nrow(X) == nrow(Xd)) # sanity

# ---- notebook cell 32 ----

dim(X)
dim(Xd)

significant_rows <- rownames(pc_df %>% filter(pvalue_1 < 0.05))
X <- X[significant_rows, , drop = FALSE]
Xd <- Xd[significant_rows, , drop = FALSE]

# Historical clustering cells use X_mat/Xd_mat names. Keep those aliases
after_filter_note <- "PCA/clustering use normalized trajectories and derivatives, not PC-score distances."
X_mat <- X
Xd_mat <- Xd
stopifnot(identical(rownames(X_mat), rownames(Xd_mat)))

dim(X_mat)
dim(Xd_mat)


# ---- notebook cell 36 ----

pca_X <- prcomp(X_mat, center = TRUE, scale. = FALSE)
cumvar <- cumsum(pca_X$sdev^2) / sum(pca_X$sdev^2)
n_pc <- which(cumvar >= 0.99)[1] # retain enough PCs for 95% variance

X_pc <- pca_X$x[ , seq_len(n_pc), drop = FALSE]
n_pc

# ---- notebook cell 37 ----

pca_Xd <- prcomp(Xd_mat, center = TRUE, scale. = FALSE)
cumvar <- cumsum(pca_Xd$sdev^2) / sum(pca_Xd$sdev^2)
n_pc <- which(cumvar >= 0.99)[1] # retain enough PCs for 95% variance

Xd_pc <- pca_Xd$x[ , seq_len(n_pc), drop = FALSE]
n_pc

# ---- notebook cell 38 ----

options(repr.plot.width = 4.8, repr.plot.height = 3.6, repr.plot.res = 300)

load_df_long <- data.frame(time = tt,
                           PC1 = pca_X$rotation[, 1],
                           PC2 = pca_X$rotation[, 2],
                           PC3 = pca_X$rotation[, 3],
                           PC4 = pca_X$rotation[, 4]) %>%
                pivot_longer(-time, names_to = "PC", values_to = "loading")


ggplot(load_df_long, aes(x = time, y = loading, color = PC)) + geom_line(linewidth = 1.2) + theme_classic()

# ---- notebook cell 39 ----

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

# ---- notebook cell 40 ----

pc_df = pc_df[rownames(pca_X$x), ]
pc_df$PC1 = pca_X$x[, 1]
pc_df$PC2 = pca_X$x[, 2]
pc_df$PC3 = pca_X$x[, 3]
pc_df$PC4 = pca_X$x[, 4]

# ---- notebook cell 42 ----

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


# ---- notebook cell 43 ----

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
    sil <- cluster::silhouette(cl, D)

    data.frame(
      K = K,
      WSS = sum(d_to_centroid^2),
      p90_centroid_distance = unname(quantile(d_to_centroid, 0.90, na.rm = TRUE)),
      mean_silhouette = mean(sil[, "sil_width"], na.rm = TRUE)
    )
  }

  metrics <- dplyr::bind_rows(lapply(K_grid, calc_for_k))

  WSS_minK <- metrics$WSS[1]
  WSS_maxK <- metrics$WSS[nrow(metrics)]

  metrics <- metrics %>%
    dplyr::mutate(
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


# ---- notebook cell 44 ----

scale_by_total_column_variance <- function(M) {
  M <- as.matrix(M)
  scale_factor <- sqrt(sum(apply(M, 2, var)))
  if (!is.finite(scale_factor) || scale_factor == 0) {
    stop("Cannot scale matrix: total column variance is zero or non-finite.")
  }
  M / scale_factor
}

# Use direct trajectory and derivative matrices for clustering.
# Each block is scaled by total column variance before alpha weighting.
# PCA scores are retained for PC1-PC4 visualization only.
X_n <- scale_by_total_column_variance(X_mat)
Xd_n <- scale_by_total_column_variance(Xd_mat)

cl.ls <- list()
sil_df.ls <- list()
D.ls <- list()


# ---- notebook cell 46 ----

lt = 0
alpha <- 0.25

subset.rownames = rownames(pc_df %>% filter(peak_log2fc > lt))
# subset.rownames = rownames(pc_df %>% filter(peak_log2fc < -lt))
# subset.rownames = rownames(pc_df)

# ---- notebook cell 47 ----

X_joint <- cbind(sqrt(1 - alpha)*X_n, sqrt(alpha)*Xd_n)[subset.rownames, ]
D_joint <- dist(X_joint, method = "euclidean")

D_x <- dist(X_mat[subset.rownames, ], method = "euclidean") 
D_pca_x = dist(X_pc[subset.rownames, ], method = "euclidean")
D_pca_dx = dist(Xd_pc[subset.rownames, ], method = "euclidean")

# ---- notebook cell 48 ----

D_use = D_joint
X_use = X_joint

hc <- hclust(D_use, method = "ward.D2")

clust_eval <- calculate_trajectory_clustering_metrics(X_mat = X_use, D = D_use, hc = hc, K_grid = 2:100, wss_fraction = 0.95, select_k_by = "elbow")
clust_eval$selection

options(repr.plot.width = 7.5, repr.plot.height = 5)
plot_trajectory_clustering_metrics(clust_eval, overlay = FALSE)

# ---- notebook cell 49 ----

cl <- cutree(hc, k = clust_eval$selection$k_silhouette)
gene_cluster_distance <- rep(NA_real_, nrow(X_use))
names(gene_cluster_distance) <- rownames(X_use)

for (k in sort(unique(cl))) {
  idx <- which(cl == k)
  Xi <- X_use[idx, , drop = FALSE]
  centroid <- colMeans(Xi)
  gene_cluster_distance[idx] <- sqrt(rowSums(sweep(Xi, 2, centroid, "-")^2))
}

disp_df <- data.frame(gene = rownames(X_use), cluster = cl, centroid_dist = gene_cluster_distance, stringsAsFactors = FALSE)

name <- "broad_induced"
cl.ls[[name]] <- cl
sil_df.ls[[name]] <- disp_df

cl <- cutree(hc, k = clust_eval$selection$k_elbow)
gene_cluster_distance <- rep(NA_real_, nrow(X_use))
names(gene_cluster_distance) <- rownames(X_use)

for (k in sort(unique(cl))) {
  idx <- which(cl == k)
  Xi <- X_use[idx, , drop = FALSE]
  centroid <- colMeans(Xi)
  gene_cluster_distance[idx] <- sqrt(rowSums(sweep(Xi, 2, centroid, "-")^2))
}

disp_df <- data.frame(gene = rownames(X_use), cluster = cl, centroid_dist = gene_cluster_distance, stringsAsFactors = FALSE)

name <- "fine_induced"
cl.ls[[name]] <- cl
sil_df.ls[[name]] <- disp_df

for (name in c("broad_induced", "fine_induced")) {
    # Add assignments back into pc_df
    pc_df[[paste0("cluster_", name)]] <- as.factor(cl.ls[[name]][rownames(pc_df)])
    pc_df[[paste0("cluster_dist_", name)]] <- sil_df.ls[[name]]$centroid_dist[match(rownames(pc_df), sil_df.ls[[name]]$gene)]
    
    # inspect cluster sizes
    print(name)
    print(sort(table(cl.ls[[name]]), decreasing = TRUE))
    
    cluster_name = paste0("cluster_", name)
    pc_df[[cluster_name]] <- factor(pc_df[[cluster_name]], levels = c(levels(pc_df[[cluster_name]]), "-1"))
    pc_df[[cluster_name]][is.na(pc_df[[cluster_name]])] <- "-1"
    pc_df[[cluster_name]] <- factor(pc_df[[cluster_name]])
}

# ---- notebook cell 50 ----

subset.rownames = rownames(pc_df %>% filter(peak_log2fc < -lt))

X_joint <- cbind(sqrt(1 - alpha)*X_n, sqrt(alpha)*Xd_n)[subset.rownames, ]
D_joint <- dist(X_joint, method = "euclidean")

D_x <- dist(X_mat[subset.rownames, ], method = "euclidean") 
D_pca_x = dist(X_pc[subset.rownames, ], method = "euclidean")
D_pca_dx = dist(Xd_pc[subset.rownames, ], method = "euclidean")

D_use = D_joint
X_use = X_joint

hc <- hclust(D_use, method = "ward.D2")

clust_eval <- calculate_trajectory_clustering_metrics(X_mat = X_use, D = D_use, hc = hc, K_grid = 2:100, wss_fraction = 0.95, select_k_by = "elbow")
clust_eval$selection

options(repr.plot.width = 7.5, repr.plot.height = 5)
plot_trajectory_clustering_metrics(clust_eval, overlay = FALSE)

# ---- notebook cell 51 ----

cl <- cutree(hc, k = clust_eval$selection$k_silhouette)
gene_cluster_distance <- rep(NA_real_, nrow(X_use))
names(gene_cluster_distance) <- rownames(X_use)

for (k in sort(unique(cl))) {
  idx <- which(cl == k)
  Xi <- X_use[idx, , drop = FALSE]
  centroid <- colMeans(Xi)
  gene_cluster_distance[idx] <- sqrt(rowSums(sweep(Xi, 2, centroid, "-")^2))
}

disp_df <- data.frame(gene = rownames(X_use), cluster = cl, centroid_dist = gene_cluster_distance, stringsAsFactors = FALSE)

name <- "broad_repressed"
cl.ls[[name]] <- cl
sil_df.ls[[name]] <- disp_df

cl <- cutree(hc, k = clust_eval$selection$k_elbow)
gene_cluster_distance <- rep(NA_real_, nrow(X_use))
names(gene_cluster_distance) <- rownames(X_use)

for (k in sort(unique(cl))) {
  idx <- which(cl == k)
  Xi <- X_use[idx, , drop = FALSE]
  centroid <- colMeans(Xi)
  gene_cluster_distance[idx] <- sqrt(rowSums(sweep(Xi, 2, centroid, "-")^2))
}

disp_df <- data.frame(gene = rownames(X_use), cluster = cl, centroid_dist = gene_cluster_distance, stringsAsFactors = FALSE)

name <- "fine_repressed"
cl.ls[[name]] <- cl
sil_df.ls[[name]] <- disp_df

for (name in c("broad_repressed", "fine_repressed")) {
    # Add assignments back into pc_df
    pc_df[[paste0("cluster_", name)]] <- as.factor(cl.ls[[name]][rownames(pc_df)])
    pc_df[[paste0("cluster_dist_", name)]] <- sil_df.ls[[name]]$centroid_dist[match(rownames(pc_df), sil_df.ls[[name]]$gene)]
    
    # inspect cluster sizes
    print(name)
    print(sort(table(cl.ls[[name]]), decreasing = TRUE))
    
    cluster_name = paste0("cluster_", name)
    pc_df[[cluster_name]] <- factor(pc_df[[cluster_name]], levels = c(levels(pc_df[[cluster_name]]), "-1"))
    pc_df[[cluster_name]][is.na(pc_df[[cluster_name]])] <- "-1"
    pc_df[[cluster_name]] <- factor(pc_df[[cluster_name]])
}

# ---- notebook cell 52 ----

RenameClusters <- function(x) {
    x <- as.character(x)
    if (x == '1') { return('impulse-up_1') }
    if (x == '2') { return('trend-up') }
    if (x == '3') { return('impulse-up_2') }
    else { return(x) }
}

pc_df$cluster_broad_induced = sapply(pc_df$cluster_broad_induced, RenameClusters, USE.NAMES = F)

# ---- notebook cell 53 ----

RenameClusters <- function(x) {
  x <- as.character(x)
  if (x == "1") return("trend-down_1")
  if (x == "2") return("trend-down_2")
  if (x == "3") return("transient-down")
  x
}

pc_df$cluster_broad_repressed <- sapply(pc_df$cluster_broad_repressed, RenameClusters, USE.NAMES = FALSE)


# ---- notebook cell 54 ----

pc_df$broad_cluster <- factor(ifelse(pc_df$cluster_broad_induced != "-1", as.character(pc_df$cluster_broad_induced), as.character(pc_df$cluster_broad_repressed)))

table(pc_df$broad_cluster)

# ---- notebook cell 55 ----

pc_df$cluster_dist_broad <- dplyr::coalesce(pc_df$cluster_dist_broad_induced, pc_df$cluster_dist_broad_repressed)
pc_df$cluster_dist_fine <- dplyr::coalesce(pc_df$cluster_dist_fine_induced, pc_df$cluster_dist_fine_repressed)

# Match the original reference row block order while preserving order within each block.
cell_type_order <- c("K562a", "K562b", "K562_puro", "THP1")
pc_df <- pc_df[order(match(as.character(pc_df$cell_type), cell_type_order)), , drop = FALSE]


# ---- notebook cell 57 ----

write.csv(pc_df, file.path(SELF_CONTAINED_OUTDIR, "pc_df.new.csv"))
saveRDS(pc_df, file.path(SELF_CONTAINED_OUTDIR, "pc_df.new.rds"))

cat("Wrote:\n")
cat("  ", file.path(SELF_CONTAINED_OUTDIR, "pc_df.new.csv"), "\n", sep = "")
cat("  ", file.path(SELF_CONTAINED_OUTDIR, "pc_df.new.rds"), "\n", sep = "")
