suppressPackageStartupMessages({
  library(fdrtool)
})

prepare_p <- function(p_raw) {
  min_pos <- min(p_raw[p_raw > 0], na.rm = TRUE)
  p <- if (is.finite(min_pos)) pmax(p_raw, min_pos / 10) else pmax(p_raw, 1e-300)
  pmin(p, 1 - 1e-16)
}

load_valid_targets <- function(valid_targets_path) {
  valid_targets <- read.csv(valid_targets_path, header = FALSE, stringsAsFactors = FALSE)[, 1]
  valid_targets[!grepl("^(MT-|ENSG)", valid_targets)]
}

load_gene_table <- function(label, path, valid_targets) {
  df <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  # Gene names: an "X" column if present, otherwise the first column (which may have an empty header).
  genes <- if ("X" %in% names(df)) df[["X"]] else df[[1]]
  keep <- genes %in% valid_targets & is.finite(df$pvalue_1)
  df <- df[keep, , drop = FALSE]
  out <- data.frame(
    dataset = label,
    gene = genes[keep],
    p_raw = df$pvalue_1,
    p = prepare_p(df$pvalue_1),
    stringsAsFactors = FALSE
  )
  # Carried along: the tradeSeq Wald statistic (tie-breaker for equally scored
  # targets in GRN edges) and an optional previously stored local FDR.
  if ("waldStat_1" %in% names(df)) out$waldStat_1 <- df$waldStat_1
  if ("lfdr" %in% names(df)) out$stored_lfdr <- df$lfdr
  out
}

safe_fdrtool <- function(p, cutoff_method) {
  tryCatch(
    fdrtool(p, statistic = "pvalue", plot = FALSE, verbose = FALSE, cutoff.method = cutoff_method),
    error = function(e) NULL
  )
}

fit_bum <- function(p) {
  p <- pmin(pmax(p, .Machine$double.xmin), 1 - 1e-12)
  negll <- function(theta) {
    pi0 <- plogis(theta[1])
    a <- plogis(theta[2])
    dens <- pi0 + (1 - pi0) * a * p^(a - 1)
    -sum(log(pmax(dens, .Machine$double.xmin)))
  }
  best <- NULL
  starts <- rbind(
    c(qlogis(0.9), qlogis(0.2)),
    c(qlogis(0.8), qlogis(0.1)),
    c(qlogis(0.95), qlogis(0.4)),
    c(qlogis(0.99), qlogis(0.05))
  )
  for (i in seq_len(nrow(starts))) {
    fit <- optim(starts[i, ], negll, method = "BFGS", control = list(maxit = 1000))
    if (is.null(best) || fit$value < best$value) {
      best <- fit
    }
  }
  pi0 <- plogis(best$par[1])
  a <- plogis(best$par[2])
  dens <- pi0 + (1 - pi0) * a * p^(a - 1)
  lfdr <- pmin(1, pi0 / pmax(dens, .Machine$double.xmin))
  list(lfdr = lfdr, pi0 = pi0, a = a, nll = best$value)
}

score_table <- function(tab) {
  p <- tab$p
  out <- tab
  methods <- list()

  for (method in c("fndr", "pct0", "locfdr")) {
    ft <- safe_fdrtool(p, method)
    col <- paste0("lfdr_fdrtool_", method)
    out[[col]] <- if (is.null(ft)) NA_real_ else ft$lfdr
    methods[[paste0("fdrtool_", method)]] <- if (is.null(ft)) {
      data.frame(method = paste0("fdrtool_", method), param = "fit_failed", value = NA_real_)
    } else {
      data.frame(method = paste0("fdrtool_", method), param = colnames(ft$param), value = as.numeric(ft$param[1, ]))
    }
  }

  out$bh_q <- p.adjust(p, method = "BH")
  out$lfdr_like_bh_q <- out$bh_q

  bum <- fit_bum(p)
  out$lfdr_bum <- bum$lfdr
  methods[["bum"]] <- data.frame(
    method = "bum",
    param = c("pi0", "a", "nll"),
    value = c(bum$pi0, bum$a, bum$nll)
  )

  score_cap <- 10
  out$gene_score_logp_cap10 <- pmin(-log10(pmax(tab$p_raw, .Machine$double.xmin)), score_cap) / score_cap
  out$lfdr_like_logp_cap10 <- 1 - out$gene_score_logp_cap10

  list(scores = out, params = do.call(rbind, methods))
}

summarise_method <- function(scores, params) {
  lfdr_cols <- grep("^(stored_lfdr|lfdr_|bh_q$)", names(scores), value = TRUE)
  rows <- lapply(lfdr_cols, function(col) {
    x <- scores[[col]]
    data.frame(
      dataset = scores$dataset[1],
      method = col,
      n = sum(is.finite(x)),
      below_0p01 = sum(x < 0.01, na.rm = TRUE),
      below_0p05 = sum(x < 0.05, na.rm = TRUE),
      below_0p10 = sum(x < 0.10, na.rm = TRUE),
      below_0p20 = sum(x < 0.20, na.rm = TRUE),
      exact_one = sum(x >= 1, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      q10 = quantile(x, 0.10, na.rm = TRUE),
      q90 = quantile(x, 0.90, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  list(summary = summary, params = transform(params, dataset = scores$dataset[1]))
}

calculate_lfdr_recalibration <- function(valid_targets_path, gene_tables) {
  valid_targets <- load_valid_targets(valid_targets_path)
  all_scores <- list()
  all_summaries <- list()
  all_params <- list()

  for (label in names(gene_tables)) {
    tab <- load_gene_table(label, gene_tables[[label]], valid_targets)
    res <- score_table(tab)
    sm <- summarise_method(res$scores, res$params)
    all_scores[[label]] <- res$scores
    all_summaries[[label]] <- sm$summary
    all_params[[label]] <- sm$params
  }

  list(
    scores = do.call(rbind, all_scores),
    summary = do.call(rbind, all_summaries),
    params = do.call(rbind, all_params)
  )
}

write_lfdr_recalibration_outputs <- function(result, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- list(
    scores = file.path(out_dir, "thp1_k562_lfdr_recalibration_scores.csv"),
    summary = file.path(out_dir, "thp1_k562_lfdr_recalibration_summary.csv"),
    params = file.path(out_dir, "thp1_k562_lfdr_recalibration_params.csv")
  )
  write.csv(result$scores, paths$scores, row.names = FALSE)
  write.csv(result$summary, paths$summary, row.names = FALSE)
  write.csv(result$params, paths$params, row.names = FALSE)
  paths
}
