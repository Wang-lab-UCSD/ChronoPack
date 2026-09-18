load_chronopack_tradeseq_libraries <- function() {
  suppressPackageStartupMessages({
    library(tradeSeq)
    library(RColorBrewer)
    library(SingleCellExperiment)
    library(slingshot)
    library(Matrix)
    library(BiocParallel)
    library(ggplot2)
  })
  RNGversion("3.5.0")
  palette(RColorBrewer::brewer.pal(8, "Dark2"))
}

report_system_ram <- function() {
  meminfo <- readLines("/proc/meminfo")
  get_kb <- function(key) {
    line <- meminfo[grep(paste0("^", key, ":"), meminfo)][1]
    as.numeric(gsub(".*:\\s+([0-9]+)\\s+kB.*", "\\1", line))
  }
  total_gb <- get_kb("MemTotal") / 1024^2
  avail_gb <- get_kb("MemAvailable") / 1024^2
  free_gb <- get_kb("MemFree") / 1024^2
  cat(sprintf(
    "System RAM: total = %.2f GB | available = %.2f GB | free = %.2f GB\n",
    total_gb, avail_gb, free_gb
  ))
}

set_tradeseq_workdir <- function(input_dir) {
  setwd(input_dir)
  invisible(normalizePath(getwd()))
}

tagged_path <- function(stem, tag = "", ext = "") {
  candidates <- if (identical(tag, "")) {
    c(
      paste0(stem, ext),
      paste0(stem, ".", ext)
    )
  } else {
    c(
      paste0(stem, ".", tag, ext),
      paste0(stem, ext, ".", tag),
      paste0(stem, ".", tag, ".", ext)
    )
  }
  candidates <- unique(candidates)
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    stop("Could not find file for ", stem, " tag='", tag, "'. Tried: ",
         paste(candidates, collapse = ", "))
  }
  hit[[1]]
}

read_tradeseq_inputs <- function(tag = "") {
  counts <- as.matrix(read.csv(
    gzfile(tagged_path("counts_genes_x_cells", tag, ".csv.gz")),
    row.names = 1,
    check.names = FALSE
  ))
  storage.mode(counts) <- "integer"

  cell_times <- read.csv(
    gzfile(tagged_path("cell_times", tag, ".csv.gz")),
    stringsAsFactors = FALSE
  )
  libs_full <- read.csv("cell_libsize_full.csv", stringsAsFactors = FALSE)
  subset_genes <- read.csv(
    tagged_path("subset_genes", tag, ".csv"),
    stringsAsFactors = FALSE
  )[, 1]

  stopifnot(
    all(colnames(counts) %in% cell_times$cell_id),
    all(colnames(counts) %in% libs_full$cell_id)
  )

  time_vec <- setNames(cell_times$time, cell_times$cell_id)[colnames(counts)]
  lib_full <- setNames(libs_full$libsize_full, libs_full$cell_id)[colnames(counts)]
  genes_keep <- intersect(subset_genes, rownames(counts))
  counts_sub <- counts[genes_keep, , drop = FALSE]
  sf <- as.numeric(lib_full)
  sf <- sf / mean(sf)

  list(
    counts = counts,
    cell_times = cell_times,
    libs_full = libs_full,
    time_vec = time_vec,
    lib_full = lib_full,
    subset_genes = subset_genes,
    genes_keep = genes_keep,
    counts_sub = counts_sub,
    sf = sf
  )
}

prepare_single_lineage <- function(counts_sub, time_vec, sf) {
  pt01 <- (time_vec - min(time_vec, na.rm = TRUE)) /
    diff(range(time_vec, na.rm = TRUE))
  PT <- matrix(pt01, ncol = 1)
  colnames(PT) <- "lineage1"
  CW <- matrix(1, nrow = ncol(counts_sub), ncol = 1)
  colnames(CW) <- "lineage1"
  offset_mat <- make_offset_matrix(counts_sub, sf)
  list(PT = PT, CW = CW, offset_mat = offset_mat)
}

prepare_hard_lineages <- function(counts_sub, time_vec, sf, cell_type_path = "cell_types.csv") {
  pt01 <- (time_vec - min(time_vec, na.rm = TRUE)) /
    diff(range(time_vec, na.rm = TRUE))
  PT <- matrix(pt01, ncol = 1)
  colnames(PT) <- "lineage1"

  cell_types <- read.csv(cell_type_path, row.names = 1, check.names = FALSE)
  types <- unique(cell_types$cell_type)
  stopifnot(length(types) == 2)

  ct <- cell_types[colnames(counts_sub), , drop = FALSE]
  stopifnot(all(rownames(ct) == colnames(counts_sub)))

  PT_mat <- cbind(K562 = PT, THP1 = PT)
  rownames(PT_mat) <- colnames(counts_sub)

  CW <- matrix(0, nrow = nrow(ct), ncol = 2, dimnames = list(rownames(ct), types))
  CW[cbind(seq_len(nrow(ct)), match(ct$cell_type, types))] <- 1
  CW <- sweep(CW, 1, rowSums(CW), "/")
  stopifnot(all(abs(rowSums(CW) - 1) < 1e-8))

  CWz <- CW
  CWz[CWz < 1e-6] <- 0
  list(
    PT = PT,
    PT_mat = PT_mat,
    CW = CW,
    CWz = CWz,
    cell_types = cell_types,
    aligned_cell_types = ct,
    offset_mat = make_offset_matrix(counts_sub, sf)
  )
}

make_offset_matrix <- function(counts_sub, sf) {
  matrix(log(sf), nrow = nrow(counts_sub), ncol = ncol(counts_sub), byrow = TRUE)
}

fit_tradeseq_sce <- function(counts_sub, pseudotime, cell_weights, offset_mat,
                             nknots = 4, workers = 4, verbose = TRUE) {
  tradeSeq::fitGAM(
    counts = counts_sub,
    pseudotime = pseudotime,
    cellWeights = cell_weights,
    offset = offset_mat,
    nknots = nknots,
    BPPARAM = BiocParallel::MulticoreParam(workers = workers),
    verbose = verbose
  )
}

fit_tradeseq_gamlist <- function(counts_sub, pseudotime, cell_weights, offset_mat,
                                 nknots = 8, workers = 4, verbose = TRUE) {
  tradeSeq::fitGAM(
    counts = counts_sub,
    pseudotime = pseudotime,
    cellWeights = cell_weights,
    offset = offset_mat,
    nknots = nknots,
    sce = FALSE,
    BPPARAM = BiocParallel::MulticoreParam(workers = workers),
    verbose = verbose
  )
}

association_df_single <- function(sce, n_points = 9) {
  asso_res <- tradeSeq::associationTest(
    sce,
    global = FALSE,
    lineages = TRUE,
    inverse = "Chol",
    nPoints = n_points
  )
  asso_res <- asso_res[order(asso_res[, "waldStat_1"], decreasing = TRUE), ]
  asso_df <- as.data.frame(asso_res)
  if (!"meanLogFC" %in% names(asso_df) && ncol(asso_df) == 4) {
    names(asso_df) <- c("waldStat_1", "df_1", "pvalue_1", "meanLogFC")
  }
  asso_df$qval_1 <- p.adjust(asso_df$pvalue_1, method = "BH")
  asso_df
}

association_df_multilineage <- function(sce, n_points = 9) {
  asso_res <- tradeSeq::associationTest(
    sce,
    global = TRUE,
    lineages = TRUE,
    inverse = "Chol",
    nPoints = n_points
  )
  asso_res <- asso_res[order(asso_res[, "waldStat_1"], decreasing = TRUE), ]
  asso_df <- as.data.frame(asso_res)
  if ("pvalue" %in% names(asso_df)) {
    asso_df$qval <- p.adjust(asso_df$pvalue, method = "BH")
  }
  if ("pvalue_1" %in% names(asso_df)) {
    asso_df$qval_1 <- p.adjust(asso_df$pvalue_1, method = "BH")
  }
  if ("pvalue_2" %in% names(asso_df)) {
    asso_df$qval_2 <- p.adjust(asso_df$pvalue_2, method = "BH")
  }
  asso_df
}

significant_single_lineage <- function(asso_df, p_threshold = 0.05,
                                      q_threshold = NULL) {
  if (!is.null(q_threshold)) {
    sig_df <- subset(asso_df, !is.na(qval_1) & qval_1 < q_threshold)
  } else {
    sig_df <- subset(asso_df, !is.na(pvalue_1) & pvalue_1 < p_threshold)
  }
  sig_df <- sig_df[!is.na(sig_df$waldStat_1), , drop = FALSE]
  rownames(sig_df)
}

significant_multilineage <- function(asso_df, q_threshold = 0.01) {
  fields <- intersect(c("qval", "qval_1", "qval_2"), names(asso_df))
  if (!length(fields)) {
    stop("No q-value columns available for multilineage significance.")
  }
  mask <- Reduce(`|`, lapply(fields, function(x) asso_df[[x]] < q_threshold))
  mask[is.na(mask)] <- FALSE
  rownames(asso_df)[mask]
}

write_csv_with_aliases <- function(x, path, aliases = character(), row.names = TRUE) {
  paths <- unique(c(path, aliases))
  for (p in paths) {
    write.csv(x, p, row.names = row.names)
  }
  invisible(normalizePath(paths, mustWork = FALSE))
}

save_rds_with_aliases <- function(x, path, aliases = character()) {
  paths <- unique(c(path, aliases))
  for (p in paths) {
    saveRDS(x, p)
  }
  invisible(normalizePath(paths, mustWork = FALSE))
}
