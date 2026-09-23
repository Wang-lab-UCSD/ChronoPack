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
  offset_vec <- make_offset_vector(sf)
  list(PT = PT, CW = CW, offset_vec = offset_vec)
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

make_offset_vector <- function(sf) {
  as.numeric(log(sf))
}

fit_tradeseq_sce <- function(counts_sub, pseudotime, cell_weights, offset,
                             nknots = 4, workers = 4, verbose = TRUE,
                             parallel = TRUE, genes = NULL) {
  if (length(offset) != ncol(counts_sub)) {
    stop("offset must be a cell-level vector with length ncol(counts_sub).")
  }
  if (nrow(pseudotime) != ncol(counts_sub)) {
    stop("pseudotime must have one row per cell.")
  }
  if (nrow(cell_weights) != ncol(counts_sub)) {
    stop("cell_weights must have one row per cell.")
  }
  tradeSeq::fitGAM(
    counts = counts_sub,
    genes = genes,
    pseudotime = pseudotime,
    cellWeights = cell_weights,
    offset = offset,
    nknots = nknots,
    sce = TRUE,
    parallel = parallel,
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

annotate_sce_association <- function(sce, asso_df,
                                     fields = c("waldStat_1", "meanLogFC", "pvalue_1", "qval_1")) {
  if (is.null(asso_df) || !nrow(asso_df)) {
    return(sce)
  }
  if (is.null(rownames(asso_df))) {
    stop("Association dataframe must have gene rownames.")
  }
  rd <- SummarizedExperiment::rowData(sce)
  gene_names <- rownames(sce)
  for (field in fields) {
    rd[[field]] <- NA_real_
    if (field %in% colnames(asso_df)) {
      idx <- match(gene_names, rownames(asso_df))
      ok <- !is.na(idx)
      rd[[field]][ok] <- as.numeric(asso_df[[field]][idx[ok]])
    }
  }
  SummarizedExperiment::rowData(sce) <- rd
  sce
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

run_tradeseq_single_lineage <- function(input_dir,
                                        n_points = 9,
                                        workers = 4,
                                        screen_nknots = 4,
                                        final_nknots = 8,
                                        p_threshold = 0.05,
                                        output_association = "asso_df.1trajectory.nknots4.csv",
                                        output_significant = "sig_df.1trajectory.nknots4.csv",
                                        output_sce = "sce.1trajectory.nknots8.rds",
                                        output_sce_screen = "sce_screen.1trajectory.nknots4.rds",
                                        output_screen_genes = "screen_genes.1trajectory.nknots4.csv",
                                        output_session_info = "sessionInfo.tradeSeqFits.txt",
                                        association_aliases = character(),
                                        significant_aliases = character(),
                                        sce_aliases = character(),
                                        sce_screen_aliases = character(),
                                        screen_genes_aliases = character(),
                                        tag = "",
                                        seed = 123,
                                        verbose = TRUE) {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  set_tradeseq_workdir(input_dir)

  inputs <- read_tradeseq_inputs(tag = tag)
  counts_sub <- inputs$counts_sub
  time_vec <- inputs$time_vec
  sf <- inputs$sf

  stopifnot(
    identical(colnames(counts_sub), names(time_vec)),
    length(sf) == ncol(counts_sub)
  )

  single_lineage <- prepare_single_lineage(counts_sub, time_vec, sf)
  PT <- single_lineage$PT
  CW <- single_lineage$CW
  offset_vec <- single_lineage$offset_vec

  set.seed(seed)
  sce_screen <- fit_tradeseq_sce(
    counts_sub = counts_sub,
    pseudotime = PT,
    cell_weights = CW,
    offset = offset_vec,
    nknots = screen_nknots,
    workers = workers,
    verbose = verbose
  )
  asso_df <- association_df_single(sce_screen, n_points = n_points)
  asso_df <- asso_df[!is.na(asso_df$waldStat_1), , drop = FALSE]
  sce_screen <- annotate_sce_association(sce_screen, asso_df)
  sce_screen_paths <- save_rds_with_aliases(
    sce_screen,
    output_sce_screen,
    sce_screen_aliases
  )

  fdr_sig_df <- asso_df[!is.na(asso_df$qval_1) & asso_df$qval_1 < 0.05, , drop = FALSE]
  fdr_sig_df <- fdr_sig_df[order(-fdr_sig_df$waldStat_1, -fdr_sig_df$meanLogFC), , drop = FALSE]

  screen_genes <- significant_single_lineage(
    asso_df,
    p_threshold = p_threshold
  )

  write_csv_with_aliases(
    asso_df,
    output_association,
    association_aliases
  )
  write_csv_with_aliases(
    fdr_sig_df,
    output_significant,
    significant_aliases
  )
  write_csv_with_aliases(
    data.frame(gene = screen_genes, stringsAsFactors = FALSE),
    output_screen_genes,
    screen_genes_aliases,
    row.names = FALSE
  )

  counts_sub2 <- counts_sub[screen_genes, , drop = FALSE]

  set.seed(seed)
  sce_final <- fit_tradeseq_sce(
    counts_sub = counts_sub2,
    pseudotime = PT,
    cell_weights = CW,
    offset = offset_vec,
    nknots = final_nknots,
    workers = workers,
    verbose = verbose
  )
  sce_final <- annotate_sce_association(sce_final, asso_df)
  sce_paths <- save_rds_with_aliases(
    sce_final,
    output_sce,
    sce_aliases
  )

  session_info_path <- normalizePath(output_session_info, mustWork = FALSE)
  writeLines(c(
    paste0("tradeSeq: ", as.character(utils::packageVersion("tradeSeq"))),
    paste0("mgcv: ", as.character(utils::packageVersion("mgcv"))),
    capture.output(utils::sessionInfo())
  ), con = output_session_info)

  list(
    input_dir = normalizePath(getwd(), mustWork = TRUE),
    counts_dim = dim(inputs$counts),
    counts_sub_dim = dim(counts_sub),
    final_counts_dim = dim(counts_sub2),
    n_subset_genes = length(inputs$subset_genes),
    n_genes_kept_initial = length(inputs$genes_keep),
    n_tested = nrow(asso_df),
    n_qval_significant = nrow(fdr_sig_df),
    n_pvalue_significant = length(screen_genes),
    genes_keep = screen_genes,
    screen_genes = screen_genes,
    association_path = normalizePath(output_association, mustWork = FALSE),
    significant_path = normalizePath(output_significant, mustWork = FALSE),
    sce_path = normalizePath(output_sce, mustWork = FALSE),
    sce_final_path = normalizePath(output_sce, mustWork = FALSE),
    sce_screen_path = normalizePath(output_sce_screen, mustWork = FALSE),
    screen_genes_path = normalizePath(output_screen_genes, mustWork = FALSE),
    session_info_path = session_info_path,
    sce_paths = sce_paths,
    sce_screen_paths = sce_screen_paths,
    pseudotime = PT,
    cell_weights = CW,
    size_factors = sf,
    association = asso_df,
    significant = fdr_sig_df
  )
}
