find_parent_with <- function(marker) {
  start <- normalizePath(getwd(), mustWork = TRUE)
  parts <- strsplit(start, .Platform$file.sep, fixed = TRUE)[[1]]
  if (parts[1] == "") parts[1] <- .Platform$file.sep
  candidates <- vapply(seq_along(parts), function(i) {
    if (i == 1 && parts[1] == .Platform$file.sep) return(.Platform$file.sep)
    do.call(file.path, as.list(parts[seq_len(i)]))
  }, character(1))
  candidates <- rev(candidates)
  hit <- candidates[file.exists(file.path(candidates, marker))]
  if (!length(hit)) {
    stop("Could not find parent containing ", marker, " from ", start)
  }
  normalizePath(hit[[1]], mustWork = TRUE)
}

repo_root <- find_parent_with("chronopack/R/tradeseq_helpers.R")
experiment_root <- find_parent_with("intmd")
source(file.path(repo_root, "chronopack", "R", "tradeseq_helpers.R"))
load_chronopack_tradeseq_libraries()
suppressPackageStartupMessages(library(dplyr))

set_tradeseq_workdir(file.path(experiment_root, "intmd", "sc-THP1_to_tradeSeq.fixed"))

TRADESEQ_CONFIG <- list(
  label = "lps_thp1",
  n_points = 6,
  workers = 4,
  screen_nknots = 4,
  final_nknots = 8,
  p_threshold = 0.05,
  outputs = list(
    association = "asso_df.1trajectory.nknots4.csv",
    significant = "sig_df.1trajectory.nknots4.csv",
    sce = "sce.1trajectory.nknots8.rds",
    gamList = "gamList.1trajectory.nknots8.rds"
  ),
  legacy_outputs = list(
    association = c("asso_df.thp1.1trajectory.nknots4.expanded.csv"),
    significant = c("sig_df.thp1.1trajectory.nknots4.expanded.csv"),
    sce = c("sce.knots8.rds"),
    gamList = c("gamList.knots8.rds")
  )
)

TRADESEQ_CONFIG


inputs <- read_tradeseq_inputs(tag = "")

counts <- inputs$counts
cell_times <- inputs$cell_times
libs_full <- inputs$libs_full
time_vec <- inputs$time_vec
lib_full <- inputs$lib_full
subset_genes <- inputs$subset_genes
genes_keep_initial <- inputs$genes_keep
counts_sub <- inputs$counts_sub
sf <- inputs$sf

stopifnot(
  identical(colnames(counts_sub), names(time_vec)),
  length(sf) == ncol(counts_sub)
)

list(
  counts_dim = dim(counts),
  counts_sub_dim = dim(counts_sub),
  n_subset_genes = length(subset_genes),
  n_genes_kept_initial = length(genes_keep_initial),
  time_range = range(time_vec)
)


single_lineage <- prepare_single_lineage(counts_sub, time_vec, sf)
PT <- single_lineage$PT
CW <- single_lineage$CW
offset_mat <- single_lineage$offset_mat

stopifnot(
  nrow(PT) == ncol(counts_sub),
  nrow(CW) == ncol(counts_sub),
  all(dim(offset_mat) == dim(counts_sub))
)


sce_screen <- fit_tradeseq_sce(
  counts_sub = counts_sub,
  pseudotime = PT,
  cell_weights = CW,
  offset_mat = offset_mat,
  nknots = TRADESEQ_CONFIG$screen_nknots,
  workers = TRADESEQ_CONFIG$workers
)

asso_df <- association_df_single(sce_screen, n_points = TRADESEQ_CONFIG$n_points)
asso_df <- asso_df[!is.na(asso_df$waldStat_1), , drop = FALSE]

sig_df <- asso_df[!is.na(asso_df$qval_1) & asso_df$qval_1 < 0.05, , drop = FALSE]
sig_df <- sig_df[order(-sig_df$waldStat_1, -sig_df$meanLogFC), , drop = FALSE]

genes_keep <- significant_single_lineage(
  asso_df,
  p_threshold = TRADESEQ_CONFIG$p_threshold
)

list(
  n_tested = nrow(asso_df),
  n_qval_significant = nrow(sig_df),
  n_pvalue_significant = length(genes_keep)
)


write_csv_with_aliases(
  asso_df,
  TRADESEQ_CONFIG$outputs$association,
  TRADESEQ_CONFIG$legacy_outputs$association
)
write_csv_with_aliases(
  sig_df,
  TRADESEQ_CONFIG$outputs$significant,
  TRADESEQ_CONFIG$legacy_outputs$significant
)


counts_sig <- counts_sub[genes_keep, , drop = FALSE]
offset_sig <- make_offset_matrix(counts_sig, sf)

stopifnot(
  nrow(counts_sig) == length(genes_keep),
  all(dim(offset_sig) == dim(counts_sig))
)

sce <- fit_tradeseq_sce(
  counts_sub = counts_sig,
  pseudotime = PT,
  cell_weights = CW,
  offset_mat = offset_sig,
  nknots = TRADESEQ_CONFIG$final_nknots,
  workers = TRADESEQ_CONFIG$workers
)
save_rds_with_aliases(
  sce,
  TRADESEQ_CONFIG$outputs$sce,
  TRADESEQ_CONFIG$legacy_outputs$sce
)

gamList <- fit_tradeseq_gamlist(
  counts_sub = counts_sig,
  pseudotime = PT,
  cell_weights = CW,
  offset_mat = offset_sig,
  nknots = TRADESEQ_CONFIG$final_nknots,
  workers = TRADESEQ_CONFIG$workers
)
save_rds_with_aliases(
  gamList,
  TRADESEQ_CONFIG$outputs$gamList,
  TRADESEQ_CONFIG$legacy_outputs$gamList
)

list(
  final_gene_count = nrow(counts_sig),
  sce_path = TRADESEQ_CONFIG$outputs$sce,
  gamList_path = TRADESEQ_CONFIG$outputs$gamList
)
