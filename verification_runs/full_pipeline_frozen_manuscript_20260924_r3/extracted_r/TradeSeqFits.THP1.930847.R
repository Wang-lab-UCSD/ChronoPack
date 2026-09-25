# Extracted from /new-stg/home/alex/tmp/ChronoPack/verification_runs/full_pipeline_frozen_manuscript_20260924_r3/thp1-k562-lps.hg38.rep/notebooks/TradeSeqFits.THP1.ipynb
options(error = function() { traceback(2); quit(status = 1) })

message('--- running TradeSeqFits.THP1.ipynb cell 2 ---')
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

message('--- running TradeSeqFits.THP1.ipynb cell 3 ---')
TRADESEQ_CONFIG <- list(
  label = "lps_thp1",
  input_dir = file.path(experiment_root, "intmd", "sc-THP1_to_tradeSeq.fixed"),
  n_points = 6,
  workers = 4,
  screen_nknots = 4,
  final_nknots = 8,
  p_threshold = 0.05,
  outputs = list(
    association = "asso_df.1trajectory.nknots4.csv",
    significant = "sig_df.1trajectory.nknots4.csv",
    sce = "sce.1trajectory.nknots8.rds",
  ),
  legacy_outputs = list(
    association = c("asso_df.thp1.1trajectory.nknots4.expanded.csv"),
    significant = c("sig_df.thp1.1trajectory.nknots4.expanded.csv"),
    sce = c("sce.knots8.rds"),
  )
)

TRADESEQ_CONFIG

message('--- running TradeSeqFits.THP1.ipynb cell 4 ---')
fit_result <- run_tradeseq_single_lineage(
  input_dir = TRADESEQ_CONFIG$input_dir,
  n_points = TRADESEQ_CONFIG$n_points,
  workers = TRADESEQ_CONFIG$workers,
  screen_nknots = TRADESEQ_CONFIG$screen_nknots,
  final_nknots = TRADESEQ_CONFIG$final_nknots,
  p_threshold = TRADESEQ_CONFIG$p_threshold,
  output_association = TRADESEQ_CONFIG$outputs$association,
  output_significant = TRADESEQ_CONFIG$outputs$significant,
  output_sce = TRADESEQ_CONFIG$outputs$sce,
  association_aliases = TRADESEQ_CONFIG$legacy_outputs$association,
  significant_aliases = TRADESEQ_CONFIG$legacy_outputs$significant,
  sce_aliases = TRADESEQ_CONFIG$legacy_outputs$sce,
  tag = "",
  verbose = TRUE
)

list(
  label = TRADESEQ_CONFIG$label,
  input_dir = fit_result$input_dir,
  counts_dim = fit_result$counts_dim,
  counts_sub_dim = fit_result$counts_sub_dim,
  n_tested = fit_result$n_tested,
  n_qval_significant = fit_result$n_qval_significant,
  n_pvalue_significant = fit_result$n_pvalue_significant,
  final_counts_dim = fit_result$final_counts_dim,
  association_path = fit_result$association_path,
  significant_path = fit_result$significant_path,
  sce_path = fit_result$sce_path,
)
