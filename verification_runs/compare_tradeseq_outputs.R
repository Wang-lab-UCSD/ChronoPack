# Compare a run's TradeSeqFits outputs with the manuscript template.
# Usage: Rscript compare_tradeseq_outputs.R RUN_ROOT [TEMPLATE_ROOT]
suppressPackageStartupMessages({
  library(SingleCellExperiment)
})

args <- commandArgs(trailingOnly = TRUE)
run_root <- args[[1]]
template_root <- if (length(args) >= 2) args[[2]] else
  "/new-stg/home/alex/tmp/ChronoPack/verification_runs/manuscript_full_pipeline_20260919_165539"

branches <- c(
  K562_TNF = "k562-tnfa.hg38.rep/intmd/sc-K562-1e4_to_tradeSeq.fixed",
  K562_puro = "k562-tnfa-puromycin.hg38/intmd/sc-K562-1e6-to-tradeSeq.fixed",
  LPS_K562 = "thp1-k562-lps.hg38.rep/intmd/sc-K562_to_tradeSeq.fixed",
  LPS_THP1 = "thp1-k562-lps.hg38.rep/intmd/sc-THP1_to_tradeSeq.fixed"
)
stat_cols <- c("waldStat_1", "pvalue_1", "meanLogFC", "qval_1")

max_rel <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  if (!any(ok)) return(NA_real_)
  max(abs(a[ok] - b[ok]) / pmax(abs(b[ok]), 1e-300))
}

get_beta <- function(sce) as.matrix(rowData(sce)$tradeSeq$beta[[1]])

for (name in names(branches)) {
  run_dir <- file.path(run_root, branches[[name]])
  ref_dir <- file.path(template_root, branches[[name]])
  cat(sprintf("\n== %s\n", name))
  run_sce_path <- file.path(run_dir, "sce.1trajectory.nknots8.rds")
  if (!file.exists(run_sce_path)) {
    cat("   not finished (no final SCE yet)\n")
    next
  }

  a <- read.csv(file.path(run_dir, "asso_df.1trajectory.nknots4.csv"), row.names = 1, check.names = FALSE)
  r <- read.csv(file.path(ref_dir, "asso_df.1trajectory.nknots4.csv"), row.names = 1, check.names = FALSE)
  common <- intersect(rownames(r), rownames(a))
  cat(sprintf("   association (nknots=4): run %d genes, ref %d, common %d, same order %s\n",
              nrow(a), nrow(r), length(common), identical(rownames(a), rownames(r))))
  for (col in stat_cols) {
    cat(sprintf("      %-10s max rel diff %.3g\n", col, max_rel(a[common, col], r[common, col])))
  }

  s <- rownames(read.csv(file.path(run_dir, "sig_df.1trajectory.nknots4.csv"), row.names = 1, check.names = FALSE))
  sr <- rownames(read.csv(file.path(ref_dir, "sig_df.1trajectory.nknots4.csv"), row.names = 1, check.names = FALSE))
  cat(sprintf("   significant (q<0.05): run %d, ref %d, identical set %s, identical order %s\n",
              length(s), length(sr), setequal(s, sr), identical(s, sr)))

  sce <- readRDS(run_sce_path)
  sce_ref <- readRDS(file.path(ref_dir, "sce.1trajectory.nknots8.rds"))
  g <- intersect(rownames(sce_ref), rownames(sce))
  cat(sprintf("   final SCE (nknots=8): run %d genes, ref %d, same genes %s, same cells %s\n",
              nrow(sce), nrow(sce_ref), setequal(rownames(sce), rownames(sce_ref)),
              identical(colnames(sce), colnames(sce_ref))))
  b <- get_beta(sce)[g, , drop = FALSE]
  br <- get_beta(sce_ref)[g, , drop = FALSE]
  cat(sprintf("      beta: same NA pattern %s, max abs diff %.3g\n",
              identical(is.na(b), is.na(br)), max(abs(b - br), na.rm = TRUE)))
  for (col in stat_cols) {
    if (col %in% colnames(rowData(sce)) && col %in% colnames(rowData(sce_ref))) {
      cat(sprintf("      rowData %-10s max rel diff %.3g\n", col,
                  max_rel(rowData(sce)[g, col], rowData(sce_ref)[g, col])))
    } else {
      cat(sprintf("      rowData %-10s missing (run %s, ref %s)\n", col,
                  col %in% colnames(rowData(sce)), col %in% colnames(rowData(sce_ref))))
    }
  }
  screen <- file.path(run_dir, "sce_screen.1trajectory.nknots4.rds")
  cat(sprintf("   screen SCE present: %s | duplicate copies present: %s\n", file.exists(screen),
              paste(basename(Filter(file.exists, file.path(run_dir, c("sce.1trajectory.nknots8.csv", "sce.knots8.rds")))), collapse = ", ")))
}
