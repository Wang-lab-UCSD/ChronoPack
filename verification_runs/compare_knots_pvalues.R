# Compare nknots=4 screen p-values with p-values recomputed from the nknots=8 fit.
# The final SCE's rowData p-values are copied from the nknots=4 screen, so the
# nknots=8 association test is recomputed here with the pipeline's settings.
# Usage: Rscript compare_knots_pvalues.R RUN_ROOT OUT_CSV
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(tradeSeq)
})

args <- commandArgs(trailingOnly = TRUE)
run_root <- args[[1]]
out_csv <- args[[2]]

branches <- list(
  K562_TNF = list(dir = "k562-tnfa.hg38.rep/intmd/sc-K562-1e4_to_tradeSeq.fixed", n_points = 9),
  K562_puro = list(dir = "k562-tnfa-puromycin.hg38/intmd/sc-K562-1e6-to-tradeSeq.fixed", n_points = 9),
  LPS_K562 = list(dir = "thp1-k562-lps.hg38.rep/intmd/sc-K562_to_tradeSeq.fixed", n_points = 6),
  LPS_THP1 = list(dir = "thp1-k562-lps.hg38.rep/intmd/sc-THP1_to_tradeSeq.fixed", n_points = 6)
)

rows <- list()
for (name in names(branches)) {
  b <- branches[[name]]
  dir <- file.path(run_root, b$dir)
  sce8_path <- file.path(dir, "sce.1trajectory.nknots8.rds")
  if (!file.exists(sce8_path)) next

  a4 <- read.csv(file.path(dir, "asso_df.1trajectory.nknots4.csv"), row.names = 1, check.names = FALSE)
  sce8 <- readRDS(sce8_path)
  a8 <- as.data.frame(associationTest(sce8, global = FALSE, lineages = TRUE,
                                      inverse = "Chol", nPoints = b$n_points))
  a8 <- a8[!is.na(a8$waldStat_1), , drop = FALSE]
  g <- intersect(rownames(a8), rownames(a4))
  p4 <- a4[g, "pvalue_1"]
  p8 <- a8[g, "pvalue_1"]
  q4 <- a4[g, "qval_1"]
  q8_within <- p.adjust(p8, method = "BH")  # BH over the 8-knot gene set only

  sig4 <- p4 < 0.05
  sig8 <- p8 < 0.05
  ok <- is.finite(p4) & is.finite(p8) & p4 > 0 & p8 > 0
  rho <- suppressWarnings(cor(-log10(p4[ok]), -log10(p8[ok]), method = "spearman"))
  both_q <- sum(q4 < 0.05 & q8_within < 0.05, na.rm = TRUE)

  cat(sprintf("\n== %s (nPoints=%d)\n", name, b$n_points))
  cat(sprintf("   genes screened at nknots=4 (all candidates): %d | p4<0.05 (= fitted at nknots=8): %d\n",
              nrow(a4), sum(a4$pvalue_1 < 0.05, na.rm = TRUE)))
  cat(sprintf("   nknots=8 fit: %d genes, %d with finite Wald statistics\n", nrow(sce8), length(g)))
  cat(sprintf("   p8<0.05: %d of %d (%.1f%%); of those, p4<0.05: %d (%.1f%%)\n",
              sum(sig8), length(g), 100 * mean(sig8), sum(sig8 & sig4), 100 * sum(sig8 & sig4) / max(sum(sig8), 1)))
  cat(sprintf("   screened genes losing p<0.05 at nknots=8: %d\n", sum(sig4 & !sig8)))
  cat(sprintf("   Spearman rho(-log10 p4, -log10 p8): %.3f | median p4 %.2e, median p8 %.2e\n",
              rho, median(p4), median(p8)))
  cat(sprintf("   q<0.05 both (q4 genome-wide BH, q8 BH within 8-knot set): %d | q4<0.05: %d | q8<0.05: %d\n",
              both_q, sum(q4 < 0.05, na.rm = TRUE), sum(q8_within < 0.05, na.rm = TRUE)))

  rows[[name]] <- data.frame(branch = name, gene = g, pvalue_nknots4 = p4, pvalue_nknots8 = p8,
                             qval_nknots4 = q4, qval_nknots8_within = q8_within,
                             waldStat_nknots4 = a4[g, "waldStat_1"], waldStat_nknots8 = a8[g, "waldStat_1"])
}
if (length(rows)) write.csv(do.call(rbind, rows), out_csv, row.names = FALSE)
