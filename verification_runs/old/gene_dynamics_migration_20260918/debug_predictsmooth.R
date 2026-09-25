suppressPackageStartupMessages({
  library(tradeSeq)
})

left <- "/new-stg/home/alex/tmp/ChronoPack/verification_runs/filename_migration_20260917_074750/k562-tnfa-puromycin.hg38/intmd/sc-K562-1e6-to-tradeSeq.fixed/gamList.1trajectory.nknots8.rds"
right <- "/new-stg/home/alex/tmp/ChronoPack/verification_runs/gene_dynamics_migration_20260918/inputs_packaged_tradeseq_ce/packaged_gamList.1trajectory.nknots8.rds"

g1 <- readRDS(left)
g2 <- readRDS(right)
cat("lengths", length(g1), length(g2), "\n")
cat("names common", length(intersect(names(g1), names(g2))), "\n")
genes <- intersect(names(g1), names(g2))
print(head(genes))

res <- try(tradeSeq::predictSmooth(g1, gene = genes, nPoints = 9, tidy = TRUE), silent = TRUE)
print(class(res))
if (inherits(res, "try-error")) {
  cat(as.character(res), "\n")
} else {
  print(dim(res))
  print(names(res))
  print(head(res))
}
