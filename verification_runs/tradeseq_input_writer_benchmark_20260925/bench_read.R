# Read both writer variants the way read_tradeseq_inputs() does and check that
# R sees identical integer count matrices. Usage: Rscript bench_read.R OUT_DIR
args <- commandArgs(trailingOnly = TRUE)
out_root <- args[[1]]

read_counts <- function(path) {
  counts <- as.matrix(read.csv(gzfile(path), row.names = 1, check.names = FALSE))
  storage.mode(counts) <- "integer"
  counts
}

ok <- TRUE
for (name in c("counts_genes_x_cells.csv.gz", "counts_genes_x_cells.full.csv.gz")) {
  t_legacy <- system.time(a <- read_counts(file.path(out_root, "legacy_float", name)))[["elapsed"]]
  t_int <- system.time(b <- read_counts(file.path(out_root, "integer", name)))[["elapsed"]]
  same <- identical(a, b)
  ok <- ok && same
  cat(sprintf("%s: dim=%s identical=%s | R read legacy=%.1fs integer=%.1fs\n",
              name, paste(dim(a), collapse = "x"), same, t_legacy, t_int))
}
cat(if (ok) "R_READ_IDENTICAL\n" else "R_READ_MISMATCH\n")
if (!ok) quit(status = 1)
