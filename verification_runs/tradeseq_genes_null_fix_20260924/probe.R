# Probe for the fit_tradeseq_sce() genes = NULL fix.
#
# Fits a small gene subset through the patched helper (cell-level vector
# offset, sce = TRUE, parallel) and checks it against a direct tradeSeq fit
# using the historical matrix offset. The matrix fit is a parity reference
# only; the pipeline representation stays cell-level.
options(error = function() { traceback(2); quit(status = 1) })

args <- commandArgs(trailingOnly = TRUE)
repo_root <- args[[1]]
input_dir <- args[[2]]
n_genes <- if (length(args) >= 3) as.integer(args[[3]]) else 30L

source(file.path(repo_root, "chronopack", "R", "tradeseq_helpers.R"))
load_chronopack_tradeseq_libraries()
cat("R:", R.version.string, "| tradeSeq:", as.character(packageVersion("tradeSeq")),
    "| mgcv:", as.character(packageVersion("mgcv")), "\n")

old_wd <- getwd()
setwd(input_dir)
inputs <- read_tradeseq_inputs(tag = "")
setwd(old_wd)

counts_small <- inputs$counts_sub[seq_len(min(n_genes, nrow(inputs$counts_sub))), , drop = FALSE]
lineage <- prepare_single_lineage(counts_small, inputs$time_vec, inputs$sf)
cat("counts_sub:", dim(inputs$counts_sub), "| probe subset:", dim(counts_small), "\n")

check <- function(ok, msg) {
  cat(if (ok) "PASS" else "FAIL", "-", msg, "\n")
  if (!ok) quit(status = 1)
}

# 1. Pipeline path: genes = NULL (default), vector offset, parallel.
set.seed(123)
sce_vec <- fit_tradeseq_sce(
  counts_sub = counts_small,
  pseudotime = lineage$PT,
  cell_weights = lineage$CW,
  offset = lineage$offset_vec,
  nknots = 4,
  workers = 4,
  verbose = FALSE
)
check(is(sce_vec, "SingleCellExperiment"), "helper returns a SingleCellExperiment")
check(nrow(sce_vec) == nrow(counts_small), sprintf("all %d genes fitted", nrow(counts_small)))
check(identical(rownames(sce_vec), rownames(counts_small)), "gene order preserved")
beta_vec <- as.matrix(rowData(sce_vec)$tradeSeq$beta[[1]])
cat("genes with non-NA coefficients:", sum(stats::complete.cases(beta_vec)), "/", nrow(beta_vec), "\n")
check(sum(stats::complete.cases(beta_vec)) > 0, "at least one gene converged to coefficients")

# 2. Explicit character gene subset still works through the helper.
sub_genes <- rownames(counts_small)[1:5]
set.seed(123)
sce_sub <- fit_tradeseq_sce(
  counts_sub = counts_small,
  pseudotime = lineage$PT,
  cell_weights = lineage$CW,
  offset = lineage$offset_vec,
  nknots = 4,
  workers = 4,
  verbose = FALSE,
  genes = sub_genes
)
check(identical(rownames(sce_sub), sub_genes), "explicit genes subset fitted")

# 3. Association test on the fitted SCE (the next pipeline step).
asso_df <- association_df_single(sce_vec, n_points = 9)
check(all(c("waldStat_1", "pvalue_1", "meanLogFC", "qval_1") %in% names(asso_df)),
      "associationTest columns present")
sce_vec <- annotate_sce_association(sce_vec, asso_df)
check(all(c("waldStat_1", "meanLogFC", "pvalue_1", "qval_1") %in% colnames(rowData(sce_vec))),
      "association fields written to rowData")

# 4. Parity reference: historical matrix offset through tradeSeq directly.
set.seed(123)
sce_mat <- tradeSeq::fitGAM(
  counts = counts_small,
  pseudotime = lineage$PT,
  cellWeights = lineage$CW,
  offset = make_offset_matrix(counts_small, inputs$sf),
  nknots = 4,
  sce = TRUE,
  parallel = TRUE,
  BPPARAM = BiocParallel::MulticoreParam(workers = 4),
  verbose = FALSE
)
beta_mat <- as.matrix(rowData(sce_mat)$tradeSeq$beta[[1]])
check(identical(is.na(beta_vec), is.na(beta_mat)), "same genes fail/converge under both offsets")
max_beta_diff <- max(abs(beta_vec - beta_mat), na.rm = TRUE)
cat("max |beta_vector - beta_matrix|:", format(max_beta_diff, digits = 6), "\n")
check(max_beta_diff < 1e-8, "vector-offset coefficients match matrix-offset coefficients")

cat("ALL CHECKS PASSED\n")
