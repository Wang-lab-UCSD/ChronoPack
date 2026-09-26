# What-if test (maintained run): DimRed gene clustering with scaling computed
# within each partition (no reference) instead of over cluster_universe; alpha = 0.1 and k rules as the
# notebook: broad k forced to 3, fine k = argmax elbow_score. Diagnostic only.
suppressPackageStartupMessages(library(matrixStats))
root <- "/new-stg/home/alex/tmp/ChronoPack"
source(file.path(root, "chronopack", "R", "dimred_clustering_helpers.R"))

vr <- file.path(root, "verification_runs")
out_root <- file.path(vr, "dimred_partition_scaling_test_20260925")
run <- file.path(vr, "full_pipeline_maintained_chronopack_20260924_r7",
                 "k562-tnfa-puromycin.hg38/self_contained_outputs/global_analysis_extended")
pc <- read.csv(file.path(run, "pc_df.new.csv"), row.names = 1, check.names = FALSE)

ari <- function(a, b) {
  t <- table(a, b); s <- function(x) sum(choose(x, 2))
  e <- s(rowSums(t)) * s(colSums(t)) / choose(sum(t), 2)
  (s(t) - e) / ((s(rowSums(t)) + s(colSums(t))) / 2 - e)
}
read_assign <- function(d) as.matrix(read.csv(file.path(d, "cluster_assignments.csv"), row.names = 1, check.names = FALSE))
fine_k <- function(d) { m <- read.csv(file.path(d, "cluster_metrics.csv")); kc <- if ("k" %in% names(m)) "k" else "K"; m[[kc]][which.max(m$elbow_score)] }
col <- function(k) paste0("k", k)

rows <- list()
for (part in c("induced", "repressed")) {
  d <- file.path(run, "cluster_candidates", part)
  out <- file.path(out_root, part)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  chronopack_cluster_trajectory_matrix_to_files(
    traj_mat_path = file.path(d, "cluster_input.fit.csv"),
    time_grid_path = file.path(d, "cluster_input.time.csv"),
    output_dir = out, alpha = 0.1, K_grid = 2:100,
    reference_traj_mat_path = NULL
  )
  a10 <- read_assign(d); a25 <- read_assign(out)
  g <- rownames(a10); stopifnot(setequal(g, rownames(a25)))
  k10 <- fine_k(d); k25 <- fine_k(out)

  # Broad: compare against the named broad labels from the alpha = 0.1 run.
  named <- pc[g, paste0("cluster_broad_", part)]
  b25 <- a25[g, col(3)]
  tab <- table(alpha_0.1 = named, per_partition_cluster = b25)
  best <- apply(tab, 2, function(x) rownames(tab)[which.max(x)])  # majority name per new cluster
  moved <- sum(named != best[as.character(b25)])
  cat(sprintf("\n==== %s (%d genes)\n", part, length(g)))
  cat("Broad k=3 contingency (rows: current named clusters, cols: per-partition scaling clusters):\n")
  print(tab)
  cat(sprintf("Broad: ARI %.4f | genes changing broad cluster (after matching clusters by majority): %d (%.1f%%)\n",
              ari(named, b25), moved, 100 * moved / length(g)))
  cat(sprintf("Broad sizes alpha=0.1: %s | per-partition: %s\n",
              paste(sort(as.integer(table(named)), decreasing = TRUE), collapse = "/"),
              paste(sort(as.integer(table(b25)), decreasing = TRUE), collapse = "/")))
  cat(sprintf("Fine k by elbow: alpha=0.1 -> %d, per-partition -> %d\n", k10, k25))
  cat(sprintf("Fine ARI: at each alpha's own k %.4f | both at k=%d %.4f | both at k=%d %.4f\n",
              ari(a10[g, col(k10)], a25[g, col(k25)]), k10, ari(a10[g, col(k10)], a25[g, col(k10)]),
              k25, ari(a10[g, col(k25)], a25[g, col(k25)])))
  rows[[part]] <- data.frame(part = part, genes = length(g), broad_ARI = ari(named, b25),
                             broad_genes_moved = moved, fine_k_alpha0.1 = k10, fine_k_per_partition = k25,
                             fine_ARI_own_k = ari(a10[g, col(k10)], a25[g, col(k25)]),
                             fine_ARI_same_k = ari(a10[g, col(k10)], a25[g, col(k10)]))
}
write.csv(do.call(rbind, rows), file.path(out_root, "summary.csv"), row.names = FALSE)
cat("TEST_OK\n")
