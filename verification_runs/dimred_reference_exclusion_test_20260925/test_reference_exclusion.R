# What-if test: DimRed gene clustering with K562_puro rows removed from the
# scaling reference (reference = cluster_universe instead of sig). Uses the
# package's own clustering function and the notebook's k-selection rules.
# Diagnostic only; nothing in the pipeline runs is modified.
suppressPackageStartupMessages(library(matrixStats))
root <- "/new-stg/home/alex/tmp/ChronoPack"
source(file.path(root, "chronopack", "R", "dimred_clustering_helpers.R"))

vr <- file.path(root, "verification_runs")
out_root <- file.path(vr, "dimred_reference_exclusion_test_20260925")
sub <- "k562-tnfa-puromycin.hg38/self_contained_outputs/global_analysis_extended/cluster_candidates"
runs <- c(frozen = "full_pipeline_frozen_manuscript_20260924_r6",
          maintained = "full_pipeline_maintained_chronopack_20260924_r7")
template <- file.path(vr, "manuscript_full_pipeline_20260919_165539", sub)

ari <- function(a, b) {
  t <- table(a, b); s <- function(x) sum(choose(x, 2))
  e <- s(rowSums(t)) * s(colSums(t)) / choose(sum(t), 2)
  (s(t) - e) / ((s(rowSums(t)) + s(colSums(t))) / 2 - e)
}
read_assign <- function(d) as.matrix(read.csv(file.path(d, "cluster_assignments.csv"), row.names = 1, check.names = FALSE))
fine_k <- function(d) { m <- read.csv(file.path(d, "cluster_metrics.csv")); kc <- if ("k" %in% names(m)) "k" else "K"; m[[kc]][which.max(m$elbow_score)] }

results <- list()
for (run in names(runs)) {
  for (part in c("induced", "repressed")) {
    d <- file.path(vr, runs[[run]], sub, part)
    ref <- read.csv(file.path(d, "cluster_reference.fit.csv"), row.names = 1, check.names = FALSE)
    ref_np <- ref[!grepl("_K562_puro$", rownames(ref)), , drop = FALSE]
    out <- file.path(out_root, run, part)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    ref_path <- file.path(out, "reference_without_puro.fit.csv")
    write.csv(ref_np, ref_path)
    chronopack_cluster_trajectory_matrix_to_files(
      traj_mat_path = file.path(d, "cluster_input.fit.csv"),
      time_grid_path = file.path(d, "cluster_input.time.csv"),
      output_dir = out, alpha = 0.1, K_grid = 2:100,
      reference_traj_mat_path = ref_path
    )
    new <- read_assign(out); old <- read_assign(file.path(template, part)); own <- read_assign(d)
    g <- rownames(old)
    stopifnot(setequal(g, rownames(new)))
    kf_new <- fine_k(out); kf_old <- fine_k(file.path(template, part)); kf_own <- fine_k(d)
    col <- function(k) paste0("k", k)
    b_new <- new[g, col(3)]; b_old <- old[g, col(3)]
    results[[paste(run, part)]] <- data.frame(
      run = run, part = part, reference_rows = nrow(ref_np), genes = length(g),
      broad_ARI_vs_template = ari(b_new, b_old),
      broad_sizes_new = paste(sort(as.integer(table(b_new)), decreasing = TRUE), collapse = "/"),
      broad_sizes_template = paste(sort(as.integer(table(b_old)), decreasing = TRUE), collapse = "/"),
      fine_k_new = kf_new, fine_k_template = kf_old, fine_k_current_run = kf_own,
      fine_ARI_vs_template = ari(new[g, col(kf_new)], old[g, col(kf_old)]),
      fine_ARI_same_k_as_template = ari(new[g, col(kf_old)], old[g, col(kf_old)]),
      fine_ARI_current_run_vs_template = ari(own[g, col(kf_own)], old[g, col(kf_old)])
    )
  }
}
res <- do.call(rbind, results)
rownames(res) <- NULL
write.csv(res, file.path(out_root, "summary.csv"), row.names = FALSE)
print(res, digits = 4)
f <- file.path(out_root, c("frozen", "maintained"), "repressed", "cluster_assignments.csv")
cat("\nfrozen vs maintained assignments identical (repressed):",
    identical(read_assign(dirname(f[1])), read_assign(dirname(f[2]))), "\n")
f <- file.path(out_root, c("frozen", "maintained"), "induced", "cluster_assignments.csv")
cat("frozen vs maintained assignments identical (induced):",
    identical(read_assign(dirname(f[1])), read_assign(dirname(f[2]))), "\n")
cat("TEST_OK\n")
