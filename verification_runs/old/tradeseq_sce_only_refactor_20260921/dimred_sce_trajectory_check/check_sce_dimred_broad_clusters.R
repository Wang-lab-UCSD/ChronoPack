suppressPackageStartupMessages({
  library(dplyr)
})
repo <- '/new-stg/home/alex/tmp/ChronoPack'
source(file.path(repo, 'src/chronopack/_backends/r/dimred_clustering_helpers.R'))
source_root <- file.path(repo, 'verification_runs', 'filename_migration_20260917_074750')
outdir <- file.path(repo, 'verification_runs', 'tradeseq_sce_only_refactor_20260921', 'dimred_sce_trajectory_check')
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cfg <- chronopack_dimred_config(
  root_puro = file.path(source_root, 'k562-tnfa-puromycin.hg38'),
  root_k562a = file.path(source_root, 'k562-tnfa.hg38.rep'),
  root_lps = file.path(source_root, 'thp1-k562-lps.hg38.rep')
)
keys_dimred <- c('K562a', 'K562_puro', 'K562b', 'THP1')
loaded <- load_dimred_fixed_tradeSeq_outputs(cfg)
feature_info <- select_dimred_features(
  asso_df.ls = loaded$asso_df.ls,
  sce.ls = loaded$sce.ls,
  counts_sub.ls = loaded$counts_sub.ls,
  keys = keys_dimred,
  pvalue_threshold = 0.05
)

t_interval_map <- c(K562a = 10, K562_puro = 10, K562b = 15, THP1 = 15)
traj.ls <- extract_dimred_smooth_trajectories(
  sce.ls = loaded$sce.ls,
  sf.ls = loaded$sf.ls,
  PT.ls = loaded$PT.ls,
  features = feature_info$usable_features.ls,
  keys = keys_dimred,
  t_interval_map = t_interval_map
)

align_res <- align_rescale_dimred_trajectories(
  traj.ls = traj.ls,
  time_map = c(K562a = 110, K562_puro = 110, K562b = 165, THP1 = 165),
  offset_map = c(K562a = 0, K562_puro = 0, K562b = -68, THP1 = -10),
  scale_map = c(THP1 = 0.65)
)
traj_res <- construct_dimred_traj_matrix(align_res$traj_trunc_rescaled.ls)
traj_mat <- traj_res$traj_mat
stt <- traj_res$tt
if (ncol(traj_mat) != 97L) {
  stop("Expected 97 one-minute trajectory timepoints, got ", ncol(traj_mat), ".")
}

pc_df <- build_dimred_pc_df(
  traj_mat = traj_mat,
  tt = stt,
  asso_df.ls = loaded$asso_df.ls,
  scale_map_use = align_res$scale_map_use,
  descriptor_args = list(
    N_minutes = 10 / length(stt),
    z_diff_threshold = 0.2,
    require_consecutive = 30,
    onset_mode = 'horizontal',
    pseudocount = 0
  )
)
matrix_res <- prepare_dimred_matrices(traj_mat = traj_mat, tt = stt, pc_df = pc_df)
pca_res <- compute_dimred_pca(matrix_res$X_mat, matrix_res$Xd_mat, variance_threshold = 0.99)
pc_df <- add_dimred_pca_columns(pc_df, pca_res$pca_X, n_columns = 4)
clustering_inputs <- prepare_dimred_clustering_inputs(matrix_res$X_mat, matrix_res$Xd_mat)
induced_res <- run_dimred_signed_clustering(
  pc_df = pc_df, X_n = clustering_inputs$X_n, Xd_n = clustering_inputs$Xd_n,
  sign = 'induced', alpha = 0.25, log2fc_threshold = 0,
  K_grid = 2:100, wss_fraction = 0.95, select_k_by = 'elbow'
)
repressed_res <- run_dimred_signed_clustering(
  pc_df = pc_df, X_n = clustering_inputs$X_n, Xd_n = clustering_inputs$Xd_n,
  sign = 'repressed', alpha = 0.25, log2fc_threshold = 0,
  K_grid = 2:100, wss_fraction = 0.95, select_k_by = 'elbow'
)
pc_df <- add_dimred_signed_clusters(pc_df, induced_res, repressed_res)
pc_df <- finalize_dimred_pc_df(pc_df)

summary <- list(
  pc_df_dim = dim(pc_df),
  traj_dim = dim(traj_mat),
  feature_counts = vapply(feature_info$usable_features.ls, length, integer(1)),
  induced_selected_k = induced_res$clust_eval$selection$k_selected,
  induced_k_elbow = induced_res$clust_eval$selection$k_elbow,
  induced_k_silhouette = induced_res$clust_eval$selection$k_silhouette,
  repressed_selected_k = repressed_res$clust_eval$selection$k_selected,
  repressed_k_elbow = repressed_res$clust_eval$selection$k_elbow,
  repressed_k_silhouette = repressed_res$clust_eval$selection$k_silhouette,
  induced_cluster_table = table(pc_df$cluster_broad_induced),
  repressed_cluster_table = table(pc_df$cluster_broad_repressed),
  broad_cluster_table = table(pc_df$broad_cluster)
)
print(summary)
write.csv(pc_df, file.path(outdir, 'pc_df.sce_trajectory.csv'))
saveRDS(pc_df, file.path(outdir, 'pc_df.sce_trajectory.rds'))
write.csv(induced_res$clust_eval$metrics, file.path(outdir, 'cluster_metrics.induced.csv'), row.names = FALSE)
write.csv(repressed_res$clust_eval$metrics, file.path(outdir, 'cluster_metrics.repressed.csv'), row.names = FALSE)
write.csv(traj_mat, file.path(outdir, 'traj_mat.sce_trajectory.csv'))
