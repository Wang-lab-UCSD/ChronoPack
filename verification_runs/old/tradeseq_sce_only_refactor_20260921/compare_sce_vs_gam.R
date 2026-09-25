suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(tradeSeq)
})

repo <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(repo, "src/chronopack/_backends/r/dimred_clustering_helpers.R"))

run_root <- file.path(repo, "verification_runs/filename_migration_20260917_074750")
output_dir <- file.path(repo, "verification_runs/tradeseq_sce_only_refactor_20260921")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cfg <- chronopack_dimred_config(
  root_puro = file.path(run_root, "k562-tnfa-puromycin.hg38"),
  root_k562a = file.path(run_root, "k562-tnfa.hg38.rep"),
  root_lps = file.path(run_root, "thp1-k562-lps.hg38.rep")
)

keys <- c("K562a", "K562_puro", "K562b", "THP1")
t_interval_map <- c(K562a = 10, K562_puro = 10, K562b = 15, THP1 = 15)

loaded <- load_dimred_fixed_tradeSeq_outputs(cfg)
feature_info <- select_dimred_features(
  asso_df.ls = loaded$asso_df.ls,
  sce.ls = loaded$sce.ls,
  counts_sub.ls = loaded$counts_sub.ls,
  keys = keys,
  pvalue_threshold = 0.05
)

prediction_rows <- list()
trajectory_rows <- list()

for (key in keys) {
  message("Checking ", key)
  sce <- loaded$sce.ls[[key]]
  gamList <- loaded$gamList.ls[[key]]
  sf <- loaded$sf.ls[[key]]
  PT <- loaded$PT.ls[[key]]
  features <- feature_info$usable_features.ls[[key]]
  nPoints <- 1 + ((length(unique(PT)) - 1) * t_interval_map[[key]])
  off_ref <- compute_tradeSeq_reference_offset(sf)

  panel_genes <- unique(c(
    head(features, 3),
    tail(features, 2)
  ))
  panel_genes <- intersect(panel_genes, intersect(rownames(sce), names(gamList)))

  dm <- SummarizedExperiment::colData(sce)$tradeSeq$dm
  df_pred <- tradeSeq:::.getPredictRangeDf(dm, lineageId = 1, nPoints = nPoints)

  for (gene in panel_genes) {
    old <- predict(
      gamList[[gene]],
      newdata = df_pred,
      type = "link",
      se.fit = TRUE,
      offset = rep(off_ref, nrow(df_pred))
    )
    new <- predict_tradeSeq_sce_with_se(
      sce = sce,
      gene = gene,
      lineageId = 1,
      nPoints = nPoints,
      off_ref = off_ref,
      se_factor = 1.96,
      return_log1p = TRUE
    )
    old_mu <- log1p(exp(old$fit))
    old_lo <- log1p(exp(old$fit - 1.96 * old$se.fit))
    old_hi <- log1p(exp(old$fit + 1.96 * old$se.fit))
    prediction_rows[[length(prediction_rows) + 1L]] <- data.frame(
      key = key,
      gene = gene,
      n_points = nPoints,
      max_abs_fit = max(abs(as.numeric(old$fit) - new$fit)),
      max_abs_se = max(abs(as.numeric(old$se.fit) - new$se.fit)),
      max_abs_mu = max(abs(old_mu - new$mu)),
      max_abs_lo = max(abs(old_lo - new$lo)),
      max_abs_hi = max(abs(old_hi - new$hi)),
      stringsAsFactors = FALSE
    )
  }

  old_traj <- .get_traj_from_predictSmooth(
    gamList = gamList,
    features = features,
    nPoints = nPoints,
    lineage_pick = 1,
    verbose = FALSE,
    sce = sce,
    sizeFactors = sf,
    prefer_predict_gam = TRUE
  )
  new_traj <- get_traj_from_tradeSeq_sce(
    sce = sce,
    features = features,
    nPoints = nPoints,
    lineage_pick = 1,
    sizeFactors = sf,
    return_log1p = TRUE,
    verbose = FALSE
  )

  common_rows <- intersect(rownames(old_traj$traj_mat), rownames(new_traj$traj_mat))
  old_mat <- old_traj$traj_mat[common_rows, , drop = FALSE]
  new_mat <- new_traj$traj_mat[common_rows, , drop = FALSE]
  trajectory_rows[[length(trajectory_rows) + 1L]] <- data.frame(
    key = key,
    n_features_requested = length(features),
    old_nrow = nrow(old_traj$traj_mat),
    new_nrow = nrow(new_traj$traj_mat),
    common_nrow = length(common_rows),
    ncol_old = ncol(old_traj$traj_mat),
    ncol_new = ncol(new_traj$traj_mat),
    same_row_order = identical(rownames(old_mat), rownames(new_mat)),
    same_time_grid = identical(as.numeric(old_traj$tt), as.numeric(new_traj$tt)),
    max_abs_traj = max(abs(old_mat - new_mat)),
    stringsAsFactors = FALSE
  )
}

prediction_df <- do.call(rbind, prediction_rows)
trajectory_df <- do.call(rbind, trajectory_rows)

write.csv(prediction_df, file.path(output_dir, "sce_vs_gam_prediction_parity.csv"), row.names = FALSE)
write.csv(trajectory_df, file.path(output_dir, "sce_vs_gam_trajectory_matrix_parity.csv"), row.names = FALSE)

print(prediction_df)
print(trajectory_df)

if (any(prediction_df$max_abs_fit > 1e-8) ||
    any(prediction_df$max_abs_se > 1e-8) ||
    any(prediction_df$max_abs_mu > 1e-8) ||
    any(prediction_df$max_abs_lo > 1e-8) ||
    any(prediction_df$max_abs_hi > 1e-8)) {
  stop("SCE-vs-gam prediction parity failed.")
}

if (any(!trajectory_df$same_row_order) ||
    any(!trajectory_df$same_time_grid) ||
    any(trajectory_df$old_nrow != trajectory_df$new_nrow) ||
    any(trajectory_df$max_abs_traj > 1e-8)) {
  stop("SCE-vs-gam trajectory matrix parity failed.")
}

message("SCE-vs-gam parity passed.")
