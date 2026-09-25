repo_root <- normalizePath("/new-stg/home/alex/tmp/ChronoPack", mustWork = TRUE)
input_dir <- file.path(repo_root, "verification_runs", "gene_dynamics_migration_20260918", "inputs_packaged_tradeseq_ce")

source(file.path(repo_root, "src", "chronopack", "_backends", "r", "tradeseq_helpers.R"))
load_chronopack_tradeseq_libraries()

result <- run_tradeseq_single_lineage(
  input_dir = input_dir,
  n_points = 9,
  workers = 4,
  screen_nknots = 4,
  final_nknots = 8,
  p_threshold = 0.05,
  output_association = "packaged_asso_df.1trajectory.nknots4.csv",
  output_significant = "packaged_sig_df.1trajectory.nknots4.csv",
  output_sce = "packaged_sce.1trajectory.nknots8.rds",
  output_gamlist = "packaged_gamList.1trajectory.nknots8.rds",
  verbose = TRUE
)

summary <- list(
  runtime = R.version.string,
  tradeSeq = as.character(utils::packageVersion("tradeSeq")),
  input_dir = input_dir,
  counts_dim = result$counts_dim,
  counts_sub_dim = result$counts_sub_dim,
  final_counts_dim = result$final_counts_dim,
  n_subset_genes = result$n_subset_genes,
  n_tested = result$n_tested,
  n_qval_significant = result$n_qval_significant,
  n_pvalue_significant = result$n_pvalue_significant,
  association_path = result$association_path,
  significant_path = result$significant_path,
  sce_path = result$sce_path,
  gamlist_path = result$gamlist_path
)

saveRDS(summary, file.path(input_dir, "packaged_run_summary.rds"))
print(summary)
