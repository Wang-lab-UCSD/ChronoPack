# Gene-Trajectory Characterization And Dynamic-Mode Migration

Date: 2026-09-18

This slice migrates the computational portion of the puromycin
`DimensionalityReduction&Clustering` workflow. It does not migrate
publication-quality plotting, mode-conservation analysis, GRN, or C2C.

## Public API

Implemented public task functions:

```python
cp.tl.characterize_gene_trajectories(...)
cp.tl.cluster_gene_trajectories(...)
```

The implementation preserves the established R computational backend and uses
Python as the orchestration / AnnData-standardization layer.

## Workflow Mapping

| Current code / stage | ChronoPack location | Public vs internal | Output |
| --- | --- | --- | --- |
| Load fixed tradeSeq sidecars | `chronopack_dimred_config()`, `load_dimred_fixed_tradeSeq_outputs()` | Internal R backend | native R lists |
| Select significant genes / available fits | `select_dimred_features()` | Internal R backend | feature sets |
| Extract fitted trajectories from gamList/SCE | `.get_traj_from_predictSmooth()` through `chronopack_characterize_gene_trajectories()` | Internal R backend | standardized trajectory matrix |
| Align/rescale time grids | `align_rescale_dimred_trajectories()` | Internal R backend | aligned fitted curves |
| Build gene-context trajectory matrix | `construct_dimred_traj_matrix()` | Internal R backend | `traj_mat`, `tt` |
| Descriptors: onset, peak, log2FC, curvature, turning, OWP, owp_ddt | `get_descriptors()` and `build_dimred_pc_df()` | Internal R backend, exposed through `cp.tl.characterize_gene_trajectories()` | `adata.var`, `pc_df.characterized.csv` |
| Z-scored trajectories and derivatives | `prepare_dimred_matrices()` | Internal R backend | `adata.varm`, CSV sidecars |
| PCA | `compute_dimred_pca()` | Internal R backend | `adata.varm[f"{key}_X_pca"]`, PCA CSV |
| Signed Ward.D2 clustering | `run_dimred_signed_clustering()` | Internal R backend, exposed through `cp.tl.cluster_gene_trajectories()` | cluster columns in `adata.var` |
| Cluster metrics | `calculate_trajectory_clustering_metrics()` | Internal R backend | result metrics / CSV sidecars |
| Broad-mode label mapping | `rename_dimred_broad_clusters()` and `finalize_dimred_pc_df()` | Internal R backend | `broad_cluster` |
| Historical `pc_df` construction/export | package result + sidecar helper files | Compatibility | `pc_df.new.csv`, `pc_df.new.rds` |

## AnnData Contract

`cp.tl.characterize_gene_trajectories()` writes:

- `adata.var[f"{key}_gene"]`
- `adata.var[f"{key}_cell_type"]`
- `adata.var[f"{key}_t_onset"]`
- `adata.var[f"{key}_t_peak"]`
- `adata.var[f"{key}_baseline_exp"]`
- `adata.var[f"{key}_peak_exp"]`
- `adata.var[f"{key}_peak_dlogexp"]`
- `adata.var[f"{key}_peak_dexp"]`
- `adata.var[f"{key}_peak_log2fc"]`
- `adata.var[f"{key}_auc_norm"]`
- `adata.var[f"{key}_owp"]`
- `adata.var[f"{key}_owp_ddt"]`
- `adata.var[f"{key}_mean_curv"]`
- `adata.var[f"{key}_max_curv"]`
- `adata.var[f"{key}_mean_turn"]`
- `adata.var[f"{key}_max_turn"]`
- `adata.varm[f"{key}_fit"]`
- `adata.varm[f"{key}_fit_ddt"]`

The time grid is recorded in ChronoPack provenance under
`adata.uns["chronopack"][key]`.

`cp.tl.cluster_gene_trajectories()` writes:

- cluster and dynamic-mode columns in `adata.var`, namespaced by `key_added`;
- `adata.varm[f"{key}_X_pca"]`;
- sidecar paths and provenance in `adata.uns["chronopack"][key]`.

The current fitted-curve representation uses gene-context rows such as
`NFKBIA_K562a`; these map directly to historical `pc_df` rows.

## Result Objects

`GeneTrajectoryCharacterizationResult` contains:

- `adata`;
- `pc_df`;
- `time_grid`;
- sidecar paths;
- the native backend result object.

`GeneTrajectoryClusteringResult` contains:

- `adata`;
- final `pc_df`;
- cluster metrics;
- sidecar paths;
- the native backend result object.

Native SCE/gamList sidecars are still required upstream because retained R
plotting functions have not yet been migrated.

## Regression Results

Verification directory:

- `verification_runs/dimred_package_migration_20260918/`

The package-backed task functions were compared against an old-style stepwise
execution of the active helper sequence on the same fixed saved tradeSeq
sidecars. No tradeSeq fits were regenerated.

Package task output versus stepwise output:

- `pc_df.new.csv`: same shape `500 x 35`;
- row order: identical;
- column order: identical;
- numeric columns: maximum absolute difference `0` for every compared column;
- categorical columns: match fraction `1.0`;
- standardized fitted trajectory matrix: same shape `500 x 8`;
- fitted trajectory values: maximum absolute difference `0`.

The current saved fixed sidecars used for this regression contain only
per-experiment fitted genes. Because the active DimRed feature-selection logic
requires fitted genes to be available across contexts, this current sidecar set
produces `500` gene-context rows. The older extended reference produced `6753`
rows because it used larger union/extension gamLists. This is a sidecar/input
universe difference, not a package wrapper difference.

## Runtime Notes

The unified R runtime did not expose the `cluster` package to the R backend.
The previous `cluster::silhouette()` dependency was replaced with an internal
implementation of the same silhouette-width definition using the distance
matrix and hard labels. This removes a fragile namespace dependency without
changing the clustering method.

## Out Of Scope

- Publication-quality R plotting wrappers.
- Mode-conservation migration.
- GRN.
- C2C.
