# Release Readiness Audit

This document records the first-release cleanup boundary after package-wide
integration testing.

## Public API Inventory

Recommended public imports:

- `chronopack.io`
  - `read_time_tagged_dge`, `load_time_tagged_anndata`
  - `filter_technical_genes`, `filter_low_count_genes`
  - `add_time_metadata`, `map_time_to_group`, `score_cell_cycle`
  - `validate`
  - (Standard QC, cell filtering and layer setup use Scanpy directly;
    `annotate_qc_metrics`, `filter_cells_by_qc`, `initialize_expression_layers`
    and `validate_preprocessed_adata` were removed on 2026-09-30.)
  - `read_h5ad`, `write_h5ad`
- `chronopack.tl`
  - `get_trajectories`
  - `cluster_cell_state_trajectories`
  - `compute_terminal_fate_probabilities`
  - `fit_gene_dynamics`
  - `characterize_gene_trajectories`
  - `cluster_gene_trajectories`
  - `mode_conservation`, `mode_switchers`
  - `CellStateTrajectoryResult`, `GeneDynamicsResult`
  - `GeneTrajectoryCharacterizationResult`
  - `GeneTrajectoryClusteringResult`
  - `ModeConservationResult`
- `chronopack.tl.grn`
  - `infer`, `rank_regulators`, `score_targets`
  - `GRNResult`, `TargetGeneScoreResult`
- `chronopack.tl.c2c`
  - `infer`
  - `C2CResult`
- `chronopack.pl`
  - `gene_trajectory`, `trajectory_features`,
    `cluster_medoids_and_top_genes`, `clustering_metrics`
  - `mode_conservation`, `mode_proportions`, `mode_sankey`, `mode_chord`,
    `cluster_pair_representative_genes`, `cluster_violin`
  - `RPlotResult`
- `chronopack.pl.c2c`
  - `interaction_scatter`, `pair_overlay`, `plot_genes_over_time`

Low-level helper functions may remain importable for compatibility, but they
are not part of the recommended public API unless listed above.

## Compatibility and Deprecation Policy

The former top-level compatibility modules (`chronopack.preprocessing`,
`chronopack.trajectories`, `chronopack.tradeseq`, `chronopack.grn_inference`,
`chronopack.omnipath_tf`, `chronopack.c2c`) and the flat source-checkout
`chronopack/` package were removed on 2026-09-30, together with the frozen
manuscript environment whose notebooks were their last users. Use
`chronopack.io`, `chronopack.tl` (including `cp.tl.grn` and `cp.tl.c2c`) and
`chronopack.pl`.

## Stale Artifact Resolution

`verification_runs/dimred_package_migration_20260918/output_pkg3/gene_trajectory_clusters.h5ad`
is a legacy verification artifact. It lacks the current gene-trajectory
contract fields:

- `adata.varm["gene_trajectories_fit"]`
- `adata.var["gene_trajectories_gene"]`
- `adata.var["gene_trajectories_cell_type"]`

Current plotting and downstream APIs should not silently reconstruct this
state. They should raise explicit validation errors and direct users to
regenerate or explicitly migrate the artifact.

The integration audit used:

`verification_runs/package_integration_audit_20260918/gene_trajectory_clusters.migrated_for_integration.h5ad`

as an explicitly migrated current-contract audit copy. This is a verification
artifact, not a package fixture.

## Sidecar and Artifact Inventory

| Artifact class | Status | Notes |
| --- | --- | --- |
| `gamList*.rds` | Compatibility-only / sidecar | Retained for historical notebooks and unmigrated custom R plotting. Core migrated plotting does not require it. |
| `sce*.rds` | Compatibility-only / sidecar | Retained for historical workflows; not stored inside AnnData. |
| association/significant CSVs | Required compatibility/export artifacts | Still consumed by some downstream restart and verification paths. |
| fixed tradeSeq inputs | Required runtime sidecars for fitting | Produced by the gene-dynamics workflow and retained for reproducibility. |
| `pc_df` CSV/RDS exports | Compatibility/export artifacts | AnnData `var`/`varm` are the primary package representation; pc_df export remains useful for older notebooks. |
| GRN compatibility CSVs | Compatibility-only artifacts | Clean GRN output is `GRNResult` plus language-neutral serialization. |
| C2C compatibility CSVs | Compatibility-only artifacts | Clean C2C output is `C2CResult`; CSV remains available for inspection. |
| `verification_runs/*` | Regression/verification artifacts | Excluded from package distributions. Do not rely on them as package fixtures. |
| external TF/LR/motif/contact resources | External scientific resources | Not bundled into ChronoPack unless separately approved and licensed. |

## Benchmarking Boundary

GRN and C2C benchmarking remain outside ChronoPack. The package provides
standardized inference outputs that external benchmarking workflows may consume,
but it does not include benchmark datasets, gold standards, comparator methods,
benchmark metrics, or benchmark-specific plotting APIs.

## Release Blockers to Review

- Confirm whether the proprietary placeholder license should be replaced before
  public distribution.
- Decide later whether clustering-sensitive runtime dependencies should be
  pinned more tightly after the known Leiden/graph drift investigation.
- Decide later when compatibility modules and native sidecars can be removed.
