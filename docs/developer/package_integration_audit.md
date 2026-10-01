# ChronoPack Package Integration / Compatibility Audit

Status: package-wide integration and compatibility audit complete. This pass
did not migrate benchmarking, remove compatibility layers wholesale, publish the
package, or redesign scientific methods.

## Scope Boundary

ChronoPack now owns reusable scientific workflows and standardized result
representations for:

- preprocessing;
- cell-state trajectories and fate;
- gene dynamics / tradeSeq;
- gene-expression trajectory characterization;
- dynamic-mode clustering;
- retained R plotting;
- cross-context mode conservation and switching;
- GRN inference;
- CellCommunication / C2C inference.

Benchmarking remains intentionally outside ChronoPack:

- no GRN benchmark namespace;
- no C2C benchmark namespace;
- no comparator methods, benchmark datasets, gold standards, benchmark metrics,
  benchmark-specific TF/LR universes, or benchmark-specific train/test logic.

## Architecture Overview

```text
AnnData
  ├── cp.io
  │     ├── data loading/conversion
  │     ├── preprocessing conventions
  │     └── preprocessing validation
  │
  ├── cp.tl cell/fate workflows
  │     ├── get_trajectories(...)
  │     ├── cluster_cell_state_trajectories(...)
  │     └── compute_terminal_fate_probabilities(...)
  │
  ├── cp.tl gene-dynamics workflows
  │     ├── fit_gene_dynamics(...)
  │     ├── characterize_gene_trajectories(...)
  │     └── cluster_gene_trajectories(...)
  │
  ├── cp.tl.grn
  │     ├── rank_regulators(...)
  │     ├── score_targets(...)
  │     └── infer(...) -> GRNResult
  │
  └── cp.tl.c2c
        └── infer(...) -> C2CResult

cross-context:
  cp.tl.mode_conservation(...) -> ModeConservationResult
  cp.tl.mode_switchers(...)    -> ModeConservationResult

visualization:
  cp.pl
    ├── gene trajectory / dynamic-mode plots
    ├── mode-conservation plots
    └── c2c plots
```

## Public API Inventory

### `cp.io`

- `read_h5ad(path, **kwargs)`
- `write_h5ad(adata, path, **kwargs)`
- preprocessing helpers:
  - `read_time_tagged_dge`
  - `add_time_metadata`
  - `filter_technical_genes`
  - `filter_low_count_genes`
  - `score_cell_cycle`
- validation:
  - `validate_adata`
- (removed 2026-09-30: `annotate_qc_metrics`, `filter_cells_by_qc`,
  `initialize_expression_layers`, `validate_preprocessed_adata`, and the
  `write_tradeseq_inputs(...)` compatibility wrapper)

### `cp.tl`

- cell-state / fate:
  - `get_trajectories`
  - `cluster_cell_state_trajectories`
  - `compute_terminal_fate_probabilities`
- gene dynamics:
  - `fit_gene_dynamics`
- gene-expression trajectories:
  - `characterize_gene_trajectories`
  - `cluster_gene_trajectories`
- cross-context mode analysis:
  - `mode_conservation`
  - `mode_switchers`
- subnamespaces:
  - `cp.tl.grn`
  - `cp.tl.c2c`

### `cp.pl`

- gene trajectory / dynamic modes:
  - `gene_trajectory`
  - `cluster_trajectories`
  - `trajectory_features`
- cross-context modes:
  - `mode_conservation`
  - `mode_sankey`
  - `mode_chord`
- C2C:
  - `cp.pl.c2c.interaction_scatter`
  - `cp.pl.c2c.pair_overlay`
  - `cp.pl.c2c.plot_genes_over_time`

### `cp.tl.grn`

- `rank_regulators(...) -> pd.DataFrame` (one context per call)
- `score_targets(...) -> TargetGeneScoreResult`
- `infer(...) -> GRNResult`
- compatibility:
  - `build_chronopack_ranked_edges`
  - `write_chronopack_ranked_edges`

### `cp.tl.c2c`

- `infer(...) -> C2CResult`
- helper exports used by the C2C notebook: `split_by_cell_type`,
  `build_pair_series` (plots in `cp.pl.c2c`).
  The `_legacy/` package was dissolved on 2026-09-30: its live code moved to
  `tl/c2c/_helpers.py`, `pl/c2c/_plots.py`, `tl/grn/_omnipath.py` and
  `tl/grn/_edges.py`.

## Standardized AnnData / Result Schema

### Preprocessing

Required by downstream workflows:

- unique `obs_names` and `var_names`;
- `adata.obs["time"]`;
- `adata.layers["counts"]`;
- workflow layers such as `cpm` and `log1p`;
- QC columns in `obs` and `var` where preprocessing validation requires them.

### Cell-State Trajectories / Fate

Conventions confirmed:

- per-cell trajectory/fate labels live in `adata.obs`;
- per-cell fate probability matrices live in `adata.obsm`;
- square cell-cell matrices live in `adata.obsp`;
- WOT/time-to-time or other non-axis-aligned state remains as result objects or
  sidecar files.

Representative post-GetTrajectories artifacts contain:

- `adata.obs["traj_clusters"]`;
- `adata.obsm["terminal_fate_probs_traj"]`;
- `adata.obsm["terminal_cell_dist"]`;
- square `adata.obsp["W_terminal_cosine_topk"]` and `adata.obsp["ot_graph"]`.

### Gene Dynamics

`cp.tl.fit_gene_dynamics(...)` writes scalar gene-level statistics to
`adata.var` under `key_added`, for example:

- `gene_dynamics_waldStat_1`
- `gene_dynamics_df_1`
- `gene_dynamics_pvalue_1`
- `gene_dynamics_meanLogFC`
- `gene_dynamics_qval_1`
- `gene_dynamics_significant`

Native tradeSeq state remains in sidecars for now.

### Gene Trajectory Characterization / Clustering

Standard keys:

- fitted trajectories:
  - `adata.varm["gene_trajectories_fit"]`
  - optional `adata.varm["gene_trajectories_fit_ddt"]`
- descriptors / association metadata:
  - `adata.var["gene_trajectories_*"]`
- PCA:
  - `adata.varm["gene_trajectory_clusters_X_pca"]`
- cluster/mode assignments:
  - `adata.var["gene_trajectory_clusters_*"]`
- provenance / time grid:
  - `adata.uns["chronopack"]["gene_trajectories"]`
  - `adata.uns["chronopack"]["gene_trajectory_clusters"]`

`cp.pl` no longer silently reconstructs missing fit matrices. It now emits a
clear error instructing users to regenerate or explicitly migrate legacy
artifacts if `adata.varm["gene_trajectories_fit"]` is absent.

### Cross-Context Modes

`ModeConservationResult` contains:

- `directional_scores`
- `conservation_scores`
- `gene_matches`
- `mode_transitions`
- `switchers`
- `parameters`

The result is not forced into either input AnnData.

### GRN

`GRNResult.edges` uses:

- `method`
- `regulator`
- `target`
- `score`
- `rank`
- `percentile`

Full networks remain in `GRNResult`, not AnnData.

### C2C

`C2CResult` contains:

- `interactions`
- `ligand_scores`
- `temporal_profiles`
- `candidate_interactions`
- detection vectors
- `parameters`
- `provenance`

The interaction network remains in `C2CResult`, not AnnData.

## Representative Integration Verification

Integration artifact:

- `verification_runs/package_integration_audit_20260918/integration_report.json`

This verification used frozen/historical artifacts and did not rerun expensive
or stochastic upstream stages.

Checks passed:

- preprocessing artifact:
  - `4034 x 18833`;
  - unique cell/gene indices;
  - `time`, `time_group`;
  - layers `counts`, `cpm`, `log1p`.
- post-GetTrajectories/fate artifacts:
  - K562+TNF and THP1+K562+LPS post-trajectory objects contain
    `traj_clusters`, terminal fate probabilities in `obsm`, and square graph
    matrices in `obsp`.
- gene dynamics artifact:
  - `4034 x 13905`;
  - expected `gene_dynamics_*` columns;
  - `.h5ad` round-trip succeeded.
- gene-trajectory / mode artifact:
  - after explicit audit-only migration of a stale artifact, plotting data
    extraction worked:
    - selected gene trajectories: 24 rows;
    - cluster trajectory summaries: 56 rows;
    - trajectory feature table: 500 rows;
  - `.h5ad` round-trip succeeded.
- mode conservation:
  - directional score table: `11 x 7`;
  - conservation score table: `6 x 3`;
  - switcher table: `96 x 10`;
  - transition table: `27 x 3`.
- GRN:
  - serialized clean public edge table has standardized columns;
  - metadata JSON contains `parameters` and `provenance`.
- C2C:
  - interactions: `189 x 11`;
  - temporal profiles: `2268 x 8`;
  - standardized context columns present.

## Inconsistencies Found / Fixed

Fixed now:

- `characterize_gene_trajectories(adata=None, ...)` had a stale variable
  reference in the standalone AnnData creation path. This was removed.
- `cp.pl.gene_trajectory(...)` now gives a more actionable missing-fit error
  for legacy artifacts lacking `adata.varm["gene_trajectories_fit"]`.

Documented compatibility issues:

- The frozen `gene_trajectory_clusters.h5ad` in
  `verification_runs/dimred_package_migration_20260918/output_pkg3/` is a
  legacy/stale artifact for the current plotting contract:
  - missing `adata.varm["gene_trajectories_fit"]`;
  - missing `adata.var["gene_trajectories_gene"]`;
  - missing `adata.var["gene_trajectories_cell_type"]`.
  The integration audit created an explicit migrated copy at
  `verification_runs/package_integration_audit_20260918/gene_trajectory_clusters.migrated_for_integration.h5ad`.
  Active release examples should use regenerated package-backed outputs or an
  explicit migration step, not silent plotting reconstruction.

- Older puromycin trajectory artifacts still contain legacy
  `traj_leiden_sr`/`fate_*_traj_last` columns. The active GetTrajectories
  notebooks use `traj_clusters`, and representative post-GetTrajectories
  artifacts for K562+TNF and THP1+K562+LPS already follow the newer contract.

## Legacy Shims / Sidecars

| Item | Classification | Keep for now? | Notes |
| --- | --- | --- | --- |
| `chronopack.{preprocessing,trajectories,tradeseq,grn_inference,omnipath_tf,c2c}` | compatibility shims | removed 2026-09-30 | last used by the frozen manuscript notebooks |
| `cp.tl.grn.build_chronopack_ranked_edges` | compatibility export | yes | reproduces old edge artifact, including explicit caller-supplied TF universe |
| `pc_df.new.csv/rds` | compatibility/export sidecar | yes | legacy downstream tables and regression comparison |
| `gene_trajectory_fit.csv`, `gene_trajectory_X_pca.csv` | generated sidecars | yes | useful for language-neutral export and artifact migration |
| native `gamList*.rds` / `sce*.rds` | generated/native sidecars | yes | retained R plotting/diagnostics and tradeSeq compatibility |
| association/significant CSVs | generated sidecars | yes | downstream GRN/C2C and legacy inspection |
| fixed tradeSeq input CSV/GZ files | generated sidecars | yes | HPC/R backend boundary and regression |
| C2C pair/ligand score CSVs | compatibility exports | yes | human inspection and existing notebook outputs |
| GRN historical ranked-edge CSV | compatibility export | yes | regression/legacy comparison only |
| BEELINE/benchmark notebooks | out of package scope | no package migration | keep outside ChronoPack APIs |

Deprecation candidates after release stabilization:

- top-level legacy imports once all active notebooks use `cp.io`, `cp.tl`, and
  `cp.pl`;
- old `pc_df` as the primary dynamic-mode representation once downstream
  notebooks use AnnData/result objects only;
- native tradeSeq sidecars if retained plotting/diagnostics no longer require
  them;
- benchmark-era compatibility edge wrappers once external benchmarking has its
  own maintained code path.

## Active Notebook Audit

All active notebooks import/use ChronoPack.

Remaining local functions are classified as notebook-local configuration or
presentation glue:

- path helpers:
  - `read_original_notebook_dir`
  - `original_project_dir`
  - `input_file`
- one-off cluster label mapping:
  - `RenameClusters`
  - `RenameCluster`
- manuscript/inspection selection helpers:
  - `top_rows_by_broad_cluster`
  - `fine_column_for_broad`
- THP1 preprocessing gene-list display helpers:
  - `_normalize_gene`
  - `_present_genes`

R TradeSeq notebooks contain R helper code such as `find_parent_with`; these are
thin R workflow clients and were not parsed as Python implementation code.

No active notebook was found defining reusable core methods equivalent to the
package-backed workflows.

## Remaining Package-Wide Issues Before Release Hardening

- Decide whether to provide an explicit one-time migration helper for legacy
  gene trajectory `.h5ad` artifacts that lack `gene_trajectories_fit` and
  standardized gene/context columns. Current behavior is intentionally strict.
- Decide when to deprecate top-level legacy modules. They are currently useful
  for old notebooks and regression scripts.
- Make release examples use regenerated package-backed artifacts rather than
  stale intermediate `.h5ad` files.
- CI/release hardening should separate lightweight tests from HPC/integration
  regression tests.
- Documentation should eventually include user-facing tutorials; current docs
  are developer/migration oriented.
