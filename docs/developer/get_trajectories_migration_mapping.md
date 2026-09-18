# GetTrajectories Migration Mapping

Date: 2026-09-18

This migration covers the cell-state / fate trajectory workflow. It does not
cover fitted gene-expression trajectories, gene-trajectory clustering, or
dynamic modes.

## Current Code To Package Mapping

| Current code / stage | ChronoPack location | Public or internal | Output location |
| --- | --- | --- | --- |
| WaddingtonOT model setup and `compute_all_transport_maps()` call | `cp.tl.get_trajectories()` via `run_wot_transport_workflow()` | Public task API plus internal helper | WOT files remain sidecar files under configured `tmap_dir` |
| WOT h5ad transport-map loading (`read_tmap`) | `chronopack.tl._trajectories.read_tmap()` | Internal helper; compatibility exported during migration | `CellStateTrajectoryResult.tmaps` |
| Long-form WOT coupling extraction (`extract_wot_couplings`) | `chronopack.tl._trajectories.extract_wot_couplings()` | Internal helper; compatibility exported during migration | `CellStateTrajectoryResult.transitions_long` |
| Rectangular time-to-time matrices (`build_transition_matrices`) | `chronopack.tl._trajectories.build_transition_matrices()` | Internal helper; compatibility exported during migration | `CellStateTrajectoryResult.tmats` |
| Square cell-by-cell transition graph construction | `chronopack.tl._trajectories.sparse_adjacency_from_tmats()` | Internal helper | `CellStateTrajectoryResult.T_full`, `T_pruned`, `T_sym`; optionally `adata.obsp['ot_graph']` |
| THP1/K562 block embedding | `chronopack.tl._trajectories.embed_block()` | Internal helper; compatibility exported during migration | Dataset-specific notebook combines blocks before storing square graph |
| UMAP from WOT graph | `cp.tl.get_trajectories()` / Scanpy | Public task API delegates to Scanpy | `adata.obsm['X_ot_umap']`, and Scanpy's `adata.obsm['X_umap']` |
| Terminal-cell fate distributions | `cp.tl.cluster_cell_state_trajectories()` through `CalculateTerminalFatesAcrossCells()` | Public task helper plus internal implementation | `adata.obsm['terminal_cell_dist']`, `adata.uns['terminal_cell_dist_idx_last']` |
| Cosine top-k graph over terminal-cell distributions | `chronopack.tl._trajectories.cosine_topk_leiden_from_terminal_dist()` | Internal helper; compatibility exported during migration | `adata.obsp['W_terminal_cosine_topk']` |
| Cell-state/fate cluster labels | `cp.tl.cluster_cell_state_trajectories()` | Public task helper | `adata.obs['traj_clusters']` |
| Terminal-cluster fate probabilities | `cp.tl.compute_terminal_fate_probabilities()` | Public task helper | `adata.obsm['terminal_fate_probs_traj']`; per-fate scalar convenience columns in `adata.obs` |
| Fate probability / entropy diagnostic plots | `chronopack.tl._trajectories.fate_umap_2x2_scanpy()` | Plotting helper retained for migration | Plot object; also convenience `adata.obs['fate_*']` columns |
| Mass-retention diagnostics | `chronopack.tl._trajectories.mass_retention_plots_and_umap()` | Diagnostic helper retained for migration | Existing `adata.obs`/`adata.uns` outputs preserved for compatibility |
| Hard cluster proportions over time plot | `chronopack.tl._trajectories.PlotHardClusterMassProportionsOverTime()` | Plotting helper retained for migration | Plot only |
| Dataset paths, tmap directories, cell-type splitting, selected plotting genes, save locations | Active notebooks | Dataset-specific | Notebook-local |

## Public API Introduced

- `cp.tl.get_trajectories()`
- `cp.tl.cluster_cell_state_trajectories()`
- `cp.tl.compute_terminal_fate_probabilities()`
- `cp.tl._trajectories.CellStateTrajectoryResult`

The result object stores WOT sidecar state that is not naturally aligned to
AnnData axes. It also supports legacy dictionary-style access such as
`result["T_pruned"]` during notebook migration.

## AnnData Contract For Cell-State Trajectory Outputs

Cell-aligned outputs:

- `adata.obsm['X_ot_umap']`: WOT graph UMAP embedding.
- `adata.obsm['terminal_cell_dist']`: per-cell distribution over terminal
  cells.
- `adata.obsm['terminal_fate_probs_traj']`: per-cell fate probabilities.
- `adata.obs['traj_clusters']`: per-cell cell-state/fate cluster label.
- `adata.obs['terminal_fate_probs_traj_<fate>']`: compatibility convenience
  columns for each fate probability.

Square cell-by-cell matrices:

- `adata.obsp['ot_graph']`: symmetrized WOT-derived cell graph.
- `adata.obsp['W_terminal_cosine_topk']`: cosine top-k graph over terminal-cell
  distributions.
- Some notebooks additionally store `adata.obsp['t_pruned']` for a square
  pruned transport graph.

Small metadata/provenance:

- `adata.uns['ot_graph']`: Scanpy neighbors metadata for the WOT graph.
- `adata.uns['terminal_cell_dist_idx_last']`: terminal-cell row positions.
- `adata.uns['terminal_fate_probs_traj_fates']` and
  `adata.uns['terminal_fate_probs_traj_names']`: fate labels.
- `adata.uns['chronopack']['cell_state_trajectories']`: ChronoPack provenance.

Non-axis-aligned outputs remain outside AnnData:

- rectangular time-to-time WOT matrices;
- long-form couplings;
- raw `tmaps` tuples;
- large WOT/model state.

These are carried by `CellStateTrajectoryResult` or by explicit sidecar files.

## Behavior Preservation Notes

- WOT parameters, transition top-k values, symmetrization choices, UMAP random
  states, terminal top-k values, Leiden resolutions, and random states remain
  notebook-configured and unchanged.
- The mixed THP1/K562 workflow still runs WOT separately per cell type and
  embeds the resulting square blocks into the full AnnData cell axis.
- Existing compatibility import paths remain temporarily:
  `chronopack.trajectories` and `chronopack._legacy.trajectories`.
- The older function name `terminal_fate_trajectory_clustering()` remains as a
  compatibility shim, but new code should prefer
  `cluster_cell_state_trajectories()` to avoid confusing this workflow with
  gene-expression trajectory clustering.

## Verification

Lightweight package tests:

- `tests/test_cell_state_trajectories.py`
- `tests/test_pre_trade_refactor.py`
- full test suite in the `chronopack` environment: `31 passed, 1 warning`.

Representative full-data check:

- Dataset: `k562-tnfa-puromycin.hg38`.
- Runtime: historical `wot2`.
- Inputs:
  - `verification_runs/upstream_preprocess_gettraj_20260917_132244/k562-tnfa-puromycin.hg38/output/adata_hvg.chk0.h5ad`
  - existing saved WOT maps in
    `verification_runs/upstream_preprocess_gettraj_20260917_132244/k562-tnfa-puromycin.hg38/notebooks/tmaps_sc`
- Comparison script:
  `verification_runs/get_trajectories_package_backed_wot2_20260918_trajectory_migration/compare_puromycin_get_trajectories.py`
- Summary output:
  `verification_runs/get_trajectories_package_backed_wot2_20260918_trajectory_migration/puromycin_get_trajectories_wot2_comparison_summary.json`

Results against the saved historical artifact:

- shape: `4034 x 3000` in both objects.
- retained cell IDs/order: identical.
- retained gene IDs/order: identical.
- `X`: allclose, maximum absolute difference
  `1.7763568394002505e-15`.
- shared layers `counts`, `cpm`, `log1p`, and `corrected`: allclose; raw
  count-derived layers are exact.
- square cell-cell graphs match exactly:
  - `adata.obsp['ot_graph']`
  - `adata.obsp['t_pruned']`
  - `adata.obsp['W_terminal_cosine_topk']`
- terminal-cell distributions match exactly:
  - `adata.obsm['terminal_cell_dist']`
- alias comparison to the older saved branch:
  - current `traj_clusters` matches historical `traj_leiden_sr` exactly.
  - current `terminal_fate_probs_traj` does not equal historical
    `fate_probs_traj_last` (`max_abs_diff=0.2618742883205414`) because that
    older artifact used a separate backward-DP fate calculation rather than
    the active notebook's terminal-cell distribution aggregation.

The saved historical artifact was produced by an older fate-labeling branch
using `traj_leiden_sr` and `fate_probs_traj_last`. The active notebooks now use
the standardized keys `traj_clusters` and `terminal_fate_probs_traj`, so those
obs-level fate probability values are not a direct reference for the active
workflow. This is a reference-artifact naming/version issue, not a discrepancy
in the migrated WOT graph construction, cell-state clustering, or terminal-cell
fate distribution machinery.
