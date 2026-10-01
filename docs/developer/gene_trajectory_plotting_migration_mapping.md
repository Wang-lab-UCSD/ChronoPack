# Gene-Trajectory Plotting Migration Mapping

Status: package-backed plotting slice complete for reusable gene-dynamics /
dynamic-mode plots. This pass is plotting architecture only; it does not change
tradeSeq fitting, trajectory characterization, PCA, clustering, or mode
assignment.

## Public API

The retained reusable plotting operations are exposed through `chronopack.pl`:

- `cp.pl.gene_trajectory(adata, genes, ...)`
  - Visualizes already evaluated fitted gene-expression trajectories.
  - Consumes standardized fitted curves from `adata.varm["gene_trajectories_fit"]`.
  - Uses gene/cell-type metadata and dynamic-mode assignments from `adata.var`.
  - Returns an `RPlotResult` containing the exact long-form plotting data.

- `cp.pl.cluster_trajectories(adata, ...)`
  - Visualizes median +/- IQR fitted trajectories by dynamic mode or fine
    cluster.
  - Consumes standardized fitted curves and cluster assignments from AnnData.
  - Performs only plotting-data summarization; it does not rerun clustering or
    descriptor calculation.

- `cp.pl.trajectory_features(adata, ...)`
  - Visualizes descriptor/PCA feature space for gene trajectories.
  - Consumes `adata.varm["gene_trajectory_clusters_X_pca"]` and descriptor /
    cluster columns from `adata.var`.

All three functions can return plotting data without starting R when no output
path is requested. When `output` is provided, the wrapper renders through the
packaged R/ggplot2 backend to `pdf`, `svg`, or `png`.

## R Backend

New explicit R backend functions were added to
`src/chronopack/_backends/r/dimred_plotting_helpers.R` and mirrored in
`chronopack/R/dimred_plotting_helpers.R`:

- `chronopack_save_ggplot()`
- `chronopack_mode_palette()`
- `chronopack_plot_gene_trajectory_csv()`
- `chronopack_plot_cluster_trajectories_csv()`
- `chronopack_plot_trajectory_features_csv()`

These functions receive explicit CSV plotting tables and output paths. They do
not depend on notebook globals, current working directory, experiment paths,
native `gamList`, or native `SCE` objects.

## Historical Function Mapping

| Historical helper | What it visualized | New package entry point | Standardized AnnData input | Native sidecar needed? |
| --- | --- | --- | --- | --- |
| `trajectory_long()` | Long-form fitted curves | internal Python plot-data extraction | `varm["gene_trajectories_fit"]`, time grid in `uns` | No |
| `median_iqr_df()` | Median/IQR curve summaries | `cp.pl.cluster_trajectories()` | fitted curves plus cluster columns in `var` | No |
| `plot_median_curves()` | Broad/fine cluster median trajectories | `cp.pl.cluster_trajectories()` | fitted curves plus dynamic-mode/fine-cluster columns | No |
| `plot_example_genes()` | Representative fitted gene panels | `cp.pl.gene_trajectory()` | fitted curves plus gene/cell-type/mode columns | No |
| `feature_plot_pc_simple()` | PCA/feature scatter | `cp.pl.trajectory_features()` for reusable feature/PCA scatter | PCA coordinates in `varm`; descriptors/modes in `var` | No |
| Sankey/chord reassignment helpers | Cross-context reassignment visualization plus analysis | Not migrated in this pass | Requires mode-conservation analysis first | Out of scope |
| Manuscript panel assembly | Figure-specific layout/lettering/output names | Notebook-local | Package plot primitives | No |

## AnnData Fields Consumed

Default keys:

- `adata.varm["gene_trajectories_fit"]`
- `adata.varm["gene_trajectory_clusters_X_pca"]`
- `adata.var["gene_trajectories_gene"]`
- `adata.var["gene_trajectories_cell_type"]`
- `adata.var["gene_trajectories_waldStat_1"]`
- `adata.var["gene_trajectories_peak_log2fc"]`
- `adata.var["gene_trajectory_clusters_broad_cluster"]`
- `adata.var["gene_trajectory_clusters_cluster_fine_induced"]`
- `adata.var["gene_trajectory_clusters_cluster_fine_repressed"]`
- `adata.uns["chronopack"]["gene_trajectories"]["parameters"]["time_grid"]`

The wrappers support explicit `fit_key` and `clustering_key` arguments so
multiple fits/clusterings can coexist.

## Active Notebook

`k562-tnfa-puromycin.hg38/DimensionalityReduction&Clustering.cluster_visualization.ipynb`
is now a thin Python workflow client. It loads
`gene_trajectory_clusters.h5ad`, performs dataset-specific representative-gene
selection, calls `cp.pl`, and writes SVG files.

The notebook no longer defines reusable ggplot helpers. It preserves
notebook-local choices such as broad-mode ordering, fine-cluster iteration, and
representative-gene selection.

## Native Sidecar Status

The migrated reusable plotting functions consume standardized AnnData only.
They do not require native `gamList*.rds` or `sce*.rds`.

Native tradeSeq sidecars are still retained because other unmigrated plotting
or diagnostic code may still depend on them, and because sidecar deprecation
should happen only after the retained plotting and mode-conservation workflows
have both been reviewed.

## Verification

Lightweight tests were added in `tests/test_gene_trajectory_plots.py`:

- standardized AnnData plot-data extraction;
- median/IQR summary agreement with the old notebook logic;
- PCA/descriptor feature extraction;
- missing-key error behavior;
- representative SVG rendering through R/ggplot2.

Visual regression focuses on the underlying plotting data rather than
pixel-perfect image identity.

Representative full-data regression used the fixed DimRed package-migration
outputs in `verification_runs/dimred_package_migration_20260918/output_pkg3`.
The standardized AnnData state was compared against the historical sidecar
tables used by the old visualization notebook:

- selected gene-trajectory long-form values: maximum absolute difference `0`;
- broad-cluster median/IQR summaries: maximum absolute difference `0`;
- rendered SVG smoke outputs:
  - `plot_smoke/cluster_trajectories.svg` (`54,796` bytes);
  - `plot_smoke/trajectory_features.svg` (`254,698` bytes).

The final clustering AnnData must retain
`adata.varm["gene_trajectories_fit"]`; `cp.tl.cluster_gene_trajectories()`
now preserves the characterization AnnData when a characterization result is
provided.


## Update 2026-10-01: TESTING.ipynb plotting port

The gene-clustering figures now follow the manuscript `TESTING.ipynb`
notebook. `cp.pl.cluster_trajectories()` (median-curve overlay) and
`chronopack_plot_cluster_trajectories_csv()` were removed: broad and fine modes
are drawn by `cp.pl.cluster_medoids_and_top_genes()`
(`plot_cluster_medoids_and_top_genes()`, TESTING version, via
`chronopack_plot_cluster_medoids_csv()`). `cp.pl.trajectory_features()` now
draws with `feature_plot_pc_simple()`. Medoids and silhouettes use the joint
clustering features, returned as `GeneTrajectoryClusteringResult.features` and
stored by the clustering notebook in `traj.obsm["X_cluster_features"]`.
