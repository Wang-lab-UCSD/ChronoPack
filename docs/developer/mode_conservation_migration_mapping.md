# Mode Conservation / Mode Switching Migration Mapping

Status: package-backed cross-context mode comparison slice complete.

This pass migrated reusable cross-context dynamic-mode analysis and its direct
visualization primitives. It did not migrate GRN, C2C, or any new scientific
methodology.

## Public API

Analysis:

- `cp.tl.mode_conservation(source, target, ...)`
- `cp.tl.mode_switchers(source, target, ...)`

Plotting:

- `cp.pl.mode_conservation(result, ...)`
- `cp.pl.mode_sankey(result, ...)`
- `cp.pl.mode_chord(result, ...)`

`mode_chord` is available as a wrapper around the existing R/circlize plotting
boundary and requires `circlize` at render time.

## Result Object

`cp.tl.mode_conservation()` returns a `ModeConservationResult` with ordinary
DataFrame payloads:

- `directional_scores`
- `conservation_scores`
- `gene_matches`
- `mode_transitions`
- `switchers`
- `parameters`

The result is intentionally not written into either input AnnData because this
is a cross-context analysis.

## Preserved Conservation Definition

For each source mode:

1. identify source high-confidence dynamic genes;
2. rank source genes by the configured NB-GAM/Wald rank column;
3. retain the current top-k source set;
4. determine the fraction of those genes that are high-confidence in the target
   and retain the same broad mode;
5. repeat in the opposite direction;
6. average the directional recovery fractions to obtain bidirectional
   conservation.

The intentional source/target asymmetry is preserved:

- source: strongest top-k mode members;
- target: broader high-confidence dynamic set, with no rank preservation
  requirement.

This was not replaced by Jaccard, ARI, correlation, or any other metric.

## Current Mapping

| Historical code / behavior | Purpose | Package location | Public? | Input / output |
| --- | --- | --- | --- | --- |
| top-k source selection embedded in Sankey/chord helpers | define source gene universe | `cp.tl.mode_conservation()` internals | public task behavior | standardized AnnData or compatibility DataFrame |
| high-confidence target recovery | directional mode conservation | `cp.tl.mode_conservation()` | public | `ModeConservationResult.directional_scores` |
| bidirectional averaging | final conservation score | `cp.tl.mode_conservation()` | public | `ModeConservationResult.conservation_scores` |
| source-target assignment table | mode switcher classification | `cp.tl.mode_switchers()` / shared implementation | public | `ModeConservationResult.switchers` |
| transition counts | plotting flow table | `cp.tl.mode_switchers()` / shared implementation | public result table | `ModeConservationResult.mode_transitions` |
| `plot_cluster_reassignment_sankey()` analysis portions | top-k union, low-significance/missing bookkeeping | `cp.tl.mode_switchers()` | public task behavior | result object |
| `plot_cluster_reassignment_sankey()` rendering | flow visualization | `cp.pl.mode_sankey()` | public | result transition table |
| `plot_cluster_reassignment_chord_pairwise()` rendering | chord-style transition visualization | `cp.pl.mode_chord()` | public, optional circlize backend | result transition table |
| manuscript-specific panel assembly / titles / filenames | figure assembly | active notebooks | notebook-local | package plot primitives |

## Standardized Inputs

The package reads gene-level outputs from the standardized gene trajectory /
clustering contract:

- gene identifiers from `adata.var["gene_trajectories_gene"]`;
- context/cell-type labels from `adata.var["gene_trajectories_cell_type"]`;
- ranks from `adata.var["gene_trajectories_waldStat_1"]`;
- significance from `adata.var["gene_trajectories_qval_1"]`;
- broad modes from `adata.var["gene_trajectory_clusters_broad_cluster"]`.

The functions also accept a `DataFrame` with historical `pc_df`-style columns
for regression and temporary notebook compatibility.

## Missing / Low-Significance Semantics

Genes that fail the significance threshold or do not belong to a canonical
dynamic mode are represented as `low_significance`.

Missing genes are represented separately as `missing` by default.
`treat_missing_as_low_significance=True` reproduces workflows that analytically
impute missing assignments into the low-significance category.

Plotting options distinguish:

- dropping low-significance flows;
- hiding low-significance-only flows;
- showing the full transition table.

Filtering remains an analysis/result-table choice, not hidden inside plotting
science.

## Active Notebook

Added `k562-tnfa-puromycin.hg38/ModeConservation.ipynb`.

The notebook loads standardized `gene_trajectory_clusters.h5ad`, selects
contexts, calls `cp.tl.mode_conservation()`, writes result tables, and renders
package-backed conservation/Sankey plots. It does not define reusable
mode-matching or Sankey helper functions.

## Regression Verification

Frozen upstream inputs were used; tradeSeq and clustering were not rerun.

Input universe:

- `verification_runs/dimred_package_migration_20260918/output_pkg3/pc_df.new.csv`
- shape: `500 x 35`

This is the current fixed-sidecar universe. The older extended
`~6753`-row universe remains an upstream sidecar/gene-universe issue and was
not reconstructed inside the conservation functions.

Representative comparisons:

- `K562a` to `THP1`
  - top-k: `20`
  - directional rows compared: `6`
  - directional recovery maximum absolute difference versus independent
    stepwise calculation: `0`
  - transition rows: `27`
  - switcher rows: `96`
  - rendered SVGs:
    - `mode_sankey.K562a_to_THP1.svg` (`81,876` bytes)
    - `mode_conservation.K562a_THP1.svg` (`41,084` bytes)

- `K562a` to `K562_puro`
  - top-k: `20`
  - directional rows compared: `6`
  - directional recovery maximum absolute difference versus independent
    stepwise calculation: `0`
  - transition rows: `26`
  - switcher rows: `99`
  - rendered SVG:
    - `mode_sankey.K562a_to_K562_puro.svg` (`81,246` bytes)

## Sidecar Status

The migrated mode-conservation and mode-switching workflow operates from
standardized AnnData or temporary `pc_df`-style compatibility tables.

It does not require:

- `gamList*.rds`;
- `sce*.rds`;
- association CSVs;
- significant-gene CSVs.

Historical `pc_df` CSV/RDS files remain useful as temporary regression and
notebook compatibility artifacts, but they are not the core package
representation.

No native sidecars were deleted.

## Tests

Added `tests/test_mode_conservation.py`.

Synthetic tests cover:

- same-mode recovery;
- mode switching;
- target low-significance;
- source/target missing genes;
- missing-as-low-significance imputation;
- AnnData standardized input extraction;
- top-k directional asymmetry;
- Sankey low-significance drop behavior;
- R-backed conservation SVG render smoke.

Full lightweight suite after this migration: `51 passed, 1 warning`.


## Update 2026-10-01: TESTING.ipynb plotting port

`03_ModeConservation` was removed; its figures are part of
`02_ClusterVisualization`. `cp.pl.mode_sankey()` and `cp.pl.mode_chord()` now
take the gene-trajectory AnnData and draw with the TESTING
`plot_cluster_reassignment_sankey()` / `plot_cluster_reassignment_chord_pairwise()`
(context pairs or a context list, top-k by `waldStat_1`), replacing the
transition-table renderers. New: `cp.pl.mode_proportions()` (pie charts),
`cp.pl.cluster_pair_representative_genes()`, `cp.pl.cluster_violin()` and
`cp.tl.topk_mode_recovery()` (`calc_topk_mode_recovery()`). `cp.tl.mode_conservation()`,
`cp.tl.mode_switchers()` and `cp.pl.mode_conservation()` are kept (per-mode,
bidirectional statistic and switcher tables; not used by the notebooks).
