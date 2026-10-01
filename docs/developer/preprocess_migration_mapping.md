# PreprocessDatasets Migration Mapping

Date: 2026-09-18

This migration preserves the existing preprocessing behavior and moves only
dataset-independent conventions into ChronoPack.

## Current Code To Package Mapping

| Current function / code block | Proposed ChronoPack location | Public or internal |
| --- | --- | --- |
| `load_time_tagged_anndata()` reading sample DGE tables, prefixing cell IDs, splitting metadata/counts, adding `sample` and time metadata | `chronopack.io.read_time_tagged_dge()` | Public |
| `map_time_to_group()` and `TIME_GROUP_CATEGORIES` | `chronopack._core._schema` constants, re-exported where needed | Internal constants, compatibility re-export |
| `add_time_metadata()` from `FINAL TIME TAG` to `obs['time']` and `obs['time_group']` | `chronopack.io.add_time_metadata()` | Public |
| `annotate_qc_flags()` adding `var['mt']`, `var['rb']`, and Scanpy QC metrics | Scanpy directly in the notebooks (`chronopack.io.annotate_qc_metrics()` removed 2026-09-30) | Not in ChronoPack |
| `filter_cells_by_qc()` explicit combined cell mask over total counts, mt percent, and detected genes | Explicit mask in the notebooks (`chronopack.io.filter_cells_by_qc()` removed 2026-09-30) | Not in ChronoPack |
| `filter_technical_genes()` excluding explicit genes, MT genes, ENSG-only genes, and ribosomal genes | `chronopack.io.filter_technical_genes()` | Public |
| `filter_low_count_genes()` using total count threshold from a count layer or `X` | `chronopack.io.filter_low_count_genes()` | Public |
| `score_cell_cycle()` Regev gene list parsing plus Scanpy scoring | `chronopack.io.score_cell_cycle()` | Public helper; delegates scoring to Scanpy |
| AnnData contract checks for unique IDs, time metadata, count layer, and matrix dimensions | `chronopack._core._schema` and `chronopack._core._validation` | Internal |
| WOT/tradeSeq handoff expectations such as `obs['time']`, `obs['sample']`, `layers['counts']`, `layers['cpm']`, `layers['log1p']`, and optional `layers['corrected']` | `chronopack._core._schema` (`chronopack.io.validate_preprocessed_adata()` removed 2026-09-30) | Internal schema |

## Notebook-Specific Code To Keep Outside The Package

- Input paths and sample names.
- QC threshold choices:
  - K562+TNFa: `min_counts=500`, `max_counts=15000`, `max_pct_mt=5`,
    `min_genes=500`.
  - Puromycin: `min_counts=200`, `max_counts=10000`, `max_pct_mt=7.5`,
    `min_genes=300`.
  - THP1+K562+LPS: first `min_counts=500`, `max_counts=10000`,
    `max_pct_mt=7.5`, then removing sample `N704`, then `min_genes=500`.
- Scrublet execution and doublet filtering calls.
- Standard Scanpy processing: normalization, log transform, HVG selection,
  regression, scaling, PCA, neighbors, UMAP, Leiden, ComBat.
- Manual cluster removals used as QC decisions.
- THP1/K562 marker lists and cell-type assignment.
- Plotting/inspection cells.
- Output file paths.

## Minimum AnnData Contract For Current Downstream Workflows

Required after preprocessing:

- `adata.obs_names` are unique cell identifiers.
- `adata.var_names` are unique gene identifiers/symbols used by downstream
  notebooks.
- `adata.obs['time']` is numeric and has one value per cell.
- `adata.obs['time_group']` is categorical with the four existing coarse time
  groups when generated from time tags.
- `adata.obs['sample']` exists for these experiments.
- `adata.layers['counts']` stores raw retained counts before normalization.
- `adata.layers['cpm']` stores normalized counts before log transform when
  produced by the active notebooks.
- `adata.layers['log1p']` stores log1p-normalized expression when produced by
  the active notebooks.
- `adata.var['mt']` and `adata.var['rb']` exist after QC annotation.
- Standard Scanpy QC columns from `calculate_qc_metrics()` exist after QC
  annotation.
- Matrix, layer, `obs`, and `var` dimensions are aligned.

Optional fields used by some downstream paths:

- `adata.layers['corrected']` for WOT input after ComBat/residualization.
- `adata.obs['cell_type']` for the mixed THP1/K562 experiment.
- `adata.obsm`/`adata.obsp` embeddings and graphs copied from the HVG object
  for downstream inspection and WOT notebooks.

## Deferred Methodology Notes

- The notebooks use manual QC cluster exclusions. These are retained as
  notebook-level choices rather than converted into a ChronoPack method.
- Scanpy operations remain direct Scanpy calls unless ChronoPack adds workflow
  semantics.
- No thresholds or biological filters should change in this migration.

## Full-Data Regression

Representative dataset:

- `k562-tnfa-puromycin.hg38`

Reference:

- old preprocessing implementation under the historical `wot2` runtime.

Validation run:

- package-backed preprocessing implementation under the same historical `wot2`
  runtime.
- verification directory:
  `verification_runs/preprocess_package_backed_wot2_20260918_134429`

Result:

- `adata_full.chk0.h5ad`: old and package-backed outputs are both
  `4034 x 18833`.
- `adata_hvg.chk0.h5ad`: old and package-backed outputs are both
  `4034 x 3000`.
- retained cell IDs and ordering match exactly.
- retained gene IDs and ordering match exactly.
- `X` matrices match exactly.
- shared layers `counts`, `cpm`, and `log1p` match exactly.
- relevant obs/var metadata match exactly except for small `doublet_score`
  numerical differences.
- `predicted_doublet` assignments match exactly.
- final retained-cell set is unchanged.

The maximum observed `doublet_score` absolute deviation was
`0.015541014255667104`. This is not treated as a refactor blocker because the
doublet calls and final preprocessing outputs are unchanged.

Conclusion:

- PreprocessDatasets package migration is scientifically equivalent under the
  historical `wot2` runtime.

## Unified Runtime Compatibility Note

The unified ChronoPack runtime produces materially different Leiden partitions
in the puromycin preprocessing notebook. This is documented as runtime/version
sensitivity, not as a package-refactor regression.

Observed versions:

- historical `wot2`: Scanpy `1.11.1`, AnnData `0.11.4`, NumPy `2.2.6`,
  SciPy `1.15.3`, scikit-learn `1.7.1`, scikit-image `0.25.2`.
- unified ChronoPack runtime: Scanpy `1.11.5`, AnnData `0.12.19`,
  NumPy `2.4.6`, SciPy `1.17.1`, scikit-learn `1.9.1`,
  scikit-image `0.26.0`.

First manual QC exclusion:

- `wot2`: 6 Leiden clusters; numeric cluster `4` contains 145 cells.
- unified runtime: 5 Leiden clusters; numeric cluster `4` contains 26 cells.
- ARI: `0.5046855055033237`.
- NMI: `0.5158418451310945`.
- unified numeric cluster `4` exactly matches `wot2` cluster `5`, indicating
  that numeric cluster IDs are not stable across runtimes.

Second manual QC exclusion:

- `wot2`: 4 Leiden clusters; numeric cluster `3` contains 66 cells.
- unified runtime: 5 Leiden clusters; numeric cluster `3` contains 948 cells.
- ARI: `0.4458315620650641`.
- NMI: `0.4938063221387045`.
- the intended `wot2` excluded cluster `3` maps best to unified cluster `4`
  with Jaccard `0.875`.

This issue is left open as a separate compatibility/methodology decision. This
migration does not repin the unified environment to `wot2`, redo manual QC, or
establish a new preprocessing baseline.

Full preprocessing execution in the unified environment also surfaced direct
runtime dependencies required by the active notebooks:

- `scikit-image` for Scrublet.
- `python-igraph` for Leiden clustering.
- `leidenalg` for Leiden clustering.
