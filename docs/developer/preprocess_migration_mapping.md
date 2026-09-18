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
| `annotate_qc_flags()` adding `var['mt']`, `var['rb']`, and Scanpy QC metrics | `chronopack.io.annotate_qc_metrics()` | Public |
| `filter_cells_by_qc()` explicit combined cell mask over total counts, mt percent, and detected genes | `chronopack.io.filter_cells_by_qc()` | Public |
| `filter_technical_genes()` excluding explicit genes, MT genes, ENSG-only genes, and ribosomal genes | `chronopack.io.filter_technical_genes()` | Public |
| `filter_low_count_genes()` using total count threshold from a count layer or `X` | `chronopack.io.filter_low_count_genes()` | Public |
| `score_cell_cycle()` Regev gene list parsing plus Scanpy scoring | `chronopack.io.score_cell_cycle()` | Public helper; delegates scoring to Scanpy |
| AnnData contract checks for unique IDs, time metadata, count layer, and matrix dimensions | `chronopack._core._schema` and `chronopack._core._validation` | Internal |
| WOT/tradeSeq handoff expectations such as `obs['time']`, `obs['sample']`, `layers['counts']`, `layers['cpm']`, `layers['log1p']`, and optional `layers['corrected']` | `chronopack._core._schema` and `chronopack.io.validate_preprocessed_adata()` | Public validation wrapper plus internal schema |

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
