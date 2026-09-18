# Gene-Dynamics Migration Mapping

Date: 2026-09-18

This migration covers `PrepareForTradeSeq` plus `TradeSeqFits` as one
gene-dynamics fitting workflow. The public task API is Python-first:

```python
cp.tl.fit_gene_dynamics(...)
```

`PrepareForTradeSeq` is treated as internal input preparation for the tradeSeq
backend, not as a separate public API unless a future workflow needs it
independently.

## Workflow Inventory

| Current code / stage | Proposed ChronoPack location | Python vs R | Public vs internal | Output location |
| --- | --- | --- | --- | --- |
| Load `chk1` / final trajectory AnnData objects | Active notebooks | Python | Notebook-local | Notebook variables |
| Split mixed LPS object into K562 and THP1 subsets | Active notebooks, using package helpers where possible | Python | Notebook-local configuration | Notebook variables |
| Negative-control logCPM loading | `chronopack.tradeseq.load_negative_control_logcpm()` initially; later internal under gene-dynamics prep | Python | Compatibility public, gene-dynamics internal | Python DataFrames |
| Temporal pseudobulk CPM/logCPM | internal gene-dynamics preparation helper | Python | Internal | Python DataFrames |
| Time-tag-bias scoring | internal gene-dynamics preparation helper; compatibility wrapper retained temporarily | Python | Internal for workflow, compatibility exported | exclusion set / diagnostic DataFrame |
| Technical-gene filtering | existing `cp.io` preprocessing helpers | Python | Public preprocessing API reused | filtered AnnData |
| Low-count filtering before tradeSeq input writing | existing `cp.io.filter_low_count_genes()` | Python | Public preprocessing API reused | filtered AnnData |
| Write counts, cell times, subset genes, and full library sizes | internal helper called by `cp.tl.fit_gene_dynamics()` or by notebook during transition | Python | Internal; compatibility wrapper retained temporarily | CSV sidecar files |
| Optional TNF soft memberships from `traj_clusters_mass_ret_cell` | active notebook for now | Python | Notebook-local / future internal if used by final workflow | `cluster_memberships.retentions.csv` |
| Read tradeSeq CSV inputs | bundled `tradeseq_helpers.R` | R backend | Internal | R matrices/vectors |
| Size factors from full library sizes | bundled `tradeseq_helpers.R::read_tradeseq_inputs()` | R backend | Internal | R numeric vector `sf = lib_full / mean(lib_full)` |
| Single-lineage pseudotime and cell weights | bundled `tradeseq_helpers.R::prepare_single_lineage()` | R backend | Internal | `PT`, `CW`, offset matrix |
| tradeSeq screen fit | bundled R backend, called by Python bridge | R backend | Internal | native SCE screen object |
| Association test | bundled R backend | R backend | Internal | association DataFrame |
| Dynamic-gene selection | bundled R backend | R backend | Internal | significant genes / `sig_df` |
| Final SCE fit on significant genes | bundled R backend | R backend | Internal | native SCE sidecar RDS |
| Final gamList fit on significant genes | bundled R backend | R backend | Internal | native gamList sidecar RDS |
| Association/significant CSV writing | R backend through explicit output config | R backend | Internal side effect | CSV sidecars |
| Native model serialization | R backend through explicit output config | R backend | Internal side effect | RDS sidecars |
| Fitted-curve extraction | not currently performed in active TradeSeqFits notebooks | Later gene-characterization migration or optional result extraction | R backend / Python | Deferred | Deferred standardized `adata.varm` |
| Gene-mode clustering / descriptors / PCA | later DimensionalityReductionClustering migration | Mixed | Out of scope | Out of scope |

## Active Experiment Configuration

Common tradeSeq settings in the active notebooks:

- input directories are the fixed PrepareForTradeSeq outputs.
- screen fit: `nknots = 4`, `sce = TRUE`.
- association test: `global = FALSE`, `lineages = TRUE`, `inverse = "Chol"`.
- significant-gene screen for the final fit: `pvalue_1 < 0.05`.
- `sig_df` reported with `qval_1 < 0.05`, ordered by decreasing
  `waldStat_1` then `meanLogFC`.
- final fit: `nknots = 8`.
- final native outputs: SCE RDS and gamList RDS.
- parallelization: `BiocParallel::MulticoreParam(workers = 4)`.

Per-experiment differences:

| Experiment / fit | Input directory | `n_points` | Legacy aliases |
| --- | --- | --- | --- |
| K562 + TNFa | `sc-K562-1e4_to_tradeSeq.fixed` | 9 | old `asso_df.1trajectory.csv`, `sig_df.1trajectory.csv`, `gamList.1trajectory.csv`, plus the historical SCE alias |
| K562 + TNFa + puromycin | `sc-K562-1e6-to-tradeSeq.fixed` | 9 | none in the active notebook |
| THP1 + K562 + LPS, K562 subset | `sc-K562_to_tradeSeq.fixed` | 6 | K562 expanded CSV aliases and `sce.knots8.rds` / `gamList.knots8.rds` |
| THP1 + K562 + LPS, THP1 subset | `sc-THP1_to_tradeSeq.fixed` | 6 | THP1 expanded CSV aliases and `sce.knots8.rds` / `gamList.knots8.rds` |

## Current Input/Output Contract

Python-side preparation writes:

- `counts_genes_x_cells.csv.gz`: genes x cells raw counts from
  `adata.layers["counts"]`.
- `cell_times.csv.gz`: `cell_id`, numeric `time`.
- `subset_genes.csv`: filtered genes used by the screen fit.
- `cell_libsize_full.csv`: full-object `total_counts`, written from the
  `tag="full"` object and used for size factors.
- optional compatibility filename variants through the current helper.

R-side fitting reads these files, aligns by cell IDs, constructs:

- `counts_sub`: rows are `subset_genes` intersected with count matrix rows.
- `time_vec`: ordered to `colnames(counts_sub)`.
- `sf`: `lib_full / mean(lib_full)`.
- `PT`: single lineage pseudotime scaled to `[0, 1]`.
- `CW`: one lineage with all weights equal to 1.
- `offset_mat`: gene x cell matrix of `log(sf)`.

R-side fitting writes:

- `asso_df.1trajectory.nknots4.csv`
- `sig_df.1trajectory.nknots4.csv`
- `sce.1trajectory.nknots8.rds`
- `gamList.1trajectory.nknots8.rds`
- configured legacy aliases where required for downstream compatibility.

## Proposed Public API

Initial public API:

```python
result = cp.tl.fit_gene_dynamics(
    adata,
    backend="tradeseq",
    key_added="gene_dynamics",
    layer="counts",
    time_key="time",
    output_dir=...,
    n_points=...,
    screen_nknots=4,
    final_nknots=8,
    p_threshold=0.05,
    workers=4,
    write_sidecars=True,
    return_model=True,
)
```

The function should validate inputs and orchestrate the tradeSeq backend. It
should not expose tradeSeq as the user's primary mental model, even though the
backend is recorded explicitly.

## AnnData Contract

For this migration, native R state remains sidecar state. ChronoPack should not
store SCE or gamList objects in `adata.uns`.

Fields to write when standardized outputs are available:

- gene-level scalar statistics in `adata.var`, namespaced by `key_added`, for
  example:
  - `{key_added}_waldStat_1`
  - `{key_added}_df_1`
  - `{key_added}_pvalue_1`
  - `{key_added}_meanLogFC`
  - `{key_added}_qval_1`
  - `{key_added}_significant`
- provenance and sidecar paths in `adata.uns["chronopack"][key_added]`.

Potential later fields:

- fitted curves in `adata.varm[f"{key_added}_fit"]`.
- matching time grid in provenance.

Fitted-curve extraction is not currently part of the active TradeSeqFits
notebooks; downstream notebooks extract trajectories from the native fitted
objects. Standardized fitted curves should be added only when the extraction
behavior is migrated from the downstream gene-trajectory characterization
workflow.

## Native Backend State

`GeneDynamicsResult` should retain lightweight references to:

- association DataFrame;
- significant-gene DataFrame;
- selected genes;
- input dimensions;
- sidecar paths;
- optional backend metadata.

Native large objects should remain as explicit sidecar files:

- SCE RDS;
- gamList RDS;
- optional screen SCE only if needed for debugging/regression.

## R Backend Refactor Target

The bundled R backend should expose explicit reusable functions, using
packaging-safe resource loading:

- source `src/chronopack/_backends/r/tradeseq_helpers.R` through
  `importlib.resources` / `source_r_resource()`;
- avoid current working directory assumptions in new backend entry points;
- pass input/output directories and configuration explicitly;
- return an explicit list containing association/significant tables, gene
  selections, dimensions, and sidecar paths.

The current helper functions already contain most scientific behavior:

- `read_tradeseq_inputs()`
- `prepare_single_lineage()`
- `make_offset_matrix()`
- `fit_tradeseq_sce()`
- `fit_tradeseq_gamlist()`
- `association_df_single()`
- `significant_single_lineage()`
- alias writing helpers

The migration should wrap these in a higher-level R function rather than
reintroducing notebook global state.

## Verification Plan

Lightweight tests:

- validation errors for missing count layer, missing time, duplicate IDs, and
  matrix misalignment;
- preparation of genes x cells counts and time table from small AnnData;
- size-factor input shape and ordering;
- key-added AnnData writes for scalar association outputs;
- packaged R resource loading;
- Python-to-R smoke function call;
- result object sidecar-path serialization/provenance.

Historical-runtime regression:

- Use `tradeSeq_CE` as the primary equivalence runtime.
- `tradeSeq_CE` currently has R/tradeSeq but does not have `rpy2`, so
  historical-runtime equivalence should call the packaged R helper directly
  from R rather than requiring the Python bridge in that environment.
- Compare historical notebook implementation versus package-backed workflow on
  a representative dataset.
- At minimum compare:
  - identical input genes and cells;
  - identical pseudotime, cell weights, size factors, and offsets;
  - `asso_df` row/column names and numeric values;
  - p-values, q-values, Wald statistics, and significant calls;
  - final selected genes;
  - final SCE/gamList dimensions and gene names;
  - downstream-required sidecar filenames.

Unified-runtime smoke:

- After `tradeSeq_CE` equivalence, run a representative package-backed fit in
  the unified ChronoPack runtime through `cp.tl.fit_gene_dynamics()`.
- Report numerical/version drift separately from refactor equivalence.

## Out Of Scope For This Slice

- DimensionalityReductionClustering.
- trajectory descriptors.
- gene-trajectory PCA/GMM/Ward/K-means dynamic-mode clustering.
- retained R plotting wrappers for gene dynamics.
- mode conservation.
- GRN.
- C2C.
