# ChronoPack Progress Summary

Date: 2026-09-16

## Directory Familiarization

- Inspected the three experiment notebook directories:
  - `k562-tnfa.hg38.rep/`
  - `k562-tnfa-puromycin.hg38/`
  - `thp1-k562-lps.hg38.rep/`
- Wrote a durable repo/task understanding note:
  - `REPO_UNDERSTANDING.md`
- Key understanding:
  - preprocessing and WOT/fate code are mostly Python
  - tradeSeq fitting, dynamic mode clustering, volcano/local-FDR-style analysis are mostly R
  - ChronoPack will likely need both Python and R components, or at least a clear cross-language boundary

## Notebook Cleanup

Created `ignore/` subdirectories and moved lower-priority notebooks there:

- `k562-tnfa.hg38.rep/ignore/GetGeneLengths.ipynb`
- `k562-tnfa.hg38.rep/ignore/VolcanoPlot.v2.ipynb`
- `k562-tnfa-puromycin.hg38/ignore/VolcanoPlot.ipynb`
- `thp1-k562-lps.hg38.rep/ignore/LocalFDR.ipynb`

## Dimensionality Reduction Notebook Comparison

Compared the three notebooks in `k562-tnfa-puromycin.hg38/`:

- `DimensionalityReductionClustering.Old.Old.ipynb`
- `DimensionalityReductionClustering.Old.ipynb`
- `DimensionalityReductionClustering.ipynb`

Current interpretation:

- `Old.Old` is the earliest exploratory version.
- `Old` is the first reorganized notebook that generated many reusable intermediate `.rds` files.
- `DimensionalityReductionClustering.ipynb` is the newest/current analysis notebook.
- The current notebook was not fully self-contained:
  - it read old checkpoint/list objects such as `asso_df_ls.rds`, `gamList_ls.rds`, `sce_ls.rds`, `counts_sub_ls.rds`, `PT_ls.rds`, `CW_ls.rds`, `sf_ls.rds`
  - it also relied on some fresh-kernel state, including clustering helpers and `X_n`/`Xd_n`

## New Self-Contained Notebook

Created:

- `k562-tnfa-puromycin.hg38/DimensionalityReduction&Clustering.self_contained.ipynb`

Goal:

- regenerate the objects needed to produce `pc_df.new.csv` and `pc_df.new.rds`
- avoid depending on outputs produced by `DimensionalityReductionClustering.Old.ipynb` or `DimensionalityReductionClustering.Old.Old.ipynb`

Main changes:

- Reads `OriginalDirectoryPath.txt` to locate the original `k562-tnfa-puromycin.hg38/notebooks` directory.
- Resolves sibling original analysis directories for:
  - `k562-tnfa.hg38.rep`
  - `thp1-k562-lps.hg38.rep`
  - `k562-tnfa-puromycin.hg38`
- Loads primary tradeSeq outputs directly from upstream tradeSeq directories.
- Regenerates extension fits for `K562a`, `K562b`, and `THP1` based on `.fixed` gene sets.
- Uses base `K562_puro` tradeSeq output directly; this matched the old checkpoint semantics for `K562_puro`.
- Reconstructs:
  - trajectory matrices
  - `pc_df`
  - filtering
  - PCA
  - finalized clustering
  - `pc_df.new.csv`
  - `pc_df.new.rds`
- Writes regenerated outputs under:
  - `k562-tnfa-puromycin.hg38/self_contained_outputs/global_analysis_extended/`

Known caveat:

- The notebook captures fixed-subset extension logic and exclusion logic.
- It does not yet explicitly add “genes dynamic in other contexts” unless those genes are already represented in the `.fixed` subset files.

## Exclusion Logic

The new notebook includes:

- low-expression exclusion via existing `exclude_low-count-genes*.csv`
- time-tag-bias exclusion via existing `exclude_time-tag-bias_genes*.csv`

These are read from the original notebook directories, not generated anew.

## Fresh-Kernel Fixes

The new notebook fixes issues found in the current notebook:

- Makes `plot_trajectory_clustering_metrics()` active rather than only present in commented code.
- Makes `calculate_trajectory_clustering_metrics()` active rather than only present in commented code.
- Defines `X_n` and `Xd_n` explicitly before final clustering.
- Removes/fixes the stale `cl_repressed` reference in the repressed-cluster renaming cell.
- Adds a file helper to handle both normal and doubled-dot fixed filenames such as:
  - `cell_times.csv.gz`
  - `cell_times..csv.gz`
  - `cell_times.full.csv.gz`

## Tests And Verification

Added structural tests:

- `tests/test_dimensionality_reduction_self_contained.py`

The tests check that:

- the new notebook exists
- it writes `pc_df.new.csv` and `pc_df.new.rds`
- it does not actively read old checkpoint objects from the old dimensionality-reduction notebooks
- all four datasets are configured as primary tradeSeq inputs
- fresh-kernel clustering helpers are active
- optional reference validation is gated by `CHRONOPACK_VALIDATE_REFERENCE=1`

`pytest` was not installed in the checked environments, so the tests were executed directly as Python assertions and passed.

R parse validation was also run in the `tradeseq_ce` environment and passed.

## Slurm Verification

Added:

- `jobs/run_dimensionality_reduction_self_contained.slurm`

This job:

- activates `tradeseq_ce`
- runs the new notebook through `jupyter nbconvert`
- requests 96GB memory
- writes outputs to `self_contained_outputs/global_analysis_extended`
- enables `CHRONOPACK_VALIDATE_REFERENCE=1`

Slurm history:

- Job `924799` failed before output generation because fixed directories used `cell_times..csv.gz`, while the first notebook version looked only for `cell_times.csv.gz`.
- The notebook was patched to handle the fixed-directory filename convention.
- Job `924802` was submitted and was running on `bigcompute` at last check.

Important status:

- Full end-to-end verification is not complete yet.
- We have not yet confirmed that regenerated `pc_df.new.csv` and `pc_df.new.rds` match the original reference object in row count, values, and cluster assignments.

## Time-Tag-Bias Filter Verification

Tracing backward from the original
`k562-tnfa-puromycin.hg38/intmd/global_analysis/extended/pc_df.new.rds`
and `pc_df.new.csv` shows that the final dimensionality-reduction object
does use the time-tag-bias and low-count exclusion lists for all four
experiment/cell-type contexts:

- `K562a`: `k562-tnfa.hg38.rep/notebooks/exclude_time-tag-bias_genes.csv`
  and `exclude_low-count-genes.k562.csv`
- `K562b`: `thp1-k562-lps.hg38.rep/notebooks/exclude_time-tag-bias_genes.k562.csv`
  and `exclude_low-count-genes.k562.csv`
- `THP1`: `thp1-k562-lps.hg38.rep/notebooks/exclude_time-tag-bias_genes.thp1.csv`
  and `exclude_low-count-genes.thp1.csv`
- `K562_puro`: `k562-tnfa-puromycin.hg38/notebooks/exclude_time-tag-bias_genes.csv`
  and `exclude_low-count-genes.k562.csv`

The time-tag-bias lists are generated by comparing each experiment's temporal
pseudobulk logCPM trajectories to the same K562 negative-control aggregate:

- `thp1-k562-lps.hg38.rep/intmd/k562_neg-ctrl_mtx-agg.pkl`

The K562 TNF and puromycin notebooks read that object via relative paths into
the LPS experiment directory. The LPS notebook reads the local copy directly.
Some save cells in the LPS notebook are now commented out, but the output CSVs
exist in the original notebook directory and are used by the final `pc_df.new`
trace.

Object-level checks on the original `pc_df.new` reference:

- RDS and CSV dimensions both equal `6724 x 38`.
- RDS and CSV row names and column names match exactly.
- Numeric columns round-trip through CSV with maximum absolute difference
  approximately `5.1e-13`.
- Strict whole-data-frame equality is false only because CSV loses R factor
  classes for cluster columns.
- No final `pc_df.new.rds` rows remain whose `gene` is present in the
  corresponding time-tag-bias or low-count exclusion list for its `cell_type`.

## Preprocessing / Trajectory / tradeSeq-Prep Refactor

Added a first local ChronoPack-style helper package:

- `chronopack/preprocessing.py`
- `chronopack/trajectories.py`
- `chronopack/tradeseq.py`

Updated the pre-tradeSeq notebooks to import shared helpers from `chronopack`
while preserving experiment-specific paths, manual QC decisions, plotting, and
threshold choices:

- `k562-tnfa.hg38.rep/PreprocessDatasets.ipynb`
- `k562-tnfa.hg38.rep/GetTrajectories.ipynb`
- `k562-tnfa.hg38.rep/PrepareForTradeSeq.ipynb`
- `thp1-k562-lps.hg38.rep/PreprocessDatasets.ipynb`
- `thp1-k562-lps.hg38.rep/GetTrajectories.ipynb`
- `thp1-k562-lps.hg38.rep/PrepareForTradeSeq.ipynb`
- `k562-tnfa-puromycin.hg38/Preprocess_GetTrajectories_MapFates.ipynb`
- `k562-tnfa-puromycin.hg38/PrepareForTradeSeq.ipynb`

Shared helpers now cover:

- standard time-group labeling and QC annotation
- technical gene filtering
- WOT transport-map loading and transition-matrix construction
- sparse adjacency assembly
- terminal-fate and trajectory-clustering support functions
- temporal pseudobulk CPM/logCPM generation
- negative-control logCPM loading
- time-tag-bias correlation scoring
- tradeSeq input CSV export

Added structural tests:

- `tests/test_pre_trade_refactor.py`

Validation run:

- Notebook JSON structure check passed for all eight refactored notebooks.
- Static test confirmed the refactored notebooks import `chronopack` and no
  longer actively redefine the targeted shared helpers.
- `chronopack/*.py` compile/import smoke test passed in the `wot2` conda
  environment.
- A tiny synthetic `ComputeSmoothCorrelations` run returned the expected output
  columns.

Known boundary of this pass:

- Exploratory/metacell-specific code was intentionally left alone.
- Manual QC cluster removal, marker checks, diagnostic plots, and final
  thresholding cells remain notebook-local.
- The separate Slurm verification job for the self-contained dimensionality
  notebook finished but failed in the R extension-fit branch with
  `offset_extra` missing dimnames. That is independent of this refactor and
  still needs a targeted patch/resubmission if we want full regenerated
  `pc_df.new` verification.

## Additional Pre-Verification Cleanup

Before running end-to-end output verification, the pre-tradeSeq notebooks were
tightened further:

- Added `load_time_tagged_anndata()` to `chronopack/preprocessing.py`.
  The preprocessing notebooks now use this shared loader for time-tagged DGE
  CSVs, sample-prefixing cell IDs, constructing AnnData, and adding `sample`,
  `time`, and `time_group` metadata.
- Added `filter_cells_by_qc()` to combine the repeated cell QC thresholds into
  one explicit mask. Scanpy's native `sc.pp.filter_cells()` allows only one of
  `min_counts`, `min_genes`, `max_counts`, or `max_genes` per call, so this
  helper is the clean combined equivalent.
- Removed metacell/KMeans-related code from the preprocessing notebooks,
  including the puromycin combined preprocessing/trajectory notebook.
- Added `run_wot_transport_workflow()` to `chronopack/trajectories.py`.
  The WOT notebooks now use this helper to run WaddingtonOT, extract WOT
  couplings, construct transition kernels, and attach the OT graph/UMAP where
  appropriate.
- Clarified the two LPS trajectory clustering passes:
  - `traj_leiden_sr`: diagnostic successor-representation clustering,
    resolution 0.5
  - `traj_clusters`: final terminal-fate-distance clustering, resolution 0.6,
    propagated to saved downstream AnnData objects
- Pruned `PrepareForTradeSeq` notebooks by removing commented-out scratch cells,
  diagnostic smoothing plots, BEELINE scratch code, and stale KMeans imports.
- Extracted BEELINE scratch cells to:
  - `k562-tnfa.hg38.rep/ignore/PrepareForTradeSeq_BEELINE_extracted.ipynb`

Validation after cleanup:

- Static syntax checks passed for all edited notebooks.
- No active production notebook contains metacell/KMeans or BEELINE scratch
  references.
- No active production notebook contains the targeted old shared helper
  definitions.
- `tests/test_pre_trade_refactor.py` and
  `tests/test_dimensionality_reduction_self_contained.py` pass when run as
  direct Python assertions.
- `chronopack/*.py` compiles in the `wot2` environment.

Still not done:

- End-to-end notebook execution and output equivalence checks remain pending.

## Dimensionality Notebook Offset Patch

Patched `DimensionalityReduction&Clustering.self_contained.ipynb` so the
tradeSeq extension-fit offset matrix receives `dimnames(counts_extra)` before
subsetting by significant gene names. This addresses the Slurm failure:

- `Error in offset_extra[sig_genes, , drop = FALSE]: no dimnames attribute for array`

Isolated verification was run in the `tradeseq_ce` environment using the real
K562a fixed count input. It confirmed:

- `counts_extra` was `665 x 6857`
- `offset_extra` dimnames matched `counts_extra`
- `offset_extra[sig_genes, , drop = FALSE]` succeeded for a real gene subset

## Interactive Pre-Trade Output Verification

Verification was run inside the current interactive Slurm allocation using the
`wot2` conda environment and sandboxed outputs under
`verification_runs/pretrade_20260916_1119/`.

Confirmed so far:

- `k562-tnfa.hg38.rep/PreprocessDatasets.ipynb` completed in the sandbox.
  The regenerated `adata_full.1e4-variant.chk0.h5ad` and
  `adata_hvg.1e4-variant.chk0.h5ad` matched the reference shapes, obs/var
  names, key obs annotations, and PCA values. UMAP coordinates differed, as
  expected for a stochastic embedding.
- `thp1-k562-lps.hg38.rep/PreprocessDatasets.ipynb` completed after fixing a
  stale diagnostic plotting variable (`genes` -> the computed marker list).
  The regenerated `adata_full.chk0.h5ad` and `adata_hvg.chk0.h5ad` matched
  reference shapes, obs/var names, key obs annotations, and PCA values. UMAP
  coordinates differed.
- `k562-tnfa-puromycin.hg38/Preprocess_GetTrajectories_MapFates.ipynb` now
  executes through WOT, terminal distributions, `traj_leiden_sr`, fate
  probabilities, retention summaries, and final `write_h5ad()` in the sandbox.
  Fixes made during this verification included:
  - keep `adata.obs['time']` numeric in `chronopack.preprocessing.add_time_metadata()`
    so WaddingtonOT's float time comparisons find cells;
  - write WOT maps as `tmaps_sc/tmaps_*.h5ad` by passing an explicit prefix to
    WOT (`tmap_dir / 'tmaps'`);
  - remove display-only cells with missing genes or notebook-only random state;
  - use the original terminal cosine-clustering settings (`topk=500`,
    `resolution=0.5`, `directed=False`, `use_weights=True`);
  - convert pandas Series values in `.uns` to one-column DataFrames before
    writing h5ad with the current AnnData version.

Remaining mismatch:

- The puromycin regenerated outputs do not yet match the reference objects.
  Current sandbox comparison:
  - `adata_full.h5ad`: regenerated `(4380, 26732)` vs reference `(4034, 18833)`
  - `adata_hvg.h5ad`: regenerated `(3873, 3000)` vs reference `(4034, 3000)`
  - HVG common cells: `3525 / 4034` reference cells
  - HVG common genes: `1844 / 3000` reference genes
  - `traj_leiden_sr` is now two clusters but counts differ: regenerated
    `{'0': 2259, '1': 1614}` vs reference `{'0': 2268, '1': 1766}`

Interpretation: the puromycin notebook is now executable but not behaviorally
verified. The remaining work is to align the saved object branch with the
reference: the saved reference appears to contain final trajectory/fate
annotations on both `adata_full` and `adata_hvg`, and it uses an 18,833-gene
filtered full feature space rather than the early raw `adata_full` snapshot.

## Puromycin Split + Downstream Verification

The puromycin combined notebook was split into the same three-stage structure as
the other experiments:

- `k562-tnfa-puromycin.hg38/PreprocessDatasets.ipynb`
- `k562-tnfa-puromycin.hg38/GetTrajectories.ipynb`
- `k562-tnfa-puromycin.hg38/PrepareForTradeSeq.ipynb`

The split was verified in the sandbox at
`verification_runs/pretrade_20260916_1119/k562-tnfa-puromycin.hg38/` using the
`wot2` environment.

Important fixes made during verification:

- Restored puromycin preprocessing thresholds from the original notebook:
  `min_genes=300` and `sc.pp.filter_genes(min_cells=3)`.
- Restored the second Leiden-based manual outlier removal that produces the
  expected 4034-cell puromycin branch.
- Changed final trajectory saving so `adata_full.h5ad` is the full-gene object
  aligned to the final trajectory/fate-annotated HVG object, rather than the
  early raw preprocessing snapshot.
- Changed ChronoPack spline smoothing to match the original notebook: when
  `s=None`, use `nanvar(y) * n * 0.5` before calling `UnivariateSpline`.
  This restored the expected time-tag-bias exclusion set.
- Removed the leftover `%store` diagnostic cell from puromycin
  `PrepareForTradeSeq.ipynb`.

Verification results:

- `PreprocessDatasets.ipynb` wrote:
  - `output/adata_full.chk0.h5ad`: `(4034, 18833)`
  - `output/adata_hvg.chk0.h5ad`: `(4034, 3000)`
- `GetTrajectories.ipynb` completed and wrote final:
  - `output/adata_full.h5ad`: `(4034, 18833)` with exact reference obs/var order
  - `output/adata_hvg.h5ad`: `(4034, 3000)` with exact reference obs/var order
- Recomputed trajectory/fate numerical annotations are close in structure but
  not byte-identical to the old reference because WOT/UMAP/graph steps were
  recomputed. This does not affect the tradeSeq inputs used downstream.
- `PrepareForTradeSeq.ipynb` completed and its downstream files match the
  reference exactly:
  - `subset_genes..csv`: `(13905, 1)`, exact
  - `subset_genes.full.csv`: `(18833, 1)`, exact
  - `cell_times..csv.gz`: `(4034, 2)`, exact
  - `cell_times.full.csv.gz`: `(4034, 2)`, exact
  - `cell_libsize_full.csv`: `(4034, 2)`, exact
  - `counts_genes_x_cells..csv.gz`: `(13905, 4034)`, exact index/columns and
    `maxabs=0`, `nnz_diff=0`
  - `counts_genes_x_cells.full.csv.gz`: `(18833, 4034)`, exact index/columns and
    `maxabs=0`, `nnz_diff=0`
- `output/adata_prefiltered.1e4-variant.fixed.h5ad` matches the existing
  reference prefiltered object on shape `(4034, 13905)`, obs/var names, and
  counts layer values (`maxabs=0`).

## GetTrajectories + PrepareForTradeSeq Verification

Verified in the sandbox at
`verification_runs/pretrade_20260916_1119/` using the `wot2` environment.

K562 + TNFa (`k562-tnfa.hg38.rep`):

- `GetTrajectories.ipynb` completed and regenerated:
  - `output/adata_full.1e4-variant.chk1.h5ad`: `(6857, 20036)`
  - `output/adata_hvg.1e4-variant.chk1.h5ad`: `(6857, 3000)`
- Both regenerated AnnData objects match the reference on shape, obs names,
  var names, sample/time/time_group/phase annotations, predicted_doublet values,
  `traj_clusters` labels, and `leiden` labels.
- Deeper numeric comparison: `counts` and `cpm` layers are exact; normalized
  `X`/`log1p`/`corrected` values differ only at floating-point roundoff
  (`maxabs <= 1.8e-14`); `terminal_cell_dist` is exact; `ot_graph` differs only
  at tiny numerical tolerance (`maxabs <= 1.9e-13`). UMAP coordinates and a few
  score/fate columns differ because those stochastic/downstream annotations were
  recomputed.
- `PrepareForTradeSeq.ipynb` completed and its downstream files match the
  reference exactly:
  - `subset_genes..csv`: `(14274, 1)`, exact
  - `subset_genes.full.csv`: `(20036, 1)`, exact
  - `cell_times..csv.gz`: `(6857, 2)`, exact
  - `cell_times.full.csv.gz`: `(6857, 2)`, exact
  - `cell_libsize_full.csv`: `(6857, 2)`, exact
  - `counts_genes_x_cells..csv.gz`: `(14274, 6857)`, exact index/columns and
    `maxabs=0`, `nnz_diff=0`
  - `counts_genes_x_cells.full.csv.gz`: `(20036, 6857)`, exact index/columns
    and `maxabs=0`, `nnz_diff=0`
- `output/adata_prefiltered.1e4-variant.fixed.h5ad` matches the reference on
  shape `(6857, 14274)`, obs/var names, and counts values (`maxabs=0`).

THP1 + K562 + LPS (`thp1-k562-lps.hg38.rep`):

- `GetTrajectories.ipynb` completed and regenerated:
  - `output/adata_full.chk1.h5ad`: `(3915, 21693)`
  - `output/adata_hvg.chk1.h5ad`: `(3915, 1994)`
- Both regenerated AnnData objects match the reference on shape, obs names,
  var names, sample/time/time_group/phase/cell_type annotations,
  predicted_doublet values, `traj_leiden_sr` labels, and `traj_clusters`
  labels.
- Deeper numeric comparison: `counts` and `cpm` layers are exact; normalized
  `X`/`log1p` values differ only at floating-point roundoff
  (`maxabs <= 8.9e-16`); `terminal_cell_dist` differs by one tiny entry
  (`maxabs <= 5.9e-11`); `ot_graph` differs only at tiny numerical tolerance
  (`maxabs <= 1.1e-13`). UMAP coordinates and a few score/fate columns differ
  because those stochastic/downstream annotations were recomputed.
- `PrepareForTradeSeq.ipynb` completed and its K562 downstream files match the
  reference exactly:
  - `subset_genes..csv`: `(15847, 1)`, exact
  - `subset_genes.full.csv`: `(21693, 1)`, exact
  - `cell_times..csv.gz`: `(3128, 2)`, exact
  - `cell_times.full.csv.gz`: `(3128, 2)`, exact
  - `cell_libsize_full.csv`: `(3128, 2)`, exact
  - `counts_genes_x_cells..csv.gz`: `(15847, 3128)`, exact index/columns and
    `maxabs=0`, `nnz_diff=0`
  - `counts_genes_x_cells.full.csv.gz`: `(21693, 3128)`, exact index/columns
    and `maxabs=0`, `nnz_diff=0`
- Its THP1 downstream files also match the reference exactly:
  - `subset_genes..csv`: `(16126, 1)`, exact
  - `subset_genes.full.csv`: `(21693, 1)`, exact
  - `cell_times..csv.gz`: `(787, 2)`, exact
  - `cell_times.full.csv.gz`: `(787, 2)`, exact
  - `cell_libsize_full.csv`: `(787, 2)`, exact
  - `counts_genes_x_cells..csv.gz`: `(16126, 787)`, exact index/columns and
    `maxabs=0`, `nnz_diff=0`
  - `counts_genes_x_cells.full.csv.gz`: `(21693, 787)`, exact index/columns
    and `maxabs=0`, `nnz_diff=0`
- `output/adata_k562_prefiltered.fixed.h5ad` matches the reference on shape
  `(3128, 15847)`, obs/var names, and counts values (`maxabs=0`).
- `output/adata_thp1_prefiltered.fixed.h5ad` matches the reference on shape
  `(787, 16126)`, obs/var names, and counts values (`maxabs=0`).

## TradeSeqFits Cleanup

Added `chronopack/R/tradeseq_helpers.R` to hold the shared R-side tradeSeq
plumbing used by the fitting notebooks:

- package loading and RNG/palette setup;
- optional RAM reporting;
- robust tradeSeq input loading, including the verified `.fixed` filenames with
  doubled dots such as `counts_genes_x_cells..csv.gz`;
- single-lineage pseudotime/weight/offset construction;
- hard two-lineage pseudotime/weight/offset construction;
- wrapper functions for `fitGAM(..., sce = TRUE)`,
  `fitGAM(..., sce = FALSE)`, and single-/multi-lineage `associationTest`;
- helper functions for selecting significant genes.

Refactored the runnable fitting setup in:

- `k562-tnfa.hg38.rep/TradeSeqFits.ipynb`
- `k562-tnfa-puromycin.hg38/TradeSeqFits.ipynb`
- `thp1-k562-lps.hg38.rep/TradeSeqFits.K562.ipynb`
- `thp1-k562-lps.hg38.rep/TradeSeqFits.THP1.ipynb`
- `thp1-k562-lps.hg38.rep/TradeSeqFits.Both.ipynb`

The single-experiment and per-cell-type LPS notebooks now point at the verified
`.fixed` PrepareForTradeSeq input directories. The LPS combined `Both` notebook
still points at the older `sc_to_tradeSeq` directory because a verified fixed
combined input directory has not yet been generated.

Light cleanup was applied to remove or neutralize stale exploratory cells that
depended on undefined notebook-local variables before the main fitting path.
Plotting and downstream diagnostic sections were left mostly intact.

Verification performed so far:

- All five notebooks are valid JSON.
- `chronopack/R/tradeseq_helpers.R` parses under the `tradeseq_ce` conda
  environment.
- Concatenated R code from each refactored TradeSeqFits notebook parses under
  `tradeseq_ce`.

The expensive tradeSeq fits were not rerun yet, so output equality has not been
verified for these notebooks.

## DimensionalityReductionClustering Self-Contained Follow-Up

Updated both self-contained puromycin dimensionality-reduction notebooks:

- `k562-tnfa-puromycin.hg38/DimensionalityReduction&Clustering.self_contained.ipynb`
- `k562-tnfa-puromycin.hg38/DimensionalityReduction&Clustering.self_contained.executed.ipynb`

Descriptor fixes:

- Kept `get_descriptors()` identical to the current
  `DimensionalityReductionClustering.ipynb` implementation.
- Changed the self-contained descriptor call to `pseudocount = 0`.
- Converted `t_onset` and `t_peak` to the aligned minute scale using the same
  per-context scale factors used when constructing the consolidated trajectory
  matrix (`THP1` uses the existing `0.65` scale factor).

Descriptor verification from saved `traj_mat.rds` and `tt.rds`:

- `t_onset`, `t_peak`, and `peak_log2fc` now match the reference to numerical
  tolerance.
- `max_curv` differs only at numerical tolerance (`maxabs ~= 1.5e-8`).
- `peak_dexp` still differs because the current `get_descriptors()` defines it
  as count-scale `expm1(y_peak) - expm1(baseline)`, while the reference column
  numerically equals `peak_dlogexp`.

Clustering fix:

- PCA remains used for `PC1`-`PC4` visualization.
- Dynamic-mode clustering now builds `X_n` and `Xd_n` from direct trajectory and
  derivative matrices (`X_mat`, `Xd_mat`) instead of from PCA scores
  (`X_pc`, `Xd_pc`).

Fast verification from saved intermediates:

- The broad-cluster match fraction improved from about `0.398` to about
  `0.836` after removing PCA from the clustering representation.
- Broad clusters now use the expected six labels only:
  `impulse-up_1`, `impulse-up_2`, `trend-up`, `transient-down`,
  `trend-down_1`, and `trend-down_2`.
- Some cluster mismatch remains, so the self-contained clustering is closer but
  not yet fully identical to the original reference object.

## Fixed TradeSeq Inputs In DimensionalityReductionClustering

Updated both puromycin self-contained dimensionality-reduction notebooks so
they now load only fixed tradeSeq outputs:

- `k562-tnfa-puromycin.hg38/DimensionalityReduction&Clustering.self_contained.ipynb`
- `k562-tnfa-puromycin.hg38/DimensionalityReduction&Clustering.self_contained.executed.ipynb`

The loader now points directly at the `.fixed` directories for all four
contexts:

- `sc-K562-1e6-to-tradeSeq.fixed`
- `sc-K562-1e4_to_tradeSeq.fixed`
- `sc-K562_to_tradeSeq.fixed`
- `sc-THP1_to_tradeSeq.fixed`

This removes the older in-notebook extension layer that treated non-fixed
tradeSeq outputs as primary and fit extra genes from fixed inputs inside
DimensionalityReductionClustering. That extension/filtering responsibility now
lives upstream in `PrepareForTradeSeq` and `TradeSeqFits`.

Also removed the explicit downstream removal layer for:

- time-tag-bias-dominated genes;
- low-count genes with otherwise significant dynamics.

Rationale: fixed `PrepareForTradeSeq` outputs are expected to have already
applied these filters, and fixed `TradeSeqFits` should only fit the corrected
gene universe. Therefore `DimensionalityReduction&Clustering` should start from
the fixed fitted objects directly and should not subtract those genes again.

Validation performed for this change:

- Both notebooks remain valid JSON.
- Concatenated R code from both notebooks parses under the existing
  `tradeseq_ce` conda environment.

The dimensionality-reduction notebook has not yet been rerun after this fixed
input simplification, so output equality has not been rechecked for this exact
version.

## Upstream PreprocessDatasets -> GetTrajectories Verification

Submitted and completed three Slurm jobs for the upstream half of the pipeline
under:

- `verification_runs/upstream_preprocess_gettraj_20260917_132244/`

Jobs:

- `925583`: `k562-tnfa.hg38.rep/PreprocessDatasets.ipynb` followed by
  `GetTrajectories.ipynb`
- `925584`: `k562-tnfa-puromycin.hg38/PreprocessDatasets.ipynb` followed by
  `GetTrajectories.ipynb`
- `925585`: `thp1-k562-lps.hg38.rep/PreprocessDatasets.ipynb` followed by
  `GetTrajectories.ipynb`

All three jobs completed successfully.

Compared the regenerated final AnnData objects against the original analysis
outputs that the `PrepareForTradeSeq` notebooks expect:

- K562 + TNFa:
  - `output/adata_full.1e4-variant.chk1.h5ad`
  - `output/adata_hvg.1e4-variant.chk1.h5ad`
- Puromycin:
  - `output/adata_full.h5ad`
  - `output/adata_hvg.h5ad`
- THP1 + K562 + LPS:
  - `output/adata_full.chk1.h5ad`
  - `output/adata_hvg.chk1.h5ad`

Core verification results:

- All six final objects match the reference shapes, obs names, and var names.
- K562 + TNFa counts and CPM layers match exactly; log-normalized/corrected
  values differ only at floating-point tolerance; `traj_clusters` and
  `traj_clusters_mass_ret_cell` match exactly.
- LPS counts, CPM, log1p, `X`, `cell_type`, `traj_clusters`, and
  `traj_leiden_sr` match exactly; OT graph values differ only at numerical
  tolerance.
- Puromycin final full/HVG objects match shape, obs names, var names, counts,
  CPM, and log1p layers exactly or at floating-point tolerance. Recomputed WOT
  graph/trajectory clustering differs from the older reference, but
  `PrepareForTradeSeq.ipynb` does not consume those puromycin trajectory
  clustering fields; it consumes the aligned cells/genes/counts/time/expression
  objects, which match as needed for downstream tradeSeq input generation.

This verifies that the current `PreprocessDatasets` and `GetTrajectories`
notebooks regenerate the objects expected by `PrepareForTradeSeq`.

## Final TradeSeqFits Verification

Verified the final remaining TradeSeqFits Slurm job:

- Job/log: `verification_runs/filename_migration_20260917_074750/slurm/lps_k562.925523.log`
- Notebook/run: `thp1-k562-lps.hg38.rep/TradeSeqFits.K562.ipynb`
- Output directory:
  `verification_runs/filename_migration_20260917_074750/thp1-k562-lps.hg38.rep/intmd/sc-K562_to_tradeSeq.fixed/`

The log completed normally and reported:

- `counts_dim`: `15847 x 3128`
- `counts_sub_dim`: `15847 x 3128`
- `n_subset_genes`: `15847`
- `n_tested`: `15822`
- `n_qval_significant`: `910`
- `n_pvalue_significant`: `1953`
- `final_gene_count`: `1953`

Verified generated outputs:

- `asso_df.1trajectory.nknots4.csv`: `15822 x 5`
- `sig_df.1trajectory.nknots4.csv`: `910 x 5`
- `sce.1trajectory.nknots8.rds`: `1953 x 3128`
- `gamList.1trajectory.nknots8.rds`: length `1953`

Canonical and legacy aliases were checked:

- association CSV aliases identical after read;
- significant CSV aliases identical after read;
- SCE RDS aliases identical after read;
- `gamList` canonical and legacy RDS files have identical MD5 checksums.

Additional alignment checks passed:

- final SCE row order exactly matches p-value-significant genes from
  `asso_df`;
- q-value-significant gene set matches `sig_df`;
- all final SCE genes are present in the input count matrix;
- final SCE cell names exactly match count matrix column names.

At this point all four submitted fixed TradeSeqFits jobs have completed and
their generated outputs have been verified:

- K562 + TNFa;
- K562 + TNFa + puromycin;
- LPS THP1;
- LPS K562.

## DimensionalityReduction&Clustering Verification

Verified the latest core dimensionality-reduction/clustering notebook:

- `k562-tnfa-puromycin.hg38/DimensionalityReduction&Clustering.self_contained.ipynb`

During verification, the notebook hit a real regression:

- normalized trajectory matrices were constructed as `X` and `Xd`;
- later PCA/clustering cells still referenced historical names
  `X_mat` and `Xd_mat`;
- execution failed with `object 'X_mat' not found`.

Fixed the notebook by assigning the historical aliases immediately after
p-value-significant row filtering:

- `X_mat <- X`
- `Xd_mat <- Xd`

This preserves the intended clustering behavior: PCA scores are used for
visualization, while clustering distances are computed from normalized
trajectory and derivative matrices rather than PC scores.

Verification run:

- `verification_runs/dimred_latest_20260917_200119/`

Because the original project tree contains only the fixed PrepareForTradeSeq
inputs, the verification run pointed the notebook to the previously verified
fixed TradeSeqFits outputs under:

- `verification_runs/filename_migration_20260917_074750/`

The failed full run successfully generated:

- `traj_mat.rds`: `6753 x 97`
- `tt.rds`: 97 time points

After applying the alias fix, verification resumed from these intermediates and
successfully generated:

- `pc_df.new.csv`
- `pc_df.new.rds`

Output checks:

- RDS and CSV row names match exactly.
- RDS and CSV column names match exactly.
- `pc_df.new.rds`: `6753 x 35`
- row set exactly matches the p-value-significant genes from the four verified
  fixed TradeSeqFits association tables;
- no duplicate rows;
- no missing values in core columns:
  - `gene`
  - `cell_type`
  - `waldStat_1`
  - `meanLogFC`
  - `pvalue_1`
  - `qval_1`
  - `PC1`-`PC4`
  - `broad_cluster`

Cell-type row counts match the verified fixed TradeSeqFits final gene counts:

- `K562a`: `2028`
- `K562_puro`: `1645`
- `K562b`: `1953`
- `THP1`: `1127`

Broad dynamic cluster counts:

- `trend-up`: `1820`
- `impulse-up_1`: `894`
- `impulse-up_2`: `519`
- `trend-down_1`: `1977`
- `trend-down_2`: `1210`
- `transient-down`: `333`

Also regenerated the rendered notebook copy:

- `k562-tnfa-puromycin.hg38/DimensionalityReduction&Clustering.self_contained.executed.ipynb`

The executed notebook was generated from the clean source notebook using the
verified fixed TradeSeqFits outputs and a capped-output executor, so plots and
tables are visible inline without embedding huge intermediate matrices. The
executed notebook contains no error outputs. Its rendered output directory
contains:

- `traj_mat.rds`
- `tt.rds`
- `pc_df.new.csv`
- `pc_df.new.rds`

The rendered `pc_df.new.rds` again has shape `6753 x 35` and the same broad
cluster counts listed above.

Added a separate visualization/review notebook:

- `k562-tnfa-puromycin.hg38/DimensionalityReduction&Clustering.cluster_visualization.ipynb`
- rendered copy:
  `k562-tnfa-puromycin.hg38/DimensionalityReduction&Clustering.cluster_visualization.executed.ipynb`

This notebook loads:

- `pc_df.new.rds`
- `traj_mat.rds`
- `tt.rds`

and visualizes:

- median +/- IQR trajectory curves for the six broad dynamic modes;
- top-Wald example genes from each broad mode;
- median +/- IQR curves for fine clusters within each broad mode;
- representative example genes from fine clusters.

The rendered visualization notebook has no error outputs.

Compared old/reference global-analysis gamList membership against the new fixed
TradeSeqFits gamLists:

- `K562a`: old main `3709`, new fixed `2028`, common `2028`;
- `K562b`: old main `3890`, new fixed `1953`, common `1953`;
- `THP1`: old main `4309`, new fixed `1127`, common `1127`.

Thus the new fixed final gamLists are strict subsets of the old large main
per-context gamLists. However, the old extension gamLists are fully represented
in the new fixed final gamLists:

- `K562a_gamList_ext.rds`: `347/347` genes present in new fixed `K562a`;
- `K562b_gamList_ext.rds`: `168/168` genes present in new fixed `K562b`;
- `THP1_gamList_ext.rds`: `54/54` genes present in new fixed `THP1`.

So the old downstream extension layer is no longer needed for those previously
extended genes, but the new fixed final gamLists do not include every gene from
the older large main gamLists.

Comparison against the old reference
`/new-stg/home/alex/projects/chrono_seq/analyses/k562-tnfa-puromycin.hg38/intmd/global_analysis/extended/pc_df.new.rds`
is not expected to be exact after the fixed-input simplification:

- old reference: `6724 x 38`
- regenerated fixed-output version: `6753 x 35`

The row-count change is explained by the new fixed TradeSeqFits outputs and the
removal of the older downstream filtering/extension layer from the
DimensionalityReduction&Clustering notebook.

## OmniPath TF Prioritization Refactor

Refactored `k562-tnfa.hg38.rep/OmnipathTFPrioritization.ipynb` into a compact,
output-free workflow backed by `chronopack/omnipath_tf.py`.

The notebook now:

- loads OmniPath signaling and transcriptional resources plus local TF and
  ligand-receptor reference files;
- builds the signed split OmniPath graph through ChronoPack helpers;
- ranks TFs for two configured contexts:
  - `TNFa`: seeded as TNF ligand through CellChat ligand-receptor links;
  - `TLR4`: seeded directly from the TLR4 signaling node;
- writes:
  - `intmd/omnipath_tf_ranks/omnipath_ranks.TNFa.csv`
  - `intmd/omnipath_tf_ranks/omnipath_ranks.TLR4.csv`

Added lightweight deterministic tests for the signed diffusion helper in
`tests/test_omnipath_tf.py`. `pytest` is not installed in the available Python
environments, so the same test functions were executed directly under the
existing `omnipath` conda environment.

Verification was run in:

- `verification_runs/omnipath_tf_20260917_160957/k562-tnfa.hg38.rep/`

The regenerated rank CSVs exactly matched the reference files at:

- `/new-stg/home/alex/tmp/grn_inference_v2/reproduction/inputs/external/projects/chrono_seq/analyses/k562-tnfa.hg38.rep/intmd/omnipath_tf_ranks/omnipath_ranks.TNFa.csv`
- `/new-stg/home/alex/tmp/grn_inference_v2/reproduction/inputs/external/projects/chrono_seq/analyses/k562-tnfa.hg38.rep/intmd/omnipath_tf_ranks/omnipath_ranks.TLR4.csv`

For both `TNFa` and `TLR4`:

- shapes matched: `501 x 9`;
- columns matched exactly;
- TF order matched exactly;
- maximum numeric difference was `0.0`;
- `pandas.DataFrame.equals()` returned `True`.

## GRN Inference Notebook Reorganization

Created a top-level `grn_inference/` analysis directory and moved
`OmnipathTFPrioritization.ipynb` there. The notebook now writes local outputs
under:

- `grn_inference/outputs/omnipath_tf_ranks/`

The moved notebook was rerun from its new location. Its generated TF-priority
tables exactly matched the original reference files:

- `omnipath_ranks.TNFa.csv`: shape `501 x 9`, identical columns and values.
- `omnipath_ranks.TLR4.csv`: shape `501 x 9`, identical columns and values.

Added modular GRN helper code:

- `chronopack/R/grn_inference_helpers.R`
  - wraps the LFDR/BUM recalibration workflow from
    `explore_thp1_lfdr_recalibration.R`;
  - produces gene-score, summary, and parameter tables.
- `chronopack/grn_inference.py`
  - wraps ChronoPack TF-gene edge scoring from
    `load_chronopack_beta_lfdr()`;
  - computes `gene_score = 1 - lfdr_bum`;
  - computes `tf_score = priority_act / sum(priority_act)`;
  - computes `EdgeWeight = outer(tf_score, gene_score)`;
  - applies the same normalization, deterministic tie breaking, reference-TF
    filtering, ranking, and percentile calculation as the reproduction script.

Added two new notebooks:

- `grn_inference/GeneScoreLfdrRecalibration.ipynb`
  - reads the protein-coding-gene reference and the K562/THP1 tradeSeq LFDR
    association tables;
  - writes local outputs under `grn_inference/outputs/tables/`;
  - verified exact equality to:
    - `thp1_k562_lfdr_recalibration_scores.csv`
    - `thp1_k562_lfdr_recalibration_summary.csv`
    - `thp1_k562_lfdr_recalibration_params.csv`
- `grn_inference/ChronoPackRegulatoryEdges.ipynb`
  - ingests local gene scores plus local TNFa/TLR4 OmniPath TF scores;
  - writes:
    - `grn_inference/outputs/tables/regulator_prioritization_chronopack_ranked_edges.csv`
  - verified that the generated ChronoPack rows (`ChronoPack_K562` and
    `ChronoPack_THP1`, `4,334,306 x 6`) exactly match the corresponding subset
    of the original combined benchmark output
    `regulator_prioritization_method_specific_ranked_edges.csv` after reading
    the generated CSV back from disk.

## Fate Mapping Modularization

The active trajectory notebooks now use the ChronoPack fate helpers rather than
hand-rolled fate-probability code:

- `chronopack/trajectories.py`
  - `CalculateTerminalFatesAcrossCells()`
  - `terminal_fates_from_terminal_clusters()`
  - `fate_umap_2x2_scanpy()`
- `k562-tnfa.hg38.rep/GetTrajectories.ipynb`
- `thp1-k562-lps.hg38.rep/GetTrajectories.ipynb`
- `k562-tnfa-puromycin.hg38/GetTrajectories.ipynb`

The older combined puromycin notebook,
`k562-tnfa-puromycin.hg38/Preprocess_GetTrajectories_MapFates.ipynb`, was also
updated to call the same helpers. For compatibility with existing exploratory
plots and saved AnnData fields, the puromycin notebooks still populate the
legacy columns:

- `fate_p0_traj_last`
- `fate_p1_traj_last`
- `fate_hard_traj_last_argmax`
- `fate_entropy_traj_last`
- `fate_maxprob_traj_last`

Lightweight verification was run by compiling all code cells in the three live
`GetTrajectories` notebooks and the older combined puromycin notebook; no syntax
errors were found. Stale outputs were cleared from the refactored puromycin fate
cells.

## Cell-Cell Communication Notebook

Added a new top-level `c2c/` analysis directory and copied the original
notebook from:

- `/new-stg/home/alex/tmp/c2c_downsample/CellCommunication.ipynb`

The copied notebook was refactored into:

- `c2c/CellCommunication.ipynb`

It now focuses on the core cell-cell communication inference workflow:

- load secreted CellChat ligand-receptor pairs;
- load THP1/K562 AnnData and tradeSeq grid exports;
- filter to ligand and receptor components expressed in the sender/receiver
  cell types;
- restrict ligands to those with temporal association rows;
- build LR activity curves against a K562 NF-kB/AP-1 response module;
- score LR pairs over the configured early causal window (`0-45` minutes);
- aggregate pair scores into ChronoPack ligand-priority scores;
- write local outputs under `c2c/outputs/`.

Reusable logic was moved into `chronopack/c2c.py`, including LR catalog
filtering, time-series construction, LR activity scoring, ligand aggregation,
and plotting helpers. The module lazy-loads Matplotlib inside plotting
functions so pure scoring helpers can be imported in minimal/broken plotting
environments.

Added focused pure-function tests in `tests/test_c2c.py`. Direct execution of
those tests passed, and `chronopack/c2c.py` plus all code cells in
`c2c/CellCommunication.ipynb` compile without syntax errors.

The full notebook was executed under the `wot2` conda environment and wrote:

- `c2c/CellCommunication.executed.ipynb`
- `c2c/outputs/chronopack_lr_pair_scores.csv`
- `c2c/outputs/chronopack_lr_pair_scores.with_rank.csv`
- `c2c/outputs/chronopack_ligand_scores.csv`
- `c2c/outputs/chronopack_lr_pair_scores.png`
- `c2c/outputs/TNF_TNFRSF1A_NFKBIA_overlay.png`
- `c2c/outputs/selected_gene_expression_over_time.png`

The executed notebook had 17 cells and no error outputs. The exported
`chronopack_lr_pair_scores.csv` was compared against the original source-table
reference:

- `/new-stg/home/alex/tmp/c2c_downsample/downsample_outputs/lr_prioritized_15min.csv`

The refactored output matched the reference shape (`189 x 8`), column order,
ordered ligand-receptor pairs, nonnumeric columns, and numeric values to
floating-point tolerance (`max_abs_diff`: `avg = 0`, `rec_det = 0`,
`logp < 5e-13`). The generated ligand-priority table contained 73 ligands, with
`TNF` ranked first.

## GetTrajectories Cleanup

Added `terminal_fate_trajectory_clustering()` to `chronopack/trajectories.py`.
This helper wraps the repeated terminal-fate workflow:

- `CalculateTerminalFatesAcrossCells()`;
- `cosine_topk_leiden_from_terminal_dist()`;
- weighted Leiden clustering into a configured trajectory-cluster column.

The live `GetTrajectories` notebooks now use this helper for terminal-fate-based
trajectory clustering:

- `k562-tnfa.hg38.rep/GetTrajectories.ipynb`
- `thp1-k562-lps.hg38.rep/GetTrajectories.ipynb`
- `k562-tnfa-puromycin.hg38/GetTrajectories.ipynb`

For `thp1-k562-lps.hg38.rep/GetTrajectories.ipynb`, the old diagnostic
`traj_leiden_sr` clustering cells and plots were removed after confirming that
`PrepareForTradeSeq.ipynb` uses `traj_clusters`, not `traj_leiden_sr`.

For `k562-tnfa-puromycin.hg38/GetTrajectories.ipynb`, the long exploratory
successor-representation and custom retention cells were removed. The notebook
now uses the shared terminal-fate clustering helper plus the modular
`mass_retention_plots_and_umap()` summary.

Static verification was run by compiling `chronopack/trajectories.py` and all
code cells in the three live `GetTrajectories` notebooks. No syntax errors were
found. A source scan confirmed that `traj_leiden_sr` no longer appears in the
THP1+K562 `GetTrajectories` or `PrepareForTradeSeq` notebooks.

Follow-up cleanup made the puromycin trajectory notebook fully analogous to the
K562+TNFa trajectory notebook:

- `k562-tnfa-puromycin.hg38/PreprocessDatasets.ipynb` now owns the
  batch-correction/residualization step and writes a `corrected` layer into
  `../output/adata_hvg.chk0.h5ad`.
- `k562-tnfa-puromycin.hg38/GetTrajectories.ipynb` now starts from
  `adata_full.chk0.h5ad` and `adata_hvg.chk0.h5ad`, sets `X` from
  `layers["corrected"]`, runs WOT, clusters terminal-fate trajectories into
  `traj_clusters`, computes `terminal_fate_probs_traj`, summarizes mass
  retention, and writes `adata_full.h5ad` / `adata_hvg.h5ad`.
- The old puromycin-specific `traj_leiden_sr`, `fate_probs_traj_last`,
  legacy fate columns, and exploratory preprocessing/PCA cells were removed
  from the active `GetTrajectories` notebook.

Outputs were cleared from the active puromycin `PreprocessDatasets`,
`GetTrajectories`, and `PrepareForTradeSeq` notebooks to avoid stale displayed
AnnData schemas. Static code-cell compilation passed for all three notebooks.

## PrepareForTradeSeq Low-Count Filtering

Reviewed the low-count gene handling across the three `PrepareForTradeSeq`
notebooks. The previous notebooks wrote `exclude_low-count-genes*.csv` sidecar
files, but those files were not consumed by `TradeSeqFits`; `TradeSeqFits`
reads `subset_genes.csv` and fits whatever genes were written by
`WriteTradeSeqInputsFromAnnData()`.

Added `filter_low_count_genes()` to `chronopack/preprocessing.py` and updated
the three active `PrepareForTradeSeq` notebooks so genes with total counts
below 10 are removed from the AnnData object before tradeSeq input files are
written:

- `k562-tnfa.hg38.rep/PrepareForTradeSeq.ipynb`
- `k562-tnfa-puromycin.hg38/PrepareForTradeSeq.ipynb`
- `thp1-k562-lps.hg38.rep/PrepareForTradeSeq.ipynb`

The old `exclude_low-count-genes*.csv` writing snippets were removed. Static
code-cell compilation passed for the three notebooks, and
`chronopack/preprocessing.py` compiles.

## Experiment Notebook Organization

Created `reference/` subdirectories in each experiment folder so old/reference
notebooks can be kept out of the active first-level workflow. The puromycin
directory had the relevant first-level legacy notebooks, which were moved to
`k562-tnfa-puromycin.hg38/reference/`:

- `DimensionalityReductionClustering.ipynb`
- `DimensionalityReductionClustering.Old.ipynb`
- `DimensionalityReductionClustering.Old.Old.ipynb`
- `Preprocess_GetTrajectories_MapFates.ipynb`
- executed/minirun DimensionalityReduction verification snapshots

The remaining first-level notebooks in the experiment directories are the
active modularized notebooks.

## DimensionalityReduction&Clustering Modularization

Refactored the active puromycin dimensionality reduction notebook,
`k562-tnfa-puromycin.hg38/DimensionalityReduction&Clustering.self_contained.ipynb`,
so the core reusable analysis logic now lives in
`chronopack/R/dimred_clustering_helpers.R`.

The helper file contains the fixed tradeSeq output loaders, feature selection,
tradeSeq smooth trajectory extraction, trajectory time alignment/rescaling,
trajectory matrix construction, descriptor calculation, normalized trajectory
and derivative matrix preparation, PCA column generation, variance-scaled joint
trajectory/derivative clustering, cluster label assignment, and final `pc_df`
formatting. The notebook now acts mostly as orchestration/configuration plus
visualization helpers.

Static verification:

- `chronopack/R/dimred_clustering_helpers.R` parses successfully in the
  `tradeseq_ce` conda environment.
- The assembled R code cells from the refactored notebook parse successfully in
  the `tradeseq_ce` conda environment.

Follow-up plotting cleanup moved generally useful dimensionality-reduction
visualizations into `chronopack/R/dimred_plotting_helpers.R`, including PC
feature plots, tradeSeq fitted gene curves, aligned trajectory traces,
cluster medoids/top genes, selected-gene overlays, gene highlights,
cluster reassignment plots, representative cluster-pair plots, cluster
confusion heatmaps, and PCA reconstruction plots.

The active notebook now sources both helper files and contains no inline
function definitions. Removed from the active notebook were unused exploratory
helpers for hinge-onset debugging, parametric shape/model fitting, sibling
cluster inspection, an old commented tradeSeq plotting implementation, and an
ad hoc `CXCL3` example plotting cell. Both helper files and the assembled
notebook code parse successfully in the `tradeseq_ce` conda environment.

## Formal Package Skeleton

Started the official package refactor using the staged-skeleton strategy.
Added a formal Python package under `src/chronopack` with public namespaces
`chronopack.io`, `chronopack.tl`, and `chronopack.pl`, plus internal
`_core`, `_utils`, `_backends`, and `_legacy` namespaces. The existing flat
`chronopack/` package now extends its source-checkout import path to the new
`src/chronopack` package so current notebooks can continue importing
ChronoPack during migration.

The first public facades establish the intended API shape while preserving
scientific behavior:

- `chronopack.io` exposes AnnData I/O, validation, preprocessing helpers, and
  tradeSeq input export through lazy compatibility imports.
- `chronopack.tl` exposes staged entry points for cell-state trajectories,
  gene dynamics, gene-expression trajectory characterization, dynamic-mode
  clustering, GRN, and C2C.
- `chronopack.pl` exposes resource-safe access to bundled R plotting helpers.

Bundled R backend resources now live under
`src/chronopack/_backends/r/` and are accessed with `importlib.resources`
rather than working-directory or repository-relative paths. Added tests for
public namespace imports, compatibility imports, and bundled R resource
availability. The dimred notebook tests were updated to reflect that core and
plotting helpers are sourced from ChronoPack rather than defined inline.

Static verification:

- Python import and `py_compile` checks passed for the new package skeleton.
- Direct lightweight test runner passed for package skeleton, dimred notebook
  structure, C2C, OmniPath TF, and pre/trade notebook-structure tests.
- Bundled R backend files parse successfully in the `tradeseq_ce` conda
  environment.

## Unified Environment Setup

Started the reproducible ChronoPack development/runtime environment pass before
migrating additional scientific workflows.

Added:

- `env/environment.yml`: human-maintained mixed Python/R/Bioconductor
  compatibility spec for the full ChronoPack workflow.
- `env/README.md`: creation, editable-install, lockfile, and smoke-test
  instructions.
- `docs/developer/dependency_audit.md`: audited imports, versions observed in
  `wot2`, `omnipath`, `tradeseq_ce`, and `r-fdr`, dependency categories, and
  installer status.
- `r/chronopackR/DESCRIPTION`: R backend dependency metadata placeholder for
  the internal R package layer.

The unified spec solves with Python 3.11, R 4.4.3, tradeSeq 1.20.0,
SingleCellExperiment 1.28.0, slingshot 2.14.0, fdrtool 1.2.18, rpy2 3.6.6,
scanpy/anndata, OmniPath, and WOT. No package-level conflict was observed.

Follow-up verification established that the earlier install failures were
filesystem/transaction issues on NFS rather than dependency conflicts. A tiny
test environment and the full ChronoPack environment both installed
successfully on `/tmp` tmpfs. The verified reference prefix is:

- `/tmp/chronopack-ref-20260918113425`

Phase E smoke tests passed:

- ChronoPack installed in editable mode.
- `import chronopack as cp`, `cp.io`, `cp.tl`, and `cp.pl` work from outside
  the repository root.
- Principal Python dependencies import successfully.
- WOT imports and computes a toy 2 x 2 transport map.
- R 4.4.3 starts from the prefix.
- tradeSeq, SingleCellExperiment, fdrtool, and retained R plotting dependencies
  load successfully.
- rpy2 initializes the prefix R when the prefix `bin` is on `PATH`, and a
  Python -> R -> Python round trip succeeds.
- Bundled R backend resources are accessible through package resources.
- `pytest` passes: 22 tests.

Generated exact snapshots:

- `env/conda-linux-64.lock`
- `env/pip-freeze.txt`

During verification, `cluster` was identified as a direct R plotting
dependency and added to `env/environment.yml` and `r/chronopackR/DESCRIPTION`.

## Packed HPC Environment Deployment

The persistent environment strategy was updated after testing the cluster
filesystem model:

- `/new-stg/home/...` is persistent shared NFS.
- `/tmp` is fast node-local tmpfs.
- `/scratch/$USER/...` is fast node-local storage within compute allocations
  but is not shared back to the login node.
- `$TMPDIR` defaults to NFS and should be overridden for environment work.

Decision:

- Do not maintain a permanent live conda prefix on NFS.
- Do not treat `/scratch` as persistent.
- Keep `env/environment.yml`, `env/conda-linux-64.lock`, and
  `env/pip-freeze.txt` as the canonical environment definition.
- Use a persistent packed conda archive as the HPC deployment artifact.
- Unpack that archive to node-local scratch per Slurm allocation/job.

Created:

- `env/activate_chronopack_hpc.sh`
- `env/smoke_test_chronopack_hpc.sh`
- `jobs/verify_chronopack_packed_env.slurm`

Packed artifact:

- `artifacts/env/chronopack-linux-64-20260918.tar.gz`
- approximately 1.2 GB
- not committed to Git

The packed archive intentionally omits an editable ChronoPack install because
`conda-pack` refuses editable packages. The activation helper unpacks the
dependency runtime, runs `conda-unpack`, activates the prefix, and installs the
current checkout in editable mode.

Packed reconstruction verification:

- Slurm job `925682`
- runtime prefix `/scratch/alex/chronopack-925682/env`
- ChronoPack import from outside the repository root passed.
- `cp.io`, `cp.tl`, and `cp.pl` imports passed.
- WOT toy transport smoke test passed.
- R 4.4.3, tradeSeq, SingleCellExperiment, fdrtool, and retained R plotting
  dependencies loaded.
- rpy2 Python -> R -> Python round trip passed.
- Bundled R resources were found through package-resource lookup.
- `pytest` passed: 22 tests.

Environment milestone status:

- Package skeleton: complete.
- Dependency audit: complete.
- Unified specification: complete.
- Solver compatibility: complete.
- Packed HPC deployment: complete.
- Editable install after unpack: complete.
- Python runtime: complete.
- R runtime: complete.
- Python/R bridge: complete.
- WOT smoke test: complete.
- tradeSeq/fdrtool: complete.
- R plotting dependencies: complete.
- pytest: complete.
- Reproducible env definition: complete.

No additional scientific workflow migration was started during this checkpoint.

## PreprocessDatasets Package Migration

Started the first scientific vertical migration after the packed HPC
environment milestone was accepted.

Added a preprocessing migration mapping:

- `docs/developer/preprocess_migration_mapping.md`

Public preprocessing API now lives under `cp.io`:

- `cp.io.read_time_tagged_dge()`
- `cp.io.add_time_metadata()`
- `cp.io.annotate_qc_metrics()`
- `cp.io.filter_cells_by_qc()`
- `cp.io.filter_technical_genes()`
- `cp.io.filter_low_count_genes()`
- `cp.io.initialize_expression_layers()`
- `cp.io.score_cell_cycle()`
- `cp.io.validate_preprocessed_adata()`

Internal schema/validation additions:

- `src/chronopack/_core/_schema.py`
- expanded `src/chronopack/_core/_validation.py`

Established the current AnnData contract for preprocessing outputs:

- unique cell IDs in `obs_names`
- unique gene IDs/symbols in `var_names`
- numeric `obs["time"]`
- coarse categorical `obs["time_group"]` where time tags are loaded
- sample metadata in `obs["sample"]` for current experiments
- raw counts in `layers["counts"]`
- normalized counts in `layers["cpm"]`
- log-normalized expression in `layers["log1p"]`
- QC flags `var["mt"]` and `var["rb"]`
- Scanpy QC columns in `obs`
- optional `layers["corrected"]` and `obs["cell_type"]` where produced
- aligned matrix/layer/obs/var dimensions

Kept dataset-specific choices in notebooks:

- input paths and sample names
- QC thresholds
- Scrublet calls
- Scanpy normalization/HVG/PCA/neighbors/UMAP/Leiden/ComBat workflow choices
- manual QC cluster exclusions
- THP1/K562 marker lists and cell-type assignment
- output paths

Updated active preprocessing notebooks to use `cp.io` for reusable behavior and
to validate final checkpoint AnnData objects before writing:

- `k562-tnfa.hg38.rep/PreprocessDatasets.ipynb`
- `k562-tnfa-puromycin.hg38/PreprocessDatasets.ipynb`
- `thp1-k562-lps.hg38.rep/PreprocessDatasets.ipynb`

Removed stale trajectory imports from the puromycin preprocessing notebook.

Compatibility shims remain temporarily:

- `chronopack/preprocessing.py`
- `src/chronopack/preprocessing.py`
- `src/chronopack/_legacy/preprocessing.py`

Regression tests:

- Added `tests/test_preprocessing_api.py` with synthetic fixtures for DGE
  loading, QC cell filtering, technical-gene filtering, low-count filtering,
  expression-layer initialization, and explicit validation errors.
- Updated notebook-structure tests to require `cp.io` preprocessing calls.
- Notebook code cells parse successfully.
- Full test suite in the verified packed runtime passed: `27 passed`.

Full-data preprocessing regeneration was not rerun in this pass; this slice
was verified with deterministic fixtures and structural notebook checks.

Full-data regression was subsequently performed for the representative
puromycin preprocessing workflow.

Historical-runtime equivalence:

- Compared old preprocessing + `wot2` against package-backed preprocessing +
  `wot2`.
- `adata_full.chk0.h5ad`: both `4034 x 18833`.
- `adata_hvg.chk0.h5ad`: both `4034 x 3000`.
- retained cell IDs/order: identical.
- retained gene IDs/order: identical.
- `X`: exact match.
- shared layers `counts`, `cpm`, and `log1p`: exact match.
- `predicted_doublet`: identical.
- final retained cells: identical.
- maximum observed `doublet_score` deviation:
  `0.015541014255667104`; acceptable because the filtering outcome is
  unchanged.

Status:

- PreprocessDatasets package migration: COMPLETE.
- Historical-output equivalence under `wot2`: CONFIRMED.
- Unified-runtime clustering compatibility: OPEN / separate issue.

Unified-runtime compatibility finding:

- The unified ChronoPack runtime produces materially different Leiden
  partitions in the manual-QC cells.
- This is documented as version-sensitive graph/Leiden drift, not a package
  refactor regression.
- First manual exclusion: ARI `0.5046855055033237`, NMI
  `0.5158418451310945`; `wot2` cluster `4` has 145 cells, unified numeric
  cluster `4` has 26 cells and instead matches `wot2` cluster `5`.
- Second manual exclusion: ARI `0.4458315620650641`, NMI
  `0.4938063221387045`; `wot2` cluster `3` has 66 cells, unified numeric
  cluster `3` has 948 cells, and the intended `wot2` cluster maps best to
  unified cluster `4` with Jaccard `0.875`.

Relevant versions:

- `wot2`: Scanpy `1.11.1`, AnnData `0.11.4`, NumPy `2.2.6`,
  SciPy `1.15.3`, scikit-learn `1.7.1`, scikit-image `0.25.2`.
- unified runtime: Scanpy `1.11.5`, AnnData `0.12.19`, NumPy `2.4.6`,
  SciPy `1.17.1`, scikit-learn `1.9.1`, scikit-image `0.26.0`.

Full preprocessing in the unified runtime exposed additional direct runtime
dependencies, now recorded in `env/environment.yml`:

- `scikit-image`
- `python-igraph`
- `leidenalg`

The provenance serialization fix is retained so ChronoPack provenance written
under `adata.uns["chronopack"]` is h5ad-safe.

## GetTrajectories Package Migration

Started the second scientific vertical migration after the preprocessing
checkpoint was accepted.

Added migration mapping:

- `docs/developer/get_trajectories_migration_mapping.md`

Terminology boundary:

- This workflow covers WOT-derived cell-state / fate trajectories.
- It does not cover fitted gene-expression trajectories, dynamic modes, or
  gene-trajectory clustering.

Public trajectory/fate API now lives under `cp.tl`:

- `cp.tl.get_trajectories()`
- `cp.tl.cluster_cell_state_trajectories()`
- `cp.tl.compute_terminal_fate_probabilities()`

Internal/result-object additions:

- `chronopack.tl._trajectories.CellStateTrajectoryResult`
- migrated the previous implementation from `_legacy.trajectories` into
  `chronopack.tl._trajectories`
- kept temporary compatibility shims:
  - `chronopack.trajectories`
  - `chronopack._legacy.trajectories`

AnnData output convention:

- `adata.obs['traj_clusters']`: per-cell cell-state/fate cluster labels.
- `adata.obsm['X_ot_umap']`: WOT graph UMAP embedding.
- `adata.obsm['terminal_cell_dist']`: per-cell distribution over terminal
  cells.
- `adata.obsm['terminal_fate_probs_traj']`: per-cell fate probabilities.
- `adata.obsp['ot_graph']`: square WOT-derived cell graph.
- `adata.obsp['W_terminal_cosine_topk']`: square terminal-distribution graph.
- rectangular WOT time-to-time matrices, long couplings, and raw tmap tuples
  stay in `CellStateTrajectoryResult`, not in AnnData.

Updated active GetTrajectories notebooks to use the new `cp.tl` public API:

- `k562-tnfa.hg38.rep/GetTrajectories.ipynb`
- `k562-tnfa-puromycin.hg38/GetTrajectories.ipynb`
- `thp1-k562-lps.hg38.rep/GetTrajectories.ipynb`

Dataset-specific configuration remains notebook-local:

- WOT directories and compute/load choices.
- WOT parameters.
- symmetrization choices.
- transition top-k values.
- UMAP random states.
- fate clustering resolutions/top-k values.
- THP1/K562 split-and-block-embed logic.
- plotting genes and output paths.

Tests added:

- `tests/test_cell_state_trajectories.py`

Verification:

- Targeted trajectory/notebook tests in the `chronopack` environment:
  `15 passed`.
- Full lightweight test suite in the `chronopack` environment:
  `31 passed, 1 warning`.

Representative full-data check:

- Dataset: `k562-tnfa-puromycin.hg38`.
- Runtime: historical `wot2`.
- The check used the historical `chk0` HVG AnnData and the existing saved
  `tmaps_sc` WOT maps, so it did not recompute WOT transport maps.
- Comparison script:
  `verification_runs/get_trajectories_package_backed_wot2_20260918_trajectory_migration/compare_puromycin_get_trajectories.py`
- Summary:
  `verification_runs/get_trajectories_package_backed_wot2_20260918_trajectory_migration/puromycin_get_trajectories_wot2_comparison_summary.json`

Full-data result:

- shape: `4034 x 3000` in package-backed and saved reference objects.
- retained cells/order: identical.
- retained genes/order: identical.
- `X`: allclose, maximum absolute difference
  `1.7763568394002505e-15`.
- shared layers `counts`, `cpm`, `log1p`, and `corrected`: allclose; raw
  count-derived layers are exact.
- square cell-cell graphs match exactly:
  - `ot_graph`
  - `t_pruned`
  - `W_terminal_cosine_topk`
- terminal-cell distribution `terminal_cell_dist` matches exactly.
- current `traj_clusters` matches the older saved `traj_leiden_sr` labels
  exactly.
- current `terminal_fate_probs_traj` does not match the older saved
  `fate_probs_traj_last` matrix (`max_abs_diff=0.2618742883205414`) because
  the older artifact used a separate backward-DP fate calculation rather than
  the active notebook's terminal-cell distribution aggregation.

Reference-artifact caveat:

- The saved historical puromycin artifact uses the older obs/obsm naming branch
  (`traj_leiden_sr`, `fate_probs_traj_last`).
- The active notebooks now use standardized names (`traj_clusters`,
  `terminal_fate_probs_traj`).
- Therefore, the representative full-data check confirms exact equivalence of
  the core WOT graph construction, cell-state clustering, and terminal-cell
  distribution machinery, but the old saved fate-probability matrix is not a
  direct reference for the active workflow.

Status:

- PreprocessDatasets: COMPLETE.
- GetTrajectories / fate migration: COMPLETE.
- Historical active-workflow equivalence: CONFIRMED.
- Lightweight regression tests: passed.
- Representative full-data core trajectory/fate comparison under `wot2`:
  passed with the naming caveat above.
- Fate probability difference versus the older `fate_probs_traj_last` artifact
  remains documented as a historical workflow-version/reference-artifact
  difference, not a package-migration regression.
- Next vertical slice: PrepareForTradeSeq + TradeSeqFits as the
  gene-dynamics fitting workflow.
- Gene-dynamics fitting migration is now in progress below.

Checkpoint note:

- A Git checkpoint was requested after accepting the GetTrajectories migration,
  but this managed session cannot write to `.git`:
  `fatal: Unable to create .../.git/index.lock: Read-only file system`.
- The working tree contains the accepted GetTrajectories changes and the
  current gene-dynamics migration work; a manual checkpoint will be needed from
  a writable Git session.

## Gene-Dynamics Fitting Migration

Started the next vertical slice:

- `PrepareForTradeSeq`
- `TradeSeqFits`

This slice treats them together as the gene-dynamics fitting workflow. The
intended public abstraction is:

- `cp.tl.fit_gene_dynamics(...)`

Added migration mapping:

- `docs/developer/gene_dynamics_migration_mapping.md`

Implemented package pieces:

- `cp.tl.fit_gene_dynamics()` now has a real tradeSeq orchestration path.
- `GeneDynamicsResult` records metadata, association/significant tables,
  selected genes, sidecar paths, and optional backend result state.
- `prepare_tradeseq_inputs()` writes the existing counts/time/subset-gene/full
  library-size CSV inputs from AnnData.
- Scalar association-test outputs can be written into `adata.var` under
  `key_added`, for example `{key_added}_waldStat_1`,
  `{key_added}_pvalue_1`, `{key_added}_qval_1`, and
  `{key_added}_significant`.
- Native SCE/gamList objects remain sidecar RDS outputs rather than being
  stored in AnnData.

Implemented R backend pieces:

- Added `run_tradeseq_single_lineage()` to the bundled R backend:
  `src/chronopack/_backends/r/tradeseq_helpers.R`.
- Synced the source-checkout compatibility helper:
  `chronopack/R/tradeseq_helpers.R`.
- The new R function wraps the existing scientific behavior:
  - read fixed tradeSeq CSV inputs;
  - construct size factors, single-lineage pseudotime, cell weights, and
    offset matrix;
  - screen fit with `nknots=4`;
  - `associationTest(..., global=FALSE, lineages=TRUE, inverse="Chol")`;
  - p-value based gene selection with `pvalue_1 < 0.05`;
  - final SCE and gamList fits with `nknots=8`;
  - write configured CSV/RDS outputs and legacy aliases.

Updated active TradeSeqFits notebooks to thin package-backed R workflow clients:

- `k562-tnfa.hg38.rep/TradeSeqFits.ipynb`
- `k562-tnfa-puromycin.hg38/TradeSeqFits.ipynb`
- `thp1-k562-lps.hg38.rep/TradeSeqFits.K562.ipynb`
- `thp1-k562-lps.hg38.rep/TradeSeqFits.THP1.ipynb`

Runtime note:

- The historical `tradeSeq_CE` environment has R/tradeSeq but does not have
  `rpy2`.
- Therefore historical-runtime equivalence should call the packaged R helper
  directly from R.
- The unified ChronoPack runtime can use the Python-first
  `cp.tl.fit_gene_dynamics()` bridge.

Verification so far:

- R helper source/entrypoint smoke under `tradeSeq_CE`: passed for packaged and
  compatibility helper paths.
- R parse validation for all four active TradeSeqFits notebooks under
  `tradeSeq_CE`: passed.
- Targeted Python tests for gene dynamics and notebook structure:
  `11 passed`.
- Full lightweight pytest suite in the `chronopack` environment:
  `35 passed, 1 warning`.

Still pending for this slice:

- Convert or further simplify active `PrepareForTradeSeq` notebooks.
- Run representative historical-runtime regression under `tradeSeq_CE`
  comparing old notebook implementation versus package-backed R helper.
- Run representative unified-runtime smoke via `cp.tl.fit_gene_dynamics()`.
- Perform full-data/HPC regression for the generated association CSVs and
  SCE/gamList sidecars.
- Do not begin DimensionalityReductionClustering until this boundary is
  reviewed after those checks.
