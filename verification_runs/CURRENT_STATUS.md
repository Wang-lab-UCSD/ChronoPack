# Current Status: ChronoPack SCE/AnnData Refactor Verification

Updated: 2026-09-24

## Update (2026-09-24, evening): tradeSeq blocker resolved, staged reruns started

- Root cause of `t.default(do.call(cbind, betaAll)) : argument is not a matrix`:
  `fit_tradeseq_sce()` always forwarded `genes = genes` with default `NULL`,
  overriding tradeSeq's `genes = seq_len(nrow(counts))`. tradeSeq then fit zero
  genes, and `t(NULL)` failed. The cell-level vector offset was not the cause.
- Fix: `fit_tradeseq_sce()` only forwards `genes` when non-NULL (both
  `chronopack/R/tradeseq_helpers.R` and `src/chronopack/_backends/r/tradeseq_helpers.R`).
  Cell-level offsets are unchanged.
- Probe: `verification_runs/tradeseq_genes_null_fix_20260924/` (jobs 931439
  frozen, 931440 maintained). 30 puromycin genes: all fitted, explicit gene
  subset works, association rowData written, vector-offset coefficients match
  the historical matrix-offset fit exactly (max |diff| = 0) in both envs.
- `verification_runs/skimage_shim/` was removed with user approval.
- Full reruns now go stage by stage (no Slurm dependency graph) via
  `verification_runs/submit_pipeline_pair_staged.py` (`setup`, then
  `stage <Name>` for each stage after checking the previous one):
  - Frozen: `verification_runs/full_pipeline_frozen_manuscript_20260924_r6`
  - Maintained: `verification_runs/full_pipeline_maintained_chronopack_20260924_r7`
- Slurm accounting is disabled (`sacct` unavailable). Staged wrappers print
  `CHRONOPACK_STAGE_OK job=<id>` as their last line; its presence in the `.out`
  log is the success signal.


### K562+TNF preprocessing cell-set drift (maintained env) and Option 2

- Stage 1 maintained K562+TNF produced 6696 cells vs the manuscript's 6857
  (397 manuscript cells removed, 236 extra kept). Frozen r6 matched the
  historical output exactly (same cells and order).
- Cause (confirmed, `verification_runs/k562_tnfa_leiden_diag_20260924/`):
  `PreprocessDatasets.ipynb` removes clusters by hard-coded Leiden label
  (`leiden != '2'` at resolution 0.5, then `leiden != '6'` at 1.0). The
  pre-filter population (7131 cells), PCA and kNN edges are identical across
  envs; each env's Leiden reproduces its own labels on either env's graph. The
  partition change comes from the Leiden implementation (leidenalg 0.10.2 ->
  0.12.0 / scanpy 1.11.1 -> 1.11.5), not the graph.
- User decision (Option 2, maintained env only): apply the manuscript's
  cluster removals by barcode. Files live in `verification_runs/maintained_overrides/`
  (notebook variant with only cells 16 and 18 changed, plus
  `exclude_leiden_clusters.manuscript.csv`: 39 `leiden_r0.5_cluster2` + 235
  `leiden_r1.0_cluster6` barcodes derived from the historical manuscript
  output). `submit_pipeline_pair_staged.py setup` applies them to the
  maintained run only. The frozen notebook is unchanged.
- The stale maintained r7 K562+TNF outputs from the label-based run were moved
  (not deleted) to `full_pipeline_maintained_chronopack_20260924_r7/k562-tnfa.hg38.rep/superseded_label_filter_20260924/`;
  GetTrajectories job 931475 was cancelled.
- GetTrajectories checks vs the manuscript reference: frozen K562+TNF,
  puromycin and THP1/K562 match exactly (K562+TNF `traj_clusters` 3793/3064).
  Maintained THP1/K562 matches exactly. Maintained puromycin differs (about 2
  preprocessing cells, ~92% `traj_clusters` agreement), which may be expected
  env drift. `cluster_cells_to_trajectories` also uses Leiden.

### Puromycin: 2 -> 3 trajectory clusters (root cause identified, not changed)

- The original analysis (`~/projects/chrono_seq/analyses/k562-tnfa-puromycin.hg38/output/adata_hvg.h5ad`)
  has 2 trajectory clusters (`traj_leiden_sr` 2268/1766) on the same 4034
  cells, in the same order. The refactor template and frozen r6 have 3 (1783/1426/825).
- Cause: the WOT input layer `corrected` lost its covariate regression. The
  original combined notebook built it as ComBat -> `regress_out(total_counts,
  pct_counts_mt, S_score, G2M_score)`; current `PreprocessDatasets.ipynb` cell
  15 is ComBat only. The original layer's gene-covariate correlations are ~1e-14;
  the current layer's are up to 0.37.
- History: `886ec17` moved ComBat -> regress_out into PreprocessDatasets
  correctly; `a60f5c2` ("Preprocessing verified") dropped `regress_out` to match
  the old standalone notebook's combat-only diagnostic path (that comparison did
  not cover `corrected`). The Sep 18 GetTrajectories check reused the saved
  original tmaps, so the first recomputed run (Sep 19 template) is where 3 clusters appear.
- Not causally confirmed by a rerun. The user chose to leave it as is: DimRed
  excludes `K562_puro`, so the final clustering is unaffected.


### Deferred fixes (user may address later)

- **APPLIED 2026-09-25 (user approved):** both `cluster_gene_trajectories`
  calls now pass `reference=cluster_universe` in the template, frozen r6,
  maintained r7 and repo-root DimRed notebooks (older failed run dirs not
  edited). Previous DimRed outputs were moved to
  `k562-tnfa-puromycin.hg38/self_contained_outputs.reference_sig_20260925/` in
  both runs; DimRed rerun as jobs 932275 (frozen) / 932276 (maintained).
  Verified: both passed; six broad counts match expected; all five cluster
  columns 100% identical to the template for the 5,091 non-puro genes in both
  envs (broad and fine, induced and repressed); K562_puro rows all -1; frozen and
  maintained identical, including centroid distances (max diff 0). Centroid
  distances differ from the template by at most 0.0012 (reference scaling changed).
  Background kept below.
- **Puromycin still influences DimRed clustering via the reference.** In
  `DimensionalityReduction&Clustering.self_contained.ipynb` (template, both run
  copies, repo root), `cluster_universe` excludes `cell_type == "K562_puro"`,
  but both `cp.tl.cluster_gene_trajectories(..., reference=sig, ...)` calls
  pass `sig`, which still includes K562_puro genes. In
  `chronopack/R/dimred_clustering_helpers.R` the reference feeds `z_row` /
  `scale_by_total_column_variance` (`prepare_dimred_clustering_inputs`), which set
  the global shape vs derivative weighting used before `hclust`, so puromycin
  genes shift the features of every clustered gene. Effect seen: maintained
  fine repressed clusters 77% identical to the template (ARI 0.89), broad
  clusters unaffected. Candidate fix: `reference=cluster_universe`. This
  changes the method, so re-verify broad/fine counts after changing it.
  - Why maintained but not frozen: the non-puro repressed inputs are identical
    in all three runs (2,667 genes, max |diff| 1.9e-12). The shape/derivative
    scaling ratio sXd/sX is template 0.954314, frozen 0.954184 (-1.3e-4),
    maintained 0.954463 (+1.5e-4). Re-clustering the maintained genes with the
    frozen scaling reproduces the template partition exactly (ARI 1.0 vs 0.89).
    So a ~1.5e-4 upward shift in derivative weight changes the k=16 Ward cut,
    likely a near-tie merge (not pinpointed); the frozen shift went the other way.
  - The fine repressed partition is therefore fragile at k=16. Changing the
    reference will shift the ratio for the frozen run too, so compare fine
    partitions (not just broad counts) against the manuscript after any fix.
  - What-if test (job 932273, `verification_runs/dimred_reference_exclusion_test_20260925/`):
    rerunning `chronopack_cluster_trajectory_matrix_to_files` with K562_puro
    rows removed from the reference (5,091 rows), using the notebook's k rules
    (broad k=3, fine k = max elbow_score), reproduces the template exactly for
    both envs: broad ARI 1.0 with identical sizes (induced 1478/474/472,
    repressed 1521/897/249), same fine k (17 induced, 16 repressed), fine ARI
    1.0; frozen and maintained assignments identical. Not compared: centroid
    distances and PCA embeddings, which use the scaled features.
- Puromycin trajectory clusters 2 -> 3: `regress_out` missing from the
  `corrected` layer (see the puromycin section above).
- Frozen puromycin tradeSeq screen differs from the template (see the
  TradeSeqFits results); user chose not to pursue it.

### DimRed what-if tests on the maintained run (not applied)

Baseline: alpha = 0.1, `reference=cluster_universe`. Scripts and outputs in
`verification_runs/dimred_alpha_test_20260925/` (job 932280) and
`verification_runs/dimred_partition_scaling_test_20260925/` (job 932282).

| Variant | Induced broad ARI / genes moved | Repressed broad ARI / genes moved | Fine k (ind/rep) | Fine ARI (ind/rep) |
|---|---|---|---|---|
| alpha = 0.25 | 0.934 / 75 (3.1%) | 0.799 / 186 (7.0%) | 15 / 16 (was 17 / 16) | 0.720 / 0.780 |
| per-partition scaling, alpha = 0.1 | 0.936 / 71 (2.9%) | 0.273 / 828 (31.0%) | 16 / 15 | 0.730 / 0.787 |

- Per-partition scaling is equivalent to shared scaling with a smaller alpha
  (induced ~0.058, repressed ~0.029), because the partitions' own slope/shape
  scale ratios (1.28, 1.82) exceed the shared reference's (0.95). So it weights
  slope less, the opposite direction from alpha = 0.25.
- Repressed per-partition: trend-down_2 splits (608 genes join trend-down_1,
  289 stay), and 219 trend-down_1 genes join transient-down.

### Fine k depends on the K grid (TESTING.ipynb comparison)

- `~/projects/chrono_seq/analyses/k562-tnfa-puromycin.hg38/notebooks/TESTING.ipynb`
  reads the template's exported `gene_trajectory_fit[_ddt]_significant_z.csv`,
  scales over all 6,723 rows (like the old `reference=sig`), and clusters with
  `K_grid = 2:50`; the pipeline uses `K_GRID = range(2, 101)`. Broad k=3 is the
  same `cutree` of the same tree, so the broad clusters match.
- Its features equal the pipeline R path to ~1e-12, the Ward trees are identical,
  and the package's metrics with K 2:50 give the same fine k as TESTING (11
  induced, 10 repressed) versus 17/16 on 2:100.
- Reason: in `calculate_trajectory_clustering_metrics`,
  `resolution = (WSS_minK - WSS) / (WSS_minK - WSS_maxK)` and
  `k_norm = (K - min K) / (max K - min K)` are both normalized to the grid
  endpoints, so `elbow_score` (and `k_resolution`) change with the grid's upper
  bound. The induced elbow at 2:50 is nearly flat (k=11 0.59579 vs k=10 0.59578).

### PrepareForTradeSeq / tradeSeq input writer

- All six PrepareForTradeSeq jobs passed (932065-932070). tradeSeq inputs
  (decompressed content) match the manuscript template for all frozen branches
  and for maintained K562+TNF and THP1/K562; maintained puromycin differs
  (expected from its drift).
- Runtime: nearly all is `fit_gene_dynamics(prepare_inputs=True)` writing
  dense float64 counts via `pandas.to_csv` into gzip (frozen K562+TNF: 21 min
  subset incl. h5ad reads, 38 min full-counts file; job 72 min).
- Refactor: `prepare_tradeseq_inputs` (`src/chronopack/tl/_gene_dynamics.py`)
  now casts integer-valued count layers to int64 before writing (non-integer
  layers unchanged). Unit tests added in `tests/test_gene_dynamics.py`
  (6 passed in both envs). Full-data benchmark + R `identical()` check:
  `verification_runs/tradeseq_input_writer_benchmark_20260925/` (jobs 932221
  frozen, 932222 maintained). The r6/r7 runs used the pre-refactor writer.
- TradeSeqFits submitted for both envs: jobs 932213-932220.

### TradeSeqFits results (all 8 jobs passed; genes = NULL fix works at full scale)

- Compared with the template (`verification_runs/compare_tradeseq_outputs.R`):
  LPS K562 and LPS THP1 match exactly in both envs (association stats,
  significant set, final SCE coefficients). K562+TNF matches to ~1e-6 (p) /
  1e-9 (beta) in both envs; same significant set and final genes; the
  degenerate gene KLF2-DT (template Wald 1.5e-7, p=1) now has NA Wald and is dropped.
- Frozen puromycin does NOT reproduce the template's nknots=4 screen: 3,121
  genes differ >1e-3 relative, 9 p<0.05 flips, 3 q<0.05 flips, final SCE
  1,635 vs 1,632 genes. Final nknots=8 coefficients on common genes agree to
  ~1e-9. Maintained puromycin differs more (upstream drift).
- Caveat: all template tradeSeq inputs were rewritten on Sep 20, after the
  template's Sep 19 TradeSeqFits outputs, so the input checksum match does not
  prove the Sep 19 fits used identical inputs. The Sep 19 helper (commit
  79d46ff) used a matrix offset and serial fitting (parallel unset); today's
  uses a cell-level vector offset and parallel = TRUE. Cause for puromycin
  not yet tested; DimRed excludes K562_puro before clustering.
- nknots=4 vs 8 p-values: `verification_runs/compare_knots_pvalues.R`,
  outputs in `verification_runs/knots_pvalue_comparison_20260925/`. The final SCE
  rowData carries nknots=4 screen p-values, and DimRed's `pvalue_1 < 0.05`
  filter uses them, so it keeps every fitted gene.
- DimRed first attempt (932261/932262) failed on launcher issues: `display`
  undefined under plain exec (frozen) and no `chronopack` import path in the
  maintained env. The staged submitter now runs DimRed with
  `run_python_notebook_with_display.py` and `PYTHONPATH=<repo>/src`, matching
  `scripts/rerun_dimred_only.sbatch`. Notebooks unchanged.

### DimRed results (jobs 932265 frozen, 932266 maintained: both passed)

- Broad clusters match the expected counts exactly in both envs:
  impulse-up_1 472, impulse-up_2 474, transient-down 249, trend-down_1 1521,
  trend-down_2 897, trend-up 1478 (K562_puro rows are `-1`).
- Per gene vs template `pc_df.new.csv`: broad clusters 100% identical in both
  envs; frozen fine induced/repressed 100% identical; maintained fine induced
  100%, fine repressed 77% (ARI 0.89, 16 clusters in both).
- Cause of the maintained fine-repressed difference: clustering passes
  `reference=sig`, which still includes K562_puro genes;
  `scale_by_total_column_variance` computes the shape/derivative scaling over
  the reference, so puromycin drift shifts the feature weighting for all
  genes. Not changed; whether puro should be in the reference is a user decision.

### Duplicate final-SCE copies removed

- TradeSeqFits no longer saves duplicate copies of the nknots=8 SCE
  (`sce.1trajectory.nknots8.csv` for K562+TNF, `sce.knots8.rds` for THP1/K562).
  Removed the `sce` legacy alias and `sce_aliases` argument from all
  TradeSeqFits notebooks (template, frozen r6 and maintained r7 copies,
  repo root). Association/significance CSV aliases are unchanged.
- DimRed config (`chronopack_dimred_config`, both R copies) now resolves only
  `sce.1trajectory.nknots8.rds`, so a stale duplicate can't be loaded silently.
- The TradeSeqFits jobs (932213-932220) were extracted before this edit and
  wrote the duplicates; `extracted_r/` records what ran. With user approval the
  six duplicates in the frozen r6 / maintained r7 runs were deleted after
  confirming each was byte-identical to its `sce.1trajectory.nknots8.rds`.
  Template and older run directories still have their copies.
- Tests: 58 passed in both envs (run with the env's `bin/` on PATH; rpy2 needs
  R). The one failure, `test_standard_scanpy_wrappers_are_not_public_io_api`,
  predates this session: commit 8df167b re-exported the cp.io QC wrappers that
  the template notebooks use, and the test still asserts they're absent.


### Code provenance (known gap, deferred by user)

- Both runs import the live working tree (`src/chronopack`, R helpers from
  `chronopack/R/`) at job start; runs do not record a git commit/diff.
- Both runs use the Sep 19 template notebooks
  (`manuscript_full_pipeline_20260919_165539/*/notebooks`), not the Sep 21
  repo-root experiment notebooks (those use plain Scanpy QC, `ot_umap`,
  `cluster_cells_to_trajectories`, explicit SCE-screen TradeSeq outputs).
  DimRed notebooks are identical in both.
- Planned (user, later): replace the repo-root experiment notebooks with the
  notebooks from the maintained run
  (`full_pipeline_maintained_chronopack_20260924_r7`). Caveat: those are
  template-based, so a straight copy would revert the Sep 21 notebook API
  cleanup; merge per notebook when doing this.

The sections below describe the state before this update.

This note is a handoff for continuing the ChronoPack manuscript-pipeline verification work.

## Goal

We are verifying the refactored ChronoPack pipeline in two full start-to-finish runs:

1. **Frozen manuscript environment**
   - Conda env: `/new-stg/home/alex/miniconda3/envs/chronopack-manuscript-20260919`
   - Purpose: reproduce the manuscript pipeline as closely as possible.

2. **Maintained ChronoPack environment**
   - Conda env: `/new-stg/home/alex/miniconda3/envs/chronopack`
   - Purpose: test the maintained package/refactored code path.

The main downstream question is whether the final `DimensionalityReduction&Clustering` stage reproduces the expected gene trajectory features and broad/dynamic cluster assignments, especially for non-puromycin experiments. Puromycin differences may be expected from environment/package drift, but K562+TNF and THP1/K562+LPS should largely match prior outputs.

## Important User Preferences / Constraints

- Be careful with scientific code and data.
- Prefer small, reversible changes.
- Ask before deleting files, modifying shell startup files, installing packages, or changing conda environments.
- Do not introduce ad hoc computational workarounds that change scientific behavior.
- Do **not** replace missing package behavior with approximate algorithmic shims.
- Preserve explicitly requested refactors:
  - fitted tradeSeq `SingleCellExperiment` objects are canonical upstream inputs for gene trajectory extraction;
  - DimensionalityReduction&Clustering uses `AnnData` trajectory objects;
  - the tradeSeq path should use **cell-level offsets** as the intended representation.

## Key Refactors Already Done

### Trajectory Data Model

ChronoPack now has a simpler trajectory-analysis model:

```text
traj: AnnData
obs = gene-context trajectories
var = aligned time points
X   = aligned fitted expression trajectories
```

Important APIs:

- `cp.tl.align_gene_trajectories(sce_dict, ...) -> AnnData`
- `cp.tl.compute_gene_trajectory_descriptors(traj, inplace=True)`
- `cp.tl.cluster_gene_trajectories(traj, k=..., alpha=..., reference=...)`
- `cp.pl.clustering_metrics(...)`

`cluster_gene_trajectories()` computes all candidate `k` solutions once. Its `assignments[k]` entries should be `pd.Series` indexed by trajectory `obs_names`.

### DimensionalityReduction&Clustering Notebook

The manuscript reproduction notebook:

`verification_runs/manuscript_full_pipeline_20260919_165539/k562-tnfa-puromycin.hg38/notebooks/DimensionalityReduction&Clustering.self_contained.ipynb`

was refactored to:

- load fitted SCEs;
- build `traj = cp.tl.align_gene_trajectories(...)`;
- add descriptors in-place;
- define significant/up/down populations explicitly in the notebook;
- exclude `cell_type == "K562_puro"` before clustering;
- use `alpha = 0.1`;
- force broad `k = 3` for induced and repressed partitions after plotting metrics.

Expected non-puromycin broad cluster counts after this DimRed refactor:

```text
Impulse-up 1: 472 genes
Impulse-up 2: 474 genes
Transient-down: 249 genes
Trend-down 1: 1,521 genes
Trend-down 2: 897 genes
Trend-up: 1,478 genes
```

A successful notebook-only run previously matched these counts.

### SCE Association Metadata

The saved/loadable SCEs should carry association fields in `rowData(sce)`:

```text
waldStat_1
meanLogFC
pvalue_1
qval_1
```

There is a migration helper:

`verification_runs/manuscript_full_pipeline_20260919_165539/scripts/migrate_sce_association_rowdata.py`

But the desired endpoint is that upstream SCE-producing workflows write these fields directly.

### Plotting

A SCE-based tradeSeq plotting function was added to the package. It reconstructs fitted smooths from SCE-stored tradeSeq coefficients/design matrices rather than loading `gamList`.

## Validated Older Results

Older K562+TNF saved outputs consistently have exactly 2 cell-trajectory clusters:

```text
traj_clusters:
0: 3793 cells
1: 3064 cells
```

Checked files included:

- `verification_runs/manuscript_full_pipeline_20260919_165539/k562-tnfa.hg38.rep/output/adata_hvg.1e4-variant.chk1.h5ad`
- old upstream/package migration outputs under `verification_runs/old/`

Therefore, if a fresh maintained-environment run gives a different number of K562+TNF trajectory clusters, treat it as a pipeline/environment issue, not as expected biology.

## Environment Notes

The maintained environment initially lacked optional Scanpy runtime dependencies:

- `scikit-image` for `sc.pp.scrublet(..., threshold=None)`
- `igraph` / `leidenalg` for `sc.tl.leiden(...)`

These were installed into:

`/new-stg/home/alex/miniconda3/envs/chronopack`

via that env's `pip`, because conda/mamba solving was killed:

```text
scikit-image 0.26.0
igraph 1.0.0
leidenalg 0.12.0
```

Verification smoke test passed:

- `skimage.filters.threshold_minimum` imports from real scikit-image;
- `scanpy.tl.leiden` uses real `scanpy.tools._leiden`;
- a tiny Leiden test produced expected two clusters.

Do **not** use `verification_runs/skimage_shim`. It was an erroneous temporary workaround and is no longer referenced by the submitter. It should eventually be removed after user approval.

## Submission Script

The script used to submit the paired full pipeline jobs is:

`verification_runs/submit_full_pipeline_pair_20260924.py`

It creates clean verification directories from:

`verification_runs/manuscript_full_pipeline_20260919_165539`

It excludes old generated outputs/intermediates but preserves the required negative-control symlink:

`thp1-k562-lps.hg38.rep/intmd/k562_neg-ctrl_mtx-agg.pkl`

Recent run directories:

- Frozen: `verification_runs/full_pipeline_frozen_manuscript_20260924_r5`
- Maintained: `verification_runs/full_pipeline_maintained_chronopack_20260924_r6`

Their `job_ids.tsv` files contain all submitted jobs and retry jobs.

## Recently Cancelled / Failed Jobs

Fresh full pipeline jobs were submitted:

Frozen `_r5`:

```text
Preprocess: 930978, 930979, 930980
GetTrajectories: 930981, 930982, 930983
PrepareForTradeSeq: 930984, 930985, 930986
TradeSeqFits original: 930987, 930988, 930989, 930990
DimRed original: 930991
```

Maintained `_r6`:

```text
Preprocess: 930992, 930993, 930994
GetTrajectories: 930995, 930996, 930997
PrepareForTradeSeq: 930998, 930999, 931000
TradeSeqFits original: 931001, 931002, 931003, 931004
DimRed original: 931005
```

The original and first retry TradeSeq jobs failed for two separate reasons:

1. R syntax issue:
   - trailing commas inside `outputs = list(...)` / `legacy_outputs = list(...)` config blocks in TradeSeq notebooks.
   - This was patched in the manuscript template and active `_r5`/`_r6` notebooks.

2. Current blocker:
   - `tradeSeq::fitGAM(..., sce = TRUE, offset = <cell-level vector>)` fails with:

```text
Error in t.default(do.call(cbind, betaAll)) : argument is not a matrix
Calls: run_tradeseq_single_lineage ... fitGAM ... .fitGAM -> data.frame -> t -> t.default
```

This happens for all single-lineage TradeSeq jobs in both environments.

An incorrect offset-matrix retry was briefly submitted:

```text
Frozen: 931065, 931066, 931067, 931068, final 931069
Maintained: 931070, 931071, 931072, 931073, final 931074
```

Those jobs were cancelled after the user clarified that the refactor intentionally uses cell-level offsets and we should not revert to matrix offsets.

The helper patch that expanded offsets to a matrix was reverted in both:

- `chronopack/R/tradeseq_helpers.R`
- `src/chronopack/_backends/r/tradeseq_helpers.R`

## Current Blocker

The next agent should investigate how to make the **cell-level offset** tradeSeq path work with SCE output.

Current function:

`fit_tradeseq_sce()` in:

- `chronopack/R/tradeseq_helpers.R`
- `src/chronopack/_backends/r/tradeseq_helpers.R`

currently calls:

```r
tradeSeq::fitGAM(
  counts = counts_sub,
  genes = genes,
  pseudotime = pseudotime,
  cellWeights = cell_weights,
  offset = offset,
  nknots = nknots,
  sce = TRUE,
  parallel = parallel,
  BPPARAM = BiocParallel::MulticoreParam(workers = workers),
  verbose = verbose
)
```

where `offset` is a cell-level vector from:

```r
make_offset_vector <- function(sf) {
  as.numeric(log(sf))
}
```

The old validated notebooks used a matrix offset:

```r
offset_mat = matrix(log(sf), nrow = nrow(counts_sub), ncol = ncol(counts_sub), byrow = TRUE)
```

But the explicit desired refactor is to keep cell-level offsets as the intended representation. Do not globally revert the pipeline to matrix-offset sidecars.

Recommended next step:

1. Run a small isolated R probe on one prepared input directory, e.g.

```text
verification_runs/full_pipeline_frozen_manuscript_20260924_r5/k562-tnfa-puromycin.hg38/intmd/sc-K562-1e6-to-tradeSeq.fixed
```

2. Test `tradeSeq::fitGAM` on a tiny subset of genes with:

- `sce = TRUE`, vector offset
- `sce = FALSE`, vector offset
- `sce = TRUE`, matrix offset
- maybe `parallel = FALSE`
- maybe `genes = rownames(counts_sub_subset)`

3. Determine whether the failure is:

- tradeSeq not supporting vector offsets with `sce = TRUE`;
- interaction with `genes = NULL`;
- parallel worker behavior;
- count matrix type/dimnames issue;
- or a bug in how we pass offsets/design data.

4. Implement the smallest fix that preserves the conceptual API:

```text
cell-level fitting offset is canonical
SCE stores the fitted tradeSeq offset/design info
DimRed can recover reference offsets from SCE
```

If tradeSeq internally requires a matrix at the `fitGAM` call boundary, document that as a private adapter only, not as a return to matrix offsets as pipeline inputs or sidecar artifacts.

## Files Modified Recently

Important modified files include:

- `src/chronopack/io/__init__.py`
  - re-exported notebook-facing preprocessing helpers:
    - `annotate_qc_metrics`
    - `filter_cells_by_qc`
    - `initialize_expression_layers`
    - `validate_preprocessed_adata`

- `src/chronopack/tl/_trajectories.py`
  - changed `get_trajectories()` to preserve `X_umap` while also writing requested `umap_key` such as `X_ot_umap`.
  - This fixed puromycin `GetTrajectories` notebook assumptions.

- `verification_runs/submit_full_pipeline_pair_20260924.py`
  - submits paired frozen/maintained full-pipeline verification runs.
  - no longer injects `skimage_shim` into `PYTHONPATH`.
  - preserves the required negative-control symlink.

- TradeSeq notebooks in:
  - `verification_runs/manuscript_full_pipeline_20260919_165539/*/notebooks/TradeSeqFits*.ipynb`
  - active `_r5` / `_r6` run copies
  - patched R trailing commas in config lists.

## Git / Workspace Notes

The working tree is dirty. Many deleted paths under `verification_runs/` are because the user moved older verification directories into `verification_runs/old/`. Do **not** revert those deletions/moves.

Before making more changes, inspect the current `git status --short` carefully.

## What Not To Do

- Do not reintroduce `gamList` loading into DimensionalityReduction&Clustering.
- Do not use approximate clustering fallbacks instead of installing/using real dependencies.
- Do not globally revert to matrix offsets as the pipeline representation.
- Do not delete old failed verification directories/logs unless the user explicitly asks.
- Do not resubmit more full runs until the tradeSeq cell-level-offset SCE blocker is understood.

