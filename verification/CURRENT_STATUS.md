# Current Status: ChronoPack SCE/AnnData Refactor Verification

Updated: 2026-10-01

## Update (2026-10-01): K562 + TNF-alpha figures (`05_Visualization`)

- New `analyses/k562-tnfa.hg38.rep/05_Visualization.ipynb` (optional stage
  `Visualization` in `run_pipeline.py`) ports the CreateFigures.ipynb figures
  and the TradeSeqFits.ipynb volcano. Inputs: `output/adata_full.1e4-variant.chk1.h5ad`
  and `intmd/sc-K562-1e4_to_tradeSeq.fixed/{asso_df,sce_screen}.1trajectory.nknots4.*`;
  figures (PDF) go to `EXP_DIR/figures/`.
- New API: `cp.tl.rotate_embedding_to_time`, `cp.tl.early_fate_correlations`,
  `cp.tl.signed_max_log2fc_from_start` (R `signed_max_log2fc_start_many`, verbatim,
  in `tradeseq_helpers.R`), `cp.pl.association_volcano` / `cp.pl.fdr_wald_break`
  (R `plot_association_volcano` in `dimred_plotting_helpers.R`; `y_break="fdr"`
  replaces the hard-coded 18.31693), and `cp.pl.cell_trajectories.{time_minutes,
  trajectory_cell_types, fate_probability_panels, two_gene_fate_scatter,
  top_fate_correlations}`. `fate_probabilities` (02 notebooks) is unchanged.
- Check `verification_runs/k562_visualization_check_20261001/` (job 934618, r9
  data via symlinks, editable env): no errors, one image per figure cell; promoted.
  Early-bin correlations identical to the legacy outputs (top/bottom 10 genes and
  r); volcano: 11,865 genes, y_break 17.3383 (RACGAP1), 712 genes q < 0.01, 22 of
  the 30 labels shared with the legacy plot (different gene universe). The rotated
  `ot_umap` is the x-mirror of the legacy layout, hence `flip_x=True`.
- Wheel rebuilt: `install_test_20260925/dist_fix15` (installed in
  chronopack-install-test). Tests: 94 passed.

## Update (2026-09-30): analyses/ layout, frozen env and puromycin removed

Everything below this section predates it: paths such as
`k562-tnfa.hg38.rep/*.ipynb`, `verification_runs/submit_*.py`, the frozen
manuscript env and the puromycin experiment refer to the old layout.

- Notebooks now live in `analyses/` (one copy, in git, keeping the outputs of
  the latest verified run): `k562-tnfa.hg38.rep/0{1..4}_*`,
  `thp1-k562-lps.hg38.rep/0{1..4}_*`, `gene_clustering/0{1..3}_*`, `c2c/`,
  `grn/`, `resources/`, `archive/` (old exploratory notebooks). Data location
  is `CHRONOPACK_DATA_DIR` (default `analyses/data/`); see `analyses/README.md`.
- `verification/run_pipeline.py setup|stage|promote RUN` runs those notebooks
  unmodified from `analyses/` against `verification_runs/RUN/data/`, saving
  executed copies to `RUN/executed/`; `promote` copies them back once checked.
  Compare with `verification/compare_runs.py STAGE NEW_DATA REF_DATA` and
  `verification/compare_tradeseq_outputs.R NEW_DATA REF_DATA` (r8's data root
  is the r8 run folder itself). Jobs request 90G (96G never fits outside
  bigcompute).
- `verification_runs/` is gitignored (run records only).
- Puromycin removed: gene clustering does not depend on it (job 934429,
  `verification_runs/dimred_without_puromycin_20260930/`: all broad/fine
  assignments identical to r8 for the 5,091 non-puromycin genes; numeric
  columns within 1.1e-9). The puromycin section of ModeConservation is gone.
- Frozen manuscript env removed from the repo: `env/manuscript-*`,
  `docs/user/manuscript_reproduction.md`, the flat `chronopack/` package
  (it shadowed `src/` when Python ran from the repo root) and the top-level
  shims `chronopack.{c2c,tradeseq,grn_inference,omnipath_tf,preprocessing,trajectories}`.
- Verification run r9 (`full_pipeline_analyses_20260930_r9`) vs r8: all stages OK
  and promoted into `analyses/`.
  - Preprocess / GetTrajectories / PrepareForTradeSeq: identical cells, counts,
    layers, neighbour graphs, traj_clusters, fate probabilities and every tradeSeq
    input table. Node-dependent only (as r8 vs r7): UMAP / ot_umap, scrublet
    doublet_score, cosine top-k `obsp/S` (brute NN with n_jobs=-1, tie order).
  - TradeSeqFits: K562+TNF (bigcompute, 160G; 90G failed with "unable to fork")
    exactly identical to r8; LPS K562/THP1 same genes, significant sets and order,
    stats within ~1e-8 relative.
  - Gene clustering: all broad/fine assignments identical to the puromycin-free
    reference (`dimred_without_puromycin_20260930`); numeric within 4e-8.
  - `02_ClusterVisualization` / `03_ModeConservation` were stale (read the old
    genes-in-var layout and `PC1`/`PC2`, never run in r8); now read `adata.obs`,
    and the feature-space plot uses `t_peak` vs `peak_log2fc`.
- After r9 (not rerun, by user decision): removed `cp.io.annotate_qc_metrics`,
  `filter_cells_by_qc`, `initialize_expression_layers` and
  `validate_preprocessed_adata` (plus its helper `validate_numeric_obs`) from the
  package; the two `01_PreprocessDatasets` notebooks no longer call the validator
  (check-only, so saved outputs are unaffected). The `chronopack-install-test`
  wheel (`dist_fix4`) predates this; rebuild it before the next verification run.
- Also after r9 (not rerun): raw inputs copied (byte-identical) into
  `analyses/inputs/` (gitignored: 8 DGE tables + `k562_neg-ctrl_mtx-agg.pkl`) and
  reference files into `analyses/resources/` (CellChatDB, `hs_tfs.txt`,
  `protein_coding_genes.csv`). Notebooks read them via `../inputs/…` /
  `../resources/…`; `CHRONOPACK_DATA_DIR` now only receives outputs, and
  `run_pipeline.py setup` no longer symlinks anything. c2c/grn still read older
  tradeSeq intermediates from `~/projects/…` and `~/tmp/grn_inference_v2/…`.
- `_legacy/` dissolved (not rerun): `c2c.py` -> `tl/c2c/_helpers.py` (compute) +
  `pl/c2c/_plots.py` (plots); `omnipath_tf.py` -> `tl/grn/_omnipath.py` (TF list /
  LR table now required args, no `/new-stg` defaults); `grn_inference.py` ->
  `tl/grn/_edges.py`. Removed unused: `_legacy/{preprocessing,trajectories,tradeseq}.py`,
  `cp.io.write_tradeseq_inputs`, `_utils/`, and the uncalled helpers
  `score_lr_pairs`, `plot_avg_vs_logp`, `load_secreted_cellchat_pairs`,
  `filter_expressed_secreted_lr_catalog`, `all_expressed`, and the duplicate
  `load_chronopack_beta_lfdr` / `build_` / `write_chronopack_ranked_edges` in the
  old edge module.
- C2C / GRN now read pipeline outputs via `CHRONOPACK_DATA_DIR` and write to
  `DATA_DIR/{c2c,grn}` (runner stages `CellCommunication`, `GRNInference`):
  - C2C: r9-style `adata_full.chk1.h5ad` (identical results to the Jan 21 file) and
    the THP1 screen-fit association (`pvalue_1`; `infer` constant changed from
    `pvalue_2`) instead of the joint two-trajectory fit: same top 10 ligands and
    labeled pair, 100 vs 189 pairs (ligands < 10 THP1 counts untested).
  - TNF overlay: `cp.tl.predict_gene_dynamics` on the THP1/K562 screen SCEs
    (user: "option b for now", expects to change it). `load_tradeseq_grid_exports`
    removed.
  - GRN: `score_targets` refits the BUM local FDR from `pvalue_1`, so r9's
    `asso_df.1trajectory.nknots4.csv` is used directly (R loader fixed for an empty
    first header / no `lfdr` column). Shared-gene ordering identical; more targets
    pass lfdr < 0.05 (smaller tested universe). Comparison:
    `verification_runs/c2c_grn_r9_inputs_20260930/`.
  - Install-test wheel rebuilt: `install_test_20260925/dist_fix5/`.
- `fit_gene_dynamics` / `prepare_tradeseq_inputs` lost `full_adata`: library sizes
  always come from `adata.obs["total_counts"]` (unchanged by gene filtering) and the
  unused `*.full.*` input copies are no longer written. `predict_gene_dynamics` fixed
  to attach SingleCellExperiment (rownames() was empty otherwise). Wheel: `dist_fix6`.
- User's simplified `03_PrepareForTradeSeq` notebooks (both experiments) verified
  against r9 GetTrajectories outputs (`verification_runs/prepare_edits_check_20260930_b/`):
  prefiltered h5ad, counts, cell times, gene lists, library sizes, memberships all
  identical; exclusion list same 1025 genes (set order differs). C2C and GRN
  notebooks run on r9 data (`r9/data/{c2c,grn}`), matching the comparison job
  exactly. These four notebooks were promoted with outputs; TradeSeqFits onward
  not rerun (inputs identical).
- 2026-10-01: soft memberships modularized (`cp.tl.trajectory_soft_memberships`);
  K562+TNF `03_PrepareForTradeSeq` writes `cell_weights.csv` (identical to the old
  `cluster_memberships.retentions.csv`); new `04_TradeSeqFits.TwoLineage.ipynb`
  (`lineage_weights=`, runner stage `TradeSeqFitsTwoLineage`, 160G bigcompute)
  verified on r9 inputs (`verification_runs/two_lineage_check_20260930/`: 11864
  genes tested, 1830 refit). GRN ties now broken explicitly by `waldStat_1` (same
  edges as before). C2C overlay uses ad hoc 8-knot fits (`fit_gene_dynamics(genes=)`)
  matching the full refit to 1e-14 (`verification_runs/c2c_adhoc_fit_check_20261001/`).
  Wheel: `dist_fix8`.

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


### 2026-09-25: repo-root notebooks refactored into public-API clients

Goal (user): repo-root experiment notebooks should look like a normal user
workflow (`import chronopack as cp`, public `cp.io` / `cp.tl` / `cp.pl` only),
keep the Sep 21 cleanups, include this session's maintained-run fixes, and
behave the same as the maintained run. Guidance followed:
`docs/developer/notes/mixed_R_python_package_architecture.md` §12,
`docs/developer/notes/ChronoPack_package_refactor_guide.md` §3/§5/§11-12,
`docs/developer/notes/additional_notes_on_mixed_python_R_architecture.md` §3-5/§26.

- Package (tests: 63 passed in both envs; only the pre-existing cp.io export
  test fails): new public `cp.tl.identify_tag_bias_candidates`,
  `cp.tl.filter_tag_bias_candidates`, `cp.tl.tag_bias_exclusion_genes`,
  `cp.tl.merge_trajectory_graphs` (per-cell-type WOT graphs -> one graph),
  `cp.tl.standardize_gene_trajectories` (fit / fit_z / ddt_z exports) and
  `cp.pl.cell_trajectories.fate_probabilities` (still writes the fate_* obs
  columns, as before).
- Notebooks (15 pipeline notebooks): removed `sys.path` repo walk-ups, private
  imports (`chronopack.tl._trajectories`, `chronopack.tl._gene_dynamics`) and the
  `chronopack.preprocessing` shim; no function definitions remain; kernelspec
  set to `python3`; restored `cp.io.validate_preprocessed_adata` checks. The 4 R
  TradeSeqFits notebooks are now Python notebooks calling
  `cp.tl.fit_gene_dynamics(input_dir=..., backend="tradeseq", ...)`, which runs
  the same bundled `run_tradeseq_single_lineage`. K562+TNF Preprocess applies
  the manuscript cluster removals by barcode
  (`k562-tnfa.hg38.rep/exclude_leiden_clusters.manuscript.csv`, whitelisted in
  `.gitignore`). DimRed uses explicit SCE paths instead of the
  `*_paths_from_dimred_roots` helpers and no association-CSV migration shim.
- Bugs fixed in the Sep 21 notebooks (never run end to end before):
  `fate_umap_2x2_scanpy(..., embedding_basis=...)` (wrong keyword) and
  `basis='X_ot_umap'` in K562+TNF PrepareForTradeSeq after the `ot_umap` rename.
- Known deviations kept on purpose: `cp.pl.cell_trajectories.mass_retention`
  and `fate_probabilities` write obs columns used downstream (plotting mutates
  state); THP1 cell_type and trajectory-cluster names still come from Leiden
  label numbers; `DimensionalityReduction&Clustering.cluster_visualization.ipynb`
  expects genes in `.var` of `gene_trajectory_clusters.h5ad` and a different
  output dir (stale vs current DimRed); only its kernelspec was changed.
- Verification: env `chronopack-install-test` = clone of `chronopack` with a
  fresh non-editable wheel (`verification_runs/install_test_20260925/`).
  New run `full_pipeline_maintained_chronopack_20260925_r8` via
  `verification_runs/submit_pipeline_user_api.py`, compared stage by stage with
  the maintained r7 run (`verification_runs/compare_runs.py`, logs in
  `verification_runs/compare_logs/`). Every stage job confirmed ChronoPack was
  imported from the installed wheel in site-packages.
  - PreprocessDatasets (933386-933388): puromycin and THP1/K562 identical to r7;
    K562+TNF same cells/genes/order, counts/cpm and kNN graphs identical, other
    values within ~1e-13. K562+TNF `doublet_score` differs from r7 by up to
    0.014 but is identical to the Sep 19 manuscript template (r6/r7 differ from
    the template; likely the in-place subset in `cp.io.filter_cells_by_qc`).
    `predicted_doublet` identical everywhere.
  - GetTrajectories (933394-933396): `traj_clusters`, fate probabilities,
    terminal-cell distributions and mass retention identical to r7 for all
    three experiments. K562+TNF `S` graph equals the template exactly; r6/r7
    differ from the template in 4 edges over 2 cells (near-tied top-k).
    `X_umap` differs by design (Sep 21 notebooks keep the Scanpy UMAP there).
  - PrepareForTradeSeq (933409-933411): all tradeSeq inputs for the four
    branches have identical genes, cells and parsed values to r7 (counts are
    now integer text; K562+TNF full counts 9.8 MB vs 11 MB).
  - TradeSeqFits as Python notebooks via `cp.tl.fit_gene_dynamics(backend="tradeseq")`
    (933442-933445): all four branches identical to r7 with max diff exactly 0
    (association stats, significant sets and order, nknots=8 SCE coefficients
    and rowData); screen SCEs present, no duplicate SCE copies. BiocParallel
    forking under rpy2 works; runtimes match r7 (THP1 ~30 min, K562+TNF ~2 h).
  - DimRed first attempt (933703) failed: `cp.tl.align_gene_trajectories` passed
    `association_paths=None` to R, which rpy2 cannot convert (this path never
    ran before because the old notebooks always passed the CSV shim). Fixed in
    `src/chronopack/tl/_gene_clustering.py` (argument omitted when None; the R
    default is NULL). Rebuilt wheel `install_test_20260925/dist_fix1/` installed
    into `chronopack-install-test`; recorded in the run's `code_provenance.txt`.
    Stages 1-4 used the first wheel (the change doesn't affect them). Rerun: job 933975.
  - DimRed (933975): passed. Six broad counts match expected; `pc_df.new.csv`
    has the same 6,725 rows, order and 31 columns as r7; every cluster column
    (broad and fine, induced and repressed) and all association statistics are
    identical; K562_puro rows all -1. Numeric descriptors and centroid distances
    differ from r7 by at most ~1e-10 (floating point); exported fit / z-score
    matrices and cluster metrics within 7e-11.
  - Conclusion: the refactored repo-root notebooks, run against a normally
    installed ChronoPack, reproduce the maintained r7 pipeline at every stage.
  - `tests/test_pre_trade_refactor.py` updated to the new notebook structure
    (public API only; Python TradeSeqFits; validation checks present) plus a new
    test that pipeline notebooks use only the installed public API. 64 passed;
    the only failure is still the pre-existing cp.io export test.

### 2026-09-28: GRN notebooks combined

- `grn_inference/GRNInference.ipynb` replaces the three GRN notebooks (moved to
  `grn_inference/old/`): one context table drives `cp.tl.grn.rank_regulators`,
  `cp.tl.grn.score_targets` and `cp.tl.grn.infer`, with results passed in
  memory; python3 kernel; outputs and parameters unchanged.
- Verified (`verification_runs/grn_inference_notebook_20260928/`, job 934103,
  installed wheel): OmniPath TF rankings and all lfdr tables identical to the
  Sep 17 `grn_inference/outputs`. Edges: 5,009,770 vs 4,334,306 in
  `regulator_prioritization_chronopack_ranked_edges.csv`; every reference edge
  is present with identical EdgeWeight, and restricting to the reference's
  regulators (176 K562 / 189 THP1) reproduces it exactly, rank order included
  (job 934104). The reference was built with a benchmark TF-universe filter,
  which the active workflow intentionally omits.
- Target-gene inputs still come from the external snapshot
  `/new-stg/home/alex/tmp/grn_inference_v2/reproduction/...`, not from the
  verified pipeline outputs.

### 2026-09-28: notebooks runnable from Jupyter with the `chronopack` env

- ChronoPack is now installed editable (`pip install --no-deps -e .`) in
  `/new-stg/home/alex/miniconda3/envs/chronopack`; no other packages changed.
- New `chronopack._backends.r.ensure_r_home()` (called when the R backend
  package is imported, before any rpy2 import) sets `R_HOME` to
  `<sys.prefix>/lib/R` when `R_HOME` is unset and no `R` is on `PATH`. Fixes
  "Unable to determine R_HOME" in Jupyter kernels that start an env's Python
  without activating the env (e.g. `~/.local/share/jupyter/kernels/chronopack`).
  Tests: `tests/test_r_home.py`; full suite 68 passed both with and without an
  activated env (only the pre-existing cp.io export test fails).
- `chronopack-install-test` reinstalled from `install_test_20260925/dist_fix2/`, then
  (2026-09-30, after the C2C/GRN API changes) from `dist_fix3/`. That env has a
  copied wheel, so it needs a rebuild after package edits; `chronopack` is editable
  (only a Jupyter kernel restart is needed).

### 2026-09-30: C2C simplified

- `cp.tl.c2c.infer(*, adata, lr_resource, association, sender, receiver, min_cells=1)`.
  Removed arguments (now fixed module constants in `tl/c2c/_inference.py`):
  cell_type/time keys, counts/log1p layers, 0-45 min window, `pvalue_2`,
  pseudocount, secreted-only, fits/tradeSeq exports, provenance. Removed
  `response_genes` (scores never used the receiver response); `temporal_profiles`
  now has LR activity only. `min_detect=1e-100` replaced by `min_cells` (each
  ligand / receptor gene detected in at least that many sender / receiver cells;
  default 1 is identical). Clear ValueError when no pairs remain.
- `cp.pl.c2c.pair_overlay` kept; takes `response=` (time-indexed Series).
- `c2c/CellCommunication.ipynb` rewritten: load, one short `infer` call with
  sender/receiver inline, `write`, three plots. No longer writes the legacy
  `chronopack_lr_pair_scores*.csv` / `chronopack_ligand_scores.csv` copies.
  With user approval, `C2CResult.legacy_pair_scores` was removed and those
  three stale CSVs were deleted from `c2c/outputs/`.
- Verified: same-node old vs new `infer` identical in every table (job 934409);
  notebook run (job 934403) gives the same 73-ligand ranking as the Sep 28
  outputs, with `logp`-derived values differing by <=1.1e-16 (different CPU).
  Tests: 69 passed (only the pre-existing cp.io export test fails).

### 2026-09-30: C2C window parameter, summed activity, reference scatter style

- `cp.tl.c2c.infer(..., causal_window=(0.0, 45.0))` replaces the
  `CAUSAL_WINDOW_MIN` constant (metadata key `causal_window`).
- Pair activity is now the **sum** over window time points (`activity_sum`),
  not the mean; `avg` and the duplicate `score` columns are gone. Scoring lives
  in `_score_pairs` in `tl/c2c/_inference.py`; legacy `score_lr_pairs` (mean
  `avg`) is untouched for the frozen notebooks.
- Every pair has 4 window points, so sum = 4 x mean: pair order, pair ranks,
  ligand order and every ligand rank identical to the mean-based outputs;
  z-scores / composite differ by <=4.4e-16 (job 934421,
  `verification_runs/c2c_sum_window_20260930/`). An old mean threshold of 0.1
  is 0.4 on the sum.
- `cp.pl.c2c.interaction_scatter` restyled to match
  `c2c_downsample/reproduction/.../candidate_lr_ligand_significance_vs_activity.pdf`
  (Reds by rec_det clipped at 0.2, size-300 bubbles, dashed p=0.05 / `x_thresh`
  lines, boxed labels for pairs past both). Notebook call:
  `interaction_scatter(c2c, x_thresh=0.4, xlim=(-0.6, 1.2))` = reference axes x4.
- Install-test wheel rebuilt as `install_test_20260925/dist_fix4/`.

### 2026-09-30: GRN API and notebook are per-context

- `cp.tl.grn.rank_regulators(source=, is_ligand=, human_tf_path=, ligand_receptor_path=)`
  returns one TF ranking DataFrame (`TFRankingResult` removed; diffusion
  history no longer returned). `score_targets(association, valid_targets=)`
  scores one tradeSeq table (no `dataset` column). `infer(targets, tfs,
  valid_targets=, name=, tf_universe=None)` builds one context's edges; the
  column-name arguments are fixed constants. `build_chronopack_ranked_edges`
  now calls `infer` per context (same output).
- `grn_inference/GRNInference.ipynb`: shared resource paths once, then a
  K562 + TNF-α section and a THP1 + LPS section (rank TFs, score targets,
  infer, write to `outputs/<context>/`). No context loop.
- Verified (job 934416, 152 s, unactivated env like a Jupyter kernel) against the
  previous combined run: TF rankings, lfdr summary/params, edges, regulators and
  target scores used by `infer` identical for both contexts. The per-gene
  `target_scores.scores.csv` differs only in CSV float text (<=3.4e-16 rel with
  pandas' default parser; the old file was R's 15-digit `write.csv`).
- Tests: 71 passed (new single-context and `score_targets` R tests); only the
  pre-existing cp.io export test fails. Docs updated.

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


### 2026-10-01: gene-clustering notebook refactor

- `01_DimensionalityReduction&Clustering` rewritten around one `traj` object:
  `compute_gene_trajectory_descriptors(..., inplace=False)`; writes only
  `gene_trajectory_characterization.h5ad`, `gene_trajectory_clusters.h5ad`, the two
  `cluster_metrics.*.csv` and metric plots (no `pc_df.*`, `gene_trajectory_fit*.csv`,
  `trajectory_payload/`, `cluster_candidates/`); separate induced / repressed sections
  (cluster -> metrics plot -> pick k -> labels into `traj.obs`), then one combine cell.
- Package: `output_dir` removed from `align_gene_trajectories`,
  `compute_gene_trajectory_descriptors`, `cluster_gene_trajectories` (R hand-off now
  always temporary); unused `_default_dimred_paths` removed; `align_gene_trajectories`
  contexts default to the `sce_dict` keys (was a fixed tuple incl. `K562_puro`, which
  rejected other context names). `compare_runs.py GeneClustering` reads
  `gene_trajectory_clusters.h5ad` obs (falls back to `pc_df.new.csv`).
- Verified on r9 fits (`verification_runs/gene_clustering_refactor_check_20261001/`):
  all broad/fine assignments identical, numeric <= 9e-10, mode conservation identical.
  Wheel: `dist_fix9`.
- Removed the old three-experiment layout API: `cp.tl.sce_paths_from_dimred_roots`,
  `cp.tl.association_paths_from_dimred_roots`, `dimred_config`, R
  `chronopack_dimred_config`, and `align_gene_trajectories(root_puro=, root_k562a=,
  root_lps=)` (`sce_dict` is now required); `K562_puro` dropped from R default maps.
- `traj.obsm["X_pca"]`: the PCA of trajectory shapes used for clustering (4 components;
  fitted on the reference = all significant trajectories, so induced and repressed
  clusterings share it). `cp.pl.trajectory_features(..., basis="X_pca")` plots it;
  `02_ClusterVisualization` uses it. Clusters unchanged vs r9. Wheel: `dist_fix10`.
- R-rendered plots now display inline: `_render` also draws a PNG preview, which
  `RPlotResult._repr_png_` returns (SVG outputs are not shown by Jupyter in untrusted,
  script-executed notebooks); `RPlotResult` repr is one line (path only); previews display at half their pixel width (`display_scale=0.5`).
  `verification/execute_notebook.py` records PNG/SVG representations.
  `02_ClusterVisualization` / `03_ModeConservation` end cells with the plot object
  (or `display(...)`) and were promoted with 10 / 2 inline plots. Wheel: `dist_fix13`.

### 2026-10-01: TESTING.ipynb plotting port (gene clustering)

- R (`dimred_plotting_helpers.R`): `plot_cluster_medoids_and_top_genes` and
  `plot_cluster_pair_representative_genes` replaced by the TESTING.ipynb versions
  (pair function calls `plot_tradeSeq_gene_sce_pub` = TESTING's SCE-only
  `plot_tradeSeq_gene`, since TESTING's `plot_tradeSeq_gene_multirun_pub` was never
  defined; silhouettes via `chronopack_silhouette`, which falls back to the same
  definition because `cluster` is not installed in the chronopack envs).
  `calc_topk_mode_recovery` added verbatim to `dimred_clustering_helpers.R`.
  New `chronopack_plot_*_csv` entry points draw to the output file and a PNG preview in
  one call; old cluster-trajectory / transition-table Sankey and chord renderers removed.
- Python: `cp.pl.trajectory_features` (now `feature_plot_pc_simple`),
  `cluster_medoids_and_top_genes`, `mode_proportions`, `mode_sankey` / `mode_chord`
  (AnnData + context pairs, TESTING renderers), `cluster_pair_representative_genes`,
  `cluster_violin`; `cp.tl.topk_mode_recovery`. Removed `cp.pl.cluster_trajectories`.
  Kept `cp.tl.mode_conservation` / `mode_switchers` / `cp.pl.mode_conservation` (unused by
  notebooks). `GeneTrajectoryClusteringResult.features` = joint clustering features; 01
  stores them in `traj.obsm["X_cluster_features"]`.
- Notebooks: 01 has prospective broad-mode medoid plots (induced / repressed); 02 rebuilt
  (PCA, broad + fine medoid panels, pies, 3 Sankeys, top-k recovery, chord, K562a x THP1
  representative genes, impulse-up violins); `03_ModeConservation` removed.
- Verified (Slurm 934596, `gene_clustering_refactor_check_20261001/`): all cluster
  columns identical to r9 (n_diff 0); no errors; every plot cell has a PNG (01: 4, 02: 14).
  Top-k recovery = TESTING saved values (induced 0.61 / 0.35 / 0.26, repressed 0.11 /
  0.06 / 0.09). Executed 01/02 promoted. Wheel: `dist_fix14`.

### 2026-10-01: gene-clustering notebooks merged; mode-conservation API removed

- `gene_clustering/01_DimensionalityReduction&Clustering` + `02_ClusterVisualization` merged
  into `gene_clustering/01_GeneTrajectoryClustering` (one setup cell, figures use the
  in-memory `traj`; same output paths). New cell after the pies: paired modes grouped
  (impulse-up, trend-up, trend-down, transient-down) as percent of responsive genes per
  context and min-max range (13-25 / 24-41 / 44-52 / 2-6%). `run_pipeline.py`
  `GeneClustering`, tests and README updated.
- Removed `cp.tl.mode_conservation` / `mode_switchers` / `ModeConservationResult`,
  `cp.pl.mode_conservation` and R `chronopack_plot_mode_conservation_csv` (with their tests);
  `tl/_mode_conservation.py` keeps `MODE_DIRECTIONS` and `topk_mode_recovery`.
- Verified (Slurm 934640, `gene_clustering_refactor_check_20261001/`): cluster columns
  identical to r9 (n_diff 0); no errors; 17 inline plots (4 + 13, as in the previous 01/02);
  text outputs of the moved cells equal the previous 02. Promoted. Wheel: `dist_fix16`.
- 2026-10-01: `cp.pl.gene_dynamics_curves` (matplotlib grid of fitted NB-GAM curves, one
  panel per gene, one curve + 95% band per group) added; K562+TNF `05_Visualization`
  ends with the two-lineage curves (`predict_gene_dynamics(..., lineage=1/2,
  sce_file="sce.2trajectory.nknots8.rds", time_max=110)`) for NFKBIA, NFKB2, ICAM1, BIRC3 /
  CXCL3, CD69, IFNGR2, DDIT4. Re-executed on r9 + two-lineage check data (7 figures, no
  errors) and promoted.

### 2026-10-01: THP-1 + K562 + LPS `05_Visualization`

- New `analyses/thp1-k562-lps.hg38.rep/05_Visualization.ipynb` (also in `run_pipeline.py`
  `Visualization` as a separate job, README, notebook tests):
  1. CreateFigures.ipynb final figure (endpoint UMAP, cells 3 + 8): new
     `cp.pl.cell_trajectories.trajectory_endpoints`; `cp.tl.rotate_embedding_to_time(groupby=)`
     rotates each cell type about its own centroid (CreateFigures cell 3, "preserve" mode).
     Input `adata_hvg.chk1.h5ad` (original: `adata_both_prefiltered.h5ad`, not produced by the
     pipeline; same cells and identical `traj_clusters` per barcode: 0 = NF-kB hi, 1 = lo, 2 = THP1).
  2. VolcanoPlot.ipynb first volcano (K562b vs THP1): a different function from the K562
     volcano, ported verbatim as R `plot_response_specificity_volcano` +
     `cp.pl.response_specificity_volcano`; `make_combined_from_pc_df` ported as
     `cp.tl.combine_context_associations`. Original inputs (`pc_df.rds` + `desc.pc_1e-1.rds` of the
     old multi-experiment global analysis) replaced by the two `asso_df.1trajectory.nknots4.csv`
     with `signed_max_log2fc_from_start(pseudocount=0.1)`; the original's gene-exclusion step was a
     no-op (0 overlaps) and is dropped; no puromycin data used. Original arguments kept (y_break
     default = 95% quantile logic). 28 of the original's 33 labels recur.
  3. NB-GAM curves (THP1 / K562b, 8 knots, 0-165 min). Ad hoc fits: THP1 IRF1
     (`fit_gene_dynamics(input_dir=THP1 fit, genes=...)`) and THP1 CD69, which is not a THP1
     tradeSeq input (8 UMIs, low-count filter): fitted from the THP1 cells of
     `adata_full.chk1.h5ad` prepared into a temporary folder (route checked: identical to the
     folder fit for IRF1/NFKBIA, max diff 0).
- Verified (Slurm 934657, `verification_runs/thp1_lps_visualization_check_20261001/`, per-file
  r9 symlinks): no errors, 3 inline PNGs; promoted. Tests 95 pass. Wheel: `dist_fix17`.

## Update (2026-10-01): NB-GAM composites in the GATA1 (`plot_tradeSeq_gene_multirun_pub`) style

- New `cp.tl.observed_gene_dynamics(fit_or_input_dir, genes, *, lineage, lineage_weights, time_max)`:
  per-cell `log1p(count * s_ref / sf)` (same reference library size as `predict_gene_dynamics`),
  `t` in minutes, `weight` (normalized lineage weight, as in the backend) and `assigned`
  (highest-weight lineage).
- `cp.pl.gene_dynamics_curves(..., observations=...)`: cell points (group color, alpha 0.3,
  rasterized), band + line, white-filled per-time-point means (weighted mean of normalized
  counts, then log1p, as in the R reference); minimal theme, x ticks every 30 min.
- Both 05_Visualization notebooks updated (K562: points colored by highest-weight lineage of
  `cell_weights.csv`, means weighted per lineage; THP1/K562b: observations from each fit folder,
  THP1 CD69 from the temporary input folder). Checks: Slurm 934686 (K562; `cell_weights.csv`
  symlinked from `two_lineage_check_20260930`, counts identical to r9) and 934685 (THP1): no
  errors; promoted. Tests: 97 passed. Wheel: `dist_fix18`.
