# CellCommunication / C2C Migration Mapping

Status: package-backed CellCommunication/C2C inference slice complete.

This pass migrated reusable temporal cell-cell communication inference and
direct reusable plotting into `cp.tl.c2c` and `cp.pl.c2c`. It did not migrate
C2C benchmarking, comparator methods, gold standards, benchmark metrics, or
benchmark-specific preprocessing/plots.

## Public API

- `cp.tl.c2c.infer(...)`
  - loads or accepts a ligand-receptor table;
  - splits sender and receiver populations from AnnData;
  - filters LR pairs by sender ligand and receiver receptor expression
    (`min_cells`, default 1: each gene detected in at least that many cells);
  - restricts ligands to genes with association rows;
  - builds temporal LR activity profiles (no receiver response module; the
    scores never used it);
  - scores early-window LR activity (summed over the time points in
    `causal_window`, default `(0, 45)` minutes) and ligand temporal
    significance;
  - returns `C2CResult`.

`cp.tl.c2c` also exports the helpers the C2C notebook uses for the TNF callout:
`split_by_cell_type` and `build_pair_series` (implementation in
`tl/c2c/_helpers.py`); `cp.pl.c2c` exports `plot_pair_overlay` and
`plot_genes_over_time` (`pl/c2c/_plots.py`). The callout's fitted curves come
from `cp.tl.predict_gene_dynamics` on the THP1 and K562 screen fits
(provisional; it replaced the joint two-trajectory grid export).
`load_secreted_cellchat_pairs`, `filter_expressed_secreted_lr_catalog`,
`score_lr_pairs`, `plot_avg_vs_logp` and `load_tradeseq_grid_exports` were
removed on 2026-09-30.

Inputs (since 2026-09-30) are pipeline outputs: the THP1 + K562 + LPS
`adata_full.chk1.h5ad` and the THP1 screen-fit association table
(`sc-THP1_to_tradeSeq.fixed/asso_df.1trajectory.nknots4.csv`, `pvalue_1`),
replacing the joint THP1+K562 two-trajectory fit (`pvalue_2`). Checked against
the joint fit: same top 10 ligands and labeled pair; 100 vs 189 candidate pairs,
because ligands with < 10 THP1 counts are not tested in the THP1 fit
(`verification_runs/c2c_grn_r9_inputs_20260930/`).

## Result Object

`C2CResult` contains:

- `interactions`: primary LR pair score table with sender/receiver context,
  ligand, receptor, pair, `activity_sum` (early-window summed activity),
  receptor detection `rec_det`, ligand `logp`, early-window point count, and
  pair rank;
- `ligand_scores`: ligand-level aggregate priority table;
- `temporal_profiles`: long-form LR activity and receiver-response profiles;
- `candidate_interactions`: filtered LR catalog used for scoring;
- `detection`: sender and receiver detection Series;
- `parameters`;
- `provenance`.

The full communication network is not stored in AnnData because LR
interactions are not aligned to either AnnData axis.

`C2CResult.write()` exports CSV tables plus JSON metadata and attempts Parquet
for the interaction table when a Parquet backend is available.

## Preserved Temporal Method

The migrated workflow preserves the current method exactly:

- LR resource: CellChat-style ligand/receptor table, optionally restricted to
  `annotation == "Secreted Signaling"`;
- ligand/receptor component parsing, including receptor shorthand expansion
  such as `TGFBR1_R2`;
- expression eligibility: every ligand / receptor gene detected (count > 0) in
  at least `min_cells` sender / receiver cells (default 1, identical to the
  former `min_detect=1e-100`);
- ligand restriction to the association table index;
- LR activity:
  - ligand time series from sender expression;
  - receptor availability from receiver receptor expression;
  - activity as the existing square-root expression/availability transform;
- score:
  - subtract baseline activity at the first time point;
  - sum over the time points within `causal_window` (default `0-45` minutes;
    time tags map to `(tag - 1) * 15` minutes, so the default window has 4
    points for every pair);
  - combine with ligand temporal significance `-log10(p + 1e-3)`;
- ligand aggregation and z-score composite ranking as in the notebook.

The pair activity is now a sum (`activity_sum`) instead of the historical mean
(`avg`). Every pair has the same number of window points, so the sum is the
mean times that number: pair ranks are unchanged, and ligand ranks are
unchanged because `max_lr_activity` enters the composite only through a
z-score. Only the activity scale (x4 for the default window) changes, so an
activity threshold of 0.1 on the old mean corresponds to 0.4 on the sum.
The historical mean-based `score_lr_pairs` helper has been removed.

## Current Workflow Mapping

| Historical code / notebook block | Purpose | Package location | Public? |
| --- | --- | --- | --- |
| CellChat CSV loading and secreted filtering | LR resource loading | `cp.tl.c2c.infer()` / helper | public task behavior |
| `split_by_cell_type()` | sender/receiver split | `cp.tl.c2c` helper | temporary public helper |
| `filter_expressed_secreted_lr_catalog()` | expression eligibility | `cp.tl.c2c.infer()` internals (`_expressed_pairs`) | public task behavior |
| `restrict_ligands_to_association()` | target ligand universe | `cp.tl.c2c.infer()` internals | public task behavior |
| `time_series()` / `module_series()` | temporal expression extraction | `cp.tl.c2c` helper | temporary public helper |
| `lr_activity()` / `build_pair_series()` | LR activity + receiver response profiles | `cp.tl.c2c.infer()` internals | public task behavior |
| `score_lr_pairs()` | early-window pair and ligand scoring | `cp.tl.c2c.infer()` internals (`_score_pairs`) | public task behavior |
| LR scatter plot | pair score visualization | `cp.pl.c2c.interaction_scatter()` | public plotting |
| LR/response overlay | manuscript callout helper | `cp.pl.c2c.pair_overlay()` for result profiles; legacy helper for custom fitted callouts | public plotting/helper |
| selected genes over time | reusable raw-expression plot | `cp.pl.c2c.plot_genes_over_time()` | public plotting |
| TNF/NFKBIA callout | specific biological example | active notebook | notebook-local |

## External Resources

External resources remain caller-provided:

- CellChatDB-style ligand-receptor table;
- AnnData object containing expression, cell type, and time metadata;
- the sender's single-lineage tradeSeq association table;
- the THP1 and K562 tradeSeq screen fits for the TNF callout curves.

No large ligand-receptor resource was bundled into ChronoPack.

## Benchmark Boundary

Cell-cell communication benchmarking is intentionally maintained outside
ChronoPack. ChronoPack provides communication inference, temporal analysis,
standardized results, and reusable visualization. Benchmark datasets,
comparator methods, evaluation metrics, gold standards, and benchmark-specific
filtering remain separate research workflows.

No `cp.tl.c2c.benchmark` namespace was added.

## Active Notebook

`c2c/CellCommunication.ipynb` is now a thin package client. It retains:

- dataset paths and resource paths;
- one `cp.tl.c2c.infer` call with sender/receiver written inline;
- package calls;
- output saving;
- manuscript-specific TNF and selected-gene plots.

Reusable inference logic now lives in `cp.tl.c2c`; reusable result-driven plots
live in `cp.pl.c2c`.

## Regression Verification

Frozen upstream inputs were used; preprocessing, WOT, tradeSeq, gene clustering,
GRN, and benchmarking workflows were not rerun.

Representative full-data regression rebuilt the active C2C workflow and
compared package-backed outputs to the existing notebook outputs:

- pair scores:
  - shape: `189 x 8`;
  - no missing LR pair identities;
  - maximum numeric absolute difference: `1.1102230246251565e-16`;
  - nonnumeric columns identical.
- pair scores with rank:
  - shape: `189 x 8`;
  - no missing LR pair identities;
  - maximum numeric absolute difference: `4.440892098500626e-16`;
  - nonnumeric columns identical.
- ligand scores:
  - shape: `73 x 11`;
  - maximum numeric absolute difference: `1.7763568394002505e-15`;
  - nonnumeric columns identical.

`C2CResult` full-data payload:

- candidate interactions: `189 x 14`;
- interactions: `189 x 11`;
- ligand scores: `73 x 11`;
- temporal profiles: `2268 x 8`.

Regression artifacts:

- `verification_runs/c2c_package_migration_20260918/`

## Tests

Updated `tests/test_c2c.py`.

Synthetic tests cover:

- ligand/receptor parsing;
- stable z-score behavior;
- pair and ligand scoring;
- package-backed `cp.tl.c2c.infer()`;
- absent receptor filtering;
- missing sender validation;
- `C2CResult` serialization;
- `cp.pl.c2c` rendering smoke.

Focused C2C test slice: `6 passed`.
