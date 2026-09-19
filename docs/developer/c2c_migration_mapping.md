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
  - filters LR pairs by sender ligand and receiver receptor expression;
  - restricts ligands to genes with association rows;
  - builds temporal LR activity profiles against a receiver response module;
  - scores early-window LR activity and ligand temporal significance;
  - returns `C2CResult`.

Reusable lower-level helpers remain available temporarily in `cp.tl.c2c` and
the top-level `chronopack.c2c` compatibility shim, including:

- `load_secreted_cellchat_pairs`;
- `load_tradeseq_grid_exports`;
- `filter_expressed_secreted_lr_catalog`;
- `build_pair_series`;
- `score_lr_pairs`;
- `plot_pair_overlay` for notebook-local callouts.

## Result Object

`C2CResult` contains:

- `interactions`: primary LR pair score table with sender/receiver context,
  ligand, receptor, pair, score, historical `avg`, receptor detection,
  ligand `logp`, early-window point count, and pair rank;
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
`C2CResult.legacy_pair_scores()` reproduces the historical pair-score CSV
schema used by the active notebook.

## Preserved Temporal Method

The migrated workflow preserves the current method exactly:

- LR resource: CellChat-style ligand/receptor table, optionally restricted to
  `annotation == "Secreted Signaling"`;
- ligand/receptor component parsing, including receptor shorthand expansion
  such as `TGFBR1_R2`;
- expression eligibility by nonzero detection in sender/receiver cells;
- ligand restriction to the association table index;
- LR activity:
  - ligand time series from sender expression;
  - receptor availability from receiver receptor expression;
  - activity as the existing square-root expression/availability transform;
- receiver response module as the existing average response-gene time series;
- score:
  - subtract baseline activity at the first time point;
  - average within the configured causal window, currently `0-45` minutes;
  - combine with ligand temporal significance `-log10(p + 1e-3)`;
- ligand aggregation and z-score composite ranking as in the notebook.

No switch was made between sum, integral, or area-under-curve definitions; the
historical early-window mean activity remains `avg`.

## Current Workflow Mapping

| Historical code / notebook block | Purpose | Package location | Public? |
| --- | --- | --- | --- |
| CellChat CSV loading and secreted filtering | LR resource loading | `cp.tl.c2c.infer()` / helper | public task behavior |
| `split_by_cell_type()` | sender/receiver split | `cp.tl.c2c` helper | temporary public helper |
| `filter_expressed_secreted_lr_catalog()` | expression eligibility | `cp.tl.c2c.infer()` internals | public task behavior |
| `restrict_ligands_to_association()` | target ligand universe | `cp.tl.c2c.infer()` internals | public task behavior |
| `time_series()` / `module_series()` | temporal expression extraction | `cp.tl.c2c` helper | temporary public helper |
| `lr_activity()` / `build_pair_series()` | LR activity + receiver response profiles | `cp.tl.c2c.infer()` internals | public task behavior |
| `score_lr_pairs()` | early-window pair and ligand scoring | `cp.tl.c2c.infer()` internals | public task behavior |
| LR scatter plot | pair score visualization | `cp.pl.c2c.interaction_scatter()` | public plotting |
| LR/response overlay | manuscript callout helper | `cp.pl.c2c.pair_overlay()` for result profiles; legacy helper for custom fitted callouts | public plotting/helper |
| selected genes over time | reusable raw-expression plot | `cp.pl.c2c.plot_genes_over_time()` | public plotting |
| TNF/NFKBIA callout | specific biological example | active notebook | notebook-local |

## External Resources

External resources remain caller-provided:

- CellChatDB-style ligand-receptor table;
- AnnData object containing expression, cell type, and time metadata;
- tradeSeq grid exports where fitted trajectories are requested;
- ligand association table.

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
- sender/receiver definitions;
- response-gene module;
- detection and causal-window parameters;
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
