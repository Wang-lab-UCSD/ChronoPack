# GRN Inference Migration Mapping

Status: package-backed GRN inference slice complete.

This pass migrated reusable ChronoPack GRN inference into `cp.tl.grn`. It did
not migrate GRN benchmarking, comparator methods, gold standards, evaluation
metrics, Perturb-seq/TRRUST/ENCODE benchmark logic, or C2C benchmarking.

## Public API

All three functions work on one context per call (2026-09-30); a workflow
runs them once for K562 + TNF-α and once for THP1 + LPS.

- `cp.tl.grn.rank_regulators(source=..., is_ligand=..., human_tf_path=..., ligand_receptor_path=...)`
  - Loads OmniPath and CellChat-style ligand-receptor resources. OmniPath is
    downloaded through the `omnipath` client (which may return its local
    cache) unless `omnipath_snapshot=` names a folder written by
    `cp.tl.grn.save_omnipath_snapshot()`; scores from a live download change
    slightly as OmniPath is updated (2026-10-19: RELA 0.727104 with the
    2025-10-21 tables vs 0.727092 with that day's download; same top 10 TFs).
    The analyses use `analyses/resources/omnipath_snapshot_2025-10-21/`, which
    reproduces the r9 TF rankings to 1e-16
    (`verification_runs/grn_omnipath_version_20261019/`).
  - Builds the signed OmniPath TF graph.
  - Ranks TFs downstream of one signal with the existing signed absorbing
    diffusion workflow.
  - Returns the TF ranking DataFrame (`TF`, ..., `priority_act`).

- `cp.tl.grn.score_targets(association, valid_targets=...)`
  - Runs the packaged R/fdrtool target-gene LFDR recalibration helper on one
    tradeSeq association table.
  - Returns `TargetGeneScoreResult` with `scores`, `summary`, and `params`
    DataFrames (no `dataset` column).

- `cp.tl.grn.infer(targets, tfs, valid_targets=..., name=...)`
  - Builds ChronoPack TF-target edge scores from recalibrated target-gene
    scores and TF priority scores; `name` fills the `method` column.
  - Returns `GRNResult`.
  - Uses the existing core calculation:
    `gene_score = 1 - lfdr_bum`,
    `tf_score = priority_act / sum(priority_act)`,
    `edge score = outer(tf_score, gene_score)`, followed by the historical
    deterministic tie-breaking/ranking behavior.

Compatibility:

- `cp.tl.grn.build_chronopack_ranked_edges(...)`
- `cp.tl.grn.write_chronopack_ranked_edges(...)`

(The top-level `chronopack.grn_inference` / `chronopack.omnipath_tf` modules
were removed on 2026-09-30.)

## Result Objects

`GRNResult` contains:

- `edges`: standardized edge DataFrame with
  `method`, `regulator`, `target`, `score`, `rank`, `percentile`;
- `regulators`: regulator priority tables used by each method/context;
- `target_scores`: target-gene scores used by each method/context;
- `parameters`;
- `provenance`.

The full TF-target network is not stored in AnnData because it is not aligned
to either AnnData axis. `GRNResult.write()` exports CSV and, when available,
Parquet plus JSON metadata. `GRNResult.to_legacy_edges()` exports the
historical `TF` / `EdgeWeight` schema for compatibility.

`TargetGeneScoreResult` contains the fdrtool/BUM score, summary, and parameter
tables. `TFRankingResult` was removed; `rank_regulators` returns the ranking
DataFrame directly (the diffusion convergence history is not returned).

## Current Workflow Mapping

| Historical code / notebook block | Purpose | Package location | Public? |
| --- | --- | --- | --- |
| `load_omnipath_tf_resources()` | load OmniPath, transcriptional, LR, and TF resources | `cp.tl.grn` / legacy compatibility | public |
| `build_omnipath_tf_graph()` | construct signed TF graph | `cp.tl.grn` / legacy compatibility | public |
| `rank_omnipath_tfs()` / diffusion cells | context-specific TF prioritization | `cp.tl.grn.rank_regulators()` | public task API |
| `calculate_lfdr_recalibration()` R helper | fdrtool/BUM target-gene scores | `cp.tl.grn.score_targets()` using packaged R resource | public task API |
| `load_chronopack_beta_lfdr()` | TF-target outer-product scoring | `cp.tl.grn.infer()` (`tl/grn/_inference.py`) | public task behavior |
| `build_chronopack_ranked_edges()` | historical ranked edge CSV | compatibility wrapper | temporary compatibility |
| TRRUST/ENCODE reference TF filter | benchmark-era TF-universe restriction | compatibility wrapper only | not public GRN inference API |
| benchmark comparisons/evaluation | method benchmarking | outside ChronoPack | out of scope |

## External Resources

The package does not bundle large scientific resources. Active notebooks pass
paths for:

- human TF list;
- CellChat ligand-receptor table;
- protein-coding target-gene annotation;
- upstream tradeSeq association tables;
- generated OmniPath TF-priority tables;
- generated target-gene score table.

OmniPath, CellChat-style ligand-receptor resources, Taiji/CIS-BP/contact priors
from older research workflows remain external resources unless explicitly
packaged later. The analyses keep a dated OmniPath snapshot (signaling and
transcriptional tables, about 4 MB) in `analyses/resources/`.

## Benchmark Boundary

Benchmarking is intentionally maintained outside the ChronoPack package.
ChronoPack provides inference outputs and standardized result representations;
benchmark datasets, gold standards, comparator methods, TF universes, and
evaluation metrics remain separate research workflows.

This applies to both GRN benchmarking and CellCommunication/C2C benchmarking.
No `cp.tl.grn.benchmark` or `cp.tl.c2c.benchmark` namespace was added.

## Active Notebooks

`grn_inference/GRNInference.ipynb` (the three earlier notebooks are in
`grn_inference/old/`) runs `rank_regulators`, `score_targets` and `infer` once
per context and writes `outputs/<context>/tf_ranking.csv`,
`target_scores.{scores,summary,params}.csv` and `grn.*`. It intentionally does
not apply benchmark-specific TF-universe filtering.

## Regression Verification

Fixed upstream inputs were used; tradeSeq fits, gene clustering, mode
conservation, and benchmark workflows were not rerun.

Target-gene scoring through `cp.tl.grn.score_targets()` reproduced the existing
local reference outputs exactly:

- scores: `23718 x 13`, maximum numeric absolute difference `0`;
- summary: `16 x 11`, maximum numeric absolute difference `0`;
- params: `30 x 4`, maximum numeric absolute difference `0`.

Historical compatibility edge regression used the old explicit TRRUST/ENCODE
TF-universe filter only to compare against the pre-existing benchmark-era edge
artifact:

- shape: `4334306 x 6`;
- no left-only or right-only TF-target rows after alignment by
  `method + TF + target`;
- maximum `EdgeWeight` absolute difference: `4.656612873077393e-10`;
- maximum rank absolute difference: `0`;
- maximum percentile absolute difference: `1.1102230246251565e-16`.

The clean public `cp.tl.grn.infer()` output, without benchmark-era TF filtering,
contains:

- shape: `5009770 x 6`;
- regulators: `263`;
- targets: `13469`.

Regression artifacts were written under:

- `verification_runs/grn_package_migration_20260918/`

## Sidecar Status

Generated sidecars/intermediates still useful for active workflows:

- OmniPath TF-rank CSVs;
- target-gene recalibration CSVs;
- clean GRN edge CSV/Parquet/metadata outputs;
- historical compatibility edge CSV only for regression/legacy comparison.

The migrated inference workflow does not require native backend objects in
AnnData and does not write full networks to `adata.uns`.

## Tests

Added `tests/test_grn_inference.py`.

Synthetic tests cover:

- TF-target alignment;
- invalid target removal;
- self-edge removal;
- explicit optional TF-universe filtering;
- legacy edge-schema export;
- CSV round-trip;
- compatibility wrapper agreement on synthetic files;
- packaged R resource discovery and `prepare_p()` smoke behavior.

Together with the existing OmniPath diffusion tests, the focused GRN test slice
passes: `7 passed`.
