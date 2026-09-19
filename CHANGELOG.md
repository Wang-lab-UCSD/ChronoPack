# Changelog

## 0.1.0 - Unreleased

Initial package-ready ChronoPack release candidate.

- Added package-backed preprocessing and AnnData validation conventions.
- Added WaddingtonOT-derived cell-state trajectory and fate infrastructure.
- Added Python-first gene-dynamics fitting orchestration backed by tradeSeq.
- Added gene-expression trajectory characterization and dynamic-mode clustering.
- Added R-backed plotting wrappers for gene trajectories and mode summaries.
- Added cross-context mode conservation and mode-switching result objects.
- Added GRN inference outputs through `chronopack.tl.grn`.
- Added CellCommunication/C2C inference outputs through `chronopack.tl.c2c`.
- Added lightweight tests, R-resource packaging checks, and release smoke tests.

Benchmark datasets, comparator methods, gold standards, and benchmarking
metrics for GRN and C2C remain outside ChronoPack for this release.
