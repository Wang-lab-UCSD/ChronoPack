# Changelog

## 0.1.0 - Unreleased

Initial package-ready ChronoPack release candidate.

- Added package-backed preprocessing and AnnData validation conventions.
- Added WaddingtonOT-derived cell-state trajectory and fate infrastructure.
- Added Python-first gene-dynamics fitting orchestration backed by tradeSeq.
- Added gene-expression trajectory characterization and dynamic-mode clustering.
- Added R-backed plotting wrappers for gene trajectories and mode summaries.
- Added cross-context comparison of dynamic modes (top-k mode recovery; mode
  proportion, Sankey and chord plots).
- Added GRN inference outputs through `chronopack.tl.grn`.
- Added CellCommunication/C2C inference outputs through `chronopack.tl.c2c`.
- Added lightweight tests, R-resource packaging checks, and release smoke tests.
- Added a K562 + TNF-alpha demo notebook (`demo/`) covering trajectories, tag-bias
  filtering, single- and per-fate gene-dynamics fits, and figures.

Changes since the previous collaborator test (2026-10-19):

- Terminal fate distributions are computed exactly (`terminal_topk=None` by
  default in `cluster_cells_to_trajectories` / `CalculateTerminalFatesAcrossCells`)
  instead of keeping each cell's 200 largest entries at every step. Trajectory
  clusters, fate probabilities, two-lineage weights and the figures built on them
  change slightly (K562 + TNF-alpha clusters ARI 0.89 vs the truncated version).
- C2C ligand-receptor pairs are ranked by `z(logp) + z(activity_sum)`
  (`composite_score`), like the ligands, instead of by `logp` then activity.
- `cp.tl.grn.rank_regulators(omnipath_snapshot=...)` and
  `cp.tl.grn.save_omnipath_snapshot()`: TF ranking from a fixed OmniPath snapshot;
  the analyses use `analyses/resources/omnipath_snapshot_2025-10-21/`.
- `cp.io.download_arrayexpress_files()` and `cp.io.sum_dge_by_time_tag()`: the
  analysis notebooks download missing ChronoSeq inputs from ArrayExpress
  (E-MTAB-15894, E-MTAB-15928) and rebuild the bulk negative control from
  E-MTAB-15966.
- `05_Visualization` (K562 + TNF-alpha) no longer needs
  `04_TradeSeqFits.TwoLineage`: without the stored fit it fits the plotted genes.
- `two_gene_fate_scatter(xlim="auto", xticks="auto")` fits the x axis to the data.

Benchmark datasets, comparator methods, gold standards, and benchmarking
metrics for GRN and C2C remain outside ChronoPack for this release.
