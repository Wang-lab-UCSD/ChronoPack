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
- Scaling to larger datasets (e.g. 6 time points x 8,000 cells):
  - `cp.io.read_time_tagged_dge` returns a sparse CSR `X` (float64, so values are
    unchanged; `sparse=False` for dense) and reads each table in chunks; the
    notebooks' `counts`/`cpm`/`log1p` layers are therefore sparse too. Results
    match the dense version to rounding (log1p differs in the last bit); the
    K562 + TNF-alpha UMAP and OT-UMAP layouts shift slightly as a result, while
    trajectory clusters, fates, tradeSeq fits and gene clusters are unchanged.
  - The WOT transition graph is built map by map with vectorized pruning instead
    of a long table, Python lists and LIL matrices. `T_full` and
    `transitions_long` on the `get_trajectories` result are built only when
    accessed. New `transition_mass=` (keep each cell's largest transitions up to
    that fraction of its mass; combine with `transition_topk=None` to use it alone).
  - The terminal-fate cosine graph is computed with dense float32 products
    instead of sparse ones. On K562 + TNF-alpha, 32 of 3.4M edges differ (ties at
    the 500th neighbour) and the trajectory clusters are identical.
- `cp.tl.prescreen_time_association()`: optional pre-screen of genes before the
  tradeSeq fits (least-squares F test on a cubic in time, keep p < 0.2). It kept
  32-41% of genes and 97-99% of the nominally significant ones in the three
  single-lineage fits (12 time tags or 6 merged). Inputs written from the
  pre-screened genes record how many were screened, and the fits' q-values count
  all of them. The `03_PrepareForTradeSeq` notebooks include it as a commented-out
  cell; the analyses do not use it.
- `cp.tl.ot_force_directed()`: force-directed layout of the WOT transition graph
  (`obsp["ot_graph"]`), an alternative to the OT-UMAP, stored in
  `obsm["ot_force_directed"]`. It runs `scanpy.tl.draw_graph` (Fruchterman-Reingold
  through python-igraph by default, weighted by the transition probabilities;
  `layout="fa"` for ForceAtlas2 needs `fa2-modified`). The analyses still use the
  OT-UMAP.
- Optional growth weighting of the terminal-fate calculation:
  `cluster_cells_to_trajectories(..., growth=...)` (and
  `CalculateTerminalFatesAcrossCells`) weights each transition by the descendant
  mass of the cell it reaches, using WOT's growth weights (each cell's coupling
  row sum relative to its map's mean), as WOT's `TransportMapModel.fates` does.
  With the unpruned transition matrix it reproduces WOT's fates. The default
  (`growth=None`) is unchanged: a Markov chain on the row-normalized transitions.
  `get_trajectories` now keeps the growth weights (`result.growth`, and
  `adata.obs["ot_growth"]` when it attaches the graph), and
  `cp.tl.coupling_from_kernel(T, growth)` rebuilds the couplings from a transition
  matrix. The terminal-cell distributions record `topk` and `growth` in
  `adata.uns["chronopack"]["terminal_cell_dist"]`. On the current data, unpruned
  growth-weighted fates match WOT's to 1e-7; with the analyses' pruned graph,
  growth weighting changes fate probabilities by 0.014 (K562 + TNF-alpha) and
  0.029 (co-culture K562 cells) on average and the dominant fate of 3% and 7% of
  cells, mostly near 50/50. The analyses keep the default.
- `cp.tl.grn.rank_regulators` (TF prioritization), after a check against canonical
  TFs for 45 stimuli (`verification_runs/grn_stimulus_panel_20261021/`):
  - Ligands (`is_ligand=True`) are seeded at the subunits of their CellChat
    receptor complexes (`IL6R_IL6ST`, `TGFbR1_R2`, ...): each receptor gets an equal
    share, split over its subunits. Before, only receptors named by a single gene
    were used, so 17 of 35 ligands tested (interferons, IL-6 family, IL-2 family,
    TGF-beta, BMP, Wnt, Hedgehog, ...) raised `ValueError` and IL1B was seeded only
    at its decoy receptor IL1R2. Ligands with single-gene receptors (TNF) are
    seeded as before.
  - Raises `ValueError` when no signal reaches any TF (e.g. IL17RA, TNFRSF11A and
    IL12RB2 have no outgoing edges in OmniPath) instead of returning every TF with
    score 0 in alphabetical order; warns when no TF gets net activating signal.
  - Edges that OmniPath flags both stimulatory and inhibitory (1,315 graph edges)
    take OmniPath's consensus sign (`conflicting_signs="consensus"`, default;
    `"stimulation"` restores the previous behaviour, which counted them all as
    activating). E.g. AKT1 -> FOXO1 is now inhibitory. Across the 45 stimuli the
    rankings are about as good as before (CREB1 for PTH rises from rank 246 to 32);
    the TNF (K562) and TLR4 (THP1) rankings keep the same top-10 TFs in a slightly
    different order (Spearman 0.95 and 0.85 over all TFs), so the stored GRN outputs
    change slightly when `GRNInference` is rerun.
  - Edges whose OmniPath consensus is tied (197 graph edges; resources disagree
    equally) carry signal of unknown sign (`tied_signs="split"`, default): half
    the mass keeps its sign and half flips, so they pass signal on without a net
    direction. Before, they counted as activating; e.g. IkBa -> RELA / NFKB1 is
    tied, and read as activating it turned IKK's inhibition of IkBa (consensus
    sign) into inhibition of NF-kB, dropping NFKB1 below STAT3 for TNF. With the
    split, NFKB1 is back above STAT3 and the 45-stimulus panel is unchanged or
    slightly better. `"stimulation"` restores the previous treatment.
  - Second-messenger steps that OmniPath's protein network lacks are added as
    activating edges (`second_messenger_bridges=True`, default): adenylyl cyclases ->
    PKA/EPAC (cAMP), phospholipase C -> calcineurin/CaMKII/IV (IP3/Ca2+) and -> PKC
    isoforms/RasGRP (DAG), 141 edges where OmniPath has none. CREB1 for PTH1R and
    ADRB2 rises from rank 32 / 29 to 10 / 3; the TNF and TLR4 rankings barely change
    (Spearman 0.99, same top 10). `False` gives the previous graph.
  - New options, off by default because they did not meet the criteria set before
    testing: `edge_weight="resources"` / `"references"` (split each node's mass by
    the evidence for each edge; helps kinase cascades but moves NFKB1 out of the top
    TFs for TNF and TLR4) and `evidence="physical"` / `"physical_curated"` (only edges
    with protein-protein, enzyme-substrate or ligand-receptor support, 18% of
    OmniPath's signaling edges and 44% of those into TFs; worse rankings overall).
  - The docstring lists the remaining weak points: TFs reached through long kinase
    cascades, the choice among STATs, and hub TFs.
  - `expressed_genes=` restricts the network to the genes expressed in the
    responding cell type (receptors, signaling proteins and TFs that are not
    expressed are removed before the diffusion); `cp.tl.grn.expressed_genes(adata)`
    lists the genes detected in at least one cell. Off by default (whole human
    network). `expression_scope="tfs"` (default when filtering) removes only
    unexpressed TFs and ligand receptors; `"network"` also removes unexpressed
    signaling proteins. For K562 (TNF) and THP1 (TLR4) the network-wide filter removes ESR1 (and IRF5, SMAD1 in
    K562) from the top 20 and keeps RELA/NFKB1 at the top (Spearman 0.91-0.93 with
    the unfiltered ranking); stricter cutoffs (>= 1-5% of cells) remove genuine
    pathway members such as IKBKG, TNFRSF1A and RELA.
  - Faster: the diffusion uses sparse matrix products (same results) and the
    OmniPath graph is built once per session and set of options (first call ~20 s,
    later calls ~0.1 s; previously ~25 s every call).
  - `docs/user/grn.md` describes how TFs are ranked, the defaults and options, how
    the ranking was checked on 45 stimuli, and its known limitations; the README,
    the quickstart and `analyses/grn/GRNInference.ipynb` link to it.
- The `wot` extra now lists `igraph` and `leidenalg`, which the Leiden clustering
  of trajectories already needed (the conda environment file already had them).

Benchmark datasets, comparator methods, gold standards, and benchmarking
metrics for GRN and C2C remain outside ChronoPack for this release.
