# GRN inference: how TFs are ranked, and known limitations

`chronopack.tl.grn` infers transcription factor (TF) → target gene edges for one
stimulus context in three steps:

1. **Rank TFs** downstream of the stimulus with `cp.tl.grn.rank_regulators`, a
   signed diffusion over the OmniPath signaling network (this page).
2. **Score target genes** with `cp.tl.grn.score_targets`: the context's tradeSeq
   association p-values become local FDRs (fdrtool / beta-uniform mixture).
3. **Infer edges** with `cp.tl.grn.infer`: each TF → target edge is weighted by
   the TF's share of the total TF priority times the target's 1 − local FDR.

The TF ranking is a prior from a generic human signaling network. It says which
TFs a signal can reach and with which sign; it does not use the cells'
expression unless asked to (see [Expression filter](#expression-filter)).

```python
tfs = cp.tl.grn.rank_regulators(
    source="TNF", is_ligand=True,                      # or e.g. source="TLR4" for a receptor
    human_tf_path="hs_tfs.txt",
    ligand_receptor_path="interaction_input_CellChatDB.csv",
    omnipath_snapshot="omnipath_snapshot_2025-10-21",  # fixed OmniPath tables
)
targets = cp.tl.grn.score_targets("asso_df.1trajectory.nknots4.csv", valid_targets="protein_coding_genes.csv")
grn = cp.tl.grn.infer(targets, tfs, valid_targets="protein_coding_genes.csv", name="ChronoPack_K562")
```

The files are in `analyses/resources/`; `analyses/grn/GRNInference.ipynb` runs
the K562 + TNF-α and THP-1 + LPS contexts.

## How the ranking works

- **Network.** Directed, signed OmniPath signaling interactions (activity flow,
  enzyme-substrate, ligand-receptor and other resources). Interactions supported
  only by resources in `DEFAULT_DROP_SOURCES` (NetPath, KEGG-MEDICUS, SPIKE,
  HPRD and a few others) are dropped, and two-gene complexes are split into their
  genes.
- **TFs.** Human TFs (`hs_tfs.txt`) that have OmniPath transcriptional targets and
  are in the network: 501 with the 2025-10-21 snapshot.
- **Start.** With `is_ligand=False`, all signal starts at the `source` receptor.
  With `is_ligand=True`, it starts at the ligand's receptors in the
  ligand-receptor table: each receptor gets an equal share, and a receptor complex
  (`IL6R_IL6ST`, `TGFbR1_R2`) is split over its subunits.
- **Diffusion.** At each step a protein passes its signal on to its targets in
  equal shares. Signal is tracked as two pools, activating and inhibiting:
  activating edges keep the pool, inhibitory edges swap them. Each step keeps 85%
  of the signal (`alpha=0.85`), and signal that reaches a TF stays there.
- **Score.** `priority_act` is the activating minus the inhibiting signal a TF
  receives, set to 0 if negative. `infer` uses each TF's share of the total.
- **Speed.** The network is built once per session and set of options (about
  20 s); each later call takes about 0.1 s.

## Defaults and options

| Parameter | Default | What it does |
|---|---|---|
| `conflicting_signs` | `"consensus"` | Edges that resources annotate as both activating and inhibiting (about 1,300) take OmniPath's consensus sign. `"stimulation"` counts them as activating, the behaviour before 2026-10-21. |
| `tied_signs` | `"split"` | Conflicting edges whose consensus is tied (about 200) pass signal of unknown sign: half keeps its sign, half flips. `"stimulation"` counts them as activating. Example: IκBα → RELA/NFKB1 is tied; read as activating, it turned IKK's (correct) inhibition of IκBα into inhibition of NF-κB. |
| `second_messenger_bridges` | `True` | Adds the steps a protein network cannot hold: adenylyl cyclases → PKA/EPAC (cAMP), phospholipase C → calcineurin/CaMKII/IV (IP3/Ca²⁺) and → PKC/RasGRP (DAG), where OmniPath has no edge (141 edges). |
| `expressed_genes`, `expression_scope` | `None`, `"tfs"` | Optional expression filter, see below. |
| `edge_weight` | `"uniform"` | `"resources"` / `"references"` split a protein's signal by how much evidence supports each edge. It helps kinase cascades, but over-weights well-studied inhibitory IκB edges and pushes NFKB1 out of the top TFs for TNF and TLR4, so it is off by default. |
| `evidence` | `"all"` | `"physical"` / `"physical_curated"` keep only edges with protein-protein, enzyme-substrate or ligand-receptor support. Only about a fifth of OmniPath's signaling edges (about 45% of edges into TFs) have such support; the rest come from pathway resources that record that A acts on B, directly or not. Restricting to them ranks canonical TFs worse. |
| `omnipath_snapshot` | `None` | Without a snapshot, OmniPath is downloaded (or read from the `omnipath` client's cache), so rankings change slightly as it is updated. Pass a folder written by `cp.tl.grn.save_omnipath_snapshot` for fixed tables. |

`rank_regulators` raises `ValueError` when no signal reaches any TF (the receptor
has no outgoing OmniPath edges), when a ligand has no receptor in the network, or
when the source or every receptor is removed by the expression filter.

## Expression filter

```python
k562 = adata[adata.obs["cell_type"] == "K562"]
tfs = cp.tl.grn.rank_regulators(
    source="TNF", is_ligand=True, ...,
    expressed_genes=cp.tl.grn.expressed_genes(k562),  # detected in >= 1 cell
)
```

- `expression_scope="tfs"` (default when filtering) removes unexpressed TFs and
  ligand receptors, the ends of the network.
- `"network"` also removes unexpressed signaling proteins in between. That risks
  dropping genuine pathway members: single-cell RNA-seq rarely detects them (in
  the K562 and THP-1 data, IKBKG is detected in 0.3% of cells, MYD88 in 2%,
  TNFRSF1A in 4%, RELA in 4–7%).
- `expressed_genes` counts a gene as expressed if it is detected in at least one
  cell. Stricter percentage cutoffs remove those same pathway members (at 5% of
  cells, TNFRSF1A and RELA drop out).
- The "detected in one cell" rule depends on the number of cells. With a few
  hundred cells, a rarely detected gene can be missed by chance; with many
  thousands, ambient RNA makes almost every gene pass.
- The filter is off by default. For K562 (TNF) and THP-1 (TLR4) it removes
  unexpressed TFs such as ESR1, IRF5 and SMAD1 from the top of the ranking.

## How well it works

The ranking was checked on 45 stimuli with textbook downstream TFs: NF-κB
pathways, interferons, gp130 and γc/βc cytokines, IL-4/13, TGF-β/BMP, Wnt,
Hedgehog, Notch, receptor tyrosine kinases, insulin, cAMP-coupled receptors and
the TCR. The figures below are for the current defaults; this check is not part
of the package.

| Seeded at | Stimuli ranked | Canonical TF in top 10 / top 3 | Median best rank | Mean average precision |
|---|---|---|---|---|
| Receptor | 42 of 45 | 39 / 32 | 2 | 0.43 |
| Ligand (when in the ligand-receptor table) | 33 of 35 | 30 / 24 | 2 | 0.41 |
| Stimulus-blind ranking (how often each TF receives signal from any receptor) | — | about 29 / 19 | 4 | about 0.25 |

- **Ranked well (canonical TF in the top 1–3):**
  - NF-κB for TNF, IL-1, CD40, LTβR and the TLRs;
  - STAT1/STAT2 for the interferons, STAT3 for the gp130 cytokines and STAT5 for
    the γc/βc cytokines;
  - SMADs for TGF-β/BMP, GLI for Hedgehog, RBPJ for Notch;
  - CREB1 for the β2-adrenergic receptor (rank 3; rank 9 for the PTH receptor).
- **How specific the top of the list is.** On average about 2 of a stimulus's top
  10 TFs are its textbook TFs. The rest are mostly network hubs or TFs that are
  plausible but not on the list.

## Known limitations

- **Long kinase cascades.** Signal is split at every step, and hub kinases
  (MAPK1, AKT1, PKA) pass it to hundreds of substrates. TFs reached only through
  RAS-RAF-MEK-ERK or PI3K-AKT rank low: ELK1 for growth-factor receptors (rank 80
  or lower), and FGF2 / NGF / insulin, whose best canonical TF ranks about 16, 24
  and 58. TFs with a direct edge from the receptor (STAT3 for EGFR) get the
  signal first.
- **Which STAT.** Every JAK links to every STAT in the network, so STAT1 and
  STAT5A lead for almost every cytokine. IL-4/IL-13 → STAT6 reaches only
  rank 7–9.
- **Hub TFs.** STAT3, STAT1, STAT5A, JUN and SMAD3 rank high for many unrelated
  signals. For the STAT3/STAT5 cytokines and EGF, the stimulus-blind ranking does
  about as well as the method.
- **Repressed TFs.** The score rewards net activation only, so TFs that a
  stimulus inhibits are not prioritized, even when their sign is right (insulin
  gives FOXO1 net inhibition, ranked low). Inhibition reached over several
  steps is less reliable: every edge sign along the path has to be right.
- **Network coverage.** Receptors without outgoing OmniPath edges cannot be
  ranked: IL17RA/IL17RC (IL-17), TNFRSF11A (RANKL) and IL12RB2 (IL-12) in the
  2025-10-21 snapshot. Signal-dependent effectors that are not in the human TF
  list (SMAD2, β-catenin) cannot be ranked.
- **Sign quality.** Edge signs come from OmniPath's resources and their
  consensus. Where these are wrong or mixed (feedback, different contexts), the
  net score can be off. Tied signs pass signal without a direction rather than
  guessing.
- **No expression or time information.** The ranking does not look at whether the
  TFs' targets actually change in the data. The optional expression filter only
  removes genes that are absent.
