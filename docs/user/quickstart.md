# Quickstart

This sketch shows the intended ChronoPack API shape. It is not a full
manuscript reproduction.

```python
import chronopack as cp

import scanpy as sc

adata = cp.io.read_time_tagged_dge(
    sample_paths=["sample_a_dge.csv", "sample_b_dge.csv"],
    sample_names=["sample_a", "sample_b"],
    time_col="TIME_TAG",
)

# Standard Scanpy QC and normalization; ChronoPack expects counts / cpm / log1p layers.
adata.var["mt"] = adata.var_names.str.startswith(("MT-", "mt-"))
adata.var["rb"] = adata.var_names.str.startswith(("RPS", "RPL"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "rb"], percent_top=None, log1p=False, inplace=True)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
adata.layers["cpm"] = adata.X.copy()
sc.pp.log1p(adata)
adata.layers["log1p"] = adata.X.copy()
```

Cell-state trajectory and fate workflows mutate AnnData where outputs are
cell-aligned and return a result object for non-axis-aligned WOT state:

```python
traj = cp.tl.get_trajectories(
    adata,
    tmap_dir="tmaps",
    key_added="cell_state_trajectories",
)

cp.tl.compute_terminal_fate_probabilities(
    adata,
    terminal_cluster_key="traj_clusters",
    out_key="terminal_fate_probs",
)
```

Gene dynamics and dynamic-mode analysis write standardized gene-level outputs
to `adata.var`, `adata.varm`, and `adata.uns["chronopack"]`:

```python
fit = cp.tl.fit_gene_dynamics(
    adata,
    pseudotime_key="pseudotime",
    cell_weights_key="cell_weights",
    counts_layer="counts",
    key_added="gene_dynamics",
)

char = cp.tl.characterize_gene_trajectories(
    adata,
    output_dir="gene_trajectories",
    key_added="gene_trajectories",
)

clusters = cp.tl.cluster_gene_trajectories(
    adata,
    output_dir="gene_trajectories",
    key_added="gene_trajectory_clusters",
)
```

Reusable plots consume standardized AnnData state and render R-backed figures
directly to files:

```python
plot = cp.pl.gene_trajectory(
    adata,
    genes=["JUN", "FOS"],
    key="gene_trajectories",
    output="jun_fos.svg",
)
```

Cross-context analyses return result objects rather than forcing paired-dataset
outputs into one AnnData:

```python
mode_res = cp.tl.mode_conservation(
    source=adata_tnf,
    target=adata_lps,
    source_label="K562_TNFa",
    target_label="THP1_LPS",
)

cp.pl.mode_sankey(mode_res, output="mode_switches.svg")
```

GRN and C2C workflows also return transparent result objects:

```python
tfs = cp.tl.grn.rank_regulators(
    source="TNF",
    is_ligand=True,
    human_tf_path="hs_tfs.txt",
    ligand_receptor_path="interaction_input_CellChatDB.csv",
)
targets = cp.tl.grn.score_targets("asso_df.1trajectory.lfdr.csv", valid_targets="protein_coding_genes.csv")
grn = cp.tl.grn.infer(targets, tfs, valid_targets="protein_coding_genes.csv", name="ChronoPack_K562")
grn.write("grn_edges")

c2c = cp.tl.c2c.infer(
    adata=adata,
    lr_resource=ligand_receptor_table,
    association=association_table,
    sender="THP1",
    receiver="K562",
    causal_window=(0, 45),  # minutes; LR activity is summed over this window
)
c2c.write("c2c_results")
fig, ax = cp.pl.c2c.interaction_scatter(c2c, x_thresh=0.4)
```

The complete analysis notebooks, with outputs from the latest verified run,
are in `analyses/` (see `analyses/README.md`).
