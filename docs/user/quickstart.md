# Quickstart

This sketch shows the intended ChronoPack API shape. It is not a full
manuscript reproduction.

```python
import chronopack as cp

adata = cp.io.read_time_tagged_dge(
    sample_paths=["sample_a_dge.csv", "sample_b_dge.csv"],
    sample_names=["sample_a", "sample_b"],
    time_col="TIME_TAG",
)

cp.io.initialize_expression_layers(adata, counts_layer="counts")
cp.io.annotate_qc_metrics(adata)
cp.io.filter_technical_genes(adata)
cp.io.validate_preprocessed_adata(adata)
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
grn = cp.tl.grn.infer(
    adata,
    regulator_scores=tf_scores,
    target_scores=target_scores,
)
grn.write("grn_edges")

c2c = cp.tl.c2c.infer(
    adata,
    lr_resource=ligand_receptor_table,
    association=association_table,
    sender="THP1",
    receiver="K562",
    response_genes=["NFKBIA", "TNFAIP3"],
    cell_type_key="cell_type",
    time_key="time",
)
c2c.write("c2c_results")
```

Legacy modules such as `chronopack.tradeseq` and `chronopack.c2c` remain for
compatibility, but new workflows should use `chronopack.io`, `chronopack.tl`,
and `chronopack.pl`.
