# ChronoPack demo

`K562_TNFa_demo.ipynb` walks through the core ChronoPack workflow on the K562 +
TNF-α time course (12 time points, 10 minutes apart, 6,857 cells), starting from
preprocessed data:

1. cell-state trajectories with Waddington-OT, trajectory clusters, fate
   probabilities and soft memberships;
2. time-tag-bias and low-count gene filtering;
3. tradeSeq gene-dynamics fits for all cells and one per trajectory, on a subset
   of genes (`N_GENES`, default 300);
4. figures: OT-UMAPs, an association volcano and fitted gene curves.

It runs in about 20 minutes on 4 cores and keeps the outputs of its last
verified run, so it can be read without running it.

## Input files

Put these in `demo/data/` (or point `CHRONOPACK_DEMO_DATA` at their folder):

| File | Contents |
|---|---|
| `k562_tnfa.hvg.h5ad` | 3,000 highly variable genes; `X` = expression corrected for QC covariates, cell cycle and sample; preprocessing UMAP |
| `k562_tnfa.counts.h5ad` | all genes; `layers["counts"]` = raw UMI counts |
| `k562_neg-ctrl_mtx-agg.pkl` | negative-control pseudobulk for the time-tag check |

`make_demo_data.py` builds them from the outputs of
`analyses/k562-tnfa.hg38.rep/01_PreprocessDatasets.ipynb`:

```bash
python demo/make_demo_data.py "$CHRONOPACK_DATA_DIR/k562-tnfa.hg38.rep/output" \
    analyses/inputs/k562_neg-ctrl_mtx-agg.pkl demo/data
```

## Running

With the ChronoPack environment (including its R packages) active:

```bash
jupyter lab demo/
```

and select the `ChronoPack` kernel. Results (transport maps, tradeSeq inputs and
fits, PDF figures) are written to `demo/work/` (or `CHRONOPACK_DEMO_WORK`).
Neither `demo/data/` nor `demo/work/` is tracked by git.
