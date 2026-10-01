# ChronoPack analyses

The notebooks behind the ChronoPack manuscript analyses. Each one uses the
installed package (`import chronopack as cp`) and keeps the outputs of the
latest verified run (see `verification/`).

## Layout and run order

| Folder | Notebooks, in order |
|---|---|
| `k562-tnfa.hg38.rep/` (K562 + TNF-α) | `01_PreprocessDatasets` → `02_GetTrajectories` → `03_PrepareForTradeSeq` → `04_TradeSeqFits` |
| `thp1-k562-lps.hg38.rep/` (THP-1 + K562 + LPS) | `01_PreprocessDatasets` → `02_GetTrajectories` → `03_PrepareForTradeSeq` → `04_TradeSeqFits.K562`, `04_TradeSeqFits.THP1` |
| `gene_clustering/` (both experiments) | `01_DimensionalityReduction&Clustering` → `02_ClusterVisualization`, `03_ModeConservation` |
| `c2c/` | `CellCommunication` |
| `grn/` | `GRNInference` |

Run both experiments through `04_TradeSeqFits*` before `gene_clustering/`.
The K562 + TNF-α `03_PrepareForTradeSeq` also reads the K562 negative-control
pseudobulk from the THP-1 + K562 + LPS folder (see below). `c2c/` and `grn/`
currently read their inputs from fixed paths set at the top of each notebook,
not from the data folder.

`resources/` holds small reference files the notebooks read (the Regev lab
cell-cycle gene list). `archive/` holds older exploratory notebooks that are
not part of the pipeline.

## Data folder

The pipeline notebooks read and write everything under one data folder, set
by the `CHRONOPACK_DATA_DIR` environment variable (default: `analyses/data/`,
which git ignores):

```
$CHRONOPACK_DATA_DIR/
  k562-tnfa.hg38.rep/
    input/            raw time-tagged DGE tables (N701…N704 *.dge.time_tags.csv.gz)
    intmd/  output/   written by the notebooks
  thp1-k562-lps.hg38.rep/
    input/            raw time-tagged DGE tables
    intmd/k562_neg-ctrl_mtx-agg.pkl   K562 negative-control pseudobulk (provided)
    intmd/  output/   written by the notebooks
  gene_clustering/    written by gene_clustering/
```

Put the raw inputs (and the negative-control pickle) in place, then start
Jupyter with the ChronoPack environment, for example:

```bash
export CHRONOPACK_DATA_DIR=/path/to/chronopack_data
jupyter lab analyses/
```

and select the `ChronoPack` kernel.
