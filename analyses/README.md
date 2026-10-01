# ChronoPack analyses

The notebooks behind the ChronoPack manuscript analyses. Each one uses the
installed package (`import chronopack as cp`) and keeps the outputs of the
latest verified run (see `verification/`).

## Layout and run order

| Folder | Notebooks, in order |
|---|---|
| `k562-tnfa.hg38.rep/` (K562 + TNF-α) | `01_PreprocessDatasets` → `02_GetTrajectories` → `03_PrepareForTradeSeq` → `04_TradeSeqFits` (and optionally `04_TradeSeqFits.TwoLineage`, not used downstream) → `05_Visualization` (figures) |
| `thp1-k562-lps.hg38.rep/` (THP-1 + K562 + LPS) | `01_PreprocessDatasets` → `02_GetTrajectories` → `03_PrepareForTradeSeq` → `04_TradeSeqFits.K562`, `04_TradeSeqFits.THP1` → `05_Visualization` (figures) |
| `gene_clustering/` (both experiments) | `01_GeneTrajectoryClustering` (clustering and figures) |
| `c2c/` | `CellCommunication` |
| `grn/` | `GRNInference` |

Run both experiments through `04_TradeSeqFits*` before `gene_clustering/` and
`grn/`; `c2c/` needs the THP-1 + K562 + LPS folder through `04_TradeSeqFits*`.

## Inputs and resources

Paths are relative to each notebook's folder, so they resolve when the
notebook runs from its own folder (as Jupyter does).

- `inputs/` (not in git): the raw data.
  - `k562-tnfa.hg38.rep/N70{1-4}_*.dge.time_tags.csv.gz`: time-tagged DGE tables
  - `thp1-k562-lps.hg38.rep/N70{1-4}_*.dge.time_tags.csv.gz`: time-tagged DGE tables
  - `k562_neg-ctrl_mtx-agg.pkl`: pre-computed K562 negative-control pseudobulk,
    read by both `03_PrepareForTradeSeq` notebooks
- `resources/` (in git): small reference files.
  - `regev_lab_cell_cycle_genes.txt`: cell-cycle genes (S then G2/M)
  - `interaction_input_CellChatDB.csv`: CellChatDB human ligand-receptor table
    (`c2c/`, `grn/`)
  - `hs_tfs.txt`: human transcription factors (`grn/`)
  - `protein_coding_genes.csv`: hg38 protein-coding genes (`grn/`)
- `archive/`: older exploratory notebooks that are not part of the pipeline.

## Data folder

The pipeline notebooks write their intermediate and final results under one
data folder, set by the `CHRONOPACK_DATA_DIR` environment variable (default:
`analyses/data/`, which git ignores). The `01_PreprocessDatasets` notebooks
create it:

```
$CHRONOPACK_DATA_DIR/
  k562-tnfa.hg38.rep/{intmd,output}/
  thp1-k562-lps.hg38.rep/{intmd,output}/
  gene_clustering/  c2c/  grn/
```

With `inputs/` in place, start Jupyter with the ChronoPack environment, for
example:

```bash
export CHRONOPACK_DATA_DIR=/path/to/chronopack_data
jupyter lab analyses/
```

and select the `ChronoPack` kernel.
