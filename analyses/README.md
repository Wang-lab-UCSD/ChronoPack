# ChronoPack analyses

The notebooks behind the ChronoPack manuscript analyses. Each one uses the
installed package (`import chronopack as cp`) and keeps the outputs of the
latest verified run (see `verification/`).

## Layout and run order

| Folder | Notebooks, in order |
|---|---|
| `k562-tnfa.hg38.rep/` (K562 + TNF-α) | `01_PreprocessDatasets` → `02_GetTrajectories` → `03_PrepareForTradeSeq` → `04_TradeSeqFits` (and optionally `04_TradeSeqFits.TwoLineage`, the genome-wide two-lineage fit, about 1.5 h) → `05_Visualization` (figures; its two-lineage gene curves use the stored two-lineage fit if present, otherwise it fits just the plotted genes, which takes a few minutes) |
| `thp1-k562-lps.hg38.rep/` (THP-1 + K562 + LPS) | `01_PreprocessDatasets` → `02_GetTrajectories` → `03_PrepareForTradeSeq` → `04_TradeSeqFits.K562`, `04_TradeSeqFits.THP1` → `05_Visualization` (figures) |
| `gene_clustering/` (both experiments) | `01_GeneTrajectoryClustering` (clustering and figures; optional variant `01_GeneTrajectoryClustering.pvalue_only`, responsive genes = p < 0.05 without the fold-change cut) |
| `c2c/` | `CellCommunication` |
| `grn/` | `GRNInference` |

Run both experiments through `04_TradeSeqFits*` before `gene_clustering/` and
`grn/`; `c2c/` needs the THP-1 + K562 + LPS folder through `04_TradeSeqFits*`.

## Data sources

The raw data are ChronoSeq time-tagged DGE tables from EMBL-EBI ArrayExpress.
The notebooks download the ones they need into `inputs/` when they are missing
(see below), or they can be downloaded by hand from each study's file list.

| Experiment | ArrayExpress | Files used | Used by |
|---|---|---|---|
| K562 + 10 ng/ml TNF-α (12 time tags, 10 min interval) | [E-MTAB-15894](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-15894) ([files](https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/894/E-MTAB-15894/Files/)) | `N70{1-4}_scK562_TNFA_Replication_RT_{1-4}.aligned.sorted.merged.tagged.dge.time_tags.csv.gz` | `k562-tnfa.hg38.rep/01_PreprocessDatasets` |
| 50:50 K562 + differentiated THP-1, 1 µg/ml LPS (12 time tags, 15 min interval) | [E-MTAB-15928](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-15928) ([files](https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/928/E-MTAB-15928/Files/)) | `N70{1-4}_scTHP1_K562_LPS_Split_RT_{1-4}_15min_Interval_Watchmaker_Replication.aligned.sorted.merged.tagged.dge.time_tags.csv.gz` | `thp1-k562-lps.hg38.rep/01_PreprocessDatasets` |
| Bulk K562 negative control, no stimulus (12 time tags) | [E-MTAB-15966](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-15966) ([files](https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/966/E-MTAB-15966/Files/)) | `N706_Bulk_K562_Control_600_Cells.aligned.sorted.merged.tagged.dge.time_tags.csv.gz` (summed per time tag into `k562_neg-ctrl_mtx-agg.pkl`) | both `03_PrepareForTradeSeq` (time-tag-bias filter) |

The studies also hold the aligned BAM files and per-sample time-tag tables, which
the notebooks do not use. OmniPath (TF ranking in `grn/`) comes from the dated
snapshot in `resources/` described below.

## Inputs and resources

Paths are relative to each notebook's folder, so they resolve when the
notebook runs from its own folder (as Jupyter does).

- `inputs/` (not in git): the raw data.
  - `k562-tnfa.hg38.rep/N70{1-4}_*.dge.time_tags.csv.gz`: time-tagged DGE tables
    (ArrayExpress E-MTAB-15894)
  - `thp1-k562-lps.hg38.rep/N70{1-4}_*.dge.time_tags.csv.gz`: time-tagged DGE tables
    (ArrayExpress E-MTAB-15928)

    The `01_PreprocessDatasets` notebooks download any of these DGE tables that
    are missing (`cp.io.download_arrayexpress_files`, about 30 MB, a few minutes);
    tables already in `inputs/` are used as they are.
  - `k562_neg-ctrl_mtx-agg.pkl`: K562 negative-control pseudobulk (the bulk K562
    control N706, ArrayExpress E-MTAB-15966, summed per time tag), read by both
    `03_PrepareForTradeSeq` notebooks. If it is missing, they download N706's DGE
    table (about 20 MB) into `inputs/` and rebuild it (about 10 minutes) with
    `cp.io.sum_dge_by_time_tag`.
- `resources/` (in git): small reference files.
  - `regev_lab_cell_cycle_genes.txt`: cell-cycle genes (S then G2/M)
  - `interaction_input_CellChatDB.csv`: CellChatDB human ligand-receptor table
    (`c2c/`, `grn/`)
  - `hs_tfs.txt`: human transcription factors (`grn/`)
  - `protein_coding_genes.csv`: hg38 protein-coding genes (`grn/`)
  - `omnipath_snapshot_2025-10-21/`: the OmniPath signaling and transcriptional
    tables (retrieved 2025-10-21) used to rank TFs in `grn/`, so the ranking does
    not change as OmniPath is updated (`cp.tl.grn.save_omnipath_snapshot` writes
    a new one)
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
