# ChronoPack Repository Reconnaissance

Date: 2026-09-16

This note records my current understanding of this workspace and the refactor target. It is intended as durable context while the Codex database is unavailable.

## Workspace layout

The root currently contains one architecture note and three experiment-specific notebook directories:

- `thp1-k562-lps.hg38.rep/`
  - Original path: `~/projects/chrono_seq/analyses/thp1-k562-lps.hg38.rep/notebooks`
  - Main notebooks: `PreprocessDatasets.ipynb`, `GetTrajectories.ipynb`, `PrepareForTradeSeq.ipynb`, `TradeSeqFits.K562.ipynb`, `TradeSeqFits.THP1.ipynb`, `TradeSeqFits.Both.ipynb`, `LocalFDR.ipynb`
- `k562-tnfa.hg38.rep/`
  - Original path: `~/projects/chrono_seq/analyses/k562-tnfa.hg38.rep/notebooks`
  - Main notebooks: `PreprocessDatasets.ipynb`, `GetTrajectories.ipynb`, `PrepareForTradeSeq.ipynb`, `TradeSeqFits.ipynb`, `OmnipathTFPrioritization.ipynb`, `GetGeneLengths.ipynb`, `VolcanoPlot.v2.ipynb`
- `k562-tnfa-puromycin.hg38/`
  - Original path: `~/projects/chrono_seq/analyses/k562-tnfa-puromycin.hg38/notebooks/`
  - Main notebooks: `Preprocess_GetTrajectories_MapFates.ipynb`, `PrepareForTradeSeq.ipynb`, `TradeSeqFits.ipynb`, `DimensionalityReductionClustering.ipynb`, older dimensionality-reduction variants, `VolcanoPlot.ipynb`

The notebooks reference data outside this workspace through relative paths from their original analysis directories. The `OriginalDirectoryPath.txt` files are the authority for resolving those paths. Per instruction, those external/original files should be treated as read-only.

There is also `mixed_R_python_package_architecture.md`, which already sketches a likely monorepo design: a formal R package plus a formal Python package, with one language serving as the primary public API and coarse-grained cross-language boundaries.

## High-level workflow

Across the three experiments, the repeated workflow appears to be:

1. Build AnnData objects from count matrices and time-tag metadata.
2. Run standard scRNA-seq QC and preprocessing with Scanpy:
   - QC metrics
   - cell/gene filtering
   - doublet/scrublet checks
   - normalization/log transform
   - highly variable genes
   - cell-cycle scoring/regression
   - PCA/neighbors/UMAP/Leiden
   - sometimes ComBat or metacell construction
3. Construct Waddington-OT transport maps and transition kernels.
4. Derive trajectory clusters and, when applicable, terminal fate probabilities/retention summaries.
5. Prepare matrices and metadata for tradeSeq:
   - counts genes x cells
   - cell times
   - full-library size factors
   - subset gene lists
   - sometimes cell-type or lineage-specific outputs
6. Fit tradeSeq NB-GAM models and run association/dynamic tests.
7. Downstream analyses:
   - dynamic mode / trajectory-shape clustering
   - volcano plots and local FDR summaries
   - fate/mass-retention summaries
   - GRN/TF prioritization using OmniPath/networkx in `k562-tnfa.hg38.rep`
   - gene length extraction from GTF annotations

## Language split

The Python-heavy portions are:

- AnnData construction and Scanpy preprocessing
- Waddington-OT transport map handling
- sparse transition-kernel utilities
- fate probability and terminal-cluster calculations
- tradeSeq input export from AnnData
- pseudobulk, smoothing, and cross-condition correlation calculations
- OmniPath/networkx GRN prioritization

The R-heavy portions are:

- tradeSeq fitting and hypothesis testing
- SingleCellExperiment construction
- tradeSeq plotting/export helpers
- dynamic trajectory/gene-mode clustering
- volcano plots
- local FDR summaries
- gene length summaries from GTF via Bioconductor tooling

This supports a mixed package strategy. A practical default is:

- Python package: data loading/preprocessing/WOT/fates/export utilities.
- R package: tradeSeq wrappers, plotting, trajectory feature clustering, gene mode utilities.
- Cross-language boundary: notebooks call each package directly, rather than tightly interleaving R -> Python -> R internals.

## Repeated code worth modularizing

Strong Python candidates:

- `map_time_to_group`
- `read_tmap`
- `build_transition_matrices`
- `_topk_csr_rows`, `topk_csr_rows`, `keep_topk`
- `CalculateTerminalFatesAcrossCells`
- `cosine_topk_leiden_from_terminal_dist`
- `terminal_fates_from_terminal_clusters`
- `fate_umap_2x2_scanpy`
- `mass_retention_plots_and_umap`
- `PlotHardClusterMassProportionsOverTime`
- `TemporalPseudobulkFromAnnData`
- `WriteTradeSeqInputsFromAnnData`
- smoothing/correlation helpers:
  - `_pearsonr_nan`
  - `_sorted_time_union`
  - `_to_float_time`
  - `smooth_series`
  - `ComputeSmoothCorrelations`
  - `variance_fraction`
  - `positive_r2`
  - `explained_fraction`

Strong R candidates:

- `plot_tradeSeq_gene_pub`
- `plot_tradeSeq_gene_pub_v2`
- `plot_tradeSeq_panel`
- `mk_lineage_vec`
- `get_pt_vec`
- `to_minutes`
- tradeSeq fit export helpers such as observed-time and grid exporters
- `plot_association_volcano`
- dynamic trajectory/gene mode helpers:
  - `z_row`
  - `scale_rows`
  - `get_descriptors`
  - onset/shape metric helpers
  - linear/logistic/impulse model helpers
  - cluster medoid/top-gene plotting helpers
  - cluster reassignment/sankey helpers

Experiment-specific code that may stay in notebooks or thin config files:

- exact input file names
- sample/replicate names and barcode conventions
- cell-type-specific cluster renaming
- ad hoc excluded-gene lists
- manually selected marker genes/query genes
- exploratory plotting cells
- one-off checks and old/legacy branches

## Candidate ChronoPack module boundaries

Python package sketch:

- `chronopack.io`
  - count/time-tag loading
  - AnnData construction
  - path/config helpers
- `chronopack.preprocessing`
  - standard Scanpy preprocessing pipeline
  - optional ComBat/batch-correction wrapper
  - metacell construction
- `chronopack.wot`
  - transport-map reading
  - transition matrix construction
  - sparse top-k pruning/normalization
- `chronopack.fates`
  - terminal fate probabilities
  - terminal cluster labeling
  - mass-retention summaries
- `chronopack.tradeseq_export`
  - AnnData-to-tradeSeq CSV export
  - library-size and subset-gene handling
- `chronopack.pseudobulk`
  - temporal pseudobulk CPM/logCPM matrices
  - smooth correlations
- `chronopack.grn`
  - OmniPath/TF prioritization helpers

R package sketch:

- `R/tradeseq_fit.R`
  - SingleCellExperiment construction
  - `fitGAM` wrappers
  - association/pattern testing wrappers
- `R/tradeseq_plot.R`
  - gene smoother plots
  - panel/grid export of fitted values
- `R/volcano.R`
  - association volcano plot helpers
- `R/dynamic_modes.R`
  - trajectory matrix construction
  - row scaling and descriptors
  - dynamic mode clustering and plotting
- `R/local_fdr.R`
  - local FDR summaries/plots
- `R/gene_lengths.R`
  - GTF-derived gene length summaries

## Conservative migration sequence

1. Create package scaffolding without changing notebook behavior.
2. Extract the most duplicated, lowest-risk pure helpers first:
   - sparse top-k matrix utilities
   - time parsing/grouping
   - pseudobulk and smoothing helpers
   - R plotting helper utilities that only transform existing objects
3. Add small tests with synthetic matrices/data frames before replacing notebook cells.
4. Replace notebook helper-definition cells with imports/library calls one notebook at a time.
5. Extract larger workflow functions only after their inputs/outputs are stable across notebooks:
   - Scanpy preprocessing pipeline
   - WOT transition/fate workflow
   - tradeSeq fitting/export wrappers
6. Keep experiment-specific paths and parameters in notebooks or YAML/JSON configs rather than hard-coding them into ChronoPack.
7. After each notebook conversion, re-run only the affected notebook or a reduced synthetic test where full data are too expensive.

## Cautions

- This is scientific analysis code. Small, reversible changes are preferable.
- Do not delete or mutate original data files. Use `OriginalDirectoryPath.txt` to resolve reads only.
- Several notebooks contain legacy/old variants. Treat the newest non-`Old` notebook as the likely source of truth unless the user says otherwise.
- The notebooks are mixed Python/R, so automated parsing can misclassify some cells. This note is based on structural inspection, headings, package calls, file literals, and repeated function names rather than full semantic execution.
- Full validation will require the original data and the relevant Python/R environments.
