# ChronoPack

ChronoPack is a Python-first package for time-resolved single-cell analysis
workflows with selected R statistical backends retained behind package APIs.
It standardizes the analysis path used by the Chrono-seq notebooks around
AnnData, result objects, and reusable plotting wrappers.

## Public API Shape

- `chronopack.io`: input loading, AnnData validation, and preprocessing
  conventions.
- `chronopack.tl`: cell-state trajectories, fate outputs, gene dynamics,
  gene-trajectory characterization, dynamic-mode clustering, and cross-context
  comparison of dynamic modes.
- `chronopack.tl.grn`: GRN inference and regulator/target scoring.
- `chronopack.tl.c2c`: cell-cell communication inference.
- `chronopack.pl`: reusable plotting wrappers, including R-backed plotting
  rendered to SVG/PDF/PNG.
- `chronopack.pl.c2c`: reusable C2C plotting helpers.

Benchmarking workflows for GRN and C2C remain separate research analyses and
are intentionally outside the ChronoPack package.

## Installation

ChronoPack has mixed Python/R dependencies. The recommended full-stack install
for the first release is the verified Conda environment:

```bash
conda env create -f env/environment.yml
conda activate chronopack
pip install -e .
pytest
```

For a non-editable package install after building a distribution:

```bash
pip install dist/chronopack-*.whl
```

The Python package does not assume that `pip` can provision the R/Bioconductor
backend. The full workflow requires the supported R stack, including
`tradeSeq`, `SingleCellExperiment`, `fdrtool`, `ggplot2`, and the packages
declared in the environment/developer documentation.

Cluster-specific packed-environment deployment is documented separately in
`env/README.md`.

## Quickstart

See `docs/user/quickstart.md` for a compact walkthrough of the intended API.

## Demo

`demo/K562_TNFa_demo.ipynb` runs the core workflow on one time course (K562
cells after TNF-α): optimal-transport cell-state trajectories, filtering of
time-tag-biased genes, and tradeSeq gene-dynamics fits for all cells and per
fate, with figures. It takes about 20 minutes on 4 cores; see `demo/README.md`
for its input files.

## Analyses

`analyses/` holds the manuscript analysis notebooks (preprocessing through
gene-trajectory clustering, cell-cell communication and GRN inference), each
using the installed package and keeping the outputs of the latest verified
run. `analyses/README.md` describes the data folder and run order.
`verification/` holds the scripts that re-run those notebooks stage by stage
on the cluster and compare the results with a reference run.

## License

ChronoPack is currently distributed under the repository license in `LICENSE`.
