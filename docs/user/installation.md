# Installation

ChronoPack combines Python workflow code with R/Bioconductor statistical
backends. For the first release, the supported full-stack installation path is
the Conda environment in `env/environment.yml`.

## Standard Full-Stack Install

```bash
conda env create -f env/environment.yml
conda activate chronopack
pip install .
```

For active development:

```bash
conda env create -f env/environment.yml
conda activate chronopack
pip install -e .
pytest
```

## Dependency Ownership

- `pyproject.toml` defines the Python package metadata and Python runtime
  dependencies.
- The R backend (`src/chronopack/_backends/r/`) needs these R packages, all
  provided by `env/environment.yml`: BiocParallel, Matrix, RColorBrewer,
  SingleCellExperiment, SummarizedExperiment, circlize, cluster, dplyr,
  fdrtool, ggplot2, ggrepel, matrixStats, patchwork, purrr, scales,
  shadowtext, slingshot, tibble, tidyr and tradeSeq (testthat for R-side
  tests).
- `env/environment.yml` and lock/snapshot files define the verified mixed
  Python/R workflow environment.

The Python package does not assume that `pip` can install the R/Bioconductor
backend.

## Checking the Environment

`env/README.md` describes the environment files and a smoke test for an
activated environment.
