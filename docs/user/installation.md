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
- `r/chronopackR/DESCRIPTION` records the R-backend dependency intent.
- `env/environment.yml` and lock/snapshot files define the verified mixed
  Python/R workflow environment.

The Python package does not assume that `pip` can install the R/Bioconductor
backend.

## HPC Deployment

The Wang-cluster packed-environment workflow is developer infrastructure, not a
package API. It is documented in `env/README.md` and hydrates a verified packed
environment onto node-local scratch inside Slurm allocations.
