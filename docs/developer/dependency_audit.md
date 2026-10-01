# ChronoPack Dependency Audit

Date: 2026-09-18

## Source Environments Inspected

- `wot2`: main Python/WOT workflow.
- `omnipath`: GRN/OmniPath workflow.
- `tradeseq_ce`: primary R/tradeSeq workflow.
- `r-fdr`: R/fdrtool target-gene-score workflow.

The goal is not to merge complete exports. The unified environment should cover
the packages actually imported by active ChronoPack code and notebooks.

## Imported Dependency Surface

Core Python:

- `anndata`, `scanpy`, `numpy`, `pandas`, `scipy`, `scikit-learn`,
  `statsmodels`, `matplotlib`, `seaborn`, `gseapy`

WOT / cell-state trajectories:

- `wot`, `scanpy`, `scipy.sparse`, `sklearn`

Python GRN/C2C:

- `networkx`, `omnipath`, `matplotlib`, `statsmodels`

R runtime / Bioconductor backend:

- `tradeSeq`, `SingleCellExperiment`, `slingshot`, `BiocParallel`,
  `SummarizedExperiment`, `Matrix`

R plotting / descriptors:

- `ggplot2`, `dplyr`, `tidyr`, `tibble`, `matrixStats`, `RColorBrewer`,
  `patchwork`, `circlize`, `ggrepel`, `shadowtext`, `purrr`, `scales`

R target-gene scoring:

- `fdrtool`

Development/test-only:

- `pytest`, `hatchling`, `ipykernel`

## Working Versions Observed

`wot2`:

- Python 3.13.2
- `anndata` 0.11.4
- `scanpy` 1.11.1
- `numpy` 2.2.6
- `pandas` 2.2.3
- `scipy` 1.15.3
- `scikit-learn` metadata showed 1.5.2, import reported 1.7.1
- `statsmodels` 0.14.4
- `matplotlib` 3.10.0
- `networkx` 3.4.2
- `wot` 1.0.8.post2

`omnipath`:

- Python 3.11.4
- `numpy` 2.3.4
- `pandas` 2.3.3
- `scipy` 1.16.2
- `matplotlib` 3.10.7
- `networkx` 3.5
- `omnipath` 1.0.12
- Does not include `scanpy`, `anndata`, `wot`, `scikit-learn`, or `statsmodels`.

`tradeseq_ce`:

- R 4.4.3
- `tradeSeq` 1.20.0
- `SingleCellExperiment` 1.28.0
- `slingshot` 2.14.0
- `BiocParallel` 1.40.0
- `SummarizedExperiment` 1.36.0
- `ggplot2` 4.0.0
- `dplyr` 1.1.4
- `tidyr` 1.3.1
- `tibble` 3.3.0
- `matrixStats` 1.5.0
- `patchwork` 1.3.2
- `circlize` 0.4.17
- `ggrepel` 0.9.8
- Missing `fdrtool` and `shadowtext`.

`r-fdr`:

- R 4.5.2
- `fdrtool` 1.2.18
- `ggplot2` 4.0.1
- `dplyr` 1.1.4
- Missing `tradeSeq`, `SingleCellExperiment`, `slingshot`,
  `BiocParallel`, `SummarizedExperiment`, `matrixStats`, `patchwork`,
  `circlize`, and `shadowtext`.

## Identified Conflicts / Risks

- The current workflows are split across R 4.4/Bioc tradeSeq and R 4.5/fdrtool
  environments. The unified environment targets R 4.4/Bioc 3.20 first because
  tradeSeq is the stricter backend constraint. The proposed unified spec solves
  with `r-fdrtool` in this stack.
- `wot2` uses Python 3.13, while the GRN `omnipath` environment uses Python
  3.11. The unified environment targets Python 3.11 because it is a safer
  common denominator for rpy2/R integration.
- `wot` is installed from PyPI in `wot2`; it is specified in the pip section of
  `environment.yml`.
- Base Python on this machine lacks `pytest`, and the default conda/mamba cache
  can be unwritable. Use the project-local cache instructions in `env/README.md`
  if needed.

## Unified Environment Verification

Human-maintained spec:

- `env/environment.yml`

Observed solve result:

- The mixed Python/R/Bioconductor spec solved successfully with Python 3.11,
  R 4.4.3, `bioconductor-tradeseq` 1.20.0,
  `bioconductor-singlecellexperiment` 1.28.0,
  `bioconductor-slingshot` 2.14.0, `r-fdrtool` 1.2.18,
  `rpy2` 3.6.6, `scanpy` 1.11.5, `anndata` 0.12.19,
  `omnipath` 1.0.12, and `wot` 1.0.8.post2 via pip.
- No package-level incompatibility was observed during solving.

Filesystem/install diagnosis:

- The repository and `/new-stg/home/alex/miniconda3/envs` are on NFS.
- `/tmp` is a local tmpfs with enough space for the reference environment.
- Conda transactions on NFS solved but stalled during extraction/linking.
- A tiny test environment on `/tmp` completed normally, confirming that the
  problem was install/extraction filesystem behavior rather than dependency
  resolution.
- The reference ChronoPack environment was successfully installed on `/tmp`.

Verified reference prefix:

- `/tmp/chronopack-ref-20260918113425`

Additional direct dependency discovered during verification:

- `cluster` is used by retained R plotting code (`cluster::silhouette`) and is
  now declared as `r-cluster` in `env/environment.yml` and listed with the
  other R backend packages in `docs/user/installation.md`.

Phase E smoke-test results:

- `pip install -e .` succeeded in the reference prefix.
- `import chronopack as cp` from outside the repository root succeeded.
- `cp.io`, `cp.tl`, and `cp.pl` import successfully.
- Principal Python imports succeeded: `numpy`, `pandas`, `scipy`,
  `scikit-learn`, `statsmodels`, `anndata`, `scanpy`, `matplotlib`, `seaborn`,
  `gseapy`, `networkx`, `omnipath`, and `wot`.
- WOT import succeeded and a toy `OTModel.compute_transport_map(0, 1)` call
  returned a 2 x 2 AnnData transport map.
- R 4.4.3 starts from the reference prefix.
- R backend package loads succeeded for `tradeSeq`, `SingleCellExperiment`,
  `fdrtool`, `ggplot2`, `dplyr`, `tidyr`, `tibble`, `matrixStats`,
  `RColorBrewer`, `patchwork`, `circlize`, `cluster`, `ggrepel`,
  `shadowtext`, `purrr`, `scales`, `slingshot`, `BiocParallel`, and
  `SummarizedExperiment`.
- `rpy2` initializes the prefix R when the prefix `bin` directory is on `PATH`;
  a Python -> R -> Python vector round trip succeeded.
- Bundled R backend resources under `src/chronopack/_backends/r/` are
  discoverable through the package-resource mechanism from outside the
  repository root.
- `pytest` passed: 22 tests.

Reproducibility artifacts:

- `env/conda-linux-64.lock`: exact conda explicit package snapshot for the
  verified prefix.
- `env/pip-freeze.txt`: exact pip package snapshot, including `wot` and `POT`.

Operational notes:

- When invoking Python by absolute path instead of activating the environment,
  prepend the prefix `bin` directory to `PATH` before using `rpy2`, so rpy2 can
  locate the matching `R`.
- Matplotlib may need a writable `MPLCONFIGDIR` on this cluster because the
  default home config directory is not writable in the sandbox.

## Packed HPC Deployment Verification

Cluster filesystem behavior:

- `/new-stg/home/...` is persistent shared NFS.
- `/tmp` is fast node-local tmpfs.
- `/scratch/$USER/...` is fast writable node-local ext4 inside compute
  allocations, but it is not visible from the login node after the allocation.
- `$TMPDIR` points to NFS by default and should be overridden inside Slurm jobs.

Deployment decision:

- Keep `env/environment.yml`, `env/conda-linux-64.lock`, and
  `env/pip-freeze.txt` as the canonical environment definition.
- (Removed 2026-09-30, unused since this audit: the packed-archive deployment
  below, `env/activate_chronopack_hpc.sh` and
  `jobs/verify_chronopack_packed_env.slurm`.)
- Use a packed conda environment as the immediate HPC deployment artifact.
- Store the packed archive on persistent NFS, but hydrate runtime prefixes onto
  node-local scratch for each allocation/job.
- Keep the live ChronoPack source checkout outside the archive; install it
  editable after unpacking.

Packed archive:

- `artifacts/env/chronopack-linux-64-20260918.tar.gz`
- Size: approximately 1.2 GB.
- The archive is not committed to Git.
- The archive contains the dependency runtime only. `conda-pack` rejected the
  original attempt with editable `chronopack` installed, so editable install is
  handled by `env/activate_chronopack_hpc.sh`.

Reusable helpers:

- `env/activate_chronopack_hpc.sh`: sourceable Slurm/runtime activation helper.
- `env/smoke_test_chronopack_hpc.sh`: smoke-test matrix for an activated
  runtime.
- `jobs/verify_chronopack_packed_env.slurm`: batch reconstruction test.

Packed reconstruction test:

- Slurm job: `925682`
- Node-local runtime prefix:
  `/scratch/alex/chronopack-925682/env`
- Archive unpack, `conda-unpack`, activation, and editable install succeeded.
- ChronoPack imported from outside the repository root.
- `cp.io`, `cp.tl`, and `cp.pl` imported successfully.
- Principal Python imports succeeded, including `scanpy`, `anndata`,
  `omnipath`, and `wot`.
- WOT toy transport map smoke test returned a 2 x 2 transport map.
- R 4.4.3 started successfully.
- `tradeSeq`, `SingleCellExperiment`, `fdrtool`, and retained R plotting
  dependencies loaded successfully.
- rpy2 completed a Python -> R -> Python round trip.
- Bundled R backend resources were found through package-resource lookup.
- `pytest` passed: 22 tests.

Conclusion:

- The packed environment reconstruction strategy is verified for this HPC
  filesystem model.
- ChronoPack should use the committed spec/snapshots as the reproducibility
  source of truth and hydrate job-local runtimes from the persistent packed
  artifact.
