# ChronoPack Environment

This directory separates the reproducible mixed-language workflow environment
from package metadata.

## Canonical Definition

The canonical ChronoPack environment definition is the committed dependency
specification and snapshots:

- `env/environment.yml`: human-maintained compatibility-level specification.
- `env/conda-linux-64.lock`: exact conda explicit package snapshot.
- `env/pip-freeze.txt`: exact pip package snapshot for the dependency runtime.

These files, not any one installed prefix, define the reference environment.
The temporary validation prefix `/tmp/chronopack-ref-20260918113425` was used
only to prove the unified stack works.

## HPC Deployment Model

On this cluster, `/new-stg/home/...` is persistent shared NFS, while `/scratch`
and `/tmp` on compute nodes are fast node-local storage. Installing or running a
large live conda prefix directly from NFS is not the preferred workflow.

The current deployment model is:

1. Build and verify the ChronoPack dependency environment on fast node-local
   storage.
2. Pack it with `conda-pack`.
3. Store the packed archive on persistent NFS:
   `artifacts/env/chronopack-linux-64-20260918.tar.gz`
4. For each Slurm allocation/job, unpack the archive to node-local scratch,
   run `conda-unpack`, activate it, and install this checkout in editable mode.

The archive is a deployment artifact and should not be committed to Git.
Because `conda-pack` refuses editable installs, the packed archive intentionally
contains the dependency runtime rather than an editable ChronoPack checkout.

## Slurm Usage

From an interactive Slurm allocation:

```bash
cd /new-stg/home/alex/tmp/ChronoPack
source env/activate_chronopack_hpc.sh
bash env/smoke_test_chronopack_hpc.sh
```

The activation helper sets node-local defaults for:

- `CHRONOPACK_LOCAL`
- `TMPDIR`
- `CONDA_PKGS_DIRS`
- `MPLCONFIGDIR`
- `CHRONOPACK_ENV_PREFIX`

It prefers `/scratch/$USER/chronopack-$SLURM_JOB_ID` when available and falls
back to `/tmp`. Set `CHRONOPACK_ENV_ARCHIVE` if the packed archive lives
somewhere other than the default path above.

Batch verification can be run with:

```bash
sbatch jobs/verify_chronopack_packed_env.slurm
```

## Standard Developer Workflow

Inside an activated ChronoPack runtime:

```bash
python -m pip install -e .
python -m pytest
```

The normal activated workflow should not require manually setting `PATH` for
`rpy2`. Cluster-specific fallbacks such as explicit `TMPDIR`, `MPLCONFIGDIR`,
or `R_HOME` should remain deployment details rather than package behavior.

## Rebuilding The Archive

Create a fresh reference prefix from the human spec on fast local storage, then
pack after verification:

```bash
conda env create --prefix /tmp/chronopack-ref -f env/environment.yml
/tmp/chronopack-ref/bin/python -m pip install -e .
env/smoke_test_chronopack_hpc.sh
/tmp/chronopack-ref/bin/python -m pip uninstall -y chronopack
conda-pack -p /tmp/chronopack-ref -o artifacts/env/chronopack-linux-64-YYYYMMDD.tar.gz
```

After a verified rebuild, refresh snapshots:

```bash
conda list --prefix /tmp/chronopack-ref --explicit > env/conda-linux-64.lock
/tmp/chronopack-ref/bin/python -m pip freeze > env/pip-freeze.txt
```

Longer term, Apptainer may be a useful HPC distribution mechanism, but it is
not implemented in this pass.
