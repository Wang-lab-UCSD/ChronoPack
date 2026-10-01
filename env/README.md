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

## Checking an Environment

After creating and activating the environment (and `pip install -e .`),
`bash env/smoke_test_chronopack_hpc.sh` checks the principal Python imports,
the R backend packages and the bundled R resources. On this cluster, run it
inside a Slurm allocation with `TMPDIR` pointed at node-local storage.

## Standard Developer Workflow

Inside an activated ChronoPack runtime:

```bash
python -m pip install -e .
python -m pytest
```

The normal activated workflow should not require manually setting `PATH` for
`rpy2`. Cluster-specific fallbacks such as explicit `TMPDIR`, `MPLCONFIGDIR`,
or `R_HOME` should remain deployment details rather than package behavior.

## Refreshing The Snapshots

Create a fresh reference prefix from the human spec on fast local storage,
verify it, then refresh the snapshots:

```bash
conda env create --prefix /tmp/chronopack-ref -f env/environment.yml
/tmp/chronopack-ref/bin/python -m pip install -e .
env/smoke_test_chronopack_hpc.sh
conda list --prefix /tmp/chronopack-ref --explicit > env/conda-linux-64.lock
/tmp/chronopack-ref/bin/python -m pip freeze --exclude chronopack > env/pip-freeze.txt
```
