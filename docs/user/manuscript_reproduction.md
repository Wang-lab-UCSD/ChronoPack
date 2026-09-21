# Manuscript Reproduction Runtime

ChronoPack has two environment concepts:

- The default ChronoPack runtime in `env/environment.yml`, intended as the
  maintained mixed Python/R development environment.
- The manuscript-reproduction runtime in `env/manuscript-reproduction.yml`,
  intended to reproduce the historical manuscript analyses as closely as the
  original dependency stack permits.

The manuscript profile intentionally preserves version-sensitive behavior from
the historical `wot2` and tradeSeq R environments. It is not a recommendation
that future ChronoPack analyses must use these versions.

## Environment Files

- `env/manuscript-reproduction.yml`: Conda environment specification for the
  mixed Python/R manuscript runtime.
- `env/manuscript-pip-overrides.txt`: post-install pip overrides needed to
  reproduce historical import-time behavior that is not resolver-consistent
  with modern package metadata.
- `env/manuscript-conda-linux-64.lock`: exact explicit Conda snapshot from the
  verified manuscript prefix.
- `env/manuscript-pip-freeze.txt`: exact pip snapshot from the verified
  manuscript prefix.

The verified prefix used in this pass was:

```text
/new-stg/home/alex/miniconda3/envs/chronopack-manuscript-20260919
```

## Creation Workflow

The environment requires two steps because the historical `wot2` runtime
contained import-time versions that conflict with declared package metadata.

```bash
conda env create --prefix /path/to/chronopack-manuscript \
  -f env/manuscript-reproduction.yml

/path/to/chronopack-manuscript/bin/python -m pip install --no-deps \
  -r env/manuscript-pip-overrides.txt

PATH=/path/to/chronopack-manuscript/bin:$PATH \
CONDA_PREFIX=/path/to/chronopack-manuscript \
/path/to/chronopack-manuscript/bin/python -m pip install -e . --no-deps
```

Use an activated environment, or otherwise put the prefix `bin` directory first
on `PATH`, before running R-backed tests or workflows. `rpy2` needs to find the
matching R executable from the same prefix.

## Recovered Version-Sensitive Python Stack

The manuscript preprocessing stack recovered from `wot2` is:

- Python `3.13.2`
- Scanpy `1.11.1`
- AnnData `0.11.4`
- NumPy `2.2.6`
- SciPy `1.15.3`
- pandas `2.2.3`
- scikit-image `0.25.2`
- python-igraph / igraph `1.0.0`
- leidenalg `0.10.2`
- texttable `1.7.0`
- networkx `3.4.2`
- WaddingtonOT `1.0.8.post2`
- umap-learn `0.5.7`
- numba `0.61.2`
- llvmlite `0.44.0`
- pynndescent `0.5.13`

Historical `wot2` has a metadata/import mismatch for scikit-learn:

- `import sklearn; sklearn.__version__` reports `1.7.1`.
- Package metadata / `pip freeze` reports `scikit-learn==1.5.2`.

The manuscript prefix preserves the import-time behavior with
`scikit-learn==1.7.1`. This requires the explicit no-deps override because
`scanpy==1.11.1` declares `scikit-learn>=1.1,<1.6`.

Similarly, `leidenalg==0.10.2` and `igraph==1.0.0` reproduce historical imports
but are not accepted together by pip's dependency resolver.

## Recovered R Stack

The manuscript R/tradeSeq stack verified in the manuscript prefix is:

- R `4.4.3`
- tradeSeq `1.20.0`
- SingleCellExperiment `1.28.0`
- BiocParallel `1.40.0`
- mgcv `1.9.3`
- fdrtool `1.2.18`
- mclust `6.1.3`
- ggplot2 `4.0.0`
- cluster `2.1.8.1`
- Matrix `1.7.4`
- MASS `7.3.65`
- slingshot `2.14.0`

## Verified Preprocessing Reproduction

The known version-sensitive preprocessing boundary was K562 TNF-alpha +
puromycin. With the full recovered stack, the package-backed preprocessing
notebook reproduced the historical reference output:

- `adata_full.chk0.h5ad`: `4034 x 18833`
- `adata_hvg.chk0.h5ad`: `4034 x 3000`
- cell IDs and order: exact match
- gene IDs and order: exact match
- `X`: exact match
- shared layers (`counts`, `cpm`, `log1p`): exact match
- common metadata: exact match except `doublet_score`
- `doublet_score` maximum absolute difference:
  `0.015541014255667104`
- `predicted_doublet`: exact match

The intermediate attempt before aligning `numba`, `llvmlite`, and
`pynndescent` produced `4052` cells instead of `4034`, confirming that the
graph/neighbor stack was part of the historical reproduction issue.

## Scope And Caveats

The package refactor equivalence result remains separate from environment
equivalence:

- package-backed preprocessing matched the old implementation under historical
  `wot2`;
- the default unified ChronoPack runtime remains suitable for package
  development, but may produce version-sensitive graph/Leiden drift;
- exact manuscript preprocessing reproduction should use this manuscript
  runtime.

Historical tradeSeq fits did not fix a random seed. Small fitted-output drift
from rerunning tradeSeq should be classified as historical stochasticity unless
it exceeds the variability already observed from rerunning the old workflow.

This pass verified the known preprocessing discrepancy. Broader downstream
manuscript-output regressions should use frozen upstream artifacts where exact
manuscript state is required.
