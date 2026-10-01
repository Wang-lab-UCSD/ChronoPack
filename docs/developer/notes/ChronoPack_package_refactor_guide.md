# ChronoPack Package Refactoring Guide

## Purpose

This document is intended to guide the next major refactoring of **ChronoPack** from a collection of analysis notebooks and utilities into a coherent, maintainable, user-facing scientific software package.

The primary design goal is to make ChronoPack feel like a natural member of the modern single-cell Python ecosystem while preserving the existing strengths of the project:

- high-resolution time-series single-cell analysis,
- NB-GAM-based gene dynamics,
- gene-trajectory clustering and characterization,
- fate/trajectory analysis,
- cross-context mode conservation,
- GRN inference and benchmarking,
- cell-cell communication analysis,
- plotting and visualization,
- mixed Python/R statistical backends.

The package should be designed around a clean public API, reproducible data contracts, modular backends, and explicit separation between computation, data storage, and plotting.

---

# 1. Recommended high-level design

Use a **Python-first, AnnData-centered package** with three primary public namespaces:

```python
import chronopack as cp

cp.io
cp.tl
cp.pl
```

These should be the main user-facing entry points.

Recommended roles:

- `chronopack.io`
  - reading/writing
  - input validation
  - AnnData schema utilities
  - R/SCE conversion
  - result export/import
  - large sidecar result handling

- `chronopack.tl`
  - scientific computation and analysis
  - gene dynamics
  - gene clustering
  - trajectory characterization
  - cell fate inference
  - driver-gene analysis
  - cross-context conservation
  - GRN analysis
  - C2C analysis

- `chronopack.pl`
  - plotting only
  - should consume already-computed results
  - should not silently rerun analysis

Avoid adding a `pp` namespace unless ChronoPack develops preprocessing operations that are genuinely specific to time-resolved analysis. Standard single-cell preprocessing should generally remain delegated to Scanpy.

---

# 2. Core architectural principle

## AnnData should be the central data container, but not the universal result object

ChronoPack should treat `AnnData` as its main data contract for anything naturally indexed by cells or genes.

Use AnnData aggressively for:

- cell metadata,
- gene metadata,
- embeddings,
- fate probabilities,
- transition matrices,
- fitted gene trajectories,
- trajectory descriptors,
- dynamic-mode assignments,
- compact provenance.

Do **not** force every possible result into AnnData.

Some outputs are inherently better represented as dedicated objects or sidecar tables:

- transport maps,
- sampled trajectories,
- GRN edge tables,
- regulator-level rankings,
- C2C ligand-receptor tables,
- cross-dataset comparison results,
- very large sparse matrices,
- benchmark outputs.

The preferred rule is:

| Result type | Storage location |
|---|---|
| one value per cell | `adata.obs` |
| multiple values per cell | `adata.obsm` |
| square cell-cell matrix | `adata.obsp` |
| one value per gene | `adata.var` |
| multiple values per gene | `adata.varm` |
| square gene-gene matrix | `adata.varp` |
| expression-like matrix | `adata.layers` |
| small metadata/provenance | `adata.uns["chronopack"]` |
| large/non-aligned table | dedicated result object and/or Parquet |
| large sparse matrix | dedicated result object and/or NPZ/Zarr |
| fitted model state | dedicated model/result object |

This distinction should be enforced throughout the package.

---

# 3. Recommended public API

The package should expose task-oriented functions rather than notebook-derived implementation names.

Example:

```python
import chronopack as cp
```

## `cp.io`

Suggested public API:

```python
cp.io.read_h5ad(...)
cp.io.write_h5ad(...)
cp.io.validate(...)
cp.io.export_results(...)
cp.io.import_results(...)
```

Potential internal responsibilities:

- checking required AnnData fields,
- validating count matrices,
- validating time variables,
- checking unique cell/gene IDs,
- converting AnnData to SingleCellExperiment when required,
- exporting large sidecar results,
- restoring analysis metadata.

---

## `cp.tl`

Suggested organization:

```text
cp.tl

CELL STATE / FATE
    transition_clustering
    FateModel
    fate_drivers

GENE DYNAMICS
    fit_gene_dynamics
    cluster_gene_trajectories
    characterize_gene_trajectories

COMPARATIVE DYNAMICS
    mode_conservation
    mode_switchers

REGULATORY ANALYSIS
    grn.*
    c2c.*
```

Possible API:

```python
cp.tl.fit_gene_dynamics(...)
cp.tl.cluster_gene_trajectories(...)
cp.tl.characterize_gene_trajectories(...)

cp.tl.FateModel(...)
cp.tl.fate_drivers(...)

cp.tl.mode_conservation(...)
cp.tl.mode_switchers(...)

cp.tl.grn.infer(...)
cp.tl.grn.rank_regulators(...)
cp.tl.grn.benchmark(...)

cp.tl.c2c.infer(...)
cp.tl.c2c.temporal_analysis(...)
```

Do not implement generic Scanpy replacements unless ChronoPack adds a genuinely time-specific method.

For example, ordinary Leiden clustering should remain:

```python
sc.tl.leiden(adata)
```

rather than:

```python
cp.tl.cluster_cells(adata)
```

unless the ChronoPack version uses temporal structure or transition kernels in a distinct way.

---

## `cp.pl`

Suggested public plotting API:

```python
cp.pl.transition_umap(...)
cp.pl.cell_fates(...)

cp.pl.volcano(...)
cp.pl.gam_fit(...)
cp.pl.mode_summary(...)

cp.pl.trajectory_features(...)
cp.pl.trajectory_feature_curve(...)

cp.pl.precision_recall(...)

cp.pl.mode_conservation(...)
cp.pl.mode_sankey(...)
cp.pl.mode_chord(...)
cp.pl.mode_switchers(...)
```

Important rule:

> `cp.pl` functions should visualize previously computed state. They should not silently rerun modeling, clustering, or testing.

For example:

```python
cp.tl.fit_gene_dynamics(adata)
cp.pl.gam_fit(adata, genes=["NFKBIA", "JUNB"])
```

not:

```python
cp.pl.gam_fit(adata, ...)
# internally refits the NB-GAM
```

---

# 4. Proposed source tree

Recommended package layout:

```text
ChronoPack/
├── pyproject.toml
├── README.md
├── LICENSE
├── CHANGELOG.md
│
├── src/
│   └── chronopack/
│       ├── __init__.py
│       │
│       ├── io/
│       │   ├── __init__.py
│       │   ├── _read.py
│       │   ├── _write.py
│       │   ├── _convert.py
│       │   ├── _export.py
│       │   └── _validate.py
│       │
│       ├── tl/
│       │   ├── __init__.py
│       │   │
│       │   ├── _gene_dynamics.py
│       │   ├── _gene_clustering.py
│       │   ├── _trajectory_features.py
│       │   ├── _conservation.py
│       │   │
│       │   ├── fates/
│       │   │   ├── __init__.py
│       │   │   ├── _model.py
│       │   │   ├── _drivers.py
│       │   │   └── _transport.py
│       │   │
│       │   ├── grn/
│       │   │   ├── __init__.py
│       │   │   ├── _infer.py
│       │   │   ├── _rank.py
│       │   │   ├── _benchmark.py
│       │   │   └── _results.py
│       │   │
│       │   └── c2c/
│       │       ├── __init__.py
│       │       ├── _infer.py
│       │       ├── _temporal.py
│       │       └── _results.py
│       │
│       ├── pl/
│       │   ├── __init__.py
│       │   ├── _fates.py
│       │   ├── _gene_dynamics.py
│       │   ├── _trajectory_features.py
│       │   ├── _conservation.py
│       │   ├── _grn.py
│       │   └── _c2c.py
│       │
│       ├── _core/
│       │   ├── __init__.py
│       │   ├── _keys.py
│       │   ├── _schema.py
│       │   ├── _validation.py
│       │   ├── _provenance.py
│       │   ├── _types.py
│       │   └── _constants.py
│       │
│       ├── _backends/
│       │   ├── __init__.py
│       │   ├── _python/
│       │   │   └── ...
│       │   └── _r/
│       │       ├── __init__.py
│       │       ├── _tradeseq.py
│       │       ├── _sce.py
│       │       └── ...
│       │
│       └── _utils/
│           ├── __init__.py
│           ├── _logging.py
│           ├── _parallel.py
│           ├── _random.py
│           └── _misc.py
│
├── tests/
│   ├── conftest.py
│   ├── test_io.py
│   ├── test_gene_dynamics.py
│   ├── test_gene_clustering.py
│   ├── test_trajectory_features.py
│   ├── test_conservation.py
│   ├── test_fates.py
│   ├── test_grn.py
│   └── test_c2c.py
│
├── docs/
│   ├── api/
│   ├── tutorials/
│   ├── concepts/
│   └── developer/
│
└── examples/
    └── ...
```

---

# 5. Public versus internal modules

Use leading underscores for implementation details that should not become part of the stable public API.

For example:

```text
chronopack/tl/_gene_dynamics.py
```

contains implementation code, while:

```python
from ._gene_dynamics import fit_gene_dynamics
```

is exposed through:

```python
chronopack.tl.fit_gene_dynamics
```

Users should generally not import internal modules directly.

Good:

```python
cp.tl.fit_gene_dynamics(...)
```

Avoid encouraging:

```python
from chronopack.tl._gene_dynamics import _fit_one_gene
```

This gives freedom to refactor internals later without breaking user code.

---

# 6. Define the AnnData schema explicitly

ChronoPack should define a documented, versioned schema for where its results live.

Do not scatter ad hoc strings throughout the source.

Centralize keys in:

```text
chronopack/_core/_keys.py
```

Example:

```python
GENE_MODE_KEY = "cp_mode"
WALD_KEY = "cp_wald"
QVALUE_KEY = "cp_qvalue"
FATE_PROBS_KEY = "cp_fate_probs"
TRANSITION_KEY = "cp_transition_matrix"
```

Better yet, support result-specific prefixes through `key_added`.

---

# 7. Recommended AnnData input schema

Typical inputs may include:

```python
adata.X
```

for normalized expression, depending on the workflow.

Raw counts should preferably live in:

```python
adata.layers["counts"]
```

Cell metadata:

```python
adata.obs["time"]
adata.obs["condition"]
adata.obs["sample"]
adata.obs["cell_type"]
```

Optional standard Scanpy outputs:

```python
adata.obsm["X_pca"]
adata.obsm["X_umap"]
```

Do not rename standard Scanpy fields unnecessarily.

---

# 8. Gene-dynamics storage conventions

Gene-level scalar results should live in `adata.var`.

Example:

```python
adata.var["cp_wald"]
adata.var["cp_pvalue"]
adata.var["cp_qvalue"]

adata.var["cp_peak_log2fc"]
adata.var["cp_onset"]
adata.var["cp_peak_time"]

adata.var["cp_mean_curvature"]
adata.var["cp_max_curvature"]
adata.var["cp_mean_turn"]
adata.var["cp_max_turn"]

adata.var["cp_owp"]
adata.var["cp_owp_ddt"]

adata.var["cp_cluster"]
adata.var["cp_mode"]
```

Fitted trajectories should be stored in `adata.varm`:

```python
adata.varm["cp_gam_fit"]
```

with shape:

```text
n_genes x n_evaluation_timepoints
```

The matching evaluation grid belongs in metadata:

```python
adata.uns["chronopack"]["gene_dynamics"]["time_grid"]
```

If derivatives are stored:

```python
adata.varm["cp_gam_derivative"]
```

If posterior mode probabilities are available:

```python
adata.varm["cp_mode_probs"]
```

---

# 9. Use `key_added` throughout the package

Follow the Scanpy-style pattern:

```python
cp.tl.fit_gene_dynamics(
    adata,
    time_key="time",
    key_added="tnf",
)
```

This should allow multiple analyses to coexist.

Possible output names:

```python
adata.var["tnf_wald"]
adata.var["tnf_qvalue"]
adata.var["tnf_mode"]
adata.varm["tnf_fit"]
```

and metadata:

```python
adata.uns["chronopack"]["tnf"]
```

This is especially important when comparing:

- K562 TNF,
- K562 TNF + puromycin,
- THP-1 LPS,
- different model settings,
- different clustering solutions.

Avoid hard-coding a single global result key.

---

# 10. Provenance must be first-class

Every major analysis step should record:

- ChronoPack version,
- backend,
- backend version,
- input keys,
- random seed,
- important parameters,
- timestamp if appropriate,
- model settings,
- result keys written.

Example:

```python
adata.uns["chronopack"]["tnf"] = {
    "analysis": "gene_dynamics",
    "chronopack_version": "0.1.0",
    "backend": "tradeseq",
    "backend_version": "...",
    "time_key": "time",
    "layer": "counts",
    "n_knots": 6,
    "parameters": {
        "...": "..."
    }
}
```

The exact schema should be standardized.

Do not allow every module to invent a different provenance format.

---

# 11. tradeSeq should be a backend, not the package identity

ChronoPack's user-facing API should describe the scientific task:

```python
cp.tl.fit_gene_dynamics(...)
```

rather than:

```python
cp.tl.tradeseq_fit(...)
```

Recommended signature:

```python
cp.tl.fit_gene_dynamics(
    adata,
    time_key="time",
    layer="counts",
    backend="tradeseq",
    n_knots=6,
    backend_kwargs=None,
    key_added="gene_dynamics",
)
```

ChronoPack should clearly document that:

- NB-GAM fitting is implemented by tradeSeq,
- ChronoPack orchestrates the workflow,
- ChronoPack adds trajectory extraction,
- ChronoPack adds trajectory characterization,
- ChronoPack adds clustering/mode assignment,
- ChronoPack adds cross-context comparison,
- ChronoPack integrates downstream fate/GRN/C2C analysis.

Do not present the tradeSeq statistical implementation as a new ChronoPack model.

---

# 12. Keep R behind a backend boundary

ChronoPack should expose a Python API even when some statistics are performed in R.

Preferred architecture:

```text
User
  ↓
ChronoPack Python API
  ↓
AnnData
  ↓
backend adapter
  ├── Python
  └── R / tradeSeq
```

Possible location:

```text
chronopack/_backends/_r/_tradeseq.py
chronopack/_backends/_r/_sce.py
```

The R-specific implementation should be isolated from higher-level business logic.

For example:

```python
cp.tl.fit_gene_dynamics(...)
```

may internally:

1. validate AnnData,
2. extract required matrices,
3. convert to `SingleCellExperiment`,
4. call tradeSeq,
5. collect fitted curves and test statistics,
6. write standardized results back to AnnData.

The rest of ChronoPack should not care how the backend was implemented.

---

# 13. Avoid exposing low-level tradeSeq arguments everywhere

The primary function should expose only common, scientifically meaningful parameters.

Example:

```python
cp.tl.fit_gene_dynamics(
    adata,
    time_key="time",
    layer="counts",
    n_knots=6,
    backend="tradeseq",
    backend_kwargs={
        ...
    },
)
```

Do not mirror every argument from `fitGAM` into the top-level ChronoPack signature.

This reduces API churn and preserves the option of supporting other backends later.

---

# 14. Fate analysis should use a model object

Cell fate and transport analyses can contain state that does not fit naturally into AnnData.

Recommended pattern:

```python
fate = cp.tl.FateModel(
    adata,
    time_key="time",
    representation="X_pca",
)

fate.fit()
fate.compute_fates()
fate.compute_driver_genes()
```

The model object may own:

```python
fate.transport_maps
fate.trajectories
fate.parameters
fate.solutions
```

while compact outputs are written back to AnnData:

```python
adata.obs["cp_fate"]
adata.obsm["cp_fate_probs"]
adata.obsp["cp_transition_matrix"]
```

This mirrors the useful pattern seen in OT-based single-cell packages: AnnData is the central biological container, while the full model state lives in a dedicated object.

---

# 15. GRN analysis should return result objects

Do not store enormous edge tables inside `adata.uns`.

Recommended:

```python
grn = cp.tl.grn.infer(adata, ...)
```

Return:

```python
GRNResult
```

Possible interface:

```python
grn.edges
grn.tf_scores
grn.metrics
grn.parameters
grn.metadata
```

Only compact summaries should be written back to AnnData.

For example:

```python
adata.var["cp_grn_target_score"]
```

and:

```python
adata.uns["chronopack"]["grn"] = {
    "method": "...",
    "n_edges": ...,
    "parameters": {...},
}
```

Large edge tables should be exported as:

```text
.parquet
```

rather than embedded into `.h5ad`.

---

# 16. C2C analysis should follow the same pattern

Recommended:

```python
c2c = cp.tl.c2c.infer(...)
```

with:

```python
C2CResult
```

containing:

```python
c2c.interactions
c2c.scores
c2c.temporal_results
c2c.metadata
```

Avoid storing millions of ligand-receptor rows in `adata.uns`.

---

# 17. Cross-context mode conservation should operate on multiple datasets

Do not force inherently multi-dataset comparisons into one AnnData.

Recommended:

```python
res = cp.tl.mode_conservation(
    adata_k562,
    adata_thp1,
    mode_key="cp_mode",
    statistic_key="cp_wald",
)
```

Return:

```python
ModeConservationResult
```

Possible fields:

```python
res.recovery
res.switch_matrix
res.shared_genes
res.switchers
res.parameters
```

Plotting then becomes:

```python
cp.pl.mode_conservation(res)
cp.pl.mode_sankey(res)
cp.pl.mode_chord(res)
cp.pl.mode_switchers(res)
```

This is preferable to stuffing cross-dataset state into either input AnnData.

---

# 18. Result objects should be lightweight and serializable

Prefer Python dataclasses or similarly transparent result containers.

For example:

```python
@dataclass
class ModeConservationResult:
    recovery: pd.DataFrame
    switch_matrix: pd.DataFrame
    shared_genes: pd.DataFrame
    switchers: pd.DataFrame
    parameters: dict
```

Desirable characteristics:

- inspectable,
- easy to serialize,
- easy to test,
- no hidden global state,
- compatible with notebooks,
- minimal magic.

Avoid large monolithic classes that own the entire analysis pipeline.

---

# 19. Storage strategy for large results

Recommended formats:

```text
.h5ad      primary AnnData
.parquet   large tabular results
.npz       sparse matrices
.zarr      very large arrays / chunked outputs
.json      small metadata / manifests
```

Avoid:

- embedding enormous DataFrames in `.uns`,
- pickling large arbitrary Python objects unless unavoidable,
- committing analysis outputs to Git.

`cp.io` should eventually provide helper functions for saving and loading these result objects.

---

# 20. Validation should be centralized

Introduce:

```python
cp.io.validate(...)
```

and internal validators such as:

```python
_validate_adata(...)
_validate_time(...)
_validate_counts(...)
_validate_unique_names(...)
_validate_result_key(...)
```

Validation should catch problems before expensive modeling starts.

Examples:

- duplicate cell IDs,
- duplicate gene IDs,
- missing `time`,
- non-numeric time,
- missing count layer,
- negative counts,
- mismatched dimensions,
- missing required fitted-curve metadata,
- invalid `key_added`,
- inconsistent time grid,
- invalid categorical labels.

Do not duplicate these checks across functions.

---

# 21. Introduce a lightweight schema system

Add:

```text
chronopack/_core/_schema.py
```

This should define expected fields for major analysis stages.

Example conceptual states:

```text
RawAnnData
GeneDynamicsReady
GeneDynamicsFitted
TrajectoryClustered
FateReady
FateComputed
```

This does not need to become a heavy framework.

The goal is simply to centralize assumptions and produce clear error messages.

---

# 22. Separate scientific computation from plotting

Each plot should depend on an explicit computed result.

Example:

```python
cp.tl.cluster_gene_trajectories(adata)
cp.pl.mode_summary(adata)
```

The plotting function should fail clearly if required data are missing.

For example:

```text
ChronoPackError:
No trajectory mode assignments were found.
Run cp.tl.cluster_gene_trajectories(...) first.
```

This is much better than silently recomputing.

---

# 23. Separate orchestration from mathematical kernels

Large functions from notebooks should be broken into:

1. public orchestration function,
2. validation,
3. low-level mathematical/statistical kernels,
4. storage,
5. plotting.

Example:

```text
fit_gene_dynamics()
    ↓
_validate_inputs()
    ↓
_prepare_backend_inputs()
    ↓
_fit_tradeseq_backend()
    ↓
_extract_fitted_curves()
    ↓
_store_gene_dynamics_results()
```

Low-level kernels should be independently testable.

---

# 24. Avoid notebook-style state

The package should not rely on:

- global variables,
- implicitly defined DataFrames,
- notebook execution order,
- environment-specific paths,
- mutable shared state,
- ad hoc lists such as `gamList.ls` unless explicitly encapsulated.

Every function should accept all required inputs explicitly.

Bad:

```python
run_analysis()
# assumes pc_df, sce_ls, gamList_ls already exist globally
```

Good:

```python
result = cp.tl.mode_conservation(
    adata1,
    adata2,
    ...
)
```

---

# 25. Make randomness explicit

Any stochastic function should accept:

```python
random_state=
```

or:

```python
seed=
```

This includes:

- K-means,
- GMM initialization,
- resampling,
- permutation testing,
- trajectory sampling,
- benchmarking,
- null models.

Store the seed in provenance.

---

# 26. Use consistent naming

Prefer verbs for functions:

```python
fit_gene_dynamics
cluster_gene_trajectories
characterize_gene_trajectories
compute_fates
compute_driver_genes
mode_conservation
infer_grn
```

Prefer nouns for result/model classes:

```python
FateModel
GRNResult
C2CResult
ModeConservationResult
```

Avoid vague names such as:

```text
run_analysis
process_data
helper
utils2
final_method
new_function
```

---

# 27. Keep domain-specific subpackages where complexity warrants them

Top-level functions are appropriate for core workflows:

```python
cp.tl.fit_gene_dynamics
cp.tl.cluster_gene_trajectories
cp.tl.mode_conservation
```

More complex areas deserve nested namespaces:

```python
cp.tl.grn
cp.tl.c2c
```

This prevents `cp.tl` from becoming an enormous flat namespace.

Fate analysis may similarly deserve:

```python
cp.tl.fates
```

internally, even if selected public helpers are promoted to `cp.tl`.

---

# 28. Suggested module mapping from current ChronoPack functionality

Current functionality should be reorganized roughly as follows.

## Gene dynamics

```text
NB-GAM fitting
dynamic gene significance
curve extraction
peak amplitude
peak time
```

→

```text
chronopack/tl/_gene_dynamics.py
```

Public:

```python
cp.tl.fit_gene_dynamics
```

---

## Gene trajectory clustering

Current:

- PCA trajectory representation,
- Ward clustering,
- K-means,
- GMM,
- silhouette,
- BIC,
- ICL,
- posterior probabilities,
- broad mode mapping.

→

```text
chronopack/tl/_gene_clustering.py
```

Public:

```python
cp.tl.cluster_gene_trajectories
```

Possible method argument:

```python
method="gmm"
method="kmeans"
method="ward"
```

Do not expose three unrelated APIs unless their inputs/outputs differ substantially.

---

## Trajectory characterization

Current:

- onset,
- peak time,
- peak response,
- curvature,
- turning metrics,
- OWP,
- OWP derivative,
- transcript features.

→

```text
chronopack/tl/_trajectory_features.py
```

Public:

```python
cp.tl.characterize_gene_trajectories
```

Optional lower-level functions may remain internal.

---

## Mode conservation

Current:

- directional recovery,
- top-k source anchoring,
- broader target threshold,
- bidirectional recovery,
- switch matrices,
- Sankey/Chord visualization,
- switcher detection.

→

```text
chronopack/tl/_conservation.py
```

Public:

```python
cp.tl.mode_conservation
cp.tl.mode_switchers
```

Plots:

```python
cp.pl.mode_sankey
cp.pl.mode_chord
cp.pl.mode_switchers
```

---

## Fate analysis

Current:

- transition kernels,
- WaddingtonOT trajectories,
- cell-state transitions,
- fates,
- driver genes.

→

```text
chronopack/tl/fates/
```

Primary API:

```python
fate = cp.tl.FateModel(...)
```

or:

```python
cp.tl.fates.FateModel(...)
```

---

## GRN analysis

Current:

- inference,
- Taiji-derived scores,
- edge ranking,
- TF prioritization,
- benchmark comparison,
- EPR,
- AUPRC ratio,
- GSEA-like regulator enrichment.

→

```text
chronopack/tl/grn/
```

Suggested public API:

```python
cp.tl.grn.infer(...)
cp.tl.grn.rank_regulators(...)
cp.tl.grn.benchmark(...)
```

---

## C2C analysis

Current:

- ligand-receptor analysis,
- temporal ordering,
- lagged comparisons,
- CellChat/LIANA/NicheNet-like integration,
- paracrine timing.

→

```text
chronopack/tl/c2c/
```

Suggested public API:

```python
cp.tl.c2c.infer(...)
cp.tl.c2c.temporal_analysis(...)
```

---

# 29. Plotting architecture

Plots should generally accept:

- AnnData,
- result objects,
- keys specifying stored results,
- standard Matplotlib axes when possible.

Recommended pattern:

```python
ax = cp.pl.gam_fit(
    adata,
    gene="NFKBIA",
    key="tnf",
    ax=ax,
)
```

Where practical:

```python
show=False
```

should allow composition into larger figures.

Avoid functions that automatically save files unless explicitly requested.

Prefer:

```python
ax = cp.pl...
fig.savefig(...)
```

over hidden saving behavior.

---

# 30. Plot function return values

Plot functions should return useful Matplotlib objects.

For simple plots:

```python
return ax
```

For multi-panel figures:

```python
return fig, axes
```

For custom composite visualizations, document the return type clearly.

This makes the plotting API composable.

---

# 31. Logging and verbosity

Use Python's `logging` module rather than raw `print()` statements.

Recommended package logger:

```python
logging.getLogger("chronopack")
```

Functions may expose:

```python
verbose=True
```

only where appropriate.

Avoid filling notebooks with low-level backend output by default.

---

# 32. Progress bars

Long-running loops may optionally use `tqdm`.

Good candidates:

- multi-K clustering evaluation,
- per-gene computations,
- permutation tests,
- GRN benchmarking,
- trajectory sampling.

Progress bars should be suppressible.

---

# 33. Optional dependencies

Do not force every user to install the entire ChronoPack stack.

Potential optional groups:

```text
chronopack
chronopack[tradeseq]
chronopack[fates]
chronopack[grn]
chronopack[c2c]
chronopack[plot]
chronopack[all]
```

R/Bioconductor dependencies may still require Conda/Bioc installation.

The Python package should fail gracefully if an optional backend is missing.

Example:

```text
ImportError:
tradeSeq backend is not available.

Install the ChronoPack tradeSeq environment or configure an R environment
containing tradeSeq.
```

Do not let imports fail globally just because one optional feature is unavailable.

---

# 34. Environment strategy

ChronoPack currently relies on a resolved Conda environment.

Preserve that reproducibility, but separate:

1. Python package metadata,
2. optional R/backend environment,
3. development environment.

Recommended files:

```text
pyproject.toml
environment.yml
environment-lock.yml
```

or equivalent lock tooling.

The package should still be importable independently of the full development environment.

---

# 35. Testing strategy

Every migrated notebook function should gain tests.

At minimum:

## Unit tests

For mathematical kernels:

- onset calculation,
- curvature,
- OWP,
- derivative calculations,
- mode mapping,
- recovery score,
- enrichment statistics.

## Integration tests

For workflows:

```python
cp.tl.fit_gene_dynamics(...)
cp.tl.cluster_gene_trajectories(...)
cp.tl.characterize_gene_trajectories(...)
```

using small synthetic datasets.

## Schema tests

Verify that functions write outputs into the expected AnnData slots.

## Regression tests

For critical existing analyses, retain a tiny fixture with known outputs and ensure refactoring does not alter results unexpectedly.

---

# 36. Use synthetic fixtures instead of real large datasets

Tests should not require:

- multi-GB `.h5ad`,
- full real datasets,
- large GRN tables,
- HPC resources.

Create tiny deterministic fixtures.

Example:

```text
50 cells
20 genes
4 time points
2 conditions
```

and a few synthetic genes with known behaviors:

- impulse-up,
- trend-up,
- trend-down,
- flat/null.

These are ideal for validating ChronoPack.

---

# 37. Backward compatibility is not yet the priority

Because ChronoPack is still being formalized, this is the right time to clean up the API aggressively.

Do not preserve awkward notebook-derived APIs solely for compatibility.

Instead:

- define the desired stable API now,
- migrate internals toward it,
- only then begin considering semantic versioning guarantees.

---

# 38. Documentation structure

Recommended docs:

```text
docs/
├── index
├── installation
├── quickstart
├── concepts/
│   ├── data_model
│   ├── gene_dynamics
│   ├── trajectory_modes
│   ├── fate_analysis
│   ├── grn
│   └── c2c
├── tutorials/
│   ├── gene_dynamics
│   ├── mode_conservation
│   ├── fate_analysis
│   ├── grn
│   └── c2c
├── api/
└── developer/
```

Documentation should clearly distinguish:

- ChronoPack methods,
- external backends,
- expected data schema,
- stored outputs,
- citations.

---

# 39. Citation strategy

ChronoPack should include a citations section describing which underlying tools should be cited.

For example:

```text
Gene dynamics:
    ChronoPack
    tradeSeq

Optimal transport:
    ChronoPack
    WaddingtonOT / relevant OT backend

GRN:
    ChronoPack
    relevant underlying dependency or method

C2C:
    ChronoPack
    relevant underlying dependency
```

Possible future helper:

```python
cp.citation()
```

or:

```python
cp.citation(feature="gene_dynamics")
```

This is useful because ChronoPack integrates several external scientific methods.

---

# 40. Use a scverse-style package template

A good starting point for the formal package skeleton is the scverse package template / cookiecutter.

The structural ideas worth adopting are:

- `src/` layout,
- `pyproject.toml`,
- pytest,
- continuous integration,
- automated documentation,
- linting,
- typing,
- version management,
- package build checks.

Do not necessarily copy every convention blindly.

The goal is to inherit mature scientific-package infrastructure rather than invent packaging machinery from scratch.

---

# 41. Do not over-engineer v0.1

The target should be a clean, maintainable scientific package, not a general-purpose framework.

Avoid prematurely building:

- complex plugin registries,
- dependency injection systems,
- extensive abstract base class hierarchies,
- elaborate configuration frameworks,
- distributed execution engines,
- database abstractions.

Prefer simple Python functions, dataclasses, and transparent state.

---

# 42. Recommended first-pass public API

A reasonable initial target is:

```python
import chronopack as cp
```

## I/O

```python
cp.io.validate(...)
cp.io.read_h5ad(...)
cp.io.write_h5ad(...)
```

## Gene dynamics

```python
cp.tl.fit_gene_dynamics(...)
cp.tl.cluster_gene_trajectories(...)
cp.tl.characterize_gene_trajectories(...)
```

## Fate analysis

```python
cp.tl.FateModel(...)
```

## Comparative dynamics

```python
cp.tl.mode_conservation(...)
cp.tl.mode_switchers(...)
```

## GRN

```python
cp.tl.grn.infer(...)
cp.tl.grn.rank_regulators(...)
cp.tl.grn.benchmark(...)
```

## C2C

```python
cp.tl.c2c.infer(...)
cp.tl.c2c.temporal_analysis(...)
```

## Plotting

```python
cp.pl.transition_umap(...)
cp.pl.cell_fates(...)
cp.pl.volcano(...)
cp.pl.gam_fit(...)
cp.pl.mode_summary(...)
cp.pl.trajectory_features(...)
cp.pl.precision_recall(...)
cp.pl.mode_sankey(...)
cp.pl.mode_chord(...)
cp.pl.mode_switchers(...)
```

---

# 43. Recommended end-user workflow

Target a workflow that looks approximately like this:

```python
import scanpy as sc
import chronopack as cp

adata = sc.read_h5ad("experiment.h5ad")

# Standard preprocessing remains in Scanpy
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# Validate inputs
cp.io.validate(
    adata,
    time_key="time",
    counts_layer="counts",
)

# Fit gene dynamics
cp.tl.fit_gene_dynamics(
    adata,
    time_key="time",
    layer="counts",
    backend="tradeseq",
    n_knots=6,
    key_added="tnf",
)

# Extract kinetic descriptors
cp.tl.characterize_gene_trajectories(
    adata,
    key="tnf",
)

# Cluster trajectories
cp.tl.cluster_gene_trajectories(
    adata,
    key="tnf",
    method="gmm",
)

# Visualize
cp.pl.gam_fit(
    adata,
    genes=["NFKBIA", "JUNB"],
    key="tnf",
)

cp.pl.mode_summary(
    adata,
    key="tnf",
)
```

This should be considered the design target.

---

# 44. Example comparative workflow

```python
res = cp.tl.mode_conservation(
    adata_k562,
    adata_thp1,
    mode_key="cp_mode",
    statistic_key="cp_wald",
    qvalue_key="cp_qvalue",
)

cp.pl.mode_sankey(res)
cp.pl.mode_chord(res)
cp.pl.mode_switchers(res)
```

---

# 45. Example fate workflow

```python
fate = cp.tl.FateModel(
    adata,
    time_key="time",
    representation="X_pca",
)

fate.fit()

fate.compute_fates()
fate.compute_driver_genes()

cp.pl.cell_fates(adata)
cp.pl.transition_umap(adata)
```

---

# 46. Refactoring strategy for Codex

Do not attempt a full rewrite in one step.

Use a staged migration.

## Phase 1 — inventory

Identify:

- notebooks,
- reusable functions,
- duplicated utilities,
- R dependencies,
- global variables,
- shared data structures,
- plotting code,
- large outputs,
- external-tool wrappers.

Create a mapping from existing functions to the target package modules.

Do not change scientific behavior yet.

---

## Phase 2 — package skeleton

Create:

```text
src/chronopack/
chronopack.io
chronopack.tl
chronopack.pl
_core
_backends
tests
docs
```

Set up:

- `pyproject.toml`,
- package installation,
- pytest,
- formatting/linting,
- minimal CI.

At this point, functionality can still be incomplete.

---

## Phase 3 — define the AnnData schema

Before migrating algorithms, define:

- required inputs,
- standard result keys,
- provenance schema,
- `key_added`,
- validation rules.

Implement central key definitions and validation.

This should happen before large-scale code migration.

---

## Phase 4 — migrate gene dynamics first

Gene dynamics should be the first complete vertical slice because it is central to ChronoPack.

Implement:

```python
cp.tl.fit_gene_dynamics
cp.tl.characterize_gene_trajectories
cp.tl.cluster_gene_trajectories
```

with:

- AnnData input,
- tradeSeq backend isolation,
- standardized results,
- provenance,
- tests.

Only after this workflow works cleanly should downstream modules migrate.

---

## Phase 5 — migrate plotting for gene dynamics

Implement:

```python
cp.pl.gam_fit
cp.pl.mode_summary
cp.pl.trajectory_features
```

Ensure plotting uses stored state only.

---

## Phase 6 — migrate comparative analysis

Implement:

```python
ModeConservationResult
cp.tl.mode_conservation
cp.tl.mode_switchers
```

Then plotting:

```python
cp.pl.mode_sankey
cp.pl.mode_chord
cp.pl.mode_switchers
```

---

## Phase 7 — migrate fate analysis

Create:

```python
FateModel
```

and separate transport/model state from compact AnnData outputs.

---

## Phase 8 — migrate GRN and C2C

Treat these as larger subpackages.

Do not force their large outputs into AnnData.

Use result objects and sidecar storage.

---

## Phase 9 — consolidate old notebooks

Once a notebook has been migrated:

- replace internal helper definitions with package imports,
- keep the notebook as a reproducible example or analysis record,
- remove duplicated code from the notebook,
- ensure its outputs match pre-refactor results.

The final notebooks should read like workflows, not source files.

---

# 47. Required behavior-preservation rule

During refactoring:

> Preserve scientific behavior before improving algorithms.

For each migrated method:

1. identify current inputs,
2. identify current outputs,
3. capture reference results,
4. refactor,
5. compare outputs,
6. only then consider changing statistical behavior.

Do not combine architectural refactoring and methodological revision unless explicitly intended.

---

# 48. Codex implementation priorities

Codex should prioritize:

1. clear module boundaries,
2. minimal duplication,
3. explicit inputs/outputs,
4. type hints,
5. concise docstrings,
6. deterministic behavior,
7. tests,
8. AnnData conventions,
9. backend isolation,
10. backwards-verifiable scientific outputs.

Codex should avoid:

- rewriting working scientific logic unnecessarily,
- changing thresholds silently,
- changing defaults without documentation,
- introducing generic abstractions with no current use,
- moving large data files into the package,
- storing analysis outputs in source control.

---

# 49. Suggested coding conventions

Use:

- Python type hints,
- `pathlib.Path`,
- `dataclasses`,
- NumPy,
- pandas,
- scipy,
- AnnData,
- Scanpy where appropriate,
- Matplotlib for plotting,
- standard logging.

Prefer:

```python
def fit_gene_dynamics(
    adata: AnnData,
    *,
    time_key: str,
    layer: str = "counts",
    backend: str = "tradeseq",
    n_knots: int = 6,
    key_added: str = "gene_dynamics",
    backend_kwargs: dict | None = None,
) -> None:
```

or return the modified AnnData only if that convention is adopted consistently.

Avoid mixed conventions where some functions mutate and return `None`, others return AnnData, and others return tuples without clear rationale.

---

# 50. Mutation convention

A Scanpy-like convention is recommended:

- functions operating on AnnData modify it in place by default,
- optional `copy=True` may return a modified copy,
- comparison/model functions return result objects.

Example:

```python
cp.tl.fit_gene_dynamics(adata)
# modifies adata
```

versus:

```python
res = cp.tl.mode_conservation(adata1, adata2)
# returns result object
```

This is intuitive and consistent with the surrounding ecosystem.

---

# 51. Error handling

Define package-specific exceptions if useful:

```python
class ChronoPackError(Exception):
    pass

class ChronoPackSchemaError(ChronoPackError):
    pass

class ChronoPackBackendError(ChronoPackError):
    pass
```

Errors should tell users what to do next.

Good:

```text
ChronoPackSchemaError:
Required counts layer 'counts' was not found.
Provide raw counts in adata.layers['counts'] or pass layer=<name>.
```

Bad:

```text
KeyError: counts
```

---

# 52. Versioning

Use semantic versioning once the public API stabilizes.

Early development can remain:

```text
0.x.y
```

Before `1.0`:

- API changes are expected,
- but should still be documented.

Record package version in analysis provenance.

---

# 53. Package data should remain small

The Python wheel should not contain:

- full `.h5ad`,
- large GRN outputs,
- full benchmark datasets,
- large notebooks,
- HPC-generated artifacts.

Small synthetic examples are fine.

Large example datasets should be externally downloadable.

---

# 54. Git repository organization

Recommended tracked content:

```text
src/
tests/
docs/
examples/
README.md
LICENSE
pyproject.toml
environment files
small notebooks if needed
```

Ignore:

```text
*.h5ad
*.h5
*.loom
large *.rds
benchmark output tables
temporary caches
.ipynb_checkpoints/
verification outputs
large generated figures
```

Use explicit `.gitignore` patterns rather than relying on manual cleanup.

---

# 55. Long-term conceptual structure

ChronoPack should ultimately read conceptually as:

```text
TIME-RESOLVED SINGLE-CELL DATA
            ↓
      Gene dynamics
            ↓
   Trajectory features/modes
            ↓
      Comparative kinetics
            ↓
   Fate / regulatory analysis
            ↓
        GRN / C2C
```

The software architecture should reinforce this scientific story.

---

# 56. Core design principle to preserve

The single most important architectural decision is:

> **AnnData is ChronoPack's lingua franca, but not ChronoPack's universal result object.**

Use AnnData for cell- and gene-indexed state.

Use dedicated objects for:

- models,
- large edge tables,
- cross-dataset comparisons,
- transport solutions,
- benchmark outputs.

This balances interoperability with maintainability.

---

# 57. Immediate recommended next action

Before refactoring more notebooks, implement the following first:

```text
1. package skeleton
2. `_core/_keys.py`
3. `_core/_schema.py`
4. `_core/_validation.py`
5. `_core/_provenance.py`
6. `cp.io.validate`
7. tradeSeq backend wrapper
8. `cp.tl.fit_gene_dynamics`
9. tests for the above
```

Then migrate:

```text
trajectory characterization
→ gene clustering
→ plotting
→ conservation
→ fates
→ GRN
→ C2C
```

This sequencing minimizes the chance of repeatedly redesigning the package interface during migration.

---

# 58. Codex task framing

When giving this document to Codex, the intended task is:

> Refactor ChronoPack incrementally into a formal Python scientific package using the architecture described above. Preserve current scientific behavior unless explicitly instructed otherwise. Make AnnData the primary public data contract, isolate R/tradeSeq behind backend adapters, separate computation from plotting, use lightweight result objects for non-AnnData outputs, add tests during migration, and avoid unnecessary abstraction or algorithmic changes during the structural refactor.

The refactor should prioritize a coherent public API and reproducibility over preserving the current notebook-level implementation structure.

