# ChronoPack Mixed Python/R Architecture Guide

## Goal

Refactor ChronoPack into a **Python-centric scientific package** with:

- a clean Python public API,
- AnnData as the primary user-facing data structure,
- R retained as an internal backend where it is already strong,
- tradeSeq and related R workflows hidden behind backend adapters,
- selected R plotting functions preserved rather than unnecessarily rewritten in Python,
- GRN and C2C analysis remaining primarily Python-based.

The implementation language should be treated as an internal detail. Users should interact with one coherent ChronoPack API.

---

# 1. Recommended overall architecture

ChronoPack should be:

```text
User
  ↓
Python ChronoPack API
  ↓
AnnData / result objects
  ↓
Implementation backends
  ├── Python
  └── R
```

The public API should remain Python-centric:

```python
import chronopack as cp

cp.io
cp.tl
cp.pl
```

Users should not normally need to know whether a particular function is implemented in Python or R.

---

# 2. Current workflow and proposed mapping

The current analysis sequence is approximately:

```text
PreprocessDatasets.py
        ↓
GetTrajectories.R
        ↓
PrepareForTradeSeq.R
        ↓
TradeSeqFits.R
        ↓
DimensionalityReductionClustering.R
```

Additional major components:

```text
GRN inference          mostly Python
C2C analysis           mostly Python
```

Recommended package-level abstraction:

```python
adata = cp.io.read(...)

cp.tl.get_trajectories(adata, ...)

cp.tl.fit_gene_dynamics(
    adata,
    backend="tradeseq",
)

cp.tl.characterize_gene_trajectories(adata, ...)

cp.tl.cluster_gene_trajectories(
    adata,
    method="gmm",
)
```

The implementation details should be hidden inside those calls.

In particular, `PrepareForTradeSeq` should probably **not** remain a public function.

It should become an internal backend step:

```text
fit_gene_dynamics()
    ↓
validate AnnData
    ↓
prepare tradeSeq inputs
    ↓
convert to R / SingleCellExperiment
    ↓
fit tradeSeq models
    ↓
run tests
    ↓
extract trajectories
    ↓
write standardized outputs back to AnnData
```

---

# 3. Public API should describe scientific tasks, not implementation details

Good:

```python
cp.tl.fit_gene_dynamics(...)
cp.tl.cluster_gene_trajectories(...)
cp.tl.characterize_gene_trajectories(...)
```

Avoid:

```python
cp.tl.prepare_for_tradeseq(...)
cp.tl.run_tradeseq(...)
cp.tl.run_r_clustering(...)
```

The public API should describe **what the user wants to accomplish**, not which internal language or package performs it.

---

# 4. tradeSeq should be treated as a backend

Recommended public interface:

```python
cp.tl.fit_gene_dynamics(
    adata,
    time_key="time",
    layer="counts",
    backend="tradeseq",
    n_knots=6,
    key_added="gene_dynamics",
)
```

Internally:

```text
chronopack Python API
        ↓
R backend adapter
        ↓
tradeSeq
```

The R-specific implementation should remain isolated from higher-level code.

Do not allow unrelated ChronoPack modules to directly call tradeSeq.

---

# 5. Recommended Python/R boundary

Use AnnData as the primary public container.

Conceptually:

```text
AnnData
   ↓
Python backend adapter
   ↓
R representation
   ↓
tradeSeq / R clustering / R plotting
   ↓
standardized outputs
   ↓
AnnData or result object
```

Potential bridge technologies:

- `rpy2`
- `anndata2ri`
- optionally `Rscript` subprocesses for large batch jobs

---

# 6. Prefer rpy2 for interactive integration

For most interactive use, especially plotting, use:

```text
Python
  ↓
rpy2
  ↓
R
```

Advantages:

- direct function calls,
- easy notebook integration,
- pandas/R conversions,
- AnnData ↔ SingleCellExperiment conversion via `anndata2ri`,
- easy access to ggplot2,
- no need for constant intermediate-file serialization.

This is a good default for:

- short R calls,
- R plotting,
- moderate tradeSeq workflows,
- interactive use.

---

# 7. Consider Rscript subprocesses for heavy compute if needed

For long-running or HPC-heavy R jobs, another valid pattern is:

```text
Python
  ↓
write temporary inputs
  ↓
Rscript backend.R ...
  ↓
write result files
  ↓
Python loads results
```

Advantages:

- better process isolation,
- R crashes do not crash Python,
- easier debugging on HPC,
- easier resource management,
- can be more robust for large batch jobs.

Possible long-term hybrid:

```text
short calls / plotting
    → rpy2

large tradeSeq jobs
    → Rscript subprocess
```

However, start with one approach first. `rpy2` is likely the simplest initial implementation.

---

# 8. Keep the existing R plotting functions

Do **not** rewrite successful R plots in Python merely for implementation-language consistency.

It is completely acceptable for:

```python
cp.pl.mode_summary(...)
```

to internally use:

```text
Python wrapper
   ↓
R backend
   ↓
ggplot2
```

Likewise:

```python
cp.pl.gam_fit(...)
cp.pl.trajectory_features(...)
```

may remain R-backed if the existing visualizations are already publication quality.

The public API does not need to expose the plotting language.

---

# 9. Do not try to convert ggplot objects into Matplotlib objects

Avoid building fragile translation layers between:

```text
ggplot2
```

and:

```text
matplotlib
```

Instead, let R render figures directly.

Recommended output formats:

```text
SVG
PDF
PNG
```

For notebooks, SVG is especially useful.

Conceptual behavior:

```python
cp.pl.mode_summary(adata)
```

may internally:

1. prepare a standardized plotting table,
2. pass the table to R,
3. build a ggplot,
4. render to a temporary SVG,
5. display the SVG inline.

For file export:

```python
cp.pl.mode_summary(
    adata,
    save="mode_summary.pdf",
)
```

the R backend can directly use `ggsave()` or an appropriate graphics device.

---

# 10. Plotting functions should consume standardized results

This is an important refactor.

Existing R plotting functions may currently depend on many analysis objects:

```text
gamList
sce
counts_sub
PT
sf
pc_df
```

Whenever possible, refactor them so they consume simpler standardized inputs.

Prefer:

```r
plot_mode_summary(
    fitted_curves,
    gene_metadata,
    time_grid,
    ...
)
```

or:

```r
plot_mode_summary(plot_df, ...)
```

instead of:

```r
plot_mode_summary(
    gamList,
    sce,
    counts,
    PT,
    sf,
    pc_df,
    ...
)
```

The plotting layer should depend on **analysis outputs**, not the entire internal model state.

---

# 11. Recommended plotting data flow

For example:

```python
cp.tl.fit_gene_dynamics(adata)
```

might populate:

```python
adata.var["cp_qvalue"]
adata.var["cp_mode"]
adata.varm["cp_fit"]
adata.uns["chronopack"]["gene_dynamics"]["time_grid"]
```

Then:

```python
cp.pl.gam_fit(
    adata,
    gene="NFKBIA",
)
```

should do:

```text
AnnData
   ↓
extract standardized fitted-curve data
   ↓
small plotting DataFrame
   ↓
R plotting function
   ↓
ggplot2 output
```

This is much cleaner than reconstructing the full tradeSeq object just to make a plot.

---

# 12. Retain backend model objects only when necessary

Some visualizations or diagnostics may genuinely require a native R model object.

In that case, do not force everything into AnnData.

Possible pattern:

```python
fit = cp.tl.fit_gene_dynamics(
    adata,
    return_model=True,
)
```

with a result object such as:

```python
GeneDynamicsResult
```

containing:

```text
fit.backend
fit.backend_model
fit.metadata
```

Alternatively, serialize the native R object:

```text
tradeseq_fit.rds
```

Use this only when needed.

Most plotting should work from standardized extracted results.

---

# 13. Mixed plotting backends are fine

Different functions under `cp.pl` may use different implementations.

Example:

```text
cp.pl.transition_umap        Python
cp.pl.cell_fates             Python

cp.pl.gam_fit                R
cp.pl.mode_summary           R
cp.pl.trajectory_features    R

cp.pl.volcano                Python or R
cp.pl.precision_recall       Python

cp.pl.mode_sankey            Python
cp.pl.mode_chord             Python or R
cp.pl.mode_switchers         Python or R
```

This is acceptable.

The implementation language should remain invisible unless a user explicitly asks about it.

Do not add a `backend=` argument to plotting functions unless multiple plotting implementations actually exist.

---

# 14. Recommended source-tree structure

A practical mixed-language layout:

```text
ChronoPack/
│
├── pyproject.toml
├── README.md
├── LICENSE
│
├── src/
│   └── chronopack/
│       │
│       ├── io/
│       │   ├── __init__.py
│       │   ├── _read.py
│       │   ├── _write.py
│       │   └── _validate.py
│       │
│       ├── tl/
│       │   ├── __init__.py
│       │   ├── _trajectories.py
│       │   ├── _gene_dynamics.py
│       │   ├── _gene_clustering.py
│       │   ├── _trajectory_features.py
│       │   ├── _conservation.py
│       │   │
│       │   ├── grn/
│       │   │   └── ...
│       │   │
│       │   └── c2c/
│       │       └── ...
│       │
│       ├── pl/
│       │   ├── __init__.py
│       │   ├── _gene_dynamics.py
│       │   ├── _trajectory_features.py
│       │   ├── _conservation.py
│       │   ├── _fates.py
│       │   └── ...
│       │
│       ├── _backends/
│       │   ├── __init__.py
│       │   └── r/
│       │       ├── __init__.py
│       │       ├── _bridge.py
│       │       ├── _tradeseq.py
│       │       ├── _clustering.py
│       │       └── _plotting.py
│       │
│       └── _core/
│           ├── _keys.py
│           ├── _schema.py
│           ├── _validation.py
│           └── _provenance.py
│
└── r/
    └── chronopackR/
        ├── DESCRIPTION
        ├── NAMESPACE
        ├── R/
        │   ├── trajectories.R
        │   ├── tradeseq.R
        │   ├── clustering.R
        │   ├── trajectory_features.R
        │   ├── plots.R
        │   └── utils.R
        │
        └── tests/
            └── testthat/
```

---

# 15. Two reasonable ways to organize the R code

## Simpler v0.1

Keep R source files bundled directly with the Python package:

```text
src/chronopack/_backends/r/
    tradeseq.R
    trajectories.R
    clustering.R
    plotting.R
```

Python wrappers use `rpy2` to source and call them.

This is sufficient for the first formal version.

---

## More mature long-term design

Turn the R implementation into a small internal R package:

```text
r/chronopackR/
    DESCRIPTION
    NAMESPACE
    R/
        tradeseq.R
        trajectories.R
        clustering.R
        plots.R
        utils.R
```

Then Python calls functions such as:

```r
chronopackR::fit_gene_dynamics(...)
chronopackR::cluster_gene_trajectories(...)
chronopackR::plot_mode_summary(...)
```

Advantages:

- R dependencies are explicitly declared,
- R functions have a proper namespace,
- R code can be tested with `testthat`,
- no uncontrolled `source()` chains,
- easier long-term maintenance,
- clearer Python/R boundary.

This is probably the best long-term architecture if the R component remains substantial.

The internal R package does not necessarily need to be published independently.

---

# 16. Suggested ownership by language

A reasonable initial division:

## Python

```text
PreprocessDatasets
AnnData handling
validation
I/O
workflow orchestration
GRN inference
GRN benchmarking
C2C
cross-dataset result objects
most infrastructure
```

## R

```text
GetTrajectories
tradeSeq preparation internals
tradeSeq fitting
selected trajectory clustering
selected trajectory characterization
selected publication-quality plotting
```

This division can evolve later.

Do not translate stable working R code unless there is a concrete benefit.

---

# 17. Python wrappers should be thin

A backend wrapper should mainly:

1. validate Python-side inputs,
2. convert data into the expected R representation,
3. invoke the R function,
4. retrieve results,
5. normalize outputs,
6. store results in AnnData or a result object,
7. record provenance.

Avoid duplicating statistical logic in both languages.

---

# 18. R functions should have clean interfaces

Refactor R notebook code into functions with explicit arguments.

Avoid functions that rely on:

- global variables,
- notebook state,
- objects created in previous cells,
- hard-coded file paths,
- working-directory assumptions.

Bad:

```r
run_tradeSeq()
# assumes sce, PT, CW, sf, counts are already defined
```

Good:

```r
fit_gene_dynamics <- function(
    counts,
    pseudotime,
    cell_weights,
    size_factors,
    nknots,
    ...
) {
    ...
}
```

Likewise for plotting.

---

# 19. AnnData should remain the primary user-facing state

After R computation, extract reusable outputs and return them to AnnData.

Examples:

```python
adata.var["cp_wald"]
adata.var["cp_qvalue"]
adata.var["cp_mode"]
adata.var["cp_onset"]
adata.var["cp_peak_time"]
adata.var["cp_peak_log2fc"]

adata.varm["cp_fit"]
adata.varm["cp_mode_probs"]

adata.uns["chronopack"]["gene_dynamics"]["time_grid"]
```

Do not require downstream Python code to understand native tradeSeq internals.

---

# 20. Keep full backend state separate from normalized outputs

The architecture should distinguish:

```text
normalized ChronoPack outputs
```

from:

```text
native backend objects
```

Normalized outputs should drive most downstream analysis.

Native backend objects should be retained only for:

- backend-specific diagnostics,
- advanced users,
- exact model re-use,
- specialized R plots that cannot be reconstructed otherwise.

---

# 21. GRN and C2C should remain Python-first

Because these modules are already mostly Python, keep them that way.

Recommended structure:

```python
cp.tl.grn.infer(...)
cp.tl.grn.rank_regulators(...)
cp.tl.grn.benchmark(...)

cp.tl.c2c.infer(...)
cp.tl.c2c.temporal_analysis(...)
```

These should consume standardized AnnData/results produced by upstream steps, not depend on R implementation details.

---

# 22. The full conceptual architecture

Recommended final structure:

```text
                         USER
                          │
                          ▼
                 ChronoPack Python API
                          │
                          ▼
                       AnnData
                          │
        ┌─────────────────┼──────────────────┐
        │                 │                  │
        ▼                 ▼                  ▼
  preprocessing       gene dynamics       GRN / C2C
      Python               │                Python
                           │
                           ▼
                      R backend
                           │
                 ┌─────────┼─────────┐
                 ▼         ▼         ▼
            trajectories tradeSeq clustering
```

Plotting:

```text
                       cp.pl
                         │
             ┌───────────┴───────────┐
             ▼                       ▼
        Python plots              R plots
                                 ggplot2
```

---

# 23. Refactoring order

Recommended sequence:

## Phase 1

Create package skeleton:

```text
cp.io
cp.tl
cp.pl
_backends/r
_core
```

## Phase 2

Refactor `PreprocessDatasets.py` into Python package functions.

## Phase 3

Create R bridge and migrate:

```text
GetTrajectories.R
PrepareForTradeSeq.R
TradeSeqFits.R
```

Expose this through:

```python
cp.tl.get_trajectories(...)
cp.tl.fit_gene_dynamics(...)
```

## Phase 4

Migrate:

```text
DimensionalityReductionClustering.R
```

Expose:

```python
cp.tl.cluster_gene_trajectories(...)
```

## Phase 5

Refactor R plotting functions so they consume standardized data.

Expose them through:

```python
cp.pl...
```

## Phase 6

Migrate Python GRN and C2C notebooks/modules.

## Phase 7

Replace notebook-local function definitions with imports from ChronoPack.

---

# 24. Preserve scientific behavior during migration

Do not change algorithms simply because the code is being reorganized.

For each migrated function:

1. capture existing notebook output,
2. isolate the existing logic,
3. turn it into a callable function,
4. add tests,
5. confirm equivalent results,
6. only then improve or generalize behavior.

Architecture refactoring and methodological changes should be separate tasks.

---

# 25. Avoid unnecessary Python rewrites

Do not rewrite working R components solely to claim that ChronoPack is “pure Python.”

The package can still be Python-centric if:

- installation begins from Python,
- users import a Python package,
- AnnData is the primary object,
- public API is Python,
- workflow orchestration is Python,
- R implementation details are internal.

Uniformity of the **public API** matters much more than uniformity of the implementation language.

---

# 26. Key design rule for Codex

The implementation should follow this principle:

> **ChronoPack is a Python-first package with an internal R statistical and plotting backend.**

Specifically:

- users interact with Python,
- AnnData is the main data contract,
- tradeSeq remains explicitly attributed but hidden behind a backend adapter,
- `PrepareForTradeSeq` becomes internal,
- selected R clustering code can remain R,
- publication-quality R plotting functions should be preserved,
- plotting wrappers should pass standardized data rather than entire notebook state,
- Python GRN and C2C modules should remain native Python,
- large backend-specific objects should remain outside AnnData unless needed,
- the language used internally should not leak into the public API.

---

# 27. Suggested Codex instruction

Use the following as the high-level implementation directive:

> Refactor ChronoPack into a Python-centric package while preserving the existing R implementation where it is scientifically useful. The public API should be exposed through `chronopack.io`, `chronopack.tl`, and `chronopack.pl`, with AnnData as the primary user-facing data structure. Isolate tradeSeq, trajectory-related R functions, R-based clustering, and selected ggplot2 plotting functions behind internal backend adapters. Do not rewrite working R code solely for language uniformity. Refactor R notebook code into explicit reusable functions, preferably evolving toward a small internal `chronopackR` package. Plotting functions should consume standardized ChronoPack outputs rather than notebook-specific objects whenever possible. Keep GRN and C2C primarily Python-based. Preserve scientific behavior during refactoring and add tests before making methodological changes.
