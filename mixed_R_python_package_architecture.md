# Converting Mixed R/Python Jupyter Analyses into a Formal Open-Source Package

## Overview

If you currently have a collection of Jupyter notebooks that depend on both **R** and **Python**, and you already have a working **Conda environment** where all dependencies are resolved, the cleanest way to formalize the project is to separate it into three layers:

1. **Software packages** — declare only the dependencies required by the public API.
2. **Cross-language interface** — decide whether R calls Python, Python calls R, or both.
3. **Reproducible development environment** — use Conda and lockfiles to recreate the environment in which the package was developed and validated.

Your existing Conda environment is an excellent starting point for the third layer, but it should **not itself become the package**.

A good target architecture is a **monorepo containing a proper R package and a proper Python package**, with one language acting as the primary public API.

If most of your current workflow logic is in R, using **R as the primary user-facing API** and calling Python functionality through `reticulate` is usually a clean design.

---

# 1. Recommended Repository Structure

A repository could look like this:

```text
chronopack/
│
├── packages/
│   ├── r/
│   │   ├── DESCRIPTION
│   │   ├── NAMESPACE
│   │   ├── R/
│   │   │   ├── clustering.R
│   │   │   ├── kinetics.R
│   │   │   ├── grn.R
│   │   │   └── python.R
│   │   ├── man/
│   │   └── tests/
│   │       └── testthat/
│   │
│   └── python/
│       ├── pyproject.toml
│       ├── src/
│       │   └── chronopack/
│       │       ├── __init__.py
│       │       ├── clustering.py
│       │       ├── ml.py
│       │       └── grn.py
│       └── tests/
│
├── notebooks/
│   ├── examples/
│   └── analyses/
│
├── env/
│   ├── environment.yml
│   └── conda-lock.yml
│
├── docs/
├── .github/
│   └── workflows/
│
├── LICENSE
└── README.md
```

The key conceptual distinction is:

> `DESCRIPTION` and `pyproject.toml` define what the software needs.  
> `environment.yml` and lockfiles define an environment known to reproduce the analyses.

These are related, but they should not be treated as the same thing.

---

# 2. Make R and Python Real Packages

## R package

R package dependencies should be declared in `DESCRIPTION`.

Example:

```text
Package: chronopack
Version: 0.1.0
Imports:
    data.table,
    dplyr,
    Matrix,
    reticulate
Suggests:
    testthat,
    ggplot2
```

Use:

- `Imports` for required runtime dependencies.
- `Suggests` for optional functionality, tests, documentation, and examples.

Your reusable notebook functions should eventually move into `R/`.

---

## Python package

Use a modern `pyproject.toml`.

Example:

```toml
[build-system]
requires = ["hatchling"]
build-backend = "hatchling.build"

[project]
name = "chronopack"
version = "0.1.0"
requires-python = ">=3.11"
dependencies = [
    "numpy>=2.0",
    "pandas>=2.2",
    "scipy>=1.13",
    "scikit-learn>=1.5",
]

[project.optional-dependencies]
ml = [
    "xgboost>=3.0",
]

dev = [
    "pytest",
]
```

Python implementation code should live under something like:

```text
packages/python/src/chronopack/
```

rather than directly inside notebooks.

Eventually, users could install the Python side independently with:

```bash
pip install chronopack
```

while R users could install the R package separately.

---

# 3. Choose One Primary Direction for R/Python Interaction

Avoid building a codebase where execution constantly alternates between languages:

```text
R → Python → R → Python → R
```

Instead, define coarse-grained boundaries.

For example:

```text
R workflow/API
      │
      ├── kinetics
      ├── visualization
      ├── tradeSeq-related analysis
      │
      └── Python backend
             ├── scikit-learn
             ├── XGBoost
             ├── graph algorithms
             └── specialized ML
```

The public API can remain R-oriented while individual computational tasks are delegated to Python internally.

For example, a user might call:

```r
fit_classifier(...)
```

while the implementation internally calls:

```python
chronopack.ml.fit_classifier(...)
```

The user should not need to care which language performs the computation.

---

# 4. If R Is the Main API, Use `reticulate`

If your current notebooks are primarily R-based, `reticulate` is likely the best interface.

It allows R to import Python packages and convert common objects such as:

- matrices
- vectors
- data frames
- lists
- NumPy arrays

A modern package setup can use `reticulate::py_require()` to declare Python requirements from the R package.

Example:

```r
.onLoad <- function(libname, pkgname) {

    reticulate::py_require(
        packages = c(
            "chronopack>=0.1",
            "numpy>=2.0",
            "scikit-learn>=1.5"
        ),
        python_version = ">=3.11"
    )
}
```

However, a cleaner design is usually to let the Python package declare its own dependencies in `pyproject.toml`.

Then the R side only needs:

```r
.onLoad <- function(libname, pkgname) {

    reticulate::py_require(
        "chronopack>=0.1"
    )
}
```

The Python module can be imported lazily:

```r
.chronopack_py <- NULL

.onLoad <- function(libname, pkgname) {

    reticulate::py_require("chronopack>=0.1")

    .chronopack_py <<-
        reticulate::import(
            "chronopack",
            delay_load = TRUE
        )
}
```

Using delayed loading prevents Python from being initialized every time the R package itself is loaded.

This is much cleaner than forcing users to manually run something like:

```r
reticulate::use_condaenv("my_exact_environment")
```

---

# 5. Keep the Conda Environment, but Change Its Role

Your current state is probably approximately:

```text
Conda environment
    ↓
everything happens to work
```

That is useful as a reference environment, but it often contains far more packages than your actual software requires.

For example:

```text
r-base
r-dplyr
r-ggplot2
r-tradeseq
numpy
scipy
pandas
jupyterlab
ipykernel
pytest
matplotlib
...
hundreds of transitive dependencies
```

Do **not** simply export the entire environment and treat that file as the package dependency specification.

Instead, maintain a deliberately small environment file.

Example:

```yaml
name: chronopack-dev

channels:
  - conda-forge
  - bioconda

dependencies:
  - python=3.12
  - r-base=4.6
  - pip
  - jupyterlab

  # system/infrastructure dependencies if needed
  - graphviz
```

Then let:

```text
DESCRIPTION
```

define R dependencies, and:

```text
pyproject.toml
```

define Python dependencies.

The Conda environment becomes the **development and reproducibility environment**, not the software API contract.

---

# 6. Maintain Both a Loose Environment and a Locked Environment

This distinction is extremely useful.

## `environment.yml`

This should be human-maintained and relatively minimal.

Example:

```yaml
dependencies:
  - python=3.12
  - r-base=4.6
  - pip
```

It says:

> These are the major components the project expects.

---

## `conda-lock.yml`

This should be machine-generated and contain exact resolved versions.

Conceptually:

```text
python 3.12.11 build xxx
numpy 2.2.6 build xxx
r-base 4.6.1 build xxx
...
```

It says:

> This exact environment was tested.

You can generate platform-specific locks with something like:

```bash
conda-lock \
    -f env/environment.yml \
    -p linux-64 \
    -p osx-arm64
```

This gives you:

```text
flexibility                     reproducibility
     │                                │
environment.yml  ──────────→  conda-lock.yml
```

---

# 7. Avoid Exact Version Pins in the Package Metadata

A common mistake is to convert the exact development environment into exact package requirements.

For example, avoid:

```toml
numpy == 2.2.6
scipy == 1.15.3
pandas == 2.3.1
```

unless those exact versions are genuinely required.

Likewise, avoid unnecessarily strict R requirements such as:

```text
Imports:
    dplyr (= 1.1.4)
```

Prefer compatibility ranges:

```toml
dependencies = [
    "numpy>=2.0",
    "scipy>=1.13",
]
```

and:

```text
Imports:
    dplyr (>= 1.1.0)
```

The package metadata defines the **supported compatibility range**.

The lockfile records the **exact versions used for testing**.

Conceptually:

```text
Package metadata
    numpy >= 2.0
       │
       │ resolver
       ↓
Reference lockfile
    numpy == 2.2.6
```

---

# 8. Do Not Add Too Many Competing Environment Managers

`renv` is excellent for R-only reproducibility.

However, if the project already has a serious mixed R/Python Conda environment, adding:

```text
Conda
+ renv
+ pip lockfile
+ another Python virtual environment
```

can become unnecessarily complicated.

Multiple systems may all try to manage the same dependencies.

For a mixed-language scientific package, a clean setup is often:

```text
DESCRIPTION       R package contract
pyproject.toml    Python package contract
environment.yml  development environment
conda-lock.yml    exact reference environment
```

That is usually sufficient.

---

# 9. Use the Existing Conda Environment to Bootstrap the Package

Before modifying your working environment, save snapshots.

## Full environment snapshot

Modern Conda:

```bash
conda export \
    --file env/current-environment.yml
```

Older approach:

```bash
conda env export > env/current-environment.yml
```

This gives you a known-working reference.

---

## Export only explicitly installed packages

A cleaner list can be generated with:

```bash
conda export \
    --from-history \
    --format=environment-yaml \
    > env/from-history.yml
```

This is useful because it excludes many transitive dependencies.

A practical migration strategy is:

1. Save the full working environment.
2. Export explicit dependencies.
3. Identify which packages are actually used by your code.
4. Move R runtime dependencies into `DESCRIPTION`.
5. Move Python runtime dependencies into `pyproject.toml`.
6. Reduce `environment.yml` to infrastructure and development dependencies.

---

# 10. Separate Runtime Dependencies from Development Dependencies

Users do not necessarily need tools such as:

```text
jupyterlab
ipykernel
pytest
ruff
devtools
testthat
roxygen2
pkgdown
```

just to run your package.

## Python

Runtime:

```toml
[project]
dependencies = [
    "numpy",
    "pandas",
]
```

Development/testing:

```toml
[dependency-groups]
test = [
    "pytest",
]

docs = [
    "mkdocs",
]
```

---

## R

Runtime:

```text
Imports:
    Matrix,
    dplyr,
    reticulate
```

Optional/development:

```text
Suggests:
    testthat,
    knitr,
    ggplot2
```

Keeping these separate makes package installation much lighter and more robust.

---

# 11. Be Careful About the R/Python Boundary

Cross-language conversion is easy for small objects, but it can become expensive for large genomics data.

Passing something like:

```text
100 genes × 5 columns
```

between R and Python is trivial.

Passing:

```text
20,000 genes × 100,000 cells
```

back and forth repeatedly is not.

Prefer interfaces based on simple, stable data structures:

```text
matrix
sparse matrix
data.frame
vector
file path
```

Avoid unnecessarily converting complex language-specific objects repeatedly.

For example, avoid pipelines like:

```text
Seurat object
    ↓ Python
AnnData
    ↓ R
SingleCellExperiment
    ↓ Python
...
```

unless that conversion is genuinely required.

Prefer APIs such as:

```r
fit_model(
    expression = X,
    metadata = meta
)
```

or, for very large datasets:

```r
fit_model(
    input = "processed_data.h5ad"
)
```

This makes the language boundary easier to maintain and debug.

---

# 12. Make Notebooks Clients of the Package

This is one of the most important architectural changes.

A notebook should not remain a library of reusable function definitions.

Current pattern:

```text
Notebook A
  function calculate_metrics(...)
  function cluster(...)
  function plot_clusters(...)

Notebook B
  slightly modified calculate_metrics(...)
  slightly modified cluster(...)
```

Target pattern:

```text
              package
          /      |       \
        /        |         \
Notebook A   Notebook B   Notebook C
```

The notebook should mostly contain:

- configuration
- dataset-specific analysis choices
- calls into the package API
- visualization
- interpretation

For example:

```r
library(chronopack)

results <- fit_trajectory_models(...)
clusters <- cluster_trajectories(results)
plot_trajectory_clusters(clusters)
```

rather than dozens of reusable function definitions.

A good rule is:

> A notebook should mostly contain analysis decisions and results, not reusable implementation code.

---

# 13. Convert Existing Notebook Functions into Exported and Internal APIs

Functions that are reusable across analyses should become package functions.

For example, helpers such as:

```r
calculate_trajectory_clustering_metrics()
get_descriptors()
plot_cluster_pair_representative_genes()
make_combined_from_pc_df()
```

can become exported functions if users are expected to call them directly.

Smaller implementation helpers should remain internal.

A useful package design distinction is:

```text
Public API
    stable
    documented
    tested
    user-facing

Internal helpers
    implementation detail
    may change freely
```

Try to keep the public API smaller than the total number of functions in the codebase.

---

# 14. Add Tests While Migrating

Do not wait until the end to add testing.

For each reusable notebook function:

```text
notebook function
      ↓
clean up the API
      ↓
move into package
      ↓
write unit test
      ↓
replace notebook code with package call
```

## R example

```r
test_that("trajectory descriptors are stable", {

    result <- get_descriptors(test_curve)

    expect_equal(
        result$onset,
        42,
        tolerance = 1e-6
    )
})
```

## Python example

```python
def test_cluster_assignment():
    result = cluster_trajectories(X)
    assert result.shape == (100,)
```

For a mixed-language package, you should ideally have three categories of tests:

```text
R unit tests
Python unit tests
R → Python integration tests
```

Integration tests are particularly important because code may work correctly in Python but fail when objects are converted through `reticulate`.

---

# 15. Use Continuous Integration

Eventually, GitHub Actions should test multiple environments.

For example:

```text
Linux
 ├── Python 3.11
 ├── Python 3.12
 └── Python 3.13

macOS
 └── Python 3.12

R
 ├── release
 └── devel
```

Run:

```bash
pytest
```

for Python tests and:

```bash
R CMD check
```

for R package checks.

Also include at least one small cross-language integration test.

A useful CI design is to have two categories of jobs.

## Locked-environment test

```text
conda-lock.yml
      ↓
Can we reproduce the exact validated environment?
```

## Fresh-install compatibility test

```text
DESCRIPTION + pyproject.toml
      ↓
Can a clean user install the currently supported dependencies?
```

Testing only the lockfile can hide missing or overly loose dependency declarations.

---

# 16. Conda Can Eventually Become a Distribution Layer

Once the R and Python packages are stable, Conda can also be used as a user-friendly installation mechanism.

You might eventually distribute:

- Python package through PyPI
- R package through CRAN/Bioconductor or GitHub
- complete environment through Conda/Bioconda/Conda-Forge

Conceptually:

```bash
mamba install -c conda-forge -c bioconda chronopack
```

could install the complete stack.

Meanwhile, Python-only users could still do:

```bash
pip install chronopack
```

and R users could install the R interface separately.

This makes Conda an **installation convenience layer**, rather than the only environment in which the software works.

---

# 17. Recommended Architecture for a Mixed R/Python Scientific Package

A strong long-term architecture is:

```text
                         GitHub repo
                             │
               ┌─────────────┴─────────────┐
               │                           │
          R package                   Python package
        chronopack R                 chronopack Python
               │                           │
        DESCRIPTION / R/            pyproject.toml /
               │                    src/chronopack/
               │                           │
               └────── reticulate ─────────┘
                             │
                             │
                     Public R workflow API
                             │
                     notebooks / users
```

Alongside:

```text
environment.yml
        +
conda-lock.yml
```

as the reproducibility layer.

---

# 18. Suggested Migration Workflow

A practical migration process would be:

## Phase 1 — Freeze the current state

Export the current environment:

```bash
conda env export > env/current-environment.yml
```

and commit it as a temporary reference.

Do not modify your known-working environment until you have a reproducible snapshot.

---

## Phase 2 — Inventory notebook functionality

Classify notebook code into:

```text
Reusable function
One-off analysis
Visualization
Data loading
Configuration
Exploratory code
```

Only reusable code should move into the package.

---

## Phase 3 — Establish package skeletons

Create:

```text
packages/r/
packages/python/
```

and initialize proper R and Python packaging metadata.

---

## Phase 4 — Move functions incrementally

For each reusable function:

1. Copy it into the package.
2. Define clean arguments and return values.
3. Remove reliance on notebook globals.
4. Document it.
5. Write a unit test.
6. Replace the notebook implementation with a package call.

---

## Phase 5 — Define language boundaries

Identify which algorithms naturally belong in:

```text
R
Python
```

Avoid rewriting code solely for language consistency.

Use R for functionality that naturally depends on:

- Bioconductor
- tradeSeq
- SingleCellExperiment
- Signac/Seurat-oriented workflows
- R-native plotting/statistics

Use Python where it naturally fits:

- scikit-learn
- XGBoost
- PyTorch
- network/graph packages
- specialized machine-learning libraries

---

## Phase 6 — Minimize dependency declarations

Populate:

```text
DESCRIPTION
pyproject.toml
```

with only actual runtime dependencies.

Keep development infrastructure in the environment configuration.

---

## Phase 7 — Add reproducible environment locking

Create:

```text
environment.yml
conda-lock.yml
```

for the platforms you intend to support.

---

## Phase 8 — Add CI

Run:

```text
R CMD check
pytest
integration tests
```

on every pull request.

---

## Phase 9 — Convert notebooks into examples/workflows

Your notebooks should become consumers of the package.

Eventually they can serve as:

- tutorials
- reproducible analyses
- documentation
- workflow examples

rather than implementation containers.

---

# 19. The Most Important Principle

Keep this distinction throughout the project:

```text
                 ┌─────────────────────────┐
                 │       PACKAGE API       │
                 │                         │
                 │ DESCRIPTION             │
                 │ pyproject.toml          │
                 │                         │
                 │ broad version ranges    │
                 │ minimal dependencies    │
                 └────────────┬────────────┘
                              │
                              │ validated in
                              ↓
                 ┌─────────────────────────┐
                 │ REFERENCE ENVIRONMENT   │
                 │                         │
                 │ environment.yml         │
                 │ conda-lock.yml          │
                 │                         │
                 │ exact versions          │
                 │ R + Python + system     │
                 └─────────────────────────┘
```

Do **not** turn your existing Conda environment directly into the package.

Instead:

> Use the current Conda environment to discover, validate, and reproduce the dependencies of the package.

That distinction is what separates a reproducible research environment from a maintainable open-source software library.
