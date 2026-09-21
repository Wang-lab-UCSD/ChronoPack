# ChronoPack Scratch Notebooks

This folder is a sandbox copy of the active experiment notebooks plus the
recent manuscript-pipeline `output/`, `intmd/`, and puromycin
`self_contained_outputs/` artifacts.

Raw `input/` directories are symlinked to the original project inputs to avoid
copying large immutable files. The copied `output/` and `intmd/` directories
are local to this scratch folder, so rerunning cells should write here rather
than into the main experiment folders.

## Recommended Environment

```bash
conda activate /new-stg/home/alex/miniconda3/envs/chronopack-manuscript-20260919
export PYTHONPATH=/new-stg/home/alex/tmp/ChronoPack/src:$PYTHONPATH
export MPLCONFIGDIR=/new-stg/home/alex/.tmp/matplotlib-chronopack
```

## Launch Pattern

Start Jupyter from the specific experiment's `notebooks/` directory. The
notebooks intentionally use paths like `../output` and `../intmd`, so the
working directory matters.

```bash
cd /new-stg/home/alex/tmp/ChronoPack/scratch_notebooks/k562-tnfa-puromycin.hg38/notebooks
jupyter lab --no-browser
```

or:

```bash
cd /new-stg/home/alex/tmp/ChronoPack/scratch_notebooks/k562-tnfa.hg38.rep/notebooks
jupyter lab --no-browser
```

or:

```bash
cd /new-stg/home/alex/tmp/ChronoPack/scratch_notebooks/thp1-k562-lps.hg38.rep/notebooks
jupyter lab --no-browser
```

## Notes

- `TradeSeqFits*` notebooks can still be slow; rerun them only intentionally.
- Puromycin `DimensionalityReduction&Clustering.*` scratch notebooks have been
  adjusted to use the scratch experiment roots and scratch
  `self_contained_outputs/`.
- If you want a clean reset, remove this `scratch_notebooks/` directory and
  recreate it from the current verification outputs.
