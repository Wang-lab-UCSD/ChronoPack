#!/usr/bin/env python
"""Benchmark the tradeSeq input writer: legacy float CSV vs integer CSV.

Writes both variants from the same frozen r6 K562+TNF objects, using the same
two calls fit_gene_dynamics(prepare_inputs=True) makes, then checks that the
parsed contents are identical. Inputs are read-only.

Usage: bench_write.py OUT_DIR
"""

from __future__ import annotations

import gzip
import json
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, "/new-stg/home/alex/tmp/ChronoPack/src")
from chronopack.tl import _gene_dynamics as gd  # noqa: E402

print("chronopack from:", gd.__file__, flush=True)

RUN = Path("/new-stg/home/alex/tmp/ChronoPack/verification_runs/full_pipeline_frozen_manuscript_20260924_r6/k562-tnfa.hg38.rep")
PIPELINE_DIR = RUN / "intmd" / "sc-K562-1e4_to_tradeSeq.fixed"

out_root = Path(sys.argv[1]).resolve()
env = sys.prefix

t0 = time.perf_counter()
hvg = ad.read_h5ad(RUN / "output" / "adata_prefiltered.1e4-variant.fixed.h5ad")
full = ad.read_h5ad(RUN / "output" / "adata_full.1e4-variant.chk1.h5ad")
load_s = time.perf_counter() - t0
print(f"env={env} load={load_s:.1f}s hvg={hvg.shape} full={full.shape}", flush=True)


def write(variant):
    out = out_root / variant
    out.mkdir(parents=True, exist_ok=True)
    timings = {}
    t = time.perf_counter()
    gd.prepare_tradeseq_inputs(hvg, out, full_adata=full, tag="")
    timings["subset"] = time.perf_counter() - t
    t = time.perf_counter()
    gd.prepare_tradeseq_inputs(full, out, full_adata=full, tag="full")
    timings["full"] = time.perf_counter() - t
    print(variant, {k: round(v, 1) for k, v in timings.items()}, flush=True)
    return out, timings


# Legacy path: identical code with the integer cast disabled.
_orig = gd._is_integer_valued
gd._is_integer_valued = lambda matrix: False
legacy_dir, legacy_t = write("legacy_float")
gd._is_integer_valued = _orig
new_dir, new_t = write("integer")


def read_counts(path):
    return pd.read_csv(path, index_col=0)


checks = {}
for name in ["counts_genes_x_cells.csv.gz", "counts_genes_x_cells.full.csv.gz"]:
    a = read_counts(legacy_dir / name)
    b = read_counts(new_dir / name)
    checks[name] = {
        "same_index": a.index.equals(b.index),
        "same_columns": a.columns.equals(b.columns),
        "values_equal": bool(np.array_equal(a.to_numpy(float), b.to_numpy(float))),
        "legacy_dtype": str(a.dtypes.iloc[0]),
        "integer_dtype": str(b.dtypes.iloc[0]),
        "size_bytes": {"legacy": (legacy_dir / name).stat().st_size, "integer": (new_dir / name).stat().st_size},
    }
    pipe = PIPELINE_DIR / name
    if pipe.exists():
        p = read_counts(pipe)
        checks[name]["integer_equals_pipeline_output"] = bool(
            p.index.equals(b.index) and p.columns.equals(b.columns)
            and np.array_equal(p.to_numpy(float), b.to_numpy(float))
        )

for name in ["cell_times.csv.gz", "cell_times.full.csv.gz", "subset_genes.csv", "subset_genes.full.csv", "cell_libsize_full.csv"]:
    opener = gzip.open if name.endswith(".gz") else open
    with opener(legacy_dir / name, "rb") as fa, opener(new_dir / name, "rb") as fb:
        checks[name] = {"identical_content": fa.read() == fb.read()}

summary = {
    "env": env,
    "pandas": pd.__version__,
    "load_s": load_s,
    "legacy_float_s": legacy_t,
    "integer_s": new_t,
    "checks": checks,
}
(out_root / "summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
