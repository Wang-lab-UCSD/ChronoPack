#!/usr/bin/env python
"""Re-execute K562+TNF PreprocessDatasets up to the first Leiden call and dump
the PCA, neighbour graph and Leiden labels that feed the hard-coded cluster
removal (cell 16: leiden != '2'). Diagnostic only; writes nothing to the run.

Usage: dump_pre_filter_graph.py NOTEBOOK OUT_NPZ
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import numpy as np  # noqa: E402
import scipy.sparse as sp  # noqa: E402

LAST_CELL = 14  # sc.tl.leiden(adata_hvg, resolution=0.5)

notebook = Path(sys.argv[1]).resolve()
out = Path(sys.argv[2]).resolve()
cells = json.loads(notebook.read_text())["cells"]
ns = {"__name__": "__main__", "__file__": str(notebook)}
for idx, cell in enumerate(cells, start=1):
    if idx > LAST_CELL:
        break
    src = "".join(cell.get("source", []))
    if cell.get("cell_type") != "code" or not src.strip():
        continue
    print(f"--- cell {idx} ---", flush=True)
    exec(compile(src, f"{notebook}:cell-{idx}", "exec"), ns)

import scanpy as sc  # noqa: E402

a = ns["adata_hvg"]
conn = sp.csr_matrix(a.obsp["connectivities"])
dist = sp.csr_matrix(a.obsp["distances"])
np.savez_compressed(
    out,
    obs_names=np.asarray(a.obs_names, dtype=str),
    X_pca=a.obsm["X_pca"],
    conn_data=conn.data, conn_indices=conn.indices, conn_indptr=conn.indptr,
    dist_data=dist.data, dist_indices=dist.indices, dist_indptr=dist.indptr,
    leiden=np.asarray(a.obs["leiden"], dtype=str),
    n_obs=a.n_obs,
)
print("scanpy", sc.__version__, "n_obs", a.n_obs, "saved", out)
