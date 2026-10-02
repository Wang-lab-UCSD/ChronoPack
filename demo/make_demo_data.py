#!/usr/bin/env python
"""Build the K562 + TNF-alpha demo inputs from the PreprocessDatasets outputs.

    make_demo_data.py EXP_OUTPUT_DIR NEGATIVE_CONTROL_PKL DEST_DIR

EXP_OUTPUT_DIR is the ``output/`` folder written by
``analyses/k562-tnfa.hg38.rep/01_PreprocessDatasets.ipynb``. DEST_DIR receives
the three files read by ``K562_TNFa_demo.ipynb``:

- ``k562_tnfa.hvg.h5ad``: the 3,000 highly variable genes, X = expression
  corrected for QC covariates, cell cycle and sample (``layers["corrected"]``),
  plus the preprocessing UMAP. Input for the cell-state trajectories.
- ``k562_tnfa.counts.h5ad``: all genes, ``layers["counts"]`` = raw UMI counts
  (no X). Input for the gene-dynamics fits.
- ``k562_neg-ctrl_mtx-agg.pkl``: the negative-control pseudobulk used to flag
  time-tag-biased genes (copied as is).

Only what the demo uses is kept, so the files are much smaller than the
pipeline objects.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

OBS_COLUMNS = ["sample", "time", "total_counts", "n_genes_by_counts", "pct_counts_mt", "phase"]
VAR_COLUMNS = ["highly_variable", "dispersions_norm"]


def _obs(adata):
    obs = adata.obs[OBS_COLUMNS].copy()
    obs["time"] = pd.to_numeric(obs["time"].astype(str)).astype(int)
    return obs


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit(__doc__)
    exp_output, negative_control, dest = (Path(p) for p in sys.argv[1:])
    dest.mkdir(parents=True, exist_ok=True)

    hvg = ad.read_h5ad(exp_output / "adata_hvg.1e4-variant.chk0.h5ad")
    demo_hvg = ad.AnnData(
        X=np.asarray(hvg.layers["corrected"], dtype=np.float32),
        obs=_obs(hvg),
        var=hvg.var[VAR_COLUMNS].copy(),
        obsm={"X_umap": hvg.obsm["X_umap"].astype(np.float32)},
    )
    demo_hvg.write_h5ad(dest / "k562_tnfa.hvg.h5ad", compression="gzip")

    full = ad.read_h5ad(exp_output / "adata_full.1e4-variant.chk0.h5ad")
    assert (full.obs_names == hvg.obs_names).all()
    counts = sp.csr_matrix(full.layers["counts"])
    assert np.array_equal(counts.data, np.round(counts.data))
    demo_full = ad.AnnData(
        X=None,
        obs=_obs(full),
        var=full.var[VAR_COLUMNS].copy(),
        layers={"counts": counts.astype(np.int32)},
    )
    demo_full.write_h5ad(dest / "k562_tnfa.counts.h5ad", compression="gzip")

    shutil.copyfile(negative_control, dest / "k562_neg-ctrl_mtx-agg.pkl")
    for path in sorted(dest.iterdir()):
        print(f"{path.name}\t{path.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
