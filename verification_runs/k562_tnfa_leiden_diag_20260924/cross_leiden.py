#!/usr/bin/env python
"""Compare the dumped pre-filter graphs and re-run this env's Leiden on both.

Usage: cross_leiden.py FROZEN_NPZ MAINTAINED_NPZ
Run it in each environment: identical labels across envs for the *same* graph
means Leiden is not the source; differing graphs point at PCA/neighbours.
"""

from __future__ import annotations

import sys

import numpy as np
import scipy.sparse as sp


def load(path):
    z = np.load(path, allow_pickle=False)
    n = int(z["n_obs"])
    conn = sp.csr_matrix((z["conn_data"], z["conn_indices"], z["conn_indptr"]), shape=(n, n))
    dist = sp.csr_matrix((z["dist_data"], z["dist_indices"], z["dist_indptr"]), shape=(n, n))
    return z["obs_names"], z["X_pca"], conn, dist, z["leiden"]


def leiden_on(conn, obs_names):
    import anndata as ad
    import scanpy as sc

    a = ad.AnnData(obs={"_": np.zeros(len(obs_names))})
    a.obs_names = obs_names
    a.obsp["connectivities"] = conn
    a.uns["neighbors"] = {"connectivities_key": "connectivities", "params": {"method": "umap"}}
    sc.tl.leiden(a, resolution=0.5)
    return np.asarray(a.obs["leiden"], dtype=str)


fo, fp, fc, fd, fl = load(sys.argv[1])
mo, mp, mc, md, ml = load(sys.argv[2])
import scanpy as sc  # noqa: E402

print("scanpy", sc.__version__)
print("same obs_names/order:", np.array_equal(fo, mo), len(fo))
print("X_pca max |diff| (abs, sign-insensitive):", float(np.abs(np.abs(fp) - np.abs(mp)).max()))
same_pattern = (fd != 0).astype(int) - (md != 0).astype(int)
print("kNN edges differing:", int(abs(same_pattern).sum()), "of", int((fd != 0).sum()))
print("connectivities max |diff|:", float(abs(fc - mc).max()))
print("stored leiden agree:", float((fl == ml).mean()), "| frozen '2':", int((fl == "2").sum()),
      "maintained '2':", int((ml == "2").sum()), "| same '2' cells:", int(((fl == "2") & (ml == "2")).sum()))
lf = leiden_on(fc, fo)
lm = leiden_on(mc, mo)
print("this-env leiden on frozen graph == frozen stored:", float((lf == fl).mean()))
print("this-env leiden on maintained graph == maintained stored:", float((lm == ml).mean()))
print("this-env leiden on maintained graph == frozen stored:", float((lm == fl).mean()))
print("this-env leiden on frozen graph == maintained stored:", float((lf == ml).mean()))
