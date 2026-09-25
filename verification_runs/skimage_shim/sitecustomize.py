"""Verification-only compatibility hooks for optional Scanpy dependencies."""

from __future__ import annotations

from importlib.util import find_spec


def _install_scanpy_leiden_fallback() -> None:
    if find_spec("igraph") is not None and find_spec("leidenalg") is not None:
        return

    try:
        import numpy as np
        import pandas as pd
        import scipy.sparse as sp
        from natsort import natsorted
        from sklearn.cluster import SpectralClustering

        import scanpy as sc
        from scanpy import _utils
    except Exception:
        return

    def _fallback_leiden(
        adata,
        resolution=1,
        *,
        restrict_to=None,
        random_state=0,
        key_added="leiden",
        adjacency=None,
        directed=None,
        use_weights=True,
        n_iterations=-1,
        partition_type=None,
        neighbors_key=None,
        obsp=None,
        copy=False,
        flavor="leidenalg",
        **clustering_args,
    ):
        if restrict_to is not None:
            raise ImportError("Fallback Leiden shim does not implement restrict_to.")
        out = adata.copy() if copy else adata
        if adjacency is None:
            adjacency = _utils._choose_graph(out, obsp, neighbors_key)
        adjacency = sp.csr_matrix(adjacency)
        n_obs = out.n_obs
        if n_obs <= 1:
            groups = np.zeros(n_obs, dtype=int)
        else:
            res = 1.0 if resolution is None else float(resolution)
            n_clusters = max(2, int(round(max(res, 0.01) * np.sqrt(n_obs) / 12.0)))
            n_clusters = min(n_clusters, n_obs)
            model = SpectralClustering(
                n_clusters=n_clusters,
                affinity="precomputed",
                assign_labels="kmeans",
                random_state=random_state,
            )
            groups = model.fit_predict(adjacency.maximum(adjacency.T))
        out.obs[key_added] = pd.Categorical(
            values=groups.astype("U"),
            categories=natsorted(map(str, np.unique(groups))),
        )
        out.uns[key_added] = {
            "params": {
                "resolution": resolution,
                "random_state": random_state,
                "n_iterations": n_iterations,
                "fallback": "spectral_clustering",
            }
        }
        return out if copy else None

    sc.tl.leiden = _fallback_leiden
    try:
        import scanpy.tools as tools

        tools.leiden = _fallback_leiden
    except Exception:
        pass


_install_scanpy_leiden_fallback()
