import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

import chronopack as cp
from chronopack.tl._trajectories import (
    fate_umap_2x2_scanpy,
    mass_retention_plots_and_umap,
)


ROOT = Path(__file__).resolve().parents[2]
REF_DIR = ROOT / "verification_runs/upstream_preprocess_gettraj_20260917_132244/k562-tnfa-puromycin.hg38"
OUT_DIR = Path(__file__).resolve().parent


def _matrix_summary_equal(a, b):
    if sp.issparse(a) or sp.issparse(b):
        a = sp.csr_matrix(a)
        b = sp.csr_matrix(b)
        diff = a - b
        max_abs = float(np.max(np.abs(diff.data))) if diff.nnz else 0.0
        return {
            "shape_equal": a.shape == b.shape,
            "nnz_a": int(a.nnz),
            "nnz_b": int(b.nnz),
            "nnz_diff": int(diff.nnz),
            "max_abs_diff": max_abs,
        }
    a = np.asarray(a)
    b = np.asarray(b)
    return {
        "shape_equal": a.shape == b.shape,
        "max_abs_diff": float(np.max(np.abs(a - b))) if a.size else 0.0,
        "allclose": bool(np.allclose(a, b, rtol=1e-10, atol=1e-12)),
    }


def _categorical_equal(a, b):
    return a.astype(str).to_numpy().tolist() == b.astype(str).to_numpy().tolist()


def main():
    bck = sc.read_h5ad(REF_DIR / "output/adata_hvg.chk0.h5ad")
    adata_hvg = bck.copy()
    adata_hvg.obs["time"] = pd.to_numeric(adata_hvg.obs["time"].astype(str))

    sc.pp.combat(adata_hvg, key="sample")
    sc.pp.regress_out(adata_hvg, keys=["total_counts", "pct_counts_mt", "S_score", "G2M_score"])
    adata_hvg.layers["corrected"] = adata_hvg.X.copy()
    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=50)
    sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata_hvg)

    adata_hvg = adata_hvg.copy()
    adata_hvg.X = adata_hvg.layers["corrected"]

    wot_out = cp.tl.get_trajectories(
        adata_hvg,
        tmap_dir=REF_DIR / "notebooks/tmaps_sc",
        ot_model_kwargs=dict(epsilon=0.05, lambda1=1, lambda2=50, growth_iters=3, day_field="time"),
        compute_transport_maps=False,
        transition_topk=30,
        symmetrize="average",
        use_unpruned_as_pruned=True,
        graph_key="ot_graph",
        umap_key="X_ot_umap",
        umap_random_state=10,
    )

    T_pruned = wot_out["T_pruned"]
    adata_hvg.obsp["t_pruned"] = T_pruned
    key_out = "X_ot_umap"
    adata_hvg.obsm["X_umap"][:, 0] *= -1
    adata_hvg.obsm[key_out] = adata_hvg.obsm["X_umap"].copy()

    cp.tl.cluster_cell_state_trajectories(
        adata_hvg,
        T_pruned,
        cluster_key="traj_clusters",
        resolution=0.5,
        terminal_topk=200,
        graph_topk=500,
        random_state=0,
    )
    cp.tl.compute_terminal_fate_probabilities(
        adata_hvg,
        terminal_cluster_key="traj_clusters",
        out_key="terminal_fate_probs_traj",
    )
    fate_umap_2x2_scanpy(
        adata_hvg,
        fate_key="terminal_fate_probs_traj",
        fate_names_key="terminal_fate_probs_traj_fates",
        umap_basis=key_out,
        size=25,
    )
    mass_retention_plots_and_umap(
        adata_hvg,
        T_pruned,
        cluster_key="traj_clusters",
        n_reps=50,
        clusters_to_plot=["0", "1"],
        umap_basis=key_out,
        size=25,
        show=False,
    )

    ref = sc.read_h5ad(REF_DIR / "output/adata_hvg.h5ad")

    summary = {
        "input_reference": str(REF_DIR),
        "shape_new": list(adata_hvg.shape),
        "shape_reference": list(ref.shape),
        "obs_names_equal": adata_hvg.obs_names.equals(ref.obs_names),
        "var_names_equal": adata_hvg.var_names.equals(ref.var_names),
        "x": _matrix_summary_equal(adata_hvg.X, ref.X),
        "layers": {
            key: _matrix_summary_equal(adata_hvg.layers[key], ref.layers[key])
            for key in sorted(set(adata_hvg.layers).intersection(ref.layers))
        },
        "obs_columns": {},
        "obsm": {},
        "obsp": {},
    }
    for key in [
        "traj_clusters",
        "terminal_fate_probs_traj_0",
        "terminal_fate_probs_traj_1",
        "fate_0",
        "fate_1",
        "fate_entropy",
        "fate_maxprob",
        "traj_clusters_mass_ret_cell",
    ]:
        if key not in adata_hvg.obs or key not in ref.obs:
            summary["obs_columns"][key] = {"present_new": key in adata_hvg.obs, "present_reference": key in ref.obs}
        elif pd.api.types.is_numeric_dtype(adata_hvg.obs[key]) and pd.api.types.is_numeric_dtype(ref.obs[key]):
            summary["obs_columns"][key] = {
                "max_abs_diff": float(np.nanmax(np.abs(adata_hvg.obs[key].to_numpy() - ref.obs[key].to_numpy()))),
                "allclose": bool(np.allclose(adata_hvg.obs[key].to_numpy(), ref.obs[key].to_numpy(), equal_nan=True, rtol=1e-10, atol=1e-12)),
            }
        else:
            summary["obs_columns"][key] = {"equal": _categorical_equal(adata_hvg.obs[key], ref.obs[key])}

    for key in ["X_umap", "X_ot_umap", "terminal_cell_dist", "terminal_fate_probs_traj"]:
        if key in adata_hvg.obsm and key in ref.obsm:
            summary["obsm"][key] = _matrix_summary_equal(adata_hvg.obsm[key], ref.obsm[key])
        else:
            summary["obsm"][key] = {"present_new": key in adata_hvg.obsm, "present_reference": key in ref.obsm}

    for key in ["ot_graph", "t_pruned", "W_terminal_cosine_topk"]:
        if key in adata_hvg.obsp and key in ref.obsp:
            summary["obsp"][key] = _matrix_summary_equal(adata_hvg.obsp[key], ref.obsp[key])
        else:
            summary["obsp"][key] = {"present_new": key in adata_hvg.obsp, "present_reference": key in ref.obsp}

    if "traj_clusters" in adata_hvg.obs and "traj_leiden_sr" in ref.obs:
        summary["legacy_alias_comparisons"] = {
            "traj_clusters_vs_traj_leiden_sr": {
                "equal": _categorical_equal(adata_hvg.obs["traj_clusters"], ref.obs["traj_leiden_sr"])
            }
        }
        if "terminal_fate_probs_traj" in adata_hvg.obsm and "fate_probs_traj_last" in ref.obsm:
            summary["legacy_alias_comparisons"]["terminal_fate_probs_traj_vs_fate_probs_traj_last"] = _matrix_summary_equal(
                adata_hvg.obsm["terminal_fate_probs_traj"],
                ref.obsm["fate_probs_traj_last"],
            )

    (OUT_DIR / "puromycin_get_trajectories_wot2_comparison_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True)
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
