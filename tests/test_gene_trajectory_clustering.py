import numpy as np
import pandas as pd
import pytest

from chronopack._core._keys import CHRONOPACK_UNS_KEY
from chronopack.tl import _gene_clustering as gc


def test_characterization_ann_data_contract():
    pc_df = pd.DataFrame(
        {
            "gene": ["g0", "g1"],
            "cell_type": ["A", "A"],
            "t_onset": [1.0, 2.0],
            "t_peak": [3.0, 4.0],
            "peak_log2fc": [1.5, -1.0],
            "owp": [0.8, 0.4],
            "owp_ddt": [0.7, 0.3],
            "mean_curv": [0.1, 0.2],
            "max_curv": [0.4, 0.5],
            "mean_turn": [0.05, 0.06],
            "max_turn": [0.2, 0.3],
        },
        index=["g0_A", "g1_A"],
    )
    fit = pd.DataFrame([[0, 1, 2], [2, 1, 0]], index=pc_df.index, columns=["0", "1", "2"])
    ddt = pd.DataFrame([[1, 1], [-1, -1]], index=pc_df.index, columns=["0.5", "1.5"])

    adata = gc._ensure_gene_adata(None, pc_df.index)
    gc._write_characterization_to_adata(
        adata,
        key_added="traj",
        pc_df=pc_df,
        fit=fit,
        time_grid=np.array([0.0, 1.0, 2.0]),
        ddt=ddt,
    )

    assert adata.shape == (0, 2)
    assert adata.var["traj_peak_log2fc"].tolist() == [1.5, -1.0]
    assert adata.varm["traj_fit"].shape == (2, 3)
    assert adata.varm["traj_fit_ddt"].shape == (2, 2)
    assert adata.uns[CHRONOPACK_UNS_KEY]["traj"]["analysis"] == "gene_expression_trajectory_characterization"


def test_clustering_ann_data_contract():
    pc_df = pd.DataFrame(
        {
            "cluster_broad_induced": ["trend-up", "-1"],
            "cluster_broad_repressed": ["-1", "trend-down_1"],
            "cluster_fine_induced": ["2", "-1"],
            "cluster_fine_repressed": ["-1", "4"],
            "broad_cluster": ["trend-up", "trend-down_1"],
            "cluster_dist_broad": [0.1, 0.2],
            "cluster_dist_fine": [0.3, 0.4],
        },
        index=["g0_A", "g1_A"],
    )
    x_pca = pd.DataFrame([[1, 0], [0, 1]], index=pc_df.index, columns=["PC1", "PC2"])

    adata = gc._ensure_gene_adata(None, pc_df.index)
    gc._write_clustering_to_adata(adata, key_added="modes", pc_df=pc_df, x_pca=x_pca)

    assert adata.var["modes_broad_cluster"].tolist() == ["trend-up", "trend-down_1"]
    assert adata.varm["modes_X_pca"].shape == (2, 2)
    assert adata.uns[CHRONOPACK_UNS_KEY]["modes"]["analysis"] == "gene_expression_trajectory_clustering"


def test_clustering_preserves_existing_standardized_fit():
    pc_df = pd.DataFrame(
        {"broad_cluster": ["trend-up", "trend-down_1"]},
        index=["g0_A", "g1_A"],
    )
    x_pca = pd.DataFrame([[1, 0], [0, 1]], index=pc_df.index, columns=["PC1", "PC2"])
    adata = gc._ensure_gene_adata(None, pc_df.index)
    adata.varm["traj_fit"] = np.array([[0.0, 1.0], [1.0, 0.0]])

    gc._write_clustering_to_adata(adata, key_added="modes", pc_df=pc_df, x_pca=x_pca)

    assert "traj_fit" in adata.varm
    np.testing.assert_array_equal(adata.varm["traj_fit"], np.array([[0.0, 1.0], [1.0, 0.0]]))


def test_r_descriptor_kernel_smoke():
    pytest.importorskip("rpy2")
    from rpy2.robjects import FloatVector, r

    from chronopack._backends.r._bridge import source_r_resource

    source_r_resource("dimred_clustering_helpers.R")
    vals = FloatVector([0, 1, 2, 2, 1, 0])
    mat = r["matrix"](vals, nrow=2, byrow=True)
    mat.rownames = r["c"]("up", "down")
    tt = FloatVector([0, 0.5, 1.0])
    desc = r["get_descriptors"](mat, tt, pseudocount=0)

    assert list(desc.rx2("gene_celltype")) == ["up", "down"]
    assert list(desc.rx2("peak_dlogexp")) == [2.0, -2.0]
