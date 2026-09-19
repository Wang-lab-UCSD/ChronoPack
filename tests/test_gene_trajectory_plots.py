import numpy as np
import pandas as pd
import pytest

from chronopack._core._keys import (
    CHRONOPACK_UNS_KEY,
    GENE_TRAJECTORY_CLUSTERING_KEY,
    GENE_TRAJECTORY_KEY,
)


def _synthetic_gene_trajectory_adata():
    ad = pytest.importorskip("anndata")

    rows = ["g0_A", "g1_A", "g2_B", "g3_B"]
    adata = ad.AnnData(
        X=np.zeros((0, len(rows)), dtype=np.float32),
        var=pd.DataFrame(index=pd.Index(rows, name="row_id")),
    )
    fit_key = GENE_TRAJECTORY_KEY
    cluster_key = GENE_TRAJECTORY_CLUSTERING_KEY
    adata.var[f"{fit_key}_gene"] = ["g0", "g1", "g2", "g3"]
    adata.var[f"{fit_key}_cell_type"] = ["A", "A", "B", "B"]
    adata.var[f"{fit_key}_waldStat_1"] = [9.0, 8.0, 7.0, 6.0]
    adata.var[f"{fit_key}_peak_log2fc"] = [2.0, 1.2, -1.1, -2.0]
    adata.var[f"{fit_key}_t_peak"] = [20.0, 30.0, 40.0, 50.0]
    adata.var[f"{cluster_key}_broad_cluster"] = [
        "trend-up",
        "trend-up",
        "trend-down_1",
        "trend-down_1",
    ]
    adata.var[f"{cluster_key}_cluster_fine_induced"] = ["1", "1", "-1", "-1"]
    adata.var[f"{cluster_key}_cluster_fine_repressed"] = ["-1", "-1", "2", "2"]
    adata.varm[f"{fit_key}_fit"] = np.array(
        [
            [0.0, 1.0, 2.0],
            [0.0, 0.5, 1.5],
            [2.0, 1.0, 0.0],
            [1.5, 0.5, -0.5],
        ]
    )
    adata.varm[f"{cluster_key}_X_pca"] = np.array(
        [
            [1.0, 0.0],
            [0.8, 0.1],
            [-1.0, 0.0],
            [-0.9, -0.2],
        ]
    )
    adata.uns[CHRONOPACK_UNS_KEY] = {
        fit_key: {"parameters": {"time_grid": [0.0, 10.0, 20.0]}}
    }
    return adata


def test_gene_trajectory_plot_data_uses_standardized_adata():
    import chronopack as cp

    adata = _synthetic_gene_trajectory_adata()
    result = cp.pl.gene_trajectory(adata, ["g0_A", "g2_B"])

    assert result.output is None
    assert result.data["row_id"].nunique() == 2
    assert sorted(result.data["time"].unique().tolist()) == [0.0, 10.0, 20.0]
    assert result.data.loc[result.data["row_id"] == "g0_A", "value"].tolist() == [0.0, 1.0, 2.0]
    assert set(result.data["broad_cluster"]) == {"trend-up", "trend-down_1"}


def test_cluster_trajectory_plot_data_matches_median_iqr_logic():
    import chronopack as cp

    adata = _synthetic_gene_trajectory_adata()
    result = cp.pl.cluster_trajectories(adata, normalize=False)
    trend_up = result.data[result.data["cluster"] == "trend-up"].sort_values("time")

    assert trend_up["n"].tolist() == [2, 2, 2]
    assert trend_up["median"].tolist() == [0.0, 0.75, 1.75]
    assert trend_up["q25"].tolist() == [0.0, 0.625, 1.625]
    assert trend_up["q75"].tolist() == [0.0, 0.875, 1.875]


def test_trajectory_features_can_read_pca_and_descriptor_columns():
    import chronopack as cp

    adata = _synthetic_gene_trajectory_adata()
    result = cp.pl.trajectory_features(adata, x="PC1", y="PC2", color="broad_cluster")

    assert result.data["PC1"].tolist() == [1.0, 0.8, -1.0, -0.9]
    assert result.data["PC2"].tolist() == [0.0, 0.1, 0.0, -0.2]
    assert result.data["broad_cluster"].tolist() == [
        "trend-up",
        "trend-up",
        "trend-down_1",
        "trend-down_1",
    ]


def test_missing_plot_key_raises_clear_error():
    import chronopack as cp

    adata = _synthetic_gene_trajectory_adata()
    del adata.varm[f"{GENE_TRAJECTORY_KEY}_fit"]

    with pytest.raises(KeyError, match="standardized fitted trajectory matrix"):
        cp.pl.gene_trajectory(adata, ["g0_A"])


def test_r_backed_plot_rendering_smoke(tmp_path):
    pytest.importorskip("rpy2")
    import chronopack as cp

    adata = _synthetic_gene_trajectory_adata()
    out = tmp_path / "cluster_trajectories.svg"
    result = cp.pl.cluster_trajectories(adata, normalize=False, output=out)

    assert result.output == str(out)
    assert out.exists()
    assert out.stat().st_size > 0
