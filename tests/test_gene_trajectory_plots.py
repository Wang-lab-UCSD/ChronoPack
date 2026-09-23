import numpy as np
import pandas as pd
import pytest


def _synthetic_gene_trajectory_adata():
    ad = pytest.importorskip("anndata")

    rows = ["g0_A", "g1_A", "g2_B", "g3_B"]
    adata = ad.AnnData(
        X=np.array(
            [
                [0.0, 1.0, 2.0],
                [0.0, 0.5, 1.5],
                [2.0, 1.0, 0.0],
                [1.5, 0.5, -0.5],
            ]
        ),
        obs=pd.DataFrame(index=pd.Index(rows, name="row_id")),
        var=pd.DataFrame({"time": [0.0, 10.0, 20.0]}, index=["t0", "t1", "t2"]),
    )
    adata.obs["gene"] = ["g0", "g1", "g2", "g3"]
    adata.obs["cell_type"] = ["A", "A", "B", "B"]
    adata.obs["waldStat_1"] = [9.0, 8.0, 7.0, 6.0]
    adata.obs["peak_log2fc"] = [2.0, 1.2, -1.1, -2.0]
    adata.obs["t_peak"] = [20.0, 30.0, 40.0, 50.0]
    adata.obs["broad_cluster"] = [
        "trend-up",
        "trend-up",
        "trend-down_1",
        "trend-down_1",
    ]
    adata.obs["cluster_fine_induced"] = ["1", "1", "-1", "-1"]
    adata.obs["cluster_fine_repressed"] = ["-1", "-1", "2", "2"]
    adata.obs["PC1"] = [1.0, 0.8, -1.0, -0.9]
    adata.obs["PC2"] = [0.0, 0.1, 0.0, -0.2]
    return adata


def test_gene_trajectory_plot_data_uses_trajectory_anndata():
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


def test_trajectory_features_reads_obs_columns():
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


def test_clustering_metrics_plot_marks_recommended_k():
    import chronopack as cp
    from chronopack.tl import GeneTrajectoryClusteringResult

    result = GeneTrajectoryClusteringResult(
        metrics=pd.DataFrame(
            {
                "k": [3, 4, 5],
                "WSS": [10.0, 7.0, 6.0],
                "p90_centroid_distance": [1.0, 0.8, 0.7],
                "resolution": [0.0, 0.75, 1.0],
                "mean_silhouette": [0.2, 0.5, 0.4],
                "elbow_score": [0.0, 0.25, 0.0],
            }
        ),
        assignments={},
    )

    plot = cp.pl.clustering_metrics(result)
    assert plot.figure is not None
    assert plot.data["k"].tolist() == [3, 4, 5]


def test_missing_plot_key_raises_clear_error():
    import chronopack as cp

    adata = _synthetic_gene_trajectory_adata()
    with pytest.raises(KeyError, match="Trajectory observation column"):
        cp.pl.trajectory_features(adata, color="missing_column")


def test_r_backed_plot_rendering_smoke(tmp_path):
    pytest.importorskip("rpy2")
    import chronopack as cp

    adata = _synthetic_gene_trajectory_adata()
    out = tmp_path / "cluster_trajectories.svg"
    result = cp.pl.cluster_trajectories(adata, normalize=False, output=out)

    assert result.output == str(out)
    assert out.exists()
    assert out.stat().st_size > 0
