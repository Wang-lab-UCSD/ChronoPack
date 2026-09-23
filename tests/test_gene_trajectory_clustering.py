import numpy as np
import pandas as pd
import pytest

from chronopack._core._keys import CHRONOPACK_UNS_KEY
from chronopack.tl import _gene_clustering as gc


def test_make_traj_adata_contract():
    pytest.importorskip("anndata")

    fit = pd.DataFrame(
        [[0.0, 1.0, 2.0], [2.0, 1.0, 0.0]],
        index=["g0_A", "g1_A"],
        columns=["t0", "t1", "t2"],
    )
    obs = pd.DataFrame(
        {
            "gene": ["g0", "g1"],
            "cell_type": ["A", "A"],
            "waldStat_1": [10.0, 9.0],
            "meanLogFC": [2.0, -2.0],
            "pvalue_1": [0.001, 0.002],
            "qval_1": [0.01, 0.02],
        },
        index=fit.index,
    )

    traj = gc._make_traj_adata(
        fit=fit,
        obs=obs,
        time_grid=np.array([0.0, 0.5, 1.0]),
        raw_payload={"A": {"fit": fit.to_numpy(), "time": np.array([0.0, 0.5, 1.0])}},
        raw_size_bytes=fit.to_numpy().nbytes,
        raw_embedded=True,
        metadata={"time_map": {"A": 1.0}},
        key_added="traj",
    )

    assert traj.shape == (2, 3)
    np.testing.assert_array_equal(traj.X, fit.to_numpy())
    assert traj.obs["gene"].tolist() == ["g0", "g1"]
    assert traj.var["time"].tolist() == [0.0, 0.5, 1.0]
    assert traj.uns[CHRONOPACK_UNS_KEY]["alignment"]["raw_trajectories_embedded"] is True
    assert "raw_trajectories" in traj.uns[CHRONOPACK_UNS_KEY]


def test_clustering_result_assignments_are_indexed_series():
    obs_names = pd.Index(["g0_A", "g1_A", "g2_B"])
    result = gc.GeneTrajectoryClusteringResult(
        metrics=pd.DataFrame({"k": [2, 3], "mean_silhouette": [0.2, 0.3]}),
        assignments={
            2: pd.Series([1, 1, 2], index=obs_names, name="k2"),
            3: pd.Series([1, 2, 3], index=obs_names, name="k3"),
        },
    )

    assert isinstance(result.assignments[2], pd.Series)
    assert result.assignments[2].index.equals(obs_names)
    assert result.assignments[3].tolist() == [1, 2, 3]
