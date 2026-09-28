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


def test_standardize_gene_trajectories_matches_notebook_math():
    import numpy as np
    import pandas as pd
    from anndata import AnnData

    import chronopack as cp

    rng = np.random.default_rng(0)
    X = rng.normal(size=(4, 6)).astype(np.float32)
    X[2] = 3.0  # constant row keeps SD 1
    traj = AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"g{i}_K562a" for i in range(4)]),
        var=pd.DataFrame({"time": [0.0, 1.0, 2.5, 4.0, 6.0, 9.0]}, index=[f"t{i}" for i in range(6)]),
    )

    fit, fit_z, ddt_z = cp.tl.standardize_gene_trajectories(traj)

    def z_row(M):  # verbatim from DimensionalityReduction&Clustering.self_contained.ipynb
        M = M.astype(float)
        mu = M.mean(axis=1, keepdims=True)
        sd = M.std(axis=1, ddof=1, keepdims=True)
        sd[(~pd.notna(sd)) | (sd == 0)] = 1
        return (M - mu) / sd

    time = traj.var["time"].to_numpy(dtype=float)
    ddt = (traj.X[:, 1:] - traj.X[:, :-1]) / (time[1:] - time[:-1])
    np.testing.assert_array_equal(fit.to_numpy(), traj.X)
    np.testing.assert_array_equal(fit_z.to_numpy(), z_row(traj.X.copy()))
    np.testing.assert_array_equal(ddt_z.to_numpy(), z_row(ddt))
    assert list(ddt_z.columns) == list(traj.var_names[1:])
    assert np.all(fit_z.iloc[2] == 0.0)
