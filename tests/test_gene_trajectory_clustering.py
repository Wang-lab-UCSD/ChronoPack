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


def test_cluster_gene_trajectories_returns_joint_features():
    pytest.importorskip("rpy2")
    ad = pytest.importorskip("anndata")

    rng = np.random.default_rng(0)
    t = np.linspace(0.1, 1, 6)
    X = np.vstack([np.sin(np.pi * t * (1 + i % 3)) + rng.normal(0, 0.05, t.size) for i in range(12)])
    traj = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"g{i}_A" for i in range(12)]),
        var=pd.DataFrame({"time": t}, index=[f"t{i}" for i in range(t.size)]),
    )
    result = gc.cluster_gene_trajectories(traj, k=[2, 3], alpha=0.1)

    # Shape (6 timepoints) and derivative (5) features, rows in obs order.
    assert result.features.shape == (12, 11)
    assert result.features.index.tolist() == traj.obs_names.tolist()
    # The centroid distances are computed on these features.
    f = result.features.to_numpy()
    labels = result.assignments[2].to_numpy()
    centroid = f[labels == labels[0]].mean(axis=0)
    assert np.linalg.norm(f[0] - centroid) == pytest.approx(result.centroid_distances[2].iloc[0], rel=1e-9)


def _toy_traj(n=20, n_t=8, seed=0):
    ad = pytest.importorskip("anndata")
    rng = np.random.default_rng(seed)
    t = np.linspace(0.1, 1, n_t)
    X = np.vstack([np.sin(np.pi * t * (1 + i % 4)) + (i % 3) * t + rng.normal(0, 0.1, t.size) for i in range(n)])
    return ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=[f"g{i}_A" for i in range(n)]),
        var=pd.DataFrame({"time": t}, index=[f"t{i}" for i in range(t.size)]),
    )


def _joint_features(traj, alpha):
    """Numpy version of the R joint features: z-scored curves and derivatives, block-scaled, alpha-weighted."""
    _, fit_z, ddt_z = gc.standardize_gene_trajectories(traj)

    def scale(M):
        return M / np.sqrt(M.var(axis=0, ddof=1).sum())

    return np.hstack([np.sqrt(1 - alpha) * scale(fit_z.to_numpy()), np.sqrt(alpha) * scale(ddt_z.to_numpy())])


def test_cluster_gene_trajectories_default_uses_full_joint_features():
    pytest.importorskip("rpy2")
    traj = _toy_traj()
    result = gc.cluster_gene_trajectories(traj, k=[2, 3, 4], alpha=0.3)

    np.testing.assert_allclose(result.features.to_numpy(), _joint_features(traj, 0.3), rtol=1e-6, atol=1e-8)
    assert result.metadata["n_pcs"] is None
    assert result.metadata["n_pcs_variance_explained"] is None


def test_cluster_gene_trajectories_full_rank_pcs_reproduce_default():
    pytest.importorskip("rpy2")
    from scipy.spatial.distance import pdist

    traj = _toy_traj()
    k = [2, 3, 4, 5, 6]
    full = gc.cluster_gene_trajectories(traj, k=k, alpha=0.3)
    n_features = full.features.shape[1]
    pcs = gc.cluster_gene_trajectories(traj, k=k, alpha=0.3, n_pcs=n_features)

    assert pcs.features.shape == (traj.n_obs, n_features)
    assert list(pcs.features.columns) == [f"PC{i}" for i in range(1, n_features + 1)]
    # PCA is a centring + rotation: same distances, same tree, same metrics.
    np.testing.assert_allclose(pdist(pcs.features.to_numpy()), pdist(full.features.to_numpy()), atol=1e-10)
    for kk in k:
        np.testing.assert_array_equal(pcs.assignments[kk].to_numpy(), full.assignments[kk].to_numpy())
        np.testing.assert_allclose(pcs.centroid_distances[kk], full.centroid_distances[kk], atol=1e-10)
    for col in ("WSS", "mean_silhouette"):
        np.testing.assert_allclose(pcs.metrics[col], full.metrics[col], atol=1e-10)
    assert pcs.metadata["n_pcs_variance_explained"] == pytest.approx(1.0)


def test_cluster_gene_trajectories_truncated_pcs():
    pytest.importorskip("rpy2")
    traj = _toy_traj()
    result = gc.cluster_gene_trajectories(traj, k=[2, 3], alpha=0.3, n_pcs=3)

    F = _joint_features(traj, 0.3)
    Fc = F - F.mean(axis=0)
    _, s, Vt = np.linalg.svd(Fc, full_matrices=False)
    scores = Fc @ Vt[:3].T
    # Same subspace up to the sign of each component.
    got = result.features.to_numpy()
    np.testing.assert_allclose(np.abs(got), np.abs(scores), rtol=1e-5, atol=1e-8)
    assert result.metadata["n_pcs"] == 3
    assert result.metadata["n_pcs_variance_explained"] == pytest.approx((s[:3] ** 2).sum() / (s**2).sum(), rel=1e-6)
    # Metrics and centroid distances are computed in PC space.
    labels = result.assignments[2].to_numpy()
    centroid = got[labels == labels[0]].mean(axis=0)
    assert np.linalg.norm(got[0] - centroid) == pytest.approx(result.centroid_distances[2].iloc[0], rel=1e-9)


def test_cluster_gene_trajectories_pcs_fit_on_reference():
    pytest.importorskip("rpy2")
    reference = _toy_traj(n=30)
    subset = reference[reference.obs_names[::2]].copy()
    on_ref = gc.cluster_gene_trajectories(reference, k=[2, 3], alpha=0.3, n_pcs=4)
    on_sub = gc.cluster_gene_trajectories(subset, k=[2, 3], alpha=0.3, n_pcs=4, reference=reference)
    alone = gc.cluster_gene_trajectories(subset, k=[2, 3], alpha=0.3, n_pcs=4)

    # The subset is projected onto the reference PCs, not onto PCs of its own.
    np.testing.assert_allclose(
        on_sub.features.to_numpy(), on_ref.features.loc[subset.obs_names].to_numpy(), rtol=1e-8, atol=1e-10
    )
    assert not np.allclose(np.abs(alone.features.to_numpy()), np.abs(on_sub.features.to_numpy()))
    assert on_sub.metadata["n_pcs_variance_explained"] == pytest.approx(on_ref.metadata["n_pcs_variance_explained"])


@pytest.mark.parametrize("bad", [0, -1, 2.5, True])
def test_cluster_gene_trajectories_rejects_invalid_n_pcs(bad):
    traj = _toy_traj()
    with pytest.raises(ValueError, match="n_pcs"):
        gc.cluster_gene_trajectories(traj, k=[2], alpha=0.3, n_pcs=bad)


def test_cluster_gene_trajectories_rejects_too_many_pcs():
    pytest.importorskip("rpy2")
    traj = _toy_traj()
    with pytest.raises(Exception, match="n_pcs"):
        gc.cluster_gene_trajectories(traj, k=[2], alpha=0.3, n_pcs=traj.n_vars * 2 + 5)

