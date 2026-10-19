import numpy as np
import pandas as pd
import pytest
from anndata import AnnData


def _fate_adata(n=120, n_genes=12, seed=0):
    rng = np.random.default_rng(seed)
    time = np.repeat([1, 2, 3, 4], n // 4)
    p1 = rng.uniform(0, 1, n)
    X = rng.poisson(1.0, size=(n, n_genes)).astype(float)
    X[:, 0] = rng.poisson(4 * p1)  # g0 tracks fate 1
    X[:, 1] = rng.poisson(4 * (1 - p1))  # g1 tracks fate 0
    X[:, 2] = 0.0
    X[:2, 2] = 1.0  # g2 rarely detected
    X = np.log1p(X)
    angle = np.deg2rad(30)
    base = np.column_stack([rng.normal(0, 1, n) + 2 * p1, time + rng.normal(0, 0.3, n)])
    R = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    adata = AnnData(
        X=X,
        obs=pd.DataFrame({"time": time, "traj_clusters": np.where(p1 > 0.5, "1", "0")}, index=[f"c{i}" for i in range(n)]),
        var=pd.DataFrame(index=[f"g{j}" for j in range(n_genes)]),
        obsm={"ot_umap": base @ R.T, "terminal_fate_probs_traj": np.column_stack([1 - p1, p1])},
        uns={"terminal_fate_probs_traj_fates": ["0", "1"]},
    )
    adata.layers["log1p"] = X.copy()
    adata.obs["fate_0"] = 1 - p1
    adata.obs["fate_1"] = p1
    return adata


def test_rotate_embedding_to_time_aligns_time_with_y():
    import chronopack as cp

    adata = _fate_adata()
    key = cp.tl.rotate_embedding_to_time(adata, basis="ot_umap")
    assert key == "ot_umap_rot"
    Xr = adata.obsm[key]
    t = adata.obs["time"].to_numpy(float)
    assert np.corrcoef(Xr[:, 1], t)[0, 1] > 0.9
    np.testing.assert_allclose(((t - t.mean())[:, None] * Xr).sum(axis=0)[0], 0.0, atol=1e-8)
    np.testing.assert_allclose(Xr.mean(axis=0), 0.0, atol=1e-12)
    X = adata.obsm["ot_umap"]
    np.testing.assert_allclose(np.linalg.norm(Xr[0] - Xr[1]), np.linalg.norm(X[0] - X[1]))

    cp.tl.rotate_embedding_to_time(adata, basis="ot_umap", key_added="flipped", flip_x=True)
    np.testing.assert_allclose(adata.obsm["flipped"], Xr * [-1, 1])


def test_early_fate_correlations_matches_scipy():
    from scipy.stats import false_discovery_control, spearmanr

    import chronopack as cp

    adata = _fate_adata()
    res = cp.tl.early_fate_correlations(adata, early_bins=(1,), layer="log1p")
    early = adata.obs["time"].to_numpy() == 1
    assert "g2" not in set(res["gene"])  # above the bin mean in 2/30 early cells < 10%
    assert res["r"].is_monotonic_decreasing
    assert res["gene"].iloc[0] == "g0" and res["gene"].iloc[-1] == "g1"

    y = (adata.obs["fate_1"] - adata.obs["fate_0"]).to_numpy()[early]
    by_gene = res.set_index("gene")
    for g in ["g0", "g1", "g5"]:
        rho = spearmanr(adata[early, g].layers["log1p"].ravel(), y).statistic
        assert by_gene.loc[g, "r"] == pytest.approx(rho, abs=1e-10)
    np.testing.assert_allclose(res["qval"], false_discovery_control(res["pval"].to_numpy()), rtol=1e-10)

    with pytest.raises(ValueError, match="Too few cells"):
        cp.tl.early_fate_correlations(adata, early_bins=(99,))


def test_cell_trajectory_figures_render_and_save(tmp_path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    import chronopack as cp

    ctp = cp.pl.cell_trajectories
    adata = _fate_adata()
    cp.tl.rotate_embedding_to_time(adata, basis="ot_umap")

    fig = ctp.time_minutes(adata, embedding_basis="ot_umap_rot", minutes_per_tag=10, output=tmp_path / "time.pdf")
    assert (tmp_path / "time.pdf").stat().st_size > 0
    assert list(adata.obs["Time (min)"].cat.categories) == [30, 20, 10, 0]
    assert adata.obs["Time (min)"].iloc[0] == 0
    assert fig.axes[0].get_xlabel() == "UMAP1"

    labels = {"0": "Inflammatory-hi", "1": "Erythroid-like"}
    ctp.trajectory_cell_types(adata, labels=labels, embedding_basis="ot_umap_rot")
    assert list(adata.obs["traj_cell_type"].cat.categories) == ["Inflammatory-hi", "Erythroid-like"]
    assert [c.upper() for c in adata.uns["traj_cell_type_colors"]] == ["#2C6DA4", "#B84A4A"]

    fig = ctp.fate_probability_panels(
        adata, embedding_basis="ot_umap_rot", which_two=(1, 0), figsize=(8, 8),
        save_panels=True, outdir=tmp_path / "panels",
    )
    assert len(fig.axes) == 8
    assert all(ax.get_title() == "" for ax in fig.axes[:4])
    assert (tmp_path / "panels" / "fate_umap__prob__1.png").exists()
    np.testing.assert_allclose(adata.obs["fate_maxprob"], adata.obsm["terminal_fate_probs_traj"].max(axis=1), rtol=1e-6)

    fig = ctp.two_gene_fate_scatter(adata, ("g0", "g1"), embedding_basis="ot_umap_rot", expression_vmax={"g0": 2})
    assert len(fig.axes) == 6
    early = adata.obs["time"].to_numpy() == 1
    from scipy.stats import spearmanr

    expected = spearmanr(adata.obs["fate_0"][early], adata[early, "g1"].layers["log1p"].ravel()).statistic
    assert fig.spearman_r["g1"] == pytest.approx(expected)
    with pytest.raises(ValueError, match="not found"):
        ctp.two_gene_fate_scatter(adata, ("g0", "missing"), embedding_basis="ot_umap_rot")

    fig = ctp.two_gene_fate_scatter(adata, ("g0", "g1"), embedding_basis="ot_umap_rot", xlim="auto", xticks="auto")
    lo, hi = adata.obs["fate_0"][early].min(), adata.obs["fate_0"][early].max()
    pad = 0.05 * (hi - lo)
    assert fig.axes[0].get_xlim() == pytest.approx((lo - pad, hi + pad))
    ticks = fig.axes[0].get_xticks()
    assert len(ticks[(ticks >= lo - pad) & (ticks <= hi + pad)]) >= 2

    res = cp.tl.early_fate_correlations(adata, layer="log1p")
    fig = ctp.top_fate_correlations(res, top_n=3, q_threshold=1.0, output=tmp_path / "corr.png")
    assert fig.top_genes["pos"][0] == "g0" and fig.top_genes["neg"][0] == "g1"
    assert (tmp_path / "corr.png").exists()
    plt.close("all")


def test_fate_probabilities_keeps_its_layout():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    import chronopack as cp

    fig, axes = cp.pl.cell_trajectories.fate_probabilities(_fate_adata(), size=5)
    assert tuple(fig.get_size_inches()) == (5.0, 4.5)
    assert axes[0, 0].get_title() == "Fate prob: 0"
    plt.close("all")


def test_rotate_embedding_to_time_by_group_keeps_each_group_in_place():
    import chronopack as cp

    adata = _fate_adata()
    shift = np.where(adata.obs["traj_clusters"] == "1", 20.0, 0.0)
    adata.obsm["ot_umap"] = adata.obsm["ot_umap"] + np.column_stack([shift, -shift])
    key = cp.tl.rotate_embedding_to_time(adata, basis="ot_umap", groupby="traj_clusters", key_added="rot_g")
    X, Xr = adata.obsm["ot_umap"], adata.obsm[key]
    t = adata.obs["time"].to_numpy(float)
    for group in ["0", "1"]:
        m = (adata.obs["traj_clusters"] == group).to_numpy()
        np.testing.assert_allclose(Xr[m].mean(axis=0), X[m].mean(axis=0), atol=1e-10)
        assert np.corrcoef(Xr[m, 1], t[m])[0, 1] > 0.8
        tc = t[m] - t[m].mean()
        np.testing.assert_allclose((tc[:, None] * (Xr[m] - Xr[m].mean(axis=0))).sum(axis=0)[0], 0.0, atol=1e-8)
        i, j = np.flatnonzero(m)[:2]
        np.testing.assert_allclose(np.linalg.norm(Xr[i] - Xr[j]), np.linalg.norm(X[i] - X[j]))


def test_trajectory_endpoints_labels_first_and_last_time_per_group(tmp_path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    import chronopack as cp

    adata = _fate_adata()
    adata.obs["cell_type"] = np.where(np.arange(adata.n_obs) % 2, "A", "B")
    out = tmp_path / "endpoints.png"
    fig = cp.pl.cell_trajectories.trajectory_endpoints(
        adata,
        labels={"1": "fate one", "0": "fate zero"},
        palette={"fate one": "#008000", "fate zero": "#E35335"},
        minutes_per_tag=15,
        output=out,
    )
    ax = fig.axes[0]
    assert out.exists()
    assert [t.get_text() for t in ax.get_legend().get_texts()] == ["fate one", "fate zero"]
    assert sorted(t.get_text() for t in ax.texts) == ["0 min", "0 min", "45 min", "45 min"]
    coords = adata.obsm["ot_umap"]
    pts = coords[(adata.obs["time"].to_numpy() == 4) & (adata.obs["cell_type"] == "A").to_numpy()]
    y_range = coords[:, 1].max() - coords[:, 1].min()
    label = ax.texts[3]  # groups in order of appearance (B, A), each first then last time
    assert label.get_text() == "45 min"
    np.testing.assert_allclose(label.get_position(), (np.median(pts[:, 0]), np.median(pts[:, 1]) + 0.015 * y_range))
    assert sum(len(c.get_offsets()) for c in ax.collections) == adata.n_obs
    with pytest.raises(KeyError, match="No colour"):
        cp.pl.cell_trajectories.trajectory_endpoints(adata, labels={"0": "a", "1": "b"}, palette={"a": "red"})
    plt.close("all")
