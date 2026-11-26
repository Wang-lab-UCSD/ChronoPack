import json

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import chronopack as cp
from chronopack._core._validation import ChronoPackError
from chronopack.tl import _gene_dynamics


def _composition_adata(induced_share_end=0.3, n_cells=150, n_genes=400, n_induced=20, seed=0, base_share=0.01):
    """Three time points; the first n_induced genes grow from base_share to induced_share_end of every library, the
    rest are unchanged in transcript numbers (each cell captures a fixed number of the other genes' UMIs)."""
    rng = np.random.default_rng(seed)
    base = rng.gamma(0.6, 1.0, n_genes)
    base /= base.sum()
    rows, times = [], []
    for k, share in enumerate([base_share, (base_share + induced_share_end) / 2, induced_share_end]):
        p_other = base.copy()
        p_other[:n_induced] = 0
        p_other /= p_other.sum()
        p_ind = np.zeros(n_genes)
        p_ind[:n_induced] = 1.0 / n_induced
        for _ in range(n_cells):
            n_other = rng.poisson(2000)
            n_ind = rng.poisson(n_other * share / (1 - share))
            rows.append(rng.multinomial(n_other, p_other) + rng.multinomial(n_ind, p_ind))
            times.append(k + 1)
    X = np.vstack(rows)
    obs = pd.DataFrame({"time": times, "total_counts": X.sum(axis=1).astype(float)},
                       index=[f"c{i}" for i in range(len(rows))])
    return AnnData(X=X.astype(float), obs=obs, var=pd.DataFrame(index=[f"g{j}" for j in range(n_genes)]),
                   layers={"counts": X})


def test_tmm_norm_factors_identical_samples_are_one():
    x = np.random.default_rng(1).poisson(20, size=(1, 300)).astype(float)
    f = cp.tl.tmm_norm_factors(np.vstack([x, 2 * x, 3 * x]))
    np.testing.assert_allclose(f, 1.0, atol=1e-9)


def test_tmm_norm_factors_ignore_an_induced_block():
    rng = np.random.default_rng(2)
    x = rng.poisson(50, size=300).astype(float)
    y = x.copy()
    y[:15] *= 20                                     # 15 genes take a large share of sample 2
    f = cp.tl.tmm_norm_factors(np.vstack([x, y]))
    share = y[:15].sum() / y.sum() - x[:15].sum() / x.sum()
    # effective library ratio = the other genes' ratio: lib_y * f_y / (lib_x * f_x) ~ 1
    ratio = y.sum() * f[1] / (x.sum() * f[0])
    assert share > 0.3
    assert ratio == pytest.approx(1.0, abs=0.03)


def test_composition_qc_measures_the_shift_and_flags_unchanged_genes():
    adata = _composition_adata()
    with pytest.warns(UserWarning, match="Library composition"):
        res = cp.tl.composition_qc(adata)
    bt = res.by_time
    X0 = adata.layers["counts"][adata.obs["time"].to_numpy() == 1]
    s0 = X0[:, :20].sum() / X0.sum()                 # the induced genes' share at baseline
    expected = np.log2(0.7 / (1 - s0))
    assert bt["composition_shift_log2"].iloc[0] == 0
    assert bt["composition_shift_log2"].iloc[-1] == pytest.approx(expected, abs=0.05)
    assert bt["induced_share_bound_log2"].iloc[-1] == pytest.approx(expected, abs=0.05)
    assert res.level["all"] == "strong"
    g = res.genes
    assert set(res.induced_genes["all"]) <= {f"g{j}" for j in range(20)}
    # unchanged genes do not move under TMM
    unchanged = g.iloc[20:][g.iloc[20:]["expressed"]]
    assert unchanged["min_log2fc_tmm"].median() > -0.1
    assert unchanged["min_log2fc_total"].median() < -0.4
    assert "composition_explainable" in adata.var.columns
    assert adata.var["composition_explainable"].iloc[20:].sum() > 0


def test_composition_qc_negligible_without_an_induced_block():
    adata = _composition_adata(induced_share_end=0.01)
    res = cp.tl.composition_qc(adata, warn_threshold=10)
    assert abs(res.max_shift_log2["all"]) < 0.1
    assert res.level["all"] == "negligible"


def test_composition_qc_per_group():
    a, b = _composition_adata(seed=3), _composition_adata(induced_share_end=0.01, seed=4)
    a.obs["cell_type"], b.obs["cell_type"] = "A", "B"
    b.obs_names = [f"b{i}" for i in range(b.n_obs)]
    import anndata as ad

    both = ad.concat([a, b])
    with pytest.warns(UserWarning):
        res = cp.tl.composition_qc(both, groupby="cell_type")
    assert res.level["A"] == "strong" and res.level["B"] == "negligible"
    assert set(res.genes["group"]) == {"A", "B"}


def test_tmm_library_sizes_are_written_and_recorded(tmp_path):
    adata = _composition_adata()
    paths = _gene_dynamics.prepare_tradeseq_inputs(adata, tmp_path, normalization="tmm")
    libs = pd.read_csv(paths["cell_libsize_full"])["libsize_full"].to_numpy()
    expected = adata.obs["total_counts"].to_numpy() * adata.obs["tmm_factor"].to_numpy()
    np.testing.assert_allclose(libs, expected)
    record = json.loads((tmp_path / "library_size.json").read_text())
    assert record["normalization"] == "tmm" and record["groupby"] == ["time"]
    # the last time point's cells get a smaller multiplier (their total over-counts the bulk of genes)
    m = adata.obs.groupby("time")["tmm_factor"].first()
    X0 = adata.layers["counts"][adata.obs["time"].to_numpy() == 1]
    assert m[3] / m[1] == pytest.approx(0.7 / (1 - X0[:, :20].sum() / X0.sum()), abs=0.04)


def test_total_normalization_is_the_default_and_recorded(tmp_path):
    adata = _composition_adata()
    paths = _gene_dynamics.prepare_tradeseq_inputs(adata, tmp_path)
    libs = pd.read_csv(paths["cell_libsize_full"])["libsize_full"].to_numpy()
    np.testing.assert_allclose(libs, adata.obs["total_counts"].to_numpy())
    assert json.loads((tmp_path / "library_size.json").read_text()) == {"normalization": "total"}
    assert "tmm_factor" not in adata.obs


def test_tmm_refuses_prescreened_gene_subsets(tmp_path):
    adata = _composition_adata()
    keep = cp.tl.prescreen_time_association(adata, df=2, p_threshold=0.2)["keep"]
    sub = adata[:, keep.to_numpy()].copy()
    with pytest.raises(ChronoPackError, match="prescreen"):
        _gene_dynamics.prepare_tradeseq_inputs(sub, tmp_path, normalization="tmm")
    # factors computed on all genes first are kept by the subset and used
    cp.tl.tmm_factors(adata)
    sub = adata[:, keep.to_numpy()].copy()
    _gene_dynamics.prepare_tradeseq_inputs(sub, tmp_path, normalization="tmm")


def test_fit_gene_dynamics_refuses_a_normalization_the_inputs_were_not_written_with(tmp_path):
    adata = _composition_adata()
    _gene_dynamics.prepare_tradeseq_inputs(adata, tmp_path)
    with pytest.raises(ChronoPackError, match="normalization"):
        cp.tl.fit_gene_dynamics(input_dir=tmp_path, normalization="tmm")


def _response_adata(mode, n_cells=400, n_genes=400, n_module=15, seed=0):
    """Baseline + three stimulated time points; module genes rise 4x at full response."""
    rng = np.random.default_rng(seed)
    mu0 = rng.gamma(2.0, 0.5, n_genes) + 0.3
    fold = np.ones(n_genes)
    fold[:n_module] = 4.0
    rows, times = [], []
    for k, level in enumerate([0.0, 0.5, 1.0, 1.0]):
        sf = rng.lognormal(0, 0.3, n_cells)
        if mode == "graded":
            mult = 1 + level * (fold - 1)[None, :] * np.ones((n_cells, 1))
        else:                                        # all-or-none: a share `level` of cells at full response
            resp = rng.random(n_cells) < level * 0.6
            mult = np.where(resp[:, None], fold[None, :], 1.0)
        mu = sf[:, None] * mu0[None, :] * mult
        rows.append(rng.negative_binomial(5, 5 / (5 + mu)))
        times += [k + 1] * n_cells
    X = np.vstack(rows)
    obs = pd.DataFrame({"time": times, "total_counts": X.sum(axis=1).astype(float)},
                       index=[f"c{i}" for i in range(len(times))])
    return AnnData(X=X.astype(float), obs=obs, var=pd.DataFrame(index=[f"g{j}" for j in range(n_genes)]),
                   layers={"counts": X})


def test_response_heterogeneity_estimates_the_share_of_non_responders():
    graded = cp.tl.response_heterogeneity(_response_adata("graded"), n_boot=100, n_sim=20)
    aon = cp.tl.response_heterogeneity(_response_adata("aon"), n_boot=100, n_sim=20)
    late_g = graded.by_time[graded.by_time.time >= 3]
    late_a = aon.by_time[aon.by_time.time >= 3]
    assert set(late_g["call"]) == {"all respond"}
    assert set(late_a["call"]) == {"non-responders"}
    # 40% of cells never respond in the all-or-none data
    assert late_a["nonresponder_share"].mean() == pytest.approx(0.4, abs=0.1)
    assert (late_a["nonresponder_lo"] < 0.4).all() and (late_a["nonresponder_hi"] > 0.3).all()
    assert late_g["nonresponder_share"].max() < 0.05
    assert set(graded.module.index) <= {f"g{j}" for j in range(15)}


def test_response_heterogeneity_writes_a_scaled_score_and_groups():
    adata = _response_adata("graded")
    adata.obs["state"] = np.where(np.arange(adata.n_obs) % 2, "x", "y")
    res = cp.tl.response_heterogeneity(adata, groupby="state", n_boot=50, n_sim=10)
    assert set(res.by_time["group"]) == {"all", "x", "y"}
    base = adata.obs["time"] == 1
    assert adata.obs.loc[base, "response_score"].median() == pytest.approx(0, abs=1e-9)
    assert adata.obs.loc[adata.obs["time"] == res.peak_time, "response_score"].median() == pytest.approx(1, abs=1e-9)


def test_response_heterogeneity_accepts_a_gene_module():
    adata = _response_adata("graded")
    res = cp.tl.response_heterogeneity(adata, genes=["g0", "g1", "g2", "g3"], n_boot=20, n_sim=5)
    assert list(res.module.index) == ["g0", "g1", "g2", "g3"]
    with pytest.raises(ChronoPackError):
        cp.tl.response_heterogeneity(adata, genes=["nope"])


def test_plots_render():
    import matplotlib

    matplotlib.use("Agg")
    adata = _composition_adata()
    res = cp.tl.composition_qc(adata, warn_threshold=10)
    assert cp.pl.composition_qc(res) is not None
    radata = _response_adata("aon")
    rres = cp.tl.response_heterogeneity(radata, n_boot=20, n_sim=5)
    assert cp.pl.response_heterogeneity(rres, adata=radata) is not None
