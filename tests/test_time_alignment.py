import anndata as ad
import numpy as np
import pandas as pd
import pytest

import chronopack as cp
from chronopack._core._validation import ChronoPackError


def _context(scale, offset, *, step, n_tags=12, n_cells=300, base=0.5, fc=20.0, repressed=False, seed=0,
             gene="R", mid=20.0, width=5.0, pulse=False):
    """Cells whose reporter follows a shared logistic programme on the reference clock, read at scale * t + offset
    (``pulse``: the rise peaks and falls back towards 40% of its height, like NFKBIA)."""
    rng = np.random.default_rng(seed)
    t = np.repeat(np.arange(n_tags) * step, n_cells)
    tau = scale * t + offset
    g = 1 / (1 + np.exp(-(tau - mid) / width))
    if pulse:
        g = g * (0.4 + 0.6 * np.exp(-np.clip(tau - mid - 2 * width, 0, None) / 30.0))
    mu = base * (1 - (1 - 1 / fc) * g) if repressed else base * (1 + (fc - 1) * g)
    sf = np.exp(rng.normal(0, 0.3, t.size))
    sf /= sf.mean()
    counts = rng.poisson(rng.gamma(3.0, sf * mu / 3.0))
    X = np.column_stack([counts, rng.poisson(2.0 * sf)]).astype(float)
    obs = pd.DataFrame({"time": t, "total_counts": 4000 * sf}, index=[f"c{i}" for i in range(t.size)])
    a = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[gene, "bg"]))
    a.layers["counts"] = a.X.copy()
    return a


@pytest.fixture(scope="module")
def pair():
    return {"ref": _context(1.0, 0.0, step=10, n_cells=500, seed=1),
            "ctx": _context(0.6, -10.0, step=15, n_cells=300, seed=2)}


@pytest.mark.parametrize("fit_scale", ["linear", "log"])
def test_recovers_known_alignment(pair, fit_scale):
    res = cp.tl.estimate_time_alignment(pair, ["R"], fit_scale=fit_scale, n_boot=150, random_state=0)
    a = res.alignment.set_index("context")
    assert a.loc["ref", "scale"] == 1.0 and a.loc["ref", "offset"] == 0.0
    r = a.loc["ctx"]
    assert r.scale_lo <= 0.6 <= r.scale_hi
    assert abs(np.log(r.scale / 0.6)) < 0.2
    assert r.offset_lo <= -10 <= r.offset_hi
    assert r.method == "registration" and r.reporters_used == "R"
    assert res.fits["resolved"].all()
    assert res.scale_map["ctx"] == pytest.approx(r.scale) and res.offset_map["ref"] == 0.0
    # the logistic rise-width estimate (secondary) is reported and, for a logistic-shaped rise, agrees
    assert r.rise_scale_lo <= 0.6 <= r.rise_scale_hi
    assert set(res.tags["context"]) == {"ref", "ctx"} and res.tags["in_window"].any()
    assert res.pairs.set_index("gene").loc["R", "used"]


def test_landmark_shifts_and_options(pair):
    res = cp.tl.estimate_time_alignment(pair, ["R"], n_boot=100, random_state=0)
    p = res.pairs.set_index("gene").loc["R"]
    # after alignment the context's onset and midpoint sit on the reference's (a logistic has no peak)
    assert p.onset_shift_lo <= 0 <= p.onset_shift_hi and p.mid_shift_lo <= 0 <= p.mid_shift_hi
    assert np.isnan(p.peak_shift)
    assert res.parameters["symmetric"] and res.parameters["interp_sd"] == 0.125
    one_way = cp.tl.estimate_time_alignment(pair, ["R"], symmetric=False, n_boot=0)
    assert abs(np.log(one_way.alignment.set_index("context").loc["ctx", "scale"] / 0.6)) < 0.2


def test_falling_context_tags_past_the_reference_window():
    """K562b-like: same clock, 15-min tags, the rise covered by ~2 tags and the context's window ending on the decline.
    The reference's whole per-tag curve is the template (not its window held flat), so the scale stays near 1."""
    data = {"ref": _context(1.0, 0.0, step=10, n_cells=500, seed=11, pulse=True),
            "ctx": _context(1.0, -65.0, step=15, n_cells=250, seed=12, pulse=True)}
    res = cp.tl.estimate_time_alignment(data, ["R"], n_boot=100, random_state=0)
    r = res.alignment.set_index("context").loc["ctx"]
    assert r.scale_lo <= 1.0 <= r.scale_hi and abs(np.log(r.scale)) < 0.25
    assert abs(r.offset - (-65.0)) < 15
    p = res.pairs.set_index("gene").loc["R"]
    assert p.onset_shift_lo <= 0 <= p.onset_shift_hi and np.isfinite(p.peak_shift)


def test_fixed_scale_and_current(pair):
    res = cp.tl.estimate_time_alignment(pair, "R", fix_scale={"ctx": 0.6}, current={"ctx": (0.9, -30.0)}, n_boot=150)
    r = res.alignment.set_index("context").loc["ctx"]
    assert r.method == "fixed scale" and r.scale == 0.6
    assert r.offset_lo <= -10 <= r.offset_hi
    assert r.current_scale == 0.9 and np.isfinite(r.offset_at_current_scale)
    assert "current_offset_outside_ci" in r["flags"]
    res2 = cp.tl.estimate_time_alignment(pair, "R", current={"ctx": (1.5, -10.0)}, n_boot=150)
    assert "current_scale_outside_ci" in res2.alignment.set_index("context").loc["ctx", "flags"]


def test_repressed_reporter():
    data = {"ref": _context(1.0, 0.0, step=10, n_cells=500, repressed=True, base=5.0, seed=3),
            "ctx": _context(0.6, -10.0, step=15, n_cells=300, repressed=True, base=5.0, seed=4)}
    res = cp.tl.estimate_time_alignment(data, ["R"], n_boot=150)
    assert (res.fits["direction"] == "repressed").all() and (res.fits["a"] < 0).all()
    r = res.alignment.set_index("context").loc["ctx"]
    assert r.scale_lo <= 0.6 <= r.scale_hi


def test_rise_between_two_tags_is_flagged():
    data = {"ref": _context(1.0, 0.0, step=10, n_cells=500, seed=5),
            "ctx": _context(1.0, 0.0, step=30, n_tags=8, n_cells=300, seed=6, width=1.0)}
    res = cp.tl.estimate_time_alignment(data, ["R"], n_boot=100)
    f = res.fits.set_index("context").loc["ctx"]
    assert not f.resolved and "width_at_lower_bound" in f["flags"] and f.n_tags >= 5
    r = res.alignment.set_index("context").loc["ctx"]
    assert "context_rise_undersampled" in r["flags"] and np.isnan(r.rise_scale)


def test_no_rise_gives_no_estimate():
    data = {"ref": _context(1.0, 0.0, step=10, n_cells=500, seed=7),
            "ctx": _context(1.0, 0.0, step=15, n_cells=300, seed=8, fc=1.0)}
    res = cp.tl.estimate_time_alignment(data, ["R"], n_boot=50)
    r = res.alignment.set_index("context").loc["ctx"]
    assert "no_reporter_with_rise" in r["flags"] and np.isnan(r.scale) and res.scale_map["ctx"] is None


def test_input_errors(pair):
    with pytest.raises(ChronoPackError, match="not found"):
        cp.tl.estimate_time_alignment(pair, ["MISSING"], n_boot=0)
    cont = pair["ctx"].copy()
    cont.obs["time"] = np.random.default_rng(0).uniform(0, 100, cont.n_obs)
    with pytest.raises(ChronoPackError, match="discrete sampling times"):
        cp.tl.estimate_time_alignment({"ref": pair["ref"], "ctx": cont}, ["R"], n_boot=0)
    with pytest.raises(ChronoPackError, match="fit_scale"):
        cp.tl.estimate_time_alignment(pair, ["R"], fit_scale="sqrt")


def test_plot_runs(pair):
    import matplotlib

    matplotlib.use("Agg")
    res = cp.tl.estimate_time_alignment(pair, ["R"], current={"ctx": (0.65, -10.0)}, n_boot=50)
    fig = cp.pl.time_alignment(res)
    assert len(fig.axes) == 4  # two own-clock panels + estimate and current overlays


def test_reads_tradeseq_sce(tmp_path):
    pytest.importorskip("rpy2")
    from rpy2.robjects import r

    paths = {}
    for name, (tags, shift) in {"A": ("seq(0, 110, by = 10)", 0), "B": ("seq(0, 165, by = 15)", 30)}.items():
        p = tmp_path / f"{name}.rds"
        r(f"""
          suppressPackageStartupMessages({{library(tradeSeq); library(SingleCellExperiment)}})
          set.seed(2)
          tags <- {tags}; per <- 150
          minutes <- rep(tags, each = per); n <- length(minutes)
          sf <- exp(rnorm(n, 0, 0.2))
          mu <- 0.5 * (1 + 19 / (1 + exp(-(minutes - 20 - {shift}) / 5)))
          counts <- rbind(R = rpois(n, sf * mu), flat = rpois(n, sf * 4), flat2 = rpois(n, sf * 2))
          colnames(counts) <- paste0("c", seq_len(n))
          sce <- fitGAM(counts, pseudotime = matrix(minutes / max(minutes), ncol = 1), cellWeights = matrix(1, n, 1),
                        offset = log(sf), nknots = 5, verbose = FALSE)
          saveRDS(sce, "{p}")
        """)
        paths[name] = str(p)
    res = cp.tl.estimate_time_alignment(paths, ["R"], time_map={"A": 110, "B": 165}, n_boot=100)
    tg = res.tags.groupby("context")["minutes"].agg(["min", "max", "nunique"])
    assert tg.loc["A"].tolist() == [0, 110, 12] and tg.loc["B"].tolist() == [0, 165, 12]
    r_ = res.alignment.set_index("context").loc["B"]
    assert r_.scale_lo <= 1.0 <= r_.scale_hi and r_.offset_lo <= -30 <= r_.offset_hi
