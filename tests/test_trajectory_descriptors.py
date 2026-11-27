import numpy as np
import pandas as pd
import pytest


def _descriptor_adata(n=41):
    ad = pytest.importorskip("anndata")

    k = np.arange(1, n + 1, dtype=float)
    ramp = (k - 1) / (n - 1)
    curves = {
        "rise_A": ramp,  # still rising at the end: censored, half change at k = 21
        "pulse_A": 2.0 * np.exp(-0.5 * ((k - 15) / 3.0) ** 2),  # peak at k = 15
        "fall_A": -1.5 * ramp,  # repressed, censored, half change at k = 21
        "flat_A": np.zeros(n),
    }
    X = 1.0 + np.vstack(list(curves.values()))
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=list(curves)),
        var=pd.DataFrame({"time": k / n}, index=[str(int(v)) for v in k]),
    )
    adata.obs["gene"] = [r.split("_")[0] for r in curves]
    adata.obs["cell_type"] = "A"
    return adata


def test_direction_half_time_and_censoring():
    pytest.importorskip("rpy2")
    import chronopack as cp

    adata = _descriptor_adata()
    cp.tl.compute_gene_trajectory_descriptors(adata)
    obs = adata.obs

    assert obs.loc["rise_A", "direction"] == "induced"
    assert obs.loc["fall_A", "direction"] == "repressed"
    assert pd.isna(obs.loc["flat_A", "direction"])

    assert bool(obs.loc["rise_A", "t_peak_censored"]) and bool(obs.loc["fall_A", "t_peak_censored"])
    assert not bool(obs.loc["pulse_A", "t_peak_censored"])
    assert not bool(obs.loc["flat_A", "t_peak_censored"])

    # Times are grid positions (no alignment scale map): t_peak of a censored gene is the window end.
    assert obs.loc["rise_A", "t_peak"] == pytest.approx(41)
    assert obs.loc["pulse_A", "t_peak"] == pytest.approx(15)
    assert obs.loc["rise_A", "t_half"] == pytest.approx(21)
    assert obs.loc["fall_A", "t_half"] == pytest.approx(21)
    assert obs.loc["pulse_A", "t_half"] < obs.loc["pulse_A", "t_peak"]
    assert pd.isna(obs.loc["flat_A", "t_half"])

    # Existing descriptors are still there.
    assert {"t_onset", "peak_log2fc", "auc_norm", "owp", "max_turn"} <= set(obs.columns)
