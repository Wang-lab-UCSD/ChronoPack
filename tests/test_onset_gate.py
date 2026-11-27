import numpy as np
import pandas as pd
import pytest


def test_onset_ignores_departures_against_the_direction():
    pytest.importorskip("rpy2")
    from rpy2.robjects import FloatVector, r

    from chronopack._backends.r._bridge import source_r_resource

    source_r_resource("dimred_clustering_helpers.R")
    k = np.arange(1, 98, dtype=float)
    # Early small hump, a dip below the start, then a larger late rise: induced, but the first long departure
    # from the early level is the dip.
    y = (0.3 + 0.05 * np.exp(-0.5 * ((k - 10) / 5) ** 2) - 0.08 * np.exp(-0.5 * ((k - 35) / 8) ** 2)
         + 0.15 / (1 + np.exp(-(k - 65) / 4)))
    r.assign("y_onset_test", FloatVector(y))
    r.assign("tt_onset_test", FloatVector(k / len(k)))
    res = r("""local({
      m <- matrix(y_onset_test, nrow = 1, dimnames = list("g_A", NULL))
      get_descriptors(m, tt_onset_test, N_minutes = 10 / length(tt_onset_test), z_diff_threshold = 0.2,
                      require_consecutive = 30, onset_mode = "horizontal", pseudocount = 0)
    })""")
    assert res.rx2("direction")[0] == "induced"
    assert 50 < res.rx2("t_onset")[0] * len(k) < 70  # the late rise, not the dip
    assert "t_onset_signed" not in list(res.names)


@pytest.fixture(scope="module")
def tiny_traj(tmp_path_factory):
    pytest.importorskip("rpy2")
    from rpy2.robjects import r

    import chronopack as cp

    d = tmp_path_factory.mktemp("onset_gate")
    sce_path = d / "sce.rds"
    r(f"""
      suppressPackageStartupMessages({{library(tradeSeq); library(SingleCellExperiment)}})
      set.seed(1)
      tags <- seq(0, 110, by = 10); per <- 80
      minutes <- rep(tags, each = per); n <- length(minutes)
      sf <- exp(rnorm(n, 0, 0.2))
      mu_up <- 4 * ifelse(minutes >= 40, 3, 1)
      counts <- rbind(up = rpois(n, sf * mu_up), flat = rpois(n, sf * 4), low = rpois(n, sf * 0.02))
      colnames(counts) <- paste0("c", seq_len(n))
      sce <- fitGAM(counts, pseudotime = matrix(minutes / 110, ncol = 1), cellWeights = matrix(1, n, 1),
                    offset = log(sf), nknots = 6, verbose = FALSE)
      for (f in c("waldStat_1", "meanLogFC", "pvalue_1", "qval_1")) rowData(sce)[[f]] <- 1
      saveRDS(sce, "{sce_path}")
    """)
    traj = cp.tl.align_gene_trajectories(
        {"A": str(sce_path), "B": str(sce_path)},
        time_map={"A": 110, "B": 110},
        offset_map={"A": 0, "B": 0},
        scale_map={"A": 1.0, "B": 1.0},
    )
    return traj


def test_onset_gate_variants(tiny_traj):
    import chronopack as cp

    base = cp.tl.compute_gene_trajectory_descriptors(tiny_traj, inplace=False)
    assert not {"onset_gate_pass", "t_onset_gated"} & set(base.obs.columns)

    v1 = cp.tl.compute_gene_trajectory_descriptors(tiny_traj, inplace=False, onset_gate="band")
    obs = v1.obs.set_index(v1.obs["gene"] + "_" + v1.obs["cell_type"])
    np.testing.assert_allclose(obs["t_onset"], base.obs.set_index(obs.index)["t_onset"].to_numpy(), equal_nan=True)
    assert not {"onset_band_pass", "onset_tag_span"} & set(obs.columns)
    assert bool(obs.loc["up_A", "onset_gate_pass"])
    assert not bool(obs.loc["flat_A", "onset_gate_pass"]) and np.isnan(obs.loc["flat_A", "t_onset_gated"])
    assert obs.loc["up_A", "t_onset_gated"] == pytest.approx(obs.loc["up_A", "t_onset"])
    assert 20 <= obs.loc["up_A", "t_onset_gated"] <= 50
    # Same fit under both context names: same answer.
    assert obs.loc["up_A", "t_onset_gated"] == obs.loc["up_B", "t_onset_gated"]

    v2 = cp.tl.compute_gene_trajectory_descriptors(tiny_traj, inplace=False, onset_gate="band_or_tags")
    obs2 = v2.obs.set_index(v2.obs["gene"] + "_" + v2.obs["cell_type"])
    assert not {"onset_band_pass", "onset_tag_span"} & set(obs2.columns)
    assert (obs2["onset_gate_pass"] | ~obs["onset_gate_pass"]).all()  # variant 2 passes everything variant 1 does
    assert bool(obs2.loc["up_A", "onset_gate_pass"]) and not bool(obs2.loc["flat_A", "onset_gate_pass"])

    with pytest.raises(ValueError, match="onset_gate must be"):
        cp.tl.compute_gene_trajectory_descriptors(tiny_traj, inplace=False, onset_gate="variant3")
    with pytest.raises(ValueError, match="Unknown onset_gate_args"):
        cp.tl.compute_gene_trajectory_descriptors(tiny_traj, inplace=False, onset_gate="band",
                                                  onset_gate_args={"deltta": 0.1})


def test_onset_gate_tag_span_and_band_from_sce(tiny_traj):
    from rpy2.robjects import FloatVector, StrVector, r

    sce = tiny_traj.uns["chronopack"]["alignment"]["source_sce"]["A"]
    res = r["chronopack_onset_gate_from_sce"](
        str(sce), StrVector(["up", "flat", "low"]), FloatVector([1, 1, 1]),
        FloatVector(list(range(1, 111))), 110.0, tags=True)
    span = dict(zip(res.rx2("gene"), res.rx2("onset_tag_span")))
    band = dict(zip(res.rx2("gene"), res.rx2("onset_band_pass")))
    assert span["up"] >= 60  # tags 40..110 are all ~3x the baseline
    assert span["flat"] < 30
    assert span["low"] == 0  # under 10 UMIs per tag
    assert band["up"] and not band["flat"]


def test_onset_gate_overrides(tiny_traj):
    import chronopack as cp

    strict = cp.tl.compute_gene_trajectory_descriptors(
        tiny_traj, inplace=False, onset_gate="band_or_tags", onset_gate_args={"delta": 5.0, "tag_span_minutes": 30})
    assert not strict.obs["onset_gate_pass"].any()
    assert pd.isna(strict.obs["t_onset_gated"]).all()
