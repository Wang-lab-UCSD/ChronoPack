import numpy as np
import pytest
import pandas as pd
from anndata import AnnData

import chronopack as cp
from chronopack._core._keys import CHRONOPACK_UNS_KEY
from chronopack.tl import _gene_dynamics


def _adata():
    return AnnData(
        X=np.array([[1, 0, 3], [0, 2, 1]], dtype=float),
        obs=pd.DataFrame({"time": [1, 2], "total_counts": [4, 3]}, index=["c0", "c1"]),
        var=pd.DataFrame(index=["g0", "g1", "g2"]),
        layers={"counts": np.array([[1, 0, 3], [0, 2, 1]], dtype=int)},
    )


def test_prepare_tradeseq_inputs_writes_expected_files(tmp_path):
    adata = _adata()

    paths = _gene_dynamics.prepare_tradeseq_inputs(adata, tmp_path, tag="")

    counts = pd.read_csv(paths["counts"], index_col=0)
    times = pd.read_csv(paths["cell_times"])
    genes = pd.read_csv(paths["subset_genes"])
    libs = pd.read_csv(paths["cell_libsize_full"])

    assert counts.shape == (3, 2)
    assert counts.index.tolist() == ["g0", "g1", "g2"]
    assert counts.columns.tolist() == ["c0", "c1"]
    assert times.to_dict(orient="list") == {"cell_id": ["c0", "c1"], "time": [1, 2]}
    assert genes["gene"].tolist() == ["g0", "g1", "g2"]
    assert libs.to_dict(orient="list") == {"cell_id": ["c0", "c1"], "libsize_full": [4, 3]}
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        "cell_libsize_full.csv", "cell_times.csv.gz", "counts_genes_x_cells.csv.gz", "subset_genes.csv",
    ]


def test_library_sizes_come_from_obs_total_counts_not_the_filtered_genes(tmp_path):
    adata = _adata()[:, ["g0"]].copy()  # gene filtering leaves obs["total_counts"] (4, 3) untouched

    paths = _gene_dynamics.prepare_tradeseq_inputs(adata, tmp_path)

    assert pd.read_csv(paths["cell_libsize_full"])["libsize_full"].tolist() == [4, 3]


def _counts_lines(path):
    import gzip

    with gzip.open(path, "rt") as handle:
        return handle.read().splitlines()


def test_prepare_tradeseq_inputs_writes_integer_valued_float_counts_as_integers(tmp_path):
    adata = _adata()
    adata.layers["counts"] = adata.layers["counts"].astype(float)

    paths = _gene_dynamics.prepare_tradeseq_inputs(adata, tmp_path, tag="")

    assert _counts_lines(paths["counts"]) == [",c0,c1", "g0,1,0", "g1,0,2", "g2,3,1"]


def test_prepare_tradeseq_inputs_keeps_non_integer_layers_as_floats(tmp_path):
    adata = _adata()
    adata.layers["cpm"] = np.array([[0.5, 0.0, 3.0], [0.0, 2.0, 1.25]])

    paths = _gene_dynamics.prepare_tradeseq_inputs(adata, tmp_path, layer="cpm", tag="")

    assert _counts_lines(paths["counts"]) == [",c0,c1", "g0,0.5,0.0", "g1,0.0,2.0", "g2,3.0,1.25"]


def test_prepare_tradeseq_inputs_writes_sparse_float_counts_as_integers(tmp_path):
    import scipy.sparse as sp

    adata = _adata()
    adata.layers["counts"] = sp.csr_matrix(adata.layers["counts"].astype(float))

    paths = _gene_dynamics.prepare_tradeseq_inputs(adata, tmp_path, tag="")

    assert _counts_lines(paths["counts"]) == [",c0,c1", "g0,1,0", "g1,0,2", "g2,3,1"]


def test_fit_gene_dynamics_can_prepare_inputs_without_running_backend(tmp_path):
    adata = _adata()

    result = cp.tl.fit_gene_dynamics(
        adata,
        output_dir=tmp_path,
        prepare_inputs=True,
        run_backend=False,
        key_added="test_fit",
    )

    assert result.backend == "tradeseq"
    assert result.key_added == "test_fit"
    assert (tmp_path / "counts_genes_x_cells.csv.gz").exists()
    assert (tmp_path / "cell_libsize_full.csv").exists()
    assert CHRONOPACK_UNS_KEY in adata.uns
    assert "test_fit" in adata.uns[CHRONOPACK_UNS_KEY]
    assert adata.uns[CHRONOPACK_UNS_KEY]["test_fit"]["backend"] == "tradeseq"


def test_fit_gene_dynamics_writes_association_stats_with_mocked_backend(tmp_path, monkeypatch):
    adata = _adata()
    association = pd.DataFrame(
        {
            "waldStat_1": [10.0, 5.0],
            "df_1": [3.0, 3.0],
            "pvalue_1": [0.001, 0.2],
            "meanLogFC": [1.5, 0.1],
            "qval_1": [0.003, 0.2],
        },
        index=["g0", "g2"],
    )
    significant = association.iloc[[0]].copy()

    def fake_run_backend(**kwargs):
        return object()

    def fake_read_outputs(result, input_dir):
        return association, significant, ["g0"], {
            "association": str(tmp_path / "asso.csv"),
            "significant": str(tmp_path / "sig.csv"),
            "sce": str(tmp_path / "sce.rds"),
            "sce_screen": str(tmp_path / "sce_screen.rds"),
            "screen_genes": str(tmp_path / "screen_genes.csv"),
        }

    monkeypatch.setattr(_gene_dynamics, "_run_tradeseq_backend", fake_run_backend)
    monkeypatch.setattr(_gene_dynamics, "_read_tradeseq_outputs", fake_read_outputs)

    result = cp.tl.fit_gene_dynamics(
        adata,
        input_dir=tmp_path,
        key_added="mock_fit",
    )

    assert result.selected_genes == ["g0"]
    assert adata.var["mock_fit_significant"].tolist() == [True, False, False]
    np.testing.assert_allclose(adata.var.loc[["g0", "g2"], "mock_fit_waldStat_1"], [10.0, 5.0])
    assert np.isnan(adata.var.loc["g1", "mock_fit_waldStat_1"])
    sidecars = adata.uns[CHRONOPACK_UNS_KEY]["mock_fit"]["result_keys"]["sidecar"]
    assert sidecars["sce"].endswith("sce.rds")
    assert sidecars["sce_screen"].endswith("sce_screen.rds")


def test_tag_bias_helpers_are_public_tl_api():
    assert cp.tl.identify_tag_bias_candidates is _gene_dynamics.identify_tag_bias_candidates
    assert cp.tl.filter_tag_bias_candidates is _gene_dynamics.filter_tag_bias_candidates
    assert cp.tl.tag_bias_exclusion_genes is _gene_dynamics.tag_bias_exclusion_genes
    for name in ["identify_tag_bias_candidates", "filter_tag_bias_candidates", "tag_bias_exclusion_genes"]:
        assert name in cp.tl.__all__


def test_prepare_tradeseq_inputs_writes_aligned_cell_weights(tmp_path):
    adata = _adata()
    weights = pd.DataFrame({"0": [0.25, 1.0], "1": [0.75, 0.0]}, index=["c1", "c0"])  # out of cell order

    paths = _gene_dynamics.prepare_tradeseq_inputs(adata, tmp_path, cell_weights=weights)

    written = pd.read_csv(paths["cell_weights"])
    assert written.to_dict(orient="list") == {"cell_id": ["c0", "c1"], "0": [1.0, 0.25], "1": [0.0, 0.75]}


def test_fit_gene_dynamics_weighted_lineages_use_multilineage_runner_and_names(tmp_path, monkeypatch):
    adata = _adata()
    weights = pd.DataFrame({"0": [1.0, 0.25], "1": [0.0, 0.75]}, index=["c0", "c1"])
    cp.tl.fit_gene_dynamics(adata, output_dir=tmp_path, cell_weights=weights, prepare_inputs=True, run_backend=False)
    seen = {}

    def fake_run_backend(**kwargs):
        seen.update(kwargs)
        return object()

    monkeypatch.setattr(_gene_dynamics, "_run_tradeseq_backend", fake_run_backend)
    monkeypatch.setattr(
        _gene_dynamics, "_read_tradeseq_outputs",
        lambda result, input_dir: (pd.DataFrame(), pd.DataFrame(), [], {}),
    )

    cp.tl.fit_gene_dynamics(input_dir=tmp_path, lineage_weights="cell_weights.csv")

    assert seen["lineage_weights"] == "cell_weights.csv"
    assert seen["outputs"]["association"] == "asso_df.2trajectory.nknots4.csv"
    assert seen["outputs"]["sce"] == "sce.2trajectory.nknots8.rds"


def test_fit_gene_dynamics_with_genes_fits_only_those_genes(tmp_path, monkeypatch):
    calls = {}

    def fake_fit_selected(**kwargs):
        calls.update(kwargs)
        association = pd.DataFrame({"waldStat_1": [3.0], "pvalue_1": [0.01]}, index=["g1"])
        return association, ["g1"], {"sce": str(tmp_path / "sce.rds")}

    def no_full_run(**kwargs):
        raise AssertionError("the full screen + refit must not run")

    monkeypatch.setattr(_gene_dynamics, "_fit_selected_genes", fake_fit_selected)
    monkeypatch.setattr(_gene_dynamics, "_run_tradeseq_backend", no_full_run)

    fit = cp.tl.fit_gene_dynamics(input_dir=tmp_path, genes=["g1"], final_nknots=8)

    assert calls["genes"] == ["g1"] and calls["nknots"] == 8
    assert fit.selected_genes == ["g1"]
    assert fit.sidecar_paths["sce"].endswith("sce.rds")
    assert list(fit.association.index) == ["g1"]


def test_signed_max_log2fc_picks_the_largest_effect_per_gene():
    import pytest

    pytest.importorskip("rpy2")
    from chronopack._backends.r._bridge import source_r_resource
    from rpy2.robjects import StrVector, r

    source_r_resource("tradeseq_helpers.R")
    # predictSmooth() stand-in: g1 rises 1 -> 4 on lineage 1 and falls 1 -> 0.3 on lineage 2
    # (|log2FC| 2 vs 1.74); g2 falls 8 -> 1 on its only lineage.
    r("""predictSmooth <- function(models, gene, nPoints, tidy) data.frame(
           lineage = c(1, 1, 1, 2, 2, 2, 1, 1, 1), time = rep(c(0, 0.5, 1), 3),
           gene = rep(c("g1", "g1", "g2"), each = 3),
           yhat = c(1, 4, 2, 1, 0.5, 0.3, 8, 4, 1))""")
    try:
        out = r["signed_max_log2fc_start_many"](r("NULL"), StrVector(["g2", "g1", "g3"]), pseudocount=0)
    finally:
        r("rm(predictSmooth, envir = globalenv())")
    np.testing.assert_allclose(list(out)[:2], [-3.0, 2.0])
    assert list(out.names) == ["g2", "g1", "g3"] and np.isnan(list(out)[2])


def test_signed_max_log2fc_from_start_reads_the_screen_fit(tmp_path):
    import pytest

    pytest.importorskip("rpy2")
    from rpy2.robjects import r

    import chronopack as cp

    sce_path = tmp_path / "sce_screen.1trajectory.nknots4.rds"
    r(f"""
      suppressPackageStartupMessages({{library(tradeSeq); library(SingleCellExperiment)}})
      set.seed(1)
      n <- 60; pt <- seq(0, 1, length.out = n)
      counts <- rbind(up = rpois(n, 2 + 20 * pt), down = rpois(n, 20 - 18 * pt), flat = rpois(n, 5))
      colnames(counts) <- paste0("c", seq_len(n))
      sce <- fitGAM(counts, pseudotime = matrix(pt, ncol = 1), cellWeights = matrix(1, n, 1),
                    nknots = 3, verbose = FALSE)
      saveRDS(sce, "{sce_path}")
      pr <- predictSmooth(sce, gene = rownames(sce), nPoints = 50, tidy = TRUE)
    """)
    pr = pd.DataFrame({k: list(r(f"pr${k}")) for k in ["gene", "time", "yhat"]})
    expected = {}
    for gene, d in pr.groupby("gene"):
        fc = np.log2((d.sort_values("time")["yhat"].to_numpy() + 0.01) / (d.sort_values("time")["yhat"].iloc[0] + 0.01))
        expected[gene] = fc[np.argmax(np.abs(fc))]

    res = cp.tl.signed_max_log2fc_from_start(tmp_path, ["flat", "up", "down"], n_points=50)
    assert res.name == "signedMaxLog2FC_start" and list(res.index) == ["flat", "up", "down"]
    np.testing.assert_allclose(res.to_numpy(), [expected[g] for g in res.index], rtol=1e-10)
    assert res["up"] > 1 and res["down"] < -1

    with pytest.raises(Exception, match="Genes not in the tradeSeq fit: nope"):
        cp.tl.signed_max_log2fc_from_start(tmp_path, ["up", "nope"])


def test_combine_context_associations_pairs_and_fills_contexts():
    import pytest

    a = pd.DataFrame(
        {"waldStat_1": [5.0, 9.0, 1.0], "pvalue_1": [0.01, 0.001, 0.5], "qval_1": [0.02, 0.003, 0.6],
         "signedMaxLog2FC_start": [1.0, -2.0, 0.1]},
        index=["g1", "g2", "g3"],
    )
    b = pd.DataFrame(
        {"waldStat_1": [7.0, 2.0], "pvalue_1": [0.002, 0.2], "qval_1": [0.004, np.nan],
         "signedMaxLog2FC_start": [0.5, 3.0]},
        index=["g2", "g4"],
    )
    both = cp.tl.combine_context_associations(a, b)
    assert list(both.index) == ["g2"]
    assert both.loc["g2", "qval"] == 0.004 and both.loc["g2", "log2FC_t2"] == 0.5
    assert list(both.columns) == ["waldStat_1", "waldStat_2", "qval", "log2FC_t1", "log2FC_t2", "pvalue_1", "pvalue_2"]

    union = cp.tl.combine_context_associations(a, b, qval_fun="min", fill_missing=True, arrange_by_wald_1=False)
    assert list(union.index) == ["g1", "g2", "g3", "g4"]
    assert union.loc["g2", "qval"] == 0.003
    assert union.loc["g4", ["waldStat_1", "log2FC_t1", "pvalue_1"]].tolist() == [0.0, 0.0, 1.0]
    assert union.loc["g4", "qval"] == 1.0  # NaN q-value filled with 1, as replace_na does
    assert union.loc["g1", ["waldStat_2", "pvalue_2"]].tolist() == [0.0, 1.0]
    assert list(cp.tl.combine_context_associations(a, b, fill_missing=True).index[:2]) == ["g2", "g1"]

    outer = cp.tl.combine_context_associations(a, b, keep_only_common=False)
    assert outer["waldStat_2"].isna().sum() == 2
    assert outer.loc["g1", "qval"] == 0.02 and np.isnan(outer.loc["g4", "qval"])
    with pytest.raises(KeyError, match="signedMaxLog2FC_start"):
        cp.tl.combine_context_associations(a.drop(columns="signedMaxLog2FC_start"), b)


def test_observed_gene_dynamics_scales_counts_to_the_reference_library_size(tmp_path):
    adata = _adata()
    weights = pd.DataFrame({"0": [2.0, 0.25], "1": [0.0, 0.75]}, index=["c0", "c1"])  # c0 not normalized
    _gene_dynamics.prepare_tradeseq_inputs(adata, tmp_path, cell_weights=weights)

    obs = cp.tl.observed_gene_dynamics(tmp_path, ["g2", "g0"], time_max=110)

    sf = np.array([4, 3]) / 3.5
    expected = np.log1p(np.array([3, 1]) * np.sqrt(sf.prod()) / sf)
    assert obs["gene"].tolist() == ["g2", "g2", "g0", "g0"]
    assert obs["cell_id"].tolist() == ["c0", "c1", "c0", "c1"]
    assert obs["t"].tolist() == [0.0, 110.0, 0.0, 110.0]
    np.testing.assert_allclose(obs["expr_log1p"].to_numpy()[:2], expected)
    assert obs["weight"].eq(1).all() and obs["assigned"].all()

    lin2 = cp.tl.observed_gene_dynamics(tmp_path, ["g2"], lineage=2, lineage_weights="cell_weights.csv")
    assert lin2["weight"].tolist() == [0.0, 0.75]
    assert lin2["assigned"].tolist() == [False, True]
    with pytest.raises(ValueError):
        cp.tl.observed_gene_dynamics(tmp_path, ["g2"], lineage=2)
    with pytest.raises(cp._core._validation.ChronoPackError):
        cp.tl.observed_gene_dynamics(tmp_path, ["MISSING"])
