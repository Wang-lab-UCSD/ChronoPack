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


def test_prepare_tradeseq_inputs_accepts_counts_only_anndata(tmp_path):
    adata = _adata()
    adata.X = None  # counts only in a layer

    paths = _gene_dynamics.prepare_tradeseq_inputs(adata, tmp_path)

    assert pd.read_csv(paths["counts"], index_col=0).shape == (3, 2)


def test_prepare_tradeseq_inputs_sparse_counts_match_dense_text_across_chunks(tmp_path):
    import scipy.sparse as sp

    rng = np.random.default_rng(0)
    counts = rng.poisson(0.3, size=(7, 11)).astype(float)
    obs = pd.DataFrame({"time": np.arange(7) % 3 + 1, "total_counts": counts.sum(axis=1)}, index=[f"c{i}" for i in range(7)])
    var = pd.DataFrame(index=[f"g{j}" for j in range(11)])
    dense = AnnData(X=counts, obs=obs, var=var, layers={"counts": counts})
    sparse = AnnData(X=counts, obs=obs, var=var, layers={"counts": sp.csr_matrix(counts)})

    expected = _counts_lines(_gene_dynamics.prepare_tradeseq_inputs(dense, tmp_path / "dense")["counts"])
    assert _counts_lines(_gene_dynamics.prepare_tradeseq_inputs(sparse, tmp_path / "sparse")["counts"]) == expected

    genes, cells = np.array(var.index), np.array(obs.index)
    _gene_dynamics._write_genes_x_cells_csv(counts.T, genes, cells, tmp_path / "dense.csv.gz")
    # Two genes per block (7 cells, at most 14 values per block), so blocks of the sparse matrix are written in turn.
    _gene_dynamics._write_genes_x_cells_csv(sp.csr_matrix(counts.T), genes, cells, tmp_path / "chunked.csv.gz", max_values_per_chunk=14)
    assert _counts_lines(tmp_path / "chunked.csv.gz") == _counts_lines(tmp_path / "dense.csv.gz")
    assert len(_counts_lines(tmp_path / "chunked.csv.gz")) == 12


def _time_course_adata(sparse_counts=False, seed=0):
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)
    n_per, times = 60, np.arange(1, 7)
    time = np.repeat(times, n_per)
    libsize = rng.integers(2000, 6000, size=time.size)
    base = rng.uniform(0.0005, 0.005, size=8)
    shape = np.ones((time.size, 8))
    shape[:, 0] = np.exp(0.5 * (time - 1))                       # rising
    shape[:, 1] = np.exp(-0.5 * (time - 3.5) ** 2) * 4 + 1       # transient
    counts = rng.poisson(libsize[:, None] * base * shape).astype(float)
    obs = pd.DataFrame({"time": time, "total_counts": libsize}, index=[f"c{i}" for i in range(time.size)])
    layer = sp.csr_matrix(counts) if sparse_counts else counts
    return AnnData(X=layer.copy(), obs=obs, var=pd.DataFrame(index=[f"g{j}" for j in range(8)]), layers={"counts": layer})


def test_prescreen_time_association_matches_a_least_squares_f_test():
    from scipy import stats

    adata = _time_course_adata()
    res = cp.tl.prescreen_time_association(adata, df=3, p_threshold=0.2, chunk_size=3)

    assert res.index.tolist() == adata.var_names.tolist()
    assert res.loc["g0", "pvalue"] < 1e-10 and res.loc["g1", "pvalue"] < 1e-10
    assert res["keep"].tolist() == (res["pvalue"] < 0.2).tolist()
    # Reference: ordinary least squares on [1, t, t^2, t^3] for one gene.
    t = adata.obs["time"].to_numpy(float)
    lib = adata.obs["total_counts"].to_numpy(float)
    y = np.log1p(adata.layers["counts"][:, 5] * (np.median(lib) / lib))
    X = np.column_stack([np.ones_like(t), t, t**2, t**3])
    rss1 = np.sum((y - X @ np.linalg.lstsq(X, y, rcond=None)[0]) ** 2)
    rss0 = np.sum((y - y.mean()) ** 2)
    F = (rss0 - rss1) / 3 / (rss1 / (len(y) - 4))
    np.testing.assert_allclose(res.loc["g5", "F"], F, rtol=1e-8)
    np.testing.assert_allclose(res.loc["g5", "pvalue"], stats.f.sf(F, 3, len(y) - 4), rtol=1e-8)

    sparse = cp.tl.prescreen_time_association(_time_course_adata(sparse_counts=True), df=3)
    np.testing.assert_allclose(sparse["F"].to_numpy(), res["F"].to_numpy(), rtol=1e-10)
    entry = adata.uns[CHRONOPACK_UNS_KEY]["time_association_prescreen"]
    assert entry["result_keys"]["n_genes_screened"] == 8
    with pytest.raises(cp._core._validation.ChronoPackError):
        cp.tl.prescreen_time_association(adata, df=6)


def test_prescreened_inputs_make_q_values_count_every_screened_gene(tmp_path, monkeypatch):
    import json

    adata = _time_course_adata()
    res = cp.tl.prescreen_time_association(adata, p_threshold=0.2)
    kept = adata[:, res["keep"]].copy()
    paths = _gene_dynamics.prepare_tradeseq_inputs(kept, tmp_path)
    record = json.loads((tmp_path / "prescreen_time_association.json").read_text())
    assert paths["prescreen"].endswith("prescreen_time_association.json")
    assert record["n_genes_screened"] == 8 and record["n_genes_written"] == kept.n_vars == record["n_genes_kept"]

    seen = {}

    def fake_run_backend(**kwargs):
        seen.update(kwargs)
        return object()

    association = pd.DataFrame({"waldStat_1": [1.0], "df_1": [3.0], "pvalue_1": [0.5], "meanLogFC": [0.1], "qval_1": [0.5]},
                               index=[kept.var_names[0]])
    monkeypatch.setattr(_gene_dynamics, "_run_tradeseq_backend", fake_run_backend)
    monkeypatch.setattr(_gene_dynamics, "_read_tradeseq_outputs",
                        lambda result, input_dir: (association, association.iloc[:0], [], {}))
    cp.tl.fit_gene_dynamics(input_dir=tmp_path, verbose=False)
    assert seen["fdr_n_tests"] == 8

    # Inputs rewritten without a pre-screen drop the record, and q-values count only the fitted genes again.
    _gene_dynamics.prepare_tradeseq_inputs(_time_course_adata(), tmp_path)
    assert not (tmp_path / "prescreen_time_association.json").exists()
    cp.tl.fit_gene_dynamics(input_dir=tmp_path, verbose=False)
    assert seen["fdr_n_tests"] is None


def test_r_bh_adjust_counts_prescreened_genes_as_p_one():
    pytest.importorskip("rpy2")
    from rpy2.robjects import FloatVector, NULL, r
    from chronopack._backends.r._bridge import source_r_resource

    source_r_resource("tradeseq_helpers.R")
    p = [0.001, 0.01, 0.02, 0.2, float("nan")]
    default = np.array(r["bh_adjust"](FloatVector(p), NULL))
    np.testing.assert_allclose(default, np.array(r["p.adjust"](FloatVector(p), method="BH")), equal_nan=True)
    padded = np.array(r["bh_adjust"](FloatVector(p), 10))
    full = np.array(r["p.adjust"](FloatVector(p[:4] + [1.0] * 6), method="BH"))[:4]
    np.testing.assert_allclose(padded[:4], full)
