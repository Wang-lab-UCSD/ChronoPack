import numpy as np
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
