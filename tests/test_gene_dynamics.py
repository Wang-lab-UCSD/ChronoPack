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
    full_paths = _gene_dynamics.prepare_tradeseq_inputs(adata, tmp_path, tag="full")

    counts = pd.read_csv(paths["counts"], index_col=0)
    times = pd.read_csv(paths["cell_times"])
    genes = pd.read_csv(paths["subset_genes"])
    libs = pd.read_csv(full_paths["cell_libsize_full"])

    assert counts.shape == (3, 2)
    assert counts.index.tolist() == ["g0", "g1", "g2"]
    assert counts.columns.tolist() == ["c0", "c1"]
    assert times.to_dict(orient="list") == {"cell_id": ["c0", "c1"], "time": [1, 2]}
    assert genes["gene"].tolist() == ["g0", "g1", "g2"]
    assert libs.to_dict(orient="list") == {"cell_id": ["c0", "c1"], "libsize_full": [4, 3]}


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
