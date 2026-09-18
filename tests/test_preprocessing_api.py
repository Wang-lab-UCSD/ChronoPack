from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

import chronopack as cp
from chronopack._core import ChronoPackError


def _write_dge(path: Path, rows: dict[str, dict[str, object]]):
    pd.DataFrame.from_dict(rows, orient="index").to_csv(path)


def test_read_time_tagged_dge_matches_notebook_loader(tmp_path):
    p1 = tmp_path / "N701_input.csv"
    p2 = tmp_path / "N702_input.csv"
    _write_dge(
        p1,
        {
            "cellA": {"FINAL TIME TAG": 1, "COUNTS": 10, "GeneA": 1, "GeneB": 0},
            "cellB": {"FINAL TIME TAG": 4, "COUNTS": 20, "GeneA": 2, "GeneB": 3},
        },
    )
    _write_dge(
        p2,
        {
            "cellA": {"FINAL TIME TAG": 10, "COUNTS": 30, "GeneA": 0, "GeneC": 5},
        },
    )

    adata = cp.io.read_time_tagged_dge([p1, p2], sample_names=["N701", "N702"])

    assert adata.obs_names.tolist() == ["N701_cellA", "N701_cellB", "N702_cellA"]
    assert adata.var_names.tolist() == ["GeneA", "GeneB", "GeneC"]
    assert adata.obs["sample"].astype(str).tolist() == ["N701", "N701", "N702"]
    assert adata.obs["time"].tolist() == [1, 4, 10]
    assert adata.obs["time_group"].astype(str).tolist() == [
        "Group 1 (T1-3)",
        "Group 2 (T4-6)",
        "Group 4 (T10-12)",
    ]
    assert np.asarray(adata.X).tolist() == [[1.0, 0.0, 0.0], [2.0, 3.0, 0.0], [0.0, 0.0, 5.0]]


def test_qc_cell_filter_matches_reference_mask():
    adata = ad.AnnData(np.ones((4, 2)))
    adata.obs_names = ["a", "b", "c", "d"]
    adata.var_names = ["g1", "g2"]
    adata.obs["total_counts"] = [100, 250, 5000, 200]
    adata.obs["pct_counts_mt"] = [1.0, 8.0, 2.0, 2.0]
    adata.obs["n_genes_by_counts"] = [50, 400, 500, 250]

    out = cp.io.filter_cells_by_qc(
        adata.copy(),
        min_counts=200,
        max_counts=1000,
        max_pct_mt=7.5,
        min_genes=200,
    )

    reference_mask = (
        (adata.obs["total_counts"].to_numpy() >= 200)
        & (adata.obs["total_counts"].to_numpy() <= 1000)
        & (adata.obs["pct_counts_mt"].to_numpy() < 7.5)
        & (adata.obs["n_genes_by_counts"].to_numpy() >= 200)
    )
    assert out.obs_names.tolist() == adata.obs_names[reference_mask].tolist()
    assert out.obs_names.tolist() == ["d"]


def test_technical_and_low_count_gene_filters_match_current_rules():
    adata = ad.AnnData(np.array([[5, 1, 0, 4, 3], [5, 0, 0, 2, 1]], dtype=float))
    adata.obs_names = ["cell1", "cell2"]
    adata.var_names = ["GeneA", "MT-ND1", "ENSG0001", "RPS3", "GeneB"]
    adata.layers["counts"] = adata.X.copy()

    technical = cp.io.filter_technical_genes(adata.copy(), exclude_genes=["GeneB"])
    assert technical.var_names.tolist() == ["GeneA"]

    low_count = cp.io.filter_low_count_genes(adata.copy(), min_total_counts=6, counts_layer="counts")
    assert low_count.var_names.tolist() == ["GeneA", "RPS3"]
    assert low_count.var["total_counts"].tolist() == [10.0, 6.0]


def test_initialize_expression_layers_preserves_counts_and_adds_normalized_layers():
    adata = ad.AnnData(np.array([[1.0, 1.0], [2.0, 0.0]]))
    adata.obs_names = ["a", "b"]
    adata.var_names = ["g1", "g2"]

    cp.io.initialize_expression_layers(adata, target_sum=10)

    assert np.asarray(adata.layers["counts"]).tolist() == [[1.0, 1.0], [2.0, 0.0]]
    assert set(["counts", "cpm", "log1p"]).issubset(adata.layers.keys())
    assert np.allclose(np.asarray(adata.layers["cpm"]).sum(axis=1), [10.0, 10.0])


def test_validate_preprocessed_adata_reports_missing_contract_fields():
    adata = ad.AnnData(np.ones((2, 2)))
    adata.obs_names = ["a", "b"]
    adata.var_names = ["g1", "g2"]

    with pytest.raises(ChronoPackError, match="Missing adata.obs\\['time'\\]"):
        cp.io.validate_preprocessed_adata(adata, require_qc=False)

    adata.obs["time"] = ["bad", "2"]
    with pytest.raises(ChronoPackError, match="numeric"):
        cp.io.validate_preprocessed_adata(adata, require_qc=False)

    adata.obs["time"] = [1, 2]
    with pytest.raises(ChronoPackError, match="Missing adata.layers\\['counts'\\]"):
        cp.io.validate_preprocessed_adata(adata, require_qc=False)
