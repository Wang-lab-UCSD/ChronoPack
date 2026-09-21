from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import chronopack as cp


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
def test_standard_scanpy_wrappers_are_not_public_io_api():
    assert not hasattr(cp.io, "annotate_qc_metrics")
    assert not hasattr(cp.io, "filter_cells_by_qc")
    assert not hasattr(cp.io, "initialize_expression_layers")
    assert not hasattr(cp.io, "validate_preprocessed_adata")
