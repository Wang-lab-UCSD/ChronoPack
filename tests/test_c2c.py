import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

import chronopack as cp
from chronopack.c2c import score_lr_pairs, split_genes, zscore


def test_split_genes_expands_receptor_shorthand():
    assert split_genes("TGFBR1_R2") == ["TGFBR1", "TGFBR2"]
    assert split_genes("IL12A_IL12B") == ["IL12A", "IL12B"]


def test_zscore_handles_constant_values():
    values = pd.Series([2.0, 2.0, 2.0])
    out = zscore(values)
    assert np.allclose(out.to_numpy(), 0.0)


def test_score_lr_pairs_ranks_pairs_and_ligands():
    pair_series = {
        ("TNF", "TNFRSF1A"): (
            pd.Series([1.0, 3.0, 5.0], index=[0, 15, 30]),
            pd.Series([0.0, 1.0, 2.0], index=[0, 15, 30]),
        ),
        ("CCL2", "CCR2"): (
            pd.Series([1.0, 1.5, 2.0], index=[0, 15, 30]),
            pd.Series([0.0, 0.2, 0.4], index=[0, 15, 30]),
        ),
    }
    association = pd.DataFrame({"pvalue_2": [1e-6, 0.05]}, index=["TNF", "CCL2"])
    detection = pd.Series({"TNFRSF1A": 0.5, "CCR2": 0.2})

    pair_scores, ligand_scores = score_lr_pairs(pair_series, association, detection)

    assert pair_scores.iloc[0]["pair"] == "TNF->TNFRSF1A"
    assert ligand_scores.iloc[0]["ligand"] == "TNF"
    assert {"composite_score", "zscore_composite_rank"}.issubset(ligand_scores.columns)


def _tiny_c2c_adata():
    obs = pd.DataFrame(
        {
            "cell_type": ["S", "S", "S", "R", "R", "R"],
            "time": ["1", "2", "3", "1", "2", "3"],
        },
        index=[f"cell{i}" for i in range(6)],
    )
    var = pd.DataFrame(index=["LIG", "REC", "RESP", "MISSING_REC"])
    X = np.array(
        [
            [1.0, 0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0, 0.0],
            [5.0, 0.0, 0.0, 0.0],
            [0.0, 4.0, 0.0, 0.0],
            [0.0, 4.0, 1.0, 0.0],
            [0.0, 4.0, 2.0, 0.0],
        ]
    )
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.layers["counts"] = sp.csr_matrix(X)
    adata.layers["log1p"] = np.log1p(X)
    return adata


def test_c2c_infer_returns_result_and_filters_missing_receptor(tmp_path):
    adata = _tiny_c2c_adata()
    lr = pd.DataFrame(
        {
            "ligand": ["LIG", "LIG"],
            "receptor": ["REC", "ABSENT"],
            "annotation": ["Secreted Signaling", "Secreted Signaling"],
        }
    )
    association = pd.DataFrame({"pvalue_2": [1e-4]}, index=["LIG"])

    result = cp.tl.c2c.infer(
        adata=adata,
        lr_resource=lr,
        association=association,
        sender="S",
        receiver="R",
        response_genes=["RESP"],
        causal_window_min=(0, 30),
        min_detect=0.0,
    )

    assert isinstance(result, cp.tl.c2c.C2CResult)
    assert result.interactions.shape[0] == 1
    row = result.interactions.iloc[0]
    assert row["sender"] == "S"
    assert row["receiver"] == "R"
    assert row["ligand"] == "LIG"
    assert row["receptor"] == "REC"
    assert row["score"] == row["avg"]
    assert result.temporal_profiles is not None
    assert set(result.temporal_profiles["time"]) == {0.0, 15.0, 30.0}

    paths = result.write(tmp_path)
    assert {"interactions_csv", "ligand_scores_csv", "metadata"}.issubset(paths)
    reread = pd.read_csv(paths["interactions_csv"])
    assert reread.shape == result.interactions.shape


def test_c2c_infer_missing_sender_has_clear_error():
    adata = _tiny_c2c_adata()
    lr = pd.DataFrame({"ligand": ["LIG"], "receptor": ["REC"]})
    association = pd.DataFrame({"pvalue_2": [1e-4]}, index=["LIG"])

    try:
        cp.tl.c2c.infer(
            adata=adata,
            lr_resource=lr,
            association=association,
            sender="ABSENT",
            receiver="R",
            response_genes=["RESP"],
        )
    except ValueError as exc:
        assert "No cells found for sender" in str(exc)
    else:
        raise AssertionError("Expected a missing-sender ValueError")


def test_c2c_plotting_smoke(tmp_path):
    adata = _tiny_c2c_adata()
    lr = pd.DataFrame({"ligand": ["LIG"], "receptor": ["REC"], "annotation": ["Secreted Signaling"]})
    association = pd.DataFrame({"pvalue_2": [1e-4]}, index=["LIG"])
    result = cp.tl.c2c.infer(
        adata=adata,
        lr_resource=lr,
        association=association,
        sender="S",
        receiver="R",
        response_genes=["RESP"],
    )

    fig, _ = cp.pl.c2c.interaction_scatter(result, max_labels=1)
    out = tmp_path / "c2c_scatter.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0

    fig, _ = cp.pl.c2c.pair_overlay(result, ("LIG", "REC"))
    out = tmp_path / "c2c_overlay.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0
