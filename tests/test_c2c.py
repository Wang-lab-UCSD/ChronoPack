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
    )

    assert isinstance(result, cp.tl.c2c.C2CResult)
    assert result.interactions.shape[0] == 1
    row = result.interactions.iloc[0]
    assert row["sender"] == "S"
    assert row["receiver"] == "R"
    assert row["ligand"] == "LIG"
    assert row["receptor"] == "REC"
    # LIG x REC activity relative to t=0, summed over the 0-45 min window.
    profile = result.temporal_profiles.set_index("time")["lr_activity"]
    assert np.isclose(row["activity_sum"], (profile - profile.loc[0.0]).sum())
    assert row["n_window_points"] == 3
    assert result.parameters["causal_window"] == [0.0, 45.0]
    assert result.temporal_profiles is not None
    assert set(result.temporal_profiles["time"]) == {0.0, 15.0, 30.0}
    assert "receiver_response" not in result.temporal_profiles.columns

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
    )

    fig, _ = cp.pl.c2c.interaction_scatter(result, x_thresh=0.0)
    out = tmp_path / "c2c_scatter.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0

    response = pd.Series([0.0, 1.0, 2.0], index=[0.0, 15.0, 30.0])
    fig, _ = cp.pl.c2c.pair_overlay(result, ("LIG", "REC"), response=response)
    out = tmp_path / "c2c_overlay.png"
    fig.savefig(out)
    assert out.exists() and out.stat().st_size > 0


def test_c2c_min_cells_requires_detection_in_sender_and_receiver_cells():
    adata = _tiny_c2c_adata()
    # LIG has counts in 3 sender cells; NOREC is a receptor gene with no counts.
    # (Underscores in receptor names are split into components, so avoid them here.)
    adata.var_names = ["LIG", "REC", "RESP", "NOREC"]
    lr = pd.DataFrame({"ligand": ["LIG", "LIG"], "receptor": ["REC", "NOREC"]})
    association = pd.DataFrame({"pvalue_2": [1e-4]}, index=["LIG"])

    def pairs(**kwargs):
        res = cp.tl.c2c.infer(adata=adata, lr_resource=lr, association=association, sender="S", receiver="R", **kwargs)
        return set(res.candidate_interactions["receptor"])

    assert pairs() == {"REC"}  # default min_cells=1 drops the undetected receptor
    assert pairs(min_cells=0) == {"REC", "NOREC"}
    assert pairs(min_cells=3) == {"REC"}
    import pytest

    with pytest.raises(ValueError, match="No ligand-receptor pairs"):
        pairs(min_cells=4)  # LIG is detected in only 3 sender cells


def test_c2c_causal_window_limits_summed_time_points():
    adata = _tiny_c2c_adata()
    lr = pd.DataFrame({"ligand": ["LIG"], "receptor": ["REC"], "annotation": ["Secreted Signaling"]})
    association = pd.DataFrame({"pvalue_2": [1e-4]}, index=["LIG"])
    kwargs = dict(adata=adata, lr_resource=lr, association=association, sender="S", receiver="R")

    full = cp.tl.c2c.infer(**kwargs).interactions.iloc[0]
    early = cp.tl.c2c.infer(**kwargs, causal_window=(0, 15)).interactions.iloc[0]

    profile = cp.tl.c2c.infer(**kwargs).temporal_profiles.set_index("time")["lr_activity"]
    delta = profile - profile.loc[0.0]
    assert early["n_window_points"] == 2
    assert np.isclose(early["activity_sum"], delta.loc[[0.0, 15.0]].sum())
    assert np.isclose(full["activity_sum"], delta.sum())
