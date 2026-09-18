import numpy as np
import pandas as pd

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
