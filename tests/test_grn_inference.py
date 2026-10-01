from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import chronopack as cp


def _synthetic_scores():
    gene_scores = pd.DataFrame({"gene": ["G1", "G2", "TF1"], "lfdr_bum": [0.2, 0.8, 0.1]})
    tf_scores = pd.DataFrame({"TF": ["TF1", "TF2", "TF_ZERO"], "priority_act": [3.0, 1.0, 0.0]})
    return gene_scores, tf_scores


def test_grn_infer_aligns_targets_and_removes_self_edges():
    gene_scores, tf_scores = _synthetic_scores()

    result = cp.tl.grn.infer(
        gene_scores,
        tf_scores,
        valid_targets=["G1", "G2", "TF1", "MT-DROP", "ENSG000001"],
        name="methodA",
    )

    edges = result.edges
    assert list(edges.columns) == ["method", "regulator", "target", "score", "rank", "percentile"]
    assert set(edges["method"]) == {"methodA"}
    assert set(edges["regulator"]) == {"TF1", "TF2"}
    assert set(edges["target"]) == {"G1", "G2", "TF1"}
    assert not ((edges["regulator"] == "TF1") & (edges["target"] == "TF1")).any()
    assert len(edges) == 5
    assert edges["rank"].max() == 5
    assert edges["score"].is_monotonic_decreasing


def test_grn_infer_optional_tf_universe_is_explicit():
    gene_scores, tf_scores = _synthetic_scores()

    result = cp.tl.grn.infer(gene_scores, tf_scores, valid_targets=["G1", "G2", "TF1"], tf_universe=["TF2"])

    assert set(result.edges["regulator"]) == {"TF2"}
    assert result.parameters["tf_universe_size"] == 1


def test_grn_infer_accepts_score_targets_result():
    gene_scores, tf_scores = _synthetic_scores()
    targets = cp.tl.grn.TargetGeneScoreResult(scores=gene_scores, summary=pd.DataFrame(), params=pd.DataFrame())

    from_result = cp.tl.grn.infer(targets, tf_scores, valid_targets=["G1", "G2", "TF1"])
    from_table = cp.tl.grn.infer(gene_scores, tf_scores, valid_targets=["G1", "G2", "TF1"])

    pd.testing.assert_frame_equal(from_result.edges, from_table.edges)


def test_grn_result_legacy_schema_and_csv_roundtrip(tmp_path):
    gene_scores, tf_scores = _synthetic_scores()
    result = cp.tl.grn.infer(gene_scores, tf_scores, valid_targets=["G1", "G2", "TF1"], name="methodA")

    legacy = result.to_legacy_edges()
    assert {"method", "TF", "target", "EdgeWeight", "rank", "percentile"}.issubset(legacy.columns)

    path = Path(tmp_path) / "edges.csv"
    result.write_edges(path, legacy_schema=True)
    reread = pd.read_csv(path)
    pd.testing.assert_frame_equal(reread, legacy)


def test_grn_legacy_wrapper_matches_per_context_infer(tmp_path):
    gene_scores, tf_scores = _synthetic_scores()
    thp1_scores = pd.DataFrame({"gene": ["G1"], "lfdr_bum": [0.5]})
    combined = pd.concat([gene_scores.assign(dataset="K562"), thp1_scores.assign(dataset="THP1")], ignore_index=True)
    scores_path = Path(tmp_path) / "scores.csv"
    k562_path = Path(tmp_path) / "k562_tf.csv"
    thp1_path = Path(tmp_path) / "thp1_tf.csv"
    valid_path = Path(tmp_path) / "valid.csv"
    combined.to_csv(scores_path, index=False)
    tf_scores.to_csv(k562_path, index=False)
    tf_scores.to_csv(thp1_path, index=False)
    pd.Series(["G1", "G2", "TF1"]).to_csv(valid_path, index=False, header=False)

    compat = cp.tl.grn.build_chronopack_ranked_edges(scores_path, k562_path, thp1_path, valid_path)
    per_context = pd.concat(
        [
            cp.tl.grn.infer(gene_scores, k562_path, valid_targets=valid_path, name="ChronoPack_K562").to_legacy_edges(),
            cp.tl.grn.infer(thp1_scores, thp1_path, valid_targets=valid_path, name="ChronoPack_THP1").to_legacy_edges(),
        ],
        ignore_index=True,
    )

    pd.testing.assert_frame_equal(compat, per_context)


def test_grn_score_targets_one_context(tmp_path):
    pytest.importorskip("rpy2")
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(300)]
    p = np.concatenate([rng.uniform(0, 1e-4, 60), rng.uniform(0, 1, 240)])
    # Same layout as R's write.csv output: row numbers, then genes in an "X" column.
    association = pd.DataFrame({"X": genes, "pvalue_1": p, "lfdr": np.clip(p * 10, 0, 1)})
    association_path = Path(tmp_path) / "asso.csv"
    association.to_csv(association_path)
    valid_path = Path(tmp_path) / "valid.csv"
    pd.Series(genes[:250]).to_csv(valid_path, index=False, header=False)

    targets = cp.tl.grn.score_targets(association_path, valid_targets=valid_path)

    assert "dataset" not in targets.scores.columns
    assert list(targets.scores["gene"]) == genes[:250]
    assert targets.scores["lfdr_bum"].between(0, 1).all()
    assert targets.scores.loc[:59, "lfdr_bum"].median() < targets.scores.loc[60:, "lfdr_bum"].median()
    assert {"method", "param", "value"}.issubset(targets.params.columns)


def test_grn_r_resource_smoke():
    pytest.importorskip("rpy2")
    from chronopack._backends.r import get_r_resource_path
    from chronopack._backends.r._bridge import source_r_resource
    from rpy2.robjects import FloatVector, r

    with get_r_resource_path("grn_inference_helpers.R") as path:
        assert path.name == "grn_inference_helpers.R"

    source_r_resource("grn_inference_helpers.R")
    prepared = np.asarray(r["prepare_p"](FloatVector([0.0, 0.5, 1.0])))
    assert np.all(np.isfinite(prepared))
    assert prepared[0] > 0
    assert prepared[-1] < 1


def test_infer_breaks_equal_target_scores_by_wald_statistic():
    import chronopack as cp

    # Three targets with the same local FDR (equal edge weights), listed out of Wald order.
    targets = pd.DataFrame(
        {"gene": ["LOW", "HIGH", "MID"], "lfdr_bum": [0.0, 0.0, 0.0], "waldStat_1": [5.0, 50.0, 20.0]}
    )
    tfs = pd.DataFrame({"TF": ["TF1"], "priority_act": [1.0]})

    edges = cp.tl.grn.infer(targets, tfs, valid_targets=["LOW", "HIGH", "MID"]).edges

    assert edges.sort_values("rank")["target"].tolist() == ["HIGH", "MID", "LOW"]
