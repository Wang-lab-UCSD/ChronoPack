from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import chronopack as cp


def _synthetic_scores():
    gene_scores = pd.DataFrame(
        {
            "dataset": ["ctxA", "ctxA", "ctxA", "ctxB"],
            "gene": ["G1", "G2", "TF1", "G1"],
            "lfdr_bum": [0.2, 0.8, 0.1, 0.5],
        }
    )
    tf_scores = {
        "methodA": pd.DataFrame(
            {
                "TF": ["TF1", "TF2", "TF_ZERO"],
                "priority_act": [3.0, 1.0, 0.0],
            }
        )
    }
    return gene_scores, tf_scores


def test_grn_infer_aligns_targets_and_removes_self_edges():
    gene_scores, tf_scores = _synthetic_scores()

    result = cp.tl.grn.infer(
        gene_scores=gene_scores,
        tf_scores=tf_scores,
        valid_targets=["G1", "G2", "TF1", "MT-DROP", "ENSG000001"],
        contexts={"methodA": "ctxA"},
    )

    edges = result.edges
    assert list(edges.columns) == ["method", "regulator", "target", "score", "rank", "percentile"]
    assert set(edges["regulator"]) == {"TF1", "TF2"}
    assert set(edges["target"]) == {"G1", "G2", "TF1"}
    assert not ((edges["regulator"] == "TF1") & (edges["target"] == "TF1")).any()
    assert len(edges) == 5
    assert edges.groupby("method")["rank"].max().item() == 5
    assert edges["score"].is_monotonic_decreasing


def test_grn_infer_optional_tf_universe_is_explicit():
    gene_scores, tf_scores = _synthetic_scores()

    result = cp.tl.grn.infer(
        gene_scores=gene_scores,
        tf_scores=tf_scores,
        valid_targets=["G1", "G2", "TF1"],
        contexts={"methodA": "ctxA"},
        tf_universe=["TF2"],
    )

    assert set(result.edges["regulator"]) == {"TF2"}
    assert result.parameters["tf_universe_size"] == 1


def test_grn_result_legacy_schema_and_csv_roundtrip(tmp_path):
    gene_scores, tf_scores = _synthetic_scores()
    result = cp.tl.grn.infer(
        gene_scores=gene_scores,
        tf_scores=tf_scores,
        valid_targets=["G1", "G2", "TF1"],
        contexts={"methodA": "ctxA"},
    )

    legacy = result.to_legacy_edges()
    assert {"method", "TF", "target", "EdgeWeight", "rank", "percentile"}.issubset(legacy.columns)

    path = Path(tmp_path) / "edges.csv"
    result.write_edges(path, legacy_schema=True)
    reread = pd.read_csv(path)
    pd.testing.assert_frame_equal(reread, legacy)


def test_grn_infer_matches_legacy_wrapper_on_synthetic_files(tmp_path):
    gene_scores, tf_scores = _synthetic_scores()
    scores_path = Path(tmp_path) / "scores.csv"
    k562_path = Path(tmp_path) / "k562_tf.csv"
    thp1_path = Path(tmp_path) / "thp1_tf.csv"
    valid_path = Path(tmp_path) / "valid.csv"
    gene_scores.replace({"ctxA": "K562", "ctxB": "THP1"}).to_csv(scores_path, index=False)
    tf_scores["methodA"].to_csv(k562_path, index=False)
    tf_scores["methodA"].to_csv(thp1_path, index=False)
    pd.Series(["G1", "G2", "TF1"]).to_csv(valid_path, index=False, header=False)

    compat = cp.tl.grn.build_chronopack_ranked_edges(scores_path, k562_path, thp1_path, valid_path)
    direct = cp.tl.grn.infer(
        gene_scores=scores_path,
        tf_scores={"ChronoPack_K562": k562_path, "ChronoPack_THP1": thp1_path},
        valid_targets=valid_path,
        contexts={"ChronoPack_K562": "K562", "ChronoPack_THP1": "THP1"},
    ).to_legacy_edges()

    pd.testing.assert_frame_equal(compat, direct)


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
