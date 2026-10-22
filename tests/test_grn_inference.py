import json
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


def test_omnipath_snapshot_roundtrip_is_used_instead_of_download(tmp_path, monkeypatch):
    from chronopack.tl.grn import _omnipath

    signaling = pd.DataFrame(
        {
            "source_genesymbol": ["NA", "B"],  # "NA" is a gene symbol, not a missing value
            "target_genesymbol": ["B", "C"],
            "is_directed": [True, True],
            "is_stimulation": [True, False],
            "is_inhibition": [False, True],
            "sources": ["SIGNOR", "SIGNOR;Reactome"],
            "n_sources": [1, 2],
            "references": [None, "SIGNOR:1"],
        }
    )
    transcriptional = pd.DataFrame({"source_genesymbol": ["C", "C"], "target_genesymbol": ["X", "Y"]})
    monkeypatch.setattr(_omnipath, "_download_omnipath_tables", lambda organism="human": (signaling, transcriptional))
    snap = cp.tl.grn.save_omnipath_snapshot(tmp_path / "snap", retrieved="2025-10-21", note="test")

    def no_download(organism="human"):
        raise AssertionError("snapshot should be read instead of downloading OmniPath")

    monkeypatch.setattr(_omnipath, "_download_omnipath_tables", no_download)
    tf_path, lr_path = tmp_path / "tfs.txt", tmp_path / "lr.csv"
    tf_path.write_text("C\n")
    pd.DataFrame({"ligand": ["L"], "receptor": ["R"]}).to_csv(lr_path)
    res = _omnipath.load_omnipath_tf_resources(human_tf_path=tf_path, ligand_receptor_path=lr_path, omnipath_snapshot=snap)

    pd.testing.assert_frame_equal(res.signaling, signaling)
    pd.testing.assert_frame_equal(res.transcriptional, transcriptional)
    meta = json.loads((snap / "metadata.json").read_text())
    assert meta["retrieved"] == "2025-10-21" and meta["n_rows"] == {"signaling": 2, "transcriptional": 2}


def _toy_tf_resources():
    from chronopack.tl.grn._omnipath import OmniPathTfResources

    edges = [  # source, target, stimulation, inhibition, consensus stimulation, consensus inhibition
        ("IL6R", "JAK1", True, False, True, False),
        ("IL6ST", "JAK1", True, False, True, False),
        ("JAK1", "STAT3", True, False, True, False),
        ("KIN", "AKT1", True, False, True, False),
        ("AKT1", "FOXO1", True, True, False, True),  # conflicting flags, consensus inhibition
        ("X", "DEAD", True, False, True, False),  # DEAD has no outgoing edges
    ]
    signaling = pd.DataFrame(edges, columns=["source_genesymbol", "target_genesymbol", "is_stimulation", "is_inhibition",
                                             "consensus_stimulation", "consensus_inhibition"])
    signaling["source"], signaling["target"] = signaling["source_genesymbol"], signaling["target_genesymbol"]
    signaling["is_directed"] = True
    signaling["sources"] = "SIGNOR"
    signaling["n_sources"] = 1
    transcriptional = pd.DataFrame({"source_genesymbol": ["STAT3", "FOXO1"], "target_genesymbol": ["T1", "T2"]})
    lr = pd.DataFrame({"ligand": ["L1", "L2"], "receptor": ["IL6R_IL6ST", "DEAD"]})
    return OmniPathTfResources(signaling=signaling, transcriptional=transcriptional, ligand_receptor=lr,
                               human_tfs={"STAT3", "FOXO1"})


def test_cellchat_receptor_complexes_are_seeded_at_their_subunits():
    from chronopack.tl.grn._omnipath import cellchat_receptor_subunits, ligand_receptor_seed

    assert cellchat_receptor_subunits("IL6R_IL6ST") == ["IL6R", "IL6ST"]
    assert cellchat_receptor_subunits("TGFbR1_R2") == ["TGFBR1", "TGFBR2"]
    assert cellchat_receptor_subunits("ACVR1_TGFbR") == ["ACVR1", "TGFBR1", "TGFBR2"]
    assert cellchat_receptor_subunits("TNFRSF1A") == ["TNFRSF1A"]
    lr = pd.DataFrame({"ligand": ["IL1B", "IL1B", "TNF", "TNF"],
                       "receptor": ["IL1R1_IL1RAP", "IL1R2", "TNFRSF1A", "TNFRSF1B"]})
    assert ligand_receptor_seed(lr, "IL1B") == {"IL1R1": 0.25, "IL1RAP": 0.25, "IL1R2": 0.5}
    assert ligand_receptor_seed(lr, "IL1B", nodes={"IL1R1", "IL1R2"}) == {"IL1R1": 0.5, "IL1R2": 0.5}
    assert ligand_receptor_seed(lr, "TNF") == {"TNFRSF1A": 0.5, "TNFRSF1B": 0.5}

    from chronopack.tl.grn._omnipath import build_omnipath_tf_graph, rank_omnipath_tfs

    res = _toy_tf_resources()
    bundle = build_omnipath_tf_graph(res)
    ranking, _ = rank_omnipath_tfs(bundle, res.ligand_receptor, source="L1", is_ligand=True)
    assert ranking.iloc[0]["TF"] == "STAT3" and ranking.iloc[0]["priority_act"] > 0


def test_tf_ranking_raises_without_signal_and_uses_consensus_signs():
    import pytest

    from chronopack.tl.grn._omnipath import build_omnipath_tf_graph, rank_omnipath_tfs

    res = _toy_tf_resources()
    bundle = build_omnipath_tf_graph(res)
    with pytest.raises(ValueError, match="DEAD"):
        rank_omnipath_tfs(bundle, res.ligand_receptor, source="L2", is_ligand=True)
    with pytest.raises(ValueError, match="No signal"):
        rank_omnipath_tfs(bundle, res.ligand_receptor, source="DEAD", is_ligand=False)

    with pytest.warns(UserWarning, match="net activating"):
        ranking, _ = rank_omnipath_tfs(bundle, res.ligand_receptor, source="KIN", is_ligand=False)
    assert ranking.set_index("TF").loc["FOXO1", "net"] < 0  # AKT1 -| FOXO1 by consensus
    old = build_omnipath_tf_graph(res, conflicting_signs="stimulation")
    ranking, _ = rank_omnipath_tfs(old, res.ligand_receptor, source="KIN", is_ligand=False, tied_signs="stimulation")
    assert ranking.set_index("TF").loc["FOXO1", "net"] > 0


def test_evidence_classes_weights_and_second_messenger_bridges():
    import numpy as np

    from chronopack.tl.grn._omnipath import (
        build_omnipath_tf_graph,
        build_signed_row_normalized_transitions,
        edge_evidence_weight,
        evidence_class,
        has_physical_evidence,
    )

    assert evidence_class("SIGNOR") == "activity flow"
    assert evidence_class("SIGNOR_ProtMapper") == "physical"
    assert evidence_class("PhosphoSite_noref") == "physical"
    assert evidence_class("REACH_ProtMapper") == "text-mined/predicted"
    assert evidence_class("STRING_talklr") == "text-mined/predicted"
    assert has_physical_evidence("SIGNOR;IntAct") and not has_physical_evidence("SIGNOR;SPIKE")
    assert has_physical_evidence("SIGNOR;NetworKIN_KEA")
    assert not has_physical_evidence("SIGNOR;NetworKIN_KEA", include_text_mined=False)
    assert not has_physical_evidence("SIGNOR;ProtMapper;REACH_ProtMapper", include_text_mined=False)
    assert has_physical_evidence("SIGNOR;ProtMapper;REACH_ProtMapper;PhosphoSite_ProtMapper", include_text_mined=False)
    attrs = {"sources": "SIGNOR|PhosphoSite|PhosphoSite_ProtMapper", "references": "SIGNOR:1;PhosphoSite:1;HPRD:3"}
    assert edge_evidence_weight(attrs, "uniform") == 1.0
    assert edge_evidence_weight(attrs, "resources") == 2.0  # PhosphoSite counted once
    assert edge_evidence_weight(attrs, "references") == 1.0 + np.log2(3)  # PMIDs 1 and 3

    res = _toy_tf_resources()
    sig = res.signaling.copy()
    sig.loc[sig["source_genesymbol"] == "IL6R", "sources"] = "SIGNOR;IntAct;HPRD"
    sig = pd.concat([sig, pd.DataFrame([{"source_genesymbol": "PLCG1", "target_genesymbol": "STAT3", "source": "PLCG1",
                                         "target": "STAT3", "is_directed": True, "is_stimulation": True,
                                         "is_inhibition": False, "consensus_stimulation": True,
                                         "consensus_inhibition": False, "sources": "SIGNOR", "n_sources": 1}])])
    res = type(res)(signaling=sig, transcriptional=res.transcriptional, ligand_receptor=res.ligand_receptor,
                    human_tfs=res.human_tfs)
    physical = build_omnipath_tf_graph(res, evidence="physical", second_messenger_bridges=False)
    assert set(physical.graph.edges) == {("IL6R", "JAK1")}
    bridged = build_omnipath_tf_graph(res)  # bridges on by default
    assert bridged.graph.has_edge("PLCG1", "PPP3CA") and bridged.graph.has_edge("ADCY1", "PRKACA")
    assert bridged.graph.edges["PLCG1", "PPP3CA"]["sources"] == "ChronoPack_second_messenger_IP3/Ca2+"
    assert not build_omnipath_tf_graph(res, second_messenger_bridges=False).graph.has_edge("PLCG1", "PPP3CA")

    # A node with two out-edges copying IL6R -> JAK1 (3 resources) and IL6ST -> JAK1 (1 resource).
    g = build_omnipath_tf_graph(res).graph
    g.add_edge("HUB", "IL6R", **g.edges["IL6R", "JAK1"])
    g.add_edge("HUB", "IL6ST", **g.edges["IL6ST", "JAK1"])
    t_uni, _, _, _ = build_signed_row_normalized_transitions(g, tfs={"STAT3", "FOXO1"})
    t_res, _, _, _ = build_signed_row_normalized_transitions(g, tfs={"STAT3", "FOXO1"}, edge_weight="resources")
    assert dict(t_uni["HUB"]) == {"IL6R": 0.5, "IL6ST": 0.5}
    assert dict(t_res["HUB"]) == {"IL6R": 0.75, "IL6ST": 0.25}


def _loop_diffusion_reference(nodes, t_pos, t_neg, *, tfs, seed, alpha=0.85, max_steps=80, min_steps=10, tol_mass=1e-12):
    """The original edge-by-edge implementation of signed_absorbing_diffusion_impulse."""
    from collections import defaultdict

    import numpy as np

    idx = {node: i for i, node in enumerate(nodes)}
    x_pos, x_neg = np.zeros(len(nodes)), np.zeros(len(nodes))
    for node, value in seed.items():
        x_pos[idx[node]] += value
    x_pos /= x_pos.sum()
    absorbed_pos, absorbed_neg, hist = defaultdict(float), defaultdict(float), []
    for step in range(max_steps):
        y_pos, y_neg = np.zeros_like(x_pos), np.zeros_like(x_neg)
        for source in nodes:
            i = idx[source]
            for target, prob in t_pos.get(source, []):
                y_pos[idx[target]] += x_pos[i] * prob
                y_neg[idx[target]] += x_neg[i] * prob
            for target, prob in t_neg.get(source, []):
                y_neg[idx[target]] += x_pos[i] * prob
                y_pos[idx[target]] += x_neg[i] * prob
        step_total = 0.0
        for tf in tfs:
            j = idx[tf]
            if y_pos[j]:
                absorbed_pos[tf] += y_pos[j]
            if y_neg[j]:
                absorbed_neg[tf] += y_neg[j]
            step_total += y_pos[j] + y_neg[j]
            y_pos[j] = y_neg[j] = 0.0
        x_pos, x_neg = alpha * y_pos, alpha * y_neg
        hist.append({"step": step, "absorbed_step": step_total,
                     "absorbed_total": sum(absorbed_pos.values()) + sum(absorbed_neg.values()),
                     "transient_mass": x_pos.sum() + x_neg.sum()})
        if step >= min_steps and x_pos.sum() + x_neg.sum() < tol_mass:
            break
    return dict(absorbed_pos), dict(absorbed_neg), x_pos, x_neg, pd.DataFrame(hist)


def test_sparse_diffusion_matches_edge_by_edge_loop():
    import networkx as nx
    import numpy as np

    from chronopack.tl.grn._omnipath import build_signed_row_normalized_transitions, signed_absorbing_diffusion_impulse

    rng = np.random.default_rng(0)
    g = nx.gnp_random_graph(60, 0.08, seed=1, directed=True)
    graph = nx.DiGraph()
    for u, v in g.edges:
        sign = rng.choice(["+", "-", "both", "none"])
        graph.add_edge(f"n{u}", f"n{v}", is_stimulation=sign in ("+", "both"), is_inhibition=sign in ("-", "both"))
    tfs = [f"n{i}" for i in range(0, 60, 6) if f"n{i}" in graph]
    t_pos, t_neg, _, _ = build_signed_row_normalized_transitions(graph, tfs=tfs)
    nodes = list(graph.nodes)
    seed = {nodes[1]: 0.7, nodes[2]: 0.3}
    new = signed_absorbing_diffusion_impulse(nodes, t_pos, t_neg, tfs=tfs, seed=seed, max_steps=80)
    ref = _loop_diffusion_reference(nodes, t_pos, t_neg, tfs=tfs, seed=seed)
    for a, b in zip(new[:2], ref[:2]):
        assert set(a) == set(b)
        np.testing.assert_allclose([a[k] for k in sorted(a)], [b[k] for k in sorted(b)], rtol=1e-12, atol=1e-15)
    np.testing.assert_allclose(new[2], ref[2], atol=1e-15)
    pd.testing.assert_frame_equal(new[4], ref[4], rtol=1e-12, atol=1e-15)


def _write_toy_inputs(tmp_path):
    res = _toy_tf_resources()
    snap = tmp_path / "snap"
    snap.mkdir()
    res.signaling.to_csv(snap / "signaling.csv.gz", index=False)
    res.transcriptional.to_csv(snap / "transcriptional.csv.gz", index=False)
    (tmp_path / "tfs.txt").write_text("STAT3\nFOXO1\n")
    res.ligand_receptor.to_csv(tmp_path / "lr.csv")
    return dict(human_tf_path=tmp_path / "tfs.txt", ligand_receptor_path=tmp_path / "lr.csv", omnipath_snapshot=snap)


def test_rank_regulators_expression_filter_and_graph_cache(tmp_path, monkeypatch):
    import anndata as ad
    import numpy as np
    import pytest
    import scipy.sparse as sp

    from chronopack.tl.grn import _inference

    counts = sp.csr_matrix(np.array([[1, 0, 2, 0], [3, 0, 0, 0], [1, 0, 1, 1]]))
    adata = ad.AnnData(X=np.zeros((3, 4)), var=pd.DataFrame(index=["IL6R", "IL6ST", "JAK1", "STAT3"]))
    adata.layers["counts"] = counts
    assert cp.tl.grn.expressed_genes(adata) == ["IL6R", "JAK1", "STAT3"]
    assert cp.tl.grn.expressed_genes(adata, min_cells=2) == ["IL6R", "JAK1"]
    assert cp.tl.grn.expressed_genes(adata, min_fraction=0.5) == ["IL6R", "JAK1"]

    kwargs = _write_toy_inputs(tmp_path)
    builds = []
    original = _inference.build_omnipath_tf_graph
    monkeypatch.setattr(_inference, "build_omnipath_tf_graph", lambda *a, **k: builds.append(1) or original(*a, **k))
    _inference._TF_GRAPH_CACHE.clear()

    full = cp.tl.grn.rank_regulators(source="L1", is_ligand=True, **kwargs)
    expressed = ["IL6R", "JAK1", "STAT3"]  # IL6ST and FOXO1 not expressed
    ranked = cp.tl.grn.rank_regulators(source="L1", is_ligand=True, expressed_genes=expressed, **kwargs)
    assert set(full["TF"]) == {"STAT3", "FOXO1"} and ranked["TF"].tolist() == ["STAT3"]
    assert ranked.iloc[0]["priority_act"] > 0
    assert len(builds) == 1  # the second call reused the cached graph
    with pytest.raises(ValueError, match="not in expressed_genes"):
        cp.tl.grn.rank_regulators(source="IL6ST", expressed_genes=expressed, **kwargs)
    # Scope: "tfs" keeps unexpressed signaling proteins (JAK1 here), "network" removes them.
    no_jak1 = ["IL6R", "IL6ST", "STAT3"]
    tfs_scope = cp.tl.grn.rank_regulators(source="L1", is_ligand=True, expressed_genes=no_jak1, **kwargs)
    assert tfs_scope["TF"].tolist() == ["STAT3"] and tfs_scope.iloc[0]["priority_act"] > 0
    with pytest.raises(ValueError, match="No signal"):
        cp.tl.grn.rank_regulators(source="L1", is_ligand=True, expressed_genes=no_jak1,
                                  expression_scope="network", **kwargs)
    with pytest.raises(ValueError, match="None of the receptors"):
        cp.tl.grn.rank_regulators(source="L1", is_ligand=True, expressed_genes=["JAK1", "STAT3"], **kwargs)
    _inference._TF_GRAPH_CACHE.clear()


def test_tied_signs_split_passes_mass_without_direction():
    import networkx as nx
    import pytest

    from chronopack.tl.grn._omnipath import build_signed_row_normalized_transitions, run_tf_diffusion_ranking

    # IKK -| IKBA, then IKBA -> NFKB with both flags (tied); IKK -> OTHER_TF activating.
    g = nx.DiGraph()
    g.add_edge("IKK", "IKBA", is_stimulation=False, is_inhibition=True)
    g.add_edge("IKBA", "NFKB", is_stimulation=True, is_inhibition=True)
    g.add_edge("IKK", "OTHER_TF", is_stimulation=True, is_inhibition=False)
    tfs = ["NFKB", "OTHER_TF"]
    t_pos, t_neg, _, _ = build_signed_row_normalized_transitions(g, tfs=tfs, tied_signs="split")
    assert t_pos["IKBA"] == [("NFKB", 0.5)] and t_neg["IKBA"] == [("NFKB", 0.5)]

    old, _ = run_tf_diffusion_ranking(g, source="IKK", tfs=tfs, tied_signs="stimulation")
    new, _ = run_tf_diffusion_ranking(g, source="IKK", tfs=tfs, tied_signs="split")
    old, new = old.set_index("TF"), new.set_index("TF")
    assert old.loc["NFKB", "net"] < 0  # the tied edge read as activating turns IKK's inhibition of IKBA into -NFKB
    assert new.loc["NFKB", "net"] == pytest.approx(0.0, abs=1e-15)
    assert new.loc["NFKB", "support"] == pytest.approx(old.loc["NFKB", "support"])  # mass still arrives
    with pytest.raises(ValueError, match="tied_signs"):
        build_signed_row_normalized_transitions(g, tfs=tfs, tied_signs="other")
