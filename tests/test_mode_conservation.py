import numpy as np
import pandas as pd
import pytest


def _context_df(context, rows):
    df = pd.DataFrame(rows)
    df["cell_type"] = context
    df.index = [f"{g}_{context}" for g in df["gene"]]
    return df


def _paired_tables():
    source = _context_df(
        "source",
        [
            {"gene": "a", "broad_cluster": "trend-up", "qval_1": 0.001, "waldStat_1": 100},
            {"gene": "b", "broad_cluster": "trend-up", "qval_1": 0.002, "waldStat_1": 90},
            {"gene": "c", "broad_cluster": "trend-up", "qval_1": 0.003, "waldStat_1": 80},
            {"gene": "d", "broad_cluster": "trend-down_1", "qval_1": 0.001, "waldStat_1": 70},
            {"gene": "e", "broad_cluster": "trend-down_1", "qval_1": 0.001, "waldStat_1": 60},
            {"gene": "f", "broad_cluster": "trend-up", "qval_1": 0.8, "waldStat_1": 50},
        ],
    )
    target = _context_df(
        "target",
        [
            {"gene": "a", "broad_cluster": "trend-up", "qval_1": 0.001, "waldStat_1": 85},
            {"gene": "b", "broad_cluster": "trend-down_1", "qval_1": 0.001, "waldStat_1": 95},
            {"gene": "c", "broad_cluster": "trend-up", "qval_1": 0.2, "waldStat_1": 75},
            {"gene": "d", "broad_cluster": "trend-down_1", "qval_1": 0.001, "waldStat_1": 65},
            {"gene": "g", "broad_cluster": "trend-up", "qval_1": 0.001, "waldStat_1": 60},
        ],
    )
    return source, target


def test_mode_conservation_preserves_directional_topk_definition():
    import chronopack as cp

    source, target = _paired_tables()
    result = cp.tl.mode_conservation(
        source,
        target,
        source_context="source",
        target_context="target",
        top_k=2,
    )

    s2t_up = result.directional_scores[
        (result.directional_scores["source_context"] == "source")
        & (result.directional_scores["mode"] == "trend-up")
    ].iloc[0]
    t2s_up = result.directional_scores[
        (result.directional_scores["source_context"] == "target")
        & (result.directional_scores["mode"] == "trend-up")
    ].iloc[0]

    assert s2t_up["source_n"] == 2
    assert s2t_up["same_mode_recovered_n"] == 1
    assert s2t_up["directional_recovery"] == 0.5
    assert t2s_up["source_n"] == 2
    assert t2s_up["same_mode_recovered_n"] == 1
    assert t2s_up["directional_recovery"] == 0.5

    up = result.conservation_scores[result.conservation_scores["mode"] == "trend-up"].iloc[0]
    assert up["bidirectional_conservation"] == 0.5


def test_mode_switchers_classify_same_switch_low_and_missing():
    import chronopack as cp

    source, target = _paired_tables()
    result = cp.tl.mode_switchers(
        source,
        target,
        source_context="source",
        target_context="target",
        top_k=3,
        treat_missing_as_low_significance=False,
    )
    by_gene = result.switchers.set_index("gene")

    assert by_gene.loc["a", "switch_class"] == "same_mode"
    assert by_gene.loc["b", "switch_class"] == "mode_switch"
    assert by_gene.loc["c", "switch_class"] == "becomes_low_significance"
    assert by_gene.loc["g", "switch_class"] == "missing_assignment"


def test_missing_can_be_imputed_as_low_significance():
    import chronopack as cp

    source, target = _paired_tables()
    result = cp.tl.mode_switchers(
        source,
        target,
        source_context="source",
        target_context="target",
        top_k=3,
        treat_missing_as_low_significance=True,
    )
    by_gene = result.switchers.set_index("gene")

    assert by_gene.loc["g", "source_mode"] == "low_significance"
    assert by_gene.loc["g", "switch_class"] == "enters_dynamic_mode"


def test_ann_data_inputs_use_standardized_var_columns():
    ad = pytest.importorskip("anndata")
    import chronopack as cp

    source, target = _paired_tables()
    cluster_key = "gene_trajectory_clusters"
    fit_key = "gene_trajectories"
    source_ad = ad.AnnData(X=np.zeros((0, len(source))), var=pd.DataFrame(index=source.index))
    target_ad = ad.AnnData(X=np.zeros((0, len(target))), var=pd.DataFrame(index=target.index))
    for adata, df in [(source_ad, source), (target_ad, target)]:
        adata.var[f"{fit_key}_gene"] = df["gene"].to_numpy()
        adata.var[f"{fit_key}_cell_type"] = df["cell_type"].to_numpy()
        adata.var[f"{fit_key}_qval_1"] = df["qval_1"].to_numpy()
        adata.var[f"{fit_key}_waldStat_1"] = df["waldStat_1"].to_numpy()
        adata.var[f"{cluster_key}_broad_cluster"] = df["broad_cluster"].to_numpy()

    result = cp.tl.mode_conservation(source_ad, target_ad, source_context="source", target_context="target", top_k=2)
    assert not result.conservation_scores.empty
    assert set(result.mode_transitions.columns) == {"source_mode", "target_mode", "n"}


def test_mode_sankey_plot_data_options():
    import chronopack as cp

    source, target = _paired_tables()
    result = cp.tl.mode_switchers(source, target, source_context="source", target_context="target", top_k=3)
    shown = cp.pl.mode_sankey(result)
    dropped = cp.pl.mode_sankey(result, drop_low_significance=True)

    assert shown.data["n"].sum() >= dropped.data["n"].sum()
    assert "low_significance" not in set(dropped.data["source_mode"]) | set(dropped.data["target_mode"])


def test_mode_conservation_render_smoke(tmp_path):
    pytest.importorskip("rpy2")
    import chronopack as cp

    source, target = _paired_tables()
    result = cp.tl.mode_conservation(source, target, source_context="source", target_context="target", top_k=2)
    out = tmp_path / "mode_conservation.svg"
    plot = cp.pl.mode_conservation(result, output=out)

    assert plot.output == str(out)
    assert out.exists()
    assert out.stat().st_size > 0
