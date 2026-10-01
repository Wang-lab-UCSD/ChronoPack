from pathlib import Path

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


MODES = ["impulse-up_1", "impulse-up_2", "trend-up", "transient-down", "trend-down_1", "trend-down_2"]


def _three_context_adata(n_genes=24):
    ad = pytest.importorskip("anndata")

    rows, obs = [], []
    for c_idx, context in enumerate(["K562a", "K562b", "THP1"]):
        for g in range(n_genes):
            mode = MODES[(g + c_idx * (g % 2)) % 6]
            rows.append(f"G{g}_{context}")
            obs.append({
                "gene": f"G{g}",
                "cell_type": context,
                "broad_cluster": mode,
                "waldStat_1": float(100 - g + c_idx),
                "qval_1": 0.001 if g % 5 else 0.2,
                "t_peak": float(10 + g),
                "t_onset": float(g),
            })
    obs = pd.DataFrame(obs, index=rows)
    return ad.AnnData(X=np.zeros((len(rows), 3)), obs=obs)


def test_mode_proportions_counts_every_mode_per_context():
    import chronopack as cp

    adata = _three_context_adata()
    result = cp.pl.mode_proportions(adata)
    assert result.output is None
    assert len(result.data) == 3 * 6
    assert result.data.groupby("cell_type")["prop"].sum().round(12).eq(1).all()
    k562a = result.data[result.data["cell_type"] == "K562a"].set_index("broad_cluster")["n"]
    assert k562a.to_dict() == {m: 4 for m in MODES}


def test_topk_mode_recovery_matches_definition():
    pytest.importorskip("rpy2")
    import chronopack as cp
    from chronopack.tl._mode_conservation import MODE_DIRECTIONS

    adata = _three_context_adata()
    out = cp.tl.topk_mode_recovery(adata, [("K562a", "K562b")], k=5)
    obs = adata.obs.assign(direction=adata.obs["broad_cluster"].map(MODE_DIRECTIONS))
    for _, row in out.iterrows():
        src = obs[(obs.cell_type == row.source) & (obs.direction == row["mode"]) & (obs.qval_1 < 0.05)]
        src = src.sort_values("waldStat_1", ascending=False).head(5)
        tgt = obs[(obs.cell_type == row.target) & (obs.direction == row["mode"]) & (obs.qval_1 < 0.05)]
        assert row.recovery == pytest.approx(src["gene"].isin(tgt["gene"]).mean())
    assert out[["mode", "source", "target", "k"]].values.tolist() == [
        ["induced", "K562a", "K562b", 5],
        ["repressed", "K562a", "K562b", 5],
    ]


def test_cross_context_plots_render(tmp_path):
    pytest.importorskip("rpy2")
    import chronopack as cp

    adata = _three_context_adata()
    sankey = cp.pl.mode_sankey(adata, [("K562a", "K562b"), ("K562a", "THP1")], top_k=10, output=tmp_path / "sankey.svg")
    assert set(sankey.data["source_cell_type"]) == {"K562a"}
    chord = cp.pl.mode_chord(adata, ["K562a", "K562b", "THP1"], top_k=10, min_flow=1, output=tmp_path / "chord.png")
    assert {"cluster_1", "cluster_2", "n"} <= set(chord.data.columns)
    pie = cp.pl.mode_proportions(adata, output=tmp_path / "pie.svg")
    adata.obs["mode2"] = adata.obs["broad_cluster"].replace({"impulse-up_1": "impulse-up", "impulse-up_2": "impulse-up"})
    violin = cp.pl.cluster_violin(
        adata, ["t_peak", "t_onset"], cluster_use="impulse-up", cluster_col="mode2", top_k=5, output=tmp_path / "violin.svg"
    )
    assert violin.data["feature"].tolist().count("t_peak") == 3
    for result in (sankey, chord, pie, violin):
        assert Path(result.output).stat().st_size > 0
        assert result.preview_png[:4] == b"\x89PNG"
