from pathlib import Path

import numpy as np
import pandas as pd
import pytest


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
