import numpy as np
import pandas as pd
import pytest


def _synthetic_gene_trajectory_adata():
    ad = pytest.importorskip("anndata")

    rows = ["g0_A", "g1_A", "g2_B", "g3_B"]
    adata = ad.AnnData(
        X=np.array(
            [
                [0.0, 1.0, 2.0],
                [0.0, 0.5, 1.5],
                [2.0, 1.0, 0.0],
                [1.5, 0.5, -0.5],
            ]
        ),
        obs=pd.DataFrame(index=pd.Index(rows, name="row_id")),
        var=pd.DataFrame({"time": [0.0, 10.0, 20.0]}, index=["t0", "t1", "t2"]),
    )
    adata.obs["gene"] = ["g0", "g1", "g2", "g3"]
    adata.obs["cell_type"] = ["A", "A", "B", "B"]
    adata.obs["waldStat_1"] = [9.0, 8.0, 7.0, 6.0]
    adata.obs["peak_log2fc"] = [2.0, 1.2, -1.1, -2.0]
    adata.obs["t_peak"] = [20.0, 30.0, 40.0, 50.0]
    adata.obs["broad_cluster"] = [
        "trend-up",
        "trend-up",
        "trend-down_1",
        "trend-down_1",
    ]
    adata.obs["cluster_fine_induced"] = ["1", "1", "-1", "-1"]
    adata.obs["cluster_fine_repressed"] = ["-1", "-1", "2", "2"]
    adata.obs["PC1"] = [1.0, 0.8, -1.0, -0.9]
    adata.obs["PC2"] = [0.0, 0.1, 0.0, -0.2]
    return adata


def test_gene_trajectory_plot_data_uses_trajectory_anndata():
    import chronopack as cp

    adata = _synthetic_gene_trajectory_adata()
    result = cp.pl.gene_trajectory(adata, ["g0_A", "g2_B"])

    assert result.output is None
    assert result.data["row_id"].nunique() == 2
    assert sorted(result.data["time"].unique().tolist()) == [0.0, 10.0, 20.0]
    assert result.data.loc[result.data["row_id"] == "g0_A", "value"].tolist() == [0.0, 1.0, 2.0]
    assert set(result.data["broad_cluster"]) == {"trend-up", "trend-down_1"}


def test_trajectory_features_reads_obs_columns():
    import chronopack as cp

    adata = _synthetic_gene_trajectory_adata()
    result = cp.pl.trajectory_features(adata, x="PC1", y="PC2", color="broad_cluster")

    assert result.data["PC1"].tolist() == [1.0, 0.8, -1.0, -0.9]
    assert result.data["PC2"].tolist() == [0.0, 0.1, 0.0, -0.2]
    assert result.data["broad_cluster"].tolist() == [
        "trend-up",
        "trend-up",
        "trend-down_1",
        "trend-down_1",
    ]


def test_clustering_metrics_plot_marks_recommended_k():
    import chronopack as cp
    from chronopack.tl import GeneTrajectoryClusteringResult

    result = GeneTrajectoryClusteringResult(
        metrics=pd.DataFrame(
            {
                "k": [3, 4, 5],
                "WSS": [10.0, 7.0, 6.0],
                "p90_centroid_distance": [1.0, 0.8, 0.7],
                "resolution": [0.0, 0.75, 1.0],
                "mean_silhouette": [0.2, 0.5, 0.4],
                "elbow_score": [0.0, 0.25, 0.0],
            }
        ),
        assignments={},
    )

    plot = cp.pl.clustering_metrics(result)
    assert plot.figure is not None
    assert plot.data["k"].tolist() == [3, 4, 5]


def test_missing_plot_key_raises_clear_error():
    import chronopack as cp

    adata = _synthetic_gene_trajectory_adata()
    with pytest.raises(KeyError, match="Trajectory observation column"):
        cp.pl.trajectory_features(adata, color="missing_column")


def test_r_backed_plot_rendering_smoke(tmp_path):
    pytest.importorskip("rpy2")
    import chronopack as cp

    adata = _synthetic_gene_trajectory_adata()
    out = tmp_path / "trajectory_features.svg"
    result = cp.pl.trajectory_features(adata, x="PC1", y="PC2", color="broad_cluster", output=out)

    assert result.output == str(out)
    assert out.exists()
    assert out.stat().st_size > 0
    assert result.preview_png[:8] == b"\x89PNG\r\n\x1a\n"


def _medoid_adata():
    ad = pytest.importorskip("anndata")

    t = np.linspace(0, 1, 8)
    curves = [np.sin(np.pi * t) + 0.05 * i for i in range(4)] + [-t ** 2 - 0.05 * i for i in range(4)]
    rows = [f"g{i}_{'A' if i % 2 else 'B'}" for i in range(8)]
    adata = ad.AnnData(
        X=np.array(curves),
        obs=pd.DataFrame(index=rows),
        var=pd.DataFrame({"time": t}, index=[f"t{i}" for i in range(t.size)]),
    )
    adata.obs["gene"] = [r.split("_")[0] for r in rows]
    adata.obs["waldStat_1"] = np.arange(8, 0, -1, dtype=float)
    adata.obs["broad_cluster"] = ["impulse-up_1"] * 4 + ["trend-down_1"] * 4
    adata.obsm["X_cluster_features"] = np.hstack([adata.X, np.diff(adata.X, axis=1)])
    return adata


def test_cluster_medoids_inputs_and_default_palette():
    import chronopack as cp
    from chronopack.pl._gene_trajectory_plots import CATEGORICAL_COLORS, MODE_COLORS, _default_colors

    adata = _medoid_adata()
    result = cp.pl.cluster_medoids_and_top_genes(adata, "broad_cluster")
    assert result.output is None
    assert result.data["cluster"].tolist() == ["impulse-up_1"] * 4 + ["trend-down_1"] * 4
    assert result.data["score"].tolist() == list(np.arange(8, 0, -1, dtype=float))
    assert _default_colors(adata.obs["broad_cluster"]) == MODE_COLORS
    assert _default_colors(["1", "2"])[:2] == list(CATEGORICAL_COLORS[:2])

    # Candidate labels (e.g. clustering assignments) as a Series; numbered labels stay numeric.
    labels = pd.Series([1, 1, 1, 1, 2, 2, 2, 2], index=adata.obs_names)
    numbered = cp.pl.cluster_medoids_and_top_genes(adata, labels, features=adata.obsm["X_cluster_features"])
    assert numbered.data["cluster"].tolist() == labels.tolist()

    del adata.obsm["X_cluster_features"]
    with pytest.raises(KeyError, match="Clustering features"):
        cp.pl.cluster_medoids_and_top_genes(adata, "broad_cluster")


def test_cluster_medoids_rendering_smoke(tmp_path):
    pytest.importorskip("rpy2")
    import chronopack as cp

    adata = _medoid_adata()
    out = tmp_path / "medoids.png"
    result = cp.pl.cluster_medoids_and_top_genes(
        adata, "broad_cluster", top_n=2, show_median_iqr=True, show_gene_titles=True,
        gene_title_mode="gene_sym", use_silhouette_alpha=False, title="", output=out,
    )
    assert out.exists() and out.stat().st_size > 0
    assert result.preview_png == out.read_bytes()


def test_trajectory_features_reads_obsm_components_and_skips_missing_rows():
    import chronopack as cp

    adata = _synthetic_gene_trajectory_adata()
    # Induced-only PCA: the two repressed rows have no coordinates.
    adata.obsm["X_pca_induced"] = np.array([[1.0, 2.0], [3.0, 4.0], [np.nan, np.nan], [np.nan, np.nan]])

    result = cp.pl.trajectory_features(adata, basis="X_pca_induced", x="PC1", y="PC2", color="broad_cluster")

    assert result.data["gene"].tolist() == ["g0", "g1"]
    assert result.data["PC1"].tolist() == [1.0, 3.0]
    assert result.data["PC2"].tolist() == [2.0, 4.0]


def test_plot_result_displays_preview_png_inline_with_short_repr(tmp_path):
    from chronopack.pl import RPlotResult

    svg = tmp_path / "plot.svg"
    svg.write_text("<svg xmlns='http://www.w3.org/2000/svg'></svg>")
    result = RPlotResult("demo", data=pd.DataFrame(), output=str(svg))

    assert result._repr_png_() is None  # SVG file and no preview: nothing to show inline
    png = b"\x89PNG\r\n\x1a\n" + b"\x00\x00\x00\rIHDR" + (2400).to_bytes(4, "big") + (1500).to_bytes(4, "big")
    with_preview = RPlotResult("demo", data=pd.DataFrame({"x": [1]}), output=str(svg), preview_png=png)
    assert with_preview._repr_png_() == (png, {"width": 1200})  # shown at half width
    assert repr(with_preview) == f"RPlotResult(plot_type='demo', output={str(svg)!r})"


def _association_table(n=40, seed=0):
    rng = np.random.default_rng(seed)
    wald = np.sort(rng.gamma(2.0, 20.0, n))[::-1]
    df = pd.DataFrame(
        {
            "waldStat_1": wald,
            "qval_1": np.linspace(0, 0.5, n) ** 2,
            "signedMaxLog2FC_start": rng.normal(0, 1.5, n),
        },
        index=[f"G{i}" for i in range(n)],
    )
    df.loc["G39", "signedMaxLog2FC_start"] = np.nan
    return df


def test_fdr_wald_break_is_the_wald_nearest_the_fdr_cutoff():
    import chronopack as cp

    df = _association_table()
    i = int(np.argmin(np.abs(df["qval_1"].to_numpy() - 0.01)))
    assert cp.pl.fdr_wald_break(df, 0.01) == df["waldStat_1"].iloc[i]

    result = cp.pl.association_volcano(df, fdr_thresh=0.01)
    assert result.output is None and result.metadata["y_break"] == df["waldStat_1"].iloc[i]
    assert cp.pl.association_volcano(df, y_break=12.5).metadata["y_break"] == 12.5
    with pytest.raises(KeyError, match="signedMaxLog2FC_start"):
        cp.pl.association_volcano(df.drop(columns="signedMaxLog2FC_start"))


def test_association_volcano_rendering_smoke(tmp_path):
    pytest.importorskip("rpy2")
    import chronopack as cp

    df = _association_table()
    out = tmp_path / "volcano.png"
    result = cp.pl.association_volcano(
        df, top_labels=5, x_lim=(-3, 3), y_lim=(1, 1000), y_lower_scale=0.1, output=out, dpi=100,
    )
    assert out.exists() and result.preview_png == out.read_bytes()
    labels = result.metadata["labels"]
    assert labels["gene"].tolist() == ["G0", "G1", "G2", "G3", "G4"]
    np.testing.assert_allclose(labels["waldStat_1"], df["waldStat_1"].iloc[:5])
    assert set(labels["side"]) <= {"other", "sig_mid", "neg", "pos"}

    pdf = cp.pl.association_volcano(df, y_break=None, output=tmp_path / "volcano.pdf")
    assert (tmp_path / "volcano.pdf").stat().st_size > 0 and pdf.preview_png[:4] == b"\x89PNG"


def test_gene_dynamics_curves_grid_follows_gene_order(tmp_path):
    import chronopack as cp

    t = np.linspace(0, 110, 5)
    curves = pd.concat(
        pd.DataFrame({"gene": g, "t": t, "mu_log1p": t / 100 + k, "lo_log1p": t / 100, "hi_log1p": t / 100 + 1, "lineage": lin})
        for g in ["A", "B", "C"] for k, lin in enumerate(["L1", "L2"])
    )
    out = tmp_path / "curves.pdf"
    fig = cp.pl.gene_dynamics_curves(curves, ["C", "A", "B"], ncols=2, colors={"L1": "#386cb0"}, output=out)

    visible = [ax for ax in fig.axes if ax.get_visible()]
    assert [ax.get_title(loc="left") for ax in visible] == ["C", "A", "B"]
    assert len(visible[0].lines) == 2 and visible[0].lines[0].get_color() == "#386cb0"
    assert out.exists()
    with pytest.raises(KeyError):
        cp.pl.gene_dynamics_curves(curves, ["A", "MISSING"])


def test_gene_dynamics_curves_draws_assigned_points_and_weighted_time_means():
    import chronopack as cp

    t = np.linspace(0, 110, 5)
    curves = pd.concat(
        pd.DataFrame({"gene": "A", "t": t, "mu_log1p": t / 100, "lineage": lin}) for lin in ["L1", "L2"]
    )
    observations = pd.concat([
        pd.DataFrame({"gene": "A", "t": [0, 0, 60], "expr_log1p": np.log1p([1.0, 3.0, 2.0]),
                      "weight": [1.0, 0.0, 0.5], "assigned": [True, False, True], "lineage": "L1"}),
        pd.DataFrame({"gene": "A", "t": [0, 0, 60], "expr_log1p": np.log1p([1.0, 3.0, 2.0]),
                      "weight": [0.0, 1.0, 0.5], "assigned": [False, True, False], "lineage": "L2"}),
    ])
    fig = cp.pl.gene_dynamics_curves(curves, ["A"], observations=observations, ncols=1,
                                     colors={"L1": "#386cb0", "L2": "#DC143C"})

    ax = fig.axes[0]
    points_l1, points_l2, means_l1, means_l2 = ax.collections[0], ax.collections[1], ax.collections[-2], ax.collections[-1]
    assert len(points_l1.get_offsets()) == 2 and len(points_l2.get_offsets()) == 1
    np.testing.assert_allclose(means_l1.get_offsets(), [[0, np.log1p(1.0)], [60, np.log1p(2.0)]])
    np.testing.assert_allclose(means_l2.get_offsets(), [[0, np.log1p(3.0)], [60, np.log1p(2.0)]])
    assert tuple(means_l1.get_facecolor()[0]) == (1.0, 1.0, 1.0, 1.0)


def _combined_table():
    df = pd.DataFrame(
        {
            "waldStat_1": [900.0, 50.0, 3.0, 40.0, 30.0, 20.0],
            "waldStat_2": [800.0, 2.0, 60.0, 0.0, 25.0, 1.0],
            "qval": [1e-6, 1e-4, 1e-3, 0.02, 0.001, 0.04],
            "log2FC_t1": [2.0, 1.5, 0.1, -1.2, 0.2, 0.3],
            "log2FC_t2": [2.5, 0.2, -2.0, 0.0, 0.1, 0.4],
        },
        index=["SHARED", "A_UP", "B_DOWN", "A_DOWN", "FLAT", "OTHER"],
    )
    return df


def test_response_specificity_volcano_checks_inputs():
    import chronopack as cp

    df = _combined_table()
    result = cp.pl.response_specificity_volcano(df, ("A", "B"))
    assert result.output is None and result.metadata["genes"] is None
    assert result.metadata["conditions"] == ("A", "B") and result.metadata["n_genes"] == 6
    assert cp.pl.response_specificity_volcano(df, {"condA": "X", "condB": "Y"}).metadata["conditions"] == ("X", "Y")
    with pytest.raises(KeyError, match="log2FC_t2"):
        cp.pl.response_specificity_volcano(df.drop(columns="log2FC_t2"), ("A", "B"))


def test_response_specificity_volcano_rendering_classifies_genes(tmp_path):
    pytest.importorskip("rpy2")
    import chronopack as cp

    out = tmp_path / "specificity.png"
    result = cp.pl.response_specificity_volcano(
        _combined_table(), ("A", "B"), require_fdr=False, fdr_thresh=0.01,
        dyn_thresh_condA=1, flat_thresh_condA=0.5, dyn_thresh_condB=1, flat_thresh_condB=0.5,
        label_top_wald=1, label_top_condA=1, label_top_condB=1, output=out, dpi=72,
    )
    assert out.exists() and result.preview_png == out.read_bytes()
    genes = result.metadata["genes"].set_index("gene")
    assert genes["side"].to_dict() == {
        "SHARED": "sig_mid", "A_UP": "condA_ind", "B_DOWN": "condB_rep",
        "A_DOWN": "condA_rep", "FLAT": "sig_mid", "OTHER": "other",
    }
    assert sorted(genes.index[genes["labelled"]]) == ["A_UP", "B_DOWN", "SHARED"]
    np.testing.assert_allclose(genes.loc["A_UP", "x"], 1.5 - 0.2)

    pdf = cp.pl.response_specificity_volcano(_combined_table(), ("A", "B"), output=tmp_path / "s.pdf")
    assert (tmp_path / "s.pdf").stat().st_size > 0 and pdf.preview_png[:4] == b"\x89PNG"
