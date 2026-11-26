"""Plots for the composition check and the response-heterogeneity analysis."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ._gene_trajectory_plots import CATEGORICAL_COLORS

_REF = "#6b6b6b"


def _groups(df):
    return list(pd.unique(df["group"])) if "group" in df else [None]


def composition_qc(result, *, time_label="Time", panel_size=(4.2, 3.4), output: str | Path | None = None):
    """Two panels per :func:`chronopack.tl.composition_qc` result, one line per group.

    Left: share of UMIs from induced genes. Right: the composition shift, i.e. the apparent log2 change of an
    unchanged gene under total-count normalisation (what ``normalization="tmm"`` removes): per time point (points)
    and its 3-point trend (line), with the drop implied by the induced share alone (dotted) and the +-0.1 log2 band
    below which composition is negligible.
    """
    import matplotlib.pyplot as plt

    df = result.by_time
    fig, axes = plt.subplots(1, 2, figsize=(2 * panel_size[0], panel_size[1]), constrained_layout=True)
    for i, g in enumerate(_groups(df)):
        d = df if g is None else df[df["group"] == g]
        col = CATEGORICAL_COLORS[i % len(CATEGORICAL_COLORS)]
        label = "all cells" if g is None else str(g)
        axes[0].plot(d["time"], 100 * d["induced_share"], "o-", color=col, lw=2, ms=5, label=label)
        axes[1].plot(d["time"], d["composition_trend_log2"], "-", color=col, lw=2, label=label)
        axes[1].plot(d["time"], d["composition_shift_log2"], "o", color=col, ms=4, alpha=0.6)
        axes[1].plot(d["time"], d["induced_share_bound_log2"], ":", color=col, lw=1.5)
    axes[1].axhspan(-0.1, 0.1, color="0.92", zorder=0)
    axes[1].axhline(0, color="0.6", lw=0.8)
    axes[0].set(xlabel=time_label, ylabel="UMIs from induced genes (%)", title="Induced genes' share of libraries")
    axes[1].set(xlabel=time_label, ylabel="log2 shift of an unchanged gene",
                title="Composition shift under total-count normalisation")
    axes[1].plot([], [], ":", color=_REF, label="implied by induced share")
    axes[1].legend(frameon=False, fontsize=8)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
    if output is not None:
        fig.savefig(output, dpi=200, bbox_inches="tight")
    return fig


def response_heterogeneity(
    result,
    *,
    adata=None,
    time_key="time",
    groups=None,
    time_label="Time",
    panel_size=(3.6, 3.2),
    output: str | Path | None = None,
):
    """Summary of a :func:`chronopack.tl.response_heterogeneity` result.

    Panels:

    - with ``adata``: the response score per time point. Cells are split by the result's ``groupby`` (one outline
      per group, same scale) when it had one, so it shows which group the response comes from; otherwise all cells.
    - response size: ``module_log2fc`` per group over time (the pooled module change, share responding x size).
    - share of non-responders per group with its 95% interval (hollow markers where the call is "indeterminate" or
      "no response").
    - one panel per group: share of cells still at baseline (points) against the graded expectation (dashed), the
      share an equal response in every cell would leave there through count noise alone.
    """
    import matplotlib.pyplot as plt

    df = result.by_time
    groups = list(pd.unique(df["group"])) if groups is None else list(groups)
    colors = {g: CATEGORICAL_COLORS[i % len(CATEGORICAL_COLORS)] for i, g in enumerate(groups)}
    n = len(groups) + 2 + (adata is not None)
    fig, axes = plt.subplots(1, n, figsize=(panel_size[0] * n, panel_size[1]), constrained_layout=True, squeeze=False)
    axes = list(axes[0])
    if adata is not None:
        ax = axes.pop(0)
        score = result.score.reindex(adata.obs_names).to_numpy(float)
        t = pd.to_numeric(adata.obs[time_key].astype(str)).to_numpy(float)
        times = np.unique(t)
        lo, hi = np.nanquantile(score, [0.005, 0.98])  # strong responders have a long tail; keep the bulk readable
        bins = np.linspace(lo, hi, 50)
        groupby = (result.parameters or {}).get("groupby")
        split = [g for g in groups if g != "all"] if groupby is not None and groupby in adata.obs else []
        if split:
            lab = adata.obs[groupby].astype(str).to_numpy()
            for j, tt in enumerate(times):
                hs = {g: np.histogram(score[(t == tt) & (lab == g)], bins=bins)[0].astype(float) for g in split}
                top = max(max(h.max() for h in hs.values()), 1.0)
                for g in split:
                    h = hs[g] / top * 0.9
                    ax.fill_between(bins[:-1], j, j + h, step="post", color=colors[g], alpha=0.25, lw=0)
                    ax.step(bins[:-1], j + h, where="post", color=colors[g], lw=1.1,
                            label=str(g) if j == 0 else None)
            ax.legend(frameon=False, fontsize=7, loc="lower right")
            title = f"Response score by {groupby}"
        else:
            cmap = plt.get_cmap("viridis")
            for j, tt in enumerate(times):
                h, e = np.histogram(score[t == tt], bins=bins, density=True)
                h = h / max(h.max(), 1e-12) * 0.9
                ax.fill_between(e[:-1], j, j + h, step="post", color=cmap(j / max(len(times) - 1, 1)), alpha=0.85,
                                lw=0)
            title = "Response score"
        ax.set_yticks(np.arange(len(times)) + 0.3, [f"{x:g}" for x in times], fontsize=8)
        ax.axvline(0, color="0.5", lw=0.8)
        ax.axvline(1, color="0.5", lw=0.8, ls="--")
        ax.set(xlabel="response score (0 = baseline, 1 = peak)", ylabel=time_label,
               title=f"{title} ({len(result.module)}-gene module)")
        ax.spines[["top", "right", "left"]].set_visible(False)

    ax = axes.pop(0)
    for g in groups:
        d = df[df["group"] == g].sort_values("time")
        ax.plot(d["time"], d["module_log2fc"], "o-", color=colors[g], lw=1.5, ms=4, label=str(g))
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set(xlabel=time_label, ylabel="module log2 fold over baseline", title="Response size (pooled)")
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes.pop(0)
    for g in groups:
        d = df[df["group"] == g].sort_values("time")
        solid = d["call"].isin(["all respond", "non-responders"]).to_numpy()
        yerr = [d["nonresponder_share"] - d["nonresponder_lo"], d["nonresponder_hi"] - d["nonresponder_share"]]
        ax.errorbar(d["time"], d["nonresponder_share"], yerr=yerr, fmt="none", ecolor=colors[g], lw=1, capsize=2)
        ax.plot(d["time"], d["nonresponder_share"], "-", color=colors[g], lw=1.5, label=str(g))
        ax.scatter(d["time"][solid], d["nonresponder_share"][solid], s=28, color=colors[g], zorder=3)
        ax.scatter(d["time"][~solid], d["nonresponder_share"][~solid], s=28, facecolor="white",
                   edgecolor=colors[g], zorder=3)
    ax.set_ylim(-0.02, 1.02)
    ax.set(xlabel=time_label, ylabel="share of non-responders", title="Non-responders beyond count noise")
    ax.spines[["top", "right"]].set_visible(False)

    for i, g in enumerate(groups):
        ax = axes[i]
        d = df[df["group"] == g].sort_values("time")
        ax.axhline(1, color="0.85", lw=0.8)
        ax.plot(d["time"], d["graded_expectation"], "--", color=_REF, lw=1.3, label="every cell responds equally")
        ax.errorbar(d["time"], d["share_at_baseline"],
                    yerr=[d["share_at_baseline"] - d["share_lo"], d["share_hi"] - d["share_at_baseline"]],
                    fmt="o", color=colors[g], ms=5, capsize=2, lw=1.2, label="data")
        ax.set_ylim(0, max(1.2, float(d["share_hi"].max()) + 0.1) if len(d) else 1.2)
        ax.set(xlabel=time_label, title=f"{g}: cells still at baseline")
        if i == 0:
            ax.set_ylabel("share below baseline quantile / quantile")
            ax.legend(frameon=False, fontsize=7, loc="upper right")
        ax.spines[["top", "right"]].set_visible(False)
    if output is not None:
        fig.savefig(output, dpi=200, bbox_inches="tight")
    return fig


__all__ = ["composition_qc", "response_heterogeneity"]
