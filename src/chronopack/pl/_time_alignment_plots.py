"""Diagnostics for :func:`chronopack.tl.estimate_time_alignment`."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ._gene_trajectory_plots import CATEGORICAL_COLORS


def _logistic(t, b, a, m, w):
    return b + a / (1.0 + np.exp(-(t - m) / w))


def time_alignment(result, *, reporters=None, contexts=None, colors=None, panel_size=(3.4, 2.5),
                   output: str | Path | None = None):
    """One row per reporter: each context on its own clock, then all contexts on the reference clock.

    Own-clock panels: per-tag values on the fit scale with 95% intervals (filled: inside the fit window); the
    registration fit (solid: the reference's per-tag curve read at scale * t + offset with the context's baseline and
    amplitude, band = the template's 95% uncertainty; on the reference, its own curve); the logistic rise fit (dashed),
    its midpoint with the bootstrap interval, and the rise fit's resolved / flag status.
    Overlay panels: every context's per-tag data mapped by the estimated alignment (and by the current one, if
    ``current`` was given), each scaled by its own rise fit (0 = baseline b, 1 = b + a; faint line: that fit), so
    only timing is compared.
    """
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    from chronopack.tl._time_alignment import _curve, _curve_eval, _irls

    fits, tags, al = result.fits, result.tags, result.alignment.set_index("context")
    reporters = list(dict.fromkeys(fits["gene"])) if reporters is None else list(reporters)
    contexts = list(al.index) if contexts is None else list(contexts)
    colors = dict(colors or {c: CATEGORICAL_COLORS[i % len(CATEGORICAL_COLORS)] for i, c in enumerate(contexts)})
    has_current = "current_scale" in al and al["current_scale"].notna().any()
    overlays = [("estimate", "scale", "offset")] + ([("current", "current_scale", "current_offset")]
                                                    if has_current else [])
    ncol = len(contexts) + len(overlays)
    fig, axes = plt.subplots(len(reporters), ncol, figsize=(panel_size[0] * ncol, panel_size[1] * len(reporters)),
                             squeeze=False, constrained_layout=True)
    F = fits.set_index(["context", "gene"])
    ref = result.parameters.get("reference", contexts[0])
    isd = result.parameters.get("interp_sd", 0.125)

    def _tagdata(c, g, window):
        d = tags[(tags["context"] == c) & (tags["gene"] == g)]
        d = d[d["in_window"]] if window else d
        return d["minutes"].to_numpy(float), d["y"].to_numpy(float), d["y_se"].to_numpy(float)

    for i, g in enumerate(reporters):
        ref_curve = _curve(*_tagdata(ref, g, False)) if (ref, g) in F.index else None
        for k, c in enumerate(contexts):
            ax = axes[i, k]
            if (c, g) not in F.index:
                ax.set_visible(False)
                continue
            f = F.loc[(c, g)]
            d = tags[(tags["context"] == c) & (tags["gene"] == g)]
            for inw, face in ((True, colors[c]), (False, "white")):
                dd = d[d["in_window"] == inw]
                ax.errorbar(dd["minutes"], dd["y"], yerr=1.96 * dd["y_se"], fmt="o", ms=3.5, mfc=face,
                            mec=colors[c], color=colors[c], elinewidth=0.8, capsize=0)
            tt = np.linspace(d["minutes"].min(), d["minutes"].max(), 400)
            s_c, o_c = (1.0, 0.0) if c == ref else (al.loc[c, "scale"], al.loc[c, "offset"]) if c in al.index \
                else (np.nan, np.nan)
            if ref_curve is not None and np.isfinite(s_c) and np.isfinite(o_c):
                if c == ref:
                    b_c, a_c = 0.0, 1.0
                else:
                    tw, yw, sw = _tagdata(c, g, True)
                    ok = np.isfinite(yw) & np.isfinite(sw) & (sw > 0)
                    H, V = _curve_eval(ref_curve, s_c * tw[ok] + o_c, isd)
                    b_c, a_c, _ = _irls(H, V, yw[ok], sw[ok])
                Hg, Vg = _curve_eval(ref_curve, s_c * tt + o_c, isd)
                band = 1.96 * abs(float(a_c)) * np.sqrt(Vg)
                yy = float(b_c) + float(a_c) * Hg
                ax.fill_between(tt, yy - band, yy + band, color=colors[c], alpha=0.15, lw=0)
                ax.plot(tt, yy, color=colors[c], lw=1.8)
            if np.isfinite(f["w"]):
                inw = tt <= f["window_end_min"]
                yy = _logistic(tt, f["b"], f["a"], f["m"], f["w"])
                ax.plot(tt[inw], yy[inw], color=colors[c], lw=1.0, ls="--")
                ax.errorbar([f["m"]], [f["b"] + f["a"] / 2], xerr=[[f["m"] - f["m_lo"]], [f["m_hi"] - f["m"]]],
                            fmt="D", ms=4, mfc="white", mec=colors[c], color=colors[c], capsize=2.5)
            status = "resolved" if f["resolved"] else (f["flags"] or "not resolved")
            ax.set_title(f"{g}, {c}: w {f['w']:.1f} ({f['w_lo']:.1f}-{f['w_hi']:.1f})\n{status}", fontsize=8)
            ax.set_xlabel(f"{c} own minutes", fontsize=7.5)
            ax.set_ylabel(f"{'log ' if f['fit_scale'] == 'log' else ''}mean count / cell", fontsize=7.5)
            ax.tick_params(labelsize=7)
        for j, (name, sk, ok) in enumerate(overlays):
            ax = axes[i, len(contexts) + j]
            for c in contexts:
                if (c, g) not in F.index or c not in al.index:
                    continue
                s, o = (1.0, 0.0) if c == ref else (al.loc[c, sk], al.loc[c, ok])
                f = F.loc[(c, g)]
                if not (np.isfinite(s) and np.isfinite(o)):
                    continue
                d = tags[(tags["context"] == c) & (tags["gene"] == g)]
                logi = np.isfinite(f["a"]) and f["a"] != 0 and np.isfinite(f["w"])
                if logi:
                    b, a = f["b"], f["a"]
                else:  # scale by the window data: first tag -> 0, largest change -> 1
                    yw = d.loc[d["in_window"], "y"].to_numpy()
                    b = yw[0]
                    a = (yw[np.nanargmax(np.abs(yw - b))] - b) or 1.0
                x = s * d["minutes"].to_numpy() + o
                ax.errorbar(x, (d["y"] - b) / a, yerr=1.96 * d["y_se"] / abs(a), fmt="o", ms=3,
                            color=colors[c], elinewidth=0.7, capsize=0)
                if logi:
                    tt = np.linspace(d["minutes"].min(), f["window_end_min"], 400)
                    ax.plot(s * tt + o, (_logistic(tt, f["b"], f["a"], f["m"], f["w"]) - b) / a, color=colors[c],
                            lw=1.0, alpha=0.6)
                else:
                    ax.plot(x, (d["y"] - b) / a, color=colors[c], lw=1.0, alpha=0.6)
            ax.axhline(0, color="#cccccc", lw=0.6)
            ax.axhline(1, color="#cccccc", lw=0.6)
            rt = tags.loc[(tags["context"] == ref) & (tags["gene"] == g), "minutes"]
            if len(rt):
                pad = 0.1 * (rt.max() - rt.min())
                ax.set_xlim(rt.min() - pad, rt.max() + pad)
            ax.set_title(f"{g}: {name} alignment", fontsize=8)
            ax.set_xlabel(f"aligned minutes ({ref} clock)", fontsize=7.5)
            ax.set_ylabel("fraction of the rise", fontsize=7.5)
            ax.tick_params(labelsize=7)
    handles = [Line2D([], [], color=colors[c], marker="o", ms=3.5, lw=1.4, label=c) for c in contexts]
    fig.legend(handles=handles, loc="outside upper center", ncol=len(contexts), fontsize=8, frameon=False)
    if output is not None:
        fig.savefig(output, dpi=150, bbox_inches="tight")
    return fig
