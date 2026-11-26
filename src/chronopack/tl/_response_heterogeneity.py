"""Response heterogeneity: what share of cells respond to the stimulus?

A pooled NB-GAM curve is the mean over cells, so at every time point it mixes the share of cells responding with how
strongly they respond. This module scores every cell on a module of induced genes and asks, per time point (and per
cell group), how many cells are still at the baseline level, compared with what count noise alone would leave there
if every cell responded equally. Method: ``verification_runs/responder_fraction_20261104`` (K562 + TNF-alpha).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import optimize, special

from chronopack._core import record_provenance, validate_adata
from chronopack._core._keys import GENE_DYNAMICS_KEY
from chronopack._core._validation import ChronoPackError


@dataclass
class ResponseHeterogeneityResult:
    """Output of :func:`response_heterogeneity`.

    ``by_time``: one row per group x time point after the baseline (see the function for the columns).
    ``module``: the genes scored, with their baseline mean, peak mean and score weight. ``score``: the per-cell
    response score (also in ``adata.obs[score_key]``).
    """

    by_time: pd.DataFrame
    module: pd.DataFrame
    score: pd.Series
    peak_time: float
    parameters: dict = field(default_factory=dict)


def _nb_ll(c, mu, r):
    return (special.gammaln(c + r) - special.gammaln(r) - special.gammaln(c + 1)
            + r * np.log(r / (r + mu)) + c * np.log(mu / (r + mu)))


def _fit_nb(c, sf):
    """NB maximum likelihood of a mean scale m (cell mean sf_i * m) and size r."""
    def f(th):
        m, r = np.exp(th)
        mu = sf * m
        g_m = np.sum(r * (c - mu) / (r + mu))
        g_r = r * np.sum(special.digamma(c + r) - special.digamma(r) + np.log(r / (r + mu)) + (mu - c) / (r + mu))
        return -np.sum(_nb_ll(c, mu, r)), -np.array([g_m, g_r])

    m0 = max(c.sum() / sf.sum(), 1e-4)
    o = optimize.minimize(f, np.log([m0, 1.0]), jac=True, method="L-BFGS-B",
                          bounds=[(np.log(1e-5), np.log(1e3)), (np.log(1e-2), np.log(1e4))])
    return np.exp(o.x)


def _tail_ratio(sc_t, sc_0, q):
    """Share of cells at t below the baseline q-quantile, over q: ~1 if nobody moved, -> 0 if everybody did."""
    return np.mean(sc_t <= np.quantile(sc_0, q)) / q


def _dense(X):
    return X.toarray() if sp.issparse(X) else np.asarray(X, dtype=float)


def _select_module(X, var_names, time, times, b, s, n_genes, min_mean_counts, min_log2fc, pvalues, p_threshold):
    """Induced genes with enough counts, ranked by how well they separate the peak from the baseline."""
    m = np.vstack([np.asarray(X[time == t].sum(axis=0)).ravel() / s[time == t].sum() for t in times])
    c = np.vstack([np.asarray(X[time == t].sum(axis=0)).ravel() for t in times])
    S = np.array([s[time == t].sum() for t in times])
    m0 = np.maximum(m[b], 0.5 / S[b])                               # half a count over the baseline cells
    others = np.flatnonzero(np.arange(len(times)) != b)
    # peaks of 3-point running means, so that one noisy time point does not select a gene
    m = m.copy()
    m[others] = pd.DataFrame(m[others]).rolling(3, center=True, min_periods=1).mean().to_numpy()
    fold = m[others] / m0
    k = others[np.argmax(fold, axis=0)]                             # each gene's peak time point
    m_pk, c_pk, S_pk = m[k, np.arange(m.shape[1])], c[k, np.arange(m.shape[1])], S[k]
    with np.errstate(divide="ignore", invalid="ignore"):
        stat = np.where(m_pk > m0, c_pk * np.log(m_pk / m0) - S_pk * (m_pk - m0), 0.0)
        peak_counts = m_pk                                           # mean counts per cell of median library size
        ok = (peak_counts >= min_mean_counts) & (np.log2(m_pk / m0) >= min_log2fc)
    if pvalues is not None:
        ok &= np.asarray(pvalues, dtype=float) < p_threshold
    order = np.argsort(-np.where(ok, stat, -np.inf))[: int(ok.sum())][:n_genes]
    return list(np.asarray(var_names)[order])


def response_heterogeneity(
    adata,
    genes=None,
    *,
    time_key="time",
    layer="counts",
    groupby=None,
    library_size_key="total_counts",
    baseline=None,
    n_genes=25,
    min_mean_counts=0.3,
    min_log2fc=1.0,
    dynamics_key=GENE_DYNAMICS_KEY,
    p_threshold=0.05,
    quantile=0.25,
    max_nonresponder=0.1,
    n_boot=200,
    n_sim=50,
    random_state=0,
    score_key="response_score",
    key_added="response_heterogeneity",
):
    """Share of cells that respond at each time point, beyond what count noise hides.

    **Module.** ``genes``, or up to ``n_genes`` induced genes chosen automatically from their profiles over time
    (3-point running means after the baseline): >= ``min_mean_counts`` counts per cell (at the median library size)
    at their peak, peak fold >= ``min_log2fc`` (log2) over the baseline, ranked by
    the Poisson log-likelihood ratio of peak vs baseline (large and well-measured changes first). If
    ``adata.var[f"{dynamics_key}_pvalue_1"]`` (or ``_pvalue``) exists, only genes with p < ``p_threshold`` are
    considered. Pass ``genes`` to score a specific programme (e.g. NF-kB targets, interferon-stimulated genes).

    **Score.** Per cell, the Poisson log-likelihood ratio of its module counts under the peak profile vs the baseline
    profile (the optimal linear detector of the response), rescaled so that the median baseline cell is 0 and the
    median cell at the peak time point is 1. Written to ``adata.obs[score_key]``.

    **Share still at baseline.** For each time point t and group (``"all"`` plus each ``groupby`` value, each against
    its own baseline cells), ``share_at_baseline`` = (share of cells below the ``quantile`` of the baseline scores) /
    ``quantile``: ~1 when no cell has moved, ~0 when every cell has moved well clear of the baseline. 95% bootstrap
    interval ``share_lo`` / ``share_hi``.

    **Graded expectation.** Count noise leaves some responding cells in the baseline tail. ``graded_expectation`` is
    the share expected if every cell responded equally: negative-binomial counts drawn for the time point's cells
    (their library sizes, the group's pooled mean at t, the baseline dispersions) against redrawn baseline cells.

    **Non-responders.** ``nonresponder_share`` = (share - graded) / (1 - graded), clipped to [0, 1], with interval
    ``nonresponder_lo`` / ``nonresponder_hi``: the share of cells left at baseline beyond what an equal response
    explains. If a share of cells never responds and the rest do, it estimates that share (slightly
    conservatively: responders sit above the pooled mean, so they leave fewer cells in the tail than the graded
    model). ``call``:

    - "all respond": the upper bound of ``nonresponder_share`` is <= ``max_nonresponder``.
    - "non-responders": its lower bound is > 0 and its upper bound is above ``max_nonresponder``.
    - "indeterminate": neither, or ``graded_expectation`` > 0.5 (the response is too small relative to count noise
      to tell responding cells from baseline ones; typically the first minutes).
    - "no response": the module has not moved above the group's baseline.

    Also per row: ``n_cells`` and ``module_log2fc`` (the score-weighted mean log2 fold of the module over the group's
    baseline: the pooled change the NB-GAM sees, i.e. share responding x response size).

    Limits: module-level only (per-gene responder shares are dominated by count noise at single-cell depth). A share
    of non-responders and a graded response that is weaker in some cells look alike when the weak cells stay within
    the baseline's spread. A two-class mixture model is deliberately not fitted: it forces non-responders to exactly
    the baseline and overcounts them (30-45% vs ~10-20% in K562 + TNF-alpha).
    """
    validate_adata(adata, time_key=time_key, layer=layer)
    rng = np.random.default_rng(random_state)
    X = adata.layers[layer] if layer is not None else adata.X
    X = X.tocsr() if sp.issparse(X) else np.asarray(X, dtype=float)
    time = pd.to_numeric(adata.obs[time_key].astype(str), errors="coerce").to_numpy(float)
    if np.isnan(time).any():
        raise ChronoPackError(f"adata.obs[{time_key!r}] contains non-numeric values.")
    times = np.unique(time)
    if times.size < 2:
        raise ChronoPackError("response_heterogeneity needs at least two time points.")
    b = 0 if baseline is None else int(np.flatnonzero(np.isclose(times, float(baseline)))[0])
    lib = (adata.obs[library_size_key].to_numpy(float) if library_size_key is not None and library_size_key in adata.obs
           else np.asarray(X.sum(axis=1)).ravel().astype(float))
    s = lib / np.median(lib)
    var_names = adata.var_names.astype(str)

    if genes is None:
        pvalues = None
        for col in (f"{dynamics_key}_pvalue_1", f"{dynamics_key}_pvalue"):
            if col in adata.var:
                pvalues = adata.var[col].fillna(1.0).to_numpy(float)
                break
        genes = _select_module(X, var_names, time, times, b, s, n_genes, min_mean_counts, min_log2fc, pvalues,
                               p_threshold)
        if len(genes) < 3:
            raise ChronoPackError(f"Only {len(genes)} genes pass the module criteria; pass genes= or relax "
                                  "min_mean_counts / min_log2fc.")
    else:
        genes = [str(g) for g in ([genes] if isinstance(genes, str) else genes)]
        missing = sorted(set(genes) - set(var_names))
        if missing:
            raise ChronoPackError(f"genes not in adata.var_names: {missing[:5]}")
    M = _dense(X[:, var_names.get_indexer(genes)])

    # pooled profiles, peak and score weights
    base = time == times[b]
    par0 = np.array([_fit_nb(M[base, j], s[base]) for j in range(len(genes))])
    r0 = par0[:, 1]
    floor = 0.5 / s[base].sum()
    mt = np.vstack([M[time == t].sum(axis=0) / s[time == t].sum() for t in times])
    m0 = np.maximum(mt[b], floor)
    others = np.flatnonzero(np.arange(len(times)) != b)
    k_pk = others[np.argmax(np.mean(np.log(np.maximum(mt[others], floor) / m0), axis=1))]
    m_pk = np.maximum(mt[k_pk], floor)
    w = np.log(m_pk / m0)
    if not np.any(w > 0):
        raise ChronoPackError("The module does not rise above the baseline at any time point.")
    raw = M @ w - s * (m_pk - m0).sum()
    lo_ref, hi_ref = np.median(raw[base]), np.median(raw[time == times[k_pk]])
    score = (raw - lo_ref) / (hi_ref - lo_ref) if hi_ref > lo_ref else raw - lo_ref

    def sim_scores(sf, means):
        """Scores of NB counts drawn for cells with size factors sf and profile means (genes or cells x genes)."""
        mu = sf[:, None] * means
        Y = rng.negative_binomial(r0[None, :], r0[None, :] / (r0[None, :] + mu))
        return Y @ w - sf * (m_pk - m0).sum()

    group_labels = [("all", np.ones(adata.n_obs, bool))]
    if groupby is not None:
        lab = adata.obs[groupby].astype(str).to_numpy()
        group_labels += [(g, lab == g) for g in pd.unique(lab)]

    rows = []
    for gname, gmask in group_labels:
        gb = gmask & base
        if gb.sum() < 10:
            continue
        g_m0 = np.maximum(M[gb].sum(axis=0) / s[gb].sum(), floor)
        g_mt = {t: M[gmask & (time == t)].sum(axis=0) / max(s[gmask & (time == t)].sum(), 1e-12) for t in times}
        shift = {t: float(w @ (g_mt[t] - g_m0)) for t in times}
        sc0 = raw[gb]
        sb = s[gb]
        for i in others:
            t = times[i]
            cm = gmask & (time == t)
            if cm.sum() < 10:
                continue
            sc_t, st = raw[cm], s[cm]
            share = _tail_ratio(sc_t, sc0, quantile)
            boot = [_tail_ratio(rng.choice(sc_t, sc_t.size), rng.choice(sc0, sc0.size), quantile)
                    for _ in range(n_boot)]
            lo, hi = np.quantile(boot, [0.025, 0.975])
            g_now = np.maximum(g_mt[t], floor)[None, :]
            graded = float(np.mean([_tail_ratio(sim_scores(st, g_now), sim_scores(sb, g_m0[None, :]), quantile)
                                    for _ in range(n_sim)]))
            nr = [float(np.clip((x - graded) / (1 - graded), 0, 1)) if graded < 1 else np.nan for x in (share, lo, hi)]
            fold = float(np.sum(w * np.log2(np.maximum(g_mt[t], floor) / g_m0)) / np.sum(np.abs(w)))
            if shift[t] <= 0:
                call = "no response"
            elif graded > 0.5:
                call = "indeterminate"
            elif nr[2] <= max_nonresponder:
                call = "all respond"
            elif nr[1] > 0:
                call = "non-responders"
            else:
                call = "indeterminate"
            rows.append({
                "group": gname, "time": t, "n_cells": int(cm.sum()), "module_log2fc": fold,
                "share_at_baseline": share, "share_lo": float(lo), "share_hi": float(hi),
                "graded_expectation": graded, "nonresponder_share": nr[0], "nonresponder_lo": nr[1],
                "nonresponder_hi": nr[2], "call": call,
            })
    by_time = pd.DataFrame(rows)
    module = pd.DataFrame({"baseline_mean": m0, "peak_mean": m_pk, "weight": w, "dispersion_size": r0},
                          index=pd.Index(genes, name="gene"))
    adata.obs[score_key] = score
    params = {"time_key": time_key, "layer": layer, "groupby": groupby, "library_size_key": library_size_key,
              "baseline": float(times[b]), "quantile": quantile, "max_nonresponder": max_nonresponder,
              "n_boot": n_boot, "n_sim": n_sim, "random_state": random_state, "genes": list(genes)}
    record_provenance(
        adata,
        key_added,
        analysis="response_heterogeneity",
        parameters=params,
        result_keys={"obs": score_key, "peak_time": float(times[k_pk]),
                     "calls": [f"{r.group}|{r.time:g}|{r.call}" for r in by_time.itertuples()]},
    )
    return ResponseHeterogeneityResult(by_time=by_time, module=module,
                                       score=pd.Series(score, index=adata.obs_names, name=score_key),
                                       peak_time=float(times[k_pk]), parameters=params)


__all__ = ["ResponseHeterogeneityResult", "response_heterogeneity"]
