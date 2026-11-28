"""Estimate the time alignment between contexts from the rise of reporter genes.

Each context's response is mapped onto a reference clock as ``aligned = scale * own minutes + offset`` (the maps
:func:`chronopack.tl.align_gene_trajectories` takes). The estimate uses the first rise of one or more reporter genes,
from the per-sampling-time (tag) means of the cells, not the NB-GAM: the two contexts' per-tag curves are registered
onto each other (no curve shape assumed), with logistic rise fits as descriptors and a secondary estimate.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import optimize, stats
from scipy.interpolate import PchipInterpolator

from chronopack._core._validation import ChronoPackError

FIT_SCALES = ("linear", "log")
_RISE_10_90 = 2.0 * np.log(9.0)  # 10-90% rise time of a logistic, in units of its width


@dataclass
class TimeAlignmentResult:
    """Output of :func:`estimate_time_alignment`.

    ``alignment``: one row per context (scale, offset, intervals, the rise-width estimate, flags). ``fits``: one row
    per context x reporter (the logistic rise fit and its flags). ``pairs``: one row per non-reference context x
    reporter (that reporter's own registration and rise-width estimates; onset / midpoint / peak shifts after the
    context's alignment). ``tags``: the per-tag data, with the fit window. ``scale_map`` / ``offset_map``: the
    estimates, ready for :func:`align_gene_trajectories` (None where a context has no estimate).
    """

    alignment: pd.DataFrame
    fits: pd.DataFrame
    pairs: pd.DataFrame
    tags: pd.DataFrame
    scale_map: dict
    offset_map: dict
    parameters: dict = field(default_factory=dict)


# ---------------------------------------------------------------------------------------------------- data loading


def _sce_cells(path, genes, lineage):
    """Counts for ``genes``, pseudotime and size factors (exp of the tradeSeq offset) from a tradeSeq SCE."""
    try:
        from rpy2.robjects import StrVector, r
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError("Reading tradeSeq SCEs requires the optional dependency 'rpy2'.") from exc
    reader = r("""
      function(path, genes, lineage) {
        sce <- if (is.character(path)) readRDS(path) else path
        dm <- SummarizedExperiment::colData(sce)$tradeSeq$dm
        if (is.null(dm)) stop("Not a tradeSeq SCE (colData(sce)$tradeSeq$dm is missing).")
        g <- intersect(genes, rownames(sce))
        cnt <- as.matrix(SummarizedExperiment::assay(sce, "counts")[g, , drop = FALSE])
        off_col <- grep("^offset", colnames(dm), value = TRUE)
        off <- if (length(off_col)) as.numeric(dm[[off_col[1]]]) else rep(0, ncol(sce))
        tcol <- paste0("t", lineage); lcol <- paste0("l", lineage)
        if (!tcol %in% colnames(dm)) stop("Lineage ", lineage, " not found in the SCE.")
        keep <- if (lcol %in% colnames(dm)) as.numeric(dm[[lcol]]) == 1 else rep(TRUE, ncol(sce))
        list(genes = g, counts = as.numeric(cnt), n_cells = ncol(cnt), t = as.numeric(dm[[tcol]]),
             offset = off, keep = keep)
      }""")
    res = reader(str(path), StrVector(list(genes)), int(lineage))
    found = [str(g) for g in res.rx2("genes")]
    n_cells = int(res.rx2("n_cells")[0])
    counts = np.asarray(res.rx2("counts"), dtype=float).reshape((len(found), n_cells), order="F").T
    t = np.asarray(res.rx2("t"), dtype=float)
    sf = np.exp(np.asarray(res.rx2("offset"), dtype=float))
    keep = np.asarray(res.rx2("keep"), dtype=bool)
    return found, counts[keep], t[keep], sf[keep]


def _adata_cells(adata, genes, time_key, layer, library_size_key):
    if time_key not in adata.obs:
        raise ChronoPackError(f"adata.obs has no {time_key!r} column.")
    X = adata.layers[layer] if layer is not None else adata.X
    var_names = pd.Index(adata.var_names.astype(str))
    found = [g for g in genes if g in var_names]
    sub = X[:, var_names.get_indexer(found)]
    counts = sub.toarray() if sp.issparse(sub) else np.asarray(sub, dtype=float)
    t = pd.to_numeric(adata.obs[time_key].astype(str), errors="coerce").to_numpy(float)
    if np.isnan(t).any():
        raise ChronoPackError(f"adata.obs[{time_key!r}] contains non-numeric values.")
    if library_size_key is not None and library_size_key in adata.obs:
        lib = adata.obs[library_size_key].to_numpy(float)
    else:
        lib = np.asarray(X.sum(axis=1)).ravel().astype(float)
    return found, np.asarray(counts, dtype=float), t, lib


def _load(data, ctx, genes, time_map, time_key, layer, library_size_key, lineage):
    """Per-cell own minutes (0 = first sampling time) and counts / size factor (mean size factor 1)."""
    obj = data[ctx]
    if isinstance(obj, (str, Path)):
        if time_map is None or ctx not in time_map:
            raise ChronoPackError(f"time_map[{ctx!r}] (minutes spanned by the SCE's pseudotime) is required for SCEs.")
        found, counts, t, sf = _sce_cells(obj, genes, lineage)
        minutes = (t - t.min()) / (t.max() - t.min()) * float(time_map[ctx])
    else:
        found, counts, t, sf = _adata_cells(obj, genes, time_key, layer, library_size_key)
        span = t.max() - t.min()
        minutes = (t - t.min()) / span * float(time_map[ctx]) if time_map and ctx in time_map else t - t.min()
    if np.any(~np.isfinite(sf)) or np.any(sf <= 0):
        raise ChronoPackError(f"{ctx}: size factors must be positive.")
    sf = sf / sf.mean()
    return found, counts / sf[:, None], np.round(minutes, 6)


# ---------------------------------------------------------------------------------------------------- per-tag data


def _tag_means(X, minutes, tags):
    mean = np.empty((len(tags), X.shape[1]))
    se = np.empty_like(mean)
    n = np.empty(len(tags), dtype=int)
    for k, tg in enumerate(tags):
        xk = X[minutes == tg]
        n[k] = xk.shape[0]
        mean[k] = xk.mean(axis=0)
        se[k] = xk.std(axis=0, ddof=1) / np.sqrt(n[k]) if n[k] > 1 else np.nan
    return mean, se, n


def _boot_tag_means(X, minutes, tags, n_boot, rng):
    """Cell bootstrap within each tag: means and SEs, shape (n_boot, tags, genes)."""
    mean = np.empty((n_boot, len(tags), X.shape[1]))
    se = np.empty_like(mean)
    for k, tg in enumerate(tags):
        xk = X[minutes == tg]
        n = xk.shape[0]
        W = rng.multinomial(n, np.full(n, 1.0 / n), size=n_boot).astype(float)
        m1 = W @ xk / n
        m2 = W @ (xk * xk) / n
        mean[:, k] = m1
        se[:, k] = np.sqrt(np.clip(m2 - m1 * m1, 0, None) * n / max(n - 1, 1) / n)
    return mean, se


def _to_fit_scale(mean, se, scale, se_floor):
    """Values and SEs on the fit scale. SE floor: ``se_floor`` on the log scale; ``se_floor`` x mean on the linear
    scale (the same relative floor)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        if scale == "log":
            y = np.where(mean > 0, np.log(np.where(mean > 0, mean, 1.0)), np.nan)
            s = np.maximum(se / mean, se_floor)
        else:
            y = mean.astype(float)
            peak = np.nanmax(np.abs(mean), axis=-2, keepdims=True)
            s = np.maximum(se, np.maximum(se_floor * np.abs(mean), 1e-3 * peak))
    return y, s


# ---------------------------------------------------------------------------------------------------- rise fits


def _sigmoid(x):
    return 0.5 * (1.0 + np.tanh(0.5 * x))


def _wls_chi2(H, y, wt):
    """Weighted least squares of y on b + a H (a >= 0; a constant fit otherwise), vectorized over H's leading axes."""
    Sw, Sy, Syy = wt.sum(), (wt * y).sum(), (wt * y * y).sum()
    Sx, Sxx, Sxy = (wt * H).sum(-1), (wt * H * H).sum(-1), (wt * H * y).sum(-1)
    den = Sw * Sxx - Sx * Sx
    with np.errstate(divide="ignore", invalid="ignore"):
        a = (Sw * Sxy - Sx * Sy) / den
        b = (Sy - a * Sx) / Sw
        chi = Syy - 2 * b * Sy - 2 * a * Sxy + b * b * Sw + 2 * a * b * Sx + a * a * Sxx
    chi_const = Syy - Sy * Sy / Sw
    ok = np.isfinite(a) & (a > 0) & (den > 1e-12 * Sw * np.maximum(Sxx, 1e-300))
    return np.where(ok, np.maximum(chi, 0.0), chi_const), np.where(ok, a, 0.0), np.where(ok, b, Sy / Sw)


def _fit_rise(t, y, s, m_bounds, w_bounds, start=None, grid=(41, 25)):
    """Logistic rise y = b + a / (1 + exp(-(t - m) / w)), a >= 0, by variable projection: (b, a) solved exactly for
    each (m, w); (m, log w) from a grid, then refined. Returns b, a, m, w, chi2."""
    ok = np.isfinite(y) & np.isfinite(s) & (s > 0)
    t, y, s = t[ok], y[ok], s[ok]
    if t.size < 5:
        return dict(b=np.nan, a=np.nan, m=np.nan, w=np.nan, chi2=np.nan, chi2_const=np.nan, n=int(t.size))
    wt = 1.0 / s**2
    lw_bounds = (np.log(w_bounds[0]), np.log(w_bounds[1]))

    def prof(theta):
        return float(_wls_chi2(_sigmoid((t - theta[0]) / np.exp(theta[1])), y, wt)[0])

    M = np.linspace(m_bounds[0], m_bounds[1], grid[0])
    LW = np.linspace(lw_bounds[0], lw_bounds[1], grid[1])
    H = _sigmoid((t[None, None, :] - M[:, None, None]) / np.exp(LW)[None, :, None])
    chi = _wls_chi2(H, y, wt)[0]
    i, j = np.unravel_index(np.argmin(chi), chi.shape)
    starts = [(M[i], LW[j])] + ([(start[0], np.log(start[1]))] if start is not None and np.all(np.isfinite(start))
                                else [])
    best = None
    for x0 in starts:
        x0 = (np.clip(x0[0], *m_bounds), np.clip(x0[1], *lw_bounds))
        o = optimize.minimize(prof, x0, method="L-BFGS-B", bounds=[m_bounds, lw_bounds])
        if best is None or o.fun < best.fun:
            best = o
    m, w = float(best.x[0]), float(np.exp(best.x[1]))
    c, a, b = _wls_chi2(_sigmoid((t - m) / w), y, wt)
    chi_const = float(_wls_chi2(np.zeros_like(t), y, wt)[0])
    return dict(b=float(b), a=float(a), m=m, w=w, chi2=float(c), chi2_const=chi_const, n=int(t.size))


def _window(mean_lin, direction, peak_frac, extra_tags):
    """Fit window: from the first tag through ``extra_tags`` past the first tag whose change from the first tag (in
    the gene's direction, linear scale) reaches ``peak_frac`` of the largest change, and at least 5 tags. Returns that
    tag's index, the window's last index and whether the largest change is at the last tag (rise may not be
    finished)."""
    z = direction * (mean_lin - mean_lin[0])
    D = np.nanmax(z)
    T = len(z)
    hits = np.flatnonzero(z >= peak_frac * D) if np.isfinite(D) and D > 0 else np.array([T - 1])
    top = int(hits[0])
    return top, min(max(top + extra_tags, 4), T - 1), bool(np.nanargmax(z) == T - 1)


def _direction(mean_lin, spec, gene):
    if spec is not None:
        d = spec.get(gene) if isinstance(spec, dict) else spec
        if d is not None:
            d = {"induced": 1, "up": 1, "repressed": -1, "down": -1}.get(d, d)
            if d not in (1, -1):
                raise ChronoPackError(f"direction for {gene!r} must be 'induced'/'repressed' (or 1/-1).")
            return int(d)
    with np.errstate(divide="ignore", invalid="ignore"):
        lg = np.log(np.maximum(mean_lin, 1e-12))
    dev = lg - lg[0]
    return 1 if dev[int(np.nanargmax(np.abs(dev)))] >= 0 else -1


# ---------------------------------------------------------------------------------------------------- registration


def _curve(t, y, s, rng=None):
    """A context's per-tag curve as a template: monotone (PCHIP) interpolation of its values, flat before the first
    and after the last tag, with the tags' SEs (for the template's own uncertainty). With ``rng`` (bootstrap
    replicates), one standard normal draw per interval between tags: the replicate's curve between those tags is
    displaced by that many interpolation SDs (see :func:`_curve_eval`), so the intervals include interpolation
    uncertainty."""
    ok = np.isfinite(t) & np.isfinite(y) & np.isfinite(s)
    t, y, s = t[ok], y[ok], s[ok]
    if t.size < 3:
        return None
    z = rng.standard_normal(t.size - 1) if rng is not None else np.zeros(t.size - 1)
    return t, y, s, PchipInterpolator(t, y, extrapolate=False), z


def _curve_eval(curve, tau, interp_sd):
    """Template value and variance at ``tau``: PCHIP value; variance = the two neighbouring tags' variances combined
    as for linear interpolation + (``interp_sd`` x 4u(1-u) x the change between them)^2 (u: position between the tags,
    so the interpolation term is 0 at the tags and ``interp_sd`` x the change at the midpoint)."""
    t, y, s, f, z = curve
    tc = np.clip(tau, t[0], t[-1])
    k = np.clip(np.searchsorted(t, tc, side="right") - 1, 0, t.size - 2)
    u = (tc - t[k]) / (t[k + 1] - t[k])
    isd = interp_sd * 4 * u * (1 - u) * (y[k + 1] - y[k])
    var = (1 - u) ** 2 * s[k] ** 2 + u**2 * s[k + 1] ** 2 + isd**2
    return f(tc) + z[k] * np.abs(isd), var


def _irls(H, VT, y, s, n_iter=3):
    """Baseline b and amplitude a >= 0 of y ~ b + a H with variance a^2 VT + s^2 (template and data uncertainty), by
    iteratively reweighted least squares; vectorized over H's leading axes. Returns b, a and the variances."""
    w = np.broadcast_to(1.0 / s**2, H.shape)
    for _ in range(n_iter + 1):
        Sw, Sx, Sy = w.sum(-1), (w * H).sum(-1), (w * y).sum(-1)
        Sxx, Sxy = (w * H * H).sum(-1), (w * H * y).sum(-1)
        den = Sw * Sxx - Sx * Sx
        with np.errstate(divide="ignore", invalid="ignore"):
            a = np.where(den > 1e-12 * Sw * np.maximum(Sxx, 1e-300), (Sw * Sxy - Sx * Sy) / den, 0.0)
        a = np.clip(np.nan_to_num(a), 0.0, None)
        b = (Sy - a * Sx) / Sw
        v = a[..., None] ** 2 * VT + s**2
        w = 1.0 / v
    return b, a, v


def _nll(H, VT, y, s):
    """Gaussian -log-likelihood of y given the template values H (variances VT) after :func:`_irls`, and the sum of
    squared standardized residuals."""
    b, a, v = _irls(H, VT, y, s)
    r2 = (y - b[..., None] - a[..., None] * H) ** 2 / v
    return 0.5 * r2.sum(-1) + 0.5 * np.log(v).sum(-1), r2.sum(-1)


def _register(obj, center, s_range, u_half, fixed_s=None, grid=(61, 61)):
    """Best (scale, offset) of ``obj(S, D)``: grid over log scale (or the fixed scale) x offset = center(scale) + u,
    |u| <= u_half, then a local refinement. Returns scale, offset, objective and whether the grid optimum is on the
    grid's edge."""
    S = np.array([float(fixed_s)]) if fixed_s is not None else np.geomspace(s_range[0], s_range[1], grid[0])
    U = np.linspace(-u_half, u_half, grid[1])
    C = center(S)
    X = obj(S[:, None] * np.ones((1, U.size)), C[:, None] + U[None, :])
    if not np.isfinite(X).any():
        return np.nan, np.nan, np.nan, True
    i, j = np.unravel_index(np.nanargmin(X), X.shape)
    edge = bool(j in (0, U.size - 1) or (fixed_s is None and i in (0, S.size - 1)))
    best = (S[i], C[i] + U[j], X[i, j])
    du = U[1] - U[0]
    if fixed_s is None:
        o = optimize.minimize(lambda th: float(obj(np.exp(th[0]), th[1])), [np.log(S[i]), C[i] + U[j]],
                              method="Nelder-Mead",
                              options=dict(xatol=1e-4, fatol=1e-7, maxiter=400,
                                           initial_simplex=[[np.log(S[i]), C[i] + U[j]],
                                                            [np.log(S[i]) + 0.05, C[i] + U[j]],
                                                            [np.log(S[i]), C[i] + U[j] + du]]))
        if o.fun < best[2]:
            best = (float(np.exp(o.x[0])), float(o.x[1]), float(o.fun))
    else:
        lo, hi = C[0] + U[max(j - 1, 0)], C[0] + U[min(j + 1, U.size - 1)]
        o = optimize.minimize_scalar(lambda d: float(obj(np.asarray(fixed_s, float), np.asarray(d, float))),
                                     bounds=(lo, hi), method="bounded")
        if o.fun < best[2]:
            best = (float(fixed_s), float(o.x), float(o.fun))
    return best[0], best[1], best[2], edge


def _crossing(t, y, frac):
    """First time the window curve (0 at its first tag, 1 at its maximum) reaches ``frac`` (linear interpolation)."""
    ok = np.isfinite(y)
    t, y = t[ok], y[ok]
    if t.size < 2 or np.max(y) <= y[0]:
        return np.nan
    z = (y - y[0]) / (np.max(y) - y[0])
    k = int(np.flatnonzero(z >= frac)[0])
    if k == 0:
        return float(t[0])
    return float(t[k - 1] + (frac - z[k - 1]) / (z[k] - z[k - 1]) * (t[k] - t[k - 1]))


def _landmarks(t, y, last, s=None):
    """Onset (10% of the rise), midpoint (50%) and peak of a per-tag curve in the gene's direction. Rise: from the
    first tag to the largest value up to the window end (index ``last``); crossings on the PCHIP curve. Peak: the
    vertex of the parabola through the first local maximum at or after that tag and its two neighbours; NaN if the
    curve is still rising at its last tag or (with SEs ``s``) does not fall significantly (> 2 SE) within the next 3
    tags (a sustained response has no peak)."""
    nan3 = (np.nan, np.nan, np.nan)
    s = np.zeros_like(y) if s is None else s
    ok = np.isfinite(y) & np.isfinite(s)
    if not ok.all() or t.size < 3:
        t, y, s = t[ok], y[ok], s[ok]
        last = min(last, t.size - 1)
    if t.size < 3:
        return nan3
    k = int(np.argmax(y[: last + 1]))
    if y[k] <= y[0] or k == 0:
        return nan3
    f = PchipInterpolator(t[: k + 1], y[: k + 1])
    g = np.linspace(t[0], t[k], 2001)
    z = (f(g) - y[0]) / (y[k] - y[0])
    onset, mid = float(g[np.argmax(z >= 0.1)]), float(g[np.argmax(z >= 0.5)])
    while k + 1 < t.size and y[k + 1] > y[k]:
        k += 1
    if k + 1 >= t.size:
        return onset, mid, np.nan
    nxt = slice(k + 1, min(k + 4, t.size))
    j = k + 1 + int(np.argmin(y[nxt]))
    if s.any() and not (y[k] - y[j] > 2 * np.hypot(s[k], s[j])):
        return onset, mid, np.nan
    c2, c1, _ = np.polyfit(t[k - 1:k + 2], y[k - 1:k + 2], 2)
    peak = float(np.clip(-c1 / (2 * c2), t[k - 1], t[k + 1])) if c2 < 0 else float(t[k])
    return onset, mid, peak


class _Registrar:
    """Registration of one context onto the reference for a set of reporters: point estimate and bootstrap.

    Objective (summed over reporters): the context's window data fitted with the reference's per-tag curve read at
    scale * t + offset, plus (``symmetric``) the reference's window data fitted with the context's per-tag curve read
    at (t - offset) / scale; each with its own baseline and amplitude, weighted by data + template uncertainty."""

    def __init__(self, ref_w, ctx_w, genes, interp_sd=0.125, symmetric=True, seed=0):
        self.genes, self.seed = genes, int(seed)
        self.phi = 1.0  # dispersion (chi2 / dof at the estimate); > 1 adds that excess tag-level noise to replicates
        self.ref_w, self.ctx_w = ref_w, ctx_w
        self.interp_sd, self.symmetric = float(interp_sd), bool(symmetric)
        t50r = np.array([_crossing(ref_w[g]["t"], ref_w[g]["y"], 0.5) for g in genes])
        t50c = np.array([_crossing(ctx_w[g]["t"], ctx_w[g]["y"], 0.5) for g in genes])
        durr = np.array([_crossing(ref_w[g]["t"], ref_w[g]["y"], 0.8) - _crossing(ref_w[g]["t"], ref_w[g]["y"], 0.2)
                         for g in genes])
        durc = np.array([_crossing(ctx_w[g]["t"], ctx_w[g]["y"], 0.8) - _crossing(ctx_w[g]["t"], ctx_w[g]["y"], 0.2)
                         for g in genes])
        ratio = durr / durc
        ratio = ratio[np.isfinite(ratio) & (ratio > 0)]
        s0 = float(np.exp(np.median(np.log(ratio)))) if ratio.size else 1.0
        self.s_range = (s0 / 4, s0 * 4)
        self.t50r, self.t50c = t50r, t50c
        self.u_half = max(float(ref_w[g]["t"][-1] - ref_w[g]["t"][0]) for g in genes) / 2

    def center(self, S):
        ok = np.isfinite(self.t50r) & np.isfinite(self.t50c)
        if not ok.any():
            return np.zeros_like(S)
        return np.array([np.mean(self.t50r[ok] - s * self.t50c[ok]) for s in S])

    def _inputs(self, r=None):
        out = []
        for gi, g in enumerate(self.genes):
            rw, cw = self.ref_w[g], self.ctx_w[g]
            rng = (lambda side: None if r is None else np.random.default_rng([self.seed, r, gi, side]))

            def pick(w, side):
                if r is None:
                    return w["fy"], w["fs"]
                y = w["fby"][r]
                if self.phi > 1:
                    y = y + rng(side + 2).standard_normal(y.size) * np.sqrt(self.phi - 1) * w["fs"]
                return y, w["fbs"][r]

            (ry, rs), (cy, cs) = pick(rw, 0), pick(cw, 1)
            rc, cc = _curve(rw["ft"], ry, rs, rng(0)), _curve(cw["ft"], cy, cs, rng(1))
            if rc is None or cc is None:
                return None
            data = []
            for w, (fy, fs) in ((cw, (cy, cs)), (rw, (ry, rs))):
                n = w["t"].size
                y, s = fy[:n], fs[:n]
                ok = np.isfinite(y) & np.isfinite(s) & (s > 0)
                if ok.sum() < 3:
                    return None
                data.append((w["t"][ok], y[ok], s[ok]))
            out.append((rc, cc, data[0], data[1]))
        return out

    def _terms(self, inputs, S, D):
        """Per direction and reporter: (nll, chi2, n) at (S, D) (broadcast)."""
        S, D = np.asarray(S, float), np.asarray(D, float)
        terms = []
        for rc, cc, (tc, yc, sc), (tr, yr, sr) in inputs:
            H, V = _curve_eval(rc, S[..., None] * tc + D[..., None], self.interp_sd)
            terms.append(_nll(H, V, yc, sc) + (tc.size,))
            if self.symmetric:
                H, V = _curve_eval(cc, (tr - D[..., None]) / S[..., None], self.interp_sd)
                terms.append(_nll(H, V, yr, sr) + (tr.size,))
        return terms

    def objective(self, inputs):
        return lambda S, D: sum(t[0] for t in self._terms(inputs, S, D))

    def fit(self, r=None, fixed_s=None, s_range=None, grid=(61, 61)):
        inputs = self._inputs(r)
        if inputs is None:
            return np.nan, np.nan, np.nan, True
        return _register(self.objective(inputs), self.center, s_range or self.s_range, self.u_half,
                         fixed_s=fixed_s, grid=grid)

    def goodness(self, s, d, fixed=False):
        """Sum of squared standardized residuals at (s, d) and its degrees of freedom."""
        inputs = self._inputs()
        terms = self._terms(inputs, s, d)
        chi = float(sum(t[1] for t in terms))
        dof = int(sum(t[2] for t in terms)) - 2 * len(terms) - (1 if fixed else 2)
        return chi, max(dof, 1)

    def fit_boot(self, n_boot, fixed_s=None):
        """Point estimate, dispersion and bootstrap: scale, offset, at-edge, replicates, chi2, dof."""
        s, d, _, e = self.fit(fixed_s=fixed_s)
        chi, dof = self.goodness(s, d, fixed=fixed_s is not None) if np.isfinite(s) else (np.nan, 1)
        self.phi = max(1.0, chi / dof) if np.isfinite(chi) else 1.0
        return s, d, e, self.boot(n_boot, (s, d), fixed_s=fixed_s), chi, dof

    def boot(self, n_boot, point, fixed_s=None):
        """Bootstrap replicates (scale, offset, at-edge), each re-registered with its own per-tag curves."""
        out = np.full((n_boot, 3), np.nan)
        if not n_boot or not np.isfinite(point[0]):
            return out
        rng_s = (max(point[0] / 3, self.s_range[0] / 2), min(point[0] * 3, self.s_range[1] * 2))
        for r in range(n_boot):
            s, d, _, e = self.fit(r, fixed_s=fixed_s, s_range=rng_s, grid=(31, 41))
            out[r] = s, d, float(e)
        return out

    def landmark_shifts(self, g, s, d, r=None):
        """Onset, midpoint and peak of the context (mapped by s, d) minus the reference's, reference minutes."""
        out = []
        for w in (self.ctx_w[g], self.ref_w[g]):
            y, se = (w["fy"], w["fs"]) if r is None else (w["fby"][r], None)
            out.append(np.array(_landmarks(w["ft"], y, w["t"].size - 1, se)))
        return s * out[0] + d - out[1]


# ---------------------------------------------------------------------------------------------------- main function


def estimate_time_alignment(
    data,
    reporters,
    *,
    reference=None,
    time_map=None,
    time_key="time",
    layer="counts",
    library_size_key="total_counts",
    lineage=1,
    fit_scale="linear",
    direction=None,
    fix_scale=None,
    current=None,
    extra_tags=1,
    peak_frac=0.9,
    se_floor=0.05,
    interp_sd=0.125,
    symmetric=True,
    amp_floor=0.05,
    n_boot=1000,
    level=0.95,
    random_state=0,
):
    """Estimate how each context's response clock maps onto a reference context's, from reporter-gene rises.

    The alignment is ``aligned minutes = scale * own minutes + offset`` (own minutes counted from each context's
    first sampling time), as used by :func:`chronopack.tl.align_gene_trajectories`. Nothing is applied: the result
    reports estimates, intervals and diagnostics, and ``scale_map`` / ``offset_map`` can be passed on if wanted.

    **Data.** ``data`` maps context -> a tradeSeq SCE (path or object, as for ``align_gene_trajectories``; counts,
    pseudotime and offsets are read from it, and ``time_map[context]`` gives the minutes its pseudotime spans) or an
    AnnData (counts in ``layer``, sampling time in ``obs[time_key]``, library size in ``obs[library_size_key]`` if
    present, else the layer's totals; with ``time_map[context]`` the times are rescaled to span that many minutes).
    Cells must come from discrete sampling times (tags). Per tag and reporter: the mean of counts / size factor over
    cells and its SE, on the linear (``fit_scale="linear"``, default) or log scale; SEs floored at ``se_floor``
    (relative). Fit window (the data each context contributes): the first rise, from the first tag through
    ``extra_tags`` past the first tag whose change from the first tag (in the gene's direction, ``direction`` or
    inferred) reaches ``peak_frac`` of the largest change (at least 5 tags).

    **Estimate (registration).** The context's window data are fitted with the reference's per-tag curve (monotone,
    PCHIP, interpolation through all its tags, held at its first value before the first tag) read at
    ``scale * t + offset`` and, with ``symmetric=True`` (default), the reference's window data with the context's
    per-tag curve read at ``(t - offset) / scale``; each with its own baseline and amplitude. The whole curve is the
    template, not only its window: a context's window can reach past the reference's (e.g. one tag into the decline
    after the peak), and holding the reference flat there mis-times the peak (K562b NFKBIA: scale 1.87 with the
    window template, 0.95 with the whole curve). Weights combine the data SEs with the template's own uncertainty: the
    neighbouring tags' SEs, interpolated, plus an interpolation term of ``interp_sd`` x the change between the two
    neighbouring tags at their midpoint (x 4u(1 - u) at position u between them, 0 at the tags). Precise tags thus
    count fully where the template is well known (e.g. a flat baseline before the onset) and less where it is
    interpolated across a steep change. Gaussian likelihood, baseline and amplitude by iteratively reweighted least
    squares, summed over reporters and both directions; the best (scale, offset) on a grid, then refined. No curve
    shape is assumed. Intervals: ``n_boot`` bootstrap resamples of cells within tags in both contexts, each
    re-registered with its own curves, which are also displaced between tags by draws of the interpolation term; if
    the fit's ``dispersion`` (sum of squared standardized residuals / degrees of freedom) exceeds 1, the replicates
    also get that excess as noise per tag. ``fix_scale`` ({context: scale}) estimates the offset only. Reporters used:
    those with a rise in both contexts. Per reporter (``pairs``) the same is done separately; ``heterogeneity_p``
    tests whether reporters agree.

    **Alignment check** (``pairs``, per reporter, at the context's estimate): the context's onset (10% of the rise
    from the first tag to the window's largest value, on the per-tag curve), midpoint (50%) and peak (vertex of the
    parabola through the first local maximum and its neighbours; only where the curve then falls by more than 2 SE
    within 3 tags), mapped onto the reference clock, minus the reference's (``onset_shift``, ``mid_shift``,
    ``peak_shift``, reference minutes, bootstrap intervals). Landmarks read from per-tag curves also carry
    interpolation error that these intervals do not include (up to a few minutes for a peak sampled every 15 min).

    **Rise fits (descriptors and a secondary estimate).** Per context and reporter, ``y = b + a / (1 + exp(-(t - m) /
    w))`` on the same window (SEs also floored at ``amp_floor`` x the rise's amplitude): midpoint ``m``, width ``w``
    (10-90% rise = 4.39 w), bootstrap intervals, in ``fits``.
    ``rise_scale`` = ``w_reference / w_context`` (``rise_offset`` lines up the midpoints), combined over reporters
    resolved in both contexts. It assumes a logistic shape, so it is a secondary check (benchmark: within ~4-12%).

    ``current`` ({context: (scale, offset)}, or the AnnData returned by ``align_gene_trajectories``) adds the offset
    implied at the current scale (registered with the scale fixed, with its interval) and whether the current values
    sit inside the intervals.

    Fit flags (``fits``): ``width_at_lower_bound`` (the rise happens between two tags), ``width_at_upper_bound``,
    ``midpoint_outside_window``, ``no_rise`` (no rise beyond a flat line, p >= 0.001), ``rise_not_finished`` (largest
    change at the last tag), ``few_tags``. Context flags (``alignment``): ``no_reporter_with_rise``,
    ``context_rise_undersampled`` / ``reference_rise_undersampled`` (every reporter's logistic width is at its lower
    bound: the rise falls between about two tags, so the scale is less certain; when the rise spans about one sampling
    interval it is also biased, see the benchmark), ``grid_edge`` (the optimum, or > 10% of bootstrap optima, at the
    search grid's edge), ``poor_fit`` (the context's curve is not a re-timed copy of the reference's within noise,
    p < 0.001, e.g. a reporter whose shape differs between contexts), ``reporters_disagree`` (``heterogeneity_p`` <
    0.05), ``current_scale_outside_ci``, ``current_offset_outside_ci``.

    Benchmark (synthetic contexts with known maps, 100 replicates per scenario; reference 10-min tags, context
    15-min tags read at scale 0.6 unless stated): on the linear scale the median estimate is within 7% of the truth
    with 96-100% interval coverage across rise shapes (logistic, sustained, pulse), context fold changes 5-80 vs 20,
    60-500 cells per tag, NB size 1-3, low expression, repressed reporters and 5-30-min context sampling. A
    K562b-like context (scale 1, rise covered by about two tags): median 0.96, the aligned onset within a median
    0.3 min of the truth and the peak within 2.8 min. THP1-like (scale 0.55, fold 9, 70 cells per tag): 0.89,
    coverage 92%. Biased: a context rise covered by about one sampling interval (flagged undersampled), -16%; a
    reference sampled more coarsely than the context (15- vs 10-min tags), +18%, so use the most finely sampled context
    as the reference; a spread of response times in one context (sd 10 min), -27%, which no per-tag method can tell
    apart from a slower response. Extra noise per tag (10% CV beyond cell sampling): unbiased, coverage 81-84%. The
    log scale is less accurate (median 0.72-1.24), so use ``fit_scale="log"`` only for comparison.

    Assumptions: the contexts share the reporter's programme, related by a change of clock; the reference samples
    its rise well (preferably the most finely sampled context); cells respond roughly in sync
    (a spread of response times in one context widens its mean rise and lowers its apparent speed, which no
    per-tag method can tell apart). Choose fast, direct targets of the pathway with short-lived mRNAs (e.g. NFKBIA for
    NF-kB); other genes may still differ in timing, which is biology the alignment leaves visible.
    """
    if fit_scale not in FIT_SCALES:
        raise ChronoPackError(f"fit_scale must be one of {FIT_SCALES}.")
    reporters = [reporters] if isinstance(reporters, str) else [str(g) for g in reporters]
    if not reporters:
        raise ChronoPackError("Give at least one reporter gene.")
    contexts = list(data)
    reference = contexts[0] if reference is None else reference
    if reference not in data or len(contexts) < 2:
        raise ChronoPackError("data needs the reference and at least one other context.")
    fix_scale = dict(fix_scale or {})
    current = _current_maps(current)
    rng = np.random.default_rng(random_state)
    reg_seed = int(np.random.default_rng(random_state).integers(2**31))  # template perturbations in the bootstrap
    alpha = (1 - level) / 2
    pct = (100 * alpha, 100 * (1 - alpha))

    fits, tag_rows, boots, windows = [], [], {}, {}
    for ctx in contexts:
        found, X, minutes = _load(data, ctx, reporters, time_map, time_key, layer, library_size_key, lineage)
        missing = sorted(set(reporters) - set(found))
        if missing:
            raise ChronoPackError(f"{ctx}: reporter genes not found: {missing}")
        X = X[:, [found.index(g) for g in reporters]]
        tags = np.unique(minutes)
        if tags.size < 5:
            raise ChronoPackError(f"{ctx}: needs at least 5 sampling times (found {tags.size}).")
        if tags.size > 60:
            raise ChronoPackError(f"{ctx}: {tags.size} distinct times; the method needs discrete sampling times "
                                  "(bin continuous pseudotime first).")
        mean, se, n = _tag_means(X, minutes, tags)
        bmean, bse = _boot_tag_means(X, minutes, tags, n_boot, rng) if n_boot else (None, None)
        step = float(np.median(np.diff(tags)))
        for gi, g in enumerate(reporters):
            d = _direction(mean[:, gi], direction, g)
            top, last, unfinished = _window(mean[:, gi], d, peak_frac, extra_tags)
            win = slice(0, last + 1)
            t = tags[win]
            fy, fs = _to_fit_scale(mean[:, gi:gi + 1], se[:, gi:gi + 1], fit_scale, se_floor)
            fy, fs = d * fy[:, 0], fs[:, 0]
            y = fy[win]
            floor = amp_floor * (np.nanmax(y) - np.nanmin(y)) if np.isfinite(y).any() else 0.0
            s = np.maximum(fs[win], floor)  # rise fits only; the registration weighs by data + template uncertainty
            span = float(t[-1] - t[0])
            m_bounds = (float(t[0] - 2 * step), float(t[-1] + 2 * step))
            w_bounds = (0.1 * step, max(span, 2 * step))
            f = _fit_rise(t, y, s, m_bounds, w_bounds)
            B = np.full((n_boot, 2), np.nan)
            fby = fbs = None
            if n_boot:
                fby, fbs = _to_fit_scale(bmean[:, :, gi:gi + 1], bse[:, :, gi:gi + 1], fit_scale, se_floor)
                fby, fbs = d * fby[:, :, 0], fbs[:, :, 0]
                by, bs = fby[:, win], np.maximum(fbs[:, win], floor)
                if np.isfinite(f["w"]):
                    for r in range(n_boot):
                        fb = _fit_rise(t, by[r], bs[r], m_bounds, w_bounds, start=(f["m"], f["w"]), grid=(21, 13))
                        B[r] = fb["m"], fb["w"]
            boots[(ctx, g)] = B
            lin = d * mean[win, gi]
            fold = float(np.max(lin) / lin[0]) if lin[0] > 0 else np.inf
            windows[(ctx, g)] = dict(t=t, y=y, ft=tags, fy=fy, fs=fs, fby=fby, fbs=fbs, top=top, direction=d, step=step,
                                     fold=fold if d > 0 else float(lin[0] / np.max(lin)) if np.max(lin) > 0 else np.inf)
            dchi = f["chi2_const"] - f["chi2"]
            rise_p = float(stats.chi2.sf(dchi, 3)) if np.isfinite(dchi) else np.nan
            m_lo, m_hi = (_pct(B[:, 0], pct[0]), _pct(B[:, 0], pct[1])) if n_boot else (np.nan, np.nan)
            w_lo, w_hi = (_pct(B[:, 1], pct[0]), _pct(B[:, 1], pct[1])) if n_boot else (np.nan, np.nan)
            flags = []
            if np.isfinite(w_lo) and w_lo <= 1.5 * w_bounds[0]:
                flags.append("width_at_lower_bound")
            if np.isfinite(w_hi) and w_hi >= w_bounds[1] / 1.5:
                flags.append("width_at_upper_bound")
            if not (t[0] <= f["m"] <= t[-1]):
                flags.append("midpoint_outside_window")
            if not (rise_p < 1e-3):
                flags.append("no_rise")
            if unfinished:
                flags.append("rise_not_finished")
            if t.size < 5:
                flags.append("few_tags")
            resolved = not ({"width_at_lower_bound", "width_at_upper_bound", "midpoint_outside_window", "no_rise",
                             "few_tags"} & set(flags)) and n_boot > 0
            fits.append({
                "context": ctx, "gene": g, "direction": "induced" if d > 0 else "repressed", "fit_scale": fit_scale,
                "b": d * f["b"], "a": d * f["a"], "m": f["m"], "w": f["w"], "m_lo": m_lo, "m_hi": m_hi,
                "w_lo": w_lo, "w_hi": w_hi, "rise_10_90": _RISE_10_90 * f["w"],
                "onset_10": f["m"] - 0.5 * _RISE_10_90 * f["w"], "fold_change": windows[(ctx, g)]["fold"],
                "top_tag_min": float(tags[top]), "window_end_min": float(t[-1]), "n_tags": int(t.size),
                "chi2": f["chi2"], "rise_p": rise_p, "fit_p": float(stats.chi2.sf(f["chi2"], max(f["n"] - 4, 1))),
                "resolved": bool(resolved), "flags": ";".join(flags),
            })
            yy, ss = _to_fit_scale(mean[:, gi:gi + 1], se[:, gi:gi + 1], fit_scale, se_floor)
            for k, tg in enumerate(tags):
                tag_rows.append({"context": ctx, "gene": g, "minutes": float(tg), "n_cells": int(n[k]),
                                 "mean": float(mean[k, gi]), "se": float(se[k, gi]), "y": float(yy[k, 0]),
                                 "y_se": float(ss[k, 0]), "in_window": bool(k <= last)})
    fits = pd.DataFrame(fits)
    F = fits.set_index(["context", "gene"])

    pair_rows, align_rows = [], []
    for ctx in contexts:
        if ctx == reference:
            continue
        ok_g = [g for g in reporters if "no_rise" not in F.loc[(reference, g)]["flags"]
                and "no_rise" not in F.loc[(ctx, g)]["flags"] and F.loc[(ctx, g)].n_tags >= 5
                and F.loc[(reference, g)].n_tags >= 5]
        fixed = fix_scale.get(ctx)
        cur = current.get(ctx)
        # per reporter
        per_boot, per_fit, ctx_pairs = {}, {}, []
        for g in reporters:
            row = {"context": ctx, "gene": g, "used": g in ok_g}
            if g in ok_g:
                reg = _Registrar({g: windows[(reference, g)]}, {g: windows[(ctx, g)]}, [g], interp_sd, symmetric,
                                 seed=reg_seed)
                s_g, d_g, e_g, Bg, chi_g, dof_g = reg.fit_boot(n_boot, fixed_s=fixed)
                per_boot[g] = Bg
                per_fit[g] = (reg, s_g, d_g, e_g, Bg, chi_g, dof_g)
                row.update(scale=s_g, scale_lo=_pct(Bg[:, 0], pct[0]) if fixed is None else np.nan,
                           scale_hi=_pct(Bg[:, 0], pct[1]) if fixed is None else np.nan, offset=d_g,
                           offset_lo=_pct(Bg[:, 1], pct[0]), offset_hi=_pct(Bg[:, 1], pct[1]), grid_edge=e_g)
            r_, c_ = F.loc[(reference, g)], F.loc[(ctx, g)]
            sb = boots[(reference, g)][:, 1] / boots[(ctx, g)][:, 1]
            row.update(rise_usable=bool(r_.resolved and c_.resolved), rise_scale=r_.w / c_.w,
                       rise_scale_lo=_pct(sb, pct[0]), rise_scale_hi=_pct(sb, pct[1]),
                       rise_offset=r_.m - (r_.w / c_.w) * c_.m)
            ctx_pairs.append(row)
        arow, reg, B = _align_row(ctx, reference, ok_g, windows, F, boots, per_boot, fixed, cur, n_boot, pct,
                                  interp_sd, symmetric, per_fit, seed=reg_seed)
        align_rows.append(arow)
        # onset / midpoint / peak of the context, mapped by the context's estimate, minus the reference's
        for row in ctx_pairs:
            if reg is None or not row["used"] or not np.isfinite(arow["scale"]):
                continue
            g = row["gene"]
            sh = reg.landmark_shifts(g, arow["scale"], arow["offset"])
            Bs = np.array([reg.landmark_shifts(g, B[r, 0], B[r, 1], r) for r in range(n_boot)
                           if np.isfinite(B[r, 0])]) if n_boot else np.empty((0, 3))
            for k, name in enumerate(("onset", "mid", "peak")):
                defined = np.isfinite(sh[k]) and len(Bs)
                row[f"{name}_shift"] = float(sh[k])
                row[f"{name}_shift_lo"] = _pct(Bs[:, k], pct[0]) if defined else np.nan
                row[f"{name}_shift_hi"] = _pct(Bs[:, k], pct[1]) if defined else np.nan
        pair_rows.extend(ctx_pairs)
    pairs = pd.DataFrame(pair_rows)
    alignment = pd.DataFrame([{"context": reference, "reference": reference, "method": "reference", "scale": 1.0,
                               "offset": 0.0, "flags": ""}] + align_rows)
    scale_map = {r.context: (None if not np.isfinite(r.scale) else float(r.scale)) for r in alignment.itertuples()}
    offset_map = {r.context: (None if not np.isfinite(r.offset) else float(r.offset)) for r in alignment.itertuples()}
    params = dict(reporters=reporters, reference=reference, fit_scale=fit_scale, extra_tags=extra_tags,
                  peak_frac=peak_frac, se_floor=se_floor, interp_sd=interp_sd, symmetric=symmetric,
                  amp_floor=amp_floor, n_boot=n_boot, level=level,
                  fix_scale=fix_scale, time_map=dict(time_map or {}), random_state=random_state)
    return TimeAlignmentResult(alignment=alignment, fits=fits, pairs=pairs, tags=pd.DataFrame(tag_rows),
                               scale_map=scale_map, offset_map=offset_map, parameters=params)


# ---------------------------------------------------------------------------------------------------- helpers


def _pct(x, q):
    x = np.asarray(x, dtype=float)
    return float(np.nanpercentile(x, q)) if np.isfinite(x).any() else np.nan


def _current_maps(current):
    if current is None:
        return {}
    if hasattr(current, "uns"):
        from chronopack._core._keys import CHRONOPACK_UNS_KEY

        meta = current.uns.get(CHRONOPACK_UNS_KEY, {}).get("alignment", {})
        sm = meta.get("scale_map_use", meta.get("scale_map", {}))
        om = meta.get("offset_map", {})
        return {str(k): (float(sm.get(k, 1.0)), float(v)) for k, v in om.items()}
    return {str(k): (float(v[0]), float(v[1])) for k, v in dict(current).items()}


def _align_row(ctx, reference, use, windows, F, boots, per_boot, fixed, cur, n_boot, pct, interp_sd, symmetric,
               per_fit=None, seed=0):
    row = {"context": ctx, "reference": reference, "method": "registration" if fixed is None else "fixed scale",
           "reporters_used": ",".join(use), "n_reporters_used": len(use)}
    flags = []
    reg, B = None, np.full((n_boot, 3), np.nan)
    nan4 = dict(scale=np.nan, scale_lo=np.nan, scale_hi=np.nan, offset=np.nan, offset_lo=np.nan, offset_hi=np.nan)
    if not use:
        row.update(nan4, fit_p=np.nan, heterogeneity_p=np.nan)
        flags.append("no_reporter_with_rise")
    else:
        if len(use) == 1 and per_fit and use[0] in per_fit:  # one reporter: the same fit as its own registration
            reg, s_hat, d_hat, edge, B, chi, dof = per_fit[use[0]]
        else:
            reg = _Registrar({g: windows[(reference, g)] for g in use}, {g: windows[(ctx, g)] for g in use}, use,
                             interp_sd, symmetric, seed=seed)
            s_hat, d_hat, edge, B, chi, dof = reg.fit_boot(n_boot, fixed_s=fixed)
        row.update(scale=s_hat, scale_lo=_pct(B[:, 0], pct[0]) if fixed is None else np.nan,
                   scale_hi=_pct(B[:, 0], pct[1]) if fixed is None else np.nan, offset=d_hat,
                   offset_lo=_pct(B[:, 1], pct[0]), offset_hi=_pct(B[:, 1], pct[1]),
                   fit_p=float(stats.chi2.sf(chi, dof)) if np.isfinite(chi) else np.nan,
                   dispersion=float(chi / dof) if np.isfinite(chi) else np.nan,
                   boot_at_grid_edge=float(np.nanmean(B[:, 2])) if n_boot else np.nan)
        if edge or (n_boot and np.nanmean(B[:, 2]) > 0.1):
            flags.append("grid_edge")
        if np.isfinite(row["fit_p"]) and row["fit_p"] < 1e-3:
            flags.append("poor_fit")
        # heterogeneity across reporters (scale; offset if the scale is fixed)
        col = 0 if fixed is None else 1
        vals, var = [], []
        for g in use:
            Bg = per_boot.get(g)
            if Bg is None:
                continue
            x = np.log(Bg[:, 0]) if col == 0 else Bg[:, 1]
            v = np.nanvar(x) if np.isfinite(x).sum() > 1 else np.nan
            if np.isfinite(v) and v > 0:
                vals.append(np.nanmedian(x))
                var.append(v)
        if len(vals) > 1:
            vals, var = np.array(vals), np.array(var)
            mu = np.sum(vals / var) / np.sum(1 / var)
            het_p = float(stats.chi2.sf(np.sum((vals - mu) ** 2 / var), len(vals) - 1))
        else:
            het_p = np.nan
        row["heterogeneity_p"] = het_p
        if np.isfinite(het_p) and het_p < 0.05:
            flags.append("reporters_disagree")
        if all("width_at_lower_bound" in F.loc[(ctx, g)]["flags"] for g in use):
            flags.append("context_rise_undersampled")
        if all("width_at_lower_bound" in F.loc[(reference, g)]["flags"] for g in use):
            flags.append("reference_rise_undersampled")
        if cur is not None:
            cs, co = cur
            oc = reg.fit(fixed_s=cs)
            Bc = reg.boot(n_boot, (cs, oc[1]), fixed_s=cs)
            row.update(current_scale=cs, current_offset=co, offset_at_current_scale=oc[1],
                       offset_at_current_scale_lo=_pct(Bc[:, 1], pct[0]),
                       offset_at_current_scale_hi=_pct(Bc[:, 1], pct[1]))
            if fixed is None and np.isfinite(row["scale_lo"]) and not (row["scale_lo"] <= cs <= row["scale_hi"]):
                flags.append("current_scale_outside_ci")
            lo, hi = row["offset_at_current_scale_lo"], row["offset_at_current_scale_hi"]
            if np.isfinite(lo) and not (lo <= co <= hi):
                flags.append("current_offset_outside_ci")
    # secondary: logistic rise widths, combined over reporters resolved in both contexts
    ru = [g for g in use if F.loc[(reference, g)].resolved and F.loc[(ctx, g)].resolved]
    row.update(_rise_estimate(ctx, reference, ru, F, boots, pct))
    row["flags"] = ";".join(flags)
    return row, reg, B


def _rise_estimate(ctx, reference, use, F, boots, pct):
    out = dict(rise_reporters=",".join(use), rise_scale=np.nan, rise_scale_lo=np.nan, rise_scale_hi=np.nan,
               rise_offset=np.nan)
    if not use:
        return out
    L = np.array([np.log(boots[(reference, g)][:, 1] / boots[(ctx, g)][:, 1]) for g in use])
    l_hat = np.array([np.log(F.loc[(reference, g)].w / F.loc[(ctx, g)].w) for g in use])
    var = np.nanvar(L, axis=1)
    wts = 1.0 / np.where(var > 0, var, np.nan)
    if not np.isfinite(wts).any():
        wts = np.ones(len(use))
    wts = np.where(np.isfinite(wts), wts, 0) / np.nansum(wts)
    s_hat = float(np.exp(np.sum(wts * l_hat)))
    s_b = np.exp(np.nansum(wts[:, None] * L, axis=0) / np.where(np.isfinite(L), wts[:, None], 0).sum(axis=0))
    mr = np.array([F.loc[(reference, g)].m for g in use])
    mc = np.array([F.loc[(ctx, g)].m for g in use])
    out.update(rise_scale=s_hat, rise_scale_lo=_pct(s_b, pct[0]), rise_scale_hi=_pct(s_b, pct[1]),
               rise_offset=float(np.sum(wts * (mr - s_hat * mc))))
    return out
