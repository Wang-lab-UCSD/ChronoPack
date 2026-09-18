"""tradeSeq input-preparation and time-tag-bias helpers."""

from __future__ import annotations

import os
import pickle

import numpy as np
import pandas as pd
import scipy.sparse as sp


def load_negative_control_logcpm(path):
    """Load a time-tag aggregate pickle and return CPM/logCPM matrices."""
    with open(path, "rb") as f:
        mtx_agg = pickle.load(f)
    counts_sum = mtx_agg.sum(axis=1).replace(0, np.nan)
    cpm_ref = mtx_agg.div(counts_sum, axis=0) * 1e6
    logcpm_ref = np.log1p(cpm_ref)
    return mtx_agg, cpm_ref, logcpm_ref


def TemporalPseudobulkFromAnnData(adata):
    """Aggregate raw counts by `adata.obs['time']` and return CPM/logCPM."""
    counts = adata.layers["counts"]
    genes = adata.var_names
    time_tags = adata.obs["time"]
    categories = time_tags.cat.categories if hasattr(time_tags, "cat") else sorted(pd.unique(time_tags))
    pseudobulk_counts = pd.DataFrame(index=sorted(categories), columns=genes)
    for t in categories:
        mask = (time_tags == t).values
        summed = counts[mask].sum(axis=0)
        if hasattr(summed, "A1"):
            summed = summed.A1
        pseudobulk_counts.loc[t] = summed
    pseudobulk_counts = pseudobulk_counts.astype(float)
    lib_sizes = pseudobulk_counts.sum(axis=1)
    cpm = pseudobulk_counts.div(lib_sizes, axis=0) * 1e6
    logcpm = np.log1p(cpm)
    return cpm, logcpm


def _sorted_time_union(df1, df2):
    def key(x):
        try:
            return 0, float(str(x))
        except Exception:
            return 1, str(x)

    return sorted(set(df1.index).union(df2.index), key=key)


def _to_float_time(times):
    vals = []
    try:
        for x in times:
            vals.append(float(str(x)))
        return np.asarray(vals, dtype=float), True
    except Exception:
        return np.arange(len(times), dtype=float), False


def smooth_series(y, t, method="spline", **kwargs):
    """Smooth a short time series using the methods used by the notebooks."""
    from scipy.interpolate import UnivariateSpline
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import savgol_filter

    y = np.asarray(y, dtype=float)
    t = np.asarray(t, dtype=float)
    method = str(method).lower()
    if method in ("moving_average", "ma"):
        window = int(kwargs.get("window", 3))
        return pd.Series(y).rolling(window, center=True, min_periods=1).mean().to_numpy()
    if method == "gaussian":
        return gaussian_filter1d(y, sigma=float(kwargs.get("sigma", 1.0)), mode="nearest")
    if method in ("savitzky_golay", "savgol"):
        window = int(kwargs.get("window_length", 5))
        if window % 2 == 0:
            window += 1
        window = min(window, len(y) if len(y) % 2 == 1 else len(y) - 1)
        if window < 3:
            return y.copy()
        return savgol_filter(y, window_length=window, polyorder=int(kwargs.get("polyorder", 2)), mode="interp")
    if method in ("loess", "lowess"):
        from statsmodels.nonparametric.smoothers_lowess import lowess

        return lowess(y, t, frac=float(kwargs.get("frac", 0.4)), it=int(kwargs.get("it", 0)), return_sorted=False)
    if method == "spline":
        finite = np.isfinite(y) & np.isfinite(t)
        k = int(kwargs.get("k", 3))
        if finite.sum() < (k + 1):
            return y * np.nan
        s = kwargs.get("s", None)
        if s is None:
            s = float(np.nanvar(y[finite]) * finite.sum() * 0.5)
        return UnivariateSpline(t[finite], y[finite], s=float(s), k=k)(t)
    raise ValueError(f"Unknown smoothing method: {method}")


def ComputeSmoothCorrelations(logcpm, logcpm_ref, method="spline", method_kwargs=None, local_r_threshold=0.8, local_fraction_threshold=0.5):
    """Score genes whose temporal signal is dominated by negative-control tag structure."""
    if method_kwargs is None:
        method_kwargs = dict(s=None, k=3)

    def pearsonr_nan(a, b):
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        m = np.isfinite(a) & np.isfinite(b)
        if m.sum() < 2:
            return np.nan
        x, y = a[m], b[m]
        if np.nanstd(x) == 0 or np.nanstd(y) == 0:
            return np.nan
        return float(np.corrcoef(x, y)[0, 1])

    def variance_fraction(component, total):
        component = np.asarray(component, dtype=float)
        total = np.asarray(total, dtype=float)
        m = np.isfinite(component) & np.isfinite(total)
        if m.sum() < 2:
            return np.nan
        vc, vt = np.var(component[m], ddof=1), np.var(total[m], ddof=1)
        if not np.isfinite(vt) or vt <= 0:
            return np.nan
        return float(vc / vt)

    def positive_r2(r):
        return np.nan if not np.isfinite(r) else float(max(r, 0.0) ** 2)

    def explained_fraction(r, component, total):
        vf, r2 = variance_fraction(component, total), positive_r2(r)
        if not np.isfinite(vf) or not np.isfinite(r2):
            return np.nan
        return float(r2 * vf)

    times = _sorted_time_union(logcpm, logcpm_ref)
    if len(times) == 0:
        raise ValueError("No timepoints found in either logcpm or logcpm_ref index.")
    t_float, _ = _to_float_time(times)
    results = []
    for gene in logcpm_ref.columns.intersection(logcpm.columns):
        y_ref = logcpm_ref.reindex(times)[gene].astype(float).to_numpy()
        y_exp = logcpm.reindex(times)[gene].astype(float).to_numpy()
        y_ref_s = smooth_series(y_ref, t_float, method=method, **method_kwargs)
        y_exp_s = smooth_series(y_exp, t_float, method=method, **method_kwargs)
        y_ref_resid = y_ref - y_ref_s
        y_exp_resid = y_exp - y_exp_s
        y_ref_diff, y_exp_diff = np.diff(y_ref), np.diff(y_exp)
        y_ref_smooth_diff, y_exp_smooth_diff = np.diff(y_ref_s), np.diff(y_exp_s)
        y_ref_resid_diff, y_exp_resid_diff = np.diff(y_ref_resid), np.diff(y_exp_resid)
        r_raw = pearsonr_nan(y_ref, y_exp)
        r_smooth = pearsonr_nan(y_ref_s, y_exp_s)
        r_resid = pearsonr_nan(y_ref_resid, y_exp_resid)
        r_diff = pearsonr_nan(y_ref_diff, y_exp_diff)
        r_smooth_diff = pearsonr_nan(y_ref_smooth_diff, y_exp_smooth_diff)
        r_resid_diff = pearsonr_nan(y_ref_resid_diff, y_exp_resid_diff)
        results.append(
            {
                "gene": gene,
                "pearson_r_raw": r_raw,
                "pearson_r_smooth": r_smooth,
                "pearson_r_residual": r_resid,
                "pearson_r_diff": r_diff,
                "pearson_r_smooth_diff": r_smooth_diff,
                "pearson_r_residual_diff": r_resid_diff,
                "exp_smooth_var_fraction": variance_fraction(y_exp_s, y_exp),
                "exp_residual_var_fraction": variance_fraction(y_exp_resid, y_exp),
                "ref_smooth_var_fraction": variance_fraction(y_ref_s, y_ref),
                "ref_residual_var_fraction": variance_fraction(y_ref_resid, y_ref),
                "smooth_explained_total_fraction": explained_fraction(r_smooth, y_exp_s, y_exp),
                "residual_explained_total_fraction": explained_fraction(r_resid, y_exp_resid, y_exp),
                "exp_smooth_diff_var_fraction": variance_fraction(y_exp_smooth_diff, y_exp_diff),
                "exp_residual_diff_var_fraction": variance_fraction(y_exp_resid_diff, y_exp_diff),
                "diff_smooth_explained_fraction": explained_fraction(r_smooth_diff, y_exp_smooth_diff, y_exp_diff),
                "diff_residual_explained_fraction": explained_fraction(r_resid_diff, y_exp_resid_diff, y_exp_diff),
                "diff_r2": positive_r2(r_diff),
            }
        )
    corr_df = pd.DataFrame(results).set_index("gene")
    corr_df["max_r_local"] = corr_df[["pearson_r_residual", "pearson_r_diff"]].max(axis=1, skipna=True)
    corr_df["tag_bias_candidate_residual"] = (
        (corr_df["pearson_r_residual"] >= local_r_threshold)
        & (corr_df["residual_explained_total_fraction"] >= local_fraction_threshold)
    )
    corr_df["tag_bias_candidate_diff"] = (
        (corr_df["pearson_r_diff"] >= local_r_threshold)
        & (corr_df["diff_residual_explained_fraction"] >= local_fraction_threshold)
    )
    corr_df["tag_bias_candidate"] = corr_df["tag_bias_candidate_residual"] | corr_df["tag_bias_candidate_diff"]
    corr_df["tag_bias_score"] = corr_df[["residual_explained_total_fraction", "diff_residual_explained_fraction"]].max(axis=1, skipna=True)
    return corr_df.sort_values("max_r_local", ascending=False)


def WriteTradeSeqInputsFromAnnData(tmp, key=None, tag=""):
    """Write counts, cell times, subset genes, and optional full library sizes."""
    if key is None:
        raise ValueError("key must be provided")
    outdir = f"../intmd/{key}"
    os.makedirs(outdir, exist_ok=True)
    M = tmp.layers["counts"].copy()
    if sp.issparse(M):
        M = M.tocsr()
    cells = tmp.obs_names.astype(str).to_numpy()
    genes = tmp.var_names.astype(str).to_numpy()
    M_gc = M.T if M.shape[0] == len(cells) else M
    if sp.issparse(M_gc):
        df_counts = pd.DataFrame.sparse.from_spmatrix(M_gc, index=genes, columns=cells)
    else:
        df_counts = pd.DataFrame(M_gc, index=genes, columns=cells)
    t = pd.to_numeric(tmp.obs["time"].astype(str), errors="coerce")
    df_time = pd.DataFrame({"cell_id": cells, "time": t.values})
    suffix = f".{tag}" if tag else ""
    df_counts.to_csv(f"{outdir}/counts_genes_x_cells{suffix}.csv.gz")
    df_time.to_csv(f"{outdir}/cell_times{suffix}.csv.gz", index=False)
    pd.Series(tmp.var_names.astype(str), name="gene").to_csv(f"{outdir}/subset_genes{suffix}.csv", index=False)
    if tag == "full":
        pd.DataFrame(
            {
                "cell_id": tmp.obs_names.astype(str).to_numpy(),
                "libsize_full": tmp.obs["total_counts"].to_numpy(),
            }
        ).to_csv(f"{outdir}/cell_libsize_full.csv", index=False)
    print(df_counts.shape, df_time.shape)
    print("Wrote subset_genes.csv and cell_libsize_full.csv")
