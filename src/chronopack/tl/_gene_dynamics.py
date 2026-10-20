"""Gene-dynamics fitting orchestration."""

from __future__ import annotations

import gzip
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp

from chronopack._backends.r._bridge import source_r_resource
from chronopack._core import record_provenance, validate_adata
from chronopack._core._keys import CHRONOPACK_UNS_KEY, GENE_DYNAMICS_KEY
from chronopack._core._validation import ChronoPackError


@dataclass
class GeneDynamicsResult:
    """Lightweight container for gene-dynamics backend outputs."""

    backend: str
    key_added: str
    metadata: dict[str, Any] = field(default_factory=dict)
    association: pd.DataFrame | None = None
    significant: pd.DataFrame | None = None
    selected_genes: list[str] = field(default_factory=list)
    sidecar_paths: dict[str, str] = field(default_factory=dict)
    backend_result: Any | None = None


def load_negative_control_logcpm(path):
    """Load the negative-control time-tag aggregate used by current workflows."""
    import pickle

    with open(path, "rb") as handle:
        counts = pickle.load(handle)
    # One float block: a table summed from libraries with different genes can hold
    # thousands of mixed int/float blocks, which makes every later column operation slow.
    counts = pd.DataFrame(counts.to_numpy(dtype=float), index=counts.index, columns=counts.columns)
    counts_sum = counts.sum(axis=1).replace(0, np.nan)
    cpm = counts.div(counts_sum, axis=0) * 1e6
    logcpm = np.log1p(cpm)
    return counts, cpm, logcpm


def temporal_pseudobulk_from_anndata(adata, *, layer="counts", time_key="time"):
    """Aggregate expression by time into CPM and logCPM matrices."""
    validate_adata(adata, time_key=time_key, layer=layer)
    counts = adata.layers[layer]
    genes = adata.var_names
    time_tags = adata.obs[time_key]
    categories = time_tags.cat.categories if hasattr(time_tags, "cat") else sorted(pd.unique(time_tags))
    pseudobulk_counts = pd.DataFrame(index=sorted(categories), columns=genes)
    for time_value in categories:
        mask = (time_tags == time_value).values
        summed = counts[mask].sum(axis=0)
        if hasattr(summed, "A1"):
            summed = summed.A1
        pseudobulk_counts.loc[time_value] = summed
    pseudobulk_counts = pseudobulk_counts.astype(float)
    lib_sizes = pseudobulk_counts.sum(axis=1)
    cpm = pseudobulk_counts.div(lib_sizes, axis=0) * 1e6
    logcpm = np.log1p(cpm)
    return cpm, logcpm


def _sorted_time_union(df1, df2):
    def key(value):
        try:
            return 0, float(str(value))
        except Exception:
            return 1, str(value)

    return sorted(set(df1.index).union(df2.index), key=key)


def _to_float_time(times):
    values = []
    try:
        for value in times:
            values.append(float(str(value)))
        return np.asarray(values, dtype=float), True
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


def compute_smooth_correlations(
    logcpm,
    logcpm_ref,
    method="spline",
    method_kwargs=None,
    local_r_threshold=0.8,
    local_fraction_threshold=0.5,
):
    """Score genes whose temporal signal is dominated by negative-control tags."""
    if method_kwargs is None:
        method_kwargs = dict(s=None, k=3)

    def pearsonr_nan(a, b):
        a = np.asarray(a, dtype=float)
        b = np.asarray(b, dtype=float)
        mask = np.isfinite(a) & np.isfinite(b)
        if mask.sum() < 2:
            return np.nan
        x, y = a[mask], b[mask]
        if np.nanstd(x) == 0 or np.nanstd(y) == 0:
            return np.nan
        return float(np.corrcoef(x, y)[0, 1])

    def variance_fraction(component, total):
        component = np.asarray(component, dtype=float)
        total = np.asarray(total, dtype=float)
        mask = np.isfinite(component) & np.isfinite(total)
        if mask.sum() < 2:
            return np.nan
        vc, vt = np.var(component[mask], ddof=1), np.var(total[mask], ddof=1)
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
    genes = logcpm_ref.columns.intersection(logcpm.columns)
    # Align both tables to the shared time points once, as float arrays (genes x times).
    ref_aligned = logcpm_ref.reindex(index=times, columns=genes).to_numpy(dtype=float).T
    exp_aligned = logcpm.reindex(index=times, columns=genes).to_numpy(dtype=float).T
    for gene, y_ref, y_exp in zip(genes, ref_aligned, exp_aligned):
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
    return filter_tag_bias_candidates(
        corr_df,
        local_r_threshold=local_r_threshold,
        local_fraction_threshold=local_fraction_threshold,
    )


def filter_tag_bias_candidates(
    corr,
    *,
    local_r_threshold=0.7,
    local_fraction_threshold=0.25,
    use_residual_only=True,
):
    """Apply the current post hoc tag-bias candidate rule to a correlation table."""
    corr = corr.copy()
    corr["max_r_local"] = corr[["pearson_r_residual", "pearson_r_diff"]].max(axis=1, skipna=True)
    corr["tag_bias_candidate_residual"] = (
        (corr["pearson_r_residual"] >= local_r_threshold)
        & (corr["residual_explained_total_fraction"] >= local_fraction_threshold)
    )
    corr["tag_bias_candidate_diff"] = (
        (corr["pearson_r_diff"] >= local_r_threshold)
        & (corr["diff_residual_explained_fraction"] >= local_fraction_threshold)
    )
    if use_residual_only:
        corr["tag_bias_candidate"] = corr["tag_bias_candidate_residual"]
    else:
        corr["tag_bias_candidate"] = corr["tag_bias_candidate_residual"] | corr["tag_bias_candidate_diff"]
    corr["tag_bias_score"] = corr[
        ["residual_explained_total_fraction", "diff_residual_explained_fraction"]
    ].max(axis=1, skipna=True)
    return corr.sort_values("max_r_local", ascending=False)


def refilter_tag_bias_candidates(*args, **kwargs):
    """Deprecated alias for :func:`filter_tag_bias_candidates`."""
    import warnings

    warnings.warn(
        "refilter_tag_bias_candidates is deprecated; use filter_tag_bias_candidates.",
        FutureWarning,
        stacklevel=2,
    )
    return filter_tag_bias_candidates(*args, **kwargs)


def identify_tag_bias_candidates(
    adata,
    negative_control_path,
    *,
    layer="counts",
    time_key="time",
    method="spline",
    method_kwargs=None,
    local_r_threshold=0.8,
    local_fraction_threshold=0.5,
):
    """Identify candidate genes dominated by time-tag bias.

    This combines the historical notebook sequence:
    load negative-control logCPM, compute temporal pseudobulk from AnnData,
    compute smoothed correlations, and apply the initial candidate threshold.
    """
    _, _, logcpm_ref = load_negative_control_logcpm(negative_control_path)
    _, logcpm = temporal_pseudobulk_from_anndata(adata, layer=layer, time_key=time_key)
    return compute_smooth_correlations(
        logcpm,
        logcpm_ref,
        method=method,
        method_kwargs=method_kwargs,
        local_r_threshold=local_r_threshold,
        local_fraction_threshold=local_fraction_threshold,
    )


def tag_bias_exclusion_genes(corr):
    """Return genes flagged by the active tag-bias rule."""
    return set(corr.index[corr["tag_bias_candidate"] == True].astype(str).tolist())


def _is_integer_valued(matrix):
    values = matrix.data if sp.issparse(matrix) else np.asarray(matrix)
    if np.issubdtype(values.dtype, np.integer):
        return True
    if not np.issubdtype(values.dtype, np.floating):
        return False
    return bool(np.all(np.isfinite(values)) and np.array_equal(values, np.round(values)))


def _write_genes_x_cells_csv(matrix_gc, genes, cells, path, max_values_per_chunk=20_000_000):
    """Write a genes x cells matrix as CSV; a sparse matrix is written a block of genes at a time.

    Both give the same text as ``pd.DataFrame(dense, index=genes, columns=cells).to_csv(path)``.
    """
    if not sp.issparse(matrix_gc):
        pd.DataFrame(np.asarray(matrix_gc), index=genes, columns=cells).to_csv(path)
        return
    matrix_gc = matrix_gc.tocsr()
    step = max(1, max_values_per_chunk // max(len(cells), 1))
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", newline="") as handle:
        for start in range(0, max(matrix_gc.shape[0], 1), step):
            block = matrix_gc[start : start + step].toarray()
            pd.DataFrame(block, index=genes[start : start + step], columns=cells).to_csv(handle, header=start == 0)


PRESCREEN_KEY = "time_association_prescreen"
PRESCREEN_RECORD = "prescreen_time_association.json"


def prescreen_time_association(
    adata,
    *,
    layer="counts",
    time_key="time",
    df=3,
    p_threshold=0.2,
    library_size_key="total_counts",
    chunk_size=2000,
):
    """Cheap per-gene test of time association, to choose genes for the tradeSeq fits.

    Each gene's log-normalized expression, ``log1p(count / library size * median
    library size)``, is regressed on a polynomial of degree ``df`` in time by least
    squares, and an F test against a constant gives ``pvalue``. Genes with
    ``pvalue < p_threshold`` are marked ``keep``.

    Use ``df`` = the tradeSeq screen's knots - 1 (3 for the default 4-knot screen).
    With ``df=3`` and ``p_threshold=0.2``, on the K562 + TNF-alpha and co-culture
    K562 / THP1 fits (12 time tags, or merged into 6), this kept 32-41% of genes
    and 97-99% of the genes that the tradeSeq screen called nominally significant
    (p < 0.05), including every gene ranked in its top 500. It was not validated
    for the per-lineage tests of a multi-lineage fit.

    The test is not independent of tradeSeq's, so the FDR must still count every
    screened gene. The screen is recorded in ``adata.uns``; tradeSeq inputs written
    from the subset (``fit_gene_dynamics(..., prepare_inputs=True)``) note the
    number of screened genes, and the fit's q-values are adjusted over all of them
    (genes left out count as p = 1).

    Returns a DataFrame indexed by ``adata.var_names`` with ``F``, ``df``,
    ``pvalue`` and ``keep``; subset with ``adata[:, result["keep"]]``.
    """
    from scipy import stats

    validate_adata(adata, time_key=time_key, layer=layer)
    X = adata.layers[layer] if layer is not None else adata.X
    time = pd.to_numeric(adata.obs[time_key].astype(str), errors="coerce").to_numpy(float)
    if np.isnan(time).any():
        raise ChronoPackError(f"adata.obs[{time_key!r}] contains non-numeric values.")
    n_times = np.unique(time).size
    if not 1 <= int(df) < n_times:
        raise ChronoPackError(f"df must be between 1 and the number of time points - 1 ({n_times - 1}), got {df}.")
    df = int(df)
    if library_size_key is not None and library_size_key in adata.obs:
        libsize = adata.obs[library_size_key].to_numpy(float)
    else:
        libsize = np.asarray(X.sum(axis=1)).ravel().astype(float)
    if (libsize <= 0).any():
        raise ChronoPackError("Every cell needs a positive library size.")

    n = adata.n_obs
    ts = (time - time.mean()) / time.std()
    Q, _ = np.linalg.qr(np.column_stack([ts**k for k in range(df + 1)]))
    scale = np.median(libsize) / libsize
    if sp.issparse(X):
        X = X.tocsc()
    rss0, rss1 = np.empty(adata.n_vars), np.empty(adata.n_vars)
    for start in range(0, adata.n_vars, chunk_size):
        block = X[:, start : start + chunk_size]
        block = block.toarray() if sp.issparse(block) else np.asarray(block, dtype=float)
        Y = np.log1p(block * scale[:, None])
        rss0[start : start + chunk_size] = ((Y - Y.mean(axis=0)) ** 2).sum(axis=0)
        rss1[start : start + chunk_size] = ((Y - Q @ (Q.T @ Y)) ** 2).sum(axis=0)
    df_resid = n - df - 1
    with np.errstate(divide="ignore", invalid="ignore"):
        F = np.clip(rss0 - rss1, 0, None) / df / (rss1 / df_resid)
    F = np.where(rss0 <= 0, 0.0, np.where(rss1 <= 0, np.inf, F))
    pvalue = stats.f.sf(F, df, df_resid)
    result = pd.DataFrame(
        {"F": F, "df": df, "pvalue": pvalue, "keep": pvalue < p_threshold},
        index=adata.var_names.copy(),
    )
    record_provenance(
        adata,
        PRESCREEN_KEY,
        analysis="prescreen_time_association",
        parameters={"layer": layer, "time_key": time_key, "df": df, "p_threshold": p_threshold,
                    "library_size_key": library_size_key},
        result_keys={"n_genes_screened": int(adata.n_vars), "n_genes_kept": int(result["keep"].sum())},
    )
    return result


def _write_prescreen_record(adata, output_dir, n_genes_written):
    """Note in the tradeSeq inputs how many genes a time-association pre-screen considered."""
    path = Path(output_dir) / PRESCREEN_RECORD
    entry = (adata.uns.get(CHRONOPACK_UNS_KEY) or {}).get(PRESCREEN_KEY)
    if entry is None:
        path.unlink(missing_ok=True)  # inputs rewritten without a pre-screen
        return None
    record = {
        "n_genes_screened": int(entry["result_keys"]["n_genes_screened"]),
        "n_genes_kept": int(entry["result_keys"]["n_genes_kept"]),
        "n_genes_written": int(n_genes_written),
        "parameters": {k: (v.item() if hasattr(v, "item") else v) for k, v in dict(entry["parameters"]).items()},
    }
    path.write_text(json.dumps(record, indent=2))
    return path


def _prescreened_gene_count(input_dir):
    path = Path(input_dir) / PRESCREEN_RECORD
    return int(json.loads(path.read_text())["n_genes_screened"]) if path.exists() else None


def prepare_tradeseq_inputs(
    adata,
    output_dir,
    *,
    layer="counts",
    time_key="time",
    gene_subset=None,
    cell_weights=None,
    tag="",
):
    """Write the tradeSeq CSV inputs from an AnnData object.

    Writes the ``layer`` counts (genes x cells), cell times and gene list, plus
    ``cell_libsize_full.csv`` from ``obs["total_counts"]``: the per-cell library
    size over all genes (set by ``sc.pp.calculate_qc_metrics`` before gene
    filtering; gene subsetting does not change it), which the backend turns into
    size-factor offsets. ``cell_weights`` (cells x lineages, e.g. from
    ``cp.tl.trajectory_soft_memberships``) is written to ``cell_weights.csv`` for a
    multi-lineage fit (``fit_gene_dynamics(lineage_weights="cell_weights.csv")``).
    """
    validate_adata(adata, time_key=time_key, layer=layer)
    if "total_counts" not in adata.obs:
        raise ChronoPackError("tradeSeq inputs need per-cell library sizes in obs['total_counts'].")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if gene_subset is None:
        gene_subset = adata.var_names.astype(str).to_list()
    else:
        gene_subset = [str(g) for g in gene_subset]
        missing = sorted(set(gene_subset) - set(adata.var_names.astype(str)))
        if missing:
            raise ChronoPackError(f"gene_subset contains genes absent from adata.var_names: {missing[:5]}")
        adata = adata[:, gene_subset].copy()

    matrix = adata.layers[layer]
    cells = adata.obs_names.astype(str).to_numpy()
    genes = adata.var_names.astype(str).to_numpy()
    matrix_gc = matrix.T if matrix.shape[0] == len(cells) else matrix
    # Counts are usually stored as float64; writing them as integers gives the
    # same values to the R reader (which casts to integer) but is much faster
    # to format and parse. Non-integer layers keep the float representation.
    if _is_integer_valued(matrix_gc):
        matrix_gc = matrix_gc.astype(np.int64)


    time = pd.to_numeric(adata.obs[time_key].astype(str), errors="coerce")
    if time.isna().any():
        raise ChronoPackError(f"adata.obs[{time_key!r}] contains non-numeric values.")
    time_df = pd.DataFrame({"cell_id": cells, "time": time.to_numpy()})

    suffix = f".{tag}" if tag else ""
    paths = {
        "counts": output_dir / f"counts_genes_x_cells{suffix}.csv.gz",
        "cell_times": output_dir / f"cell_times{suffix}.csv.gz",
        "subset_genes": output_dir / f"subset_genes{suffix}.csv",
    }
    _write_genes_x_cells_csv(matrix_gc, genes, cells, paths["counts"])
    time_df.to_csv(paths["cell_times"], index=False)
    pd.Series(genes, name="gene").to_csv(paths["subset_genes"], index=False)
    prescreen_path = _write_prescreen_record(adata, output_dir, len(genes))
    if prescreen_path is not None:
        paths["prescreen"] = prescreen_path

    paths["cell_libsize_full"] = output_dir / "cell_libsize_full.csv"
    pd.DataFrame({"cell_id": cells, "libsize_full": adata.obs["total_counts"].to_numpy()}).to_csv(
        paths["cell_libsize_full"], index=False
    )

    if cell_weights is not None:
        weights = pd.DataFrame(cell_weights)
        weights.index = weights.index.astype(str)
        missing = pd.Index(cells).difference(weights.index)
        if len(missing):
            raise ChronoPackError(f"cell_weights is missing {len(missing)} cells, e.g. {list(missing[:3])}.")
        weights = weights.loc[cells]
        if weights.shape[1] < 2 or (weights.to_numpy(float) < 0).any() or (weights.sum(axis=1) <= 0).any():
            raise ChronoPackError("cell_weights needs >= 2 lineage columns of non-negative weights per cell.")
        paths["cell_weights"] = output_dir / "cell_weights.csv"
        weights.to_csv(paths["cell_weights"], index_label="cell_id")

    return {key: str(path) for key, path in paths.items()}


def _read_tradeseq_outputs(result, input_dir):
    input_dir = Path(input_dir)
    def _optional_path(name):
        try:
            value = str(result.rx2(name)[0])
        except Exception:
            return None
        if value in {"NA", "NA_character_", "nan", "None"}:
            return None
        return value

    paths = {
        "association": str(result.rx2("association_path")[0]),
        "significant": str(result.rx2("significant_path")[0]),
        "sce": str(result.rx2("sce_path")[0]),
    }
    optional = {
        "sce_screen": _optional_path("sce_screen_path"),
        "sce_final": _optional_path("sce_final_path"),
        "screen_genes": _optional_path("screen_genes_path"),
        "session_info": _optional_path("session_info_path"),
    }
    paths.update({key: value for key, value in optional.items() if value})
    association = pd.read_csv(paths["association"], index_col=0)
    significant = pd.read_csv(paths["significant"], index_col=0)
    selected_genes = [str(x) for x in list(result.rx2("genes_keep"))]
    return association, significant, selected_genes, paths


def _run_tradeseq_backend(
    *,
    input_dir,
    n_points,
    workers,
    screen_nknots,
    final_nknots,
    p_threshold,
    outputs,
    legacy_outputs,
    tag,
    verbose,
    lineage_weights=None,
    fdr_n_tests=None,
):
    try:
        from rpy2.robjects import NULL, StrVector, r
    except ImportError as exc:
        raise ImportError("tradeSeq backend execution requires the optional 'rpy2' dependency.") from exc

    source_r_resource("tradeseq_helpers.R")
    r["load_chronopack_tradeseq_libraries"]()
    if lineage_weights is None:
        run, extra = r["run_tradeseq_single_lineage"], {}
    else:
        run, extra = r["run_tradeseq_weighted_lineages"], {"weights_file": str(lineage_weights)}
    return run(
        **extra,
        input_dir=str(input_dir),
        n_points=int(n_points),
        workers=int(workers),
        screen_nknots=int(screen_nknots),
        final_nknots=int(final_nknots),
        p_threshold=float(p_threshold),
        output_association=str(outputs["association"]),
        output_significant=str(outputs["significant"]),
        output_sce=str(outputs["sce"]),
        output_sce_screen=str(outputs.get("sce_screen", "sce_screen.1trajectory.nknots4.rds")),
        output_screen_genes=str(outputs.get("screen_genes", "screen_genes.1trajectory.nknots4.csv")),
        output_session_info=str(outputs.get("session_info", "sessionInfo.tradeSeqFits.txt")),
        association_aliases=StrVector(legacy_outputs.get("association", [])),
        significant_aliases=StrVector(legacy_outputs.get("significant", [])),
        sce_aliases=StrVector(legacy_outputs.get("sce", [])),
        sce_screen_aliases=StrVector(legacy_outputs.get("sce_screen", [])),
        screen_genes_aliases=StrVector(legacy_outputs.get("screen_genes", [])),
        tag=str(tag),
        fdr_n_tests=NULL if fdr_n_tests is None else int(fdr_n_tests),
        verbose=bool(verbose),
    )


def _fit_selected_genes(*, input_dir, genes, nknots, n_points, lineage_weights, workers, tag, verbose):
    import tempfile

    try:
        from rpy2.robjects import NULL, StrVector, r
    except ImportError as exc:
        raise ImportError("tradeSeq backend execution requires the optional 'rpy2' dependency.") from exc

    genes = [str(g) for g in ([genes] if isinstance(genes, str) else genes)]
    out_dir = tempfile.mkdtemp(prefix="chronopack-gene-fit-")
    source_r_resource("tradeseq_helpers.R")
    r["load_chronopack_tradeseq_libraries"]()
    result = r["fit_tradeseq_genes"](
        input_dir=str(input_dir),
        genes=StrVector(genes),
        out_dir=out_dir,
        nknots=int(nknots),
        weights_file=NULL if lineage_weights is None else str(lineage_weights),
        n_points=int(n_points),
        workers=int(workers),
        tag=str(tag),
        verbose=bool(verbose),
    )
    association_path = str(result.rx2("association_path")[0])
    association = pd.read_csv(association_path, index_col=0)
    return association, genes, {"sce": str(result.rx2("sce_path")[0]), "association": association_path}


def _write_association_to_var(adata, association, selected_genes, key_added):
    selected = set(map(str, selected_genes))
    for col in association.columns:
        out_col = f"{key_added}_{col}"
        adata.var[out_col] = np.nan
        common = adata.var_names.intersection(association.index.astype(str))
        if len(common):
            vals = association.loc[common, col]
            adata.var.loc[common, out_col] = pd.to_numeric(vals, errors="coerce").to_numpy()
    sig_col = f"{key_added}_significant"
    adata.var[sig_col] = adata.var_names.astype(str).isin(selected)


def fit_gene_dynamics(
    adata=None,
    *,
    time_key="time",
    layer="counts",
    backend="tradeseq",
    key_added=GENE_DYNAMICS_KEY,
    input_dir=None,
    output_dir=None,
    prepare_inputs=False,
    cell_weights=None,
    lineage_weights=None,
    genes=None,
    run_backend=True,
    n_points=9,
    workers=4,
    screen_nknots=4,
    final_nknots=8,
    p_threshold=0.05,
    outputs=None,
    legacy_outputs=None,
    tag="",
    verbose=True,
    copy=False,
    return_model=True,
):
    """Fit gene dynamics with the tradeSeq backend.

    By default one lineage is fitted along time: a ``screen_nknots`` screen of
    all genes with a per-gene association test, then a ``final_nknots`` refit of
    genes with p < ``p_threshold``. With ``lineage_weights`` (a cell-weights CSV in
    ``input_dir``, written by ``prepare_inputs=True`` with ``cell_weights``), one
    lineage per weight column is fitted along the same time axis, with global and
    per-lineage association tests; genes with p < ``p_threshold`` in any of them
    are refitted, and outputs are named ``*.{n}trajectory.*``.

    With ``genes``, only those genes are fitted, at ``final_nknots`` and without a
    screen (e.g. for a few genes to plot). The fit uses the same inputs, offsets and
    lineages as a full run (a gene's fit equals its fit there), is written to a
    temporary folder rather than ``input_dir``, and can be passed to
    :func:`predict_gene_dynamics`. Its association table covers only those genes,
    so q-values are adjusted over them alone. ``tag`` selects tagged inputs
    (``counts_genes_x_cells.{tag}.csv.gz``, ...); ``tag="full"`` reads the
    unfiltered count matrix, so genes removed from the tradeSeq inputs (e.g. as
    low-count) can be fitted too, on the same cells and offsets.

    If the inputs were written from genes kept by
    :func:`prescreen_time_association`, ``input_dir`` records how many genes were
    screened and the q-values are adjusted over all of them (``fdr_n_tests`` in
    the provenance), not only over the genes fitted.

    tradeSeq SingleCellExperiment objects are retained as sidecar RDS files.
    Scalar association-test outputs are written to ``adata.var`` when
    ``adata`` is supplied.
    """
    if backend != "tradeseq":
        raise ValueError("Only backend='tradeseq' is currently implemented.")

    target = adata.copy() if (copy and adata is not None) else adata
    if target is not None:
        validate_adata(target, time_key=time_key, layer=layer)

    if output_dir is not None and input_dir is None:
        input_dir = output_dir
    if input_dir is None:
        raise ValueError("input_dir or output_dir must be provided.")
    input_dir = Path(input_dir)

    if outputs is None:
        n = "1trajectory"
        if lineage_weights is not None and run_backend:
            weights_path = input_dir / lineage_weights
            if not weights_path.exists() and not prepare_inputs:
                raise FileNotFoundError(f"Missing lineage weights file {weights_path}.")
            n_lineages = (
                pd.DataFrame(cell_weights).shape[1] if cell_weights is not None
                else len(pd.read_csv(weights_path, nrows=0).columns) - 1
            )
            n = f"{n_lineages}trajectory"
        outputs = {
            "association": f"asso_df.{n}.nknots4.csv",
            "significant": f"sig_df.{n}.nknots4.csv",
            "sce_screen": f"sce_screen.{n}.nknots4.rds",
            "sce": f"sce.{n}.nknots8.rds",
            "screen_genes": f"screen_genes.{n}.nknots4.csv",
            "session_info": "sessionInfo.tradeSeqFits.txt" if n == "1trajectory" else f"sessionInfo.tradeSeqFits.{n}.txt",
        }
    if legacy_outputs is None:
        legacy_outputs = {
            "association": [],
            "significant": [],
            "sce_screen": [],
            "sce": [],
            "screen_genes": [],
        }

    prepared_paths = {}
    if prepare_inputs:
        if target is None:
            raise ValueError("prepare_inputs=True requires adata.")
        prepared_paths.update(
            prepare_tradeseq_inputs(
                target,
                input_dir,
                layer=layer,
                time_key=time_key,
                cell_weights=cell_weights,
                tag="",
            )
        )

    metadata = {
        "time_key": time_key,
        "layer": layer,
        "input_dir": str(input_dir),
        "n_points": n_points,
        "workers": workers,
        "screen_nknots": screen_nknots,
        "final_nknots": final_nknots,
        "p_threshold": p_threshold,
        "outputs": outputs,
        "legacy_outputs": legacy_outputs,
        "lineage_weights": lineage_weights,
        "prepared_paths": prepared_paths,
    }

    association = None
    significant = None
    selected_genes = []
    sidecar_paths = {}
    backend_result = None

    if run_backend and genes is not None:
        association, selected_genes, sidecar_paths = _fit_selected_genes(
            input_dir=input_dir,
            genes=genes,
            nknots=final_nknots,
            n_points=n_points,
            lineage_weights=lineage_weights,
            workers=workers,
            tag=tag,
            verbose=verbose,
        )
        metadata["genes"] = selected_genes
        if target is not None:
            _write_association_to_var(target, association, selected_genes, key_added)
    elif run_backend:
        fdr_n_tests = _prescreened_gene_count(input_dir)
        metadata["fdr_n_tests"] = fdr_n_tests
        if fdr_n_tests is not None and verbose:
            print(f"q-values are adjusted over the {fdr_n_tests} genes of the time-association pre-screen.")
        backend_result = _run_tradeseq_backend(
            input_dir=input_dir,
            n_points=n_points,
            workers=workers,
            screen_nknots=screen_nknots,
            final_nknots=final_nknots,
            p_threshold=p_threshold,
            outputs=outputs,
            legacy_outputs=legacy_outputs,
            tag=tag,
            verbose=verbose,
            lineage_weights=lineage_weights,
            fdr_n_tests=fdr_n_tests,
        )
        association, significant, selected_genes, sidecar_paths = _read_tradeseq_outputs(
            backend_result,
            input_dir,
        )
        if target is not None:
            _write_association_to_var(target, association, selected_genes, key_added)

    if target is not None:
        result_keys = {"var": []}
        if association is not None:
            result_keys["var"] = [f"{key_added}_{col}" for col in association.columns] + [
                f"{key_added}_significant"
            ]
        if sidecar_paths:
            result_keys["sidecar"] = sidecar_paths
        record_provenance(
            target,
            key_added,
            analysis="gene_dynamics",
            backend=backend,
            parameters=metadata,
            result_keys=result_keys,
        )

    result = GeneDynamicsResult(
        backend=backend,
        key_added=key_added,
        metadata=metadata,
        association=association,
        significant=significant,
        selected_genes=selected_genes,
        sidecar_paths=sidecar_paths,
        backend_result=backend_result,
    )

    if copy and target is not None:
        return (target, result) if return_model else target
    return result if return_model else None


def predict_gene_dynamics(
    fit,
    genes=None,
    *,
    lineage: int = 1,
    n_points: int = 100,
    time_max: float = 165.0,
    sce_file: str = "sce_screen.1trajectory.nknots4.rds",
) -> pd.DataFrame:
    """Fitted tradeSeq curves for selected genes on an even time grid.

    ``fit`` is either a :class:`GeneDynamicsResult` from
    ``fit_gene_dynamics(..., genes=[...])`` (``genes`` defaults to its genes) or a
    tradeSeq input directory, whose ``sce_file`` is read (by default the 4-knot
    screen fit, which covers every tested gene). Each gene is predicted for
    ``lineage`` at the reference library size, the geometric mean of the
    normalized full-library sizes. Pseudotime 0-1 is mapped to 0-``time_max``
    minutes. Returns one row per gene and time point with ``gene``, ``t``
    (minutes) and ``mu_log1p`` / ``lo_log1p`` / ``hi_log1p`` (log1p of the
    predicted mean count and its 95% band).
    """
    import tempfile

    if isinstance(fit, GeneDynamicsResult):
        input_dir = Path(fit.metadata["input_dir"])
        sce_path = Path(fit.sidecar_paths["sce"])
        genes = fit.selected_genes if genes is None else genes
    else:
        input_dir = Path(fit)
        sce_path = input_dir / sce_file
        if genes is None:
            raise ValueError("genes is required when fit is an input directory.")
    genes = [str(g) for g in ([genes] if isinstance(genes, str) else genes)]
    libsize = pd.read_csv(input_dir / "cell_libsize_full.csv")["libsize_full"].to_numpy(float)
    off_ref = float(np.mean(np.log(libsize / libsize.mean())))

    source_r_resource("dimred_plotting_helpers.R")  # predict_tradeSeq_sce_with_se()
    try:
        from rpy2.robjects import StrVector, r
    except ImportError as exc:
        raise ImportError("Predicting tradeSeq curves requires the optional 'rpy2' dependency.") from exc
    predict = r(
        """
        function(sce_path, genes, lineage, n_points, off_ref, out) {
          # rownames()/rowData() on the SCE only dispatch correctly with the package attached.
          suppressPackageStartupMessages(library(SingleCellExperiment))
          sce <- readRDS(sce_path)
          missing <- setdiff(genes, rownames(sce))
          if (length(missing)) stop("Genes not in the tradeSeq fit: ", paste(missing, collapse = ", "))
          res <- do.call(rbind, lapply(genes, function(g)
            predict_tradeSeq_sce_with_se(sce, g, lineageId = lineage, nPoints = n_points, off_ref = off_ref)))
          write.csv(res, out, row.names = FALSE)
        }
        """
    )
    with tempfile.TemporaryDirectory(prefix="chronopack-curves-") as tmp:
        out = Path(tmp) / "curves.csv"
        predict(str(sce_path), StrVector(genes), int(lineage), int(n_points), off_ref, str(out))
        curves = pd.read_csv(out)
    return pd.DataFrame(
        {
            "gene": curves["gene"],
            "t": curves["pt"] * float(time_max),
            "mu_log1p": curves["mu"],
            "lo_log1p": curves["lo"],
            "hi_log1p": curves["hi"],
        }
    )


def observed_gene_dynamics(
    fit,
    genes=None,
    *,
    lineage: int = 1,
    lineage_weights=None,
    time_max: float = 165.0,
) -> pd.DataFrame:
    """Per-cell normalized expression of selected genes, to plot under fitted curves.

    ``fit`` is a :class:`GeneDynamicsResult` (``genes`` defaults to its genes and
    ``lineage_weights`` to the fit's) or a tradeSeq input directory. Counts are
    scaled to the reference library size used by
    :func:`predict_gene_dynamics` (``count * s_ref / sf``, with ``sf`` the
    normalized full-library sizes and ``s_ref`` their geometric mean), so points
    and curves share one scale. Pseudotime 0-1 is mapped to 0-``time_max``
    minutes. Returns one row per gene and cell with ``gene``, ``cell_id``, ``t``,
    ``expr_log1p``, ``weight`` (the cell's weight in ``lineage``: 1 for a
    single-lineage fit, else its normalized ``lineage_weights`` column, a CSV in
    the input directory as in :func:`fit_gene_dynamics`) and ``assigned`` (True
    when ``lineage`` is the cell's highest-weight lineage).
    """
    if isinstance(fit, GeneDynamicsResult):
        input_dir = Path(fit.metadata["input_dir"])
        genes = fit.selected_genes if genes is None else genes
        if lineage_weights is None:
            lineage_weights = fit.metadata.get("lineage_weights")
    else:
        input_dir = Path(fit)
        if genes is None:
            raise ValueError("genes is required when fit is an input directory.")
    genes = [str(g) for g in ([genes] if isinstance(genes, str) else genes)]

    wanted = set(genes)
    chunks = []
    for chunk in pd.read_csv(input_dir / "counts_genes_x_cells.csv.gz", index_col=0, chunksize=2000):
        chunk.index = chunk.index.astype(str)
        chunks.append(chunk[chunk.index.isin(wanted)])
    counts = pd.concat(chunks)
    missing = [g for g in genes if g not in counts.index]
    if missing:
        raise ChronoPackError(f"Genes not in the tradeSeq inputs: {', '.join(missing)}")
    counts = counts.loc[genes]
    cells = counts.columns.astype(str)

    times = pd.read_csv(input_dir / "cell_times.csv.gz", index_col="cell_id")["time"]
    times.index = times.index.astype(str)
    time = times.reindex(cells).to_numpy(float)
    pt = (time - np.nanmin(time)) / (np.nanmax(time) - np.nanmin(time))
    libsize = pd.read_csv(input_dir / "cell_libsize_full.csv", index_col="cell_id")["libsize_full"]
    libsize.index = libsize.index.astype(str)
    sf = libsize.reindex(cells).to_numpy(float)
    sf = sf / sf.mean()
    scale = np.exp(np.mean(np.log(sf))) / sf

    if lineage_weights is None:
        if lineage != 1:
            raise ValueError("lineage must be 1 without lineage_weights.")
        weight, assigned = np.ones(len(cells)), np.ones(len(cells), dtype=bool)
    else:
        if isinstance(lineage_weights, pd.DataFrame):
            cw = lineage_weights.copy()
        else:
            cw = pd.read_csv(input_dir / lineage_weights, index_col=0)
        cw.index = cw.index.astype(str)
        cw = cw.reindex(cells).to_numpy(float)
        if np.isnan(cw).any():
            raise ChronoPackError("lineage_weights is missing cells of the tradeSeq inputs.")
        if not 1 <= lineage <= cw.shape[1]:
            raise ValueError(f"lineage must be between 1 and {cw.shape[1]}.")
        # As in the backend: rows normalized to 1, tiny weights treated as "not in the lineage".
        cw = cw / cw.sum(axis=1, keepdims=True)
        cw[cw < 1e-6] = 0
        weight, assigned = cw[:, lineage - 1], cw.argmax(axis=1) == lineage - 1

    expr = np.log1p(counts.to_numpy(float) * scale)
    n_genes, n_cells = expr.shape
    return pd.DataFrame(
        {
            "gene": np.repeat(genes, n_cells),
            "cell_id": np.tile(cells.to_numpy(), n_genes),
            "t": np.tile(pt * float(time_max), n_genes),
            "expr_log1p": expr.ravel(),
            "weight": np.tile(weight, n_genes),
            "assigned": np.tile(assigned, n_genes),
        }
    )


def signed_max_log2fc_from_start(
    fit,
    genes=None,
    *,
    n_points: int = 100,
    pseudocount: float = 0.01,
    sce_file: str = "sce_screen.1trajectory.nknots4.rds",
) -> pd.Series:
    """Signed maximum log2 fold change of each fitted curve versus its start.

    The tradeSeq means (``predictSmooth``, ``n_points`` grid points) are compared
    with the first grid point as log2((mu + pseudocount) / (mu_start + pseudocount));
    the value of largest magnitude is kept, with its sign (for several lineages,
    the lineage with the largest magnitude). ``fit`` is a
    :class:`GeneDynamicsResult` or a tradeSeq input directory whose ``sce_file``
    is read (by default the 4-knot screen fit, which covers every tested gene,
    matching the association table ``asso_df.1trajectory.nknots4.csv``).
    Returns a Series named ``signedMaxLog2FC_start`` indexed by gene, the x-axis
    of :func:`chronopack.pl.association_volcano`.
    """
    import tempfile

    if isinstance(fit, GeneDynamicsResult):
        sce_path = Path(fit.metadata["input_dir"]) / sce_file
        if genes is None and fit.association is not None:
            genes = fit.association.index
    else:
        sce_path = Path(fit) / sce_file
    if genes is None:
        raise ValueError("genes is required when fit is an input directory.")
    genes = [str(g) for g in ([genes] if isinstance(genes, str) else genes)]

    source_r_resource("tradeseq_helpers.R")  # signed_max_log2fc_start_many()
    try:
        from rpy2.robjects import StrVector, r
    except ImportError as exc:
        raise ImportError("Computing tradeSeq fold changes requires the optional 'rpy2' dependency.") from exc
    with tempfile.TemporaryDirectory(prefix="chronopack-log2fc-") as tmp:
        out = Path(tmp) / "log2fc.csv"
        r["signed_max_log2fc_start_rds"](str(sce_path), StrVector(genes), str(out),
                                         n_points=int(n_points), pseudocount=float(pseudocount))
        res = pd.read_csv(out, keep_default_na=False, na_values=["NA"])
    return pd.Series(res["value"].to_numpy(float), index=pd.Index(res["gene"].astype(str), name="gene"),
                     name="signedMaxLog2FC_start")



def combine_context_associations(
    asso_1: pd.DataFrame,
    asso_2: pd.DataFrame,
    *,
    wald_col: str = "waldStat_1",
    log2fc_col: str = "signedMaxLog2FC_start",
    pvalue_col: str = "pvalue_1",
    qval_col: str = "qval_1",
    qval_fun: str = "max",
    keep_only_common: bool = True,
    fill_missing: bool = False,
    arrange_by_wald_1: bool = True,
) -> pd.DataFrame:
    """Pair two contexts' single-lineage association tables by gene (VolcanoPlot ``make_combined_from_pc_df``).

    ``asso_1`` / ``asso_2`` are indexed by gene, e.g. ``asso_df.1trajectory.nknots4.csv``
    with ``signedMaxLog2FC_start`` added. Returns ``waldStat_1``, ``waldStat_2``,
    ``qval`` (``qval_fun`` = ``"max"`` or ``"min"`` of the two q-values, ignoring
    missing ones), ``log2FC_t1``, ``log2FC_t2``, ``pvalue_1``, ``pvalue_2``, indexed by
    gene, the input of :func:`chronopack.pl.response_specificity_volcano`. Genes in
    both tables are kept, or all genes with ``keep_only_common=False``;
    ``fill_missing=True`` keeps all genes and gives a gene missing from one context
    Wald 0, log2FC 0 and p- and q-value 1 there (as do missing values). Sorted by ``waldStat_1``
    (descending) unless ``arrange_by_wald_1=False``.
    """
    if qval_fun not in ("max", "min"):
        raise ValueError("qval_fun must be 'max' or 'min'.")
    cols = {pvalue_col: "pvalue", wald_col: "waldStat", log2fc_col: "log2FC_t", qval_col: "qval_tmp"}
    parts = []
    for i, asso in enumerate((asso_1, asso_2), start=1):
        missing = [c for c in cols if c not in asso]
        if missing:
            raise KeyError(f"Association table {i} is missing column(s): {', '.join(missing)}.")
        part = asso.loc[:, list(cols)].apply(pd.to_numeric, errors="coerce")
        part.columns = [f"{name}{'' if name == 'log2FC_t' else '_'}{i}" for name in cols.values()]
        part.index = part.index.astype(str)
        parts.append(part)
    how = "outer" if (fill_missing or not keep_only_common) else "inner"
    combined = parts[0].join(parts[1], how=how, sort=False)
    if how == "outer":
        combined = combined.reindex(list(dict.fromkeys([*parts[0].index, *parts[1].index])))
    if fill_missing:
        fill = {}
        for i in (1, 2):
            fill.update({f"waldStat_{i}": 0.0, f"log2FC_t{i}": 0.0, f"pvalue_{i}": 1.0, f"qval_tmp_{i}": 1.0})
        combined = combined.fillna(fill)
    q = combined[["qval_tmp_1", "qval_tmp_2"]]
    combined["qval"] = q.max(axis=1) if qval_fun == "max" else q.min(axis=1)
    combined = combined[["waldStat_1", "waldStat_2", "qval", "log2FC_t1", "log2FC_t2", "pvalue_1", "pvalue_2"]]
    if arrange_by_wald_1:
        combined = combined.sort_values("waldStat_1", ascending=False, kind="stable")
    combined.index.name = "gene"
    return combined

__all__ = [
    "GeneDynamicsResult",
    "combine_context_associations",
    "compute_smooth_correlations",
    "filter_tag_bias_candidates",
    "fit_gene_dynamics",
    "identify_tag_bias_candidates",
    "load_negative_control_logcpm",
    "predict_gene_dynamics",
    "prepare_tradeseq_inputs",
    "refilter_tag_bias_candidates",
    "signed_max_log2fc_from_start",
    "smooth_series",
    "tag_bias_exclusion_genes",
    "temporal_pseudobulk_from_anndata",
]
