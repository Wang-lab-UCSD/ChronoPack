"""Gene-dynamics fitting orchestration."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp

from chronopack._backends.r._bridge import source_r_resource
from chronopack._core import record_provenance, validate_adata
from chronopack._core._keys import GENE_DYNAMICS_KEY
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
    return refilter_tag_bias_candidates(
        corr_df,
        local_r_threshold=local_r_threshold,
        local_fraction_threshold=local_fraction_threshold,
    )


def refilter_tag_bias_candidates(
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


def tag_bias_exclusion_genes(corr):
    """Return genes flagged by the active tag-bias rule."""
    return set(corr.index[corr["tag_bias_candidate"] == True].astype(str).tolist())


def prepare_tradeseq_inputs(
    adata,
    output_dir,
    *,
    full_adata=None,
    layer="counts",
    time_key="time",
    gene_subset=None,
    tag="",
):
    """Write the current tradeSeq CSV inputs from an AnnData object."""
    validate_adata(adata, time_key=time_key, layer=layer)
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

    if sp.issparse(matrix_gc):
        counts_df = pd.DataFrame.sparse.from_spmatrix(matrix_gc, index=genes, columns=cells)
    else:
        counts_df = pd.DataFrame(np.asarray(matrix_gc), index=genes, columns=cells)

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
    counts_df.to_csv(paths["counts"])
    time_df.to_csv(paths["cell_times"], index=False)
    pd.Series(genes, name="gene").to_csv(paths["subset_genes"], index=False)

    if tag == "full":
        lib_source = full_adata if full_adata is not None else adata
        if "total_counts" not in lib_source.obs:
            raise ChronoPackError("Writing full tradeSeq inputs requires obs['total_counts'].")
        lib_obs = lib_source.obs.loc[cells]
        lib_df = pd.DataFrame(
            {
                "cell_id": cells,
                "libsize_full": lib_obs["total_counts"].to_numpy(),
            }
        )
        paths["cell_libsize_full"] = output_dir / "cell_libsize_full.csv"
        lib_df.to_csv(paths["cell_libsize_full"], index=False)

    return {key: str(path) for key, path in paths.items()}


def _read_tradeseq_outputs(result, input_dir):
    input_dir = Path(input_dir)
    paths = {
        "association": str(result.rx2("association_path")[0]),
        "significant": str(result.rx2("significant_path")[0]),
        "sce": str(result.rx2("sce_path")[0]),
        "gamList": str(result.rx2("gamlist_path")[0]),
    }
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
):
    try:
        from rpy2.robjects import StrVector, r
    except ImportError as exc:
        raise ImportError("tradeSeq backend execution requires the optional 'rpy2' dependency.") from exc

    source_r_resource("tradeseq_helpers.R")
    r["load_chronopack_tradeseq_libraries"]()
    run = r["run_tradeseq_single_lineage"]
    return run(
        input_dir=str(input_dir),
        n_points=int(n_points),
        workers=int(workers),
        screen_nknots=int(screen_nknots),
        final_nknots=int(final_nknots),
        p_threshold=float(p_threshold),
        output_association=str(outputs["association"]),
        output_significant=str(outputs["significant"]),
        output_sce=str(outputs["sce"]),
        output_gamlist=str(outputs["gamList"]),
        association_aliases=StrVector(legacy_outputs.get("association", [])),
        significant_aliases=StrVector(legacy_outputs.get("significant", [])),
        sce_aliases=StrVector(legacy_outputs.get("sce", [])),
        gamlist_aliases=StrVector(legacy_outputs.get("gamList", [])),
        tag=str(tag),
        verbose=bool(verbose),
    )


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
    full_adata=None,
    prepare_inputs=False,
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

    Native tradeSeq objects are retained as sidecar RDS files. Scalar
    association-test outputs are written to ``adata.var`` when ``adata`` is
    supplied.
    """
    if backend != "tradeseq":
        raise ValueError("Only backend='tradeseq' is currently implemented.")

    target = adata.copy() if (copy and adata is not None) else adata
    if target is not None:
        validate_adata(target, time_key=time_key, layer=layer)

    if outputs is None:
        outputs = {
            "association": "asso_df.1trajectory.nknots4.csv",
            "significant": "sig_df.1trajectory.nknots4.csv",
            "sce": "sce.1trajectory.nknots8.rds",
            "gamList": "gamList.1trajectory.nknots8.rds",
        }
    if legacy_outputs is None:
        legacy_outputs = {
            "association": [],
            "significant": [],
            "sce": [],
            "gamList": [],
        }

    if output_dir is not None and input_dir is None:
        input_dir = output_dir
    if input_dir is None:
        raise ValueError("input_dir or output_dir must be provided.")
    input_dir = Path(input_dir)

    prepared_paths = {}
    if prepare_inputs:
        if target is None:
            raise ValueError("prepare_inputs=True requires adata.")
        prepared_paths.update(
            prepare_tradeseq_inputs(
                target,
                input_dir,
                full_adata=full_adata,
                layer=layer,
                time_key=time_key,
                tag="",
            )
        )
        full_source = full_adata if full_adata is not None else target
        prepared_paths.update(
            {
                f"full_{key}": value
                for key, value in prepare_tradeseq_inputs(
                    full_source,
                    input_dir,
                    full_adata=full_source,
                    layer=layer,
                    time_key=time_key,
                    tag="full",
                ).items()
            }
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
        "prepared_paths": prepared_paths,
    }

    association = None
    significant = None
    selected_genes = []
    sidecar_paths = {}
    backend_result = None

    if run_backend:
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


__all__ = [
    "GeneDynamicsResult",
    "compute_smooth_correlations",
    "fit_gene_dynamics",
    "load_negative_control_logcpm",
    "prepare_tradeseq_inputs",
    "refilter_tag_bias_candidates",
    "smooth_series",
    "tag_bias_exclusion_genes",
    "temporal_pseudobulk_from_anndata",
]
