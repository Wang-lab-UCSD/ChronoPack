from pathlib import Path
import json

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


ROOT = Path("/new-stg/home/alex/tmp/ChronoPack/verification_runs/manuscript_reproduction_20260919")
OLD_ROOT = Path(
    "/new-stg/home/alex/tmp/ChronoPack/verification_runs/upstream_preprocess_gettraj_20260917_132244/"
    "k562-tnfa-puromycin.hg38/output"
)
NEW_ROOT = Path(
    __import__("os").environ.get("CHRONOPACK_PREPROCESS_NEW_ROOT", str(ROOT / "preprocess_puromycin/output"))
)
OUT_PATH = Path(
    __import__("os").environ.get(
        "CHRONOPACK_PREPROCESS_COMPARE_OUT", str(ROOT / "preprocess_puromycin_comparison_summary.json")
    )
)


def matrix_stats(a, b):
    if sp.issparse(a):
        a = a.toarray()
    else:
        a = np.asarray(a)
    if sp.issparse(b):
        b = b.toarray()
    else:
        b = np.asarray(b)
    if a.shape != b.shape:
        return {"same_shape": False, "shape_old": list(a.shape), "shape_new": list(b.shape)}
    diff = np.abs(a - b)
    return {
        "same_shape": True,
        "exact_equal": bool(np.array_equal(a, b)),
        "allclose_1e-8": bool(np.allclose(a, b, rtol=1e-8, atol=1e-8, equal_nan=True)),
        "max_abs_diff": float(np.nanmax(diff)) if diff.size else 0.0,
        "n_diff_gt_1e-8": int(np.sum(diff > 1e-8)),
    }


def series_equal(a, b):
    a_s = pd.Series(a).astype(object).where(pd.notna(a), "<NA>").astype(str)
    b_s = pd.Series(b).astype(object).where(pd.notna(b), "<NA>").astype(str)
    return bool(a_s.equals(b_s))


def compare_pair(name):
    old = ad.read_h5ad(OLD_ROOT / name)
    new = ad.read_h5ad(NEW_ROOT / name)
    res = {
        "shape_old": list(old.shape),
        "shape_new": list(new.shape),
        "shape_equal": old.shape == new.shape,
        "obs_names_equal": bool(old.obs_names.equals(new.obs_names)),
        "var_names_equal": bool(old.var_names.equals(new.var_names)),
        "obs_names_unique_new": bool(new.obs_names.is_unique),
        "var_names_unique_new": bool(new.var_names.is_unique),
        "layers_old": sorted(old.layers.keys()),
        "layers_new": sorted(new.layers.keys()),
        "obs_cols_old": sorted(old.obs.columns),
        "obs_cols_new": sorted(new.obs.columns),
        "var_cols_old": sorted(old.var.columns),
        "var_cols_new": sorted(new.var.columns),
        "obsm_old": sorted(old.obsm.keys()),
        "obsm_new": sorted(new.obsm.keys()),
        "obsp_old": sorted(old.obsp.keys()),
        "obsp_new": sorted(new.obsp.keys()),
        "X": matrix_stats(old.X, new.X),
    }
    res["layer_stats"] = {
        layer: matrix_stats(old.layers[layer], new.layers[layer])
        for layer in sorted(set(old.layers.keys()) & set(new.layers.keys()))
    }
    obs_stats = {}
    for col in sorted(set(old.obs.columns) & set(new.obs.columns)):
        obs_stats[col] = {"values_equal_as_strings": series_equal(old.obs[col], new.obs[col])}
        if pd.api.types.is_numeric_dtype(old.obs[col]) and pd.api.types.is_numeric_dtype(new.obs[col]):
            a = pd.to_numeric(old.obs[col], errors="coerce").to_numpy()
            b = pd.to_numeric(new.obs[col], errors="coerce").to_numpy()
            if a.dtype == bool or b.dtype == bool:
                diff = np.logical_xor(a.astype(bool), b.astype(bool))
                obs_stats[col].update(
                    {
                        "allclose_1e-8": bool(np.array_equal(a, b)),
                        "max_abs_diff": int(np.sum(diff)),
                    }
                )
                continue
            obs_stats[col].update(
                {
                    "allclose_1e-8": bool(np.allclose(a, b, rtol=1e-8, atol=1e-8, equal_nan=True)),
                    "max_abs_diff": float(np.nanmax(np.abs(a - b))) if len(a) else 0.0,
                }
            )
    res["obs_common_stats"] = obs_stats
    var_stats = {}
    for col in sorted(set(old.var.columns) & set(new.var.columns)):
        var_stats[col] = {"values_equal_as_strings": series_equal(old.var[col], new.var[col])}
        if pd.api.types.is_numeric_dtype(old.var[col]) and pd.api.types.is_numeric_dtype(new.var[col]):
            a = pd.to_numeric(old.var[col], errors="coerce").to_numpy()
            b = pd.to_numeric(new.var[col], errors="coerce").to_numpy()
            if a.dtype == bool or b.dtype == bool:
                diff = np.logical_xor(a.astype(bool), b.astype(bool))
                var_stats[col].update(
                    {
                        "allclose_1e-8": bool(np.array_equal(a, b)),
                        "max_abs_diff": int(np.sum(diff)),
                    }
                )
                continue
            var_stats[col].update(
                {
                    "allclose_1e-8": bool(np.allclose(a, b, rtol=1e-8, atol=1e-8, equal_nan=True)),
                    "max_abs_diff": float(np.nanmax(np.abs(a - b))) if len(a) else 0.0,
                }
            )
    res["var_common_stats"] = var_stats
    res["all_shared_layers_exact"] = all(v.get("exact_equal", False) for v in res["layer_stats"].values())
    res["all_common_obs_string_equal_except_doublet_score"] = all(
        v["values_equal_as_strings"] for k, v in obs_stats.items() if k != "doublet_score"
    )
    res["all_common_var_string_equal"] = all(v["values_equal_as_strings"] for v in var_stats.values())
    return res


summary = {name: compare_pair(name) for name in ["adata_full.chk0.h5ad", "adata_hvg.chk0.h5ad"]}
OUT_PATH.write_text(json.dumps(summary, indent=2, sort_keys=True))
compact = {
    k: {
        "shape_old": v["shape_old"],
        "shape_new": v["shape_new"],
        "shape_equal": v["shape_equal"],
        "obs_names_equal": v["obs_names_equal"],
        "var_names_equal": v["var_names_equal"],
        "X": v["X"],
        "layers_exact": v["all_shared_layers_exact"],
        "obs_except_doublet_equal": v["all_common_obs_string_equal_except_doublet_score"],
        "var_equal": v["all_common_var_string_equal"],
        "doublet_score_max_abs_diff": v["obs_common_stats"].get("doublet_score", {}).get("max_abs_diff"),
        "predicted_doublet_equal": v["obs_common_stats"].get("predicted_doublet", {}).get("values_equal_as_strings"),
    }
    for k, v in summary.items()
}
print(json.dumps(compact, indent=2, sort_keys=True))
