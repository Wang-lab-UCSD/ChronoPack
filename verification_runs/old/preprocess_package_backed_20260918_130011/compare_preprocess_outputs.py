from pathlib import Path
import json
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

old_root = Path('/new-stg/home/alex/tmp/ChronoPack/verification_runs/upstream_preprocess_gettraj_20260917_132244/k562-tnfa-puromycin.hg38/output')
new_root = Path('/new-stg/home/alex/tmp/ChronoPack/verification_runs/preprocess_package_backed_20260918_130011/k562-tnfa-puromycin.hg38/output')
out_path = Path('/new-stg/home/alex/tmp/ChronoPack/verification_runs/preprocess_package_backed_20260918_130011/preprocess_comparison_summary.json')


def matrix_stats(a, b):
    if sp.issparse(a):
        a = a.toarray()
    else:
        a = np.asarray(a)
    if sp.issparse(b):
        b = b.toarray()
    else:
        b = np.asarray(b)
    same_shape = a.shape == b.shape
    if not same_shape:
        return {"same_shape": False, "shape_old": a.shape, "shape_new": b.shape}
    diff = np.abs(a - b)
    return {
        "same_shape": True,
        "exact_equal": bool(np.array_equal(a, b)),
        "allclose_1e-8": bool(np.allclose(a, b, rtol=1e-8, atol=1e-8, equal_nan=True)),
        "max_abs_diff": float(np.nanmax(diff)) if diff.size else 0.0,
        "n_diff_gt_1e-8": int(np.sum(diff > 1e-8)),
    }


def series_equal(a, b):
    a_s = pd.Series(a).astype(object).where(pd.notna(a), '<NA>').astype(str)
    b_s = pd.Series(b).astype(object).where(pd.notna(b), '<NA>').astype(str)
    return bool(a_s.equals(b_s))


def compare_pair(name):
    old = ad.read_h5ad(old_root / name)
    new = ad.read_h5ad(new_root / name)
    res = {
        "shape_old": old.shape,
        "shape_new": new.shape,
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
    }
    res["X"] = matrix_stats(old.X, new.X)
    layer_stats = {}
    for layer in sorted(set(old.layers.keys()).intersection(new.layers.keys())):
        layer_stats[layer] = matrix_stats(old.layers[layer], new.layers[layer])
    res["layer_stats"] = layer_stats

    obs_stats = {}
    for col in sorted(set(old.obs.columns).intersection(new.obs.columns)):
        obs_stats[col] = {"values_equal_as_strings": series_equal(old.obs[col], new.obs[col])}
        if pd.api.types.is_numeric_dtype(old.obs[col]) and pd.api.types.is_numeric_dtype(new.obs[col]):
            a = pd.to_numeric(old.obs[col], errors='coerce').to_numpy()
            b = pd.to_numeric(new.obs[col], errors='coerce').to_numpy()
            obs_stats[col].update({
                "allclose_1e-8": bool(np.allclose(a, b, rtol=1e-8, atol=1e-8, equal_nan=True)),
                "max_abs_diff": float(np.nanmax(np.abs(a - b))) if len(a) else 0.0,
            })
    res["obs_common_stats"] = obs_stats

    var_stats = {}
    for col in sorted(set(old.var.columns).intersection(new.var.columns)):
        var_stats[col] = {"values_equal_as_strings": series_equal(old.var[col], new.var[col])}
        if pd.api.types.is_numeric_dtype(old.var[col]) and pd.api.types.is_numeric_dtype(new.var[col]):
            a = pd.to_numeric(old.var[col], errors='coerce').to_numpy()
            b = pd.to_numeric(new.var[col], errors='coerce').to_numpy()
            var_stats[col].update({
                "allclose_1e-8": bool(np.allclose(a, b, rtol=1e-8, atol=1e-8, equal_nan=True)),
                "max_abs_diff": float(np.nanmax(np.abs(a - b))) if len(a) else 0.0,
            })
    res["var_common_stats"] = var_stats
    return res

summary = {
    "adata_full.chk0.h5ad": compare_pair('adata_full.chk0.h5ad'),
    "adata_hvg.chk0.h5ad": compare_pair('adata_hvg.chk0.h5ad'),
}
out_path.write_text(json.dumps(summary, indent=2, sort_keys=True))
print(json.dumps(summary, indent=2, sort_keys=True))
