from pathlib import Path
import json

import anndata as ad


ROOT = Path("/new-stg/home/alex/tmp/ChronoPack/verification_runs/manuscript_reproduction_20260919")
REF_ROOT = Path(
    "/new-stg/home/alex/tmp/ChronoPack/verification_runs/upstream_preprocess_gettraj_20260917_132244/"
    "k562-tnfa-puromycin.hg38/output"
)
NEW_ROOT = ROOT / "preprocess_puromycin/output"
OUT = ROOT / "preprocess_puromycin_shape_overlap.json"

summary = {}
for name in ["adata_full.chk0.h5ad", "adata_hvg.chk0.h5ad"]:
    old = ad.read_h5ad(REF_ROOT / name, backed="r")
    new = ad.read_h5ad(NEW_ROOT / name, backed="r")
    old_obs = set(old.obs_names)
    new_obs = set(new.obs_names)
    old_var = set(old.var_names)
    new_var = set(new.var_names)
    summary[name] = {
        "old_shape": list(old.shape),
        "new_shape": list(new.shape),
        "obs_overlap": len(old_obs & new_obs),
        "obs_old_only": len(old_obs - new_obs),
        "obs_new_only": len(new_obs - old_obs),
        "var_overlap": len(old_var & new_var),
        "var_old_only": len(old_var - new_var),
        "var_new_only": len(new_var - old_var),
        "first_old_only_obs": sorted(old_obs - new_obs)[:20],
        "first_new_only_obs": sorted(new_obs - old_obs)[:20],
        "obs_cols_old": sorted(old.obs.columns),
        "obs_cols_new": sorted(new.obs.columns),
        "var_cols_old": sorted(old.var.columns),
        "var_cols_new": sorted(new.var.columns),
    }
    old.file.close()
    new.file.close()

OUT.write_text(json.dumps(summary, indent=2, sort_keys=True))
print(json.dumps(summary, indent=2, sort_keys=True))
