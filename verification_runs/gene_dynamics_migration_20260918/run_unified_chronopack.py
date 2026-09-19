from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

import chronopack as cp
from chronopack._core._keys import CHRONOPACK_UNS_KEY


repo_root = Path("/new-stg/home/alex/tmp/ChronoPack")
input_dir = repo_root / "verification_runs" / "gene_dynamics_migration_20260918" / "inputs_unified_chronopack"

counts = pd.read_csv(input_dir / "counts_genes_x_cells.csv.gz", index_col=0)
times = pd.read_csv(input_dir / "cell_times.csv.gz")
adata = ad.AnnData(
    X=counts.T.to_numpy(dtype=np.float32),
    obs=pd.DataFrame(index=counts.columns.astype(str)),
    var=pd.DataFrame(index=counts.index.astype(str)),
)
adata.obs["time"] = pd.to_numeric(
    times.set_index("cell_id").loc[adata.obs_names, "time"],
    errors="raise",
).to_numpy()
adata.layers["counts"] = counts.T.to_numpy(dtype=np.int32)

result = cp.tl.fit_gene_dynamics(
    adata,
    input_dir=input_dir,
    n_points=9,
    workers=4,
    screen_nknots=4,
    final_nknots=8,
    p_threshold=0.05,
    outputs={
        "association": "unified_asso_df.1trajectory.nknots4.csv",
        "significant": "unified_sig_df.1trajectory.nknots4.csv",
        "sce": "unified_sce.1trajectory.nknots8.rds",
        "gamList": "unified_gamList.1trajectory.nknots8.rds",
    },
    verbose=True,
)

adata.write_h5ad(input_dir / "unified_gene_dynamics_smoke.h5ad")
roundtrip = ad.read_h5ad(input_dir / "unified_gene_dynamics_smoke.h5ad")

summary = {
    "backend": result.backend,
    "n_association_rows": None if result.association is None else int(result.association.shape[0]),
    "n_significant_rows": None if result.significant is None else int(result.significant.shape[0]),
    "n_selected_genes": len(result.selected_genes),
    "sidecar_paths": result.sidecar_paths,
    "h5ad_shape": list(roundtrip.shape),
    "h5ad_has_chronopack_uns": CHRONOPACK_UNS_KEY in roundtrip.uns,
    "var_columns": [col for col in adata.var.columns if col.startswith("gene_dynamics_")],
}
pd.Series(summary, dtype=object).to_json(input_dir / "unified_run_summary.json", indent=2)
print(summary)
