from pathlib import Path

import json

import chronopack as cp


ROOT = Path("/new-stg/home/alex/tmp/ChronoPack")
SOURCE = ROOT / "verification_runs" / "filename_migration_20260917_074750"
OUT = ROOT / "verification_runs" / "dimred_package_migration_20260918" / "output_pkg3"
OUT.mkdir(parents=True, exist_ok=True)

char = cp.tl.characterize_gene_trajectories(
    root_puro=SOURCE / "k562-tnfa-puromycin.hg38",
    root_k562a=SOURCE / "k562-tnfa.hg38.rep",
    root_lps=SOURCE / "thp1-k562-lps.hg38.rep",
    output_dir=OUT,
)
cluster = cp.tl.cluster_gene_trajectories(
    characterization=char,
    output_dir=OUT,
)

char.adata.write_h5ad(OUT / "gene_trajectory_characterization.h5ad")
cluster.adata.write_h5ad(OUT / "gene_trajectory_clusters.h5ad")

summary = {
    "characterized_shape": list(char.adata.shape),
    "clustered_shape": list(cluster.adata.shape),
    "time_grid": char.metadata["time_grid"],
    "characterization_sidecars": char.sidecar_paths,
    "clustering_sidecars": cluster.sidecar_paths,
}
(OUT / "package_dimred_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
