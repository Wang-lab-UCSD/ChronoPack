"""Run cp.tl.c2c.infer (summed window activity) and compare rankings to the mean-based reference outputs."""
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

import chronopack as cp

HERE = Path(__file__).resolve().parent
DATA = Path("/new-stg/home/alex/projects/chrono_seq/analyses/thp1-k562-lps.hg38.rep")
adata = sc.read_h5ad(DATA / "output/adata_full.chk1.h5ad")

c2c = cp.tl.c2c.infer(
    adata=adata,
    lr_resource="/new-stg/home/alex/ref/cellchatdb/hs/interaction_input_CellChatDB.csv",
    association=DATA / "intmd/sc_to_tradeSeq/asso_df.both.2trajectory.nknots4.expanded.csv",
    sender="THP1",
    receiver="K562",
    causal_window=(0, 45),
)
c2c.write(HERE / "outputs/chronopack_c2c")
fig, ax = cp.pl.c2c.interaction_scatter(c2c, x_thresh=0.4, xlim=(-0.6, 1.2))
fig.savefig(HERE / "outputs/chronopack_lr_pair_scores.png", dpi=300, bbox_inches="tight")
fig.savefig(HERE / "outputs/chronopack_lr_pair_scores.pdf", bbox_inches="tight")

ref_pairs = pd.read_csv(HERE / "outputs/reference_mean.interactions.csv")
ref_lig = pd.read_csv(HERE / "outputs/reference_mean.ligand_scores.csv")
new_pairs, new_lig = c2c.interactions, c2c.ligand_scores

print("n pairs ref/new:", len(ref_pairs), len(new_pairs))
print("n_window_points values:", sorted(new_pairs["n_window_points"].unique()))
print("pair order identical:", ref_pairs["pair"].tolist() == new_pairs["pair"].tolist())
m = ref_pairs.merge(new_pairs, on="pair", suffixes=("_ref", "_new"))
print("max |sum - 4*mean|:", float(np.nanmax(np.abs(m["activity_sum"] - m["n_window_points_new"] * m["avg"]))))
print("pair_rank identical:", bool((m["pair_rank_ref"] == m["pair_rank_new"]).all()))
print("logp identical:", bool(np.allclose(m["logp_ref"], m["logp_new"])))
print("ligand order identical:", ref_lig["ligand"].tolist() == new_lig["ligand"].tolist())
ml = ref_lig.merge(new_lig, on="ligand", suffixes=("_ref", "_new"))
for col in ["zscore_composite_rank", "current_ligand_rank", "best_pair_rank"]:
    print(f"{col} identical:", bool((ml[f"{col}_ref"] == ml[f"{col}_new"]).all()))
for col in ["z_max_lr_activity", "composite_score"]:
    print(f"max |{col} diff|:", float(np.nanmax(np.abs(ml[f"{col}_ref"] - ml[f"{col}_new"]))))
print("top 10 ligands:", new_lig["ligand"].head(10).tolist())
print("pairs labeled (p<0.05 & sum>=0.4):", new_pairs.loc[(new_pairs.logp > -np.log10(0.05)) & (new_pairs.activity_sum >= 0.4), "pair"].tolist())
print("CHRONOPACK_STAGE_OK")
