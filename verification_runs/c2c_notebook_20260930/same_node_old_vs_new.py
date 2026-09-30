"""Run the previous (git HEAD) and current cp.tl.c2c.infer on the same node and compare in memory."""
import importlib.util, platform, subprocess
from pathlib import Path
import scanpy as sc
import chronopack as cp

here = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("old_c2c_inference", here / "old_c2c_inference.py")
old = importlib.util.module_from_spec(spec)
import sys; sys.modules[spec.name] = old
spec.loader.exec_module(old)

DATA = Path("/new-stg/home/alex/projects/chrono_seq/analyses/thp1-k562-lps.hg38.rep")
adata = sc.read_h5ad(DATA / "output/adata_full.chk1.h5ad")
common = dict(adata=adata, lr_resource="/new-stg/home/alex/ref/cellchatdb/hs/interaction_input_CellChatDB.csv",
              association=DATA / "intmd/sc_to_tradeSeq/asso_df.both.2trajectory.nknots4.expanded.csv", sender="THP1", receiver="K562")
res_old = old.infer(**common, response_genes=["NFKBIA", "NFKBIZ", "TNFAIP3", "BIRC3", "IER3", "CXCL1", "CXCL2", "CXCL3"],
                    min_detect=1e-100, tradeseq_export_dir=DATA / "intmd/sc_to_tradeSeq/ts_export_grid_100pts")
res_new = cp.tl.c2c.infer(**common)
print("node:", platform.node(), "| cpu:", subprocess.run("grep -m1 'model name' /proc/cpuinfo", shell=True, capture_output=True, text=True).stdout.split(":")[-1].strip())
for t in ["interactions", "ligand_scores", "candidate_interactions"]:
    print(f"{t:24s} identical old vs new (same node): {getattr(res_old, t).equals(getattr(res_new, t))}")
print("temporal_profiles (minus receiver_response):", res_old.temporal_profiles.drop(columns="receiver_response").equals(res_new.temporal_profiles))
