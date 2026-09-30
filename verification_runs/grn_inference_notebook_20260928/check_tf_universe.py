import numpy as np, pandas as pd
D = "/new-stg/home/alex/tmp/ChronoPack/verification_runs/grn_inference_notebook_20260928"
new = pd.read_csv(f"{D}/outputs/tables/chronopack_grn/chronopack_grn.edges.csv").rename(columns={"regulator": "TF", "score": "EdgeWeight"})
ref = pd.read_csv("/new-stg/home/alex/tmp/ChronoPack/grn_inference/outputs/tables/regulator_prioritization_chronopack_ranked_edges.csv")
m = ref.merge(new, on=["method", "TF", "target"], suffixes=("_ref", "_new"))
print("shared edges:", len(m), "of ref", len(ref))
print("EdgeWeight max rel diff on shared edges:", float((np.abs(m.EdgeWeight_new - m.EdgeWeight_ref) / np.abs(m.EdgeWeight_ref)).max()))
for meth, g in new.groupby("method"):
    tfs = set(ref.loc[ref.method == meth, "TF"])
    sub = g[g.TF.isin(tfs)].sort_values("EdgeWeight", ascending=False, kind="mergesort").reset_index(drop=True)
    sub["rank"] = np.arange(1, len(sub) + 1)
    r = ref[ref.method == meth].reset_index(drop=True)
    same_order = (sub[["TF", "target"]].values == r[["TF", "target"]].values).all() if len(sub) == len(r) else False
    print(f"{meth}: new restricted to the {len(tfs)} reference regulators -> {len(sub)} edges vs ref {len(r)} | same edges in same rank order {same_order}")
    print(f"   extra regulators in new: {sorted(set(g.TF) - tfs)[:12]} ...")
