"""Compare per-context GRN outputs with the previous combined run (itself identical to Sep 17)."""
from pathlib import Path
import pandas as pd

NEW = Path(__file__).resolve().parent / "outputs"
OLD = Path("/new-stg/home/alex/tmp/ChronoPack/verification_runs/grn_inference_notebook_20260928/outputs")

def same(label, a, b):
    a = a.reset_index(drop=True); b = b.reset_index(drop=True)
    print(f"  {label:34s} {str(a.shape):14s} identical: {a.equals(b)}")
    if not a.equals(b) and a.shape == b.shape:
        for c in a.columns:
            if not a[c].equals(b[c]): print(f"      differs: {c}")

for ctx, stim in [("K562", "TNFa"), ("THP1", "TLR4")]:
    print(f"== {ctx}")
    same("tf_ranking vs omnipath_ranks." + stim, pd.read_csv(NEW / ctx / "tf_ranking.csv"),
         pd.read_csv(OLD / "omnipath_tf_ranks" / f"omnipath_ranks.{stim}.csv"))
    for part in ["scores", "summary", "params"]:
        old = pd.read_csv(OLD / "tables" / f"thp1_k562_lfdr_recalibration_{part}.csv")
        old = old[old["dataset"] == ctx].drop(columns="dataset")
        same(f"target_scores.{part}", pd.read_csv(NEW / ctx / f"target_scores.{part}.csv"), old)
    for part in ["edges", "regulators", "target_scores"]:
        old = pd.read_csv(OLD / "tables" / "chronopack_grn" / f"chronopack_grn.{part}.csv")
        old = old[old["method"] == f"ChronoPack_{ctx}"].drop(columns=[c for c in ["dataset"] if c in old.columns])
        same(f"grn.{part}", pd.read_csv(NEW / ctx / f"grn.{part}.csv"), old)
print("GRN_COMPARE_OK")
