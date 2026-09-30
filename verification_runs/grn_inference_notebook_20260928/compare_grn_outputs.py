#!/usr/bin/env python
"""Compare the combined GRN notebook's outputs with the Sep 17 outputs in grn_inference/outputs."""
from pathlib import Path

import numpy as np
import pandas as pd

NEW = Path(__file__).resolve().parent / "outputs"
REF = Path("/new-stg/home/alex/tmp/ChronoPack/grn_inference/outputs")


def compare(name, a, b, key):
    a = a.drop(columns=[c for c in a.columns if c.startswith("Unnamed")])
    b = b.drop(columns=[c for c in b.columns if c.startswith("Unnamed")])
    same_cols = list(a.columns) == list(b.columns)
    same_keys = a[key].astype(str).agg("|".join, axis=1).equals(b[key].astype(str).agg("|".join, axis=1)) if len(a) == len(b) else False
    diffs = {}
    for c in a.columns.intersection(b.columns):
        if a[c].dtype.kind in "fi" and b[c].dtype.kind in "fi" and len(a) == len(b):
            x, y = a[c].to_numpy(float), b[c].to_numpy(float)
            d = np.abs(x - y) / np.maximum(np.abs(y), 1e-300)
            d[np.isnan(x) & np.isnan(y)] = 0
            m = float(np.nanmax(d)) if d.size else 0.0
            if m:
                diffs[c] = m
        elif len(a) == len(b) and not (a[c].astype(str).values == b[c].astype(str).values).all():
            diffs[c] = "values differ"
    print(f"{name:52s} rows {len(a)} vs {len(b)} | same columns {same_cols} | same {key} in same order {same_keys} | max rel diff: {diffs or 'none'}")


for label in ["TNFa", "TLR4"]:
    compare(f"omnipath_ranks.{label}.csv",
            pd.read_csv(NEW / "omnipath_tf_ranks" / f"omnipath_ranks.{label}.csv"),
            pd.read_csv(REF / "omnipath_tf_ranks" / f"omnipath_ranks.{label}.csv"), ["TF"])
for part in ["scores", "summary", "params"]:
    name = f"thp1_k562_lfdr_recalibration_{part}.csv"
    a, b = pd.read_csv(NEW / "tables" / name), pd.read_csv(REF / "tables" / name)
    key = ["dataset", "gene"] if part == "scores" else [a.columns[0]]
    compare(name, a, b, key)

new_edges = sorted((NEW / "tables" / "chronopack_grn").glob("*edges*.csv"))
print("new edge files:", [p.name for p in new_edges])
edges = pd.read_csv(new_edges[0])
legacy = edges.rename(columns={"regulator": "TF", "score": "EdgeWeight"})[["method", "TF", "target", "EdgeWeight", "rank", "percentile"]]
ref_edges = pd.read_csv(REF / "tables" / "regulator_prioritization_chronopack_ranked_edges.csv")
compare("edges vs regulator_prioritization_chronopack_ranked_edges.csv", legacy, ref_edges, ["method", "TF", "target"])
for m in sorted(set(legacy["method"]) | set(ref_edges["method"])):
    a = legacy[legacy["method"] == m]; b = ref_edges[ref_edges["method"] == m]
    ka = set(zip(a["TF"], a["target"])); kb = set(zip(b["TF"], b["target"]))
    print(f"   {m}: edges new {len(a)} ref {len(b)} | pairs only new {len(ka - kb)} only ref {len(kb - ka)} | regulators new {a['TF'].nunique()} ref {b['TF'].nunique()}")
print("COMPARE_GRN_OK")
