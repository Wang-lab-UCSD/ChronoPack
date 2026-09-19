from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/new-stg/home/alex/tmp/ChronoPack")
RUN = ROOT / "verification_runs" / "dimred_package_migration_20260918" / "output_pkg3"
REF = ROOT / "verification_runs" / "dimred_latest_20260917_200119" / "output"


def numeric_summary(left: pd.DataFrame, right: pd.DataFrame):
    common_rows = left.index.intersection(right.index)
    common_cols = left.columns.intersection(right.columns)
    out = {
        "shape_left": list(left.shape),
        "shape_right": list(right.shape),
        "same_index": left.index.tolist() == right.index.tolist(),
        "same_columns": left.columns.tolist() == right.columns.tolist(),
        "row_overlap": int(len(common_rows)),
        "column_overlap": int(len(common_cols)),
        "columns": {},
    }
    for col in common_cols:
        l = pd.to_numeric(left.loc[common_rows, col], errors="coerce")
        r = pd.to_numeric(right.loc[common_rows, col], errors="coerce")
        both = l.notna() & r.notna()
        if both.any():
            delta = (l[both] - r[both]).abs()
            out["columns"][col] = {
                "max_abs_diff": float(delta.max()),
                "median_abs_diff": float(delta.median()),
            }
    return out


def categorical_match(left: pd.Series, right: pd.Series):
    common = left.index.intersection(right.index)
    return {
        "n_common": int(len(common)),
        "match_fraction": float((left.loc[common].astype(str) == right.loc[common].astype(str)).mean()),
        "left_counts": left.loc[common].astype(str).value_counts().to_dict(),
        "right_counts": right.loc[common].astype(str).value_counts().to_dict(),
    }


def main():
    ref = pd.read_csv(REF / "pc_df.new.csv", index_col=0)
    new = pd.read_csv(RUN / "pc_df.new.csv", index_col=0)
    ref_traj = pd.read_csv(REF / "traj_mat.csv", index_col=0) if (REF / "traj_mat.csv").exists() else None
    new_traj = pd.read_csv(RUN / "gene_trajectory_fit.csv", index_col=0)

    summary = {
        "pc_df": numeric_summary(ref, new),
        "broad_cluster": categorical_match(ref["broad_cluster"], new["broad_cluster"]),
    }
    for col in [
        "cluster_broad_induced",
        "cluster_broad_repressed",
        "cluster_fine_induced",
        "cluster_fine_repressed",
    ]:
        if col in ref and col in new:
            summary[col] = categorical_match(ref[col], new[col])

    # Reference traj_mat is RDS-only in this run; compare against the CSV when
    # present and otherwise report the new standardized shape.
    summary["standardized_fit_shape"] = list(new_traj.shape)

    out = RUN / "dimred_package_comparison_summary.json"
    out.write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
