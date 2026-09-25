from pathlib import Path
import json
import pandas as pd

ROOT = Path("/new-stg/home/alex/tmp/ChronoPack/verification_runs/dimred_package_migration_20260918")
PKG = ROOT / "output_pkg3"
STEP = ROOT / "output_stepwise"
OLD_REF = Path("/new-stg/home/alex/tmp/ChronoPack/verification_runs/dimred_latest_20260917_200119/output")


def compare_df(a, b):
    left = pd.read_csv(a, index_col=0)
    right = pd.read_csv(b, index_col=0)
    common_rows = left.index.intersection(right.index)
    common_cols = left.columns.intersection(right.columns)
    out = {
        "shape_left": list(left.shape),
        "shape_right": list(right.shape),
        "same_index": left.index.tolist() == right.index.tolist(),
        "same_columns": left.columns.tolist() == right.columns.tolist(),
        "row_overlap": len(common_rows),
        "column_overlap": len(common_cols),
        "numeric": {},
        "categorical": {},
    }
    for col in common_cols:
        lnum = pd.to_numeric(left.loc[common_rows, col], errors="coerce")
        rnum = pd.to_numeric(right.loc[common_rows, col], errors="coerce")
        both = lnum.notna() & rnum.notna()
        if both.any():
            delta = (lnum[both] - rnum[both]).abs()
            out["numeric"][col] = {
                "max_abs_diff": float(delta.max()),
                "median_abs_diff": float(delta.median()),
            }
        else:
            out["categorical"][col] = float(
                (left.loc[common_rows, col].astype(str) == right.loc[common_rows, col].astype(str)).mean()
            )
    return out


summary = {
    "package_vs_stepwise_pc_df": compare_df(PKG / "pc_df.new.csv", STEP / "pc_df.new.csv"),
    "package_vs_stepwise_fit": compare_df(PKG / "gene_trajectory_fit.csv", STEP / "gene_trajectory_fit.csv"),
    "package_rows_vs_old_extended_reference": {
        "package_rows": int(pd.read_csv(PKG / "pc_df.new.csv", index_col=0).shape[0]),
        "old_extended_reference_rows": int(pd.read_csv(OLD_REF / "pc_df.new.csv", index_col=0).shape[0]),
    },
}
(ROOT / "output_pkg3" / "package_vs_stepwise_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2))
