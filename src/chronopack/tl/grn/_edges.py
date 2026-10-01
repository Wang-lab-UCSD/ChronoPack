"""GRN inference helpers for ChronoPack-derived TF-gene predictions."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def load_valid_targets(path: str | Path) -> set[str]:
    """Load valid protein-coding target genes, excluding MT and ENSG identifiers."""

    genes = pd.read_csv(path, header=None)[0].dropna().astype(str)
    return {gene for gene in genes if not (gene.startswith("MT-") or gene.startswith("ENSG"))}


def load_reference_tfs(
    valid_targets: set[str],
    *,
    trrust_path: str | Path,
    chip_encode_gt_path: str | Path,
) -> set[str]:
    """Load reference TFs used for regulator-prioritization universe filtering."""

    trrust = pd.read_csv(trrust_path, sep="\t", header=None)
    trrust_tfs = set(
        trrust[trrust[1].astype(str).isin(valid_targets)][0]
        .dropna()
        .astype(str)
        .str.strip()
    )

    chip_mat = pd.read_csv(chip_encode_gt_path, sep="\t", index_col=0)
    chip_tfs = set(chip_mat.columns.astype(str).str.strip())
    return trrust_tfs | chip_tfs


def normalize_edges(df: pd.DataFrame, method: str, valid_targets: set[str]) -> pd.DataFrame:
    """Normalize edge table schema and remove invalid/self edges."""

    out = df[["TF", "target", "EdgeWeight"]].copy()
    out["TF"] = out["TF"].astype(str).str.strip()
    out["target"] = out["target"].astype(str).str.strip()
    out["EdgeWeight"] = pd.to_numeric(out["EdgeWeight"], errors="coerce")
    out = out.dropna(subset=["TF", "target", "EdgeWeight"])
    out = out[out["TF"] != out["target"]]
    out = out[out["target"].isin(valid_targets)]
    out = (
        out.sort_values("EdgeWeight", ascending=False, kind="mergesort")
        .drop_duplicates(["TF", "target"], keep="first")
        .reset_index(drop=True)
    )
    out["method"] = method
    return out


def break_ties(df: pd.DataFrame) -> pd.DataFrame:
    """Apply the deterministic dense-rank tie breaker used in the benchmark script."""

    out = df.copy()
    score = out["EdgeWeight"].astype(float)
    dense = score.rank(method="dense", ascending=False)
    n_unique = int(dense.max()) if len(dense) else 0
    base = (n_unique - (dense - 1)).to_numpy(float)
    within = out.groupby("EdgeWeight", sort=False).cumcount().to_numpy(float)
    group_size = out.groupby("EdgeWeight", sort=False)["EdgeWeight"].transform("size").to_numpy(float)
    out["EdgeWeight"] = base + 1e-6 * ((group_size - 1 - within) / group_size)
    return out
