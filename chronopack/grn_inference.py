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


def load_chronopack_beta_lfdr(
    *,
    name: str,
    dataset: str,
    gene_scores_path: str | Path,
    omnipath_tf_scores_path: str | Path,
    valid_targets: set[str],
) -> pd.DataFrame:
    """Construct ChronoPack TF-gene edges from BUM LFDR gene scores and TF priorities."""

    scores = pd.read_csv(gene_scores_path)
    genes = scores[scores["dataset"] == dataset][["gene", "lfdr_bum"]].copy()
    genes["target"] = genes["gene"].astype(str)
    genes["gene_score"] = 1 - pd.to_numeric(genes["lfdr_bum"], errors="coerce")
    genes = genes[genes["target"].isin(valid_targets)]
    genes = genes[genes["gene_score"] > 0].dropna(subset=["gene_score"])

    omni = pd.read_csv(omnipath_tf_scores_path)
    omni["TF"] = omni["TF"].astype(str).str.strip()
    omni["priority_act"] = pd.to_numeric(omni["priority_act"], errors="coerce")
    omni = omni[omni["priority_act"] > 0].dropna(subset=["priority_act"])
    omni["tf_score"] = omni["priority_act"] / omni["priority_act"].sum()

    edge_weight = np.outer(omni["tf_score"].to_numpy(float), genes["gene_score"].to_numpy(float))
    pred = (
        pd.DataFrame(edge_weight, index=omni["TF"], columns=genes["target"])
        .stack()
        .rename("EdgeWeight")
        .reset_index()
        .rename(columns={"level_0": "TF", "level_1": "target"})
    )
    return break_ties(normalize_edges(pred, name, valid_targets))


def build_chronopack_ranked_edges(
    *,
    gene_scores_path: str | Path,
    k562_tf_scores_path: str | Path,
    thp1_tf_scores_path: str | Path,
    valid_targets_path: str | Path,
    trrust_path: str | Path | None = None,
    chip_encode_gt_path: str | Path | None = None,
) -> pd.DataFrame:
    """Build method-specific ranked ChronoPack edges for K562 and THP1 contexts."""

    valid_targets = load_valid_targets(valid_targets_path)
    methods = {
        "ChronoPack_K562": load_chronopack_beta_lfdr(
            name="ChronoPack_K562",
            dataset="K562",
            gene_scores_path=gene_scores_path,
            omnipath_tf_scores_path=k562_tf_scores_path,
            valid_targets=valid_targets,
        ),
        "ChronoPack_THP1": load_chronopack_beta_lfdr(
            name="ChronoPack_THP1",
            dataset="THP1",
            gene_scores_path=gene_scores_path,
            omnipath_tf_scores_path=thp1_tf_scores_path,
            valid_targets=valid_targets,
        ),
    }

    reference_tfs = None
    if trrust_path is not None and chip_encode_gt_path is not None:
        reference_tfs = load_reference_tfs(
            valid_targets,
            trrust_path=trrust_path,
            chip_encode_gt_path=chip_encode_gt_path,
        )

    ranked = []
    for method, edges in methods.items():
        pred = edges
        if reference_tfs is not None:
            pred = pred[pred["TF"].isin(reference_tfs)].copy()
        pred = (
            pred.sort_values("EdgeWeight", ascending=False, kind="mergesort")
            .drop_duplicates(["TF", "target"], keep="first")
            .reset_index(drop=True)
        )
        pred["rank"] = np.arange(1, len(pred) + 1)
        pred["percentile"] = 1 - (pred["rank"] - 1) / max(len(pred), 1)
        ranked.append(pred[["method", "TF", "target", "EdgeWeight", "rank", "percentile"]])

    return pd.concat(ranked, ignore_index=True)


def write_chronopack_ranked_edges(ranked_edges: pd.DataFrame, path: str | Path) -> Path:
    """Write ChronoPack method-specific ranked edges."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ranked_edges.to_csv(path, index=False)
    return path
