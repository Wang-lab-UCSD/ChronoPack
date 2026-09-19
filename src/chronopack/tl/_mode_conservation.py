"""Cross-context dynamic-mode conservation and switching analysis."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from chronopack._core._keys import GENE_TRAJECTORY_CLUSTERING_KEY, GENE_TRAJECTORY_KEY

CANONICAL_DYNAMIC_MODES = (
    "trend-up",
    "impulse-up_1",
    "impulse-up_2",
    "trend-down_1",
    "trend-down_2",
    "transient-down",
)
LOW_SIGNIFICANCE_LABEL = "low_significance"
MISSING_LABEL = "missing"


@dataclass
class ModeConservationResult:
    """Transparent result object for cross-context dynamic-mode comparisons."""

    directional_scores: pd.DataFrame
    conservation_scores: pd.DataFrame
    gene_matches: pd.DataFrame
    mode_transitions: pd.DataFrame
    switchers: pd.DataFrame
    parameters: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        """Return serializable DataFrame payloads keyed by result component."""
        return {
            "directional_scores": self.directional_scores.copy(),
            "conservation_scores": self.conservation_scores.copy(),
            "gene_matches": self.gene_matches.copy(),
            "mode_transitions": self.mode_transitions.copy(),
            "switchers": self.switchers.copy(),
            "parameters": dict(self.parameters),
        }

    def write(self, output_dir: str | Path, *, prefix: str = "mode_conservation") -> dict[str, str]:
        """Write result tables as CSV files and return paths."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        paths = {
            "directional_scores": output_dir / f"{prefix}.directional_scores.csv",
            "conservation_scores": output_dir / f"{prefix}.conservation_scores.csv",
            "gene_matches": output_dir / f"{prefix}.gene_matches.csv",
            "mode_transitions": output_dir / f"{prefix}.mode_transitions.csv",
            "switchers": output_dir / f"{prefix}.switchers.csv",
        }
        self.directional_scores.to_csv(paths["directional_scores"], index=False)
        self.conservation_scores.to_csv(paths["conservation_scores"], index=False)
        self.gene_matches.to_csv(paths["gene_matches"], index=False)
        self.mode_transitions.to_csv(paths["mode_transitions"], index=False)
        self.switchers.to_csv(paths["switchers"], index=False)
        return {key: str(path) for key, path in paths.items()}


def _prefixed_col(columns, key: str, name: str) -> str:
    for candidate in (f"{key}_{name}", name):
        if candidate in columns:
            return candidate
    raise KeyError(f"Could not find column '{key}_{name}' or '{name}'.")


def mode_annotation_table(
    obj,
    *,
    context: str | None = None,
    context_key: str = "cell_type",
    gene_key: str = "gene",
    mode_key: str = "broad_cluster",
    significance_key: str = "qval_1",
    rank_key: str = "waldStat_1",
    clustering_key: str = GENE_TRAJECTORY_CLUSTERING_KEY,
    fit_key: str = GENE_TRAJECTORY_KEY,
    significance_threshold: float = 0.05,
    low_significance_label: str = LOW_SIGNIFICANCE_LABEL,
    dynamic_modes: tuple[str, ...] | list[str] | None = CANONICAL_DYNAMIC_MODES,
) -> pd.DataFrame:
    """Extract standardized per-gene dynamic-mode annotations from AnnData or a DataFrame."""
    if isinstance(obj, pd.DataFrame):
        var = obj.copy()
        var.index = var.index.astype(str)
    else:
        var = obj.var.copy()
        var.index = obj.var_names.astype(str)

    gene_col = _prefixed_col(var.columns, fit_key, gene_key)
    context_col = _prefixed_col(var.columns, fit_key, context_key)
    mode_col = _prefixed_col(var.columns, clustering_key, mode_key)
    try:
        sig_col = _prefixed_col(var.columns, fit_key, significance_key)
    except KeyError:
        sig_col = _prefixed_col(var.columns, clustering_key, significance_key)
    try:
        rank_col = _prefixed_col(var.columns, fit_key, rank_key)
    except KeyError:
        rank_col = _prefixed_col(var.columns, clustering_key, rank_key)

    out = pd.DataFrame(
        {
            "row_id": var.index.astype(str),
            "gene": var[gene_col].astype(str).to_numpy(),
            "context": var[context_col].astype(str).to_numpy(),
            "mode": var[mode_col].astype(str).to_numpy(),
            "significance": pd.to_numeric(var[sig_col], errors="coerce").to_numpy(),
            "rank": pd.to_numeric(var[rank_col], errors="coerce").to_numpy(),
        }
    )
    if context is not None:
        out = out[out["context"].astype(str) == str(context)].copy()
    if out.empty:
        raise ValueError(f"No gene-mode annotations remain for context={context!r}.")

    dynamic_modes = tuple(dynamic_modes or sorted(out["mode"].dropna().astype(str).unique()))
    out["is_significant"] = out["significance"].le(float(significance_threshold))
    out["is_dynamic_mode"] = out["mode"].isin(dynamic_modes)
    out["is_high_confidence"] = out["is_significant"] & out["is_dynamic_mode"] & np.isfinite(out["rank"])
    out["mode_with_low"] = np.where(out["is_high_confidence"], out["mode"], low_significance_label)
    out = out.sort_values(["gene", "is_high_confidence", "rank", "row_id"], ascending=[True, False, False, True])
    out = out.drop_duplicates("gene", keep="first").reset_index(drop=True)
    return out


def _top_genes_by_mode(df: pd.DataFrame, *, top_k: int, rank_desc: bool = True) -> dict[str, pd.DataFrame]:
    out = {}
    for mode, sub in df[df["is_high_confidence"]].groupby("mode", sort=True):
        sub = sub.sort_values(["rank", "gene"], ascending=[not rank_desc, True])
        out[str(mode)] = sub.head(int(top_k)).copy()
    return out


def _directional_recovery(
    source: pd.DataFrame,
    target: pd.DataFrame,
    *,
    source_label: str,
    target_label: str,
    top_k: int,
    rank_desc: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    target_by_gene = target.set_index("gene", drop=False)
    rows = []
    gene_rows = []
    for mode, src_top in _top_genes_by_mode(source, top_k=top_k, rank_desc=rank_desc).items():
        genes = src_top["gene"].astype(str).tolist()
        recovered = []
        missing = []
        target_modes = []
        for gene in genes:
            if gene not in target_by_gene.index:
                recovered.append(False)
                missing.append(True)
                target_modes.append(MISSING_LABEL)
                target_high = False
            else:
                tgt = target_by_gene.loc[gene]
                target_mode = str(tgt["mode"])
                target_high = bool(tgt["is_high_confidence"])
                recovered.append(target_high and target_mode == mode)
                missing.append(False)
                target_modes.append(target_mode if target_high else LOW_SIGNIFICANCE_LABEL)
            gene_rows.append(
                {
                    "source_context": source_label,
                    "target_context": target_label,
                    "source_mode": mode,
                    "gene": gene,
                    "source_rank": float(src_top.loc[src_top["gene"] == gene, "rank"].iloc[0]),
                    "target_mode": target_modes[-1],
                    "target_high_confidence": target_high,
                    "same_mode_recovered": recovered[-1],
                    "target_missing": missing[-1],
                }
            )
        denom = len(genes)
        rows.append(
            {
                "source_context": source_label,
                "target_context": target_label,
                "mode": mode,
                "source_top_k": int(top_k),
                "source_n": int(denom),
                "same_mode_recovered_n": int(sum(recovered)),
                "directional_recovery": float(sum(recovered) / denom) if denom else np.nan,
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(gene_rows)


def _transition_table(
    source: pd.DataFrame,
    target: pd.DataFrame,
    *,
    source_label: str,
    target_label: str,
    top_k: int,
    rank_desc: bool,
    low_significance_label: str,
    missing_label: str,
    treat_missing_as_low_significance: bool,
    topk_set_mode: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    source_top_parts = list(_top_genes_by_mode(source, top_k=top_k, rank_desc=rank_desc).values())
    target_top_parts = list(_top_genes_by_mode(target, top_k=top_k, rank_desc=rank_desc).values())
    source_top = pd.concat(source_top_parts, ignore_index=True) if source_top_parts else source.iloc[0:0].copy()
    target_top = pd.concat(target_top_parts, ignore_index=True) if target_top_parts else target.iloc[0:0].copy()
    source_genes = set(source_top["gene"].astype(str))
    target_genes = set(target_top["gene"].astype(str))
    if topk_set_mode == "union":
        genes = sorted(source_genes | target_genes)
    elif topk_set_mode == "intersection":
        genes = sorted(source_genes & target_genes)
    else:
        raise ValueError("topk_set_mode must be 'union' or 'intersection'.")

    source_by_gene = source.set_index("gene", drop=False)
    target_by_gene = target.set_index("gene", drop=False)
    rows = []
    for gene in genes:
        src_missing = gene not in source_by_gene.index
        tgt_missing = gene not in target_by_gene.index
        if src_missing:
            src_mode = low_significance_label if treat_missing_as_low_significance else missing_label
            src_high = False
        else:
            src = source_by_gene.loc[gene]
            src_high = bool(src["is_high_confidence"])
            src_mode = str(src["mode"]) if src_high else low_significance_label
        if tgt_missing:
            tgt_mode = low_significance_label if treat_missing_as_low_significance else missing_label
            tgt_high = False
        else:
            tgt = target_by_gene.loc[gene]
            tgt_high = bool(tgt["is_high_confidence"])
            tgt_mode = str(tgt["mode"]) if tgt_high else low_significance_label

        if src_mode == missing_label or tgt_mode == missing_label:
            switch_class = "missing_assignment"
        elif src_mode == low_significance_label and tgt_mode == low_significance_label:
            switch_class = "low_to_low"
        elif src_mode == low_significance_label:
            switch_class = "enters_dynamic_mode"
        elif tgt_mode == low_significance_label:
            switch_class = "becomes_low_significance"
        elif src_mode == tgt_mode:
            switch_class = "same_mode"
        else:
            switch_class = "mode_switch"
        rows.append(
            {
                "gene": gene,
                "source_context": source_label,
                "target_context": target_label,
                "source_mode": src_mode,
                "target_mode": tgt_mode,
                "source_high_confidence": src_high,
                "target_high_confidence": tgt_high,
                "source_missing": src_missing,
                "target_missing": tgt_missing,
                "switch_class": switch_class,
            }
        )
    switchers = pd.DataFrame(rows)
    if switchers.empty:
        transitions = pd.DataFrame(columns=["source_mode", "target_mode", "n"])
    else:
        transitions = (
            switchers.groupby(["source_mode", "target_mode"], dropna=False)
            .size()
            .reset_index(name="n")
            .sort_values(["source_mode", "target_mode"])
            .reset_index(drop=True)
        )
    return switchers, transitions


def mode_conservation(
    source,
    target,
    *,
    source_context: str | None = None,
    target_context: str | None = None,
    source_label: str | None = None,
    target_label: str | None = None,
    gene_key: str = "gene",
    context_key: str = "cell_type",
    mode_key: str = "broad_cluster",
    significance_key: str = "qval_1",
    significance_threshold: float = 0.05,
    rank_key: str = "waldStat_1",
    top_k: int = 300,
    rank_desc: bool = True,
    clustering_key: str = GENE_TRAJECTORY_CLUSTERING_KEY,
    fit_key: str = GENE_TRAJECTORY_KEY,
    dynamic_modes: tuple[str, ...] | list[str] | None = CANONICAL_DYNAMIC_MODES,
    low_significance_label: str = LOW_SIGNIFICANCE_LABEL,
    missing_label: str = MISSING_LABEL,
    treat_missing_as_low_significance: bool = False,
    topk_set_mode: str = "union",
) -> ModeConservationResult:
    """Calculate bidirectional dynamic-mode conservation between two contexts.

    The conservation score preserves the established source/target asymmetry:
    source sets are top-ranked high-confidence genes within each source mode,
    while target recovery only requires high-confidence membership in the same
    broad mode.
    """
    source_label = source_label or source_context or "source"
    target_label = target_label or target_context or "target"
    src = mode_annotation_table(
        source,
        context=source_context,
        context_key=context_key,
        gene_key=gene_key,
        mode_key=mode_key,
        significance_key=significance_key,
        significance_threshold=significance_threshold,
        rank_key=rank_key,
        clustering_key=clustering_key,
        fit_key=fit_key,
        low_significance_label=low_significance_label,
        dynamic_modes=dynamic_modes,
    )
    tgt = mode_annotation_table(
        target,
        context=target_context,
        context_key=context_key,
        gene_key=gene_key,
        mode_key=mode_key,
        significance_key=significance_key,
        significance_threshold=significance_threshold,
        rank_key=rank_key,
        clustering_key=clustering_key,
        fit_key=fit_key,
        low_significance_label=low_significance_label,
        dynamic_modes=dynamic_modes,
    )
    s2t, s2t_genes = _directional_recovery(
        src,
        tgt,
        source_label=source_label,
        target_label=target_label,
        top_k=top_k,
        rank_desc=rank_desc,
    )
    t2s, t2s_genes = _directional_recovery(
        tgt,
        src,
        source_label=target_label,
        target_label=source_label,
        top_k=top_k,
        rank_desc=rank_desc,
    )
    directional = pd.concat([s2t, t2s], ignore_index=True)
    gene_matches = pd.concat([s2t_genes, t2s_genes], ignore_index=True)
    if directional.empty:
        conservation = pd.DataFrame(columns=["mode", "source_to_target", "target_to_source", "bidirectional_conservation"])
    else:
        tmp = directional.pivot_table(
            index="mode",
            columns=["source_context", "target_context"],
            values="directional_recovery",
            aggfunc="first",
        )
        rows = []
        for mode in sorted(directional["mode"].astype(str).unique()):
            vals = directional.loc[directional["mode"].astype(str) == mode, "directional_recovery"].dropna().to_numpy()
            rows.append(
                {
                    "mode": mode,
                    "bidirectional_conservation": float(vals.mean()) if len(vals) else np.nan,
                    "n_directions": int(len(vals)),
                }
            )
        conservation = pd.DataFrame(rows)

    switchers, transitions = _transition_table(
        src,
        tgt,
        source_label=source_label,
        target_label=target_label,
        top_k=top_k,
        rank_desc=rank_desc,
        low_significance_label=low_significance_label,
        missing_label=missing_label,
        treat_missing_as_low_significance=treat_missing_as_low_significance,
        topk_set_mode=topk_set_mode,
    )
    params = {
        "source_label": source_label,
        "target_label": target_label,
        "source_context": source_context,
        "target_context": target_context,
        "mode_key": mode_key,
        "significance_key": significance_key,
        "significance_threshold": significance_threshold,
        "rank_key": rank_key,
        "top_k": top_k,
        "rank_desc": rank_desc,
        "dynamic_modes": list(dynamic_modes or []),
        "low_significance_label": low_significance_label,
        "missing_label": missing_label,
        "treat_missing_as_low_significance": treat_missing_as_low_significance,
        "topk_set_mode": topk_set_mode,
    }
    return ModeConservationResult(
        directional_scores=directional,
        conservation_scores=conservation,
        gene_matches=gene_matches,
        mode_transitions=transitions,
        switchers=switchers,
        parameters=params,
    )


def mode_switchers(*args, **kwargs) -> ModeConservationResult:
    """Return the mode-switching tables generated by ``mode_conservation``."""
    return mode_conservation(*args, **kwargs)


__all__ = [
    "CANONICAL_DYNAMIC_MODES",
    "LOW_SIGNIFICANCE_LABEL",
    "MISSING_LABEL",
    "ModeConservationResult",
    "mode_annotation_table",
    "mode_conservation",
    "mode_switchers",
]
