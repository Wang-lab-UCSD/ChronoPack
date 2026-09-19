"""Cell-cell communication inference orchestration."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from chronopack._legacy import c2c as _legacy


@dataclass
class C2CResult:
    """Transparent result object for temporal cell-cell communication inference."""

    interactions: pd.DataFrame
    ligand_scores: pd.DataFrame
    temporal_profiles: pd.DataFrame | None = None
    candidate_interactions: pd.DataFrame | None = None
    detection: dict[str, pd.Series] = field(default_factory=dict)
    parameters: dict[str, Any] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)

    def legacy_pair_scores(self, *, interval_min: float | None = None, include_rank: bool = False) -> pd.DataFrame:
        """Return the historical LR pair score table schema."""
        cols = ["ligand", "receptor", "pair", "avg", "rec_det", "logp", "n_window_points"]
        if include_rank and "pair_rank" in self.interactions.columns:
            cols.append("pair_rank")
        out = self.interactions[cols].copy()
        if interval_min is not None:
            out.insert(0, "interval_min", interval_min)
        return out

    def write(self, output_dir: str | Path, *, prefix: str = "chronopack_c2c") -> dict[str, str]:
        """Write result tables and small JSON metadata."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        paths = {
            "interactions_csv": output_dir / f"{prefix}.interactions.csv",
            "ligand_scores_csv": output_dir / f"{prefix}.ligand_scores.csv",
            "metadata": output_dir / f"{prefix}.metadata.json",
        }
        self.interactions.to_csv(paths["interactions_csv"], index=False)
        self.ligand_scores.to_csv(paths["ligand_scores_csv"], index=False)
        try:
            parquet_path = output_dir / f"{prefix}.interactions.parquet"
            self.interactions.to_parquet(parquet_path, index=False)
            paths["interactions_parquet"] = parquet_path
        except Exception:
            pass
        if self.temporal_profiles is not None:
            paths["temporal_profiles_csv"] = output_dir / f"{prefix}.temporal_profiles.csv"
            self.temporal_profiles.to_csv(paths["temporal_profiles_csv"], index=False)
        if self.candidate_interactions is not None:
            paths["candidate_interactions_csv"] = output_dir / f"{prefix}.candidate_interactions.csv"
            self.candidate_interactions.to_csv(paths["candidate_interactions_csv"], index=False)
        metadata = {"parameters": self.parameters, "provenance": self.provenance}
        paths["metadata"].write_text(json.dumps(metadata, indent=2, sort_keys=True))
        return {key: str(path) for key, path in paths.items()}


def _read_table(obj: str | Path | pd.DataFrame) -> pd.DataFrame:
    if isinstance(obj, pd.DataFrame):
        return obj.copy()
    path = Path(obj)
    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path)


def _read_association(obj: str | Path | pd.DataFrame) -> pd.DataFrame:
    if isinstance(obj, pd.DataFrame):
        return obj.copy()
    return pd.read_csv(obj, index_col=0)


def _profiles_from_pair_series(pair_series: Mapping[tuple[str, str], tuple[pd.Series, pd.Series]]) -> pd.DataFrame:
    rows = []
    for (ligand, receptor), (activity, response) in pair_series.items():
        common = activity.index.intersection(response.index)
        for t in common:
            rows.append(
                {
                    "ligand": ligand,
                    "receptor": receptor,
                    "pair": f"{ligand}->{receptor}",
                    "time": float(t),
                    "lr_activity": float(activity.loc[t]),
                    "receiver_response": float(response.loc[t]),
                }
            )
    return pd.DataFrame(rows)


def infer(
    *,
    adata,
    lr_resource: str | Path | pd.DataFrame,
    association: str | Path | pd.DataFrame,
    sender: str,
    receiver: str,
    response_genes: Sequence[str],
    cell_type_key: str = "cell_type",
    time_key: str = "time",
    counts_layer: str = "counts",
    expression_layer: str = "log1p",
    min_detect: float = 0.0,
    causal_window_min: tuple[float, float] = (0.0, 45.0),
    association_p_col: str = "pvalue_2",
    p_pseudocount: float = 1e-3,
    prefer_fits: bool = False,
    fits_df: str | Path | pd.DataFrame | None = None,
    wmeans_df: str | Path | pd.DataFrame | None = None,
    tradeseq_export_dir: str | Path | None = None,
    secreted_only: bool = True,
    provenance: Mapping[str, Any] | None = None,
) -> C2CResult:
    """Infer temporal sender-to-receiver ligand-receptor communication.

    This function preserves the current ChronoPack C2C workflow: secreted LR
    catalog filtering, sender/receiver expression eligibility, ligand temporal
    association restriction, LR activity construction, and early-window pair
    scoring.
    """
    if cell_type_key not in adata.obs:
        raise KeyError(f"AnnData.obs is missing cell type key {cell_type_key!r}.")
    if time_key not in adata.obs:
        raise KeyError(f"AnnData.obs is missing time key {time_key!r}.")

    interactions = _read_table(lr_resource)
    if secreted_only and "annotation" in interactions.columns:
        interactions = interactions[interactions["annotation"] == "Secreted Signaling"].copy()
    required = {"ligand", "receptor"}
    missing = required.difference(interactions.columns)
    if missing:
        raise ValueError(f"LR resource is missing required columns: {sorted(missing)}")

    association_df = _read_association(association)
    sender_adata, receiver_adata = _legacy.split_by_cell_type(
        adata,
        cell_type_key=cell_type_key,
        sender=sender,
        receiver=receiver,
    )
    if sender_adata.n_obs == 0:
        raise ValueError(f"No cells found for sender={sender!r}.")
    if receiver_adata.n_obs == 0:
        raise ValueError(f"No cells found for receiver={receiver!r}.")

    if tradeseq_export_dir is not None:
        exports = _legacy.load_tradeseq_grid_exports(tradeseq_export_dir)
        fits_df = exports["fits"] if fits_df is None else fits_df
        wmeans_df = exports["wmeans"] if wmeans_df is None else wmeans_df
    fits = _read_table(fits_df) if fits_df is not None else None
    wmeans = _read_table(wmeans_df) if wmeans_df is not None else None

    lr_expressed, sender_detection, receiver_detection = _legacy.filter_expressed_secreted_lr_catalog(
        interactions,
        sender_adata=sender_adata,
        receiver_adata=receiver_adata,
        sender_layer=counts_layer,
        receiver_layer=counts_layer,
        min_detect=min_detect,
    )
    lr_scored = _legacy.restrict_ligands_to_association(lr_expressed, association_df)
    pair_series = _legacy.build_pair_series(
        lr_scored,
        sender_adata=sender_adata,
        receiver_adata=receiver_adata,
        response_genes=response_genes,
        time_key=time_key,
        layer=expression_layer,
        prefer_fits=prefer_fits,
        fits_df=fits,
        wmeans_df=wmeans,
    )
    pair_scores, ligand_scores = _legacy.score_lr_pairs(
        pair_series,
        association_df=association_df,
        receiver_detection=receiver_detection,
        causal_window_min=causal_window_min,
        association_p_col=association_p_col,
        p_pseudocount=p_pseudocount,
    )
    pair_scores = pair_scores.assign(sender=sender, receiver=receiver)
    pair_scores["score"] = pair_scores["avg"]
    front = ["sender", "receiver", "ligand", "receptor", "pair", "score"]
    pair_scores = pair_scores[front + [col for col in pair_scores.columns if col not in front]]
    temporal_profiles = _profiles_from_pair_series(pair_series)
    if not temporal_profiles.empty:
        temporal_profiles = temporal_profiles.assign(sender=sender, receiver=receiver)
        temporal_profiles = temporal_profiles[
            ["sender", "receiver", "ligand", "receptor", "pair", "time", "lr_activity", "receiver_response"]
        ]

    return C2CResult(
        interactions=pair_scores,
        ligand_scores=ligand_scores,
        temporal_profiles=temporal_profiles,
        candidate_interactions=lr_scored.reset_index(drop=True),
        detection={"sender": sender_detection, "receiver": receiver_detection},
        parameters={
            "sender": sender,
            "receiver": receiver,
            "cell_type_key": cell_type_key,
            "time_key": time_key,
            "counts_layer": counts_layer,
            "expression_layer": expression_layer,
            "min_detect": min_detect,
            "causal_window_min": tuple(causal_window_min),
            "association_p_col": association_p_col,
            "p_pseudocount": p_pseudocount,
            "prefer_fits": prefer_fits,
            "secreted_only": secreted_only,
            "response_genes": list(response_genes),
        },
        provenance=dict(provenance or {}),
    )


__all__ = ["C2CResult", "infer"]
