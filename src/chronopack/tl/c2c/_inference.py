"""Cell-cell communication inference orchestration."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from . import _helpers


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


# Fixed conventions of the ChronoPack C2C workflow.
CELL_TYPE_KEY = "cell_type"
TIME_KEY = "time"
COUNTS_LAYER = "counts"  # detection (expressed / not expressed)
EXPRESSION_LAYER = "log1p"  # LR activity time series
ASSOCIATION_P_COL = "pvalue_1"  # time-association p-value of the ligand in the sender's own tradeSeq fit
P_PSEUDOCOUNT = 1e-3


def _detected_cells(adata) -> pd.Series:
    """Number of cells with a nonzero count for each gene."""
    X = adata.layers[COUNTS_LAYER] if COUNTS_LAYER in adata.layers else adata.X
    counts = np.asarray((X > 0).sum(axis=0)).ravel()
    return pd.Series(counts, index=adata.var_names.astype(str))


def _expressed_pairs(interactions: pd.DataFrame, sender_adata, receiver_adata, min_cells: int) -> pd.DataFrame:
    """Keep LR pairs whose ligand and receptor genes are detected in enough sender / receiver cells."""
    sender_cells = _detected_cells(sender_adata)
    receiver_cells = _detected_cells(receiver_adata)

    def detected(genes, cells):
        return bool(genes) and all(gene in cells.index and cells[gene] >= min_cells for gene in genes)

    ligand_components = interactions["ligand"].apply(_helpers.split_genes)
    receptor_components = interactions["receptor"].apply(_helpers.split_genes)
    keep = ligand_components.apply(lambda g: detected(g, sender_cells)) & receptor_components.apply(
        lambda g: detected(g, receiver_cells)
    )
    out = interactions.loc[keep].copy()
    out["lig_components"] = ligand_components.loc[keep]
    out["rec_components"] = receptor_components.loc[keep]
    return out


def _score_pairs(activities, association_df, receiver_detection, causal_window):
    """Score LR pairs by summed early activity and aggregate them into ligand scores."""
    rows = []
    for (ligand, receptor), activity_raw in activities.items():
        activity = activity_raw - activity_raw.iloc[0]
        in_window = (activity.index >= causal_window[0]) & (activity.index <= causal_window[1])
        receptor_genes = receptor.split("_")
        rows.append(
            {
                "ligand": ligand,
                "receptor": receptor,
                "pair": f"{ligand}->{receptor}",
                "activity_sum": activity.loc[in_window].sum() if in_window.any() else np.nan,
                "rec_det": float(np.nanmean([receiver_detection.get(gene, np.nan) for gene in receptor_genes])),
                "logp": -np.log10(float(association_df.loc[ligand, ASSOCIATION_P_COL]) + P_PSEUDOCOUNT),
                "n_window_points": int(in_window.sum()),
            }
        )
    pair_scores = (
        pd.DataFrame(rows).sort_values(["logp", "activity_sum"], ascending=[False, False]).reset_index(drop=True)
    )
    pair_scores["pair_rank"] = np.arange(1, len(pair_scores) + 1)

    current_ligand_rank = pair_scores.drop_duplicates("ligand", keep="first").copy()
    current_ligand_rank["current_ligand_rank"] = np.arange(1, len(current_ligand_rank) + 1)
    ligand_scores = (
        pair_scores.groupby("ligand", as_index=False)
        .agg(
            ligand_significance=("logp", "max"),
            max_lr_activity=("activity_sum", "max"),
            best_pair=("pair", "first"),
            best_pair_rank=("pair_rank", "first"),
            n_receptor_pairs=("pair", "size"),
        )
        .merge(current_ligand_rank[["ligand", "current_ligand_rank"]], on="ligand", how="left")
    )
    ligand_scores["z_ligand_significance"] = _helpers.zscore(ligand_scores["ligand_significance"])
    ligand_scores["z_max_lr_activity"] = _helpers.zscore(ligand_scores["max_lr_activity"])
    ligand_scores["composite_score"] = ligand_scores["z_ligand_significance"] + ligand_scores["z_max_lr_activity"]
    ligand_scores = ligand_scores.sort_values("composite_score", ascending=False).reset_index(drop=True)
    ligand_scores["zscore_composite_rank"] = np.arange(1, len(ligand_scores) + 1)
    return pair_scores, ligand_scores


def infer(
    *,
    adata,
    lr_resource: str | Path | pd.DataFrame,
    association: str | Path | pd.DataFrame,
    sender: str,
    receiver: str,
    min_cells: int = 1,
    causal_window: tuple[float, float] = (0.0, 45.0),
) -> C2CResult:
    """Infer temporal sender-to-receiver ligand-receptor communication.

    ``adata`` needs ``obs["cell_type"]`` and ``obs["time"]`` plus ``counts`` and
    ``log1p`` layers. ``lr_resource`` is a CellChat-style table (``ligand``,
    ``receptor``, optionally ``annotation``; only ``Secreted Signaling`` pairs are
    used). ``association`` is the sender cell type's single-lineage tradeSeq
    association table (indexed by gene, ``pvalue_1`` column), e.g. the
    ``asso_df.1trajectory.nknots4.csv`` written by ``fit_gene_dynamics``.

    LR pairs are kept when every ligand gene is detected in at least ``min_cells``
    sender cells and every receptor gene in at least ``min_cells`` receiver
    cells, and the ligand has an association row. Each pair is scored by its
    LR activity (ligand expression x receptor availability, relative to the
    first time point) summed over the time points within ``causal_window``
    (minutes, inclusive; time tags map to ``(tag - 1) * 15`` minutes), and ranked
    by ligand temporal association (-log10 p) then activity.
    """
    for key in (CELL_TYPE_KEY, TIME_KEY):
        if key not in adata.obs:
            raise KeyError(f"AnnData.obs is missing {key!r}.")

    interactions = _read_table(lr_resource)
    if "annotation" in interactions.columns:
        interactions = interactions[interactions["annotation"] == "Secreted Signaling"].copy()
    missing = {"ligand", "receptor"}.difference(interactions.columns)
    if missing:
        raise ValueError(f"LR resource is missing required columns: {sorted(missing)}")

    association_df = _read_association(association)
    sender_adata, receiver_adata = _helpers.split_by_cell_type(
        adata, cell_type_key=CELL_TYPE_KEY, sender=sender, receiver=receiver
    )
    if sender_adata.n_obs == 0:
        raise ValueError(f"No cells found for sender={sender!r}.")
    if receiver_adata.n_obs == 0:
        raise ValueError(f"No cells found for receiver={receiver!r}.")

    lr_expressed = _expressed_pairs(interactions, sender_adata, receiver_adata, min_cells)
    lr_scored = _helpers.restrict_ligands_to_association(lr_expressed, association_df)

    activities = {}
    for ligand, receptor in zip(lr_scored["ligand"], lr_scored["receptor"]):
        activity = _helpers.lr_activity(
            sender_adata, receiver_adata, ligand, receptor, time_key=TIME_KEY, layer=EXPRESSION_LAYER
        )
        if activity is not None:
            activities[(ligand, receptor)] = activity
    if not activities:
        raise ValueError(
            f"No ligand-receptor pairs to score: none pass min_cells={min_cells} detection "
            "and the ligand association filter."
        )

    receiver_detection = _helpers.detection_series(receiver_adata, layer=COUNTS_LAYER)
    pair_scores, ligand_scores = _score_pairs(activities, association_df, receiver_detection, causal_window)
    pair_scores = pair_scores.assign(sender=sender, receiver=receiver)
    front = ["sender", "receiver", "ligand", "receptor", "pair"]
    pair_scores = pair_scores[front + [col for col in pair_scores.columns if col not in front]]

    temporal_profiles = pd.DataFrame(
        [
            {
                "sender": sender,
                "receiver": receiver,
                "ligand": ligand,
                "receptor": receptor,
                "pair": f"{ligand}->{receptor}",
                "time": float(t),
                "lr_activity": float(value),
            }
            for (ligand, receptor), activity in activities.items()
            for t, value in activity.items()
        ],
        columns=["sender", "receiver", "ligand", "receptor", "pair", "time", "lr_activity"],
    )

    return C2CResult(
        interactions=pair_scores,
        ligand_scores=ligand_scores,
        temporal_profiles=temporal_profiles,
        candidate_interactions=lr_scored.reset_index(drop=True),
        detection={
            "sender": _helpers.detection_series(sender_adata, layer=COUNTS_LAYER),
            "receiver": receiver_detection,
        },
        parameters={
            "sender": sender,
            "receiver": receiver,
            "min_cells": min_cells,
            "causal_window": [float(causal_window[0]), float(causal_window[1])],
            "association_p_col": ASSOCIATION_P_COL,
            "p_pseudocount": P_PSEUDOCOUNT,
            "expression_layer": EXPRESSION_LAYER,
            "counts_layer": COUNTS_LAYER,
        },
    )


__all__ = ["C2CResult", "infer"]
