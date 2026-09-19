"""ChronoPack GRN inference helpers."""

from __future__ import annotations

import json
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import pandas as pd

from chronopack._backends.r._bridge import source_r_resource
from chronopack._legacy import grn_inference as _legacy_edges
from chronopack._legacy.omnipath_tf import (
    OmniPathTfGraph,
    OmniPathTfResources,
    build_omnipath_tf_graph,
    load_omnipath_tf_resources,
    rank_omnipath_tfs,
    run_tf_diffusion_ranking,
)


@dataclass
class TargetGeneScoreResult:
    """Target-gene scores produced by the fdrtool/BUM recalibration workflow."""

    scores: pd.DataFrame
    summary: pd.DataFrame
    params: pd.DataFrame
    parameters: dict[str, Any] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)
    sidecar_paths: dict[str, str] = field(default_factory=dict)

    def write(self, output_dir: str | Path, *, prefix: str = "target_gene_scores") -> dict[str, str]:
        """Write result tables as CSV files."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        paths = {
            "scores": output_dir / f"{prefix}.scores.csv",
            "summary": output_dir / f"{prefix}.summary.csv",
            "params": output_dir / f"{prefix}.params.csv",
        }
        self.scores.to_csv(paths["scores"], index=False)
        self.summary.to_csv(paths["summary"], index=False)
        self.params.to_csv(paths["params"], index=False)
        return {key: str(path) for key, path in paths.items()}


@dataclass
class TFRankingResult:
    """OmniPath TF-prioritization result."""

    rankings: dict[str, pd.DataFrame]
    histories: dict[str, pd.DataFrame] = field(default_factory=dict)
    parameters: dict[str, Any] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)

    def write(self, output_dir: str | Path, *, prefix: str = "omnipath_ranks") -> dict[str, str]:
        """Write one ranking CSV per context."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        paths = {}
        for label, df in self.rankings.items():
            path = output_dir / f"{prefix}.{label}.csv"
            df.to_csv(path, index=False)
            paths[label] = str(path)
        return paths


@dataclass
class GRNResult:
    """Transparent result object for ChronoPack TF-target GRN inference."""

    edges: pd.DataFrame
    regulators: pd.DataFrame | None = None
    target_scores: pd.DataFrame | None = None
    parameters: dict[str, Any] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)

    def to_legacy_edges(self) -> pd.DataFrame:
        """Return the historical edge-table schema."""
        edges = self.edges.copy()
        rename = {}
        if "regulator" in edges.columns:
            rename["regulator"] = "TF"
        if "score" in edges.columns:
            rename["score"] = "EdgeWeight"
        return edges.rename(columns=rename)

    def write_edges(self, path: str | Path, *, legacy_schema: bool = False) -> str:
        """Write the edge table to CSV or Parquet based on the filename suffix."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        edges = self.to_legacy_edges() if legacy_schema else self.edges
        if path.suffix.lower() == ".parquet":
            edges.to_parquet(path, index=False)
        else:
            edges.to_csv(path, index=False)
        return str(path)

    def write(self, output_dir: str | Path, *, prefix: str = "chronopack_grn") -> dict[str, str]:
        """Write language-neutral result tables and small JSON metadata."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        paths = {
            "edges_csv": output_dir / f"{prefix}.edges.csv",
            "metadata": output_dir / f"{prefix}.metadata.json",
        }
        self.edges.to_csv(paths["edges_csv"], index=False)
        try:
            parquet_path = output_dir / f"{prefix}.edges.parquet"
            self.edges.to_parquet(parquet_path, index=False)
            paths["edges_parquet"] = parquet_path
        except Exception:
            pass
        if self.regulators is not None:
            paths["regulators"] = output_dir / f"{prefix}.regulators.csv"
            self.regulators.to_csv(paths["regulators"], index=False)
        if self.target_scores is not None:
            paths["target_scores"] = output_dir / f"{prefix}.target_scores.csv"
            self.target_scores.to_csv(paths["target_scores"], index=False)
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


def _load_valid_targets(valid_targets: str | Path | Iterable[str] | pd.DataFrame | None) -> set[str] | None:
    if valid_targets is None:
        return None
    if isinstance(valid_targets, (str, Path)):
        return _legacy_edges.load_valid_targets(valid_targets)
    if isinstance(valid_targets, pd.DataFrame):
        values = valid_targets.iloc[:, 0].astype(str)
    else:
        values = pd.Series(list(valid_targets), dtype=str)
    return {gene for gene in values if not gene.startswith("MT-") and not gene.startswith("ENSG")}


def _normalize_contexts(
    tf_scores: Mapping[str, str | Path | pd.DataFrame],
    contexts: Mapping[str, str | Mapping[str, Any]] | None,
) -> dict[str, dict[str, Any]]:
    if contexts is None:
        return {label: {"dataset": label, "tf_scores": scores} for label, scores in tf_scores.items()}
    out = {}
    for label, spec in contexts.items():
        if isinstance(spec, Mapping):
            dataset = spec.get("dataset", label)
            score_obj = spec.get("tf_scores", tf_scores.get(label))
        else:
            dataset = str(spec)
            score_obj = tf_scores.get(label)
        if score_obj is None:
            raise ValueError(f"No TF score table supplied for context {label!r}.")
        out[str(label)] = {"dataset": str(dataset), "tf_scores": score_obj}
    return out


def _context_edges(
    *,
    method: str,
    dataset: str,
    gene_scores: pd.DataFrame,
    tf_scores: pd.DataFrame,
    valid_targets: set[str] | None,
    tf_universe: set[str] | None,
    gene_column: str,
    dataset_column: str,
    lfdr_column: str,
    tf_column: str,
    tf_score_column: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    missing = {gene_column, dataset_column, lfdr_column}.difference(gene_scores.columns)
    if missing:
        raise ValueError(f"gene_scores is missing required columns: {sorted(missing)}")
    missing = {tf_column, tf_score_column}.difference(tf_scores.columns)
    if missing:
        raise ValueError(f"tf_scores for {method!r} is missing required columns: {sorted(missing)}")

    genes = gene_scores.loc[gene_scores[dataset_column].astype(str) == str(dataset), [gene_column, lfdr_column]].copy()
    genes = genes.rename(columns={gene_column: "target", lfdr_column: "lfdr"})
    genes["gene_score"] = 1.0 - pd.to_numeric(genes["lfdr"], errors="coerce")
    genes = genes[np.isfinite(genes["gene_score"]) & (genes["gene_score"] > 0)].copy()
    if valid_targets is not None:
        genes = genes[genes["target"].astype(str).isin(valid_targets)].copy()
    genes = genes.drop_duplicates("target", keep="first").reset_index(drop=True)
    if genes.empty:
        raise ValueError(f"No positive target-gene scores remain for dataset {dataset!r}.")

    tfs = tf_scores[[tf_column, tf_score_column]].copy()
    tfs = tfs.rename(columns={tf_column: "TF", tf_score_column: "priority_act"})
    tfs["priority_act"] = pd.to_numeric(tfs["priority_act"], errors="coerce")
    tfs = tfs[np.isfinite(tfs["priority_act"]) & (tfs["priority_act"] > 0)].copy()
    if tf_universe is not None:
        tfs = tfs[tfs["TF"].astype(str).isin(tf_universe)].copy()
    tfs = tfs.drop_duplicates("TF", keep="first").reset_index(drop=True)
    tf_total = float(tfs["priority_act"].sum())
    if not np.isfinite(tf_total) or tf_total <= 0:
        raise ValueError(f"No positive TF priority scores remain for context {method!r}.")
    tfs["tf_score"] = tfs["priority_act"] / tf_total

    edge_weight = np.outer(tfs["tf_score"].to_numpy(dtype=float), genes["gene_score"].to_numpy(dtype=float))
    rows = pd.DataFrame(
        {
            "TF": np.repeat(tfs["TF"].astype(str).to_numpy(), len(genes)),
            "target": np.tile(genes["target"].astype(str).to_numpy(), len(tfs)),
            "EdgeWeight": edge_weight.ravel(),
        }
    )
    rows = _legacy_edges.normalize_edges(rows, method, valid_targets or set(genes["target"].astype(str)))
    rows = _legacy_edges.break_ties(rows)
    return rows, tfs, genes


def infer(
    *,
    gene_scores: str | Path | pd.DataFrame,
    tf_scores: Mapping[str, str | Path | pd.DataFrame],
    valid_targets: str | Path | Iterable[str] | pd.DataFrame | None = None,
    contexts: Mapping[str, str | Mapping[str, Any]] | None = None,
    tf_universe: Iterable[str] | None = None,
    gene_column: str = "gene",
    dataset_column: str = "dataset",
    lfdr_column: str = "lfdr_bum",
    tf_column: str = "TF",
    tf_score_column: str = "priority_act",
    provenance: Mapping[str, Any] | None = None,
) -> GRNResult:
    """Infer ChronoPack TF-target edge scores from target and TF priorities."""
    gene_score_df = _read_table(gene_scores)
    valid = _load_valid_targets(valid_targets)
    tf_filter = {str(tf) for tf in tf_universe} if tf_universe is not None else None
    context_map = _normalize_contexts(tf_scores, contexts)

    edge_parts = []
    regulator_parts = []
    target_parts = []
    for method, spec in context_map.items():
        edge_df, tf_df, gene_df = _context_edges(
            method=method,
            dataset=str(spec["dataset"]),
            gene_scores=gene_score_df,
            tf_scores=_read_table(spec["tf_scores"]),
            valid_targets=valid,
            tf_universe=tf_filter,
            gene_column=gene_column,
            dataset_column=dataset_column,
            lfdr_column=lfdr_column,
            tf_column=tf_column,
            tf_score_column=tf_score_column,
        )
        edge_parts.append(edge_df)
        regulator_parts.append(tf_df.assign(method=method))
        target_parts.append(gene_df.assign(method=method, dataset=str(spec["dataset"])))

    legacy_edges = pd.concat(edge_parts, ignore_index=True)
    legacy_edges = legacy_edges.sort_values(["method", "EdgeWeight"], ascending=[True, False], kind="mergesort")
    legacy_edges = legacy_edges.drop_duplicates(["method", "TF", "target"], keep="first").reset_index(drop=True)
    legacy_edges["rank"] = legacy_edges.groupby("method").cumcount() + 1
    group_n = legacy_edges.groupby("method")["rank"].transform("max")
    legacy_edges["percentile"] = 1 - (legacy_edges["rank"] - 1) / group_n

    edges = legacy_edges.rename(columns={"TF": "regulator", "EdgeWeight": "score"})
    edges = edges[["method", "regulator", "target", "score", "rank", "percentile"]]
    return GRNResult(
        edges=edges,
        regulators=pd.concat(regulator_parts, ignore_index=True),
        target_scores=pd.concat(target_parts, ignore_index=True),
        parameters={
            "gene_column": gene_column,
            "dataset_column": dataset_column,
            "lfdr_column": lfdr_column,
            "tf_column": tf_column,
            "tf_score_column": tf_score_column,
            "tf_universe_size": None if tf_filter is None else len(tf_filter),
        },
        provenance=dict(provenance or {}),
    )


def score_targets(
    *,
    valid_targets_path: str | Path,
    gene_tables: Mapping[str, str | Path],
    output_dir: str | Path | None = None,
) -> TargetGeneScoreResult:
    """Run the packaged R/fdrtool target-gene score recalibration workflow."""
    source_r_resource("grn_inference_helpers.R")
    try:
        from rpy2.robjects import ListVector, r
    except ImportError as exc:
        raise ImportError("GRN target scoring requires the optional 'rpy2' dependency.") from exc

    r_gene_tables = ListVector({str(label): str(Path(path)) for label, path in gene_tables.items()})
    result = r["calculate_lfdr_recalibration"](str(Path(valid_targets_path)), r_gene_tables)

    if output_dir is None:
        with tempfile.TemporaryDirectory(prefix="chronopack-grn-targets-") as tmp:
            paths = r["write_lfdr_recalibration_outputs"](result, tmp)
            sidecars = {}
            scores = pd.read_csv(str(paths.rx2("scores")[0]))
            summary = pd.read_csv(str(paths.rx2("summary")[0]))
            params = pd.read_csv(str(paths.rx2("params")[0]))
    else:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        paths = r["write_lfdr_recalibration_outputs"](result, str(out))
        sidecars = {name: str(paths.rx2(name)[0]) for name in ("scores", "summary", "params")}
        scores = pd.read_csv(sidecars["scores"])
        summary = pd.read_csv(sidecars["summary"])
        params = pd.read_csv(sidecars["params"])

    return TargetGeneScoreResult(
        scores=scores,
        summary=summary,
        params=params,
        parameters={"valid_targets_path": str(valid_targets_path), "datasets": list(gene_tables.keys())},
        provenance={"backend": "fdrtool/BUM R helper", "r_resource": "grn_inference_helpers.R"},
        sidecar_paths=sidecars,
    )


def rank_regulators(
    *,
    human_tf_path: str | Path,
    ligand_receptor_path: str | Path,
    contexts: Mapping[str, Mapping[str, Any]],
    organism: str = "human",
    output_dir: str | Path | None = None,
    w_signed: float = 1.0,
    w_unsigned: float = 0.5,
    alpha: float = 0.85,
    max_steps: int = 80,
    min_steps: int = 10,
    tol_mass: float = 1e-12,
) -> TFRankingResult:
    """Rank context-specific regulators with the OmniPath signed diffusion workflow."""
    resources = load_omnipath_tf_resources(
        human_tf_path=human_tf_path,
        ligand_receptor_path=ligand_receptor_path,
        organism=organism,
    )
    graph_bundle = build_omnipath_tf_graph(resources)
    rankings = {}
    histories = {}
    for label, cfg in contexts.items():
        ranking, history = rank_omnipath_tfs(
            graph_bundle,
            resources.ligand_receptor,
            source=str(cfg["source"]),
            is_ligand=bool(cfg.get("is_ligand", False)),
            w_signed=w_signed,
            w_unsigned=w_unsigned,
            alpha=alpha,
            max_steps=max_steps,
            min_steps=min_steps,
            tol_mass=tol_mass,
        )
        rankings[str(label)] = ranking
        histories[str(label)] = history

    out = TFRankingResult(
        rankings=rankings,
        histories=histories,
        parameters={
            "organism": organism,
            "w_signed": w_signed,
            "w_unsigned": w_unsigned,
            "alpha": alpha,
            "max_steps": max_steps,
            "min_steps": min_steps,
            "tol_mass": tol_mass,
        },
        provenance={"human_tf_path": str(human_tf_path), "ligand_receptor_path": str(ligand_receptor_path)},
    )
    if output_dir is not None:
        out.write(output_dir)
    return out


def build_chronopack_ranked_edges(
    gene_scores_path: str | Path,
    k562_tf_scores_path: str | Path,
    thp1_tf_scores_path: str | Path,
    valid_targets_path: str | Path,
    trrust_path: str | Path | None = None,
    chip_encode_gt_path: str | Path | None = None,
) -> pd.DataFrame:
    """Compatibility wrapper returning the historical ranked-edge schema."""
    valid_targets = _legacy_edges.load_valid_targets(valid_targets_path)
    result = infer(
        gene_scores=gene_scores_path,
        tf_scores={
            "ChronoPack_K562": k562_tf_scores_path,
            "ChronoPack_THP1": thp1_tf_scores_path,
        },
        valid_targets=valid_targets,
        contexts={"ChronoPack_K562": "K562", "ChronoPack_THP1": "THP1"},
    )
    edges = result.to_legacy_edges()
    if trrust_path is not None and chip_encode_gt_path is not None:
        tf_universe = _legacy_edges.load_reference_tfs(
            valid_targets,
            trrust_path=trrust_path,
            chip_encode_gt_path=chip_encode_gt_path,
        )
        edges = edges[edges["TF"].isin(tf_universe)].copy()
    edges = edges.sort_values(["method", "EdgeWeight"], ascending=[True, False], kind="mergesort")
    edges = edges.drop_duplicates(["method", "TF", "target"], keep="first").reset_index(drop=True)
    edges["rank"] = edges.groupby("method").cumcount() + 1
    group_n = edges.groupby("method")["rank"].transform("max")
    edges["percentile"] = 1 - (edges["rank"] - 1) / group_n
    return edges


def write_chronopack_ranked_edges(ranked_edges: GRNResult | pd.DataFrame, path: str | Path) -> str:
    """Write a historical ranked-edge CSV."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(ranked_edges, GRNResult):
        df = ranked_edges.to_legacy_edges()
    else:
        df = ranked_edges.copy()
    df.to_csv(path, index=False)
    return str(path)


__all__ = [
    "GRNResult",
    "TargetGeneScoreResult",
    "TFRankingResult",
    "OmniPathTfResources",
    "OmniPathTfGraph",
    "build_omnipath_tf_graph",
    "load_omnipath_tf_resources",
    "rank_omnipath_tfs",
    "run_tf_diffusion_ranking",
    "infer",
    "rank_regulators",
    "score_targets",
    "build_chronopack_ranked_edges",
    "write_chronopack_ranked_edges",
]
