"""ChronoPack GRN inference helpers."""

from __future__ import annotations

import json
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from chronopack._backends.r._bridge import source_r_resource
from . import _edges
from ._omnipath import (
    OmniPathTfGraph,
    OmniPathTfResources,
    build_omnipath_tf_graph,
    cellchat_receptor_subunits,
    load_omnipath_tf_resources,
    rank_omnipath_tfs,
    ligand_receptor_seed,
    read_omnipath_snapshot,
    restrict_tf_graph,
    run_tf_diffusion_ranking,
    save_omnipath_snapshot,
)

# OmniPath graphs built by rank_regulators, reused across calls (building one takes ~30 s).
_TF_GRAPH_CACHE: dict = {}
_TF_GRAPH_CACHE_SIZE = 4


@dataclass
class TargetGeneScoreResult:
    """Target-gene scores produced by the fdrtool/BUM recalibration workflow."""

    scores: pd.DataFrame
    summary: pd.DataFrame
    params: pd.DataFrame
    parameters: dict[str, Any] = field(default_factory=dict)
    provenance: dict[str, Any] = field(default_factory=dict)

    def write(self, output_dir: str | Path, *, prefix: str = "target_scores") -> dict[str, str]:
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
        return _edges.load_valid_targets(valid_targets)
    if isinstance(valid_targets, pd.DataFrame):
        values = valid_targets.iloc[:, 0].astype(str)
    else:
        values = pd.Series(list(valid_targets), dtype=str)
    return {gene for gene in values if not gene.startswith("MT-") and not gene.startswith("ENSG")}


# Column conventions of the inputs to infer().
TARGET_GENE_COLUMN = "gene"
TARGET_LFDR_COLUMN = "lfdr_bum"  # beta-uniform-mixture local FDR from score_targets()
TARGET_TIEBREAK_COLUMN = "waldStat_1"  # orders equally scored targets (higher first)
TF_COLUMN = "TF"
TF_SCORE_COLUMN = "priority_act"  # net activating mass from rank_regulators()


def _context_edges(
    *,
    method: str,
    gene_scores: pd.DataFrame,
    tf_scores: pd.DataFrame,
    valid_targets: set[str] | None,
    tf_universe: set[str] | None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    missing = {TARGET_GENE_COLUMN, TARGET_LFDR_COLUMN}.difference(gene_scores.columns)
    if missing:
        raise ValueError(f"Target-gene scores are missing required columns: {sorted(missing)}")
    missing = {TF_COLUMN, TF_SCORE_COLUMN}.difference(tf_scores.columns)
    if missing:
        raise ValueError(f"TF scores are missing required columns: {sorted(missing)}")

    genes = gene_scores.copy()
    if TARGET_TIEBREAK_COLUMN in genes.columns:
        # Edges with equal weights keep the target order (stable sorts below),
        # so ordering targets by the tradeSeq Wald statistic breaks those ties.
        genes = genes.sort_values(TARGET_TIEBREAK_COLUMN, ascending=False, kind="mergesort", na_position="last")
    genes = genes[[TARGET_GENE_COLUMN, TARGET_LFDR_COLUMN]]
    genes = genes.rename(columns={TARGET_GENE_COLUMN: "target", TARGET_LFDR_COLUMN: "lfdr"})
    genes["gene_score"] = 1.0 - pd.to_numeric(genes["lfdr"], errors="coerce")
    genes = genes[np.isfinite(genes["gene_score"]) & (genes["gene_score"] > 0)].copy()
    if valid_targets is not None:
        genes = genes[genes["target"].astype(str).isin(valid_targets)].copy()
    genes = genes.drop_duplicates("target", keep="first").reset_index(drop=True)
    if genes.empty:
        raise ValueError(f"No positive target-gene scores remain for {method!r}.")

    tfs = tf_scores[[TF_COLUMN, TF_SCORE_COLUMN]].copy()
    tfs = tfs.rename(columns={TF_COLUMN: "TF", TF_SCORE_COLUMN: "priority_act"})
    tfs["priority_act"] = pd.to_numeric(tfs["priority_act"], errors="coerce")
    tfs = tfs[np.isfinite(tfs["priority_act"]) & (tfs["priority_act"] > 0)].copy()
    if tf_universe is not None:
        tfs = tfs[tfs["TF"].astype(str).isin(tf_universe)].copy()
    tfs = tfs.drop_duplicates("TF", keep="first").reset_index(drop=True)
    tf_total = float(tfs["priority_act"].sum())
    if not np.isfinite(tf_total) or tf_total <= 0:
        raise ValueError(f"No positive TF priority scores remain for {method!r}.")
    tfs["tf_score"] = tfs["priority_act"] / tf_total

    edge_weight = np.outer(tfs["tf_score"].to_numpy(dtype=float), genes["gene_score"].to_numpy(dtype=float))
    rows = pd.DataFrame(
        {
            "TF": np.repeat(tfs["TF"].astype(str).to_numpy(), len(genes)),
            "target": np.tile(genes["target"].astype(str).to_numpy(), len(tfs)),
            "EdgeWeight": edge_weight.ravel(),
        }
    )
    rows = _edges.normalize_edges(rows, method, valid_targets or set(genes["target"].astype(str)))
    rows = _edges.break_ties(rows)
    return rows, tfs, genes


def infer(
    targets: TargetGeneScoreResult | str | Path | pd.DataFrame,
    tfs: str | Path | pd.DataFrame,
    *,
    valid_targets: str | Path | Iterable[str] | pd.DataFrame | None = None,
    name: str = "ChronoPack",
    tf_universe: Iterable[str] | None = None,
) -> GRNResult:
    """Infer ranked TF -> target edges for one context.

    ``targets`` is the result of :func:`score_targets` (or a table / CSV with
    ``gene`` and ``lfdr_bum`` columns); ``tfs`` is the ranking from
    :func:`rank_regulators` (``TF`` and ``priority_act`` columns). Each edge is
    weighted by the TF's share of the total priority times the target's
    ``1 - lfdr``, restricted to ``valid_targets``, then ranked. Equal weights
    are ordered by TF priority, then by the target's tradeSeq ``waldStat_1``
    (higher first) when the target scores carry it. ``name`` fills the edges' ``method`` column.
    ``tf_universe`` optionally restricts the regulators (benchmark use).
    """
    gene_scores = targets.scores if isinstance(targets, TargetGeneScoreResult) else _read_table(targets)
    valid = _load_valid_targets(valid_targets)
    tf_filter = {str(tf) for tf in tf_universe} if tf_universe is not None else None
    edge_df, tf_df, gene_df = _context_edges(
        method=name,
        gene_scores=gene_scores,
        tf_scores=_read_table(tfs),
        valid_targets=valid,
        tf_universe=tf_filter,
    )

    edge_df = edge_df.sort_values(["method", "EdgeWeight"], ascending=[True, False], kind="mergesort")
    edge_df = edge_df.drop_duplicates(["method", "TF", "target"], keep="first").reset_index(drop=True)
    edge_df["rank"] = edge_df.groupby("method").cumcount() + 1
    group_n = edge_df.groupby("method")["rank"].transform("max")
    edge_df["percentile"] = 1 - (edge_df["rank"] - 1) / group_n

    edges = edge_df.rename(columns={"TF": "regulator", "EdgeWeight": "score"})
    edges = edges[["method", "regulator", "target", "score", "rank", "percentile"]]
    return GRNResult(
        edges=edges,
        regulators=tf_df.assign(method=name),
        target_scores=gene_df.assign(method=name),
        parameters={"name": name, "tf_universe_size": None if tf_filter is None else len(tf_filter)},
    )


def score_targets(association: str | Path, *, valid_targets: str | Path) -> TargetGeneScoreResult:
    """Score target genes for one context by local FDR (fdrtool / beta-uniform mixture, R).

    ``association`` is a tradeSeq association CSV with ``pvalue_1`` and ``lfdr``
    columns and gene names in an ``X`` column (as written by R's ``write.csv``)
    or otherwise the first column; ``valid_targets`` is a one-column CSV
    (no header) of allowed target genes. Returns per-gene scores (``lfdr_bum``
    is used by :func:`infer`) plus a summary and the fitted parameters.
    """
    source_r_resource("grn_inference_helpers.R")
    try:
        from rpy2.robjects import ListVector, r
    except ImportError as exc:
        raise ImportError("GRN target scoring requires the optional 'rpy2' dependency.") from exc

    # The R helper labels rows with a dataset name; one context per call, so it is dropped.
    result = r["calculate_lfdr_recalibration"](
        str(Path(valid_targets)), ListVector({"context": str(Path(association))})
    )
    with tempfile.TemporaryDirectory(prefix="chronopack-grn-targets-") as tmp:
        paths = r["write_lfdr_recalibration_outputs"](result, tmp)
        tables = {
            name: pd.read_csv(str(paths.rx2(name)[0])).drop(columns="dataset")
            for name in ("scores", "summary", "params")
        }

    return TargetGeneScoreResult(
        scores=tables["scores"],
        summary=tables["summary"],
        params=tables["params"],
        parameters={"association": str(association), "valid_targets": str(valid_targets)},
        provenance={"backend": "fdrtool/BUM R helper", "r_resource": "grn_inference_helpers.R"},
    )


def expressed_genes(adata, *, min_cells: int = 1, min_fraction: float = 0.0, layer: str | None = "counts") -> list[str]:
    """Genes detected (count > 0) in at least ``min_cells`` cells and ``min_fraction`` of the cells of ``adata``.

    Pass the cells of one cell type (e.g. ``adata[adata.obs["cell_type"] == "THP1"]``)
    and give the result to :func:`rank_regulators` as ``expressed_genes``. Uses
    ``layer`` (``counts`` by default, ``adata.X`` if absent).

    The default (detected in any cell) is deliberate: signaling proteins are
    rarely captured by single-cell RNA-seq (in the K562 and THP1 data RELA is
    detected in 4-7% of cells, TNFRSF1A in 4%, IKBKG in 0.3%), so a percentage
    cutoff removes genuine pathway members; at 5% the TNF receptor TNFRSF1A and
    RELA drop out.
    """
    X = adata.layers[layer] if layer is not None and layer in adata.layers else adata.X
    detected = np.asarray((X > 0).sum(axis=0)).ravel()
    keep = (detected >= min_cells) & (detected >= min_fraction * adata.n_obs)
    return [str(gene) for gene in np.asarray(adata.var_names)[keep]]


def _path_key(path) -> tuple:
    path = Path(path).resolve()
    files = sorted(path.glob("*")) if path.is_dir() else [path]
    return (str(path), tuple((f.name, f.stat().st_mtime_ns) for f in files if f.is_file()))


def _cached_tf_graph(*, human_tf_path, ligand_receptor_path, organism, omnipath_snapshot, **graph_kwargs):
    import datetime

    key = (_path_key(human_tf_path), _path_key(ligand_receptor_path), organism,
           _path_key(omnipath_snapshot) if omnipath_snapshot is not None else ("download", datetime.date.today()),
           tuple(sorted(graph_kwargs.items())))
    if key not in _TF_GRAPH_CACHE:
        resources = load_omnipath_tf_resources(
            human_tf_path=human_tf_path,
            ligand_receptor_path=ligand_receptor_path,
            organism=organism,
            omnipath_snapshot=omnipath_snapshot,
        )
        _TF_GRAPH_CACHE[key] = (resources, build_omnipath_tf_graph(resources, **graph_kwargs))
        while len(_TF_GRAPH_CACHE) > _TF_GRAPH_CACHE_SIZE:
            _TF_GRAPH_CACHE.pop(next(iter(_TF_GRAPH_CACHE)))
    return _TF_GRAPH_CACHE[key]


def rank_regulators(
    *,
    source: str,
    human_tf_path: str | Path,
    ligand_receptor_path: str | Path,
    is_ligand: bool = False,
    organism: str = "human",
    w_signed: float = 1.0,
    w_unsigned: float = 0.5,
    alpha: float = 0.85,
    max_steps: int = 80,
    min_steps: int = 10,
    tol_mass: float = 1e-12,
    omnipath_snapshot: str | Path | None = None,
    conflicting_signs: str = "consensus",
    edge_weight: str = "uniform",
    evidence: str = "all",
    second_messenger_bridges: bool = True,
    expressed_genes: Iterable[str] | None = None,
    expression_scope: str = "tfs",
    tied_signs: str = "split",
) -> pd.DataFrame:
    """Rank TFs downstream of one upstream signal by signed OmniPath diffusion.

    Mass starts at ``source`` (for a ligand, ``is_ligand=True``, at its receptors
    from the ligand-receptor table: each receptor gets an equal share, and a
    receptor complex such as ``IL6R_IL6ST`` is split over its subunits), spreads
    over OmniPath signaling edges, flipping sign across inhibitory edges, and is
    absorbed by TFs (human TFs with OmniPath transcriptional targets). Returns one
    row per TF sorted by ``priority_act`` (net activating mass absorbed). Raises
    ``ValueError`` if no mass reaches any TF (e.g. the receptor has no outgoing
    OmniPath edges, as for IL17RA, TNFRSF11A or IL12RB2 in the 2025-10-21 snapshot).

    Edges that OmniPath flags both stimulatory and inhibitory take its consensus
    sign (``conflicting_signs="consensus"``, default); ``"stimulation"`` treats
    them as activating, as ChronoPack did before 2026-10-21. Edges whose
    consensus is tied keep both flags; ``tied_signs="split"`` (default) treats
    them as of unknown sign (half the mass keeps its sign, half flips; see
    :func:`build_signed_row_normalized_transitions`), ``"stimulation"`` as
    activating (ChronoPack before 2026-10-22).

    ``second_messenger_bridges=True`` (default) adds the cAMP, IP3/Ca2+ and DAG
    steps that OmniPath's protein network lacks (adenylyl cyclase -> PKA/EPAC,
    phospholipase C -> calcineurin/CaMK/PKC/RasGRP); ``False`` gives the graph used
    before 2026-10-21. ``edge_weight="resources"`` / ``"references"`` splits each
    node's mass over its edges by the evidence for them, and
    ``evidence="physical"`` / ``"physical_curated"`` keeps only edges with
    physical-interaction support; both are off by default (see
    :func:`build_omnipath_tf_graph` and :func:`edge_evidence_weight`).

    ``expressed_genes`` (e.g. from :func:`expressed_genes` on the responding cell
    type) removes genes that are not expressed before the diffusion; ``None``
    (default) uses the whole human network. ``expression_scope="tfs"`` (default)
    removes only unexpressed TFs and ligand receptors, the ends of the network;
    ``"network"`` also removes unexpressed signaling proteins, which risks dropping
    genuine pathway members that single-cell RNA-seq rarely detects (e.g. IKBKG).
    The source receptor (``is_ligand=False``) must be expressed in either case.

    The OmniPath graph is built once per set of inputs and options and reused by
    later calls in the same session.

    OmniPath is downloaded (or taken from the ``omnipath`` client's cache), so the
    ranking can change slightly as the database is updated. Pass
    ``omnipath_snapshot``, a folder written by :func:`save_omnipath_snapshot`, to
    use fixed OmniPath tables instead.

    The ranking is a prior from a generic signaling network. It finds TFs one or
    two steps from the receptor well, and is weaker for TFs reached through long
    kinase cascades, for the choice among STATs, and against hub TFs that rank high
    for many signals. Why the defaults were chosen, how the ranking was checked and
    its known limitations: ``docs/user/grn.md`` in the ChronoPack repository.
    """
    resources, graph_bundle = _cached_tf_graph(
        human_tf_path=human_tf_path,
        ligand_receptor_path=ligand_receptor_path,
        organism=organism,
        omnipath_snapshot=omnipath_snapshot,
        conflicting_signs=conflicting_signs,
        evidence=evidence,
        second_messenger_bridges=second_messenger_bridges,
    )
    source = str(source)
    if expressed_genes is not None:
        full = graph_bundle
        receptors = {gene for receptor in resources.ligand_receptor["receptor"].dropna().astype(str)
                     for gene in cellchat_receptor_subunits(receptor)}
        expressed_set = {str(gene) for gene in expressed_genes}
        if not is_ligand:
            receptors.add(source)
        graph_bundle = restrict_tf_graph(full, expressed_set, scope=expression_scope, receptor_genes=receptors)
        if is_ligand:
            seed = ligand_receptor_seed(resources.ligand_receptor, source, nodes=full.graph)
            if seed and not ligand_receptor_seed(resources.ligand_receptor, source, nodes=graph_bundle.graph):
                raise ValueError(f"None of the receptors of ligand '{source}' ({sorted(seed)}) is in expressed_genes.")
        elif source in full.graph and source not in graph_bundle.graph:
            raise ValueError(f"Source '{source}' is not in expressed_genes.")
    ranking, _history = rank_omnipath_tfs(
        graph_bundle,
        resources.ligand_receptor,
        source=source,
        is_ligand=bool(is_ligand),
        w_signed=w_signed,
        w_unsigned=w_unsigned,
        alpha=alpha,
        max_steps=max_steps,
        min_steps=min_steps,
        tol_mass=tol_mass,
        edge_weight=edge_weight,
        tied_signs=tied_signs,
    )
    return ranking


def build_chronopack_ranked_edges(
    gene_scores_path: str | Path,
    k562_tf_scores_path: str | Path,
    thp1_tf_scores_path: str | Path,
    valid_targets_path: str | Path,
    trrust_path: str | Path | None = None,
    chip_encode_gt_path: str | Path | None = None,
) -> pd.DataFrame:
    """Compatibility wrapper returning the historical ranked-edge schema.

    ``gene_scores_path`` is the combined K562/THP1 target-score table with a
    ``dataset`` column, as written by the former multi-context workflow.
    """
    valid_targets = _edges.load_valid_targets(valid_targets_path)
    gene_scores = _read_table(gene_scores_path)
    parts = []
    for dataset, tf_path in (("K562", k562_tf_scores_path), ("THP1", thp1_tf_scores_path)):
        result = infer(
            gene_scores[gene_scores["dataset"].astype(str) == dataset],
            tf_path,
            valid_targets=valid_targets,
            name=f"ChronoPack_{dataset}",
        )
        parts.append(result.to_legacy_edges())
    edges = pd.concat(parts, ignore_index=True)
    if trrust_path is not None and chip_encode_gt_path is not None:
        tf_universe = _edges.load_reference_tfs(
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
    "OmniPathTfResources",
    "OmniPathTfGraph",
    "build_omnipath_tf_graph",
    "load_omnipath_tf_resources",
    "rank_omnipath_tfs",
    "read_omnipath_snapshot",
    "save_omnipath_snapshot",
    "run_tf_diffusion_ranking",
    "infer",
    "expressed_genes",
    "rank_regulators",
    "score_targets",
    "build_chronopack_ranked_edges",
    "write_chronopack_ranked_edges",
]
