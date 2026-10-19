"""OmniPath TF prioritization helpers."""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd


DEFAULT_DROP_SOURCES = {
    "NetPath",
    "Wang",
    "InnateDB",
    "KEGG-MEDICUS",
    "SPIKE",
    "SPIKE_LC",
    "Lit-BM-17",
    "HINT",
    "HPRD",
    "Macrophage",
}


@dataclass(frozen=True)
class OmniPathTfResources:
    """Raw resources used to construct an OmniPath TF ranking graph."""

    signaling: pd.DataFrame
    transcriptional: pd.DataFrame
    ligand_receptor: pd.DataFrame
    human_tfs: set[str]


@dataclass(frozen=True)
class OmniPathTfGraph:
    """Split signed graph plus the TF universe used for diffusion ranking."""

    graph: nx.DiGraph
    tfs: list[str]
    interactions: pd.DataFrame
    split_interactions: pd.DataFrame
    merged_interactions: pd.DataFrame


SNAPSHOT_FILES = {"signaling": "signaling.csv.gz", "transcriptional": "transcriptional.csv.gz"}
SNAPSHOT_METADATA = "metadata.json"


def _download_omnipath_tables(organism: str = "human") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Signed directed signaling and transcriptional interactions from the OmniPath web service (or the client's cache)."""

    import omnipath as op

    signaling = op.interactions.OmniPath.get(
        genesymbols=1,
        directed=1,
        signed=1,
        organisms=organism,
    )
    transcriptional = op.interactions.Transcriptional.get(
        genesymbols=1,
        organisms=organism,
    )
    return signaling, transcriptional


def save_omnipath_snapshot(
    path: str | Path,
    *,
    organism: str = "human",
    retrieved: str | None = None,
    note: str | None = None,
) -> Path:
    """Save the OmniPath tables used by :func:`rank_regulators` to a snapshot folder.

    OmniPath is updated over time, so TF rankings computed from a live download change slightly
    between runs. The snapshot (``signaling.csv.gz``, ``transcriptional.csv.gz`` and
    ``metadata.json``) can be passed as ``omnipath_snapshot`` to :func:`rank_regulators` to
    reproduce a ranking exactly. The tables come from the ``omnipath`` client, which may return
    its local cache rather than a fresh download; ``retrieved`` (default: today) and ``note``
    record where they came from.
    """
    import datetime
    import json
    import omnipath as op

    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    tables = dict(zip(("signaling", "transcriptional"), _download_omnipath_tables(organism)))
    for name, table in tables.items():
        table.to_csv(path / SNAPSHOT_FILES[name], index=False)
    metadata = {
        "retrieved": retrieved or datetime.date.today().isoformat(),
        "saved": datetime.date.today().isoformat(),
        "organism": organism,
        "omnipath_client_version": getattr(op, "__version__", None),
        "queries": {
            "signaling": "omnipath.interactions.OmniPath.get(genesymbols=1, directed=1, signed=1, organisms=organism)",
            "transcriptional": "omnipath.interactions.Transcriptional.get(genesymbols=1, organisms=organism)",
        },
        "n_rows": {name: int(len(table)) for name, table in tables.items()},
        "note": note,
    }
    (path / SNAPSHOT_METADATA).write_text(json.dumps(metadata, indent=2) + "\n")
    return path


def read_omnipath_snapshot(path: str | Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Signaling and transcriptional tables of a snapshot written by :func:`save_omnipath_snapshot`."""

    path = Path(path)
    missing = [f for f in SNAPSHOT_FILES.values() if not (path / f).exists()]
    if missing:
        raise FileNotFoundError(f"OmniPath snapshot {path} is missing {missing}.")
    # Only empty fields are missing values, so gene symbols such as "NA" stay strings.
    read = lambda f: pd.read_csv(path / f, keep_default_na=False, na_values=[""], low_memory=False)
    return read(SNAPSHOT_FILES["signaling"]), read(SNAPSHOT_FILES["transcriptional"])


def load_omnipath_tf_resources(
    *,
    human_tf_path: str | Path,
    ligand_receptor_path: str | Path,
    organism: str = "human",
    omnipath_snapshot: str | Path | None = None,
) -> OmniPathTfResources:
    """Load OmniPath signaling/TF tables and local TF/LR reference files.

    The OmniPath tables are downloaded (or taken from the ``omnipath`` client's cache) unless
    ``omnipath_snapshot`` names a folder written by :func:`save_omnipath_snapshot`.
    """

    if omnipath_snapshot is None:
        signaling, transcriptional = _download_omnipath_tables(organism)
    else:
        signaling, transcriptional = read_omnipath_snapshot(omnipath_snapshot)
    ligand_receptor = pd.read_csv(ligand_receptor_path, index_col=0)
    human_tfs = set(pd.read_csv(human_tf_path, header=None).iloc[:, 0].values)
    return OmniPathTfResources(
        signaling=signaling,
        transcriptional=transcriptional,
        ligand_receptor=ligand_receptor,
        human_tfs=human_tfs,
    )


def unique_sources(df: pd.DataFrame, col: str = "sources", sep: str = ";") -> set[str]:
    """Return unique source database names encoded in OmniPath's sources column."""

    src = (
        df[col]
        .dropna()
        .astype(str)
        .str.split(sep)
        .explode()
        .str.strip()
    )
    return set(src[src.ne("")].unique())


def build_overlap_with_complexes(
    human_tfs: set[str],
    omnipath_tfs: set[str],
    *,
    sep: str = "_",
) -> dict[str, set[str]]:
    """Build the TF universe, retaining OmniPath complexes whose parts are known TFs."""

    hs = {s.strip().upper() for s in human_tfs if isinstance(s, str)}
    op_tfs = {s.strip().upper() for s in omnipath_tfs if isinstance(s, str)}

    singles_overlap = hs & op_tfs
    complex_strings = {s for s in op_tfs if sep in s}
    complex_ok = set()
    complex_components = set()

    for complex_name in complex_strings:
        parts = [p for p in complex_name.split(sep) if p]
        if all(p in hs for p in parts):
            complex_ok.add(complex_name)
            complex_components.update(parts)

    return {
        "overlap_tfs": singles_overlap | complex_ok,
        "singles_overlap": singles_overlap,
        "complexes_included": complex_ok,
        "components_from_included_complexes": complex_components,
    }


def filter_signaling_interactions(
    signaling: pd.DataFrame,
    *,
    drop_sources: set[str] = DEFAULT_DROP_SOURCES,
) -> pd.DataFrame:
    """Apply the source and directed-edge filters used by the original notebook."""

    keep_sources = unique_sources(signaling, "sources") - set(drop_sources)
    source_pattern = "|".join(keep_sources)

    interactions = signaling[signaling["is_directed"] == True].copy()
    interactions = interactions[interactions["source"] != interactions["target"]]
    interactions = interactions[interactions["sources"].str.contains(source_pattern)]
    interactions = interactions[interactions["n_sources"] >= 1]
    return interactions


def _is_valid_symbol(symbol: str) -> bool:
    return bool(re.fullmatch(r"[A-Z0-9]+", symbol))


def explode_symbol(
    symbol: str,
    *,
    sep: str = "_",
    known_symbols: set[str] | None = None,
    only_two_parts: bool = True,
) -> list[str]:
    """Split simple compound symbols such as NFKB1_RELA into component symbols."""

    if not isinstance(symbol, str) or sep not in symbol:
        return [symbol]
    parts = symbol.split(sep)
    if only_two_parts and len(parts) != 2:
        return [symbol]
    left, right = parts[0].strip(), parts[1].strip()
    if not left or not right:
        return [symbol]
    if not (_is_valid_symbol(left) and _is_valid_symbol(right)):
        return [symbol]
    if known_symbols is not None and ((left not in known_symbols) or (right not in known_symbols)):
        return [symbol]
    return [left, right]


def explode_edges(
    interactions: pd.DataFrame,
    *,
    src_col: str = "source_genesymbol",
    tgt_col: str = "target_genesymbol",
    known_symbols: set[str] | None = None,
    only_two_parts: bool = True,
) -> pd.DataFrame:
    """Expand compound source/target symbols and keep the original edge metadata."""

    missing = [col for col in [src_col, tgt_col] if col not in interactions.columns]
    if missing:
        raise ValueError(f"interactions missing columns: {missing}")

    rows = []
    for _, row in interactions.iterrows():
        sources = explode_symbol(row[src_col], known_symbols=known_symbols, only_two_parts=only_two_parts)
        targets = explode_symbol(row[tgt_col], known_symbols=known_symbols, only_two_parts=only_two_parts)
        for source in sources:
            for target in targets:
                if source == target:
                    continue
                expanded = row.to_dict()
                expanded["u"] = str(source)
                expanded["v"] = str(target)
                rows.append(expanded)

    expanded = pd.DataFrame(rows)
    expanded = expanded.dropna(subset=["u", "v"])
    return expanded[expanded["u"] != expanded["v"]].copy()


def _merge_sources(series: pd.Series) -> str | float:
    values = []
    for value in series.dropna().tolist():
        if isinstance(value, (list, tuple, set)):
            values.extend(list(value))
        elif isinstance(value, str):
            if "|" in value:
                values.extend([token for token in value.split("|") if token])
            elif ";" in value:
                values.extend([token for token in value.split(";") if token])
            else:
                values.append(value)
        else:
            values.append(str(value))
    values = sorted({str(value).strip() for value in values if value and str(value).strip()})
    return "|".join(values) if values else np.nan


def _merge_refs(series: pd.Series) -> str | float:
    values = []
    for value in series.dropna().tolist():
        if isinstance(value, str):
            if ";" in value:
                values.extend([token for token in value.split(";") if token])
            elif "|" in value:
                values.extend([token for token in value.split("|") if token])
            else:
                values.append(value)
        else:
            values.append(str(value))
    values = sorted({str(value).strip() for value in values if value and str(value).strip()})
    return ";".join(values) if values else np.nan


def build_original_edge_sign_map(
    interactions: pd.DataFrame,
    tfs: set[str],
    *,
    src_col: str = "source_genesymbol",
    tgt_col: str = "target_genesymbol",
) -> dict[tuple[str, str], tuple[bool, bool, bool]]:
    """Map original unambiguous source-to-TF signs before compound splitting."""

    df = interactions.copy()
    for col in ["is_stimulation", "is_inhibition"]:
        df[col] = df[col].fillna(False).astype(bool) if col in df.columns else False

    df["u0"] = df[src_col].astype(str)
    df["v0"] = df[tgt_col].astype(str)
    df = df[df["v0"].isin(set(map(str, tfs)))].copy()
    if df.empty:
        return {}

    merged = (
        df.groupby(["u0", "v0"], as_index=False)
        .agg({"is_stimulation": "any", "is_inhibition": "any"})
    )
    unambiguous = merged[merged["is_stimulation"] ^ merged["is_inhibition"]].copy()
    return {
        (row["u0"], row["v0"]): (bool(row["is_stimulation"]), bool(row["is_inhibition"]), True)
        for _, row in unambiguous.iterrows()
    }


def merge_duplicate_edges(
    expanded: pd.DataFrame,
    *,
    prefer_sign_map: dict[tuple[str, str], tuple[bool, bool, bool]] | None = None,
) -> pd.DataFrame:
    """Collapse duplicate expanded edges, preserving signs and selected metadata."""

    df = expanded.copy()
    for col in ["is_stimulation", "is_inhibition"]:
        df[col] = df[col].fillna(False).astype(bool) if col in df.columns else False

    agg_dict = {"is_stimulation": "any", "is_inhibition": "any"}
    if "curation_effort" in df.columns:
        agg_dict["curation_effort"] = "max"
    if "sources" in df.columns:
        agg_dict["sources"] = _merge_sources
    if "references" in df.columns:
        agg_dict["references"] = _merge_refs

    merged = df.groupby(["u", "v"], as_index=False).agg(agg_dict)

    if prefer_sign_map:
        preferences = pd.DataFrame(
            [
                (u, v, stimulation, inhibition)
                for (u, v), (stimulation, inhibition, ok) in prefer_sign_map.items()
                if ok
            ],
            columns=["u", "v", "_pref_stim", "_pref_inhib"],
        )
        if not preferences.empty:
            merged = merged.merge(preferences, on=["u", "v"], how="left")
            mask = merged["_pref_stim"].notna() & merged["_pref_inhib"].notna()
            merged.loc[mask, "is_stimulation"] = merged.loc[mask, "_pref_stim"].astype(bool)
            merged.loc[mask, "is_inhibition"] = merged.loc[mask, "_pref_inhib"].astype(bool)
            merged = merged.drop(columns=["_pref_stim", "_pref_inhib"])

    merged["unsigned_for_diffusion"] = merged["is_stimulation"] & merged["is_inhibition"]
    return merged


def build_digraph_from_merged(merged: pd.DataFrame) -> nx.DiGraph:
    """Create a NetworkX DiGraph from merged edge metadata."""

    graph = nx.DiGraph()
    for _, row in merged.iterrows():
        source, target = row["u"], row["v"]
        attrs = {key: row[key] for key in merged.columns if key not in ("u", "v")}
        graph.add_edge(source, target, **attrs)
    return graph


def explode_tf_set(tfs: set[str], *, known_symbols: set[str] | None = None) -> set[str]:
    """Apply the same compound-splitting rules to a TF set."""

    out = set()
    for tf in tfs:
        out.update(explode_symbol(tf, known_symbols=known_symbols, only_two_parts=True))
    return out


def build_omnipath_tf_graph(
    resources: OmniPathTfResources,
    *,
    known_symbols: set[str] | None = None,
    drop_sources: set[str] = DEFAULT_DROP_SOURCES,
) -> OmniPathTfGraph:
    """Build the split signed OmniPath graph and corresponding TF list."""

    omnipath_tfs = set(resources.transcriptional["source_genesymbol"].unique())
    overlap = build_overlap_with_complexes(resources.human_tfs, omnipath_tfs)
    overlap_tfs = overlap["overlap_tfs"]
    interactions = filter_signaling_interactions(resources.signaling, drop_sources=drop_sources)

    base_tfs = set()
    for tf in overlap_tfs:
        parts = explode_symbol(tf, known_symbols=known_symbols, only_two_parts=True)
        if len(parts) == 1 and str(parts[0]) == str(tf):
            base_tfs.add(str(tf))

    prefer_sign_map = build_original_edge_sign_map(interactions, base_tfs)
    expanded = explode_edges(
        interactions,
        src_col="source_genesymbol",
        tgt_col="target_genesymbol",
        known_symbols=known_symbols,
        only_two_parts=True,
    )
    merged = merge_duplicate_edges(expanded, prefer_sign_map=prefer_sign_map)
    graph = build_digraph_from_merged(merged)
    tfs = sorted(explode_tf_set(overlap_tfs, known_symbols=known_symbols) & set(graph.nodes()))

    return OmniPathTfGraph(
        graph=graph,
        tfs=tfs,
        interactions=interactions,
        split_interactions=expanded,
        merged_interactions=merged,
    )


def edge_sign_from_attrs(attrs: dict) -> int:
    """Return +1/-1/0 from OmniPath stimulation/inhibition edge attributes."""

    if not attrs:
        return 0
    if attrs.get("is_stimulation"):
        return +1
    if attrs.get("is_inhibition"):
        return -1
    return 0


def build_lr_graph(
    lr_df: pd.DataFrame,
    ligand_col: str = "ligand",
    receptor_col: str = "receptor",
    normalize: bool = True,
) -> nx.DiGraph:
    """Construct a ligand-receptor graph from a CellChat-style table."""

    graph = nx.DiGraph()
    if lr_df is None or len(lr_df) == 0:
        return graph
    df = lr_df[[ligand_col, receptor_col]].dropna().copy()
    if normalize:
        df[ligand_col] = df[ligand_col].astype(str).str.strip()
        df[receptor_col] = df[receptor_col].astype(str).str.strip()
    for ligand, receptor in df.itertuples(index=False, name=None):
        if ligand and receptor:
            graph.add_edge(ligand, receptor)
    return graph


def build_signed_row_normalized_transitions(
    graph: nx.DiGraph,
    *,
    tfs: set[str] | list[str],
    is_ligand: bool = False,
    lr_df: pd.DataFrame | None = None,
    ligand_col: str = "ligand",
    receptor_col: str = "receptor",
    restrict_lr_to_nodes_in_G: bool = True,
    w_signed: float = 1.0,
    w_unsigned: float = 0.5,
    treat_unsigned_as_positive: bool = True,
) -> tuple[dict, dict, dict, dict]:
    """Build positive/negative transition maps for signed absorbing diffusion."""

    tfs = set(tfs)
    lr_graph = build_lr_graph(lr_df, ligand_col, receptor_col) if is_ligand else None
    lr_receptors = {}
    if is_ligand and lr_graph is not None:
        for ligand in lr_graph.nodes:
            receptors = list(lr_graph.successors(ligand))
            if restrict_lr_to_nodes_in_G:
                receptors = [receptor for receptor in receptors if receptor in graph]
            if receptors:
                lr_receptors[ligand] = receptors

    t_pos = {}
    t_neg = {}
    outdeg_eff = {}

    for source in graph.nodes:
        if source in tfs:
            t_pos[source] = []
            t_neg[source] = []
            outdeg_eff[source] = 0
            continue

        pos_edges = []
        neg_edges = []
        for target in graph.successors(source):
            sign = edge_sign_from_attrs(graph.get_edge_data(source, target, {}) or {})
            if sign == -1:
                neg_edges.append((target, w_signed))
            elif sign == +1:
                pos_edges.append((target, w_signed))
            elif treat_unsigned_as_positive:
                pos_edges.append((target, w_unsigned))
            else:
                pos_edges.append((target, w_unsigned))

        denom = sum(weight for _, weight in pos_edges) + sum(weight for _, weight in neg_edges)
        if denom <= 0:
            t_pos[source] = []
            t_neg[source] = []
            outdeg_eff[source] = 0
        else:
            t_pos[source] = [(target, weight / denom) for target, weight in pos_edges]
            t_neg[source] = [(target, weight / denom) for target, weight in neg_edges]
            outdeg_eff[source] = len({target for target, _ in pos_edges + neg_edges})

    return t_pos, t_neg, outdeg_eff, lr_receptors


def signed_absorbing_diffusion_impulse(
    nodes: list[str],
    t_pos: dict,
    t_neg: dict,
    *,
    tfs: set[str] | list[str],
    seed: dict[str, float],
    alpha: float = 0.85,
    max_steps: int = 200,
    min_steps: int = 10,
    tol_mass: float = 1e-12,
) -> tuple[dict, dict, np.ndarray, np.ndarray, pd.DataFrame]:
    """Run signed impulse diffusion with TF nodes as absorbing states."""

    idx = {node: i for i, node in enumerate(nodes)}
    x_pos = np.zeros(len(nodes), dtype=float)
    x_neg = np.zeros(len(nodes), dtype=float)

    for node, value in seed.items():
        if node in idx:
            x_pos[idx[node]] += float(value)

    total = x_pos.sum()
    if total <= 0:
        raise ValueError("Seed has zero mass; seed nodes are not in graph.")
    x_pos /= total

    tfs = list(tfs)
    absorbed_pos = defaultdict(float)
    absorbed_neg = defaultdict(float)
    hist = []

    for step in range(max_steps):
        y_pos = np.zeros_like(x_pos)
        y_neg = np.zeros_like(x_neg)

        for source in nodes:
            source_idx = idx[source]
            pos_mass = x_pos[source_idx]
            neg_mass = x_neg[source_idx]
            if pos_mass == 0 and neg_mass == 0:
                continue

            for target, prob in t_pos.get(source, []):
                target_idx = idx[target]
                y_pos[target_idx] += pos_mass * prob
                y_neg[target_idx] += neg_mass * prob

            for target, prob in t_neg.get(source, []):
                target_idx = idx[target]
                y_neg[target_idx] += pos_mass * prob
                y_pos[target_idx] += neg_mass * prob

        absorbed_step = 0.0
        for tf in tfs:
            tf_idx = idx.get(tf)
            if tf_idx is None:
                continue
            ap = y_pos[tf_idx]
            an = y_neg[tf_idx]
            if ap:
                absorbed_pos[tf] += float(ap)
            if an:
                absorbed_neg[tf] += float(an)
            absorbed_step += float(ap + an)
            y_pos[tf_idx] = 0.0
            y_neg[tf_idx] = 0.0

        x_pos = alpha * y_pos
        x_neg = alpha * y_neg
        transient_mass = float(x_pos.sum() + x_neg.sum())
        absorbed_total = float(sum(absorbed_pos.values()) + sum(absorbed_neg.values()))
        hist.append(
            {
                "step": step,
                "absorbed_step": absorbed_step,
                "absorbed_total": absorbed_total,
                "transient_mass": transient_mass,
            }
        )

        if step >= min_steps and transient_mass < tol_mass:
            break

    return dict(absorbed_pos), dict(absorbed_neg), x_pos, x_neg, pd.DataFrame(hist)


def run_tf_diffusion_ranking(
    graph: nx.DiGraph,
    *,
    source: str,
    tfs: set[str] | list[str],
    is_ligand: bool = False,
    lr_df: pd.DataFrame | None = None,
    ligand_col: str = "ligand",
    receptor_col: str = "receptor",
    restrict_lr_to_nodes_in_G: bool = True,
    w_signed: float = 1.0,
    w_unsigned: float = 0.5,
    alpha: float = 0.85,
    max_steps: int = 80,
    min_steps: int = 10,
    tol_mass: float = 1e-12,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rank TFs by signed absorbing diffusion from a source node or ligand."""

    tfs = list(tfs)
    t_pos, t_neg, _, lr_map = build_signed_row_normalized_transitions(
        graph,
        tfs=tfs,
        is_ligand=is_ligand,
        lr_df=lr_df,
        ligand_col=ligand_col,
        receptor_col=receptor_col,
        restrict_lr_to_nodes_in_G=restrict_lr_to_nodes_in_G,
        w_signed=w_signed,
        w_unsigned=w_unsigned,
        treat_unsigned_as_positive=True,
    )

    nodes = list(graph.nodes)
    seed = {}
    if is_ligand:
        receptors = lr_map.get(source, [])
        if len(receptors) == 0:
            raise ValueError(f"No ligand-receptor targets found for ligand '{source}'.")
        weight = 1.0 / len(receptors)
        for receptor in receptors:
            seed[receptor] = seed.get(receptor, 0.0) + weight
    else:
        if source not in graph:
            raise ValueError(f"Source '{source}' not in graph.")
        seed[source] = 1.0

    absorbed_pos, absorbed_neg, _, _, hist = signed_absorbing_diffusion_impulse(
        nodes,
        t_pos,
        t_neg,
        tfs=tfs,
        seed=seed,
        alpha=alpha,
        max_steps=max_steps,
        min_steps=min_steps,
        tol_mass=tol_mass,
    )

    rows = []
    for tf in tfs:
        abs_pos = absorbed_pos.get(tf, 0.0)
        abs_neg = absorbed_neg.get(tf, 0.0)
        support = abs_pos + abs_neg
        net = abs_pos - abs_neg
        score = (net / support) if support > 0 else 0.0
        rows.append(
            {
                "TF": tf,
                "abs_pos": abs_pos,
                "abs_neg": abs_neg,
                "support": support,
                "net": net,
                "score": score,
            }
        )

    ranking = pd.DataFrame(rows)
    ranking["priority_act"] = ranking["net"].clip(lower=0)
    ranking["priority_abs"] = ranking["support"] * ranking["score"].abs()
    ranking = ranking.sort_values(["priority_act", "support"], ascending=[False, False]).reset_index(drop=True)
    return ranking, hist


def rank_omnipath_tfs(
    graph_bundle: OmniPathTfGraph,
    ligand_receptor: pd.DataFrame,
    *,
    source: str,
    is_ligand: bool,
    ligand_col: str = "ligand",
    receptor_col: str = "receptor",
    w_signed: float = 1.0,
    w_unsigned: float = 0.5,
    alpha: float = 0.85,
    max_steps: int = 80,
    min_steps: int = 10,
    tol_mass: float = 1e-12,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Convenience wrapper to rank TFs for a source in an OmniPath TF graph."""

    return run_tf_diffusion_ranking(
        graph_bundle.graph,
        source=source,
        tfs=graph_bundle.tfs,
        is_ligand=is_ligand,
        lr_df=ligand_receptor,
        ligand_col=ligand_col,
        receptor_col=receptor_col,
        restrict_lr_to_nodes_in_G=True,
        w_signed=w_signed,
        w_unsigned=w_unsigned,
        alpha=alpha,
        max_steps=max_steps,
        min_steps=min_steps,
        tol_mass=tol_mass,
    )
