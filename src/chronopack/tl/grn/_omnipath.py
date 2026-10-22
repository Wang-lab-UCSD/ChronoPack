"""OmniPath TF prioritization helpers."""

from __future__ import annotations

import re
import warnings
from dataclasses import dataclass
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from scipy import sparse


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


# OmniPath resources by the kind of evidence they hold (base names; "X_Y" tokens are X's data processed by Y).
PPI_RESOURCES = {"BioGRID", "DIP", "DOMINO", "DeathDomain", "ELM", "HINT", "HPRD", "HuRI", "InnateDB", "IntAct",
                 "LMPID", "Lit-BM-17", "MPPI", "PDZBase"}
ENZYME_SUBSTRATE_RESOURCES = {"DEPOD", "HPRD-phos", "KEA", "Li2012", "MIMP", "PhosphoNetworks", "PhosphoPoint",
                              "PhosphoSite", "ProtMapper", "dbPTM", "iPTMnet", "phosphoELM"}
LIGAND_RECEPTOR_RESOURCES = {"Baccin2019", "CellCall", "CellChatDB", "CellPhoneDB", "CellTalkDB", "Cellinker", "DLRP",
                             "EMBRACE", "Fantom5", "Guide2Pharma", "HPMR", "ICELLNET", "Kirouac2010", "LRdb",
                             "Ramilowski2015", "UniProt", "Wojtowicz2020", "connectomeDB2020", "iTALK", "scConnect",
                             "talklr"}
TEXT_MINED_OR_PREDICTED_RESOURCES = {"REACH", "Sparser", "RLIMS-P", "NetworKIN", "Kinexus", "STRING"}
PHYSICAL_RESOURCES = PPI_RESOURCES | ENZYME_SUBSTRATE_RESOURCES | LIGAND_RECEPTOR_RESOURCES

# Second-messenger steps that OmniPath's protein-protein signaling tables cannot hold (cAMP, IP3/Ca2+, DAG):
# adenylyl cyclases -> cAMP -> PKA catalytic subunits and EPAC; phospholipase C -> IP3/Ca2+ -> calcineurin and
# CaMKII/IV, and -> DAG -> classical/novel PKC isoforms and RasGRPs. Added as activating edges when absent.
_ADCY = [f"ADCY{i}" for i in range(1, 11)]
_PLC = ["PLCB1", "PLCB2", "PLCB3", "PLCB4", "PLCG1", "PLCG2"]
SECOND_MESSENGER_BRIDGES = (
    [(a, t, "cAMP") for a in _ADCY for t in ("PRKACA", "PRKACB", "RAPGEF3", "RAPGEF4")]
    + [(p, t, "IP3/Ca2+") for p in _PLC for t in ("PPP3CA", "PPP3CB", "PPP3CC", "CAMK2A", "CAMK2B", "CAMK2D",
                                                  "CAMK2G", "CAMK4")]
    + [(p, t, "DAG") for p in _PLC for t in ("PRKCA", "PRKCB", "PRKCG", "PRKCD", "PRKCE", "PRKCH", "PRKCQ",
                                             "RASGRP1", "RASGRP3")]
)

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


def resolve_conflicting_signs(interactions: pd.DataFrame) -> pd.DataFrame:
    """Use OmniPath's consensus sign for edges flagged both stimulatory and inhibitory.

    Some resources call an edge activating and others inhibiting; OmniPath then sets
    both ``is_stimulation`` and ``is_inhibition`` and records the majority in
    ``consensus_stimulation`` / ``consensus_inhibition``. Where that consensus is one
    sign, the edge gets that sign; where it is tied (or the consensus columns are
    absent) the edge keeps both flags.
    """

    df = interactions.copy()
    if not {"consensus_stimulation", "consensus_inhibition"}.issubset(df.columns):
        return df
    stim = df["is_stimulation"].fillna(False).astype(bool)
    inhib = df["is_inhibition"].fillna(False).astype(bool)
    c_stim = df["consensus_stimulation"].fillna(False).astype(bool)
    c_inhib = df["consensus_inhibition"].fillna(False).astype(bool)
    resolve = stim & inhib & (c_stim ^ c_inhib)
    df.loc[resolve, "is_stimulation"] = c_stim[resolve]
    df.loc[resolve, "is_inhibition"] = c_inhib[resolve]
    return df


def _source_tokens(sources) -> list[str]:
    if not isinstance(sources, str):
        return []
    return [token.strip() for token in re.split(r"[;|]", sources) if token.strip()]


def evidence_class(resource: str) -> str:
    """``"physical"``, ``"text-mined/predicted"`` or ``"activity flow"`` for an OmniPath resource token.

    Physical: protein-protein interaction, enzyme-substrate (incl. sites mapped by
    ProtMapper, MIMP or KEA) and ligand-receptor resources. Activity flow: pathway
    resources (SIGNOR, SPIKE, SignaLink, ...) that record that A acts on B, directly
    or not.
    """

    base, _, via = resource.partition("_")
    if base in TEXT_MINED_OR_PREDICTED_RESOURCES:
        return "text-mined/predicted"
    if base in PHYSICAL_RESOURCES or via in {"ProtMapper", "MIMP", "KEA"} or via in LIGAND_RECEPTOR_RESOURCES:
        return "physical"
    return "activity flow"


# Tokens naming an aggregator whose records also appear as "<origin>_<aggregator>" tokens (e.g. REACH_ProtMapper).
AGGREGATOR_TOKENS = {"ProtMapper", "MIMP", "KEA"}


def has_physical_evidence(sources, *, include_text_mined: bool = True) -> bool:
    """Whether an OmniPath ``sources`` string includes a physical-interaction resource.

    A bare aggregator token (``ProtMapper``, ``MIMP``, ``KEA``) defers to the
    ``<origin>_<aggregator>`` tokens that come with it, so text-mined sites
    (``REACH_ProtMapper``) are not counted as curated through ``ProtMapper``.
    """

    tokens = _source_tokens(sources)
    specific = [token for token in tokens if token not in AGGREGATOR_TOKENS]
    classes = {evidence_class(token) for token in specific}
    for aggregator in AGGREGATOR_TOKENS.intersection(tokens):
        if not any(token.endswith("_" + aggregator) for token in specific):
            classes.add("physical")
    return "physical" in classes or (include_text_mined and "text-mined/predicted" in classes)


def edge_evidence_weight(attrs: dict, mode: str = "uniform") -> float:
    """Evidence weight of a graph edge: 1 (``"uniform"``), the number of distinct
    primary resources (``"resources"``) or ``1 + log2(1 + n references)``
    (``"references"``)."""

    if mode == "uniform":
        return 1.0
    if mode == "resources":
        return float(max(1, len({token.partition("_")[0] for token in _source_tokens(attrs.get("sources"))})))
    if mode == "references":
        refs = attrs.get("references")
        n = len({token.split(":")[-1] for token in re.split(r"[;|]", refs) if token.strip()}) if isinstance(refs, str) else 0
        return 1.0 + float(np.log2(1 + n))
    raise ValueError(f"edge_weight must be 'uniform', 'resources' or 'references', got {mode!r}.")


def cellchat_receptor_subunits(receptor: str) -> list[str]:
    """Gene symbols of a CellChat receptor name.

    CellChat joins the subunits of a receptor complex with ``_`` (``IL6R_IL6ST``,
    ``IL2RA_IL2RB_IL2RG``) and abbreviates some (``TGFbR1_R2`` is TGFBR1 + TGFBR2,
    ``ACVR1_TGFbR`` is ACVR1 + TGFBR1 + TGFBR2).
    """

    out, prev = [], None
    for token in str(receptor).split("_"):
        token = token.strip().upper()
        if not token:
            continue
        if re.fullmatch(r"R\d+", token) and prev is not None:
            token = re.sub(r"R\d+$", "", prev) + token
        if token == "TGFBR":
            out += ["TGFBR1", "TGFBR2"]
        else:
            out.append(token)
        prev = token
    return list(dict.fromkeys(out))


def ligand_receptor_seed(
    lr_df: pd.DataFrame,
    ligand: str,
    *,
    nodes=None,
    ligand_col: str = "ligand",
    receptor_col: str = "receptor",
) -> dict[str, float]:
    """Seed mass on a ligand's receptors: each receptor (complex) gets an equal share,
    split equally over its subunits (restricted to ``nodes`` when given)."""

    if lr_df is None or len(lr_df) == 0:
        return {}
    df = lr_df[[ligand_col, receptor_col]].dropna()
    receptors = df.loc[df[ligand_col].astype(str).str.strip() == str(ligand), receptor_col].astype(str).str.strip()
    complexes = []
    for receptor in dict.fromkeys(receptors):
        if not receptor:
            continue
        parts = [receptor] if nodes is not None and receptor in nodes else cellchat_receptor_subunits(receptor)
        if nodes is not None:
            parts = [part for part in parts if part in nodes]
        if parts:
            complexes.append(parts)
    seed = {}
    for parts in complexes:
        for part in parts:
            seed[part] = seed.get(part, 0.0) + 1.0 / len(complexes) / len(parts)
    return seed


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
    conflicting_signs: str = "consensus",
    evidence: str = "all",
    second_messenger_bridges: bool = True,
) -> OmniPathTfGraph:
    """Build the split signed OmniPath graph and corresponding TF list.

    ``conflicting_signs="consensus"`` (default) gives edges flagged both
    stimulatory and inhibitory OmniPath's consensus sign
    (:func:`resolve_conflicting_signs`); ``"stimulation"`` keeps both flags, which
    the diffusion treats as activating (ChronoPack before 2026-10-21).

    ``evidence="all"`` (default) keeps every directed signed edge. OmniPath's
    signaling edges are not all physical: many come only from pathway (activity
    flow) resources that record that A acts on B, directly or not.
    ``"physical"`` keeps edges supported by at least one protein-protein,
    enzyme-substrate or ligand-receptor resource (text-mined or predicted ones
    included), ``"physical_curated"`` only by curated or experimental ones
    (:func:`evidence_class`).

    ``second_messenger_bridges=True`` (default) adds the activating steps through
    cAMP, IP3/Ca2+ and DAG that a protein network lacks
    (:data:`SECOND_MESSENGER_BRIDGES`: adenylyl cyclases -> PKA/EPAC, phospholipase C
    -> calcineurin/CaMK and -> PKC/RasGRP), where OmniPath has no edge; ``False``
    gives the graph used before 2026-10-21.
    """

    if conflicting_signs not in ("consensus", "stimulation"):
        raise ValueError(f"conflicting_signs must be 'consensus' or 'stimulation', got {conflicting_signs!r}.")
    if evidence not in ("all", "physical", "physical_curated"):
        raise ValueError(f"evidence must be 'all', 'physical' or 'physical_curated', got {evidence!r}.")
    omnipath_tfs = set(resources.transcriptional["source_genesymbol"].unique())
    overlap = build_overlap_with_complexes(resources.human_tfs, omnipath_tfs)
    overlap_tfs = overlap["overlap_tfs"]
    interactions = filter_signaling_interactions(resources.signaling, drop_sources=drop_sources)
    if conflicting_signs == "consensus":
        interactions = resolve_conflicting_signs(interactions)
    if evidence != "all":
        keep = interactions["sources"].map(
            lambda value: has_physical_evidence(value, include_text_mined=evidence == "physical"))
        interactions = interactions[keep.astype(bool)].copy()

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
    if second_messenger_bridges:
        present = set(zip(merged["u"], merged["v"]))
        bridges = pd.DataFrame(
            [(u, v, True, False, 1, f"ChronoPack_second_messenger_{via}", np.nan, False)
             for u, v, via in SECOND_MESSENGER_BRIDGES if (u, v) not in present],
            columns=["u", "v", "is_stimulation", "is_inhibition", "curation_effort", "sources", "references",
                     "unsigned_for_diffusion"],
        )
        merged = pd.concat([merged, bridges[[c for c in bridges.columns if c in merged.columns or c in ("u", "v")]]],
                           ignore_index=True)
    graph = build_digraph_from_merged(merged)
    tfs = sorted(explode_tf_set(overlap_tfs, known_symbols=known_symbols) & set(graph.nodes()))

    return OmniPathTfGraph(
        graph=graph,
        tfs=tfs,
        interactions=interactions,
        split_interactions=expanded,
        merged_interactions=merged,
    )


def restrict_tf_graph(
    graph_bundle: OmniPathTfGraph,
    genes,
    *,
    scope: str = "network",
    receptor_genes=None,
) -> OmniPathTfGraph:
    """Remove the nodes whose genes are not all in ``genes`` (e.g. those expressed in a cell type).

    ``scope="network"`` applies this to every node (receptors, signaling proteins
    and TFs). ``scope="tfs"`` applies it only to the ends of the network: the TFs
    and the receptors in ``receptor_genes``; signaling proteins in between stay,
    since many are barely detected by single-cell RNA-seq. A complex node such as
    ``IRF9_STAT1_STAT2`` counts as expressed when every subunit is in ``genes``.
    """

    if scope not in ("network", "tfs"):
        raise ValueError(f"scope must be 'network' or 'tfs', got {scope!r}.")
    genes = {str(gene) for gene in genes}
    ends = set(graph_bundle.tfs) | {str(gene) for gene in (receptor_genes or ())}

    def expressed(node):
        return all(part in genes for part in str(node).split("_") if part)

    keep = [node for node in graph_bundle.graph.nodes
            if (scope == "tfs" and node not in ends) or expressed(node)]
    graph = graph_bundle.graph.subgraph(keep).copy()
    kept = set(keep)
    merged = graph_bundle.merged_interactions
    merged = merged[merged["u"].isin(kept) & merged["v"].isin(kept)].copy()
    return OmniPathTfGraph(
        graph=graph,
        tfs=[tf for tf in graph_bundle.tfs if tf in kept],
        interactions=graph_bundle.interactions,
        split_interactions=graph_bundle.split_interactions,
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
    edge_weight: str = "uniform",
    tied_signs: str = "split",
) -> tuple[dict, dict, dict, dict]:
    """Build positive/negative transition maps for signed absorbing diffusion.

    Each node's outgoing mass is split over its edges in proportion to
    ``w_signed`` / ``w_unsigned`` times the edge's evidence weight
    (:func:`edge_evidence_weight`; ``"uniform"`` by default).

    Edges flagged both stimulatory and inhibitory (resources disagree and
    OmniPath's consensus is tied, or duplicates merged with both signs) carry
    signal of unknown sign with ``tied_signs="split"`` (default): half of the mass
    keeps its sign and half flips, so the edge passes mass on but adds no net
    direction. ``"stimulation"`` treats them as activating (ChronoPack before
    2026-10-22).
    """

    if tied_signs not in ("stimulation", "split"):
        raise ValueError(f"tied_signs must be 'stimulation' or 'split', got {tied_signs!r}.")

    tfs = set(tfs)
    lr_receptors = {}
    if is_ligand and lr_df is not None and len(lr_df):
        # Receptor complexes (e.g. IL6R_IL6ST) are seeded at their subunits; see ligand_receptor_seed.
        for ligand in lr_df[ligand_col].dropna().astype(str).str.strip().unique():
            seed = ligand_receptor_seed(lr_df, ligand, nodes=graph if restrict_lr_to_nodes_in_G else None,
                                        ligand_col=ligand_col, receptor_col=receptor_col)
            if seed:
                lr_receptors[ligand] = list(seed)

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
            attrs = graph.get_edge_data(source, target, {}) or {}
            sign = edge_sign_from_attrs(attrs)
            evidence = edge_evidence_weight(attrs, edge_weight)
            if tied_signs == "split" and attrs.get("is_stimulation") and attrs.get("is_inhibition"):
                pos_edges.append((target, w_signed * evidence / 2))
                neg_edges.append((target, w_signed * evidence / 2))
            elif sign == -1:
                neg_edges.append((target, w_signed * evidence))
            elif sign == +1:
                pos_edges.append((target, w_signed * evidence))
            elif treat_unsigned_as_positive:
                pos_edges.append((target, w_unsigned * evidence))
            else:
                pos_edges.append((target, w_unsigned * evidence))

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
    """Run signed impulse diffusion with TF nodes as absorbing states.

    Positive mass moves along activating edges and switches sign along
    inhibitory ones; mass reaching a TF is absorbed. Each step is two sparse
    matrix-vector products (the same result as following every edge in turn).
    """

    idx = {node: i for i, node in enumerate(nodes)}
    n = len(nodes)
    x_pos = np.zeros(n, dtype=float)
    x_neg = np.zeros(n, dtype=float)

    for node, value in seed.items():
        if node in idx:
            x_pos[idx[node]] += float(value)

    total = x_pos.sum()
    if total <= 0:
        raise ValueError("Seed has zero mass; seed nodes are not in graph.")
    x_pos /= total

    def transposed(transitions):
        rows, cols, vals = [], [], []
        for source, edges in transitions.items():
            for target, prob in edges:
                rows.append(idx[target])
                cols.append(idx[source])
                vals.append(prob)
        return sparse.csr_matrix((vals, (rows, cols)), shape=(n, n))

    pos_t, neg_t = transposed(t_pos), transposed(t_neg)
    tf_names = list(dict.fromkeys(tf for tf in tfs if tf in idx))
    tf_idx = np.array([idx[tf] for tf in tf_names], dtype=int)
    absorbed_pos_arr = np.zeros(len(tf_names))
    absorbed_neg_arr = np.zeros(len(tf_names))
    reached_pos = np.zeros(len(tf_names), dtype=bool)
    reached_neg = np.zeros(len(tf_names), dtype=bool)
    hist = []

    for step in range(max_steps):
        y_pos = pos_t @ x_pos + neg_t @ x_neg
        y_neg = pos_t @ x_neg + neg_t @ x_pos

        ap = y_pos[tf_idx]
        an = y_neg[tf_idx]
        absorbed_pos_arr += ap
        absorbed_neg_arr += an
        reached_pos |= ap != 0
        reached_neg |= an != 0
        absorbed_step = float(ap.sum() + an.sum())
        y_pos[tf_idx] = 0.0
        y_neg[tf_idx] = 0.0

        x_pos = alpha * y_pos
        x_neg = alpha * y_neg
        transient_mass = float(x_pos.sum() + x_neg.sum())
        absorbed_total = float(absorbed_pos_arr.sum() + absorbed_neg_arr.sum())
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

    absorbed_pos = {tf: float(v) for tf, v, hit in zip(tf_names, absorbed_pos_arr, reached_pos) if hit}
    absorbed_neg = {tf: float(v) for tf, v, hit in zip(tf_names, absorbed_neg_arr, reached_neg) if hit}
    return absorbed_pos, absorbed_neg, x_pos, x_neg, pd.DataFrame(hist)


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
    edge_weight: str = "uniform",
    tied_signs: str = "split",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rank TFs by signed absorbing diffusion from a source node or ligand."""

    tfs = list(tfs)
    t_pos, t_neg, _, _ = build_signed_row_normalized_transitions(
        graph,
        tfs=tfs,
        is_ligand=False,  # ligand seeds come from ligand_receptor_seed below
        lr_df=lr_df,
        ligand_col=ligand_col,
        receptor_col=receptor_col,
        restrict_lr_to_nodes_in_G=restrict_lr_to_nodes_in_G,
        w_signed=w_signed,
        w_unsigned=w_unsigned,
        treat_unsigned_as_positive=True,
        edge_weight=edge_weight,
        tied_signs=tied_signs,
    )

    nodes = list(graph.nodes)
    seed = {}
    if is_ligand:
        seed = ligand_receptor_seed(lr_df, source, nodes=graph if restrict_lr_to_nodes_in_G else None,
                                    ligand_col=ligand_col, receptor_col=receptor_col)
        if not seed:
            raise ValueError(f"No ligand-receptor targets found for ligand '{source}'.")
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

    if sum(absorbed_pos.values()) + sum(absorbed_neg.values()) <= 0:
        dead = [node for node in seed if node in graph and graph.out_degree(node) == 0]
        raise ValueError(
            f"No signal from '{source}' reaches any TF (seeded at {sorted(seed)}"
            + (f"; without outgoing OmniPath edges: {sorted(dead)}" if dead else "")
            + "), so every TF would score 0."
        )
    if not any(absorbed_pos.get(tf, 0.0) > absorbed_neg.get(tf, 0.0) for tf in tfs):
        warnings.warn(f"No TF receives net activating signal from '{source}'; priority_act is 0 for every TF.",
                      stacklevel=2)

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
    edge_weight: str = "uniform",
    tied_signs: str = "split",
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
        edge_weight=edge_weight,
        tied_signs=tied_signs,
    )
