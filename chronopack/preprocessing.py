"""Preprocessing helpers shared across Chrono-seq notebooks."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad


def map_time_to_group(time_tag):
    """Map 1..12 time tags onto the four coarse groups used by the notebooks."""
    time_tag = int(time_tag)
    if 1 <= time_tag <= 3:
        return "Group 1 (T1-3)"
    if 4 <= time_tag <= 6:
        return "Group 2 (T4-6)"
    if 7 <= time_tag <= 9:
        return "Group 3 (T7-9)"
    if 10 <= time_tag <= 12:
        return "Group 4 (T10-12)"
    return "Unknown Group"


TIME_GROUP_CATEGORIES = [
    "Group 1 (T1-3)",
    "Group 2 (T4-6)",
    "Group 3 (T7-9)",
    "Group 4 (T10-12)",
]


def annotate_qc_flags(adata, mt_prefixes=("MT-", "mt-"), rb_prefixes=("RPS", "RPL")):
    """Add mitochondrial/ribosomal flags and Scanpy QC metrics in place."""
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefixes)
    adata.var["rb"] = adata.var_names.str.startswith(rb_prefixes)
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "rb"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    return adata


def score_cell_cycle(adata, gene_set_path):
    """Score cell cycle using the Regev-lab gene list convention used here."""
    cell_cycle_genes = [x.strip() for x in open(gene_set_path)]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    s_genes = [g for g in s_genes if g in adata.var_names]
    g2m_genes = [g for g in g2m_genes if g in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    return adata


def add_time_metadata(adata, time_col="FINAL TIME TAG"):
    """Standardize per-cell time labels used by WOT and tradeSeq prep."""
    adata.obs["time"] = pd.to_numeric(adata.obs[time_col]).astype(int)
    adata.obs["time_group"] = pd.Categorical(
        adata.obs[time_col].apply(map_time_to_group),
        categories=TIME_GROUP_CATEGORIES,
        ordered=True,
    )
    return adata


def load_time_tagged_anndata(sample_paths, sample_names=None, metadata_cols=("FINAL TIME TAG", "COUNTS")):
    """Load time-tagged DGE tables into one AnnData object with sample-prefixed cells."""
    if sample_names is None:
        sample_names = [Path(path).name.split("_", 1)[0] for path in sample_paths]
    if len(sample_paths) != len(sample_names):
        raise ValueError("sample_paths and sample_names must have the same length")

    frames = []
    for path, sample in zip(sample_paths, sample_names):
        df = pd.read_csv(path, sep=",", index_col=0)
        df.index = f"{sample}_" + df.index.astype(str)
        frames.append(df)

    total_rows = sum(frame.shape[0] for frame in frames)
    unique_rows = len(set().union(*[set(frame.index) for frame in frames]))
    if total_rows != unique_rows:
        raise ValueError("Sample-prefixed cell IDs are not unique")

    df = pd.concat(frames).fillna(0)
    metadata = df.loc[:, list(metadata_cols)].copy()
    mtx = df.drop(columns=list(metadata_cols))
    adata = ad.AnnData(mtx)
    adata.obs = metadata

    source_map = {}
    for frame, sample in zip(frames, sample_names):
        source_map.update({bc: sample for bc in frame.index})
    adata.obs["sample"] = adata.obs_names.map(source_map).astype("category")
    add_time_metadata(adata, time_col="FINAL TIME TAG")
    adata.obs.pop("FINAL TIME TAG")
    return adata


def filter_cells_by_qc(adata, min_counts=None, max_counts=None, max_pct_mt=None, min_genes=None):
    """Apply the common cell-level QC thresholds with one explicit mask."""
    mask = np.ones(adata.n_obs, dtype=bool)
    if min_counts is not None:
        mask &= adata.obs["total_counts"].to_numpy() >= min_counts
    if max_counts is not None:
        mask &= adata.obs["total_counts"].to_numpy() <= max_counts
    if max_pct_mt is not None:
        mask &= adata.obs["pct_counts_mt"].to_numpy() < max_pct_mt
    if min_genes is not None:
        mask &= adata.obs["n_genes_by_counts"].to_numpy() >= min_genes
    return adata[mask, :].copy()


def filter_technical_genes(adata, exclude_genes=()):
    """Remove excluded genes plus MT, ENSG-only, and ribosomal features."""
    mask = ~np.isin(adata.var_names, list(exclude_genes))
    mask &= ~adata.var_names.str.contains("^MT-")
    mask &= ~adata.var_names.str.contains("^ENSG")
    mask &= ~adata.var_names.str.contains(r"^RP[SL][0-9]")
    return adata[:, mask].copy()


def filter_low_count_genes(adata, min_total_counts=10, counts_layer="counts"):
    """Remove genes whose total counts across retained cells are below threshold."""
    if counts_layer is not None and counts_layer in adata.layers:
        counts = adata.layers[counts_layer]
    else:
        counts = adata.X
    total_counts = np.asarray(counts.sum(axis=0)).ravel()
    adata.var["total_counts"] = total_counts
    keep = total_counts >= min_total_counts
    return adata[:, keep].copy()
