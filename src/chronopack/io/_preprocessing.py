"""ChronoPack preprocessing I/O and AnnData preparation helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from chronopack._core._keys import DEFAULT_COUNTS_LAYER, DEFAULT_TIME_KEY
from chronopack._core._provenance import record_provenance
from chronopack._core._schema import (
    CPM_LAYER,
    LOG1P_LAYER,
    PREPROCESSING_SCHEMA,
    SAMPLE_KEY,
    TIME_GROUP_CATEGORIES,
    TIME_GROUP_KEY,
    TIME_TAG_KEY,
    map_time_to_group,
)
from chronopack._core._validation import ChronoPackError, validate_preprocessed_adata


def _target(adata, copy: bool):
    return adata.copy() if copy else adata


def read_time_tagged_dge(
    sample_paths: Sequence[str | Path],
    *,
    sample_names: Sequence[str] | None = None,
    metadata_cols: Sequence[str] = (TIME_TAG_KEY, "COUNTS"),
    time_col: str = TIME_TAG_KEY,
    sample_key: str = SAMPLE_KEY,
) -> ad.AnnData:
    """Read time-tagged DGE CSV tables into one AnnData object.

    The loader preserves the existing notebook convention: cell IDs are
    prefixed with the sample name, metadata columns are moved into ``obs``,
    expression columns become ``X``, missing genes across samples are filled
    with zero, and ``obs['time']``/``obs['time_group']`` are added from the time
    tag column.
    """
    sample_paths = [Path(path) for path in sample_paths]
    if sample_names is None:
        sample_names = [path.name.split("_", 1)[0] for path in sample_paths]
    if len(sample_paths) != len(sample_names):
        raise ChronoPackError("sample_paths and sample_names must have the same length.")
    if time_col not in metadata_cols:
        raise ChronoPackError(f"time_col {time_col!r} must be included in metadata_cols.")

    frames: list[pd.DataFrame] = []
    for path, sample in zip(sample_paths, sample_names):
        df = pd.read_csv(path, sep=",", index_col=0)
        missing = [col for col in metadata_cols if col not in df.columns]
        if missing:
            raise ChronoPackError(f"{path} is missing metadata columns: {missing}.")
        df.index = f"{sample}_" + df.index.astype(str)
        frames.append(df)

    total_rows = sum(frame.shape[0] for frame in frames)
    unique_rows = len(set().union(*(set(frame.index) for frame in frames)))
    if total_rows != unique_rows:
        raise ChronoPackError("Sample-prefixed cell IDs are not unique.")

    df = pd.concat(frames).fillna(0)
    metadata = df.loc[:, list(metadata_cols)].copy()
    matrix = df.drop(columns=list(metadata_cols))
    if matrix.columns.duplicated().any():
        dupes = matrix.columns[matrix.columns.duplicated()].unique().tolist()
        raise ChronoPackError(f"Duplicate gene identifiers after loading: {dupes[:5]}.")

    adata = ad.AnnData(matrix)
    adata.obs = metadata

    source_map = {}
    for frame, sample in zip(frames, sample_names):
        source_map.update({barcode: sample for barcode in frame.index})
    adata.obs[sample_key] = adata.obs_names.map(source_map).astype("category")
    add_time_metadata(adata, time_col=time_col)
    if time_col in adata.obs:
        adata.obs.pop(time_col)

    record_provenance(
        adata,
        "read_time_tagged_dge",
        analysis="preprocessing.read_time_tagged_dge",
        parameters={
            "n_samples": len(sample_paths),
            "metadata_cols": list(metadata_cols),
            "time_col": time_col,
            "sample_key": sample_key,
        },
        result_keys={"obs": [DEFAULT_TIME_KEY, TIME_GROUP_KEY, sample_key]},
    )
    validate_preprocessed_adata(adata, require_layers=(), require_qc=False, require_sample=True)
    return adata


def add_time_metadata(
    adata,
    *,
    time_col: str = TIME_TAG_KEY,
    time_key: str = DEFAULT_TIME_KEY,
    time_group_key: str = TIME_GROUP_KEY,
    copy: bool = False,
):
    """Add numeric time and coarse time-group metadata."""
    out = _target(adata, copy)
    if time_col not in out.obs:
        raise ChronoPackError(f"Missing adata.obs[{time_col!r}].")
    time = pd.to_numeric(out.obs[time_col], errors="raise").astype(int)
    out.obs[time_key] = time
    out.obs[time_group_key] = pd.Categorical(
        time.apply(map_time_to_group),
        categories=list(TIME_GROUP_CATEGORIES),
        ordered=True,
    )
    record_provenance(
        out,
        "add_time_metadata",
        analysis="preprocessing.add_time_metadata",
        parameters={"time_col": time_col, "time_key": time_key, "time_group_key": time_group_key},
        result_keys={"obs": [time_key, time_group_key]},
    )
    return out


def annotate_qc_metrics(
    adata,
    *,
    mt_prefixes: tuple[str, ...] = ("MT-", "mt-"),
    rb_prefixes: tuple[str, ...] = ("RPS", "RPL"),
    copy: bool = False,
):
    """Add mitochondrial/ribosomal flags and Scanpy QC metrics."""
    out = _target(adata, copy)
    out.var[PREPROCESSING_SCHEMA.mito_flag_key] = out.var_names.str.startswith(mt_prefixes)
    out.var[PREPROCESSING_SCHEMA.ribosomal_flag_key] = out.var_names.str.startswith(rb_prefixes)
    sc.pp.calculate_qc_metrics(
        out,
        qc_vars=[PREPROCESSING_SCHEMA.mito_flag_key, PREPROCESSING_SCHEMA.ribosomal_flag_key],
        percent_top=None,
        log1p=False,
        inplace=True,
    )
    record_provenance(
        out,
        "annotate_qc_metrics",
        analysis="preprocessing.annotate_qc_metrics",
        parameters={"mt_prefixes": mt_prefixes, "rb_prefixes": rb_prefixes},
        result_keys={
            "var": [PREPROCESSING_SCHEMA.mito_flag_key, PREPROCESSING_SCHEMA.ribosomal_flag_key],
            "obs": ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_rb"],
        },
    )
    return out


def filter_cells_by_qc(
    adata,
    *,
    min_counts: float | None = None,
    max_counts: float | None = None,
    max_pct_mt: float | None = None,
    min_genes: float | None = None,
    copy: bool = False,
):
    """Filter cells with the existing combined QC mask."""
    out = _target(adata, copy)
    required = ["total_counts", "pct_counts_mt", "n_genes_by_counts"]
    missing = [col for col in required if col not in out.obs]
    if missing:
        raise ChronoPackError(f"Missing QC columns in adata.obs: {missing}.")

    mask = np.ones(out.n_obs, dtype=bool)
    if min_counts is not None:
        mask &= out.obs["total_counts"].to_numpy() >= min_counts
    if max_counts is not None:
        mask &= out.obs["total_counts"].to_numpy() <= max_counts
    if max_pct_mt is not None:
        mask &= out.obs["pct_counts_mt"].to_numpy() < max_pct_mt
    if min_genes is not None:
        mask &= out.obs["n_genes_by_counts"].to_numpy() >= min_genes

    out._inplace_subset_obs(mask)
    record_provenance(
        out,
        "filter_cells_by_qc",
        analysis="preprocessing.filter_cells_by_qc",
        parameters={
            "min_counts": min_counts,
            "max_counts": max_counts,
            "max_pct_mt": max_pct_mt,
            "min_genes": min_genes,
        },
        result_keys={"n_obs": out.n_obs},
    )
    return out


def filter_technical_genes(
    adata,
    *,
    exclude_genes: Iterable[str] = (),
    copy: bool = False,
):
    """Remove the technical gene classes used by the current notebooks."""
    out = _target(adata, copy)
    exclude = set(exclude_genes)
    mask = ~np.isin(out.var_names, list(exclude))
    mask &= ~out.var_names.str.contains("^MT-")
    mask &= ~out.var_names.str.contains("^ENSG")
    mask &= ~out.var_names.str.contains(r"^RP[SL][0-9]")
    out._inplace_subset_var(mask)
    record_provenance(
        out,
        "filter_technical_genes",
        analysis="preprocessing.filter_technical_genes",
        parameters={"n_explicit_exclude_genes": len(exclude)},
        result_keys={"n_vars": out.n_vars},
    )
    return out


def filter_low_count_genes(
    adata,
    *,
    min_total_counts: float = 10,
    counts_layer: str | None = DEFAULT_COUNTS_LAYER,
    copy: bool = False,
):
    """Remove genes whose total retained-cell counts are below threshold."""
    out = _target(adata, copy)
    if counts_layer is not None:
        if counts_layer not in out.layers:
            raise ChronoPackError(f"Missing adata.layers[{counts_layer!r}].")
        counts = out.layers[counts_layer]
    else:
        counts = out.X
    total_counts = np.asarray(counts.sum(axis=0)).ravel()
    out.var["total_counts"] = total_counts
    out._inplace_subset_var(total_counts >= min_total_counts)
    record_provenance(
        out,
        "filter_low_count_genes",
        analysis="preprocessing.filter_low_count_genes",
        parameters={"min_total_counts": min_total_counts, "counts_layer": counts_layer},
        result_keys={"n_vars": out.n_vars, "var": ["total_counts"]},
    )
    return out


def score_cell_cycle(
    adata,
    gene_set_path: str | Path,
    *,
    copy: bool = False,
):
    """Score cell cycle using the Regev-lab gene-list split used here."""
    out = _target(adata, copy)
    with open(gene_set_path) as handle:
        cell_cycle_genes = [line.strip() for line in handle if line.strip()]
    s_genes = [gene for gene in cell_cycle_genes[:43] if gene in out.var_names]
    g2m_genes = [gene for gene in cell_cycle_genes[43:] if gene in out.var_names]
    sc.tl.score_genes_cell_cycle(out, s_genes=s_genes, g2m_genes=g2m_genes)
    record_provenance(
        out,
        "score_cell_cycle",
        analysis="preprocessing.score_cell_cycle",
        parameters={"gene_set_path": str(gene_set_path), "n_s_genes": len(s_genes), "n_g2m_genes": len(g2m_genes)},
        result_keys={"obs": ["S_score", "G2M_score", "phase"]},
    )
    return out


def initialize_expression_layers(
    adata,
    *,
    counts_layer: str = DEFAULT_COUNTS_LAYER,
    cpm_layer: str = CPM_LAYER,
    log1p_layer: str = LOG1P_LAYER,
    target_sum: float = 1e4,
    copy: bool = False,
):
    """Create the standard counts, normalized, and log1p layers.

    This wraps the recurring notebook idiom while preserving Scanpy for the
    actual normalization and log transform.
    """
    out = _target(adata, copy)
    out.layers[counts_layer] = out.X.copy()
    sc.pp.normalize_total(out, target_sum=target_sum)
    out.layers[cpm_layer] = out.X.copy()
    sc.pp.log1p(out)
    out.layers[log1p_layer] = out.X.copy()
    record_provenance(
        out,
        "initialize_expression_layers",
        analysis="preprocessing.initialize_expression_layers",
        parameters={"target_sum": target_sum},
        result_keys={"layers": [counts_layer, cpm_layer, log1p_layer]},
    )
    return out


# Backward-compatible names used by active notebooks before the formal API.
load_time_tagged_anndata = read_time_tagged_dge
annotate_qc_flags = annotate_qc_metrics


__all__ = [
    "TIME_GROUP_CATEGORIES",
    "add_time_metadata",
    "annotate_qc_flags",
    "annotate_qc_metrics",
    "filter_cells_by_qc",
    "filter_low_count_genes",
    "filter_technical_genes",
    "initialize_expression_layers",
    "load_time_tagged_anndata",
    "map_time_to_group",
    "read_time_tagged_dge",
    "score_cell_cycle",
    "validate_preprocessed_adata",
]
