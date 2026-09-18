"""Central AnnData validation helpers."""

from __future__ import annotations

import numpy as np

from ._keys import DEFAULT_COUNTS_LAYER, DEFAULT_TIME_KEY
from ._schema import PREPROCESSING_SCHEMA, SAMPLE_KEY, TIME_GROUP_KEY


class ChronoPackError(ValueError):
    """Raised when ChronoPack input validation fails."""


def validate_adata(
    adata,
    *,
    time_key: str | None = None,
    layer: str | None = None,
    require_unique: bool = True,
):
    """Validate common AnnData assumptions and return ``adata``."""
    if require_unique:
        if not adata.obs_names.is_unique:
            raise ChronoPackError("AnnData obs_names must be unique.")
        if not adata.var_names.is_unique:
            raise ChronoPackError("AnnData var_names must be unique.")
    if time_key is not None:
        if time_key not in adata.obs:
            raise ChronoPackError(f"Missing adata.obs[{time_key!r}].")
        vals = np.asarray(adata.obs[time_key])
        if len(vals) != adata.n_obs:
            raise ChronoPackError(f"adata.obs[{time_key!r}] has invalid length.")
    if layer is not None and layer not in adata.layers:
        raise ChronoPackError(f"Missing adata.layers[{layer!r}].")
    if adata.X.shape != (adata.n_obs, adata.n_vars):
        raise ChronoPackError("adata.X shape is inconsistent with obs/var dimensions.")
    for key, layer_value in adata.layers.items():
        if layer_value.shape != (adata.n_obs, adata.n_vars):
            raise ChronoPackError(f"adata.layers[{key!r}] shape is inconsistent with adata.X.")
    return adata


def validate_numeric_obs(adata, key: str):
    """Validate that an obs column exists and is numeric/coercible."""
    if key not in adata.obs:
        raise ChronoPackError(f"Missing adata.obs[{key!r}].")
    vals = np.asarray(adata.obs[key])
    numeric = np.asarray(np.asarray(vals), dtype=object)
    try:
        numeric = numeric.astype(float)
    except (TypeError, ValueError) as exc:
        raise ChronoPackError(f"adata.obs[{key!r}] must be numeric or numeric-coercible.") from exc
    if len(numeric) != adata.n_obs:
        raise ChronoPackError(f"adata.obs[{key!r}] has invalid length.")
    if not np.all(np.isfinite(numeric)):
        raise ChronoPackError(f"adata.obs[{key!r}] contains non-finite values.")
    return adata


def validate_preprocessed_adata(
    adata,
    *,
    time_key: str = DEFAULT_TIME_KEY,
    counts_layer: str = DEFAULT_COUNTS_LAYER,
    require_layers: tuple[str, ...] | list[str] | None = None,
    require_qc: bool = True,
    require_sample: bool = False,
    require_time_group: bool = False,
    require_cell_type: bool = False,
):
    """Validate the minimal AnnData contract for current preprocessing outputs."""
    validate_adata(adata, require_unique=True)
    validate_numeric_obs(adata, time_key)

    if require_sample and SAMPLE_KEY not in adata.obs:
        raise ChronoPackError(f"Missing adata.obs[{SAMPLE_KEY!r}].")
    if require_time_group and TIME_GROUP_KEY not in adata.obs:
        raise ChronoPackError(f"Missing adata.obs[{TIME_GROUP_KEY!r}].")
    if require_cell_type and PREPROCESSING_SCHEMA.cell_type_key not in adata.obs:
        raise ChronoPackError(f"Missing adata.obs[{PREPROCESSING_SCHEMA.cell_type_key!r}].")

    if require_layers is None:
        require_layers = (counts_layer,)
    for layer in require_layers:
        if layer not in adata.layers:
            raise ChronoPackError(f"Missing adata.layers[{layer!r}].")
        if adata.layers[layer].shape != (adata.n_obs, adata.n_vars):
            raise ChronoPackError(f"adata.layers[{layer!r}] shape is inconsistent with adata.X.")

    if require_qc:
        obs_cols = ["total_counts", "n_genes_by_counts", "pct_counts_mt", "pct_counts_rb"]
        var_cols = [PREPROCESSING_SCHEMA.mito_flag_key, PREPROCESSING_SCHEMA.ribosomal_flag_key]
        missing_obs = [col for col in obs_cols if col not in adata.obs]
        missing_var = [col for col in var_cols if col not in adata.var]
        if missing_obs:
            raise ChronoPackError(f"Missing QC columns in adata.obs: {missing_obs}.")
        if missing_var:
            raise ChronoPackError(f"Missing QC columns in adata.var: {missing_var}.")

    return adata
