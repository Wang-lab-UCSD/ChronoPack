"""Central AnnData validation helpers."""

from __future__ import annotations

import numpy as np


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
    return adata
