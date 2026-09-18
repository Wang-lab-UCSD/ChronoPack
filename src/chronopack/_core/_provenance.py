"""Small provenance helpers for AnnData-mutating ChronoPack functions."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any

from ._keys import CHRONOPACK_UNS_KEY


def _h5ad_safe(value: Any) -> Any:
    """Convert small provenance payloads to AnnData/HDF5-writable values."""
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, tuple):
        return [_h5ad_safe(item) for item in value]
    if isinstance(value, list):
        return [_h5ad_safe(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _h5ad_safe(item) for key, item in value.items()}
    return str(value)


def record_provenance(
    adata,
    key: str,
    *,
    analysis: str,
    parameters: dict[str, Any] | None = None,
    result_keys: dict[str, Any] | None = None,
    backend: str | None = None,
    version: str | None = None,
) -> dict[str, Any]:
    """Record compact ChronoPack provenance under ``adata.uns``."""
    if CHRONOPACK_UNS_KEY not in adata.uns:
        adata.uns[CHRONOPACK_UNS_KEY] = {}
    entry = {
        "analysis": analysis,
        "chronopack_version": version,
        "backend": backend,
        "parameters": _h5ad_safe(dict(parameters or {})),
        "result_keys": _h5ad_safe(dict(result_keys or {})),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    adata.uns[CHRONOPACK_UNS_KEY][key] = entry
    return entry
