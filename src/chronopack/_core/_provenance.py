"""Small provenance helpers for AnnData-mutating ChronoPack functions."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any

from ._keys import CHRONOPACK_UNS_KEY


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
        "parameters": dict(parameters or {}),
        "result_keys": dict(result_keys or {}),
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    adata.uns[CHRONOPACK_UNS_KEY][key] = entry
    return entry
