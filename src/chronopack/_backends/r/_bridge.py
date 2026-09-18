"""Optional rpy2 bridge for bundled R backend functions."""

from __future__ import annotations

from ._resources import get_r_resource_path


def source_r_resource(name: str):
    """Source a bundled R file into the active rpy2 global environment."""
    try:
        from rpy2.robjects import r
    except ImportError as exc:
        raise ImportError("R backend calls require the optional 'rpy2' dependency.") from exc
    with get_r_resource_path(name) as path:
        return r["source"](str(path))
