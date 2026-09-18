"""Plotting wrappers for ChronoPack results."""

from __future__ import annotations

from chronopack._backends.r import get_r_resource_path, read_r_resource_text


def r_plotting_resource_path(name="dimred_plotting_helpers.R"):
    """Yield a path to a bundled R plotting helper resource."""
    return get_r_resource_path(name)


def read_r_plotting_resource(name="dimred_plotting_helpers.R"):
    """Read a bundled R plotting helper resource."""
    return read_r_resource_text(name)


__all__ = ["r_plotting_resource_path", "read_r_plotting_resource"]
