"""Deprecation helpers for compatibility modules."""

from __future__ import annotations

import warnings


def warn_legacy_module(module_name: str, replacement: str) -> None:
    """Warn that a top-level compatibility module is not the preferred API."""
    warnings.warn(
        (
            f"{module_name} is a legacy compatibility interface and is not "
            f"recommended for new code. Use {replacement} instead."
        ),
        FutureWarning,
        stacklevel=3,
    )
