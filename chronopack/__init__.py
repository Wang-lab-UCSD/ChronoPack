"""Compatibility package for source-checkout ChronoPack imports.

The formal package now lives under ``src/chronopack``. This flat package is
kept temporarily so existing notebooks that import ``chronopack.*`` from a
source checkout continue to work during the staged migration.
"""

from __future__ import annotations

from pathlib import Path

_SRC_PACKAGE = Path(__file__).resolve().parents[1] / "src" / "chronopack"
if _SRC_PACKAGE.exists():
    src_path = str(_SRC_PACKAGE)
    if src_path not in __path__:
        __path__.append(src_path)

__version__ = "0.1.0"

__all__ = ["io", "tl", "pl", "__version__"]


def __getattr__(name):
    if name in {"io", "tl", "pl"}:
        import importlib

        return importlib.import_module(f"{__name__}.{name}")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
