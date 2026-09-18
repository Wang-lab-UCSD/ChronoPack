"""ChronoPack public Python API."""

__version__ = "0.1.0"

__all__ = ["io", "tl", "pl", "__version__"]


def __getattr__(name):
    if name in {"io", "tl", "pl"}:
        import importlib

        return importlib.import_module(f"{__name__}.{name}")
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
