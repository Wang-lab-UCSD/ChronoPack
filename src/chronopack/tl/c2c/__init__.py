"""Cell-cell communication analysis tools."""

from chronopack._legacy.c2c import *  # noqa: F401,F403

infer = score_lr_pairs  # type: ignore[name-defined]

__all__ = [name for name in globals() if not name.startswith("_")]
