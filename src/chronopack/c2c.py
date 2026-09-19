"""Compatibility imports for C2C helpers."""

from chronopack._utils._deprecation import warn_legacy_module

warn_legacy_module("chronopack.c2c", "chronopack.tl.c2c and chronopack.pl.c2c")

from chronopack._legacy.c2c import *  # noqa: F401,F403
from chronopack.tl.c2c import C2CResult, infer  # noqa: F401
