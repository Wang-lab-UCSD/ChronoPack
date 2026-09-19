"""Compatibility imports for legacy cell-state trajectory helpers."""

from chronopack._utils._deprecation import warn_legacy_module

warn_legacy_module("chronopack.trajectories", "chronopack.tl")

from chronopack._legacy.trajectories import *  # noqa: F401,F403
