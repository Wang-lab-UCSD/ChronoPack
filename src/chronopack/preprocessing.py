"""Temporary compatibility imports for preprocessing helpers."""

from chronopack._utils._deprecation import warn_legacy_module

warn_legacy_module("chronopack.preprocessing", "chronopack.io")

from chronopack.io import *  # noqa: F401,F403
