"""Compatibility imports for OmniPath TF helpers."""

from chronopack._utils._deprecation import warn_legacy_module

warn_legacy_module("chronopack.omnipath_tf", "chronopack.tl.grn.rank_regulators")

from chronopack._legacy.omnipath_tf import *  # noqa: F401,F403
from chronopack.tl.grn import rank_regulators  # noqa: F401
