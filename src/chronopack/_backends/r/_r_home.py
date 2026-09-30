"""Locate the R installation for rpy2 when the environment is not activated."""

from __future__ import annotations

import os
import shutil
import sys
from pathlib import Path


def ensure_r_home() -> str | None:
    """Set ``R_HOME`` to the running environment's R if rpy2 could not find R.

    rpy2 reads ``R_HOME`` or runs ``R RHOME`` from ``PATH``. A Jupyter kernel
    that starts an environment's Python without activating the environment has
    neither, so R calls fail with "Unable to determine R_HOME". When ``R_HOME``
    is unset and no ``R`` is on ``PATH``, point it at ``<sys.prefix>/lib/R``
    (``Lib/R`` on Windows) if that exists. Returns the ``R_HOME`` in effect, or
    ``None`` if it is left for rpy2 to resolve.
    """
    if os.environ.get("R_HOME"):
        return os.environ["R_HOME"]
    if shutil.which("R"):
        return None
    for candidate in (Path(sys.prefix) / "lib" / "R", Path(sys.prefix) / "Lib" / "R"):
        if (candidate / "bin").is_dir():
            os.environ["R_HOME"] = str(candidate)
            return str(candidate)
    return None
