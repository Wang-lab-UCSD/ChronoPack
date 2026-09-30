"""Bundled R backend resources."""

from ._r_home import ensure_r_home
from ._resources import available_r_resources, get_r_resource_path, read_r_resource_text

# Every R-backed module imports this package before rpy2, so R_HOME is in place
# before rpy2 initializes (see ensure_r_home).
ensure_r_home()

__all__ = ["available_r_resources", "ensure_r_home", "get_r_resource_path", "read_r_resource_text"]
