"""Package-resource access for bundled R backend files."""

from __future__ import annotations

from contextlib import contextmanager
from importlib import resources
from pathlib import Path
from typing import Iterator

_PACKAGE = __package__


def available_r_resources() -> tuple[str, ...]:
    """Return bundled R resource filenames."""
    root = resources.files(_PACKAGE)
    return tuple(sorted(path.name for path in root.iterdir() if path.name.endswith(".R")))


@contextmanager
def get_r_resource_path(name: str) -> Iterator[Path]:
    """Yield a filesystem path for a bundled R resource."""
    if not name.endswith(".R"):
        raise ValueError("R resource names must end with '.R'.")
    resource = resources.files(_PACKAGE).joinpath(name)
    if not resource.is_file():
        raise FileNotFoundError(f"Bundled R resource not found: {name}")
    with resources.as_file(resource) as path:
        yield path


def read_r_resource_text(name: str) -> str:
    """Read a bundled R resource as text."""
    resource = resources.files(_PACKAGE).joinpath(name)
    if not resource.is_file():
        raise FileNotFoundError(f"Bundled R resource not found: {name}")
    return resource.read_text()
