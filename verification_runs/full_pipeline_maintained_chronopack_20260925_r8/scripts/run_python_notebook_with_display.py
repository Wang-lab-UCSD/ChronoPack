#!/usr/bin/env python
"""Execute code cells from a notebook with a minimal IPython-like display."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def _display(*objects):
    for obj in objects:
        print(obj)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("notebook", type=Path)
    args = parser.parse_args()

    notebook = args.notebook.resolve()
    data = json.loads(notebook.read_text())
    namespace = {
        "__name__": "__main__",
        "__file__": str(notebook),
        "display": _display,
    }

    for idx, cell in enumerate(data.get("cells", []), start=1):
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        if not source.strip():
            continue
        print(f"\n--- running {notebook.name} cell {idx} ---", flush=True)
        exec(compile(source, f"{notebook}:cell-{idx}", "exec"), namespace)


if __name__ == "__main__":
    main()
