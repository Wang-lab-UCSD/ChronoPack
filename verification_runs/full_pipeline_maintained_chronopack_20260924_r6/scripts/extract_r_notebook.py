#!/usr/bin/env python
"""Extract R code cells from a notebook to a plain R script."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("notebook", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    notebook = args.notebook.resolve()
    output = args.output.resolve()
    data = json.loads(notebook.read_text())
    chunks: list[str] = [
        f"# Extracted from {notebook}\n",
        "options(error = function() { traceback(2); quit(status = 1) })\n",
    ]
    for idx, cell in enumerate(data.get("cells", []), start=1):
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        if not source.strip():
            continue
        chunks.append(f"\nmessage('--- running {notebook.name} cell {idx} ---')\n")
        chunks.append(source)
        if not source.endswith("\n"):
            chunks.append("\n")

    output.write_text("".join(chunks))


if __name__ == "__main__":
    main()
