#!/usr/bin/env python
"""Execute notebook code cells sequentially with simple timing output."""

from __future__ import annotations

import json
import runpy
import sys
import time
import traceback
from pathlib import Path

from IPython.display import display


def main() -> int:
    if len(sys.argv) != 2:
        print("usage: run_notebook_cells.py NOTEBOOK.ipynb", file=sys.stderr)
        return 2

    notebook_path = Path(sys.argv[1]).resolve()
    workdir = notebook_path.parent
    sys.path.insert(0, str(workdir))

    nb = json.loads(notebook_path.read_text())
    namespace = {"__name__": "__main__", "display": display}

    for idx, cell in enumerate(nb["cells"]):
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        if not source.strip():
            continue

        first = source.strip().splitlines()[0][:100]
        print(f"START cell {idx}: {first}", flush=True)
        started = time.time()
        try:
            exec(compile(source, f"{notebook_path.name}:cell_{idx}", "exec"), namespace)
        except Exception:
            elapsed = time.time() - started
            print(f"ERROR cell {idx} after {elapsed:.1f}s", flush=True)
            traceback.print_exc()
            return 1
        elapsed = time.time() - started
        print(f"END cell {idx}: {elapsed:.1f}s", flush=True)

    print("SCRIPT COMPLETE", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
