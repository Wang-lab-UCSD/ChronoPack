#!/usr/bin/env python
"""Execute a notebook and save a copy with its outputs.

    execute_notebook.py SOURCE.ipynb EXECUTED.ipynb

Runs the code cells of SOURCE in one namespace, from SOURCE's folder (so
relative paths resolve as they do in Jupyter), and writes EXECUTED with each
cell's stdout/stderr, displayed objects, final expression value and open
matplotlib figures. SOURCE is never modified. Kept dependency-free (no
nbclient / papermill) so it runs in the plain ChronoPack environment.
"""

from __future__ import annotations

import ast
import base64
import contextlib
import io
import json
import os
import sys
import traceback
from pathlib import Path


def _stream(name: str, text: str) -> dict | None:
    if not text:
        return None
    return {"output_type": "stream", "name": name, "text": text.splitlines(keepends=True)}


def _repr(obj) -> str:
    try:
        return repr(obj)
    except Exception:
        return f"<unrepresentable {type(obj).__name__}>"


def _mime(obj) -> tuple[dict, dict]:
    """(data, metadata) for an output, as Jupyter would record it."""
    data, metadata = {"text/plain": _repr(obj)}, {}
    # Rich representations Jupyter would display (HTML tables, SVG / PNG plots). A
    # _repr_*_ method may return (value, metadata), e.g. a PNG with a display width.
    for method, mime in (("_repr_html_", "text/html"), ("_repr_svg_", "image/svg+xml"), ("_repr_png_", "image/png")):
        if hasattr(obj, method):
            try:
                value = getattr(obj, method)()
            except Exception:
                value = None
            if isinstance(value, tuple):
                value, md = value
                if md:
                    metadata[mime] = md
            if value is not None:
                data[mime] = base64.b64encode(value).decode("ascii") if isinstance(value, bytes) else value
    return data, metadata


def _figures() -> list[dict]:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return []
    outs = []
    for num in list(plt.get_fignums()):
        buf = io.BytesIO()
        plt.figure(num).savefig(buf, format="png", bbox_inches="tight", dpi=100)
        outs.append({
            "output_type": "display_data",
            "data": {"image/png": base64.b64encode(buf.getvalue()).decode("ascii")},
            "metadata": {},
        })
    plt.close("all")
    return outs


def _split_last_expr(source: str):
    tree = ast.parse(source, mode="exec")
    if not tree.body or not isinstance(tree.body[-1], ast.Expr):
        return tree, None
    body = ast.Module(body=tree.body[:-1], type_ignores=[])
    return body, ast.Expression(body=tree.body[-1].value)


def main() -> None:
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    source = Path(sys.argv[1]).resolve()
    executed = Path(sys.argv[2]).resolve()
    executed.parent.mkdir(parents=True, exist_ok=True)
    nb = json.loads(source.read_text())
    os.chdir(source.parent)

    displayed: list[dict] = []
    namespace = {
        "__name__": "__main__",
        "display": lambda *objs: displayed.extend(
            dict(zip(("data", "metadata"), _mime(o)), output_type="display_data") for o in objs
        ),
    }

    def save():
        executed.write_text(json.dumps(nb, indent=1, ensure_ascii=False) + "\n")

    count = 0
    for idx, cell in enumerate(nb["cells"]):
        if cell["cell_type"] != "code":
            continue
        code = "".join(cell["source"])
        cell["outputs"] = []
        cell["execution_count"] = None
        if not code.strip():
            continue
        count += 1
        cell["execution_count"] = count
        displayed.clear()
        out, err = io.StringIO(), io.StringIO()
        label = f"{source.name}:cell-{idx}"
        print(f"--- {source.name} cell {idx}", flush=True)
        try:
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                body, last = _split_last_expr(code)
                exec(compile(body, label, "exec"), namespace)
                result = eval(compile(last, label, "eval"), namespace) if last is not None else None
        except Exception as exc:
            cell["outputs"] += [o for o in (_stream("stdout", out.getvalue()), _stream("stderr", err.getvalue())) if o]
            cell["outputs"].append({
                "output_type": "error",
                "ename": type(exc).__name__,
                "evalue": str(exc),
                "traceback": traceback.format_exception(type(exc), exc, exc.__traceback__),
            })
            save()
            raise
        cell["outputs"] += [o for o in (_stream("stdout", out.getvalue()), _stream("stderr", err.getvalue())) if o]
        cell["outputs"] += displayed
        if result is not None:
            cell["outputs"].append({
                "output_type": "execute_result", "execution_count": count, **dict(zip(("data", "metadata"), _mime(result))),
            })
        cell["outputs"] += _figures()
    save()


if __name__ == "__main__":
    main()
