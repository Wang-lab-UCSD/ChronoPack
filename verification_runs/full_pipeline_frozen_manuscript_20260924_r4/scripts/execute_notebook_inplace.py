#!/usr/bin/env python
"""Execute a Python notebook in place without nbclient/nbformat.

This runner is intentionally small and is only meant for the ChronoPack
verification notebooks in this repository. It records stdout/stderr, explicit
display() calls, the final expression value, and open matplotlib figures.
"""

from __future__ import annotations

import ast
import base64
import contextlib
import io
import json
import sys
import traceback
from pathlib import Path


def _stream_output(name: str, text: str) -> dict | None:
    if not text:
        return None
    return {"output_type": "stream", "name": name, "text": text.splitlines(keepends=True)}


def _format_text(obj) -> str:
    try:
        return repr(obj)
    except Exception:
        return f"<unrepresentable {type(obj).__name__}>"


def _execute_result(obj, execution_count: int) -> dict:
    return {
        "output_type": "execute_result",
        "execution_count": execution_count,
        "data": {"text/plain": _format_text(obj)},
        "metadata": {},
    }


def _display_data(obj) -> dict:
    data = {"text/plain": _format_text(obj)}
    if hasattr(obj, "_repr_html_"):
        try:
            html = obj._repr_html_()
            if html is not None:
                data["text/html"] = html
        except Exception:
            pass
    return {"output_type": "display_data", "data": data, "metadata": {}}


def _figure_outputs() -> list[dict]:
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return []

    outs = []
    for num in list(plt.get_fignums()):
        fig = plt.figure(num)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight", dpi=120)
        outs.append(
            {
                "output_type": "display_data",
                "data": {"image/png": base64.b64encode(buf.getvalue()).decode("ascii")},
                "metadata": {},
            }
        )
    plt.close("all")
    return outs


def _split_final_expr(source: str):
    tree = ast.parse(source, mode="exec")
    if not tree.body or not isinstance(tree.body[-1], ast.Expr):
        return source, None
    prefix = ast.Module(body=tree.body[:-1], type_ignores=[])
    ast.fix_missing_locations(prefix)
    expr = ast.Expression(body=tree.body[-1].value)
    ast.fix_missing_locations(expr)
    return prefix, expr


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: execute_notebook_inplace.py NOTEBOOK.ipynb")

    notebook = Path(sys.argv[1]).resolve()
    data = json.loads(notebook.read_text())

    outputs: list[dict] = []

    def display(*objects):
        for obj in objects:
            outputs.append(_display_data(obj))

    namespace = {
        "__name__": "__main__",
        "__file__": str(notebook),
        "display": display,
    }

    execution_count = 0
    for idx, cell in enumerate(data.get("cells", []), start=1):
        if cell.get("cell_type") != "code":
            continue
        source = "".join(cell.get("source", []))
        cell["outputs"] = []
        if not source.strip():
            cell["execution_count"] = None
            continue

        execution_count += 1
        cell["execution_count"] = execution_count
        outputs = []
        stdout = io.StringIO()
        stderr = io.StringIO()

        print(f"--- executing {notebook.name} cell {idx} ---", flush=True)
        try:
            with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
                prefix, expr = _split_final_expr(source)
                if expr is None:
                    exec(compile(source, f"{notebook}:cell-{idx}", "exec"), namespace)
                else:
                    exec(compile(prefix, f"{notebook}:cell-{idx}", "exec"), namespace)
                    result = eval(compile(expr, f"{notebook}:cell-{idx}", "eval"), namespace)
                    if result is not None:
                        outputs.append(_execute_result(result, execution_count))
        except Exception as exc:
            tb = traceback.format_exception(type(exc), exc, exc.__traceback__)
            out = _stream_output("stdout", stdout.getvalue())
            if out:
                cell["outputs"].append(out)
            err = _stream_output("stderr", stderr.getvalue())
            if err:
                cell["outputs"].append(err)
            cell["outputs"].append(
                {
                    "output_type": "error",
                    "ename": type(exc).__name__,
                    "evalue": str(exc),
                    "traceback": tb,
                }
            )
            notebook.write_text(json.dumps(data, indent=1))
            raise

        out = _stream_output("stdout", stdout.getvalue())
        if out:
            cell["outputs"].append(out)
        err = _stream_output("stderr", stderr.getvalue())
        if err:
            cell["outputs"].append(err)
        cell["outputs"].extend(outputs)
        cell["outputs"].extend(_figure_outputs())

    notebook.write_text(json.dumps(data, indent=1))


if __name__ == "__main__":
    main()
