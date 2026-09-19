"""Plotting wrappers for cross-context dynamic-mode comparison results."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from chronopack._backends.r._bridge import source_r_resource
from chronopack.pl._gene_trajectory_plots import RPlotResult


def _result_table(result_or_table, name: str) -> pd.DataFrame:
    if isinstance(result_or_table, pd.DataFrame):
        return result_or_table.copy()
    if hasattr(result_or_table, name):
        return getattr(result_or_table, name).copy()
    raise TypeError(f"Expected a DataFrame or an object with attribute '{name}'.")


def _render(plot_function: str, data: pd.DataFrame, output, **kwargs) -> str | None:
    if output is None:
        return None
    output = str(Path(output))
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    try:
        from rpy2.robjects import r
    except ImportError as exc:
        raise ImportError("R-backed plotting requires optional dependency 'rpy2'.") from exc
    source_r_resource("dimred_plotting_helpers.R")
    with TemporaryDirectory() as tmpdir:
        data_path = Path(tmpdir) / "plot_data.csv"
        data.to_csv(data_path, index=False)
        r[plot_function](str(data_path), output, **kwargs)
    return output


def mode_conservation(
    result,
    *,
    output: str | Path | None = None,
    width: float = 7,
    height: float = 4,
    title: str | None = None,
) -> RPlotResult:
    """Plot bidirectional dynamic-mode conservation scores."""
    data = _result_table(result, "conservation_scores")
    rendered = _render(
        "chronopack_plot_mode_conservation_csv",
        data,
        output,
        width=float(width),
        height=float(height),
        title="" if title is None else str(title),
    )
    return RPlotResult("mode_conservation", data=data, output=rendered)


def mode_sankey(
    result,
    *,
    output: str | Path | None = None,
    width: float = 7,
    height: float = 5,
    title: str | None = None,
    drop_low_significance: bool = False,
    hide_low_significance: bool = False,
    low_significance_label: str = "low_significance",
) -> RPlotResult:
    """Render source-target dynamic-mode flow counts from a result object."""
    data = _result_table(result, "mode_transitions")
    if drop_low_significance:
        data = data[
            (data["source_mode"].astype(str) != low_significance_label)
            & (data["target_mode"].astype(str) != low_significance_label)
        ].copy()
    elif hide_low_significance:
        data = data[
            ~(
                (data["source_mode"].astype(str) == low_significance_label)
                & (data["target_mode"].astype(str) == low_significance_label)
            )
        ].copy()
    rendered = _render(
        "chronopack_plot_mode_sankey_csv",
        data,
        output,
        width=float(width),
        height=float(height),
        title="" if title is None else str(title),
    )
    return RPlotResult(
        "mode_sankey",
        data=data,
        output=rendered,
        metadata={
            "drop_low_significance": drop_low_significance,
            "hide_low_significance": hide_low_significance,
            "low_significance_label": low_significance_label,
        },
    )


def mode_chord(
    result,
    *,
    output: str | Path | None = None,
    width: float = 6,
    height: float = 6,
    title: str | None = None,
) -> RPlotResult:
    """Render a chord-style source-target mode transition display when circlize is available."""
    data = _result_table(result, "mode_transitions")
    rendered = _render(
        "chronopack_plot_mode_chord_csv",
        data,
        output,
        width=float(width),
        height=float(height),
        title="" if title is None else str(title),
    )
    return RPlotResult("mode_chord", data=data, output=rendered)


__all__ = ["mode_chord", "mode_conservation", "mode_sankey"]
