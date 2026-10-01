"""Cross-context recovery of dynamic response modes."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd


# TESTING.ipynb `response_direction` of the named broad modes.
MODE_DIRECTIONS = {
    "impulse-up_1": "induced",
    "impulse-up_2": "induced",
    "trend-up": "induced",
    "trend-down_1": "repressed",
    "trend-down_2": "repressed",
    "transient-down": "repressed",
}


def topk_mode_recovery(
    adata,
    pairs,
    *,
    modes=("induced", "repressed"),
    k: int = 100,
    q_cut: float = 0.05,
    mode_key: str = "response_direction",
    groupby: str = "broad_cluster",
    rank_key: str = "waldStat_1",
    q_key: str = "qval_1",
    gene_key: str = "gene",
) -> pd.DataFrame:
    """Top-k recovery of response modes between contexts (TESTING ``calc_topk_mode_recovery``).

    For each mode and ``(source, target)`` context pair: the fraction of the top
    ``k`` source genes of that mode (by ``rank_key``, ``q_key`` < ``q_cut``) whose
    gene has the same mode with ``q_key`` < ``q_cut`` in the target context. Modes
    come from ``obs[mode_key]``; when that column is absent they are the response
    direction of ``obs[groupby]`` (``MODE_DIRECTIONS``). ``adata`` is a
    gene-trajectory AnnData (or its ``obs`` DataFrame) with a ``cell_type`` column.
    """
    try:
        from rpy2.robjects import StrVector, r
    except ImportError as exc:
        raise ImportError("topk_mode_recovery requires optional dependency 'rpy2'.") from exc
    from chronopack._backends.r._bridge import source_r_resource
    from chronopack.pl._gene_trajectory_plots import _r_to_pandas

    obs = adata.copy() if isinstance(adata, pd.DataFrame) else adata.obs.copy()
    if mode_key not in obs:
        obs[mode_key] = obs[groupby].astype(str).map(MODE_DIRECTIONS)
    table = pd.DataFrame({
        "cell_type": obs["cell_type"].astype(str).to_numpy(),
        mode_key: obs[mode_key].astype(object).to_numpy(),
        rank_key: pd.to_numeric(obs[rank_key], errors="coerce").to_numpy(),
        q_key: pd.to_numeric(obs[q_key], errors="coerce").to_numpy(),
        gene_key: obs[gene_key].astype(str).to_numpy(),
    })
    combos = [(str(m), str(s), str(t)) for m in modes for s, t in pairs]
    source_r_resource("dimred_clustering_helpers.R")
    with TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "pc_df.csv"
        table.to_csv(path, index=False)
        value = r["chronopack_topk_mode_recovery_csv"](
            str(path),
            sources=StrVector([c[1] for c in combos]),
            targets=StrVector([c[2] for c in combos]),
            modes=StrVector([c[0] for c in combos]),
            mode_col=mode_key,
            k=int(k),
            q_cut=float(q_cut),
            rank_col=rank_key,
            q_col=q_key,
            gene_col=gene_key,
        )
    out = _r_to_pandas(value)
    out.insert(3, "k", int(k))
    return out


__all__ = [
    "MODE_DIRECTIONS",
    "topk_mode_recovery",
]
