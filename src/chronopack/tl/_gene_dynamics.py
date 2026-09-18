"""Gene-dynamics public wrappers and backend boundaries."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from chronopack._core import record_provenance, validate_adata


@dataclass
class GeneDynamicsResult:
    """Lightweight container for gene-dynamics backend outputs."""

    backend: str
    key_added: str
    metadata: dict[str, Any] = field(default_factory=dict)
    backend_model: Any | None = None


def fit_gene_dynamics(
    adata,
    *,
    time_key="time",
    layer="counts",
    backend="tradeseq",
    key_added="gene_dynamics",
    copy=False,
    return_model=False,
    **backend_kwargs,
):
    """Validate and record a gene-dynamics analysis request.

    The full tradeSeq fitting implementation remains in the R backend/notebooks
    during this staged migration. This wrapper establishes the public API,
    mutation convention, and provenance contract without changing scientific
    behavior.
    """
    target = adata.copy() if copy else adata
    validate_adata(target, time_key=time_key, layer=layer)
    record_provenance(
        target,
        key_added,
        analysis="gene_dynamics",
        backend=backend,
        parameters={"time_key": time_key, "layer": layer, **backend_kwargs},
    )
    result = GeneDynamicsResult(
        backend=backend,
        key_added=key_added,
        metadata={"time_key": time_key, "layer": layer, **backend_kwargs},
    )
    if copy:
        return (target, result) if return_model else target
    return result if return_model else None


__all__ = ["GeneDynamicsResult", "fit_gene_dynamics"]
