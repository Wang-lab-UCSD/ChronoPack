"""Gene-dynamics fitting orchestration."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy.sparse as sp

from chronopack._backends.r._bridge import source_r_resource
from chronopack._core import record_provenance, validate_adata
from chronopack._core._keys import GENE_DYNAMICS_KEY
from chronopack._core._validation import ChronoPackError


@dataclass
class GeneDynamicsResult:
    """Lightweight container for gene-dynamics backend outputs."""

    backend: str
    key_added: str
    metadata: dict[str, Any] = field(default_factory=dict)
    association: pd.DataFrame | None = None
    significant: pd.DataFrame | None = None
    selected_genes: list[str] = field(default_factory=list)
    sidecar_paths: dict[str, str] = field(default_factory=dict)
    backend_result: Any | None = None


def prepare_tradeseq_inputs(
    adata,
    output_dir,
    *,
    full_adata=None,
    layer="counts",
    time_key="time",
    gene_subset=None,
    tag="",
):
    """Write the current tradeSeq CSV inputs from an AnnData object."""
    validate_adata(adata, time_key=time_key, layer=layer)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if gene_subset is None:
        gene_subset = adata.var_names.astype(str).to_list()
    else:
        gene_subset = [str(g) for g in gene_subset]
        missing = sorted(set(gene_subset) - set(adata.var_names.astype(str)))
        if missing:
            raise ChronoPackError(f"gene_subset contains genes absent from adata.var_names: {missing[:5]}")
        adata = adata[:, gene_subset].copy()

    matrix = adata.layers[layer]
    cells = adata.obs_names.astype(str).to_numpy()
    genes = adata.var_names.astype(str).to_numpy()
    matrix_gc = matrix.T if matrix.shape[0] == len(cells) else matrix

    if sp.issparse(matrix_gc):
        counts_df = pd.DataFrame.sparse.from_spmatrix(matrix_gc, index=genes, columns=cells)
    else:
        counts_df = pd.DataFrame(np.asarray(matrix_gc), index=genes, columns=cells)

    time = pd.to_numeric(adata.obs[time_key].astype(str), errors="coerce")
    if time.isna().any():
        raise ChronoPackError(f"adata.obs[{time_key!r}] contains non-numeric values.")
    time_df = pd.DataFrame({"cell_id": cells, "time": time.to_numpy()})

    suffix = f".{tag}" if tag else ""
    paths = {
        "counts": output_dir / f"counts_genes_x_cells{suffix}.csv.gz",
        "cell_times": output_dir / f"cell_times{suffix}.csv.gz",
        "subset_genes": output_dir / f"subset_genes{suffix}.csv",
    }
    counts_df.to_csv(paths["counts"])
    time_df.to_csv(paths["cell_times"], index=False)
    pd.Series(genes, name="gene").to_csv(paths["subset_genes"], index=False)

    if tag == "full":
        lib_source = full_adata if full_adata is not None else adata
        if "total_counts" not in lib_source.obs:
            raise ChronoPackError("Writing full tradeSeq inputs requires obs['total_counts'].")
        lib_obs = lib_source.obs.loc[cells]
        lib_df = pd.DataFrame(
            {
                "cell_id": cells,
                "libsize_full": lib_obs["total_counts"].to_numpy(),
            }
        )
        paths["cell_libsize_full"] = output_dir / "cell_libsize_full.csv"
        lib_df.to_csv(paths["cell_libsize_full"], index=False)

    return {key: str(path) for key, path in paths.items()}


def _read_tradeseq_outputs(result, input_dir):
    input_dir = Path(input_dir)
    paths = {
        "association": str(result.rx2("association_path")[0]),
        "significant": str(result.rx2("significant_path")[0]),
        "sce": str(result.rx2("sce_path")[0]),
        "gamList": str(result.rx2("gamlist_path")[0]),
    }
    association = pd.read_csv(paths["association"], index_col=0)
    significant = pd.read_csv(paths["significant"], index_col=0)
    selected_genes = [str(x) for x in list(result.rx2("genes_keep"))]
    return association, significant, selected_genes, paths


def _run_tradeseq_backend(
    *,
    input_dir,
    n_points,
    workers,
    screen_nknots,
    final_nknots,
    p_threshold,
    outputs,
    legacy_outputs,
    tag,
    verbose,
):
    try:
        from rpy2.robjects import StrVector, r
    except ImportError as exc:
        raise ImportError("tradeSeq backend execution requires the optional 'rpy2' dependency.") from exc

    source_r_resource("tradeseq_helpers.R")
    r["load_chronopack_tradeseq_libraries"]()
    run = r["run_tradeseq_single_lineage"]
    return run(
        input_dir=str(input_dir),
        n_points=int(n_points),
        workers=int(workers),
        screen_nknots=int(screen_nknots),
        final_nknots=int(final_nknots),
        p_threshold=float(p_threshold),
        output_association=str(outputs["association"]),
        output_significant=str(outputs["significant"]),
        output_sce=str(outputs["sce"]),
        output_gamlist=str(outputs["gamList"]),
        association_aliases=StrVector(legacy_outputs.get("association", [])),
        significant_aliases=StrVector(legacy_outputs.get("significant", [])),
        sce_aliases=StrVector(legacy_outputs.get("sce", [])),
        gamlist_aliases=StrVector(legacy_outputs.get("gamList", [])),
        tag=str(tag),
        verbose=bool(verbose),
    )


def _write_association_to_var(adata, association, selected_genes, key_added):
    selected = set(map(str, selected_genes))
    for col in association.columns:
        out_col = f"{key_added}_{col}"
        adata.var[out_col] = np.nan
        common = adata.var_names.intersection(association.index.astype(str))
        if len(common):
            vals = association.loc[common, col]
            adata.var.loc[common, out_col] = pd.to_numeric(vals, errors="coerce").to_numpy()
    sig_col = f"{key_added}_significant"
    adata.var[sig_col] = adata.var_names.astype(str).isin(selected)


def fit_gene_dynamics(
    adata=None,
    *,
    time_key="time",
    layer="counts",
    backend="tradeseq",
    key_added=GENE_DYNAMICS_KEY,
    input_dir=None,
    output_dir=None,
    full_adata=None,
    prepare_inputs=False,
    run_backend=True,
    n_points=9,
    workers=4,
    screen_nknots=4,
    final_nknots=8,
    p_threshold=0.05,
    outputs=None,
    legacy_outputs=None,
    tag="",
    verbose=True,
    copy=False,
    return_model=True,
):
    """Fit gene dynamics with the tradeSeq backend.

    Native tradeSeq objects are retained as sidecar RDS files. Scalar
    association-test outputs are written to ``adata.var`` when ``adata`` is
    supplied.
    """
    if backend != "tradeseq":
        raise ValueError("Only backend='tradeseq' is currently implemented.")

    target = adata.copy() if (copy and adata is not None) else adata
    if target is not None:
        validate_adata(target, time_key=time_key, layer=layer)

    if outputs is None:
        outputs = {
            "association": "asso_df.1trajectory.nknots4.csv",
            "significant": "sig_df.1trajectory.nknots4.csv",
            "sce": "sce.1trajectory.nknots8.rds",
            "gamList": "gamList.1trajectory.nknots8.rds",
        }
    if legacy_outputs is None:
        legacy_outputs = {
            "association": [],
            "significant": [],
            "sce": [],
            "gamList": [],
        }

    if output_dir is not None and input_dir is None:
        input_dir = output_dir
    if input_dir is None:
        raise ValueError("input_dir or output_dir must be provided.")
    input_dir = Path(input_dir)

    prepared_paths = {}
    if prepare_inputs:
        if target is None:
            raise ValueError("prepare_inputs=True requires adata.")
        prepared_paths.update(
            prepare_tradeseq_inputs(
                target,
                input_dir,
                full_adata=full_adata,
                layer=layer,
                time_key=time_key,
                tag="",
            )
        )
        full_source = full_adata if full_adata is not None else target
        prepared_paths.update(
            {
                f"full_{key}": value
                for key, value in prepare_tradeseq_inputs(
                    full_source,
                    input_dir,
                    full_adata=full_source,
                    layer=layer,
                    time_key=time_key,
                    tag="full",
                ).items()
            }
        )

    metadata = {
        "time_key": time_key,
        "layer": layer,
        "input_dir": str(input_dir),
        "n_points": n_points,
        "workers": workers,
        "screen_nknots": screen_nknots,
        "final_nknots": final_nknots,
        "p_threshold": p_threshold,
        "outputs": outputs,
        "legacy_outputs": legacy_outputs,
        "prepared_paths": prepared_paths,
    }

    association = None
    significant = None
    selected_genes = []
    sidecar_paths = {}
    backend_result = None

    if run_backend:
        backend_result = _run_tradeseq_backend(
            input_dir=input_dir,
            n_points=n_points,
            workers=workers,
            screen_nknots=screen_nknots,
            final_nknots=final_nknots,
            p_threshold=p_threshold,
            outputs=outputs,
            legacy_outputs=legacy_outputs,
            tag=tag,
            verbose=verbose,
        )
        association, significant, selected_genes, sidecar_paths = _read_tradeseq_outputs(
            backend_result,
            input_dir,
        )
        if target is not None:
            _write_association_to_var(target, association, selected_genes, key_added)

    if target is not None:
        result_keys = {"var": []}
        if association is not None:
            result_keys["var"] = [f"{key_added}_{col}" for col in association.columns] + [
                f"{key_added}_significant"
            ]
        if sidecar_paths:
            result_keys["sidecar"] = sidecar_paths
        record_provenance(
            target,
            key_added,
            analysis="gene_dynamics",
            backend=backend,
            parameters=metadata,
            result_keys=result_keys,
        )

    result = GeneDynamicsResult(
        backend=backend,
        key_added=key_added,
        metadata=metadata,
        association=association,
        significant=significant,
        selected_genes=selected_genes,
        sidecar_paths=sidecar_paths,
        backend_result=backend_result,
    )

    if copy and target is not None:
        return (target, result) if return_model else target
    return result if return_model else None


__all__ = ["GeneDynamicsResult", "fit_gene_dynamics", "prepare_tradeseq_inputs"]
