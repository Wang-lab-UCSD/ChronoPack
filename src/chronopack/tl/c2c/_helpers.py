"""Cell-cell communication helpers used by `cp.tl.c2c.infer` and the C2C notebook."""

from __future__ import annotations

import re
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd


def time_to_minute(values: pd.Series, base_interval_min: float = 15.0) -> pd.Series:
    """Convert integer-like time tags to elapsed minutes."""
    vals = pd.to_numeric(pd.Series(values).astype(str), errors="coerce")
    return (vals - 1) * float(base_interval_min)


def to_dense(matrix):
    """Return a dense numpy array from dense or sparse matrix-like input."""
    return matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)


def split_by_cell_type(adata, cell_type_key: str = "cell_type", sender: str = "THP1", receiver: str = "K562"):
    """Split an AnnData object into sender and receiver cell-type subsets."""
    sender_adata = adata[adata.obs[cell_type_key] == sender, :]
    receiver_adata = adata[adata.obs[cell_type_key] == receiver, :]
    return sender_adata, receiver_adata


def detection_series(adata, layer: Optional[str] = None) -> pd.Series:
    """Calculate fraction of cells with nonzero expression for each gene."""
    X = adata.layers[layer] if layer is not None and layer in adata.layers else adata.X
    det = (to_dense(X) > 0).mean(axis=0)
    return pd.Series(np.asarray(det).ravel(), index=adata.var_names.astype(str))


def split_genes(value: str) -> list[str]:
    """Split CellChat ligand/receptor strings into component genes."""
    if pd.isna(value):
        return []
    tokens = [tok for tok in re.split(r"[|,/;+\s]+|_", str(value).strip()) if tok]
    out = []
    last_full = None
    for token in tokens:
        match = re.fullmatch(r"R(\d+)", token)
        if match and last_full is not None:
            base = re.sub(r"\d+$", "", last_full)
            out.append(base + match.group(1))
        else:
            out.append(token)
            last_full = token
    return out


def restrict_ligands_to_association(lr_pairs: pd.DataFrame, association_df: pd.DataFrame) -> pd.DataFrame:
    """Keep LR pairs whose ligand has a temporal association row."""
    ligands = set(lr_pairs["ligand"].astype(str)) & set(association_df.index.astype(str))
    return lr_pairs[lr_pairs["ligand"].astype(str).isin(ligands)].copy()


def moving_average(values, window: int = 3):
    """Centered moving average with edge correction."""
    x = np.asarray(values, dtype=float)
    if window <= 1:
        return x
    kernel = np.ones(window) / window
    y = np.convolve(x, kernel, mode="same")
    for i in range(window // 2):
        y[i] = np.mean(x[: i + 1])
        y[-i - 1] = np.mean(x[-(i + 1) :])
    return y


def time_series(
    adata,
    genes: Sequence[str],
    time_key: str = "time",
    layer: str = "log1p",
    agg: str = "mean",
    min_detect: float = 0.0,
    prefer_fits: bool = True,
    lineage: Union[int, str] = 1,
    use: str = "fit",
    return_ci: bool = False,
    fits_df: Optional[pd.DataFrame] = None,
    wmeans_df: Optional[pd.DataFrame] = None,
    base_interval_min: float = 15.0,
    smooth: bool = False,
    smooth_window: int = 1,
):
    """Build a gene-set time series from tradeSeq fits or raw AnnData expression."""
    lineage_label = f"L{int(lineage)}" if not (isinstance(lineage, str) and lineage.startswith("L")) else lineage
    genes = [gene for gene in genes if gene in adata.var_names]
    if not genes:
        return None

    X = adata.layers[layer] if layer in adata.layers else adata.X
    expr = pd.DataFrame(
        to_dense(X[:, [adata.var_names.get_loc(gene) for gene in genes]]),
        index=adata.obs_names,
        columns=genes,
    )
    expr[time_key] = time_to_minute(adata.obs[time_key], base_interval_min=base_interval_min).to_numpy()

    detect = (expr[genes] > 0).mean()
    keep = [gene for gene in genes if detect[gene] >= min_detect]
    if not keep:
        return None

    if prefer_fits and fits_df is not None and "gene" in fits_df.columns:
        if use == "fit":
            sub = fits_df[(fits_df["gene"].isin(keep)) & (fits_df["lineage"] == lineage_label)].copy()
            if not sub.empty:
                pivot_mu = sub.pivot_table(index="t", columns="gene", values="mu_log1p", aggfunc="mean")
                mu = pivot_mu.mean(axis=1, skipna=True).sort_index()
                if return_ci:
                    pivot_lo = sub.pivot_table(index="t", columns="gene", values="lo_log1p", aggfunc="mean")
                    pivot_hi = sub.pivot_table(index="t", columns="gene", values="hi_log1p", aggfunc="mean")
                    return pd.DataFrame(
                        {
                            "mu": mu,
                            "lo": pivot_lo.mean(axis=1, skipna=True),
                            "hi": pivot_hi.mean(axis=1, skipna=True),
                        }
                    ).sort_index()
                return mu
        elif use == "weighted" and wmeans_df is not None:
            required = {"gene", "lineage", "pt_mid", "mean_norm_log1p"}
            if required.issubset(wmeans_df.columns):
                sub = wmeans_df[(wmeans_df["gene"].isin(keep)) & (wmeans_df["lineage"] == lineage_label)].copy()
                if not sub.empty:
                    pivot = sub.pivot_table(index="pt_mid", columns="gene", values="mean_norm_log1p", aggfunc="mean")
                    return pivot.mean(axis=1, skipna=True).sort_index()

    grouped = expr.groupby(time_key)[keep]
    if agg == "mean":
        ts = grouped.mean().mean(axis=1)
    else:
        ts = grouped.sum().sum(axis=1)
    values = ts.sort_index().to_numpy(dtype=float)
    if smooth:
        values = moving_average(values, window=smooth_window)
    return pd.Series(values, index=np.sort(expr[time_key].unique()))


def module_series(adata, genes, time_key="time", layer="log1p", **kwargs):
    """Return a module-average time series."""
    return time_series(adata, genes, time_key=time_key, layer=layer, agg="mean", **kwargs)


def lr_activity(
    sender_adata,
    receiver_adata,
    ligand: str,
    receptor: str,
    time_key: str = "time",
    layer: str = "log1p",
    prefer_fits: bool = False,
    fits_df: Optional[pd.DataFrame] = None,
    wmeans_df: Optional[pd.DataFrame] = None,
    sender_lineage: Union[int, str] = 2,
    receiver_lineage: Union[int, str] = 1,
) -> Optional[pd.Series]:
    """Compute LR activity as ligand expression times mean receptor availability."""
    ligand_series = time_series(
        sender_adata,
        [ligand],
        time_key=time_key,
        layer=layer,
        lineage=sender_lineage,
        prefer_fits=prefer_fits,
        fits_df=fits_df,
        wmeans_df=wmeans_df,
    )
    if ligand_series is None:
        return None

    receptor_genes = receptor.split("_")
    receptor_series = [
        time_series(
            receiver_adata,
            [gene],
            time_key=time_key,
            layer=layer,
            lineage=receiver_lineage,
            prefer_fits=prefer_fits,
            fits_df=fits_df,
            wmeans_df=wmeans_df,
        )
        for gene in receptor_genes
    ]
    if any(series is None for series in receptor_series):
        return None

    common = ligand_series.index
    for series in receptor_series:
        common = common.intersection(series.index)
    if len(common) == 0:
        return None

    ligand_series = ligand_series.loc[common]
    receptor_values = [np.clip(series.loc[common].to_numpy(dtype=float), 0, None) for series in receptor_series]
    if len(receptor_values) == 1:
        receptor_availability = receptor_values[0].mean()
    else:
        receptor_availability = np.sqrt(np.prod(np.vstack(receptor_values), axis=0)).mean()
    activity = np.sqrt(np.clip(ligand_series.to_numpy(dtype=float), 0, None) * receptor_availability)
    return pd.Series(activity, index=common)


def build_pair_series(
    lr_pairs: pd.DataFrame,
    sender_adata,
    receiver_adata,
    response_genes: Sequence[str],
    time_key: str = "time",
    layer: str = "log1p",
    prefer_fits: bool = False,
    fits_df: Optional[pd.DataFrame] = None,
    wmeans_df: Optional[pd.DataFrame] = None,
) -> dict:
    """Build aligned LR activity and receiver response series for all pairs."""
    response = module_series(
        receiver_adata,
        response_genes,
        time_key=time_key,
        layer=layer,
        lineage=1,
        prefer_fits=prefer_fits,
        fits_df=fits_df,
        wmeans_df=wmeans_df,
    )
    out = {}
    if response is None:
        return out
    for ligand, receptor in zip(lr_pairs["ligand"], lr_pairs["receptor"]):
        activity = lr_activity(
            sender_adata,
            receiver_adata,
            ligand,
            receptor,
            time_key=time_key,
            layer=layer,
            prefer_fits=prefer_fits,
            fits_df=fits_df,
            wmeans_df=wmeans_df,
        )
        if activity is None:
            continue
        common = activity.index.intersection(response.index)
        out[(ligand, receptor)] = (activity.loc[common], response.loc[common])
    return out


def zscore(values: pd.Series) -> pd.Series:
    """Population z-score with stable zero-variance handling."""
    sd = values.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return values * 0.0
    return (values - values.mean()) / sd


