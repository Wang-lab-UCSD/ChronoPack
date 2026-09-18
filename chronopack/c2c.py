"""Cell-cell communication helpers for ChronoPack analyses."""

from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Optional, Sequence, Union

import numpy as np
import pandas as pd
import scipy.sparse as sp


def _load_matplotlib():
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize, TwoSlopeNorm

    return plt, Normalize, TwoSlopeNorm


def time_to_minute(values: pd.Series, base_interval_min: float = 15.0) -> pd.Series:
    """Convert integer-like time tags to elapsed minutes."""
    vals = pd.to_numeric(pd.Series(values).astype(str), errors="coerce")
    return (vals - 1) * float(base_interval_min)


def to_dense(matrix):
    """Return a dense numpy array from dense or sparse matrix-like input."""
    return matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)


def load_secreted_cellchat_pairs(cellchat_path: Union[str, Path]) -> pd.DataFrame:
    """Load CellChatDB interactions and keep secreted signaling pairs."""
    df = pd.read_csv(cellchat_path)
    if "annotation" in df.columns:
        df = df[df["annotation"] == "Secreted Signaling"].copy()
    return df


def split_by_cell_type(adata, cell_type_key: str = "cell_type", sender: str = "THP1", receiver: str = "K562"):
    """Split an AnnData object into sender and receiver cell-type subsets."""
    sender_adata = adata[adata.obs[cell_type_key] == sender, :]
    receiver_adata = adata[adata.obs[cell_type_key] == receiver, :]
    return sender_adata, receiver_adata


def load_tradeseq_grid_exports(export_dir: Union[str, Path], time_max_min: float = 165.0) -> dict:
    """Load tradeSeq grid export tables and normalize their time columns."""
    export_dir = Path(export_dir).expanduser()
    fits = pd.read_csv(export_dir / "fits_on_grid.csv")
    if "time" in fits.columns:
        fits = fits.rename(columns={"time": "t"})
    if "t" in fits.columns:
        fits["t"] = pd.to_numeric(fits["t"], errors="coerce")
        if fits["t"].max(skipna=True) <= 1.5:
            fits["t"] = fits["t"] * float(time_max_min)
        else:
            fits["t"] = (fits["t"] / fits["t"].max(skipna=True)) * float(time_max_min)

    wmeans_path = export_dir / "obs_binned_means.csv"
    wmeans = pd.read_csv(wmeans_path) if wmeans_path.exists() else None
    if wmeans is not None:
        if "time_mid" in wmeans.columns:
            wmeans["time_mid"] = pd.to_numeric(wmeans["time_mid"], errors="coerce")
            if wmeans["time_mid"].max(skipna=True) <= 1.5:
                wmeans["time_mid"] = wmeans["time_mid"] * float(time_max_min)
            else:
                wmeans["time_mid"] = (wmeans["time_mid"] / wmeans["time_mid"].max(skipna=True)) * float(time_max_min)
            wmeans["pt_mid"] = wmeans["time_mid"]

    meta_path = export_dir / "meta.csv"
    meta = pd.read_csv(meta_path) if meta_path.exists() else None

    return {"fits": fits, "wmeans": wmeans, "meta": meta}


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


def all_expressed(genes: Sequence[str], det: pd.Series, min_detect: float = 0.0) -> bool:
    """Return True when every gene is present and detected above threshold."""
    if not genes:
        return False
    return all(gene in det.index and float(det[gene]) >= min_detect for gene in genes)


def filter_expressed_secreted_lr_catalog(
    interactions: pd.DataFrame,
    sender_adata,
    receiver_adata,
    sender_layer: Optional[str] = "counts",
    receiver_layer: Optional[str] = "counts",
    min_detect: float = 0.0,
) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    """Keep LR pairs whose ligand and receptor components are expressed."""
    sender_det = detection_series(sender_adata, layer=sender_layer)
    receiver_det = detection_series(receiver_adata, layer=receiver_layer)

    ligand_components = interactions["ligand"].apply(split_genes)
    receptor_components = interactions["receptor"].apply(split_genes)
    ligand_ok = ligand_components.apply(lambda genes: all_expressed(genes, sender_det, min_detect))
    receptor_ok = receptor_components.apply(lambda genes: all_expressed(genes, receiver_det, min_detect))

    keep = ligand_ok & receptor_ok
    out = interactions.loc[keep].copy()
    out["lig_components"] = ligand_components.loc[keep]
    out["rec_components"] = receptor_components.loc[keep]
    return out, sender_det, receiver_det


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


def score_lr_pairs(
    pair_series: dict,
    association_df: pd.DataFrame,
    receiver_detection: pd.Series,
    causal_window_min: tuple[float, float] = (0.0, 45.0),
    association_p_col: str = "pvalue_2",
    p_pseudocount: float = 1e-3,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Score LR pairs and aggregate them into ChronoPack ligand scores."""
    rows = []
    for (ligand, receptor), (activity_raw, _) in pair_series.items():
        if ligand not in association_df.index:
            continue
        activity = activity_raw - activity_raw.iloc[0]
        mask = (activity.index >= causal_window_min[0]) & (activity.index <= causal_window_min[1])
        receptor_genes = receptor.split("_")
        rec_det = float(np.nanmean([receiver_detection.get(gene, np.nan) for gene in receptor_genes]))
        rows.append(
            {
                "ligand": ligand,
                "receptor": receptor,
                "pair": f"{ligand}->{receptor}",
                "avg": activity.loc[mask].mean() if mask.any() else np.nan,
                "rec_det": rec_det,
                "logp": -np.log10(float(association_df.loc[ligand, association_p_col]) + p_pseudocount),
                "n_window_points": int(mask.sum()),
            }
        )

    pair_scores = pd.DataFrame(rows).sort_values(["logp", "avg"], ascending=[False, False]).reset_index(drop=True)
    pair_scores["pair_rank"] = np.arange(1, len(pair_scores) + 1)

    current_ligand_rank = pair_scores.drop_duplicates("ligand", keep="first").copy()
    current_ligand_rank["current_ligand_rank"] = np.arange(1, len(current_ligand_rank) + 1)
    ligand_scores = (
        pair_scores.groupby("ligand", as_index=False)
        .agg(
            ligand_significance=("logp", "max"),
            max_lr_activity=("avg", "max"),
            best_pair=("pair", "first"),
            best_pair_rank=("pair_rank", "first"),
            n_receptor_pairs=("pair", "size"),
        )
        .merge(current_ligand_rank[["ligand", "current_ligand_rank"]], on="ligand", how="left")
    )
    ligand_scores["z_ligand_significance"] = zscore(ligand_scores["ligand_significance"])
    ligand_scores["z_max_lr_activity"] = zscore(ligand_scores["max_lr_activity"])
    ligand_scores["composite_score"] = ligand_scores["z_ligand_significance"] + ligand_scores["z_max_lr_activity"]
    ligand_scores = ligand_scores.sort_values("composite_score", ascending=False).reset_index(drop=True)
    ligand_scores["zscore_composite_rank"] = np.arange(1, len(ligand_scores) + 1)
    return pair_scores, ligand_scores


def plot_avg_vs_logp(
    plot_df: pd.DataFrame,
    x_col: str = "avg",
    y_col: str = "logp",
    pair_col: str = "pair",
    color_col: Optional[str] = None,
    color_diverging: bool = False,
    color_center: float = 0.0,
    p_thresh: float = 0.1,
    x_thresh: Optional[float] = None,
    max_labels: Optional[int] = 20,
    figsize=(4.8, 4.0),
    dpi: int = 300,
    cmap: str = "viridis",
    title: str = "Ligand-receptor pairs",
    xlabel: str = "early LR activity",
    ylabel: str = "-log10(p)",
):
    """Scatter LR activity against ligand temporal significance."""
    plt, Normalize, TwoSlopeNorm = _load_matplotlib()
    df = plot_df.replace([np.inf, -np.inf], np.nan).dropna(subset=[x_col, y_col, pair_col]).copy()
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    sig_y = df[y_col].to_numpy(float) > -np.log10(p_thresh)
    ax.axhline(-np.log10(p_thresh), color="0.35", lw=1.5, ls="--", zorder=0)
    if x_thresh is not None:
        ax.axvline(x_thresh, color="0.35", lw=1.5, ls="--", zorder=0)

    if color_col is None:
        ax.scatter(df[x_col], df[y_col], s=np.where(sig_y, 90, 65), color="#4C78A8", edgecolors="0.25", linewidths=0.6)
    else:
        values = df[color_col].to_numpy(float)
        if color_diverging:
            vmax = np.nanmax(np.abs(values))
            norm = TwoSlopeNorm(vmin=-vmax, vcenter=color_center, vmax=vmax)
            cmap = "coolwarm"
        else:
            norm = Normalize(vmin=np.nanmin(values), vmax=np.nanmax(values))
        scatter = ax.scatter(df[x_col], df[y_col], s=np.where(sig_y, 90, 65), c=values, cmap=cmap, norm=norm, edgecolors="0.25", linewidths=0.6)
        fig.colorbar(scatter, ax=ax, pad=0.012, fraction=0.045).ax.set_title(color_col)

    label_df = df.sort_values(y_col, ascending=False).head(max_labels or len(df))
    for _, row in label_df.iterrows():
        ax.annotate(str(row[pair_col]).replace("->", "->"), (row[x_col], row[y_col]), xytext=(0, 5), textcoords="offset points", ha="center", fontsize=7)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(alpha=0.18)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    return fig, ax


def plot_pair_overlay(
    pair_series: dict,
    pair=("TNF", "TNFRSF1A"),
    norm_mode: str = "zscore",
    nfkb_onset_time: Optional[float] = None,
    tmin_min: float = 0.0,
    tmax_min: float = 165.0,
    xtick_step: float = 30.0,
    figsize=(5, 2),
    dpi: int = 300,
    title: Optional[str] = None,
):
    """Overlay one LR activity curve with its receiver response curve."""
    plt, _, _ = _load_matplotlib()
    if pair not in pair_series:
        raise KeyError(f"{pair} not found in pair_series")

    def normalize(values):
        values = np.asarray(values, dtype=float)
        if norm_mode == "zscore":
            sd = np.nanstd(values)
            return np.zeros_like(values) if not np.isfinite(sd) or sd == 0 else (values - np.nanmean(values)) / sd
        if norm_mode == "minmax":
            rng = np.nanmax(values) - np.nanmin(values)
            return np.zeros_like(values) if not np.isfinite(rng) or rng == 0 else (values - np.nanmin(values)) / rng
        if norm_mode in {"none", "raw"}:
            return values
        raise ValueError("norm_mode must be 'zscore', 'minmax', or 'none'")

    def scale_index(index):
        arr = np.asarray(index, dtype=float)
        lo, hi = np.nanmin(arr), np.nanmax(arr)
        return np.full_like(arr, tmin_min, dtype=float) if hi == lo else (arr - lo) / (hi - lo) * (tmax_min - tmin_min) + tmin_min

    activity, response = pair_series[pair]
    x_activity = scale_index(activity.index)
    x_response = scale_index(response.index)
    y_activity = normalize(activity.to_numpy())
    y_response = normalize(response.to_numpy())

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ligand, receptor = pair
    ax.plot(x_activity, y_activity, linewidth=2.2, label=f"{ligand}->{receptor}")
    response_line, = ax.plot(x_response, y_response, linewidth=2.2, linestyle="--", label="receiver response")

    if nfkb_onset_time is not None:
        order = np.argsort(x_response)
        y_onset = np.interp(float(nfkb_onset_time), x_response[order], y_response[order])
        ax.plot(float(nfkb_onset_time), y_onset, marker="D", markersize=7, markerfacecolor="white", markeredgecolor=response_line.get_color(), markeredgewidth=1.8, linestyle="None", zorder=10)

    ax.set_xlim(tmin_min, tmax_min)
    ax.set_xticks(np.arange(tmin_min, tmax_min + 1, xtick_step))
    ax.set_xlabel("time (min)")
    ax.set_ylabel({"zscore": "z-score", "minmax": "min-max normalized"}.get(norm_mode, "value"))
    ax.set_title(title or f"{ligand}->{receptor}")
    ax.legend(frameon=False, fontsize=7)
    ax.grid(alpha=0.15)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    return fig, ax


def plot_genes_over_time(
    adata,
    gene_list: Sequence[str],
    time_key: str = "time",
    cluster_key: str = "cell_type",
    layer: Optional[str] = None,
    ncols: int = 3,
    figsize_per_panel=(4.0, 3.2),
):
    """Plot raw expression trajectories for a short list of genes."""
    plt, _, _ = _load_matplotlib()
    time_values = pd.to_numeric(pd.Series(adata.obs[time_key]).astype(str), errors="coerce")
    cluster_values = pd.Series(adata.obs[cluster_key]).astype(str)
    mask_valid = time_values.notna()
    clusters = sorted(cluster_values[mask_valid].unique())
    colors = dict(zip(clusters, plt.rcParams["axes.prop_cycle"].by_key().get("color", [])))
    nrows = math.ceil(len(gene_list) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows))
    axes = np.asarray(axes).ravel()
    matrix = adata.layers[layer] if layer is not None and layer in adata.layers else adata.X
    for ax, gene in zip(axes, gene_list):
        if gene not in adata.var_names:
            ax.text(0.5, 0.5, f"{gene}\nnot found", ha="center", va="center")
            ax.axis("off")
            continue
        expr = pd.Series(np.asarray(to_dense(matrix[:, adata.var_names.get_loc(gene)])).ravel())[mask_valid.to_numpy()].reset_index(drop=True)
        t = time_values[mask_valid].reset_index(drop=True)
        c = cluster_values[mask_valid].reset_index(drop=True)
        for cluster in clusters:
            sub = pd.DataFrame({"time": t[c == cluster], "expr": expr[c == cluster]})
            stats = sub.groupby("time").agg(mean=("expr", "mean"), n=("expr", "size"), sd=("expr", "std")).reset_index()
            stats["se"] = stats["sd"].fillna(0.0) / np.sqrt(stats["n"].clip(lower=1))
            ax.plot(stats["time"], stats["mean"], label=cluster, color=colors.get(cluster))
            ax.fill_between(stats["time"], stats["mean"] - 1.96 * stats["se"], stats["mean"] + 1.96 * stats["se"], alpha=0.2, color=colors.get(cluster))
        ax.set_title(gene)
        ax.set_xlabel("time")
        ax.set_ylabel("expression")
    for ax in axes[len(gene_list) :]:
        ax.axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=max(1, len(labels)), frameon=False)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    return fig, axes
