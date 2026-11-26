"""Composition-robust library sizes (TMM) and a check of library-composition effects.

Per-cell total-count normalisation assumes the genes a cell expresses keep a constant share of its library. When a
block of genes is strongly induced it takes a growing share of every library, so every other gene's normalised
expression falls without any change in its own transcript numbers. TMM (trimmed mean of M-values; Robinson & Oshlack
2010, edgeR's default) estimates the library-size correction from the bulk of genes after trimming the most changed
ones, so an induced block does not move it.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

from chronopack._core import record_provenance, validate_adata
from chronopack._core._keys import CHRONOPACK_UNS_KEY
from chronopack._core._validation import ChronoPackError

TMM_FACTOR_KEY = "tmm_factor"
COMPOSITION_LEVELS = ((0.3, "strong"), (0.1, "moderate"), (0.0, "negligible"))


# --------------------------------------------------------------------------------------------------------------------
# pseudobulk and TMM
# --------------------------------------------------------------------------------------------------------------------
def _group_codes(obs, groupby):
    keys = [groupby] if isinstance(groupby, str) else list(groupby)
    missing = [k for k in keys if k not in obs]
    if missing:
        raise ChronoPackError(f"adata.obs has no column(s) {missing}.")
    labels = obs[keys].astype(str).agg("|".join, axis=1) if len(keys) > 1 else obs[keys[0]]
    codes, uniques = pd.factorize(labels, sort=True)
    return codes, pd.Index(uniques), keys


def _group_sums(X, codes, n_groups):
    """Sum the rows (cells) of X within each group: dense array groups x genes."""
    G = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(n_groups, len(codes)))
    S = G @ (X if sp.issparse(X) else np.asarray(X, dtype=float))
    return np.asarray(S.toarray() if sp.issparse(S) else S, dtype=float)


def _tmm_one(obs, ref, n_o, n_r, logratio_trim, sum_trim):
    """edgeR's ``.calcFactorTMM`` (weighted, no A cut-off) for one sample against the reference sample."""
    with np.errstate(divide="ignore", invalid="ignore"):
        log_r = np.log2((obs / n_o) / (ref / n_r))
        abs_e = (np.log2(obs / n_o) + np.log2(ref / n_r)) / 2
        v = (n_o - obs) / n_o / obs + (n_r - ref) / n_r / ref
    fin = np.isfinite(log_r) & np.isfinite(abs_e)
    log_r, abs_e, v = log_r[fin], abs_e[fin], v[fin]
    if log_r.size == 0 or np.max(np.abs(log_r)) < 1e-6:
        return 1.0
    n = log_r.size
    lo_l, lo_s = np.floor(n * logratio_trim) + 1, np.floor(n * sum_trim) + 1
    hi_l, hi_s = n + 1 - lo_l, n + 1 - lo_s
    rl, ra = stats.rankdata(log_r), stats.rankdata(abs_e)
    keep = (rl >= lo_l) & (rl <= hi_l) & (ra >= lo_s) & (ra <= hi_s)
    f = np.nansum(log_r[keep] / v[keep]) / np.nansum(1 / v[keep])
    return float(2 ** (0.0 if np.isnan(f) else f))


def tmm_norm_factors(counts, lib_size=None, *, logratio_trim=0.3, sum_trim=0.05, ref_sample=None):
    """TMM normalisation factors of pseudobulk samples, as ``edgeR::calcNormFactors(DGEList(t(counts),
    lib.size=lib_size), method="TMM")``.

    ``counts`` is samples x genes; ``lib_size`` (default: the row sums) the samples' library sizes, e.g. summed
    per-cell totals over all genes. The reference is ``ref_sample`` (a row position) or, as in edgeR, the sample
    whose upper-quartile CPM is closest to the mean upper quartile. Factors are scaled to a geometric mean of 1; the
    effective library size of a sample is ``lib_size`` x its factor.
    """
    counts = np.asarray(counts, dtype=float)
    lib = counts.sum(axis=1) if lib_size is None else np.asarray(lib_size, dtype=float)
    counts = counts[:, (counts > 0).any(axis=0)]
    if (lib <= 0).any():
        raise ChronoPackError("Every pseudobulk sample needs a positive total count.")
    if ref_sample is None:
        f75 = np.quantile(counts / lib[:, None], 0.75, axis=1)
        ref_sample = (int(np.argmax(np.sqrt(counts).sum(axis=1))) if np.median(f75) < 1e-20
                      else int(np.argmin(np.abs(f75 - f75.mean()))))
    f = np.array([_tmm_one(counts[i], counts[ref_sample], lib[i], lib[ref_sample], logratio_trim, sum_trim)
                  for i in range(len(counts))])
    return f / np.exp(np.mean(np.log(f)))


def _check_not_prescreened(adata):
    entry = (adata.uns.get(CHRONOPACK_UNS_KEY) or {}).get("time_association_prescreen")
    if entry is not None and adata.n_vars < int(entry["result_keys"]["n_genes_screened"]):
        raise ChronoPackError(
            "adata holds only the genes kept by prescreen_time_association, which are enriched for changing genes; "
            "TMM needs the bulk of genes. Run cp.tl.tmm_factors on the unscreened AnnData first (it writes "
            "obs['tmm_factor'], which the subset keeps)."
        )


def tmm_factors(
    adata,
    *,
    groupby="time",
    layer="counts",
    library_size_key="total_counts",
    key_added=TMM_FACTOR_KEY,
    logratio_trim=0.3,
    sum_trim=0.05,
):
    """Composition-robust library-size multipliers from TMM on pseudobulk samples.

    Cells are summed per ``groupby`` sample (a time column, or e.g. ``["library", "time"]`` when libraries span the
    time points), TMM factors are computed between those pseudobulks (``tmm_norm_factors``), and every cell gets its
    sample's multiplier in ``adata.obs[key_added]``. A cell's composition-robust library size is
    ``obs[library_size_key] * obs[key_added]``; ``prepare_tradeseq_inputs(normalization="tmm")`` writes that as the
    tradeSeq library size.

    The multiplier is the sample's TMM factor computed with the samples' summed ``library_size_key`` as their
    library sizes (``total_counts`` is computed before gene filtering, so UMIs in genes filtered out of ``adata``
    count too), as ``edgeR::calcNormFactors(DGEList(counts, lib.size=...))``. With ``library_size_key=None`` the
    library size is the cell's total over ``adata``'s genes.

    Run it on the cells of one cell type or context (the factors compare samples of the same cells), and on the
    bulk of the genes: after technical and low-count filtering is fine, but not on the genes kept by
    ``prescreen_time_association``, which are enriched for changing genes.

    Assumption: most genes do not change per unit of RNA. If stimulated cells make more mRNA overall, changes become
    relative to the bulk of genes rather than absolute.

    Returns a DataFrame per sample with ``n_cells``, ``counts`` (over ``adata``'s genes), ``library_size`` and
    ``multiplier``.
    """
    validate_adata(adata, layer=layer)
    _check_not_prescreened(adata)
    X = adata.layers[layer] if layer is not None else adata.X
    codes, samples, keys = _group_codes(adata.obs, groupby)
    S = _group_sums(X, codes, len(samples))
    counts = S.sum(axis=1)
    if library_size_key is not None:
        if library_size_key not in adata.obs:
            raise ChronoPackError(f"adata.obs has no {library_size_key!r}; pass library_size_key=None to use the "
                                  "totals over adata's genes.")
        lib = np.bincount(codes, weights=adata.obs[library_size_key].to_numpy(float), minlength=len(samples))
    else:
        lib = counts
    mult = tmm_norm_factors(S, lib, logratio_trim=logratio_trim, sum_trim=sum_trim)
    table = pd.DataFrame(
        {"n_cells": np.bincount(codes, minlength=len(samples)), "counts": counts, "library_size": lib,
         "multiplier": mult},
        index=pd.Index(samples.astype(str), name="|".join(keys)),
    )
    adata.obs[key_added] = mult[codes]
    record_provenance(
        adata,
        key_added,
        analysis="tmm_factors",
        parameters={"groupby": keys, "layer": layer, "library_size_key": library_size_key,
                    "logratio_trim": logratio_trim, "sum_trim": sum_trim, "n_genes": int(adata.n_vars)},
        result_keys={"obs": key_added, "multiplier": table["multiplier"].to_dict()},
    )
    return table


# --------------------------------------------------------------------------------------------------------------------
# composition check
# --------------------------------------------------------------------------------------------------------------------
@dataclass
class CompositionQCResult:
    """Output of :func:`composition_qc`.

    ``by_time``: one row per (group,) time point; ``genes``: per-gene pseudobulk log2 fold changes against the
    baseline under both normalisations, with the composition flags; ``max_shift_log2`` and ``level`` per group.
    """

    by_time: pd.DataFrame
    genes: pd.DataFrame
    induced_genes: dict = field(default_factory=dict)
    max_shift_log2: dict = field(default_factory=dict)
    level: dict = field(default_factory=dict)

    def summary(self):
        return pd.DataFrame({"max_shift_log2": self.max_shift_log2, "level": self.level,
                             "n_induced_genes": {k: len(v) for k, v in self.induced_genes.items()}})


def _level(shift):
    return next(name for cut, name in COMPOSITION_LEVELS if abs(shift) >= cut)


def _composition_one(X, obs, time_key, library_size_key, baseline, induced_genes, genes, min_log2fc,
                     min_counts, gene_log2fc, prior_cpm, logratio_trim, sum_trim):
    time = pd.to_numeric(obs[time_key].astype(str), errors="coerce")
    if time.isna().any():
        raise ChronoPackError(f"adata.obs[{time_key!r}] contains non-numeric values.")
    codes, times = pd.factorize(time.to_numpy(float), sort=True)
    times = np.asarray(times, dtype=float)
    if times.size < 2:
        raise ChronoPackError("composition_qc needs at least two time points.")
    b = 0 if baseline is None else int(np.flatnonzero(np.isclose(times, float(baseline)))[0])
    S = _group_sums(X, codes, len(times))
    n_cells = np.bincount(codes, minlength=len(times))
    counts = S.sum(axis=1)
    lib = (np.bincount(codes, weights=obs[library_size_key].to_numpy(float), minlength=len(times))
           if library_size_key is not None and library_size_key in obs else counts)
    mult = tmm_norm_factors(S, lib, logratio_trim=logratio_trim, sum_trim=sum_trim)
    mult = mult / mult[b]
    shift = np.log2(mult)                                   # apparent log2 change of an unchanged gene, total vs TMM
    trend = pd.Series(shift).rolling(3, center=True, min_periods=1).mean().to_numpy()

    cpm_total = S / lib[:, None] * 1e6
    cpm_tmm = S / (lib * mult)[:, None] * 1e6
    lfc_total = np.log2((cpm_total + prior_cpm) / (cpm_total[b] + prior_cpm))
    lfc_tmm = np.log2((cpm_tmm + prior_cpm) / (cpm_tmm[b] + prior_cpm))
    expressed = S.sum(axis=0) >= (min_counts if min_counts is not None else 20 * len(times))

    if induced_genes is None:
        induced = expressed & (np.delete(lfc_tmm, b, axis=0).max(axis=0) >= min_log2fc)
    else:
        induced = genes.isin(pd.Index(induced_genes).astype(str))
    share = S[:, induced].sum(axis=1) / counts
    by_time = pd.DataFrame({
        "time": times, "n_cells": n_cells, "umis_per_cell": lib / n_cells,
        "induced_share": share,
        "induced_share_bound_log2": np.log2((1 - share) / (1 - share[b])),
        "tmm_multiplier": mult,
        "composition_shift_log2": shift,
        "composition_trend_log2": trend,
    })

    others = np.ones(len(times), bool)
    others[b] = False
    g = pd.DataFrame({
        "min_log2fc_total": lfc_total[others].min(axis=0), "min_log2fc_tmm": lfc_tmm[others].min(axis=0),
        "max_log2fc_total": lfc_total[others].max(axis=0), "max_log2fc_tmm": lfc_tmm[others].max(axis=0),
        "expressed": expressed, "induced": induced,
    }, index=genes)
    g["composition_explainable"] = (expressed & (g.min_log2fc_total <= -gene_log2fc)
                                    & (g.min_log2fc_tmm > -gene_log2fc))
    g["composition_masked"] = (expressed & (g.max_log2fc_tmm >= gene_log2fc) & (g.max_log2fc_total < gene_log2fc))
    return by_time, g, list(genes[induced])


def composition_qc(
    adata,
    *,
    time_key="time",
    layer="counts",
    groupby=None,
    library_size_key="total_counts",
    baseline=None,
    induced_genes=None,
    min_log2fc=1.0,
    min_counts=None,
    gene_log2fc=0.5,
    prior_cpm=1.0,
    warn_threshold=0.1,
    logratio_trim=0.3,
    sum_trim=0.05,
    key_added="composition",
):
    """Check how much library composition moves total-count-normalised expression over time.

    Cells are summed per time point (per ``groupby`` value too, e.g. a cell-type column; each group is checked on
    its own). For every time point:

    - ``composition_shift_log2``: the apparent log2 change, against ``baseline`` (default: the first time point), of
      a gene whose transcript numbers do not change, under per-cell total-count normalisation (what
      ``prepare_tradeseq_inputs`` writes by default). It is log2 of the TMM multiplier relative to the baseline:
      negative when an induced block takes over the libraries. This is what ``normalization="tmm"`` would remove.
    - ``induced_share``: the share of UMIs from induced genes (pseudobulk log2 fold change >= ``min_log2fc`` at
      some time point under TMM, >= ``min_counts`` UMIs in total; or ``induced_genes``), and
      ``induced_share_bound_log2``, the drop that share alone causes in every other gene,
      log2((1 - share_t) / (1 - share_baseline)). The two estimates usually agree; the TMM one also includes genes
      below the fold-change cut and repressed genes.
    - ``composition_trend_log2``: the shift averaged over each time point and its neighbours (3-point running
      mean). Composition driven by an induced programme changes smoothly; a single time point off the trend is a
      sample-level difference (e.g. one tag's libraries), not the response. ``level`` uses the trend.
    - ``umis_per_cell``.

    Per gene (expressed: >= ``min_counts`` UMIs in total, default 20 per time point), the largest pseudobulk fall
    and rise against the baseline under both normalisations, and two flags (read them together with the
    dynamic-gene calls: for noisy genes near the cut-off a small shift flips them either way):

    - ``composition_explainable``: falls by >= ``gene_log2fc`` log2 under total-count normalisation but not under
      TMM. A "repressed" call for such a gene can be the composition shift alone.
    - ``composition_masked``: rises by >= ``gene_log2fc`` under TMM but not under total-count normalisation.

    Without ``groupby`` the flags are written to ``adata.var[f"{key_added}_explainable"]`` and
    ``[f"{key_added}_masked"]``. ``level`` per group is "negligible" (largest |trend| < 0.1 log2), "moderate"
    (< 0.3) or "strong"; a warning is raised when it reaches ``warn_threshold``. ``max_shift_log2`` is that
    largest trend value.

    These are pseudobulk checks; they need no tradeSeq fits and can run before them, to choose the normalisation.
    Run on the bulk of genes (not after ``prescreen_time_association``).
    """
    validate_adata(adata, time_key=time_key, layer=layer)
    _check_not_prescreened(adata)
    X = adata.layers[layer] if layer is not None else adata.X
    if sp.issparse(X):
        X = X.tocsr()
    genes = adata.var_names.astype(str)
    groups = [(None, np.ones(adata.n_obs, bool))]
    if groupby is not None:
        lab = adata.obs[groupby].astype(str).to_numpy()
        groups = [(gname, lab == gname) for gname in pd.unique(lab)]

    rows, gene_tables, induced, max_shift, level = [], [], {}, {}, {}
    for gname, mask in groups:
        bt, gt, ind = _composition_one(
            X[mask], adata.obs.loc[mask], time_key, library_size_key, baseline, induced_genes, genes, min_log2fc,
            min_counts, gene_log2fc, prior_cpm, logratio_trim, sum_trim)
        name = "all" if gname is None else gname
        if gname is not None:
            bt.insert(0, "group", gname)
            gt.insert(0, "group", gname)
        rows.append(bt)
        gene_tables.append(gt)
        induced[name] = ind
        s = bt["composition_trend_log2"].to_numpy()
        max_shift[name] = float(s[np.argmax(np.abs(s))])
        level[name] = _level(max_shift[name])
        if abs(max_shift[name]) >= warn_threshold:
            label = "" if gname is None else f" in {gname}"
            warnings.warn(
                f"Library composition shifts total-count-normalised expression by up to {max_shift[name]:+.2f} "
                f"log2{label} ({level[name]}); {int(gt['composition_explainable'].sum())} genes fall by >= "
                f"{gene_log2fc} log2 only under total-count normalisation. Consider "
                "fit_gene_dynamics(..., normalization='tmm').",
                stacklevel=2,
            )
    by_time = pd.concat(rows, ignore_index=True)
    gene_table = pd.concat(gene_tables) if groupby is None else pd.concat(gene_tables).reset_index(names="gene")
    if groupby is None:
        adata.var[f"{key_added}_explainable"] = gene_table["composition_explainable"].to_numpy()
        adata.var[f"{key_added}_masked"] = gene_table["composition_masked"].to_numpy()
    record_provenance(
        adata,
        key_added,
        analysis="composition_qc",
        parameters={"time_key": time_key, "layer": layer, "groupby": groupby, "library_size_key": library_size_key,
                    "baseline": baseline, "min_log2fc": min_log2fc, "min_counts": min_counts,
                    "gene_log2fc": gene_log2fc, "prior_cpm": prior_cpm},
        result_keys={"max_shift_log2": max_shift, "level": level,
                     "var": [] if groupby is not None else [f"{key_added}_explainable", f"{key_added}_masked"]},
    )
    return CompositionQCResult(by_time=by_time, genes=gene_table, induced_genes=induced, max_shift_log2=max_shift,
                               level=level)


__all__ = ["CompositionQCResult", "composition_qc", "tmm_factors", "tmm_norm_factors"]
