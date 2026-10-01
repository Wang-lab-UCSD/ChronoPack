#!/usr/bin/env python
"""Compare the outputs of two pipeline runs, stage by stage.

Usage: compare_runs.py STAGE NEW_DATA_DIR REF_DATA_DIR

Each DATA_DIR holds one folder per experiment (a run's ``data/`` folder, or the
run root for runs made before the analyses/ layout). STAGE is
PreprocessDatasets, GetTrajectories, PrepareForTradeSeq (AnnData outputs,
file by file), TradeSeqInputs (parsed tradeSeq input tables) or
GeneClustering (cluster assignments; REF_DATA_DIR is then the reference
gene-clustering output folder). Layers and X are compared in row chunks so
large dense files stay within memory. TradeSeqFits outputs are compared by
compare_tradeseq_outputs.R.
"""

from __future__ import annotations

import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from anndata.io import read_elem

FILES = {
    "PreprocessDatasets": {
        "k562-tnfa.hg38.rep": ["adata_full.1e4-variant.chk0.h5ad", "adata_hvg.1e4-variant.chk0.h5ad"],
        "thp1-k562-lps.hg38.rep": ["adata_full.chk0.h5ad", "adata_hvg.chk0.h5ad"],
    },
    "GetTrajectories": {
        "k562-tnfa.hg38.rep": ["adata_full.1e4-variant.chk1.h5ad", "adata_hvg.1e4-variant.chk1.h5ad"],
        "thp1-k562-lps.hg38.rep": ["adata_full.chk1.h5ad", "adata_hvg.chk1.h5ad"],
    },
    "PrepareForTradeSeq": {
        "k562-tnfa.hg38.rep": ["adata_prefiltered.1e4-variant.fixed.h5ad"],
        "thp1-k562-lps.hg38.rep": ["adata_k562_prefiltered.fixed.h5ad", "adata_thp1_prefiltered.fixed.h5ad"],
    },
}
# Keys renamed by the Sep 21 notebook cleanup (new name -> maintained r7 name).
RENAMED = {"ot_umap": "X_ot_umap", "S": "W_terminal_cosine_topk"}
CLUSTER_COLUMNS = ["cluster_broad_induced", "cluster_fine_induced", "cluster_broad_repressed",
                   "cluster_fine_repressed", "broad_cluster"]


def max_abs_diff(a, b):
    if a.dtype.kind in "OUSb" or b.dtype.kind in "OUSb":
        return 0.0 if np.array_equal(np.asarray(a, dtype=str), np.asarray(b, dtype=str)) else np.inf
    a = np.asarray(a, dtype=float); b = np.asarray(b, dtype=float)
    both_nan = np.isnan(a) & np.isnan(b)
    d = np.abs(a - b)
    d[both_nan] = 0
    return float(np.nanmax(d)) if d.size else 0.0


def matrix_diff(ga, gb, rows_b=None):
    """Max |a - b| for a dense dataset or CSR/CSC group, row-chunked."""
    def dense_rows(g, idx):
        if isinstance(g, h5py.Dataset):
            return g[np.sort(idx)] if idx is not None else g[:]
        m = read_elem(g)
        m = m[idx] if idx is not None else m
        return m.toarray() if hasattr(m, "toarray") else np.asarray(m)
    if isinstance(ga, h5py.Dataset) and isinstance(gb, h5py.Dataset) and rows_b is None:
        if ga.shape != gb.shape:
            return f"shape {ga.shape} vs {gb.shape}"
        out = 0.0
        for s in range(0, ga.shape[0], 512):
            out = max(out, max_abs_diff(ga[s:s + 512], gb[s:s + 512]))
        return out
    a = dense_rows(ga, None)
    b = dense_rows(gb, rows_b)
    if a.shape != b.shape:
        return f"shape {a.shape} vs {b.shape}"
    return max_abs_diff(a, b)


def compare_file(new, ref):
    lines = []
    with h5py.File(new, "r") as hn, h5py.File(ref, "r") as hr:
        on, orf = read_elem(hn["obs"]), read_elem(hr["obs"])
        vn, vr = read_elem(hn["var"]), read_elem(hr["var"])
        same_cells = on.index.equals(orf.index)
        same_genes = vn.index.equals(vr.index)
        lines.append(f"    cells {len(on)} vs {len(orf)} (same order {same_cells}) | genes {len(vn)} vs {len(vr)} (same order {same_genes})")
        if not same_cells:
            lines.append(f"      cells only in new {len(on.index.difference(orf.index))}, only in ref {len(orf.index.difference(on.index))}")
        common = on.index.intersection(orf.index)
        diffs = {c: max_abs_diff(on.loc[common, c].to_numpy(), orf.loc[common, c].to_numpy())
                 for c in on.columns.intersection(orf.columns)}
        bad = {c: d for c, d in diffs.items() if d != 0}
        lines.append(f"    obs: {len(diffs)} shared columns, {len(diffs) - len(bad)} identical; differing: {bad}")
        extra = sorted(set(on.columns) ^ set(orf.columns))
        if extra:
            lines.append(f"      obs columns only in new: {sorted(set(on.columns) - set(orf.columns))} | only in ref: {sorted(set(orf.columns) - set(on.columns))}")
        vdiff = {c: max_abs_diff(vn[c].to_numpy(), vr.loc[vn.index, c].to_numpy())
                 for c in vn.columns.intersection(vr.columns)} if same_genes else {}
        vbad = {c: d for c, d in vdiff.items() if d != 0}
        lines.append(f"    var: {len(vdiff)} shared columns, differing: {vbad}")
        if same_cells and same_genes:
            for key in ["X"] + [f"layers/{k}" for k in hn["layers"]]:
                if key in hr:
                    lines.append(f"    {key}: max |diff| {matrix_diff(hn[key], hr[key])}")
                else:
                    lines.append(f"    {key}: not in ref")
        for grp in ["obsm", "obsp"]:
            if grp not in hn:
                continue
            for k in hn[grp]:
                rk = k if k in hr.get(grp, {}) else RENAMED.get(k)
                if rk is None or rk not in hr[grp]:
                    lines.append(f"    {grp}/{k}: not in ref")
                    continue
                d = matrix_diff(hn[grp][k], hr[grp][rk]) if same_cells else "cells differ"
                lines.append(f"    {grp}/{k}{' (ref ' + rk + ')' if rk != k else ''}: max |diff| {d}")
    return lines


TRADESEQ_INPUT_DIRS = [
    "k562-tnfa.hg38.rep/intmd/sc-K562-1e4_to_tradeSeq.fixed",
    "thp1-k562-lps.hg38.rep/intmd/sc-K562_to_tradeSeq.fixed",
    "thp1-k562-lps.hg38.rep/intmd/sc-THP1_to_tradeSeq.fixed",
]


def compare_tradeseq_inputs(new_root, ref_root):
    """Parsed-value comparison; counts are integers in the new writer, floats before."""
    for d in TRADESEQ_INPUT_DIRS:
        print(f"\n== {d}")
        new_dir, ref_dir = new_root / d, ref_root / d
        # ``*.full.*`` copies were unused and are no longer written (2026-09-30).
        names = sorted(p.name for p in ref_dir.iterdir() if p.suffix in {".csv", ".gz"} and ".full." not in p.name)
        for name in names:
            new, ref = new_dir / name, ref_dir / name
            if not new.exists():
                print(f"    {name}: MISSING in new run")
                continue
            a = pd.read_csv(new, index_col=0)
            b = pd.read_csv(ref, index_col=0)
            same_index = a.index.equals(b.index)
            same_cols = a.columns.equals(b.columns)
            vals = (np.array_equal(a.to_numpy(float), b.to_numpy(float), equal_nan=True)
                    if same_index and same_cols and a.shape == b.shape and all(a.dtypes.map(lambda t: t.kind in "iuf"))
                    else a.equals(b))
            print(f"    {name:36s} shape {a.shape} vs {b.shape} | same rows {same_index} | same cols {same_cols} | values equal {vals}")
        extra = sorted(p.name for p in new_dir.iterdir() if p.name not in names and p.suffix in {".csv", ".gz"})
        if extra:
            print(f"    only in new run: {extra}")


def _trajectory_obs(folder):
    """Per-trajectory table: obs of gene_trajectory_clusters.h5ad, or pc_df.new.csv for older runs."""
    h5ad = folder / "gene_trajectory_clusters.h5ad"
    if h5ad.exists():
        with h5py.File(h5ad, "r") as handle:
            return read_elem(handle["obs"])
    return pd.read_csv(folder / "pc_df.new.csv", index_col=0)


def compare_gene_clustering(new_dir, ref_dir):
    new = _trajectory_obs(new_dir)
    ref = _trajectory_obs(ref_dir)
    print(f"rows by cell_type: new {new['cell_type'].value_counts().to_dict()} | ref {ref['cell_type'].value_counts().to_dict()}")
    common = new.index.intersection(ref.index)
    print(f"rows: new {len(new)}, ref {len(ref)}, common {len(common)}, same set {set(new.index) == set(ref.index)}")
    new, ref = new.loc[common], ref.loc[common]
    for col in CLUSTER_COLUMNS:
        n_diff = int((new[col].astype(str) != ref[col].astype(str)).sum())
        print(f"    {col:26s} identical {n_diff == 0} (n_diff {n_diff})")
    numeric = [c for c in new.columns.intersection(ref.columns) if pd.api.types.is_numeric_dtype(ref[c])]
    diffs = {c: max_abs_diff(new[c].to_numpy(), ref[c].to_numpy()) for c in numeric}
    print(f"    numeric columns: {len(diffs)}; max |diff| > 1e-10: {dict((c, d) for c, d in diffs.items() if d > 1e-10)}")


def main():
    stage, new_root, ref_root = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    if stage == "TradeSeqInputs":
        compare_tradeseq_inputs(new_root, ref_root)
        return
    if stage == "GeneClustering":
        compare_gene_clustering(new_root, ref_root)
        return
    for exp, files in FILES[stage].items():
        for name in files:
            new = new_root / exp / "output" / name
            ref = ref_root / exp / "output" / name
            print(f"\n== {exp}/{name}")
            if not new.exists() or not ref.exists():
                print(f"    missing: new {new.exists()} ref {ref.exists()}")
                continue
            print("\n".join(compare_file(new, ref)), flush=True)


if __name__ == "__main__":
    main()
