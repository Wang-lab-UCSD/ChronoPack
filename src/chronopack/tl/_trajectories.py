"""Cell-state trajectory and fate analysis helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import sparse

from chronopack._core import record_provenance, validate_adata
from chronopack._core._keys import (
    CELL_STATE_TRAJECTORIES_KEY,
    CELL_STATE_TRAJECTORY_CLUSTER_KEY,
    OT_FORCE_DIRECTED_KEY,
    OT_GRAPH_KEY,
    OT_GROWTH_KEY,
    OT_UMAP_KEY,
    TERMINAL_CELL_DIST_INDEX_KEY,
    TERMINAL_CELL_DIST_KEY,
    TERMINAL_COSINE_GRAPH_KEY,
    TERMINAL_FATE_PROBS_KEY,
)


def _scanpy():
    import scanpy as sc

    return sc


def _pyplot():
    import matplotlib.pyplot as plt

    return plt


def _sklearn_neighbors_and_normalize():
    from sklearn.neighbors import NearestNeighbors
    from sklearn.preprocessing import normalize

    return NearestNeighbors, normalize


class CellStateTrajectoryResult:
    """Sidecar outputs from WOT-derived cell-state trajectory analysis.

    Rectangular time-to-time WOT couplings and sampled/model state are not
    AnnData-axis-aligned, so ChronoPack keeps them in this result object.
    Square cell-cell graphs may also be attached to ``adata.obsp`` by the
    workflow for Scanpy-compatible downstream use.

    ``transitions_long`` (one row per cell pair) and ``T_full`` (every coupling
    in one cell x cell matrix) hold all n_t x n_t+1 couplings, which is large for
    big time points. Either may be given as a zero-argument function; it is then
    built the first time it is accessed.

    ``growth`` holds WOT's growth weight of each cell (its coupling row sum
    relative to its map's mean, see :func:`tmap_growth`) as a Series indexed like
    ``adata.obs_names``; cells at the last time point are NaN.
    """

    _FIELDS = ("tmaps", "transitions_long", "tmats", "T_full", "T_pruned", "T_sym", "growth")

    def __init__(self, tmaps, transitions_long, tmats, T_full, T_pruned, T_sym, growth=None):
        self.tmaps = tmaps
        self._transitions_long = transitions_long
        self.tmats = tmats
        self._T_full = T_full
        self.T_pruned = T_pruned
        self.T_sym = T_sym
        self.growth = growth

    @property
    def transitions_long(self) -> pd.DataFrame:
        if callable(self._transitions_long):
            self._transitions_long = self._transitions_long()
        return self._transitions_long

    @property
    def T_full(self) -> sp.spmatrix:
        if callable(self._T_full):
            self._T_full = self._T_full()
        return self._T_full

    def __repr__(self) -> str:
        return f"CellStateTrajectoryResult(n_tmaps={len(self.tmaps)}, T_pruned={getattr(self.T_pruned, 'shape', None)})"

    def as_dict(self) -> dict:
        """Return the legacy dictionary representation (builds every lazy field)."""
        return {key: getattr(self, key) for key in self._FIELDS}

    def __getitem__(self, key):
        if key not in self._FIELDS:
            raise KeyError(key)
        return getattr(self, key)

    def keys(self):
        return list(self._FIELDS)

    def items(self):
        return self.as_dict().items()

    def values(self):
        return self.as_dict().values()

    def __iter__(self) -> Iterator[str]:
        return iter(self._FIELDS)


def _read_tmap_with_growth(fpath):
    """:func:`read_tmap` plus each source cell's coupling row sum divided by the map's mean row sum."""
    fpath = Path(fpath)
    ad = _scanpy().read_h5ad(fpath)
    t_src, t_tgt = map(float, fpath.stem.split("_")[1:])
    T = sparse.csr_matrix(ad.X)
    row_sums = np.asarray(T.sum(axis=1)).ravel()
    inv_row = 1.0 / np.maximum(row_sums, 1e-12)
    T_norm = T.multiply(inv_row[:, None])
    mean = row_sums.mean() if row_sums.size else 0.0
    growth = row_sums / mean if mean > 0 else row_sums
    return (t_src, t_tgt, T_norm, ad.obs_names.to_numpy(), ad.var_names.to_numpy()), growth


def read_tmap(fpath):
    """Read a WOT transport map and return a row-normalized sparse matrix."""
    return _read_tmap_with_growth(fpath)[0]


def _renormalize_rows(T):
    row_sums = np.asarray(T.sum(axis=1)).flatten()
    nz = row_sums != 0
    row_sums[nz] = 1.0 / row_sums[nz]
    return sparse.diags(row_sums) @ T


def _transition_entry(T, src_ids, tgt_ids, t_next, make_dense=False):
    return {
        "csr": T,
        "dense": T.toarray() if make_dense else None,
        "src_ids": src_ids,
        "tgt_ids": tgt_ids,
        "t_next": t_next,
        "src_lookup": {cid: i for i, cid in enumerate(src_ids)},
    }


def build_transition_matrices(df, make_dense=False, renorm_rows=True):
    """Build per-time transport CSR matrices from long-form WOT couplings."""
    matrices = {}
    for (t_src, t_tgt), g in df.groupby(["t_source", "t_target"]):
        src_ids = g["cell_i"].unique()
        tgt_ids = g["cell_j"].unique()
        src_lookup = {cid: i for i, cid in enumerate(src_ids)}
        tgt_lookup = {cid: i for i, cid in enumerate(tgt_ids)}
        rows = g["cell_i"].map(src_lookup).to_numpy()
        cols = g["cell_j"].map(tgt_lookup).to_numpy()
        data = g["p_ij"].to_numpy()
        T = sparse.csr_matrix((data, (rows, cols)), shape=(len(src_ids), len(tgt_ids)))
        if renorm_rows:
            T = _renormalize_rows(T)
        matrices[t_src] = _transition_entry(T, src_ids, tgt_ids, t_tgt, make_dense)
    print(f"Built {len(matrices)} transport matrices.")
    return matrices


def transition_matrices_from_tmaps(tmaps, make_dense=False, renorm_rows=True):
    """Build per-time transport CSR matrices directly from :func:`read_tmap` outputs.

    Gives the same matrices as ``build_transition_matrices`` on the long table of
    the same maps, without building that table (one row per cell pair).
    """
    matrices = {}
    for t_src, t_tgt, T, src_ids, tgt_ids in sorted(tmaps, key=lambda m: (m[0], m[1])):
        coo = T.tocoo()
        # Cells are ordered as they first appear among the nonzero couplings, as in the long table.
        row_order = pd.unique(coo.row)
        col_order = pd.unique(coo.col)
        row_pos = np.empty(T.shape[0], dtype=np.int64)
        row_pos[row_order] = np.arange(row_order.size)
        col_pos = np.empty(T.shape[1], dtype=np.int64)
        col_pos[col_order] = np.arange(col_order.size)
        M = sparse.csr_matrix((coo.data, (row_pos[coo.row], col_pos[coo.col])), shape=(row_order.size, col_order.size))
        if renorm_rows:
            M = _renormalize_rows(M)
        matrices[t_src] = _transition_entry(M, np.asarray(src_ids)[row_order], np.asarray(tgt_ids)[col_order], t_tgt, make_dense)
    print(f"Built {len(matrices)} transport matrices.")
    return matrices


def _transitions_long_from_tmaps(tmaps):
    tidy_frames = []
    for t_src, t_tgt, T, src_ids, tgt_ids in tmaps:
        coo = T.tocoo()
        tidy_frames.append(
            pd.DataFrame(
                {
                    "t_source": t_src,
                    "t_target": t_tgt,
                    "cell_i": src_ids[coo.row],
                    "cell_j": tgt_ids[coo.col],
                    "p_ij": coo.data,
                }
            )
        )
    return pd.concat(tidy_frames, ignore_index=True)


def read_tmaps(tmap_path):
    """Read every WOT map (``*.h5ad``) in a directory with :func:`read_tmap`."""
    return [read_tmap(fp) for fp in sorted(Path(tmap_path).glob("*.h5ad"))]


def _read_tmaps_with_growth(tmap_path):
    tmaps, parts = [], []
    for fp in sorted(Path(tmap_path).glob("*.h5ad")):
        tmap, growth = _read_tmap_with_growth(fp)
        tmaps.append(tmap)
        parts.append(pd.Series(growth, index=pd.Index(tmap[3]).astype(str)))
    growth = pd.concat(parts) if parts else pd.Series(dtype=float)
    growth.name = OT_GROWTH_KEY
    return tmaps, growth


def tmap_growth(tmap_path):
    """WOT growth weight of every source cell of the maps in a directory.

    The weight is the cell's row sum in its WOT coupling divided by the mean row
    sum of that map, so the cells of each time point average 1. It is the growth
    that WOT's ``TransportMapModel.fates`` and ``pull_back`` carry along (see the
    ``growth`` option of :func:`CalculateTerminalFatesAcrossCells`). Returns a
    Series indexed by cell; cells of the last time point have no outgoing map
    and are not included.
    """
    return _read_tmaps_with_growth(tmap_path)[1]


def coupling_from_kernel(T, growth):
    """Scale each row of a row-normalized transition matrix by the cell's growth weight.

    Returns ``diag(growth) @ T``: from the transition matrices ChronoPack keeps
    (e.g. ``result.T_full``) and ``result.growth``, this gives WOT's couplings, each
    map divided by its mean row sum. ``growth`` is ordered like the rows of ``T``.
    Cells without a growth weight (NaN, e.g. the last time point) must have empty
    rows; they stay empty. For a pruned ``T`` each row sums to the cell's growth
    weight, i.e. the kept transitions carry the cell's full coupling mass.
    """
    T = sp.csr_matrix(T)
    g = np.asarray(growth, dtype=float).ravel()
    if g.shape != (T.shape[0],):
        raise ValueError(f"growth has {g.size} values for a matrix with {T.shape[0]} rows.")
    missing = ~np.isfinite(g)
    if (missing & (np.diff(T.indptr) > 0)).any():
        raise ValueError("growth is missing (NaN/inf) for cells that have outgoing transitions.")
    return (sp.diags(np.where(missing, 0.0, g), format="csr") @ T).tocsr()

def extract_wot_couplings(tmap_path, make_dense=False):
    """Read all WOT maps from a directory into tmaps, long table, and tmats."""
    tmaps = read_tmaps(tmap_path)
    transitions_long = _transitions_long_from_tmaps(tmaps)
    tmats = transition_matrices_from_tmaps(tmaps, make_dense=make_dense)
    return tmaps, transitions_long, tmats


def run_wot_transport_workflow(
    adata,
    tmap_dir,
    ot_model_kwargs=None,
    compute_transport_maps=True,
    transition_topk=30,
    transition_mass=None,
    symmetrize="maximum",
    use_unpruned_as_pruned=False,
    attach_graph=True,
    graph_key="ot_graph",
    umap_key=OT_UMAP_KEY,
    umap_random_state=10,
    umap_kwargs=None,
):
    """Run WaddingtonOT maps, extract couplings, and attach transition graph outputs.

    Each cell keeps its largest forward transitions: at most ``transition_topk``
    (``None`` for no limit) and, with ``transition_mass`` (a fraction in (0, 1]),
    only as many as needed to reach that fraction of its transition mass. The
    pruned graph is built map by map, so the full n_t x n_t+1 couplings are only
    held as the maps themselves; ``T_full`` and ``transitions_long`` on the result
    are built when first accessed.

    The transition matrices are row-normalized; each cell's WOT growth weight
    (:func:`tmap_growth`) is kept as ``result.growth`` and, with ``attach_graph``,
    in ``adata.obs["ot_growth"]`` (NaN at the last time point).
    """
    tmap_dir = Path(tmap_dir)
    tmap_dir.mkdir(parents=True, exist_ok=True)
    if compute_transport_maps:
        import wot

        kwargs = dict(epsilon=0.05, lambda1=1, lambda2=50, growth_iters=3, day_field="time")
        if ot_model_kwargs:
            kwargs.update(ot_model_kwargs)
        ot_model = wot.ot.OTModel(adata, **kwargs)
        ot_model.compute_all_transport_maps(tmap_out=str(tmap_dir / "tmaps"))

    tmaps, growth = _read_tmaps_with_growth(tmap_dir)
    tmats = transition_matrices_from_tmaps(tmaps)
    obs_names = adata.obs_names.copy()
    growth = growth.reindex(obs_names.astype(str))
    growth.index = obs_names
    if use_unpruned_as_pruned:
        T_pruned = embed_transition_matrices(obs_names, tmats)
    else:
        T_pruned = embed_transition_matrices(obs_names, tmats, topk=transition_topk, mass=transition_mass)
    T_sym = symmetrize_graph(T_pruned, symmetrize)
    if attach_graph:
        adata.obs[OT_GROWTH_KEY] = growth.to_numpy()
        adata.obsp[graph_key] = T_sym
        adata.uns[graph_key] = {
            "connectivities_key": graph_key,
            "distances_key": graph_key,
            "params": {
                "method": "WOT_tmats",
                "top_k": transition_topk,
                "symmetric": symmetrize,
                **({"mass": transition_mass} if transition_mass is not None else {}),
            },
        }
        kwargs = dict(init_pos="spectral", copy=False, random_state=umap_random_state)
        if umap_kwargs:
            kwargs.update(umap_kwargs)
        original_umap_uns = adata.uns["umap"].copy() if "umap" in adata.uns and hasattr(adata.uns["umap"], "copy") else adata.uns.get("umap")
        _scanpy().tl.umap(adata, neighbors_key=graph_key, **kwargs)
        adata.obsm[umap_key] = adata.obsm["X_umap"].copy()
        if original_umap_uns is not None:
            adata.uns["umap"] = original_umap_uns
    return CellStateTrajectoryResult(
        tmaps=tmaps,
        transitions_long=lambda: _transitions_long_from_tmaps(tmaps),
        tmats=tmats,
        T_full=lambda: embed_transition_matrices(obs_names, tmats),
        T_pruned=T_pruned,
        T_sym=T_sym,
        growth=growth,
    )


def get_trajectories(
    adata,
    *,
    tmap_dir,
    time_key="time",
    key_added=CELL_STATE_TRAJECTORIES_KEY,
    graph_key=OT_GRAPH_KEY,
    umap_key=OT_UMAP_KEY,
    copy=False,
    **kwargs,
):
    """Compute WOT-derived cell-state trajectory outputs.

    The function follows Scanpy-style AnnData mutation semantics. It modifies
    ``adata`` in place by default and returns a ``CellStateTrajectoryResult``
    containing non-axis-aligned WOT sidecar objects. With ``copy=True``, it
    returns ``(adata_copy, result)``.
    """
    target = adata.copy() if copy else adata
    validate_adata(target, time_key=time_key)
    result = run_wot_transport_workflow(
        target,
        tmap_dir=tmap_dir,
        graph_key=graph_key,
        umap_key=umap_key,
        **kwargs,
    )
    record_provenance(
        target,
        key_added,
        analysis="cell_state_trajectories",
        backend="waddington_ot",
        parameters={
            "time_key": time_key,
            "graph_key": graph_key,
            "umap_key": umap_key,
            **kwargs,
        },
        result_keys={
            "obs": [OT_GROWTH_KEY],
            "obsp": [graph_key],
            "obsm": [umap_key],
            "sidecar": ["tmaps", "transitions_long", "tmats", "T_full", "T_pruned", "T_sym", "growth"],
        },
    )
    return (target, result) if copy else result


def ot_force_directed(
    adata,
    *,
    graph_key=OT_GRAPH_KEY,
    key_added=OT_FORCE_DIRECTED_KEY,
    layout="fr",
    weighted=True,
    init_pos=None,
    random_state=0,
    copy=False,
    **layout_kwargs,
):
    """Force-directed layout of the WOT transition graph, an alternative to the OT-UMAP.

    Lays out ``adata.obsp[graph_key]`` (the graph written by :func:`get_trajectories`)
    with :func:`scanpy.tl.draw_graph` and stores the coordinates in
    ``adata.obsm[key_added]``. Cells are linked only to their transitions at the
    neighbouring time points, so the layout spreads time along one axis and fates
    along the other.

    ``layout="fr"`` (Fruchterman-Reingold, through python-igraph) needs no
    dependency beyond Leiden clustering's; with ``weighted=True`` stronger
    transitions pull harder, and ``niter`` defaults to 1000. ``layout="fa"``
    (ForceAtlas2) needs the ``fa2-modified`` package. ``init_pos`` is an ``obsm``
    key to start from (e.g. ``"ot_umap"``); by default cells start at random
    positions drawn with ``random_state``. Other keyword arguments go to the
    layout function. Modifies ``adata`` in place, or returns a copy with
    ``copy=True``.
    """
    from importlib.util import find_spec

    sc = _scanpy()
    target = adata.copy() if copy else adata
    if graph_key not in target.obsp:
        raise KeyError(f"adata.obsp[{graph_key!r}] not found; run get_trajectories first.")
    if layout == "fa":
        if find_spec("fa2_modified") is None:
            raise ImportError("layout='fa' (ForceAtlas2) needs the fa2-modified package: pip install fa2-modified")
    elif weighted:
        layout_kwargs.setdefault("weights", "weight")
    if layout == "fr":
        layout_kwargs.setdefault("niter", 1000)

    # draw_graph turns every stored entry into an undirected edge, so a symmetric matrix would give each pair two
    # parallel edges (twice the pull); pass each pair once, with the larger of its two weights.
    G = sp.csr_matrix(target.obsp[graph_key])
    adjacency = sp.triu(G.maximum(G.T), k=1).tocsr()
    original_uns = target.uns.get("draw_graph")
    tmp_ext = "_chronopack_tmp"
    sc.tl.draw_graph(
        target,
        layout=layout,
        adjacency=adjacency,
        init_pos=init_pos if init_pos is not None else False,
        random_state=random_state,
        key_added_ext=tmp_ext,
        **layout_kwargs,
    )
    target.obsm[key_added] = target.obsm.pop(f"X_draw_graph_{tmp_ext}")
    if original_uns is None:
        target.uns.pop("draw_graph", None)
    else:
        target.uns["draw_graph"] = original_uns
    record_provenance(
        target,
        key_added,
        analysis="ot_force_directed_layout",
        backend="scanpy.tl.draw_graph",
        parameters={
            "graph_key": graph_key,
            "layout": layout,
            "weighted": bool(weighted and layout != "fa"),
            "init_pos": init_pos,
            "random_state": random_state,
            **{k: v for k, v in layout_kwargs.items() if k != "weights"},
        },
        result_keys={"obsm": [key_added]},
    )
    return target if copy else None


def prune_transition_rows(M, topk=None, mass=None):
    """Keep only each row's largest entries of a sparse matrix (values unchanged).

    ``topk`` keeps at most that many entries per row. ``mass`` (a fraction in
    (0, 1]) keeps the fewest largest entries whose sum reaches that fraction of
    the row's total; the largest entry is always kept. With both, a row stops at
    whichever limit comes first. Ties are broken by column order.
    """
    M = sp.csr_matrix(M, copy=True)
    M.sum_duplicates()
    if topk is None and mass is None:
        return M
    if topk is not None and (int(topk) != topk or topk < 1):
        raise ValueError(f"topk must be a positive integer or None, got {topk!r}.")
    if mass is not None and not 0 < mass <= 1:
        raise ValueError(f"mass must be in (0, 1] or None, got {mass!r}.")
    n_rows = M.shape[0]
    counts = np.diff(M.indptr)
    rows = np.repeat(np.arange(n_rows), counts)
    order = np.lexsort((-M.data, rows))  # by row, then by decreasing value
    keep_sorted = np.ones(M.nnz, dtype=bool)
    if topk is not None:
        rank = np.arange(M.nnz) - np.repeat(M.indptr[:-1], counts)
        keep_sorted &= rank < topk
    if mass is not None:
        vals = M.data[order]
        csum = np.concatenate([[0.0], np.cumsum(vals)])
        start = np.repeat(csum[M.indptr[:-1]], counts)
        total = np.repeat(csum[M.indptr[1:]] - csum[M.indptr[:-1]], counts)
        before = csum[:-1] - start  # row mass of the larger entries
        keep_sorted &= before < mass * total
    keep = np.zeros(M.nnz, dtype=bool)
    keep[order[keep_sorted]] = True
    indptr = np.concatenate([[0], np.cumsum(np.bincount(rows[keep], minlength=n_rows))]).astype(M.indptr.dtype)
    return sp.csr_matrix((M.data[keep], M.indices[keep], indptr), shape=M.shape)


def keep_topk(mat, k):
    """For each sparse row, keep only the k largest elements."""
    return prune_transition_rows(mat, topk=k)


def embed_transition_matrices(obs_names, tmats, *, topk=None, mass=None):
    """Place per-time transition matrices into one cell x cell matrix over ``obs_names``.

    With ``topk`` / ``mass``, each map's rows are pruned first
    (:func:`prune_transition_rows`). A cell is the source of a single map, so
    this equals pruning the rows of the full matrix.
    """
    obs_names = pd.Index(obs_names)
    n = len(obs_names)
    rows, cols, data = [], [], []
    for tdict in tmats.values():
        M = tdict["csr"]
        if topk is not None or mass is not None:
            M = prune_transition_rows(M, topk=topk, mass=mass)
        M = M.tocoo()
        rows.append(obs_names.get_indexer(np.asarray(tdict["src_ids"]))[M.row])
        cols.append(obs_names.get_indexer(np.asarray(tdict["tgt_ids"]))[M.col])
        data.append(M.data)
    if not rows:
        return sp.csr_matrix((n, n))
    return sp.csr_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))


def symmetrize_graph(T, symmetrize="maximum"):
    """Symmetrize a directed transition graph: elementwise ``"maximum"`` or ``"average"``."""
    if symmetrize == "average":
        return 0.5 * (T + T.T)
    return T.maximum(T.T)


def sparse_adjacency_from_tmats(bdata, tmats, topk=30, symmetrize="maximum", mass=None):
    """Embed per-time transition matrices into one cell-by-cell adjacency."""
    T_full = embed_transition_matrices(bdata.obs_names, tmats)
    T_pruned = embed_transition_matrices(bdata.obs_names, tmats, topk=topk, mass=mass)
    return T_full, T_pruned, symmetrize_graph(T_pruned, symmetrize)


def embed_block(A, idx, n):
    """Place square block A at rows/cols given by idx inside an n x n sparse matrix."""
    A = A.tocoo()
    return sparse.csr_matrix((A.data, (idx[A.row], idx[A.col])), shape=(n, n))


def merge_trajectory_graphs(adata, graphs):
    """Combine per-subset cell-cell graphs into one graph over ``adata.obs_names``.

    ``graphs`` is an iterable of ``(obs_names, graph)`` pairs, where each square
    sparse ``graph`` is ordered like its ``obs_names``. Used when transport maps
    are computed separately per cell type. Subsets must not overlap.
    """
    obs_index = pd.Index(adata.obs_names)
    n = adata.n_obs
    merged = None
    seen = pd.Index([])
    for names, graph in graphs:
        names = pd.Index(names)
        if graph.shape != (len(names), len(names)):
            raise ValueError(f"Graph shape {graph.shape} does not match {len(names)} obs_names.")
        missing = names.difference(obs_index)
        if len(missing):
            raise ValueError(f"{len(missing)} obs_names are not in adata.obs_names.")
        if len(seen.intersection(names)):
            raise ValueError("Cell subsets passed to merge_trajectory_graphs overlap.")
        seen = seen.append(names)
        block = embed_block(graph, obs_index.get_indexer(names), n)
        merged = block if merged is None else merged + block
    if merged is None:
        raise ValueError("merge_trajectory_graphs needs at least one graph.")
    return merged.tocsr()


def topk_csr_rows(M_csr, k=200, renorm=True, eps=1e-15):
    """Keep only top-k values per CSR row, optionally row-normalizing."""
    M = M_csr.tocsr()
    indptr, indices, data = [0], [], []
    for i in range(M.shape[0]):
        start, end = M.indptr[i], M.indptr[i + 1]
        cols = M.indices[start:end]
        vals = M.data[start:end]
        if vals.size == 0:
            indptr.append(indptr[-1])
            continue
        if vals.size > k:
            keep = np.argpartition(vals, -k)[-k:]
            order = keep[np.argsort(-vals[keep])]
        else:
            order = np.argsort(-vals)
        cols = cols[order]
        vals = vals[order]
        if renorm and vals.sum() > 0:
            vals = vals / (vals.sum() + eps)
        indices.extend(cols.tolist())
        data.extend(vals.tolist())
        indptr.append(len(indices))
    out = sp.csr_matrix(
        (np.asarray(data, dtype=np.float64), np.asarray(indices, dtype=np.int32), np.asarray(indptr, dtype=np.int64)),
        shape=M.shape,
    )
    out.eliminate_zeros()
    return out


def _topk_csr_rows(M_csr, k=100, eps=1e-12):
    return topk_csr_rows(M_csr, k=k, renorm=False, eps=eps)


def _resolve_growth(adata, growth):
    """Per-cell growth weights in ``adata.obs_names`` order from an obs key, a Series or an array."""
    if growth is None:
        return None
    if isinstance(growth, str):
        if growth not in adata.obs:
            raise KeyError(f"adata.obs[{growth!r}] not found; get_trajectories stores the WOT growth weights there.")
        growth = adata.obs[growth]
    if isinstance(growth, pd.Series):
        growth = growth.reindex(adata.obs_names)
    g = np.asarray(growth, dtype=float).ravel()
    if g.shape != (adata.n_obs,):
        raise ValueError(f"growth has {g.size} values for {adata.n_obs} cells.")
    return g


def CalculateTerminalFatesAcrossCells(adata, T, topk=None, growth=None):
    """Compute terminal-cell fate distributions for every cell.

    Each cell's distribution over the cells of the last time point is propagated
    backwards through the row-normalized forward transitions of ``T``, which equals
    the corresponding rows of ``P^k`` when ``T`` links consecutive time points.
    ``topk=None`` (default) keeps every terminal cell, so the distributions are
    exact; an integer keeps only each cell's ``topk`` largest terminal-cell
    probabilities at every step (renormalized), an approximation that sharpens
    early cells' fate probabilities.

    ``growth`` (default ``None``: every cell counts equally) weights the
    transitions by the cells' growth, as WOT's ``TransportMapModel.fates`` does: an
    obs key (``"ot_growth"``, stored by :func:`get_trajectories`), a Series indexed
    by cell, or an array in ``adata.obs_names`` order. Each cell's descendant mass
    at the last time point is ``M_i = growth_i * sum_j P_ij M_j`` (``M = 1`` at the
    last time point), and the step from cell ``i`` follows ``P_ij M_j``
    (renormalized), so a cell's distribution is that of its descendants' mass
    rather than of a single random walk. With the unpruned transition matrix and
    ``growth=result.growth`` this reproduces WOT's fates. Growth weights are only
    needed before the last time point; there they must be finite and non-negative.
    Since ``T`` is row-normalized at every step, it may also be a coupling matrix
    (pass its row sums as ``growth`` to keep them).
    """
    assert sp.issparse(T), "T must be a scipy sparse matrix."
    assert T.shape == (adata.n_obs, adata.n_obs), "T must be n_cells x n_cells."
    T = T.tocsr().astype(np.float64).copy()
    t = adata.obs["time"]
    time_levels = list(t.cat.categories) if hasattr(t, "cat") else sorted(pd.unique(t))
    tt = t.to_numpy().astype(object)
    idx_by_time_all = [np.where(tt == lvl)[0] for lvl in time_levels]
    idx_by_time = [ii for ii in idx_by_time_all if ii.size > 0]
    time_levels_nonempty = [lvl for lvl, ii in zip(time_levels, idx_by_time_all) if ii.size > 0]
    assert len(idx_by_time) >= 2, "Need at least 2 timepoints."
    idx_last = idx_by_time[-1]
    n_last = idx_last.size
    n = adata.n_obs
    print(f"Terminal timepoint: {time_levels_nonempty[-1]} | n_terminal={n_last}")
    g = _resolve_growth(adata, growth)
    if g is not None:
        before_last = np.concatenate(idx_by_time[:-1])
        bad = ~np.isfinite(g[before_last]) | (g[before_last] < 0)
        if bad.any():
            raise ValueError(f"growth must be finite and non-negative before the last time point; {int(bad.sum())} cells are not.")
    D_by_time = [None] * len(idx_by_time)
    D_by_time[-1] = sp.identity(n_last, format="csr", dtype=np.float64)
    mass_by_time = [None] * len(idx_by_time)  # descendant mass at the last time point (with growth)
    mass_by_time[-1] = np.ones(n_last)
    for ti in range(len(idx_by_time) - 2, -1, -1):
        ii = idx_by_time[ti]
        future_idx = np.concatenate(idx_by_time[ti + 1 :])
        D_future = sp.vstack([D_by_time[j] for j in range(ti + 1, len(idx_by_time))], format="csr")
        A = T[ii][:, future_idx].tocsr()
        rs = np.asarray(A.sum(axis=1)).ravel()
        inv = np.zeros_like(rs, dtype=np.float64)
        nz = rs > 0
        inv[nz] = 1.0 / (rs[nz] + 1e-15)
        P = sp.diags(inv, format="csr") @ A
        if g is not None:
            W = (P @ sp.diags(np.concatenate(mass_by_time[ti + 1 :]), format="csr")).tocsr()
            s = np.asarray(W.sum(axis=1)).ravel()
            inv_s = np.zeros_like(s)
            inv_s[s > 0] = 1.0 / s[s > 0]
            P = sp.diags(inv_s, format="csr") @ W
            mass_by_time[ti] = g[ii] * s
        D_t = (P @ D_future).tocsr()
        if topk is not None:
            D_t = topk_csr_rows(D_t, k=topk, renorm=True)
        else:
            D_t.eliminate_zeros()
        zero_rows = np.where(np.asarray(D_t.sum(axis=1)).ravel() == 0)[0]
        if zero_rows.size:
            D_t = D_t.tolil()
            for r in zero_rows:
                D_t.rows[r] = [0]
                D_t.data[r] = [1.0]
            D_t = D_t.tocsr()
        D_by_time[ti] = D_t
        print(f"Computed terminal distributions for time bin {ti + 1}/{len(idx_by_time)} (n={ii.size})")
    rows_all, cols_all, data_all = [], [], []
    for ti, ii in enumerate(idx_by_time):
        D_t = D_by_time[ti].tocoo()
        rows_all.append(ii[D_t.row])
        cols_all.append(D_t.col)
        data_all.append(D_t.data)
    D_full = sp.coo_matrix(
        (np.concatenate(data_all), (np.concatenate(rows_all), np.concatenate(cols_all))),
        shape=(n, n_last),
    ).tocsr()
    rs = np.asarray(D_full.sum(axis=1)).ravel()
    D_full = sp.diags(1.0 / np.where(rs > 0, rs, 1.0), format="csr") @ D_full
    adata.obsm[TERMINAL_CELL_DIST_KEY] = D_full.astype(np.float32)
    adata.uns[TERMINAL_CELL_DIST_INDEX_KEY] = idx_last
    record_provenance(
        adata,
        TERMINAL_CELL_DIST_KEY,
        analysis="terminal_cell_distributions",
        backend="chronopack",
        parameters={"topk": topk, "growth": growth if growth is None or isinstance(growth, str) else "per-cell values"},
        result_keys={"obsm": [TERMINAL_CELL_DIST_KEY], "uns": [TERMINAL_CELL_DIST_INDEX_KEY]},
    )
    print(f"Stored adata.obsm['{TERMINAL_CELL_DIST_KEY}'] with shape {D_full.shape} (sparse nnz={D_full.nnz})")
    return D_full


def _cosine_topk_neighbors(X, k, max_values_per_block=50_000_000):
    """Sparse matrix of each row's ``k`` largest cosine similarities (itself included).

    ``X`` is dense with L2-normalized rows; similarities are computed a block of
    rows at a time as dense products, so memory stays at a block of n values.
    """
    n = X.shape[0]
    k = min(k, n)
    step = max(1, max_values_per_block // max(n, 1))
    rows, cols, vals = [], [], []
    for start in range(0, n, step):
        sim = X[start : start + step] @ X.T
        idx = np.argpartition(-sim, k - 1, axis=1)[:, :k] if k < n else np.broadcast_to(np.arange(n), sim.shape)
        rows.append(np.repeat(np.arange(start, start + sim.shape[0]), k))
        cols.append(idx.ravel())
        vals.append(np.take_along_axis(sim, idx, axis=1).ravel())
    vals = np.concatenate(vals)
    # As sklearn's cosine distance: 1 - similarity clipped to [0, 2], then back to a similarity.
    vals = 1.0 - np.clip(1.0 - vals, 0.0, 2.0)
    return sp.csr_matrix((vals.astype(np.float32), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))


def cosine_topk_leiden_from_terminal_dist(adata, dist_key=TERMINAL_CELL_DIST_KEY, topk=100, random_state=0):
    """Build a cosine top-k graph from terminal-cell distributions.

    The distributions are mostly dense (every terminal cell is reachable from
    early cells), so they are compared as a dense float32 matrix.
    """
    assert dist_key in adata.obsm, f"Need adata.obsm['{dist_key}']"
    _, normalize = _sklearn_neighbors_and_normalize()
    X = adata.obsm[dist_key]
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    n = adata.n_obs
    Xn = normalize(X.astype(np.float32, copy=False), norm="l2", axis=1)
    S = _cosine_topk_neighbors(Xn, topk + 1)
    S.setdiag(0.0)
    S.eliminate_zeros()
    S = _topk_csr_rows(S.maximum(S.T).tocsr(), k=topk)
    adata.obsp[TERMINAL_COSINE_GRAPH_KEY] = S
    return S


def terminal_fate_trajectory_clustering(
    adata,
    T,
    cluster_key=CELL_STATE_TRAJECTORY_CLUSTER_KEY,
    resolution=0.3,
    terminal_topk=None,
    graph_topk=500,
    random_state=0,
    directed=False,
    growth=None,
):
    """Compute terminal fate distributions and Leiden trajectory clusters.

    ``terminal_topk`` is passed to :func:`CalculateTerminalFatesAcrossCells` as
    ``topk``: ``None`` (default) for exact distributions. ``growth`` (default
    ``None``) is passed on as well, to weight the transitions by WOT's growth.
    """
    CalculateTerminalFatesAcrossCells(adata, T, topk=terminal_topk, growth=growth)
    S = cosine_topk_leiden_from_terminal_dist(
        adata,
        dist_key=TERMINAL_CELL_DIST_KEY,
        topk=graph_topk,
        random_state=random_state,
    )
    _scanpy().tl.leiden(
        adata,
        adjacency=S.tocsr(),
        directed=directed,
        use_weights=True,
        key_added=cluster_key,
        resolution=resolution,
        random_state=random_state,
    )
    return S


def cluster_cells_to_trajectories(
    adata,
    S=None,
    *,
    cluster_key=CELL_STATE_TRAJECTORY_CLUSTER_KEY,
    resolution=0.3,
    terminal_topk=None,
    graph_topk=500,
    random_state=0,
    directed=False,
    growth=None,
):
    """Cluster cells into cell-state trajectories from a square cell-state graph.

    ``terminal_topk=None`` (default) uses exact terminal-cell distributions; an
    integer keeps each cell's ``terminal_topk`` largest at every step.

    ``S`` is usually the pruned transition matrix (``result.T_pruned`` from
    :func:`get_trajectories`); passing ``result.T_full`` instead uses every WOT
    transition (and, with ``growth``, gives exactly WOT's fates) while the OT
    graph and OT-UMAP stay pruned, at the cost of memory for large time points.

    ``growth=None`` (default) treats the transitions as a Markov chain, every cell
    counting equally. ``growth="ot_growth"`` (or a Series of
    ``result.growth``, e.g. concatenated across per-cell-type runs) weights them by
    WOT's growth weights, so terminal distributions, and with them the trajectory
    clusters and fate probabilities, follow descendant mass as in WOT's
    ``TransportMapModel.fates``; see :func:`CalculateTerminalFatesAcrossCells`.
    """
    if S is None:
        if OT_GRAPH_KEY not in adata.obsp:
            raise KeyError(f"Missing adata.obsp[{OT_GRAPH_KEY!r}]; pass S explicitly or run get_trajectories first.")
        S = adata.obsp[OT_GRAPH_KEY]
    return terminal_fate_trajectory_clustering(
        adata,
        S,
        cluster_key=cluster_key,
        resolution=resolution,
        terminal_topk=terminal_topk,
        graph_topk=graph_topk,
        random_state=random_state,
        directed=directed,
        growth=growth,
    )


def cluster_cell_state_trajectories(*args, **kwargs):
    """Deprecated alias for :func:`cluster_cells_to_trajectories`."""
    import warnings

    warnings.warn(
        "cluster_cell_state_trajectories is deprecated; use cluster_cells_to_trajectories.",
        FutureWarning,
        stacklevel=2,
    )
    return cluster_cells_to_trajectories(*args, **kwargs)


def terminal_fates_from_terminal_clusters(
    adata,
    terminal_cluster_key=CELL_STATE_TRAJECTORY_CLUSTER_KEY,
    out_key=TERMINAL_FATE_PROBS_KEY,
):
    """Aggregate terminal-cell distributions by terminal trajectory cluster."""
    D = adata.obsm[TERMINAL_CELL_DIST_KEY].tocsr()
    idx_last = adata.uns[TERMINAL_CELL_DIST_INDEX_KEY]
    term_clusters = adata.obs.iloc[idx_last][terminal_cluster_key].astype(str).to_numpy()
    fate_names = sorted(pd.unique(term_clusters))
    F_cols = []
    for fate in fate_names:
        cols = np.where(term_clusters == fate)[0]
        F_cols.append(np.asarray(D[:, cols].sum(axis=1)).ravel())
    F = np.vstack(F_cols).T
    F = F / np.maximum(F.sum(axis=1, keepdims=True), 1e-12)
    adata.obsm[out_key] = F
    adata.uns[f"{out_key}_fates"] = fate_names
    adata.uns[f"{out_key}_names"] = fate_names
    for j, fate in enumerate(fate_names):
        adata.obs[f"{out_key}_{fate}"] = F[:, j]
    return F, fate_names


def compute_terminal_fate_probabilities(
    adata,
    *,
    terminal_cluster_key=CELL_STATE_TRAJECTORY_CLUSTER_KEY,
    out_key=TERMINAL_FATE_PROBS_KEY,
):
    """Aggregate terminal-cell distributions into per-cell fate probabilities."""
    return terminal_fates_from_terminal_clusters(
        adata,
        terminal_cluster_key=terminal_cluster_key,
        out_key=out_key,
    )


def fate_umap_2x2_scanpy(
    adata,
    fate_key="terminal_fate_probs_traj",
    fate_names_key="terminal_fate_probs_traj_fates",
    umap_basis="umap",
    fate_obs_prefix="fate",
    which_two=(0, 1),
    max_entropy_norm=True,
    size=15,
):
    """Plot two fate probabilities plus entropy and max fate probability."""
    plt = _pyplot()
    assert fate_key in adata.obsm, f"Need adata.obsm['{fate_key}']"
    assert fate_names_key in adata.uns, f"Need adata.uns['{fate_names_key}']"
    F = np.asarray(adata.obsm[fate_key])
    fates = list(adata.uns[fate_names_key])
    i0, i1 = which_two
    for j, name in enumerate(fates):
        key = f"{fate_obs_prefix}_{name}"
        if key not in adata.obs:
            adata.obs[key] = F[:, j].astype(np.float32)
    eps = 1e-12
    P = np.clip(F, eps, 1.0)
    ent = -np.sum(P * np.log(P), axis=1)
    if max_entropy_norm and len(fates) > 1:
        ent = ent / np.log(len(fates))
    adata.obs[f"{fate_obs_prefix}_entropy"] = ent.astype(np.float32)
    adata.obs[f"{fate_obs_prefix}_maxprob"] = np.max(F, axis=1).astype(np.float32)
    fig, axes = plt.subplots(2, 2, figsize=(5, 4.5), constrained_layout=True)
    sc = _scanpy()
    sc.pl.embedding(adata, basis=umap_basis, color=f"{fate_obs_prefix}_{fates[i0]}", cmap="Reds", vmin=0, vmax=1, size=size, ax=axes[0, 0], show=False, title=f"Fate prob: {fates[i0]}")
    sc.pl.embedding(adata, basis=umap_basis, color=f"{fate_obs_prefix}_{fates[i1]}", cmap="Blues", vmin=0, vmax=1, size=size, ax=axes[0, 1], show=False, title=f"Fate prob: {fates[i1]}")
    sc.pl.embedding(adata, basis=umap_basis, color=f"{fate_obs_prefix}_entropy", cmap="viridis", vmin=0, vmax=1 if max_entropy_norm else None, size=size, ax=axes[1, 0], show=False, title="Fate entropy" + (" (norm)" if max_entropy_norm else ""))
    sc.pl.embedding(adata, basis=umap_basis, color=f"{fate_obs_prefix}_maxprob", cmap="viridis", vmin=0, vmax=1, size=size, ax=axes[1, 1], show=False, title="Fate maxprob")
    plt.show()
    return fig, axes


def mass_retention_plots_and_umap(
    adata,
    T,
    cluster_key,
    time_key="time",
    umap_basis="umap",
    n_reps=50,
    clusters_to_plot=None,
    size=15,
    eps=1e-12,
    random_state=0,
    store_prefix=None,
    show=True,
):
    """Compute and plot cluster mass retention across adjacent time bins."""
    assert sp.issparse(T), "T must be scipy sparse"
    assert T.shape == (adata.n_obs, adata.n_obs)
    assert cluster_key in adata.obs
    assert time_key in adata.obs
    if store_prefix is None:
        store_prefix = cluster_key
    T = T.tocsr().astype(np.float64).copy()
    rs = np.asarray(T.sum(axis=1)).ravel()
    zero = rs == 0
    if np.any(zero):
        T[zero, zero] = 1.0
        rs = np.asarray(T.sum(axis=1)).ravel()
    T = sp.diags(1.0 / (rs + 1e-15), format="csr") @ T
    t = adata.obs[time_key]
    time_levels = list(t.cat.categories) if hasattr(t, "cat") else sorted(pd.unique(t))
    tt = t.to_numpy().astype(object)
    idx_by_time_all = [np.where(tt == lvl)[0] for lvl in time_levels]
    idx_by_time = [ii for ii in idx_by_time_all if ii.size > 0]
    time_levels_nonempty = [lvl for lvl, ii in zip(time_levels, idx_by_time_all) if ii.size > 0]
    n_steps = len(idx_by_time) - 1
    step_labels = [f"{time_levels_nonempty[i]}->{time_levels_nonempty[i + 1]}" for i in range(n_steps)]
    lab = adata.obs[cluster_key].astype(str).to_numpy()
    clusters_all = pd.unique(lab)
    rng = np.random.default_rng(random_state)
    cell_obs = np.full(adata.n_obs, np.nan, dtype=np.float32)
    rows, rows_overall = [], []
    for ti in range(n_steps):
        idx_src, idx_dst = idx_by_time[ti], idx_by_time[ti + 1]
        A = T[idx_src, :][:, idx_dst].tocsr()
        denom = np.asarray(A.sum(axis=1)).ravel()
        src_lab, dst_lab = lab[idx_src], lab[idx_dst]
        dst_frac = (pd.Series(dst_lab).value_counts() / max(1, len(idx_dst))).to_dict()
        src_rows_by_c = {c: np.where(src_lab == c)[0] for c in clusters_all}
        total_out_by_c = {c: float(denom[src_rows_by_c[c]].sum()) if src_rows_by_c[c].size else 0.0 for c in clusters_all}
        ret_obs, base_analytic = {}, {}
        for c in clusters_all:
            rr = src_rows_by_c[c]
            total_out = total_out_by_c[c]
            if total_out <= 0:
                ret_obs[c] = np.nan
                base_analytic[c] = np.nan
                continue
            base_analytic[c] = float(dst_frac.get(c, 0.0))
            dst_cols = np.where(dst_lab == c)[0]
            same_per_row = np.asarray(A[:, dst_cols].sum(axis=1)).ravel() if dst_cols.size else np.zeros(A.shape[0])
            ret_obs[c] = float(same_per_row[rr].sum()) / total_out
            d = denom[rr]
            obs_i = np.zeros_like(d, dtype=np.float64)
            good = d > 0
            obs_i[good] = same_per_row[rr][good] / d[good]
            cell_obs[idx_src[rr]] = obs_i.astype(np.float32)
        total_out_all = float(denom.sum())
        same_mass_all = sum(float(A[src_rows_by_c[c], :][:, np.where(dst_lab == c)[0]].sum()) for c in clusters_all if total_out_by_c[c] > 0)
        base_mass_all = sum(total_out_by_c[c] * float(dst_frac.get(c, 0.0)) for c in clusters_all if total_out_by_c[c] > 0)
        overall_perm = []
        ret_perm = {c: [] for c in clusters_all}
        for _ in range(max(0, int(n_reps or 0))):
            dst_perm = np.asarray(dst_lab, dtype=object)[rng.permutation(len(dst_lab))]
            same_all = 0.0
            for c in clusters_all:
                rr = src_rows_by_c[c]
                total_out = total_out_by_c[c]
                if total_out <= 0:
                    continue
                dst_cols = np.where(dst_perm == c)[0]
                same_mass = float(A[rr, :][:, dst_cols].sum()) if dst_cols.size else 0.0
                ret_perm[c].append(same_mass / total_out)
                same_all += same_mass
            overall_perm.append(same_all / max(total_out_all, eps))
        for c in clusters_all:
            vals = ret_perm[c]
            rows.append(
                {
                    "step": ti,
                    "step_label": step_labels[ti],
                    "time_src": time_levels_nonempty[ti],
                    "time_dst": time_levels_nonempty[ti + 1],
                    "cluster": c,
                    "ret_obs": ret_obs.get(c, np.nan),
                    "base_analytic": base_analytic.get(c, np.nan),
                    "base_perm_mean": float(np.mean(vals)) if vals else np.nan,
                    "base_perm_sd": float(np.std(vals)) if vals else np.nan,
                    "total_out_mass": total_out_by_c.get(c, 0.0),
                    "n_src_cells": int((src_lab == c).sum()),
                    "n_dst_cells": int((dst_lab == c).sum()),
                }
            )
        rows_overall.append(
            {
                "step": ti,
                "step_label": step_labels[ti],
                "ret_obs": same_mass_all / max(total_out_all, eps),
                "base_analytic": base_mass_all / max(total_out_all, eps),
                "base_perm_mean": float(np.mean(overall_perm)) if overall_perm else np.nan,
                "base_perm_sd": float(np.std(overall_perm)) if overall_perm else np.nan,
                "total_out_mass": total_out_all,
            }
        )
    df = pd.DataFrame(rows)
    df_overall = pd.DataFrame(rows_overall)
    if clusters_to_plot is None:
        clusters_to_plot = list(df.groupby("cluster")["total_out_mass"].sum().sort_values(ascending=False).head(2).index)
    clusters_to_plot = [str(x) for x in clusters_to_plot]
    adata.obs[f"{store_prefix}_mass_ret_cell"] = cell_obs
    adata.uns[f"{store_prefix}_mass_retention_table"] = df
    adata.uns[f"{store_prefix}_mass_retention_overall"] = df_overall
    adata.uns[f"{store_prefix}_mass_retention_clusters_plotted"] = clusters_to_plot
    if show:
        plt = _pyplot()
        x = df_overall["step"].to_numpy()
        plt.figure(figsize=(9, 3.5))
        plt.plot(x, df_overall["ret_obs"], marker="o", label="observed")
        plt.plot(x, df_overall["base_analytic"], marker="o", linestyle="--", label="analytic baseline")
        if np.isfinite(df_overall["base_perm_mean"]).any():
            m = df_overall["base_perm_mean"].to_numpy()
            s = df_overall["base_perm_sd"].to_numpy()
            plt.plot(x, m, marker="o", linestyle="--", label="permutation baseline (mean)")
            plt.fill_between(x, m - s, m + s, alpha=0.2, label="perm +/-1 sd")
        plt.xticks(x, df_overall["step_label"], rotation=45, ha="right")
        plt.xlabel("time step")
        plt.ylabel("mass retention")
        plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
        plt.tight_layout()
        plt.show()
        _scanpy().pl.embedding(
            adata,
            basis=umap_basis,
            color=f"{store_prefix}_mass_ret_cell",
            cmap="viridis",
            vmin=0,
            vmax=1,
            size=size,
            na_color="lightgray",
            title=f"Per-cell mass retention to own cluster (k={cluster_key})",
        )
    return {"per_cluster": df, "overall": df_overall, "cell_key": f"{store_prefix}_mass_ret_cell", "clusters_plotted": clusters_to_plot}


def PlotHardClusterMassProportionsOverTime(adata, cluster_key="traj_clusters", time_key="time"):
    """Plot hard cluster proportions by time as a stacked area plot."""
    plt = _pyplot()
    assert time_key in adata.obs, f"Need adata.obs['{time_key}']"
    assert cluster_key in adata.obs, f"Need adata.obs['{cluster_key}']"
    t = adata.obs[time_key]
    c = adata.obs[cluster_key].astype(str)
    if hasattr(t, "cat"):
        time_levels = list(t.cat.categories)
    else:
        uu = pd.unique(t.astype(str))
        num = pd.to_numeric(uu, errors="coerce")
        time_levels = [str(x) for x in uu[np.argsort(num)]] if np.all(~np.isnan(num)) else sorted([str(x) for x in uu])
    tab = pd.crosstab(
        pd.Categorical(t.astype(str), categories=[str(x) for x in time_levels], ordered=True),
        pd.Categorical(c),
    )
    props = tab.div(tab.sum(axis=1), axis=0).fillna(0.0)
    x = np.arange(props.shape[0])
    plt.figure(figsize=(8, 4))
    plt.stackplot(x, props.T.to_numpy(), labels=[str(cl) for cl in props.columns])
    plt.xticks(x, props.index.to_list(), rotation=45, ha="right")
    plt.xlabel("time")
    plt.ylabel("hard proportion")
    plt.title(f"Hard mass proportions over time (stacked): {cluster_key}")
    plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
    plt.tight_layout()
    plt.show()


def trajectory_soft_memberships(
    adata,
    *,
    cluster_key: str = CELL_STATE_TRAJECTORY_CLUSTER_KEY,
    retention_key: str | None = None,
    key_added: str = "softmem",
) -> pd.DataFrame:
    """Soft memberships of two trajectory clusters from per-cell mass retention.

    ``obs[retention_key]`` (default ``f"{cluster_key}_mass_ret_cell"``, written by
    ``cp.pl.cell_trajectories.mass_retention``) is the fraction of a cell's forward
    transport mass that stays in its own cluster. It becomes the cell's membership
    of its own cluster and ``1 - retention`` that of the other cluster; cells without
    a retention value (the last time point) get full membership of their own
    cluster. Adds ``{key_added}_{cluster}`` for each cluster plus
    ``{key_added}_maxprob``, ``{key_added}_argmax`` and ``{key_added}_entropy``
    (bits) to ``adata.obs``, and returns the memberships (cells x clusters, in
    cluster order), e.g. as cell weights for a two-lineage tradeSeq fit.
    """
    retention_key = retention_key or f"{cluster_key}_mass_ret_cell"
    for key in (cluster_key, retention_key):
        if key not in adata.obs:
            raise KeyError(f"Missing adata.obs[{key!r}].")
    labels = adata.obs[cluster_key].astype(str).to_numpy()
    clusters = pd.Categorical(adata.obs[cluster_key].astype(str)).categories.tolist()
    if len(clusters) != 2:
        raise ValueError(f"Expected exactly 2 clusters in obs[{cluster_key!r}], got {len(clusters)}: {clusters}")

    p_self = adata.obs[retention_key].astype(float).to_numpy()
    p_self = np.clip(np.where(np.isfinite(p_self), p_self, 1.0), 0.0, 1.0)
    p_other = 1.0 - p_self
    P = np.column_stack([np.where(labels == cluster, p_self, p_other) for cluster in clusters])

    for j, cluster in enumerate(clusters):
        adata.obs[f"{key_added}_{cluster}"] = P[:, j].astype(np.float32)
    adata.obs[f"{key_added}_maxprob"] = P.max(axis=1).astype(np.float32)
    adata.obs[f"{key_added}_argmax"] = pd.Categorical(
        np.where(P[:, 1] > P[:, 0], clusters[1], clusters[0]), categories=clusters
    )
    entropy = -(P * np.log(P + 1e-12)).sum(axis=1) / np.log(2.0)
    adata.obs[f"{key_added}_entropy"] = entropy.astype(np.float32)

    return pd.DataFrame(
        {cluster: adata.obs[f"{key_added}_{cluster}"].to_numpy() for cluster in clusters},
        index=adata.obs_names.astype(str),
    )


def _rotate_to_time(Xc, t):
    """Rotate centred 2-D coordinates so that ``t`` increases along +y."""
    tc = t - np.nanmean(t)
    tc = np.where(np.isfinite(tc), tc, 0.0)
    v = (tc[:, None] * Xc).sum(axis=0)
    norm = np.linalg.norm(v)
    if not np.isfinite(norm) or norm < 1e-12:
        return Xc.copy()
    ux, uy = v / (norm + 1e-12)
    R = np.array([[uy, -ux],
                  [ux, uy]], dtype=float)
    Xr = Xc @ R.T
    ok = np.isfinite(Xr[:, 1]) & np.isfinite(t)
    if ok.sum() >= 3 and np.corrcoef(Xr[ok, 1], t[ok])[0, 1] < 0:
        Xr[:, 1] *= -1
    return Xr


def rotate_embedding_to_time(
    adata,
    *,
    basis=OT_UMAP_KEY,
    time_key="time",
    key_added=None,
    flip_x=False,
    groupby=None,
):
    """Rotate a 2-D embedding so that time increases along +y.

    The centred coordinates are rotated so the direction most correlated with
    ``obs[time_key]`` becomes the y-axis (time increasing upward) and recentred.
    A rotation cannot undo a mirror image; ``flip_x=True`` also negates x (e.g.
    to match an earlier layout). With ``groupby`` (an ``obs`` column, e.g. the
    cell type), each group is rotated about its own centroid, which stays where
    it was, so separate groups keep their layout (groups of fewer than 3 cells
    are not rotated; ``flip_x`` mirrors each group about its centroid). Stores
    the result in ``adata.obsm[key_added]`` (default ``f"{basis}_rot"``) and
    returns the key.
    """
    if basis not in adata.obsm:
        raise KeyError(f"Missing adata.obsm[{basis!r}].")
    key_added = f"{basis}_rot" if key_added is None else key_added
    X = np.asarray(adata.obsm[basis], dtype=float)[:, :2].copy()
    t = np.asarray(adata.obs[time_key]).astype(float)

    if groupby is None:
        Xr = _rotate_to_time(X - X.mean(axis=0, keepdims=True), t)
        if flip_x:
            Xr[:, 0] *= -1
        Xr -= Xr.mean(axis=0, keepdims=True)
    else:
        groups = np.asarray(adata.obs[groupby]).astype(str)
        Xr = np.full_like(X, np.nan)
        for group in np.unique(groups):
            m = groups == group
            centroid = np.nanmean(X[m], axis=0, keepdims=True)
            Xg = X[m] - centroid
            if m.sum() >= 3:
                Xg = _rotate_to_time(Xg, t[m])
            if flip_x:
                Xg[:, 0] *= -1
            Xr[m] = Xg - np.nanmean(Xg, axis=0, keepdims=True) + centroid

    adata.obsm[key_added] = Xr
    return key_added


def early_fate_correlations(
    adata,
    fate0_key="fate_0",
    fate1_key="fate_1",
    time_key="time",
    early_bins=(1,),
    layer=None,
    use_raw_if_available=True,
    min_detect_rate=0.10,
    method="spearman",
):
    """Correlate early-bin gene expression with the fate-probability difference.

    For cells whose ``obs[time_key]`` is in ``early_bins``, correlates each gene
    (``layers[layer]``, else ``.raw`` if present and ``use_raw_if_available``,
    else ``.X``) with ``fate1 - fate0`` after removing per-bin means. Genes
    detected in fewer than ``min_detect_rate`` of those cells are dropped.
    Returns ``gene, r, pval, qval, detect_rate_early`` (t-approximation p,
    Benjamini-Hochberg q) sorted by ``r`` descending; r > 0 predicts fate 1.
    """
    from scipy.stats import rankdata, t as tdist

    # --- pick matrix and gene names ---
    if use_raw_if_available and (adata.raw is not None) and (layer is None):
        X = adata.raw.X
        genes = adata.raw.var_names.to_numpy()
    else:
        X = adata.layers[layer] if (layer is not None) else adata.X
        genes = adata.var_names.to_numpy()

    # --- subset to early bins ---
    early_bins = list(early_bins) if isinstance(early_bins, (list, tuple, np.ndarray)) else [early_bins]
    early_mask = adata.obs[time_key].isin(early_bins).to_numpy()

    if early_mask.sum() < 3:
        raise ValueError(f"Too few cells in early_bins={early_bins}: n={early_mask.sum()}")

    if fate0_key not in adata.obs or fate1_key not in adata.obs:
        raise KeyError(f"Need adata.obs['{fate0_key}'] and adata.obs['{fate1_key}']")

    y = (adata.obs.loc[early_mask, fate1_key].to_numpy()
         - adata.obs.loc[early_mask, fate0_key].to_numpy()).astype(float)

    X_e = X[early_mask]
    X_e = X_e.toarray() if sp.issparse(X_e) else np.array(X_e)

    # --- remove per-bin means (controls for bin effect across early bins) ---
    bins = adata.obs.loc[early_mask, time_key].astype(str).to_numpy()
    for b in np.unique(bins):
        idx = np.where(bins == b)[0]
        if len(idx) == 0:
            continue
        X_e[idx, :] -= X_e[idx, :].mean(axis=0, keepdims=True)
        y[idx] -= y[idx].mean()

    # --- filter low-detection genes ---
    detect = (X_e > 0).mean(axis=0)
    keep = detect >= float(min_detect_rate)
    X_e = X_e[:, keep]
    genes_keep = genes[keep]
    detect_keep = detect[keep]

    if X_e.shape[1] == 0:
        raise ValueError(f"No genes left after detect-rate filter (min_detect_rate={min_detect_rate}).")

    # --- correlations ---
    method = method.lower()
    if method == "spearman":
        Xr = np.apply_along_axis(rankdata, 0, X_e)
        yr = rankdata(y)
        Xc = Xr - Xr.mean(axis=0)
        yc = yr - yr.mean()
    elif method == "pearson":
        Xc = X_e - X_e.mean(axis=0)
        yc = y - y.mean()
    else:
        raise ValueError("method must be 'spearman' or 'pearson'")

    num = Xc.T @ yc
    den = np.sqrt((Xc**2).sum(axis=0) * (yc**2).sum())
    r = (num / np.maximum(den, 1e-12)).astype(float)

    # p-values (t approximation)
    n = X_e.shape[0]
    tstat = r * np.sqrt((n - 2) / np.maximum(1e-12, 1 - r**2))
    p = 2 * tdist.sf(np.abs(tstat), df=n - 2)

    # Benjamini-Hochberg FDR
    m = len(p)
    order = np.argsort(p)
    ranks = np.arange(1, m + 1)
    q = np.empty_like(p)
    q[order] = p[order] * m / ranks
    for i in range(m - 2, -1, -1):
        q[order[i]] = min(q[order[i]], q[order[i + 1]])

    return pd.DataFrame({
        "gene": genes_keep,
        "r": r,
        "pval": p,
        "qval": q,
        "detect_rate_early": detect_keep,
    }).sort_values("r", ascending=False).reset_index(drop=True)


__all__ = [
    "CalculateTerminalFatesAcrossCells",
    "CellStateTrajectoryResult",
    "PlotHardClusterMassProportionsOverTime",
    "build_transition_matrices",
    "cluster_cell_state_trajectories",
    "cluster_cells_to_trajectories",
    "compute_terminal_fate_probabilities",
    "cosine_topk_leiden_from_terminal_dist",
    "coupling_from_kernel",
    "early_fate_correlations",
    "embed_block",
    "extract_wot_couplings",
    "fate_umap_2x2_scanpy",
    "get_trajectories",
    "keep_topk",
    "mass_retention_plots_and_umap",
    "merge_trajectory_graphs",
    "ot_force_directed",
    "read_tmap",
    "rotate_embedding_to_time",
    "run_wot_transport_workflow",
    "sparse_adjacency_from_tmats",
    "terminal_fate_trajectory_clustering",
    "terminal_fates_from_terminal_clusters",
    "tmap_growth",
    "topk_csr_rows",
    "trajectory_soft_memberships",
]
