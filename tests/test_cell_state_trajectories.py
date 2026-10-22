import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from chronopack.tl._trajectories import (
    CellStateTrajectoryResult,
    CalculateTerminalFatesAcrossCells,
    _transitions_long_from_tmaps,
    build_transition_matrices,
    compute_terminal_fate_probabilities,
    cosine_topk_leiden_from_terminal_dist,
    embed_transition_matrices,
    prune_transition_rows,
    sparse_adjacency_from_tmats,
    transition_matrices_from_tmaps,
)


def test_build_transition_matrices_renormalizes_rows():
    df = pd.DataFrame(
        {
            "t_source": [0.0, 0.0, 0.0],
            "t_target": [1.0, 1.0, 1.0],
            "cell_i": ["c0", "c0", "c1"],
            "cell_j": ["c2", "c3", "c2"],
            "p_ij": [2.0, 2.0, 5.0],
        }
    )

    tmats = build_transition_matrices(df)
    mat = tmats[0.0]["csr"]

    np.testing.assert_allclose(np.asarray(mat.sum(axis=1)).ravel(), [1.0, 1.0])
    np.testing.assert_allclose(mat.toarray(), [[0.5, 0.5], [1.0, 0.0]])
    assert tmats[0.0]["src_ids"].tolist() == ["c0", "c1"]
    assert tmats[0.0]["tgt_ids"].tolist() == ["c2", "c3"]


def test_sparse_adjacency_from_tmats_returns_square_cell_graphs():
    adata = AnnData(
        X=np.zeros((4, 2)),
        obs=pd.DataFrame(index=["c0", "c1", "c2", "c3"]),
        var=pd.DataFrame(index=["g0", "g1"]),
    )
    tmats = {
        0.0: {
            "csr": sp.csr_matrix([[0.5, 0.5], [1.0, 0.0]]),
            "src_ids": np.array(["c0", "c1"]),
            "tgt_ids": np.array(["c2", "c3"]),
        }
    }

    T_full, T_pruned, T_sym = sparse_adjacency_from_tmats(adata, tmats, topk=1, symmetrize="maximum")

    assert T_full.shape == (4, 4)
    assert T_pruned.shape == (4, 4)
    assert T_sym.shape == (4, 4)
    assert T_pruned[1, 2] == 1.0
    assert T_sym[2, 1] == 1.0


def test_terminal_fate_probabilities_are_cell_aligned_and_normalized():
    adata = AnnData(
        X=np.zeros((4, 2)),
        obs=pd.DataFrame(
            {"time": [0, 0, 1, 1], "traj_clusters": ["early", "early", "A", "B"]},
            index=["c0", "c1", "c2", "c3"],
        ),
        var=pd.DataFrame(index=["g0", "g1"]),
    )
    T = sp.csr_matrix(
        [
            [0.0, 0.0, 0.25, 0.75],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0],
        ]
    )

    CalculateTerminalFatesAcrossCells(adata, T, topk=2)
    F, fate_names = compute_terminal_fate_probabilities(adata)

    assert fate_names == ["A", "B"]
    assert adata.obsm["terminal_cell_dist"].shape == (4, 2)
    assert adata.obsm["terminal_fate_probs_traj"].shape == (4, 2)
    np.testing.assert_allclose(F.sum(axis=1), np.ones(4))
    np.testing.assert_allclose(F[0], [0.25, 0.75])
    np.testing.assert_allclose(F[1], [1.0, 0.0])
    np.testing.assert_allclose(F[2], [1.0, 0.0])
    np.testing.assert_allclose(F[3], [0.0, 1.0])


def test_terminal_cell_distributions_are_exact_by_default():
    rng = np.random.default_rng(0)
    sizes = [3, 4, 5]
    time = np.repeat([0, 1, 2], sizes)
    n = time.size
    P = np.zeros((n, n))
    starts = np.cumsum([0] + sizes)
    for t in range(2):
        rows, cols = slice(starts[t], starts[t + 1]), slice(starts[t + 1], starts[t + 2])
        block = rng.random((sizes[t], sizes[t + 1]))
        P[rows, cols] = block / block.sum(axis=1, keepdims=True)
    adata = AnnData(X=np.zeros((n, 1)), obs=pd.DataFrame({"time": time}, index=[f"c{i}" for i in range(n)]),
                    var=pd.DataFrame(index=["g0"]))

    D = CalculateTerminalFatesAcrossCells(adata, sp.csr_matrix(P)).toarray()

    last = slice(starts[2], starts[3])
    np.testing.assert_allclose(D[: starts[1]], (P @ P)[: starts[1], last], atol=1e-6)
    np.testing.assert_allclose(D[starts[1]: starts[2]], P[starts[1]: starts[2], last], atol=1e-6)
    np.testing.assert_allclose(D[last], np.eye(sizes[2]), atol=1e-6)

    D_top = CalculateTerminalFatesAcrossCells(adata, sp.csr_matrix(P), topk=2).toarray()
    assert (np.count_nonzero(D_top, axis=1) <= 2).all()
    assert not np.allclose(D_top, D)


def test_cell_state_trajectory_result_supports_legacy_mapping_access():
    matrix = sp.eye(2, format="csr")
    result = CellStateTrajectoryResult(
        tmaps=[],
        transitions_long=pd.DataFrame(),
        tmats={},
        T_full=matrix,
        T_pruned=matrix,
        T_sym=matrix,
    )

    assert result["T_pruned"] is matrix
    assert result.as_dict()["T_sym"] is matrix
    assert set(result.keys()) == {"tmaps", "transitions_long", "tmats", "T_full", "T_pruned", "T_sym", "growth"}


def test_merge_trajectory_graphs_matches_block_embedding():
    import chronopack as cp
    from chronopack.tl._trajectories import embed_block

    adata = AnnData(
        X=np.zeros((5, 1)),
        obs=pd.DataFrame(index=["a0", "b0", "a1", "b1", "b2"]),
        var=pd.DataFrame(index=["g0"]),
    )
    a_cells, b_cells = ["a0", "a1"], ["b0", "b1", "b2"]
    A = sp.csr_matrix(np.array([[0.0, 1.0], [0.5, 0.5]], dtype=np.float32))
    B = sp.csr_matrix(np.arange(9, dtype=np.float32).reshape(3, 3))

    merged = cp.tl.merge_trajectory_graphs(adata, [(a_cells, A), (b_cells, B)])

    idx = pd.Index(adata.obs_names)
    expected = embed_block(A, idx.get_indexer(a_cells), 5) + embed_block(B, idx.get_indexer(b_cells), 5)
    assert merged.dtype == expected.dtype
    np.testing.assert_array_equal(merged.toarray(), expected.toarray())
    assert merged[2, 0] == 0.5  # a1 -> a0
    assert merged[4, 3] == 7.0  # b2 -> b1


def test_merge_trajectory_graphs_rejects_overlapping_or_unknown_cells():
    import pytest

    import chronopack as cp

    adata = AnnData(X=np.zeros((2, 1)), obs=pd.DataFrame(index=["c0", "c1"]), var=pd.DataFrame(index=["g0"]))
    one = sp.csr_matrix([[1.0]])
    with pytest.raises(ValueError, match="overlap"):
        cp.tl.merge_trajectory_graphs(adata, [(["c0"], one), (["c0"], one)])
    with pytest.raises(ValueError, match="not in adata.obs_names"):
        cp.tl.merge_trajectory_graphs(adata, [(["cX"], one)])


def test_fate_probabilities_plot_writes_fate_summaries(tmp_path, monkeypatch):
    import matplotlib

    matplotlib.use("Agg")
    import chronopack as cp

    adata = AnnData(
        X=np.zeros((3, 1)),
        obs=pd.DataFrame(index=["c0", "c1", "c2"]),
        var=pd.DataFrame(index=["g0"]),
        obsm={"ot_umap": np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
              "terminal_fate_probs_traj": np.array([[1.0, 0.0], [0.5, 0.5], [0.2, 0.8]])},
        uns={"terminal_fate_probs_traj_fates": ["0", "1"]},
    )

    cp.pl.cell_trajectories.fate_probabilities(adata, size=5)

    np.testing.assert_allclose(adata.obs["fate_1"], [0.0, 0.5, 0.8])
    np.testing.assert_allclose(adata.obs["fate_maxprob"], [1.0, 0.5, 0.8])
    assert adata.obs["fate_entropy"].iloc[1] == np.float32(1.0)


def test_trajectory_soft_memberships_from_mass_retention():
    import anndata as ad
    import numpy as np
    import pandas as pd

    import chronopack as cp

    obs = pd.DataFrame(
        {
            "traj_clusters": pd.Categorical(["0", "0", "1", "1"]),
            "traj_clusters_mass_ret_cell": [0.9, np.nan, 0.6, 1.2],  # NaN = last time point; > 1 is clipped
        },
        index=["a", "b", "c", "d"],
    )
    adata = ad.AnnData(np.zeros((4, 1)), obs=obs)

    weights = cp.tl.trajectory_soft_memberships(adata)

    np.testing.assert_allclose(weights["0"], [0.9, 1.0, 0.4, 0.0], rtol=1e-6)
    np.testing.assert_allclose(weights["1"], [0.1, 0.0, 0.6, 1.0], rtol=1e-6)
    assert list(weights.columns) == ["0", "1"]
    assert adata.obs["softmem_argmax"].tolist() == ["0", "0", "1", "1"]
    np.testing.assert_allclose(adata.obs["softmem_entropy"], [0.469, 0.0, 0.971, 0.0], atol=1e-3)


def _random_tmaps(rng, sizes):
    """Dense WOT-like maps between consecutive time points, as read_tmap returns them."""
    names = [[f"t{t}_c{i}" for i in range(n)] for t, n in enumerate(sizes)]
    tmaps = []
    for t in range(len(sizes) - 1):
        T = sp.csr_matrix(rng.random((sizes[t], sizes[t + 1])))
        T = T.multiply(1.0 / np.asarray(T.sum(axis=1)).ravel()[:, None])
        tmaps.append((float(t), float(t + 1), T, np.array(names[t], dtype=object), np.array(names[t + 1], dtype=object)))
    obs_names = [name for block in names for name in block]
    return tmaps, obs_names


def _legacy_keep_topk(mat, k):
    """The previous LIL implementation of keep_topk, as a reference."""
    mat = mat.tolil()
    for i in range(mat.shape[0]):
        if len(mat.data[i]) > k:
            keep = np.argsort(mat.data[i])[-k:]
            mat.rows[i] = list(np.array(mat.rows[i])[keep])
            mat.data[i] = list(np.array(mat.data[i])[keep])
    return mat.tocsr()


def test_transition_matrices_from_tmaps_match_the_long_table_route():
    rng = np.random.default_rng(1)
    tmaps, _ = _random_tmaps(rng, [6, 7, 5])
    tmaps = tmaps[::-1]  # file order (tmaps_10_11 before tmaps_1_2) is not time order

    direct = transition_matrices_from_tmaps(tmaps)
    via_long = build_transition_matrices(_transitions_long_from_tmaps(tmaps))

    assert list(direct) == list(via_long)
    for t in direct:
        assert direct[t]["src_ids"].tolist() == via_long[t]["src_ids"].tolist()
        assert direct[t]["tgt_ids"].tolist() == via_long[t]["tgt_ids"].tolist()
        assert direct[t]["t_next"] == via_long[t]["t_next"]
        np.testing.assert_array_equal(direct[t]["csr"].toarray(), via_long[t]["csr"].toarray())


def test_pruned_graph_matches_pruning_the_full_matrix():
    rng = np.random.default_rng(2)
    tmaps, obs_names = _random_tmaps(rng, [40, 50, 45, 35])
    tmats = transition_matrices_from_tmaps(tmaps)
    adata = AnnData(X=np.zeros((len(obs_names), 1)), obs=pd.DataFrame(index=obs_names), var=pd.DataFrame(index=["g0"]))

    T_full, T_pruned, T_sym = sparse_adjacency_from_tmats(adata, tmats, topk=30, symmetrize="maximum")

    np.testing.assert_array_equal(T_pruned.toarray(), _legacy_keep_topk(T_full, 30).toarray())
    np.testing.assert_array_equal(T_sym.toarray(), np.maximum(T_pruned.toarray(), T_pruned.toarray().T))
    np.testing.assert_array_equal(embed_transition_matrices(obs_names, tmats).toarray(), T_full.toarray())


def test_prune_transition_rows_by_topk_and_mass():
    M = sp.csr_matrix(np.array([
        [0.1, 0.4, 0.0, 0.3, 0.2],
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0],
        [0.25, 0.25, 0.25, 0.25, 0.0],
    ]))

    np.testing.assert_array_equal(prune_transition_rows(M, topk=2).toarray()[0], [0, 0.4, 0, 0.3, 0])
    # 0.4 + 0.3 = 70% of the row: 0.7 keeps two entries, 0.71 needs a third.
    np.testing.assert_array_equal(prune_transition_rows(M, mass=0.7).toarray()[0], [0, 0.4, 0, 0.3, 0])
    np.testing.assert_array_equal(prune_transition_rows(M, mass=0.71).toarray()[0], [0, 0.4, 0, 0.3, 0.2])
    np.testing.assert_array_equal(prune_transition_rows(M, mass=0.01).toarray()[2], [0, 0, 1.0, 0, 0])
    assert prune_transition_rows(M, mass=0.9, topk=1).toarray()[0].tolist() == [0, 0.4, 0, 0, 0]
    assert prune_transition_rows(M, topk=2).getrow(1).nnz == 0
    assert prune_transition_rows(M, topk=2).getrow(3).nnz == 2  # ties keep the earlier columns
    np.testing.assert_array_equal(prune_transition_rows(M).toarray(), M.toarray())

    rng = np.random.default_rng(3)
    R = sp.random(30, 40, density=0.5, random_state=4, format="csr")
    for k in (1, 5, 50):
        np.testing.assert_array_equal(prune_transition_rows(R, topk=k).toarray(), _legacy_keep_topk(R, k).toarray())
    pruned = prune_transition_rows(R, mass=0.8).toarray()
    for dense_row, kept in zip(R.toarray(), pruned):
        if dense_row.sum() == 0:
            continue
        values = np.sort(dense_row)[::-1]
        n_needed = int(np.searchsorted(np.cumsum(values), 0.8 * dense_row.sum() - 1e-12) + 1)
        assert np.count_nonzero(kept) == n_needed
        assert kept.sum() >= 0.8 * dense_row.sum() - 1e-12


def test_result_builds_full_matrix_and_long_table_only_when_accessed():
    calls = []

    def full():
        calls.append("T_full")
        return sp.eye(2, format="csr")

    result = CellStateTrajectoryResult(tmaps=[], transitions_long=lambda: calls.append("long") or pd.DataFrame(),
                                       tmats={}, T_full=full, T_pruned=None, T_sym=None)
    assert result["T_pruned"] is None and calls == []
    assert set(result.keys()) == {"tmaps", "transitions_long", "tmats", "T_full", "T_pruned", "T_sym", "growth"}
    assert result["T_full"].shape == (2, 2) and calls == ["T_full"]
    result.T_full
    assert calls == ["T_full"]
    result.as_dict()
    assert calls == ["T_full", "long"]


def test_dense_cosine_graph_matches_sklearn_sparse_neighbors():
    from sklearn.neighbors import NearestNeighbors
    from sklearn.preprocessing import normalize

    from chronopack.tl._trajectories import _topk_csr_rows

    rng = np.random.default_rng(5)
    n, m, k = 120, 25, 10
    D = rng.random((n, m)) ** 4
    D[D < 0.05] = 0
    D[-m:] = np.eye(m)  # terminal cells: one-hot rows, as in the exact fate distributions
    D = sp.csr_matrix(D / D.sum(axis=1, keepdims=True), dtype=np.float32)
    adata = AnnData(X=np.zeros((n, 1)), obs=pd.DataFrame(index=[f"c{i}" for i in range(n)]), var=pd.DataFrame(index=["g0"]))
    adata.obsm["terminal_cell_dist"] = D

    S = cosine_topk_leiden_from_terminal_dist(adata, topk=k)

    # The previous implementation: sklearn brute-force cosine neighbours on the sparse matrix.
    Xn = normalize(D, norm="l2", axis=1, copy=True)
    nn = NearestNeighbors(n_neighbors=k + 1, metric="cosine", algorithm="brute").fit(Xn)
    ref = nn.kneighbors_graph(Xn, mode="distance").tocsr()
    ref.data = (1.0 - ref.data).astype(np.float32)
    ref.setdiag(0.0)
    ref.eliminate_zeros()
    ref = _topk_csr_rows(ref.maximum(ref.T).tocsr(), k=k)

    assert S.shape == ref.shape
    np.testing.assert_allclose(S.toarray(), ref.toarray(), atol=1e-6)
    assert (S.toarray() > 0).sum() == (ref.toarray() > 0).sum()


def _chain_graph_adata(n_per_time=20, n_times=3, seed=0):
    """Symmetric graph linking each time slice to the next, with some explicit zero entries."""
    rng = np.random.default_rng(seed)
    n = n_per_time * n_times
    rows, cols, vals = [], [], []
    for t in range(n_times - 1):
        for i in range(t * n_per_time, (t + 1) * n_per_time):
            js = rng.choice(np.arange((t + 1) * n_per_time, (t + 2) * n_per_time), size=5, replace=False)
            rows += [i] * 5
            cols += js.tolist()
            vals += rng.random(5).tolist()
    vals[0] = 0.0  # stored zero, as in the OT graph
    G = sp.csr_matrix((vals, (rows, cols)), shape=(n, n))
    G = G.maximum(G.T).tocsr()
    adata = AnnData(obs=pd.DataFrame({"time": np.repeat(np.arange(n_times), n_per_time)}, index=[f"c{i}" for i in range(n)]))
    adata.obsp["ot_graph"] = G
    return adata


def test_ot_force_directed_matches_weighted_igraph_layout():
    import random

    import igraph

    from chronopack.tl import ot_force_directed

    adata = _chain_graph_adata()
    adata.uns["draw_graph"] = {"params": {"layout": "existing"}}
    ot_force_directed(adata, random_state=3)
    Y = adata.obsm["ot_force_directed"]
    assert Y.shape == (adata.n_obs, 2) and np.isfinite(Y).all()
    assert not [k for k in adata.obsm if k.startswith("X_draw_graph")]
    assert adata.uns["draw_graph"] == {"params": {"layout": "existing"}}
    prov = adata.uns["chronopack"]["ot_force_directed"]
    assert prov["parameters"]["layout"] == "fr" and prov["parameters"]["weighted"]

    # Same as a weighted Fruchterman-Reingold layout in igraph from scanpy's random start.
    G = sp.triu(adata.obsp["ot_graph"], k=1).tocoo()
    keep = G.data != 0
    g = igraph.Graph(n=adata.n_obs, edges=list(zip(G.row[keep].tolist(), G.col[keep].tolist())), directed=False)
    g.es["weight"] = G.data[keep].tolist()
    np.random.seed(3)
    init = np.random.random((adata.n_obs, 2))
    random.seed(3)
    expected = np.asarray(g.layout_fruchterman_reingold(weights="weight", niter=1000, seed=init.tolist()).coords)
    np.testing.assert_allclose(Y, expected, rtol=1e-6, atol=1e-8)

    # Deterministic, and the weights matter.
    again = ot_force_directed(adata, random_state=3, copy=True)
    np.testing.assert_array_equal(again.obsm["ot_force_directed"], Y)
    unweighted = ot_force_directed(adata, random_state=3, weighted=False, copy=True)
    assert not np.allclose(unweighted.obsm["ot_force_directed"], Y)


def test_ot_force_directed_errors():
    import pytest
    from importlib.util import find_spec

    from chronopack.tl import ot_force_directed

    adata = _chain_graph_adata()
    with pytest.raises(KeyError):
        ot_force_directed(adata, graph_key="missing")
    if find_spec("fa2_modified") is None:
        with pytest.raises(ImportError):
            ot_force_directed(adata, layout="fa")


def _write_raw_tmaps(path, rng, sizes):
    """WOT-style map files (tmaps_<t>_<t+1>.h5ad) holding couplings with uneven row sums (growth)."""
    names = [[f"t{t}_c{i}" for i in range(n)] for t, n in enumerate(sizes, start=1)]
    couplings = []
    for t in range(len(sizes) - 1):
        C = rng.random((sizes[t], sizes[t + 1])) ** 3 * rng.uniform(0.2, 3.0, size=(sizes[t], 1))
        AnnData(X=C, obs=pd.DataFrame(index=names[t]), var=pd.DataFrame(index=names[t + 1])).write_h5ad(
            path / f"tmaps_{t + 1}_{t + 2}.h5ad")
        couplings.append(C)
    obs = pd.DataFrame({"time": np.repeat(np.arange(1, len(sizes) + 1), sizes)},
                       index=[name for block in names for name in block])
    return couplings, obs


def test_workflow_keeps_wot_growth_weights(tmp_path):
    import pytest

    from chronopack.tl import coupling_from_kernel
    from chronopack.tl._trajectories import run_wot_transport_workflow, tmap_growth

    couplings, obs = _write_raw_tmaps(tmp_path, np.random.default_rng(7), [5, 6, 4])
    adata = AnnData(X=np.zeros((len(obs), 1)), obs=obs, var=pd.DataFrame(index=["g0"]))
    result = run_wot_transport_workflow(adata, tmp_path, compute_transport_maps=False, transition_topk=None,
                                        attach_graph=False)

    expected = np.concatenate([C.sum(axis=1) / C.sum(axis=1).mean() for C in couplings])
    assert result.growth.index.equals(adata.obs_names)
    np.testing.assert_allclose(result.growth.iloc[:11], expected)
    assert result.growth.iloc[11:].isna().all()  # the last time point has no outgoing map
    pd.testing.assert_series_equal(tmap_growth(tmp_path).sort_index(), result.growth.dropna().sort_index(),
                                   check_names=False)
    assert "ot_growth" not in adata.obs  # only written with attach_graph

    # The couplings come back from the row-normalized matrix, each map divided by its mean row sum.
    C = coupling_from_kernel(result.T_full, result.growth).toarray()
    np.testing.assert_allclose(C[:5, 5:11], couplings[0] / couplings[0].sum(axis=1).mean())
    np.testing.assert_allclose(C[5:11, 11:], couplings[1] / couplings[1].sum(axis=1).mean())
    assert not C[11:].any()
    with pytest.raises(ValueError):
        coupling_from_kernel(result.T_full, np.full(adata.n_obs, np.nan))


def test_growth_weighted_fates_match_wot(tmp_path):
    import pytest

    wot = pytest.importorskip("wot")
    from chronopack.tl._trajectories import run_wot_transport_workflow, terminal_fates_from_terminal_clusters

    sizes = [12, 15, 10, 14]
    _, obs = _write_raw_tmaps(tmp_path, np.random.default_rng(8), sizes)
    last = obs.index[obs["time"] == len(sizes)]
    obs["traj_clusters"] = "0"
    obs.loc[last, "traj_clusters"] = np.array(["a", "b"])[np.arange(len(last)) % 2]
    adata = AnnData(X=np.zeros((len(obs), 1)), obs=obs, var=pd.DataFrame(index=["g0"]))
    result = run_wot_transport_workflow(adata, tmp_path, compute_transport_maps=False, transition_topk=None,
                                        attach_graph=False)

    model = wot.tmap.TransportMapModel.from_directory(str(tmp_path / "tmaps"))
    sets = {f: list(last[obs.loc[last, "traj_clusters"] == f]) for f in ("a", "b")}
    res = model.fates(model.population_from_cell_sets(sets, at_time=float(len(sizes))))
    expected = pd.DataFrame(res.X, index=res.obs_names, columns=res.var_names).loc[adata.obs_names, ["a", "b"]]

    def fates(growth):
        tmp = adata.copy()
        CalculateTerminalFatesAcrossCells(tmp, result.T_full, growth=growth)
        F, names = terminal_fates_from_terminal_clusters(tmp)
        return pd.DataFrame(F, index=tmp.obs_names, columns=names)[["a", "b"]].to_numpy(), tmp

    weighted, tmp = fates(result.growth)
    np.testing.assert_allclose(weighted, expected.to_numpy(), atol=1e-6)  # terminal_cell_dist is stored as float32
    assert tmp.uns["chronopack"]["terminal_cell_dist"]["parameters"]["growth"] == "per-cell values"
    default, tmp = fates(None)
    assert np.abs(default - expected.to_numpy()).max() > 1e-2  # the Markov-chain default differs
    assert tmp.uns["chronopack"]["terminal_cell_dist"]["parameters"]["growth"] is None

    adata.obs["ot_growth"] = result.growth.to_numpy()
    by_key, tmp = fates("ot_growth")
    np.testing.assert_array_equal(by_key, weighted)
    assert tmp.uns["chronopack"]["terminal_cell_dist"]["parameters"]["growth"] == "ot_growth"

    with pytest.raises(KeyError):
        CalculateTerminalFatesAcrossCells(adata.copy(), result.T_full, growth="missing")
    bad = result.growth.copy()
    bad.iloc[0] = np.nan
    with pytest.raises(ValueError):
        CalculateTerminalFatesAcrossCells(adata.copy(), result.T_full, growth=bad)
