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
    assert set(result.keys()) == {"tmaps", "transitions_long", "tmats", "T_full", "T_pruned", "T_sym"}


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
    assert set(result.keys()) == {"tmaps", "transitions_long", "tmats", "T_full", "T_pruned", "T_sym"}
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
