import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData

from chronopack.tl._trajectories import (
    CellStateTrajectoryResult,
    CalculateTerminalFatesAcrossCells,
    build_transition_matrices,
    compute_terminal_fate_probabilities,
    sparse_adjacency_from_tmats,
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
