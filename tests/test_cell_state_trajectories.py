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
