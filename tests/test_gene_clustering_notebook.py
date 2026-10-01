import json
import re
from pathlib import Path


NOTEBOOK = Path(__file__).resolve().parents[1] / "analyses/gene_clustering/01_GeneTrajectoryClustering.ipynb"


def _code_cells():
    nb = json.loads(NOTEBOOK.read_text())
    return ["".join(cell.get("source", "")) for cell in nb["cells"] if cell["cell_type"] == "code"]


def _active_code():
    lines = []
    for src in _code_cells():
        for line in src.splitlines():
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                lines.append(stripped)
    return "\n".join(lines)


def test_notebook_exists_and_writes_expected_outputs():
    assert NOTEBOOK.exists()
    active = _active_code()
    assert "pc_df" not in active  # traj.obs holds every per-trajectory column
    assert "gene_trajectory_characterization.h5ad" in active
    assert "gene_trajectory_clusters.h5ad" in active


def test_checkpoint_reads_from_old_dimensionality_notebooks_are_not_active_inputs():
    active = _active_code()
    banned_reads = [
        "readRDS('traj_mat.rds')",
        'readRDS("traj_mat.rds")',
        "readRDS('X.rds')",
        "readRDS('pc.rds')",
        "readRDS('tt.rds')",
        "readRDS('asso_df_ls.rds')",
        "readRDS('gamList_ls.rds')",
        "readRDS('sce_ls.rds')",
        "readRDS('counts_sub_ls.rds')",
        "readRDS('PT_ls.rds')",
        "readRDS('CW_ls.rds')",
        "readRDS('sf_ls.rds')",
        "readRDS('K562a_asso-df.rds')",
        "readRDS('K562b_asso-df.rds')",
        "readRDS('THP1_asso-df.rds')",
    ]
    for banned in banned_reads:
        assert banned not in active


def test_all_three_contexts_are_loaded_from_primary_tradeseq_fits():
    active = _active_code()
    assert "cp.tl.align_gene_trajectories(" in active
    assert "cp.tl.compute_gene_trajectory_descriptors(" in active
    assert "cp.tl.cluster_gene_trajectories(" in active
    assert "cp.pl.clustering_metrics(" in active
    for context in ['"K562a"', '"K562b"', '"THP1"']:
        assert context in active
    assert "puro" not in active.lower()


def test_dimred_notebook_uses_public_chronopack_api():
    active = _active_code()
    assert "import chronopack as cp" in active
    assert "run_dimred_signed_clustering(" not in active
    assert "prepare_dimred_clustering_inputs(" not in active
    assert "finalize_gene_trajectory_clusters(" not in active
    assert "cl_repressed_renamed" not in active


def test_notebook_has_no_inline_function_definitions():
    for source in _code_cells():
        stripped = "\n".join(
            line for line in source.splitlines() if not line.lstrip().startswith("#")
        )
        assert "<- function" not in stripped
