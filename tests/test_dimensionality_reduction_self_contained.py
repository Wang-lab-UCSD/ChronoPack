import json
import re
from pathlib import Path


NOTEBOOK = Path(
    "k562-tnfa-puromycin.hg38/DimensionalityReduction&Clustering.self_contained.ipynb"
)


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
    assert 'pc_df.new.csv' in active
    assert 'pc_df.new.rds' in active
    assert 'saveRDS(pc_df, file.path(SELF_CONTAINED_OUTDIR, "pc_df.new.rds"))' in active


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


def test_all_four_datasets_are_loaded_from_primary_tradeseq_configs():
    active = _active_code()
    assert "chronopack_dimred_config(" in active
    assert "load_dimred_fixed_tradeSeq_outputs(cfg)" in active
    for key in ["K562_puro", "K562a", "K562b", "THP1"]:
        assert re.search(rf'"{key}"', active)


def test_dimred_helpers_are_sourced_from_chronopack():
    active = _active_code()
    assert "dimred_clustering_helpers.R" in active
    assert "dimred_plotting_helpers.R" in active
    assert "run_dimred_signed_clustering(" in active
    assert "prepare_dimred_clustering_inputs(X_mat, Xd_mat)" in active
    assert "cl_repressed_renamed" not in active


def test_notebook_has_no_inline_function_definitions():
    for source in _code_cells():
        stripped = "\n".join(
            line for line in source.splitlines() if not line.lstrip().startswith("#")
        )
        assert "<- function" not in stripped


def test_optional_reference_validation_is_gated_by_environment_variable():
    active = _active_code()
    assert 'Sys.getenv("CHRONOPACK_VALIDATE_REFERENCE")' in active
    assert 'VALIDATE_AGAINST_REFERENCE <- identical(Sys.getenv("CHRONOPACK_VALIDATE_REFERENCE"), "1")' in active
