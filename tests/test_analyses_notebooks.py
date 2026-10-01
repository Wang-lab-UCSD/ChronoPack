import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1] / "analyses"
EXPERIMENTS = ["k562-tnfa.hg38.rep", "thp1-k562-lps.hg38.rep"]
PREPROCESS = [f"{e}/01_PreprocessDatasets.ipynb" for e in EXPERIMENTS]
GET_TRAJECTORIES = [f"{e}/02_GetTrajectories.ipynb" for e in EXPERIMENTS]
PREPARE = [f"{e}/03_PrepareForTradeSeq.ipynb" for e in EXPERIMENTS]
TRADESEQ_FITS = [
    "k562-tnfa.hg38.rep/04_TradeSeqFits.ipynb",
    "k562-tnfa.hg38.rep/04_TradeSeqFits.TwoLineage.ipynb",
    "thp1-k562-lps.hg38.rep/04_TradeSeqFits.K562.ipynb",
    "thp1-k562-lps.hg38.rep/04_TradeSeqFits.THP1.ipynb",
]
GENE_CLUSTERING = ["gene_clustering/01_GeneTrajectoryClustering.ipynb"]
NOTEBOOKS = PREPROCESS + GET_TRAJECTORIES + PREPARE
VISUALIZATION = ["k562-tnfa.hg38.rep/05_Visualization.ipynb", "thp1-k562-lps.hg38.rep/05_Visualization.ipynb"]

SHARED_DEFS = [
    "def read_tmap",
    "def build_transition_matrices",
    "def TemporalPseudobulkFromAnnData",
    "def ComputeSmoothCorrelations",
    "def WriteTradeSeqInputsFromAnnData",
    "def map_time_to_group",
    "def SparseAdjacencyFromTmats",
    "def ExtractWOTCouplings",
]


def _code_sources(path):
    nb = json.loads((ROOT / path).read_text())
    return [
        "".join(cell.get("source", []))
        for cell in nb["cells"]
        if cell.get("cell_type") == "code"
    ]


def test_refactored_notebooks_import_chronopack():
    for path in NOTEBOOKS:
        joined = "\n".join(_code_sources(path))
        assert ("from chronopack." in joined) or ("import chronopack" in joined), path


def test_shared_helpers_are_not_redefined_in_active_cells():
    for path in NOTEBOOKS:
        for source in _code_sources(path):
            stripped = "\n".join(
                line for line in source.splitlines() if not line.lstrip().startswith("#")
            )
            for marker in SHARED_DEFS:
                assert marker not in stripped, (path, marker)


def test_preprocess_notebooks_use_loader_and_qc_helper():
    for path in PREPROCESS:
        joined = "\n".join(_code_sources(path))
        assert "cp.io.read_time_tagged_dge" in joined, path
        assert "sc.pp.calculate_qc_metrics" in joined, path
        assert "sc.pp.normalize_total" in joined, path
        assert "sc.pp.log1p" in joined, path
        assert "cp.io.annotate_qc_metrics" not in joined, path
        assert "cp.io.filter_cells_by_qc" not in joined, path
        assert "cp.io.initialize_expression_layers" not in joined, path
        assert "cp.io.validate_preprocessed_adata" not in joined, path


def test_wot_notebooks_use_workflow_helper():
    for path in GET_TRAJECTORIES:
        joined = "\n".join(_code_sources(path))
        assert "cp.tl.get_trajectories" in joined, path
        assert "cp.tl.cluster_cells_to_trajectories" in joined, path
        assert "cp.tl.compute_terminal_fate_probabilities" in joined, path
        assert "cp.pl.cell_trajectories.proportions_over_time" in joined, path
        assert "cp.pl.cell_trajectories.mass_retention" in joined, path
        assert "cluster_cell_state_trajectories" not in joined, path
        assert "from chronopack.trajectories import" not in joined, path


def test_production_notebooks_do_not_contain_metacell_or_beeline_scratch_code():
    forbidden = ["metacell", "Metacell", "KMeans", "Kmeans", "kmeans", "BEELINE", "beeline", "to_beeline"]
    for path in NOTEBOOKS:
        joined = "\n".join(_code_sources(path))
        for marker in forbidden:
            assert marker not in joined, (path, marker)


def test_beeline_scratch_cells_were_extracted():
    assert (ROOT / "archive/k562-tnfa.hg38.rep/PrepareForTradeSeq_BEELINE_extracted.ipynb").exists()


def test_tradeseq_fit_notebooks_use_backend_entrypoint():
    for path in TRADESEQ_FITS:
        joined = "\n".join(_code_sources(path))
        assert "cp.tl.fit_gene_dynamics(" in joined, path
        assert 'backend="tradeseq"' in joined, path
        assert "run_tradeseq_single_lineage" not in joined, path
        assert "source(" not in joined, path
        assert "find_parent_with" not in joined, path


def test_prepare_tradeseq_notebooks_use_gene_dynamics_backend():
    for path in PREPARE:
        joined = "\n".join(_code_sources(path))
        assert "cp.tl.fit_gene_dynamics" in joined, path
        assert "cp.tl.identify_tag_bias_candidates" in joined, path
        assert "cp.tl.filter_tag_bias_candidates" in joined, path
        assert "cp.tl.tag_bias_exclusion_genes" in joined, path
        assert "gd." not in joined, path
        assert "compute_smooth_correlations" not in joined, path
        assert "refilter_tag_bias_candidates" not in joined, path
        assert "WriteTradeSeqInputsFromAnnData" not in joined, path


PIPELINE_NOTEBOOKS = NOTEBOOKS + TRADESEQ_FITS + VISUALIZATION + GENE_CLUSTERING + ["c2c/CellCommunication.ipynb", "grn/GRNInference.ipynb"]
# Notebooks whose analysis logic must come from the package (no `def`; c2c/ and
# grn/ may define small row-selection helpers).
ANALYSIS_NOTEBOOKS = NOTEBOOKS + TRADESEQ_FITS + VISUALIZATION + GENE_CLUSTERING


def test_pipeline_notebooks_use_only_the_installed_public_api():
    import re

    for path in PIPELINE_NOTEBOOKS:
        nb = json.loads((ROOT / path).read_text())
        assert nb["metadata"]["kernelspec"]["name"] == "chronopack", path
        joined = "\n".join(_code_sources(path))
        assert "import chronopack as cp" in joined, path
        assert "sys.path" not in joined, path
        assert not re.search(r"(from|import) chronopack(\.\w+)*\._\w+", joined), path
        assert not re.search(r"from chronopack(\.\w+)* import _", joined), path
        assert not re.search(r"(from|import) chronopack\.(preprocessing|trajectories|tradeseq|c2c|grn_inference|omnipath_tf)\b", joined), path
        if path in ANALYSIS_NOTEBOOKS:
            assert not re.search(r"^\s*def \w+", joined, re.M), path


def test_pipeline_notebooks_read_and_write_under_the_data_dir():
    import re

    for path in PIPELINE_NOTEBOOKS:
        joined = "\n".join(_code_sources(path))
        assert 'DATA_DIR = Path(os.environ.get("CHRONOPACK_DATA_DIR", "../data"))' in joined, path
        assert not re.search(r"""['"]\.\./(input/|intmd/|output/|\.\./)""", joined), path
        assert "/new-stg/" not in joined, path
        assert "puro" not in joined.lower(), path


def test_notebook_inputs_and_resources_resolve_inside_the_repo():
    import re

    for path in PIPELINE_NOTEBOOKS:
        joined = "\n".join(_code_sources(path))
        for rel in re.findall(r"""['"](\.\./resources/[^'"]+)['"]""", joined):
            assert (ROOT / path).parent.joinpath(rel).exists(), (path, rel)
        # analyses/inputs/ is not in git; only check the layout of the references.
        for rel in re.findall(r"""['"](\.\./inputs/[^'"]+)['"]""", joined):
            assert rel.count("/") >= 2, (path, rel)


def test_visualization_notebook_uses_the_figure_api():
    for path in VISUALIZATION:
        joined = "\n".join(_code_sources(path))
        assert "cp.tl.signed_max_log2fc_from_start(" in joined, path
        assert "cp.tl.predict_gene_dynamics(" in joined, path
        assert "cp.pl.gene_dynamics_curves(" in joined, path
        assert 'FIG_DIR = EXP_DIR / "figures"' in joined, path
    k562, thp1 = ("\n".join(_code_sources(path)) for path in VISUALIZATION)
    assert "cp.pl.association_volcano(" in k562
    assert "cp.tl.early_fate_correlations(" in k562
    assert "cp.pl.cell_trajectories.trajectory_endpoints(" in thp1
    assert "cp.tl.combine_context_associations(" in thp1
    assert "cp.pl.response_specificity_volcano(" in thp1
    assert "cp.tl.fit_gene_dynamics(" in thp1
