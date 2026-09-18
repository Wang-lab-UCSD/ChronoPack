import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

NOTEBOOKS = [
    "k562-tnfa.hg38.rep/PreprocessDatasets.ipynb",
    "k562-tnfa.hg38.rep/GetTrajectories.ipynb",
    "k562-tnfa.hg38.rep/PrepareForTradeSeq.ipynb",
    "thp1-k562-lps.hg38.rep/PreprocessDatasets.ipynb",
    "thp1-k562-lps.hg38.rep/GetTrajectories.ipynb",
    "thp1-k562-lps.hg38.rep/PrepareForTradeSeq.ipynb",
    "k562-tnfa-puromycin.hg38/PreprocessDatasets.ipynb",
    "k562-tnfa-puromycin.hg38/GetTrajectories.ipynb",
    "k562-tnfa-puromycin.hg38/PrepareForTradeSeq.ipynb",
]

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
        assert "from chronopack." in joined, path


def test_shared_helpers_are_not_redefined_in_active_cells():
    for path in NOTEBOOKS:
        for source in _code_sources(path):
            stripped = "\n".join(
                line for line in source.splitlines() if not line.lstrip().startswith("#")
            )
            for marker in SHARED_DEFS:
                assert marker not in stripped, (path, marker)


def test_chronopack_modules_exist():
    for relpath in [
        "chronopack/__init__.py",
        "chronopack/preprocessing.py",
        "chronopack/trajectories.py",
        "chronopack/tradeseq.py",
    ]:
        assert (ROOT / relpath).exists()


def test_preprocess_notebooks_use_loader_and_qc_helper():
    for path in [
        "k562-tnfa.hg38.rep/PreprocessDatasets.ipynb",
        "thp1-k562-lps.hg38.rep/PreprocessDatasets.ipynb",
        "k562-tnfa-puromycin.hg38/PreprocessDatasets.ipynb",
    ]:
        joined = "\n".join(_code_sources(path))
        assert "load_time_tagged_anndata" in joined, path
        assert "filter_cells_by_qc" in joined, path
        assert "sc.pp.filter_cells" not in joined, path


def test_wot_notebooks_use_workflow_helper():
    for path in [
        "k562-tnfa.hg38.rep/GetTrajectories.ipynb",
        "thp1-k562-lps.hg38.rep/GetTrajectories.ipynb",
        "k562-tnfa-puromycin.hg38/GetTrajectories.ipynb",
    ]:
        joined = "\n".join(_code_sources(path))
        assert "run_wot_transport_workflow" in joined, path


def test_production_notebooks_do_not_contain_metacell_or_beeline_scratch_code():
    forbidden = ["metacell", "Metacell", "KMeans", "Kmeans", "kmeans", "BEELINE", "beeline", "to_beeline"]
    for path in NOTEBOOKS:
        joined = "\n".join(_code_sources(path))
        for marker in forbidden:
            assert marker not in joined, (path, marker)


def test_beeline_scratch_cells_were_extracted():
    assert (ROOT / "k562-tnfa.hg38.rep/ignore/PrepareForTradeSeq_BEELINE_extracted.ipynb").exists()
