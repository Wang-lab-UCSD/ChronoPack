import json
import re
from pathlib import Path

DEMO = Path(__file__).resolve().parents[1] / "demo"
NOTEBOOK = DEMO / "K562_TNFa_demo.ipynb"


def _nb():
    return json.loads(NOTEBOOK.read_text())


def _code(nb):
    return ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code"]


def test_demo_notebook_uses_the_public_api_and_relative_paths():
    nb = _nb()
    joined = "\n".join(_code(nb))
    assert nb["metadata"]["kernelspec"]["name"] == "chronopack"
    assert "import chronopack as cp" in joined
    assert 'os.environ.get("CHRONOPACK_DEMO_DATA", "data")' in joined
    assert "/new-stg" not in joined
    assert not re.search(r"^\s*def ", joined, flags=re.M)
    assert not re.search(r"chronopack\.(_|\w+\._)", joined)  # no private modules
    for call in (
        "cp.tl.get_trajectories(",
        "cp.tl.cluster_cells_to_trajectories(",
        "cp.tl.compute_terminal_fate_probabilities(",
        "cp.tl.trajectory_soft_memberships(",
        "cp.tl.identify_tag_bias_candidates(",
        "cp.tl.fit_gene_dynamics(",
        'lineage_weights="cell_weights.csv"',
        "cp.pl.association_volcano(",
        "cp.pl.gene_dynamics_curves(",
    ):
        assert call in joined, call


def test_demo_notebook_keeps_executed_outputs():
    nb = _nb()
    code_cells = [c for c in nb["cells"] if c["cell_type"] == "code" and "".join(c["source"]).strip()]
    assert all(c["execution_count"] for c in code_cells)
    assert not any(o["output_type"] == "error" for c in code_cells for o in c["outputs"])
    assert sum("image/png" in o.get("data", {}) for c in code_cells for o in c["outputs"]) >= 8


def test_demo_data_builder_writes_the_files_the_notebook_reads():
    script = (DEMO / "make_demo_data.py").read_text()
    joined = "\n".join(_code(_nb()))
    for name in ("k562_tnfa.hvg.h5ad", "k562_tnfa.counts.h5ad", "k562_neg-ctrl_mtx-agg.pkl"):
        assert name in script and name in joined, name
