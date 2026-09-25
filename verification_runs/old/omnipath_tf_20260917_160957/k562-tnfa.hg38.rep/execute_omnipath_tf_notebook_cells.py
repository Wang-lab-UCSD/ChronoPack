try:
    from IPython.display import display
except Exception:
    display = print


# %% cell 1

from pathlib import Path
import sys

import pandas as pd

NOTEBOOK_DIR = Path.cwd().resolve()
ANALYSIS_DIR = NOTEBOOK_DIR.parent if NOTEBOOK_DIR.name == "notebooks" else NOTEBOOK_DIR
REPO_ROOT = ANALYSIS_DIR.parent

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chronopack.omnipath_tf import (
    build_omnipath_tf_graph,
    load_omnipath_tf_resources,
    rank_omnipath_tfs,
)


# %% cell 2

HUMAN_TF_PATH = Path("/new-stg/home/alex/ref/gene_sets/hs_tfs.txt")
LIGAND_RECEPTOR_PATH = Path("/new-stg/home/alex/ref/cellchatdb/hs/interaction_input_CellChatDB.csv")
OUTPUT_DIR = ANALYSIS_DIR / "intmd" / "omnipath_tf_ranks"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

RANKING_CONFIGS = {
    "TNFa": {"source": "TNF", "is_ligand": True},
    "TLR4": {"source": "TLR4", "is_ligand": False},
}


# %% cell 4

resources = load_omnipath_tf_resources(
    human_tf_path=HUMAN_TF_PATH,
    ligand_receptor_path=LIGAND_RECEPTOR_PATH,
)
graph_bundle = build_omnipath_tf_graph(resources)

print(f"Graph: {graph_bundle.graph.number_of_nodes():,} nodes, {graph_bundle.graph.number_of_edges():,} edges")
print(f"TFs: {len(graph_bundle.tfs):,}")


# %% cell 6

rankings = {}
histories = {}

for label, config in RANKING_CONFIGS.items():
    df, hist = rank_omnipath_tfs(
        graph_bundle,
        resources.ligand_receptor,
        source=config["source"],
        is_ligand=config["is_ligand"],
        w_signed=1.0,
        w_unsigned=0.5,
        alpha=0.85,
        max_steps=80,
        min_steps=10,
        tol_mass=1e-12,
    )
    rankings[label] = df
    histories[label] = hist
    print(label, df.shape, "steps", len(hist))
    display(df.head(10))


# %% cell 8

for label, df in rankings.items():
    path = OUTPUT_DIR / f"omnipath_ranks.{label}.csv"
    df.to_csv(path)
    print(path)
