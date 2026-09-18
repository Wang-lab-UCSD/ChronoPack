try:
    from IPython.display import display
except Exception:
    display = print


# %% cell 1

from pathlib import Path
import sys

import numpy as np
import pandas as pd

NOTEBOOK_DIR = Path.cwd().resolve()
REPO_ROOT = NOTEBOOK_DIR.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from chronopack.grn_inference import build_chronopack_ranked_edges, write_chronopack_ranked_edges


# %% cell 2

VALID_TARGETS_PATH = Path("/new-stg/home/alex/tmp/grn_inference_v2/reproduction/inputs/external/ref/genome/hg38/annotations/protein_coding_genes.csv")
TRRUST_PATH = Path("/new-stg/home/alex/ref/grns/ttrust/trrust_rawdata.human.tsv")
CHIP_ENCODE_GT_PATH = Path("/new-stg/home/alex/downloads/k562-chip_encode/results/tf_gene.1kbp_200bp.stringent.tsv")

GENE_SCORES_PATH = NOTEBOOK_DIR / "outputs" / "tables" / "thp1_k562_lfdr_recalibration_scores.csv"
K562_TF_SCORES_PATH = NOTEBOOK_DIR / "outputs" / "omnipath_tf_ranks" / "omnipath_ranks.TNFa.csv"
THP1_TF_SCORES_PATH = NOTEBOOK_DIR / "outputs" / "omnipath_tf_ranks" / "omnipath_ranks.TLR4.csv"
OUTPUT_PATH = NOTEBOOK_DIR / "outputs" / "tables" / "regulator_prioritization_chronopack_ranked_edges.csv"
REFERENCE_RANKED_EDGES_PATH = Path("/new-stg/home/alex/tmp/grn_inference_v2/reproduction/benchmark_outputs/tables/regulator_prioritization_method_specific_ranked_edges.csv")


# %% cell 3

ranked_edges = build_chronopack_ranked_edges(
    gene_scores_path=GENE_SCORES_PATH,
    k562_tf_scores_path=K562_TF_SCORES_PATH,
    thp1_tf_scores_path=THP1_TF_SCORES_PATH,
    valid_targets_path=VALID_TARGETS_PATH,
    trrust_path=TRRUST_PATH,
    chip_encode_gt_path=CHIP_ENCODE_GT_PATH,
)
write_chronopack_ranked_edges(ranked_edges, OUTPUT_PATH)
print(OUTPUT_PATH)
print(ranked_edges.shape)
ranked_edges.head()


# %% cell 4

reference_chunks = []
for chunk in pd.read_csv(REFERENCE_RANKED_EDGES_PATH, chunksize=1_000_000):
    chunk = chunk[chunk["method"].isin(["ChronoPack_K562", "ChronoPack_THP1"])]
    if len(chunk):
        reference_chunks.append(chunk)
reference = pd.concat(reference_chunks, ignore_index=True)

print("generated", ranked_edges.shape)
print("reference", reference.shape)
print("columns equal", list(ranked_edges.columns) == list(reference.columns))
print("exact equals", ranked_edges.equals(reference))
if ranked_edges.shape == reference.shape and list(ranked_edges.columns) == list(reference.columns):
    numeric_cols = ["EdgeWeight", "rank", "percentile"]
    max_diff = max(float(np.nanmax(np.abs(ranked_edges[c].to_numpy() - reference[c].to_numpy()))) for c in numeric_cols)
    print("max numeric diff", max_diff)
