"""Centralized AnnData keys used by ChronoPack."""

CHRONOPACK_UNS_KEY = "chronopack"

DEFAULT_TIME_KEY = "time"
DEFAULT_COUNTS_LAYER = "counts"

FATE_LABEL_KEY = "cp_fate"
FATE_PROBS_KEY = "cp_fate_probs"
TRANSITION_MATRIX_KEY = "cp_transition_matrix"

CELL_STATE_TRAJECTORIES_KEY = "cell_state_trajectories"
OT_GRAPH_KEY = "ot_graph"
OT_UMAP_KEY = "ot_umap"
OT_FORCE_DIRECTED_KEY = "ot_force_directed"
OT_GROWTH_KEY = "ot_growth"
TERMINAL_CELL_DIST_KEY = "terminal_cell_dist"
TERMINAL_CELL_DIST_INDEX_KEY = "terminal_cell_dist_idx_last"
TERMINAL_FATE_PROBS_KEY = "terminal_fate_probs_traj"
CELL_STATE_TRAJECTORY_CLUSTER_KEY = "traj_clusters"
TERMINAL_COSINE_GRAPH_KEY = "S"

GENE_DYNAMICS_KEY = "gene_dynamics"
GENE_TRAJECTORY_KEY = "gene_trajectories"
GENE_TRAJECTORY_CLUSTERING_KEY = "gene_trajectory_clusters"
GENE_MODE_KEY = "cp_mode"
GENE_FIT_KEY = "cp_fit"
