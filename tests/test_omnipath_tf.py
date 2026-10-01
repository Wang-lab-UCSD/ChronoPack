import networkx as nx
import pandas as pd

from chronopack._legacy.omnipath_tf import run_tf_diffusion_ranking


def test_signed_diffusion_absorbs_positive_and_negative_tf_mass():
    graph = nx.DiGraph()
    graph.add_edge("SRC", "TF_POS", is_stimulation=True, is_inhibition=False)
    graph.add_edge("SRC", "TF_NEG", is_stimulation=False, is_inhibition=True)

    ranking, hist = run_tf_diffusion_ranking(
        graph,
        source="SRC",
        tfs=["TF_NEG", "TF_POS"],
        alpha=0.85,
        max_steps=3,
        min_steps=0,
    )

    by_tf = ranking.set_index("TF")
    assert by_tf.loc["TF_POS", "abs_pos"] == 0.5
    assert by_tf.loc["TF_POS", "abs_neg"] == 0.0
    assert by_tf.loc["TF_NEG", "abs_pos"] == 0.0
    assert by_tf.loc["TF_NEG", "abs_neg"] == 0.5
    assert hist.loc[0, "absorbed_total"] == 1.0


def test_ligand_seed_uses_receptors_in_graph_only():
    graph = nx.DiGraph()
    graph.add_edge("REC_A", "TF1", is_stimulation=True, is_inhibition=False)
    graph.add_edge("REC_B", "TF2", is_stimulation=True, is_inhibition=False)
    lr_df = pd.DataFrame(
        {
            "ligand": ["LIG", "LIG", "LIG"],
            "receptor": ["REC_A", "REC_B", "MISSING"],
        }
    )

    ranking, _ = run_tf_diffusion_ranking(
        graph,
        source="LIG",
        tfs=["TF1", "TF2"],
        is_ligand=True,
        lr_df=lr_df,
        alpha=0.85,
        max_steps=3,
        min_steps=0,
    )

    by_tf = ranking.set_index("TF")
    assert by_tf.loc["TF1", "abs_pos"] == 0.5
    assert by_tf.loc["TF2", "abs_pos"] == 0.5
