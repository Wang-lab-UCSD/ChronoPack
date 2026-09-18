from importlib import resources

import chronopack as cp


def test_public_namespaces_importable():
    assert cp.io is not None
    assert cp.tl is not None
    assert cp.pl is not None
    assert hasattr(cp.tl, "get_trajectories")
    assert hasattr(cp.tl, "fit_gene_dynamics")
    assert hasattr(cp.tl, "cluster_gene_trajectories")


def test_legacy_imports_still_available():
    from chronopack.tradeseq import WriteTradeSeqInputsFromAnnData
    from chronopack.c2c import score_lr_pairs

    assert callable(WriteTradeSeqInputsFromAnnData)
    assert callable(score_lr_pairs)


def test_bundled_r_resources_accessible():
    from chronopack._backends.r import available_r_resources, read_r_resource_text

    resources_found = available_r_resources()
    assert "tradeseq_helpers.R" in resources_found
    assert "dimred_clustering_helpers.R" in resources_found
    assert "dimred_plotting_helpers.R" in resources_found
    assert "fit_tradeseq_sce <- function" in read_r_resource_text("tradeseq_helpers.R")


def test_bundled_r_resources_are_package_data():
    root = resources.files("chronopack._backends.r")
    assert root.joinpath("tradeseq_helpers.R").is_file()
