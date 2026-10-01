from importlib import resources

import chronopack as cp


def test_public_namespaces_importable():
    assert cp.io is not None
    assert cp.tl is not None
    assert cp.pl is not None
    assert hasattr(cp.tl, "get_trajectories")
    assert hasattr(cp.tl, "fit_gene_dynamics")
    assert hasattr(cp.tl, "cluster_gene_trajectories")


def test_legacy_compatibility_modules_are_gone():
    import importlib

    import pytest

    for name in ["c2c", "tradeseq", "grn_inference", "omnipath_tf", "preprocessing", "trajectories"]:
        with pytest.raises(ModuleNotFoundError):
            importlib.import_module(f"chronopack.{name}")


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
