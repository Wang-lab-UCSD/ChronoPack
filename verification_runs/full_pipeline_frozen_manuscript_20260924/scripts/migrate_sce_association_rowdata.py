#!/usr/bin/env python
"""Annotate fitted tradeSeq SCE RDS files with association statistics in rowData."""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import chronopack as cp
from chronopack._backends.r._bridge import source_r_resource
from rpy2.robjects import r


FIELDS = ("waldStat_1", "meanLogFC", "pvalue_1", "qval_1")


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: migrate_sce_association_rowdata.py RUN_ROOT")
    run_root = Path(sys.argv[1]).resolve()
    root_puro = run_root / "k562-tnfa-puromycin.hg38"
    root_k562a = run_root / "k562-tnfa.hg38.rep"
    root_lps = run_root / "thp1-k562-lps.hg38.rep"

    sce_paths = cp.tl.sce_paths_from_dimred_roots(root_puro, root_k562a, root_lps)
    association_paths = cp.tl.association_paths_from_dimred_roots(root_puro, root_k562a, root_lps)
    source_r_resource("dimred_clustering_helpers.R")

    for key, sce_path_str in sce_paths.items():
        sce_path = Path(sce_path_str)
        asso_path = Path(association_paths[key])
        backup_path = sce_path.with_suffix(sce_path.suffix + ".pre_rowdata_assoc.bak")
        if not backup_path.exists():
            shutil.copy2(sce_path, backup_path)

        r_sce = r["readRDS"](str(sce_path))
        r_asso = r["read.csv"](
            str(asso_path),
            **{"row.names": 1, "check.names": False},
        )
        r_sce = r["annotate_sce_association"](r_sce, r_asso)
        r["saveRDS"](r_sce, str(sce_path))

        rd = r["as.data.frame"](r("SummarizedExperiment::rowData")(r_sce))
        names = set(map(str, r["colnames"](rd)))
        missing = [field for field in FIELDS if field not in names]
        if missing:
            raise RuntimeError(f"{key}: missing rowData fields after migration: {missing}")
        print(f"{key}: annotated {sce_path} (backup: {backup_path})")


if __name__ == "__main__":
    main()
