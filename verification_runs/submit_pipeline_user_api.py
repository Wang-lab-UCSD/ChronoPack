#!/usr/bin/env python
"""Stage-by-stage verification of the refactored repo-root notebooks.

Runs the public-API notebooks from the repo root against a normally installed
ChronoPack (a non-editable wheel in the ``chronopack-install-test`` env), one
stage at a time, like ``submit_pipeline_pair_staged.py``.

Usage::

    python submit_pipeline_user_api.py setup
    python submit_pipeline_user_api.py stage PreprocessDatasets [--experiment NAME]
    ... GetTrajectories, PrepareForTradeSeq, TradeSeqFits,
        DimensionalityReductionClustering
"""

from __future__ import annotations

import argparse
import hashlib
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from submit_full_pipeline_pair_20260924 import (  # noqa: E402
    EXPERIMENTS,
    ROOT,
    VERIFY_ROOT,
    rsync_template,
    sbatch,
)

RUN_NAME = "full_pipeline_maintained_chronopack_20260925_r8"
PREFIX = Path("/new-stg/home/alex/miniconda3/envs/chronopack-install-test")
GIT = "/new-stg/home/alex/miniconda3/envs/metacellR/bin/git"
WRAPPER = "run_notebook_installed_chronopack.sbatch"
SUCCESS_MARKER = "CHRONOPACK_STAGE_OK"

# Notebooks copied from the repo root into each experiment's notebooks/ dir.
NOTEBOOKS = {
    "k562-tnfa.hg38.rep": [
        "PreprocessDatasets.ipynb", "GetTrajectories.ipynb", "PrepareForTradeSeq.ipynb",
        "TradeSeqFits.ipynb", "exclude_leiden_clusters.manuscript.csv",
    ],
    "k562-tnfa-puromycin.hg38": [
        "PreprocessDatasets.ipynb", "GetTrajectories.ipynb", "PrepareForTradeSeq.ipynb",
        "TradeSeqFits.ipynb", "DimensionalityReduction&Clustering.self_contained.ipynb",
        "DimensionalityReduction&Clustering.cluster_visualization.ipynb",
    ],
    "thp1-k562-lps.hg38.rep": [
        "PreprocessDatasets.ipynb", "GetTrajectories.ipynb", "PrepareForTradeSeq.ipynb",
        "TradeSeqFits.K562.ipynb", "TradeSeqFits.THP1.ipynb",
    ],
}

STAGES = {
    "PreprocessDatasets": [("PreprocessDatasets", e, "PreprocessDatasets.ipynb") for e in EXPERIMENTS],
    "GetTrajectories": [("GetTrajectories", e, "GetTrajectories.ipynb") for e in EXPERIMENTS],
    "PrepareForTradeSeq": [("PrepareForTradeSeq", e, "PrepareForTradeSeq.ipynb") for e in EXPERIMENTS],
    "TradeSeqFits": [
        ("TradeSeqFits", "k562-tnfa.hg38.rep", "TradeSeqFits.ipynb"),
        ("TradeSeqFits", "k562-tnfa-puromycin.hg38", "TradeSeqFits.ipynb"),
        ("TradeSeqFits.K562", "thp1-k562-lps.hg38.rep", "TradeSeqFits.K562.ipynb"),
        ("TradeSeqFits.THP1", "thp1-k562-lps.hg38.rep", "TradeSeqFits.THP1.ipynb"),
    ],
    "DimensionalityReductionClustering": [
        ("DimensionalityReductionClustering", "k562-tnfa-puromycin.hg38",
         "DimensionalityReduction&Clustering.self_contained.ipynb"),
    ],
}


def write_wrapper(run_root: Path) -> None:
    logs = run_root / "logs"
    logs.mkdir(exist_ok=True)
    (run_root / "scripts" / WRAPPER).write_text(f"""#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output={logs}/%x.%j.out
#SBATCH --error={logs}/%x.%j.err

set -euo pipefail
: "${{NOTEBOOK:?NOTEBOOK must be set}}"
: "${{WORKDIR:?WORKDIR must be set}}"
: "${{RUN_ROOT:?RUN_ROOT must be set}}"

PREFIX="{PREFIX}"
export PATH="$PREFIX/bin:$PATH"
export CONDA_PREFIX="$PREFIX"
export PYTHONNOUSERSITE=1
unset PYTHONPATH
export MPLCONFIGDIR="$RUN_ROOT/tmp/mpl-${{SLURM_JOB_ID:-manual}}"
export TMPDIR="$RUN_ROOT/tmp/tmp-${{SLURM_JOB_ID:-manual}}"
mkdir -p "$MPLCONFIGDIR" "$TMPDIR"

cd "$WORKDIR"
# ChronoPack must come from the installed wheel, never from the source tree.
python -c 'import sys, chronopack; f = chronopack.__file__; print("chronopack from", f); sys.exit(0 if "/site-packages/" in f else 1)'
python "$RUN_ROOT/scripts/run_python_notebook_with_display.py" "$NOTEBOOK"
echo "{SUCCESS_MARKER} job=${{SLURM_JOB_ID:-manual}}"
""")


def write_provenance(run_root: Path) -> None:
    def git(*args):
        return subprocess.run([GIT, *args], cwd=ROOT, text=True, capture_output=True, check=True).stdout
    wheels = sorted((VERIFY_ROOT / "install_test_20260925" / "dist").glob("chronopack-*.whl"))
    lines = [f"run: {RUN_NAME}", f"env: {PREFIX}", f"git HEAD: {git('rev-parse', 'HEAD').strip()}", ""]
    for w in wheels:
        lines.append(f"wheel: {w} sha256={hashlib.sha256(w.read_bytes()).hexdigest()}")
    lines += ["", "uncommitted changes (git status --short):", git("status", "--short")]
    (run_root / "code_provenance.txt").write_text("\n".join(lines))
    (run_root / "code_provenance.diff").write_text(git("diff", "HEAD"))


def setup() -> None:
    run_root = VERIFY_ROOT / RUN_NAME
    print(f"Creating {run_root} with env {PREFIX}")
    rsync_template(run_root)
    for exp, files in NOTEBOOKS.items():
        dst = run_root / exp / "notebooks"
        for name in files:
            shutil.copy2(ROOT / exp / name, dst / name)
    write_wrapper(run_root)
    write_provenance(run_root)


def submit_stage(stage: str, experiments: list[str] | None) -> None:
    run_root = VERIFY_ROOT / RUN_NAME
    if not (run_root / "scripts" / WRAPPER).exists():
        raise SystemExit(f"{run_root} is not set up; run `setup` first.")
    for label, experiment, notebook in STAGES[stage]:
        if experiments and experiment not in experiments:
            continue
        job_id = sbatch(run_root, stage=label, experiment=experiment, notebook=notebook, wrapper=WRAPPER)
        print(f"{RUN_NAME}\t{label}\t{experiment}\t{job_id}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("setup")
    stage_parser = sub.add_parser("stage")
    stage_parser.add_argument("name", choices=list(STAGES))
    stage_parser.add_argument("--experiment", action="append", choices=EXPERIMENTS)
    args = parser.parse_args()
    if args.command == "setup":
        setup()
    else:
        submit_stage(args.name, args.experiment)


if __name__ == "__main__":
    main()
