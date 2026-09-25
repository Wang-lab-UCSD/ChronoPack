#!/usr/bin/env python
"""Stage-by-stage paired pipeline verification.

Unlike ``submit_full_pipeline_pair_20260924.py``, which submits the whole
dependency graph at once, this script submits one stage at a time (for both
environments) with no Slurm dependencies. The next stage is submitted only
after the previous one has been checked.

Usage::

    python submit_pipeline_pair_staged.py setup
    python submit_pipeline_pair_staged.py stage PreprocessDatasets
    python submit_pipeline_pair_staged.py stage GetTrajectories
    python submit_pipeline_pair_staged.py stage PrepareForTradeSeq
    python submit_pipeline_pair_staged.py stage TradeSeqFits
    python submit_pipeline_pair_staged.py stage DimensionalityReductionClustering
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from submit_full_pipeline_pair_20260924 import (  # noqa: E402
    EXPERIMENTS,
    VERIFY_ROOT,
    rsync_template,
    sbatch,
    write_wrappers,
)


RUNS = {
    "full_pipeline_frozen_manuscript_20260924_r6": Path(
        "/new-stg/home/alex/miniconda3/envs/chronopack-manuscript-20260919"
    ),
    "full_pipeline_maintained_chronopack_20260924_r7": Path(
        "/new-stg/home/alex/miniconda3/envs/chronopack"
    ),
}

# Per-run files copied over the template after setup. Only the maintained run
# gets the barcode-based Leiden cluster removal for K562+TNF preprocessing; the
# frozen run keeps the manuscript's label-based removal unchanged.
OVERRIDES = {
    "full_pipeline_maintained_chronopack_20260924_r7": VERIFY_ROOT / "maintained_overrides",
}

PY = "run_python_notebook.env.sbatch"
R = "run_r_notebook.env.sbatch"
# DimRed uses IPython's display() and imports chronopack without adding the
# repo to sys.path, so it runs like the earlier DimRed launchers
# (scripts/rerun_dimred_only.sbatch): the display-aware runner plus src on
# PYTHONPATH. The maintained env has no editable ChronoPack install.
PY_DISPLAY = "run_python_notebook_display.env.sbatch"

# Slurm accounting is disabled on this cluster, so sacct cannot report job
# outcomes. The wrappers run under `set -e`, so this line is only reached when
# the notebook finished successfully.
SUCCESS_MARKER = "CHRONOPACK_STAGE_OK"

# stage -> [(job stage label, experiment, notebook, wrapper)]
STAGES = {
    "PreprocessDatasets": [
        ("PreprocessDatasets", exp, "PreprocessDatasets.ipynb", PY) for exp in EXPERIMENTS
    ],
    "GetTrajectories": [
        ("GetTrajectories", exp, "GetTrajectories.ipynb", PY) for exp in EXPERIMENTS
    ],
    "PrepareForTradeSeq": [
        ("PrepareForTradeSeq", exp, "PrepareForTradeSeq.ipynb", PY) for exp in EXPERIMENTS
    ],
    "TradeSeqFits": [
        ("TradeSeqFits", "k562-tnfa.hg38.rep", "TradeSeqFits.ipynb", R),
        ("TradeSeqFits", "k562-tnfa-puromycin.hg38", "TradeSeqFits.ipynb", R),
        ("TradeSeqFits.K562", "thp1-k562-lps.hg38.rep", "TradeSeqFits.K562.ipynb", R),
        ("TradeSeqFits.THP1", "thp1-k562-lps.hg38.rep", "TradeSeqFits.THP1.ipynb", R),
    ],
    "DimensionalityReductionClustering": [
        (
            "DimensionalityReductionClustering",
            "k562-tnfa-puromycin.hg38",
            "DimensionalityReduction&Clustering.self_contained.ipynb",
            PY_DISPLAY,
        ),
    ],
}


def setup() -> None:
    for run_name, prefix in RUNS.items():
        run_root = VERIFY_ROOT / run_name
        print(f"Creating {run_root} with env {prefix}")
        rsync_template(run_root)
        write_wrappers(run_root, prefix)
        for wrapper in (PY, R):
            path = run_root / "scripts" / wrapper
            path.write_text(path.read_text() + f'echo "{SUCCESS_MARKER} job=${{SLURM_JOB_ID:-manual}}"\n')
        write_display_wrapper(run_root)
        apply_overrides(run_name)


def write_display_wrapper(run_root: Path) -> None:
    scripts = run_root / "scripts"
    text = (scripts / PY).read_text()
    old = 'python "$RUN_ROOT/scripts/run_python_notebook.py" "$NOTEBOOK"'
    if old not in text:
        raise SystemExit(f"Unexpected wrapper contents in {scripts / PY}")
    new = (
        'export PYTHONPATH="/new-stg/home/alex/tmp/ChronoPack/src${PYTHONPATH:+:$PYTHONPATH}"\n'
        'python "$RUN_ROOT/scripts/run_python_notebook_with_display.py" "$NOTEBOOK"'
    )
    (scripts / PY_DISPLAY).write_text(text.replace(old, new))


def apply_overrides(run_name: str) -> None:
    overrides = OVERRIDES.get(run_name)
    if overrides is None:
        return
    run_root = VERIFY_ROOT / run_name
    print(f"Applying {overrides} to {run_root}")
    subprocess.check_call(["rsync", "-a", "--itemize-changes", f"{overrides}/", f"{run_root}/"])


def submit_stage(stage: str, runs: list[str], experiments: list[str] | None) -> None:
    for run_name in runs:
        run_root = VERIFY_ROOT / run_name
        if not (run_root / "scripts" / PY).exists():
            raise SystemExit(f"{run_root} is not set up; run `setup` first.")
        if not (run_root / "scripts" / PY_DISPLAY).exists():
            write_display_wrapper(run_root)
        for label, experiment, notebook, wrapper in STAGES[stage]:
            if experiments and experiment not in experiments:
                continue
            job_id = sbatch(
                run_root,
                stage=label,
                experiment=experiment,
                notebook=notebook,
                wrapper=wrapper,
            )
            print(f"{run_name}\t{label}\t{experiment}\t{job_id}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("setup")
    overrides_parser = sub.add_parser("apply-overrides")
    overrides_parser.add_argument("run", choices=list(OVERRIDES))
    stage_parser = sub.add_parser("stage")
    stage_parser.add_argument("name", choices=list(STAGES))
    stage_parser.add_argument("--run", action="append", choices=list(RUNS),
                              help="limit to these runs (default: both)")
    stage_parser.add_argument("--experiment", action="append", choices=EXPERIMENTS,
                              help="limit to these experiments (default: all)")
    args = parser.parse_args()
    if args.command == "setup":
        setup()
    elif args.command == "apply-overrides":
        apply_overrides(args.run)
    else:
        submit_stage(args.name, args.run or list(RUNS), args.experiment)


if __name__ == "__main__":
    main()
