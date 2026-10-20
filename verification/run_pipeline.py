#!/usr/bin/env python
"""Run the analyses/ notebooks as a pipeline verification, one stage at a time.

    python verification/run_pipeline.py setup RUN_NAME
    python verification/run_pipeline.py stage STAGE RUN_NAME [--experiment NAME] [--notebook PATH]
                                              [--mem 90G] [--allow-bigcompute]
    python verification/run_pipeline.py promote RUN_NAME

``setup`` creates ``verification_runs/RUN_NAME/`` with an empty data folder
(raw inputs are read from ``analyses/inputs/``) and records code provenance.
``stage`` submits one Slurm job per notebook of a stage. Each job executes the
notebook straight from ``analyses/`` (``CHRONOPACK_DATA_DIR`` pointed at the
run's data folder) and saves the executed copy, with outputs, under
``RUN_NAME/executed/``. A job succeeded when its log ends with SUCCESS_MARKER
(Slurm accounting is disabled on this cluster). ``promote`` copies the executed
notebooks back into ``analyses/`` once a run has been checked.

ChronoPack is imported from a normal (non-editable) wheel install in ENV, never
from the source tree; rebuild and reinstall the wheel before ``setup`` after
package edits.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ANALYSES = ROOT / "analyses"
RUNS = ROOT / "verification_runs"
ENV = Path("/new-stg/home/alex/miniconda3/envs/chronopack-install-test")
GIT = "/new-stg/home/alex/miniconda3/envs/metacellR/bin/git"
SUCCESS_MARKER = "CHRONOPACK_STAGE_OK"

EXPERIMENTS = ["k562-tnfa.hg38.rep", "thp1-k562-lps.hg38.rep"]

K, T = EXPERIMENTS
# stage -> notebooks (relative to analyses/). Notebooks listed in one tuple run
# sequentially in one job; separate tuples are separate jobs.
STAGES = {
    "PreprocessDatasets": [(f"{K}/01_PreprocessDatasets.ipynb",), (f"{T}/01_PreprocessDatasets.ipynb",)],
    "GetTrajectories": [(f"{K}/02_GetTrajectories.ipynb",), (f"{T}/02_GetTrajectories.ipynb",)],
    "PrepareForTradeSeq": [(f"{K}/03_PrepareForTradeSeq.ipynb",), (f"{T}/03_PrepareForTradeSeq.ipynb",)],
    "TradeSeqFits": [
        (f"{K}/04_TradeSeqFits.ipynb",),
        (f"{T}/04_TradeSeqFits.K562.ipynb",),
        (f"{T}/04_TradeSeqFits.THP1.ipynb",),
    ],
    # Optional: the two-lineage K562 + TNF-alpha fit (not used downstream). Fits in the
    # default 90G on a regular node (~85 min, 2026-10-20). Avoid --allow-bigcompute for
    # tradeSeq stages: the all-cells K562 fit took ~3 h there vs 18 min on a regular node.
    "TradeSeqFitsTwoLineage": [(f"{K}/04_TradeSeqFits.TwoLineage.ipynb",)],
    # Optional: figures of both experiments (need GetTrajectories and TradeSeqFits).
    "Visualization": [(f"{K}/05_Visualization.ipynb",), (f"{T}/05_Visualization.ipynb",)],
    "GeneClustering": [("gene_clustering/01_GeneTrajectoryClustering.ipynb",)],
    # Optional: the same clustering with responsive genes = p < 0.05 only (no |peak log2FC| cut).
    "GeneClusteringPvalueOnly": [("gene_clustering/01_GeneTrajectoryClustering.pvalue_only.ipynb",)],
    "CellCommunication": [("c2c/CellCommunication.ipynb",)],
    "GRNInference": [("grn/GRNInference.ipynb",)],
}


def git(*args: str) -> str:
    return subprocess.run([GIT, *args], cwd=ROOT, text=True, capture_output=True, check=True).stdout


def setup(run: Path) -> None:
    if run.exists():
        raise SystemExit(f"Refusing to overwrite existing run directory: {run}")
    for sub in ("data", "executed", "logs", "tmp"):
        (run / sub).mkdir(parents=True)

    pip = subprocess.run([str(ENV / "bin/python"), "-m", "pip", "show", "-f", "chronopack"],
                         text=True, capture_output=True, check=True).stdout
    location = next(l.split(": ", 1)[1] for l in pip.splitlines() if l.startswith("Location: "))
    direct_url = Path(location).glob("chronopack-*.dist-info/direct_url.json")
    wheel = json.loads(next(direct_url).read_text())["url"].removeprefix("file://")
    lines = [
        f"run: {run.name}",
        f"env: {ENV}",
        f"wheel: {wheel} sha256={hashlib.sha256(Path(wheel).read_bytes()).hexdigest()}",
        f"git HEAD: {git('rev-parse', 'HEAD').strip()}",
        "",
        "uncommitted changes (git status --short):",
        git("status", "--short"),
    ]
    (run / "code_provenance.txt").write_text("\n".join(lines))
    (run / "code_provenance.diff").write_text(git("diff", "HEAD"))
    (run / "job_ids.tsv").write_text("stage\tnotebooks\tjob_id\n")
    print(f"Created {run}")


def job_script(run: Path, stage: str, notebooks: tuple[str, ...], mem: str, allow_bigcompute: bool) -> str:
    cmds = "\n".join(
        f'"$ENV/bin/python" "{ROOT}/verification/execute_notebook.py" "{ANALYSES / nb}" "{run / "executed" / nb}"'
        for nb in notebooks
    )
    return f"""#!/bin/bash
#SBATCH --job-name={stage}
#SBATCH --cpus-per-task=8
#SBATCH --mem={mem}
#SBATCH --time=48:00:00
{'' if allow_bigcompute else '#SBATCH --exclude=bigcompute'}
#SBATCH --output={run}/logs/%x.%j.out
#SBATCH --error={run}/logs/%x.%j.err
set -euo pipefail
ENV="{ENV}"
export PATH="$ENV/bin:$PATH" CONDA_PREFIX="$ENV" PYTHONNOUSERSITE=1 MPLBACKEND=Agg
unset PYTHONPATH
export CHRONOPACK_DATA_DIR="{run / 'data'}"
export MPLCONFIGDIR="{run}/tmp/mpl-$SLURM_JOB_ID" TMPDIR="{run}/tmp/tmp-$SLURM_JOB_ID"
mkdir -p "$MPLCONFIGDIR" "$TMPDIR"
cd "{run}"
# ChronoPack must come from the installed wheel, never from the source tree.
"$ENV/bin/python" -c 'import sys, chronopack; f = chronopack.__file__; print("chronopack from", f); sys.exit("/site-packages/" not in f)'
{cmds}
echo "{SUCCESS_MARKER} job=$SLURM_JOB_ID"
"""


def submit(run: Path, stage: str, experiment: str | None, mem: str, allow_bigcompute: bool, notebook: str | None) -> None:
    if not (run / "code_provenance.txt").exists():
        raise SystemExit(f"{run} is not set up; run `setup` first.")
    for notebooks in STAGES[stage]:
        if experiment and not notebooks[0].startswith(experiment + "/"):
            continue
        if notebook and notebook not in notebooks:
            continue
        script = run / "tmp" / f"{Path(notebooks[0]).parent.name}.{Path(notebooks[0]).stem}.sbatch"
        script.write_text(job_script(run, stage, notebooks, mem, allow_bigcompute))
        job_id = subprocess.run(["sbatch", "--parsable", str(script)], text=True,
                                capture_output=True, check=True).stdout.strip()
        with open(run / "job_ids.tsv", "a") as fh:
            fh.write(f"{stage}\t{','.join(notebooks)}\t{job_id}\n")
        print(f"{stage}\t{job_id}\t{', '.join(notebooks)}")


def promote(run: Path) -> None:
    executed = sorted((run / "executed").rglob("*.ipynb"))
    if not executed:
        raise SystemExit(f"No executed notebooks under {run / 'executed'}")
    for nb in executed:
        rel = nb.relative_to(run / "executed")
        dest = ANALYSES / rel
        src_cells = [c["source"] for c in json.loads(dest.read_text())["cells"]]
        exe_cells = [c["source"] for c in json.loads(nb.read_text())["cells"]]
        if src_cells != exe_cells:
            raise SystemExit(f"{rel}: source changed since {run.name} ran; rerun before promoting.")
        shutil.copy2(nb, dest)
        print(f"promoted {rel}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("setup").add_argument("run")
    stage = sub.add_parser("stage")
    stage.add_argument("name", choices=list(STAGES))
    stage.add_argument("run")
    stage.add_argument("--experiment", choices=EXPERIMENTS + ["gene_clustering", "c2c", "grn"])
    stage.add_argument("--notebook", help="only this notebook (path relative to analyses/)")
    # 96G+ requests only fit bigcompute on this cluster; the largest other nodes allow ~94G.
    stage.add_argument("--mem", default="90G")
    stage.add_argument("--allow-bigcompute", action="store_true")
    sub.add_parser("promote").add_argument("run")
    args = parser.parse_args()
    run = RUNS / args.run
    if args.command == "setup":
        setup(run)
    elif args.command == "stage":
        submit(run, args.name, args.experiment, args.mem, args.allow_bigcompute, args.notebook)
    else:
        promote(run)


if __name__ == "__main__":
    main()
