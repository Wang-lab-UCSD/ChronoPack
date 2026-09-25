#!/usr/bin/env python
"""Submit paired full-pipeline verification runs.

Creates two fresh verification directories from the current manuscript pipeline
template and submits the full notebook dependency graph:

1. frozen manuscript environment
2. maintained ChronoPack environment
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path


ROOT = Path("/new-stg/home/alex/tmp/ChronoPack")
VERIFY_ROOT = ROOT / "verification_runs"
TEMPLATE = VERIFY_ROOT / "manuscript_full_pipeline_20260919_165539"

RUNS = {
    "full_pipeline_frozen_manuscript_20260924_r5": Path(
        "/new-stg/home/alex/miniconda3/envs/chronopack-manuscript-20260919"
    ),
    "full_pipeline_maintained_chronopack_20260924_r6": Path(
        "/new-stg/home/alex/miniconda3/envs/chronopack"
    ),
}

EXPERIMENTS = [
    "k562-tnfa.hg38.rep",
    "k562-tnfa-puromycin.hg38",
    "thp1-k562-lps.hg38.rep",
]


def run(cmd: list[str], *, cwd: Path | None = None, env: dict[str, str] | None = None) -> str:
    out = subprocess.check_output(cmd, cwd=str(cwd or ROOT), env=env, text=True)
    return out.strip()


def rsync_template(run_root: Path) -> None:
    if run_root.exists():
        raise SystemExit(f"Refusing to overwrite existing run directory: {run_root}")
    run_root.mkdir(parents=True)

    for name in EXPERIMENTS:
        src = TEMPLATE / name
        dst = run_root / name
        cmd = [
            "rsync",
            "-a",
            "--exclude=output/",
            "--exclude=intmd/",
            "--exclude=self_contained_outputs/",
            "--exclude=notebooks/.ipynb_checkpoints/",
            f"{src}/",
            f"{dst}/",
        ]
        run(cmd)
        (dst / "output").mkdir(exist_ok=True)
        (dst / "intmd").mkdir(exist_ok=True)

    neg_ctrl = TEMPLATE / "thp1-k562-lps.hg38.rep" / "intmd" / "k562_neg-ctrl_mtx-agg.pkl"
    if neg_ctrl.exists():
        dst = run_root / "thp1-k562-lps.hg38.rep" / "intmd" / neg_ctrl.name
        run(["rsync", "-a", str(neg_ctrl), str(dst)])

    for name in ["scripts", "extracted_r", "logs", "tmp"]:
        (run_root / name).mkdir(exist_ok=True)
    run(["rsync", "-a", f"{TEMPLATE / 'scripts'}/", f"{run_root / 'scripts'}/"])
    (run_root / "job_ids.tsv").write_text("stage\texperiment\tjob_id\tdependency\n")


def write_wrappers(run_root: Path, prefix: Path) -> None:
    scripts = run_root / "scripts"
    logs = run_root / "logs"
    logs.mkdir(exist_ok=True)

    common = f"""set -euo pipefail

: "${{NOTEBOOK:?NOTEBOOK must be set}}"
: "${{WORKDIR:?WORKDIR must be set}}"
: "${{RUN_ROOT:?RUN_ROOT must be set}}"

PREFIX="{prefix}"
export PATH="$PREFIX/bin:$PATH"
export CONDA_PREFIX="$PREFIX"
export PYTHONNOUSERSITE=1
export MPLCONFIGDIR="$RUN_ROOT/tmp/mpl-${{SLURM_JOB_ID:-manual}}"
export TMPDIR="$RUN_ROOT/tmp/tmp-${{SLURM_JOB_ID:-manual}}"
mkdir -p "$MPLCONFIGDIR" "$TMPDIR" "$RUN_ROOT/extracted_r"
"""

    (scripts / "run_python_notebook.env.sbatch").write_text(
        f"""#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=48:00:00
#SBATCH --output={logs}/%x.%j.out
#SBATCH --error={logs}/%x.%j.err

{common}
cd "$WORKDIR"
python "$RUN_ROOT/scripts/run_python_notebook.py" "$NOTEBOOK"
"""
    )
    (scripts / "run_r_notebook.env.sbatch").write_text(
        f"""#!/bin/bash
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=96:00:00
#SBATCH --output={logs}/%x.%j.out
#SBATCH --error={logs}/%x.%j.err

{common}
SCRIPT="$RUN_ROOT/extracted_r/$(basename "$NOTEBOOK" .ipynb).${{SLURM_JOB_ID:-manual}}.R"
python "$RUN_ROOT/scripts/extract_r_notebook.py" "$NOTEBOOK" "$SCRIPT"
cd "$WORKDIR"
Rscript "$SCRIPT"
"""
    )


def sbatch(
    run_root: Path,
    *,
    stage: str,
    experiment: str,
    notebook: str,
    wrapper: str,
    dependency: str | None = None,
) -> str:
    env = os.environ.copy()
    exp_root = run_root / experiment
    notebook_path = exp_root / "notebooks" / notebook
    env.update(
        {
            "RUN_ROOT": str(run_root),
            "WORKDIR": str(exp_root / "notebooks"),
            "NOTEBOOK": str(notebook_path),
        }
    )
    cmd = [
        "sbatch",
        f"--job-name=fp-{run_root.name[:10]}-{stage[:8]}-{experiment[:8]}",
    ]
    if dependency:
        cmd.append(f"--dependency=afterok:{dependency}")
    cmd.append(str(run_root / "scripts" / wrapper))
    out = run(cmd, env=env)
    job_id = out.split()[-1]
    with (run_root / "job_ids.tsv").open("a") as handle:
        handle.write(f"{stage}\t{experiment}\t{job_id}\t{dependency or ''}\n")
    return job_id


def submit_graph(run_root: Path) -> None:
    jobs: dict[tuple[str, str], str] = {}
    for exp in EXPERIMENTS:
        jobs[("PreprocessDatasets", exp)] = sbatch(
            run_root,
            stage="PreprocessDatasets",
            experiment=exp,
            notebook="PreprocessDatasets.ipynb",
            wrapper="run_python_notebook.env.sbatch",
        )
    for exp in EXPERIMENTS:
        jobs[("GetTrajectories", exp)] = sbatch(
            run_root,
            stage="GetTrajectories",
            experiment=exp,
            notebook="GetTrajectories.ipynb",
            wrapper="run_python_notebook.env.sbatch",
            dependency=jobs[("PreprocessDatasets", exp)],
        )
    for exp in EXPERIMENTS:
        jobs[("PrepareForTradeSeq", exp)] = sbatch(
            run_root,
            stage="PrepareForTradeSeq",
            experiment=exp,
            notebook="PrepareForTradeSeq.ipynb",
            wrapper="run_python_notebook.env.sbatch",
            dependency=jobs[("GetTrajectories", exp)],
        )

    jobs[("TradeSeqFits", "k562-tnfa.hg38.rep")] = sbatch(
        run_root,
        stage="TradeSeqFits",
        experiment="k562-tnfa.hg38.rep",
        notebook="TradeSeqFits.ipynb",
        wrapper="run_r_notebook.env.sbatch",
        dependency=jobs[("PrepareForTradeSeq", "k562-tnfa.hg38.rep")],
    )
    jobs[("TradeSeqFits", "k562-tnfa-puromycin.hg38")] = sbatch(
        run_root,
        stage="TradeSeqFits",
        experiment="k562-tnfa-puromycin.hg38",
        notebook="TradeSeqFits.ipynb",
        wrapper="run_r_notebook.env.sbatch",
        dependency=jobs[("PrepareForTradeSeq", "k562-tnfa-puromycin.hg38")],
    )
    jobs[("TradeSeqFits.K562", "thp1-k562-lps.hg38.rep")] = sbatch(
        run_root,
        stage="TradeSeqFits.K562",
        experiment="thp1-k562-lps.hg38.rep",
        notebook="TradeSeqFits.K562.ipynb",
        wrapper="run_r_notebook.env.sbatch",
        dependency=jobs[("PrepareForTradeSeq", "thp1-k562-lps.hg38.rep")],
    )
    jobs[("TradeSeqFits.THP1", "thp1-k562-lps.hg38.rep")] = sbatch(
        run_root,
        stage="TradeSeqFits.THP1",
        experiment="thp1-k562-lps.hg38.rep",
        notebook="TradeSeqFits.THP1.ipynb",
        wrapper="run_r_notebook.env.sbatch",
        dependency=jobs[("PrepareForTradeSeq", "thp1-k562-lps.hg38.rep")],
    )

    dimred_dep = ":".join(
        [
            jobs[("TradeSeqFits", "k562-tnfa.hg38.rep")],
            jobs[("TradeSeqFits", "k562-tnfa-puromycin.hg38")],
            jobs[("TradeSeqFits.K562", "thp1-k562-lps.hg38.rep")],
            jobs[("TradeSeqFits.THP1", "thp1-k562-lps.hg38.rep")],
        ]
    )
    sbatch(
        run_root,
        stage="DimensionalityReductionClustering",
        experiment="k562-tnfa-puromycin.hg38",
        notebook="DimensionalityReduction&Clustering.self_contained.ipynb",
        wrapper="run_python_notebook.env.sbatch",
        dependency=dimred_dep,
    )


def main() -> None:
    for run_name, prefix in RUNS.items():
        run_root = VERIFY_ROOT / run_name
        print(f"Creating {run_root} with env {prefix}")
        rsync_template(run_root)
        write_wrappers(run_root, prefix)
        submit_graph(run_root)
        print((run_root / "job_ids.tsv").read_text())


if __name__ == "__main__":
    main()
