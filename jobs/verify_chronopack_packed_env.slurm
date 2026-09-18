#!/usr/bin/env bash
#SBATCH --job-name=chronopack-env
#SBATCH --output=jobs/logs/chronopack-env-%j.out
#SBATCH --error=jobs/logs/chronopack-env-%j.err
#SBATCH --time=00:45:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

set -euo pipefail

cd /new-stg/home/alex/tmp/ChronoPack

if [[ -n "${SLURM_JOB_ID:-}" && -n "${USER:-}" && -d "/scratch/${USER}" && -w "/scratch/${USER}" ]]; then
  export CHRONOPACK_LOCAL="/scratch/${USER}/chronopack-${SLURM_JOB_ID}"
else
  export CHRONOPACK_LOCAL="/tmp/chronopack-${USER:-alex}-${SLURM_JOB_ID:-local}"
fi

export TMPDIR="${CHRONOPACK_LOCAL}/tmp"
export CONDA_PKGS_DIRS="${CHRONOPACK_LOCAL}/conda-pkgs"
export MPLCONFIGDIR="${TMPDIR}/matplotlib"

source env/activate_chronopack_hpc.sh
bash env/smoke_test_chronopack_hpc.sh
