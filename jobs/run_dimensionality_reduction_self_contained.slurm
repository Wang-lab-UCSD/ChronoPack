#!/usr/bin/env bash
#SBATCH --job-name=chronopack_dimred
#SBATCH --output=/new-stg/home/alex/tmp/ChronoPack/k562-tnfa-puromycin.hg38/self_contained_outputs/slurm-%j.out
#SBATCH --error=/new-stg/home/alex/tmp/ChronoPack/k562-tnfa-puromycin.hg38/self_contained_outputs/slurm-%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00

set -euo pipefail

source /new-stg/home/alex/miniconda3/etc/profile.d/conda.sh
conda activate tradeseq_ce

WORKDIR=/new-stg/home/alex/tmp/ChronoPack/k562-tnfa-puromycin.hg38
OUTDIR="${WORKDIR}/self_contained_outputs/global_analysis_extended"

mkdir -p "${WORKDIR}/self_contained_outputs" "${OUTDIR}"

export CHRONOPACK_TRADESEQ_WORKERS="${SLURM_CPUS_PER_TASK:-4}"
export CHRONOPACK_DIMRED_OUTDIR="${OUTDIR}"
export CHRONOPACK_VALIDATE_REFERENCE=1

cd "${WORKDIR}"

jupyter nbconvert \
  --to notebook \
  --execute "DimensionalityReduction&Clustering.self_contained.ipynb" \
  --output "DimensionalityReduction&Clustering.self_contained.executed.ipynb" \
  --ExecutePreprocessor.timeout=-1 \
  --ExecutePreprocessor.kernel_name=ir
