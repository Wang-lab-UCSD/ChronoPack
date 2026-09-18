#!/usr/bin/env bash
# Source this file from an interactive Slurm allocation or a Slurm job.
#
# It hydrates the packed ChronoPack dependency environment onto node-local
# storage, runs conda-unpack, activates the environment, and installs the
# current checkout in editable mode.

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  echo "This script must be sourced so it can activate the environment:" >&2
  echo "  source env/activate_chronopack_hpc.sh" >&2
  exit 2
fi

_chronopack_script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" || return 1
_chronopack_repo="${CHRONOPACK_REPO:-$(cd "${_chronopack_script_dir}/.." && pwd)}"
_chronopack_archive="${CHRONOPACK_ENV_ARCHIVE:-${_chronopack_repo}/artifacts/env/chronopack-linux-64-20260918.tar.gz}"

if [[ ! -f "${_chronopack_archive}" ]]; then
  echo "ChronoPack environment archive is missing:" >&2
  echo "  ${_chronopack_archive}" >&2
  echo "Set CHRONOPACK_ENV_ARCHIVE to a valid conda-pack archive path." >&2
  return 1
fi

_chronopack_job_id="${SLURM_JOB_ID:-local}"
if [[ -n "${CHRONOPACK_LOCAL:-}" ]]; then
  _chronopack_local="${CHRONOPACK_LOCAL}"
elif [[ -n "${USER:-}" && -d "/scratch/${USER}" && -w "/scratch/${USER}" ]]; then
  _chronopack_local="/scratch/${USER}/chronopack-${_chronopack_job_id}"
else
  _chronopack_user="${USER:-$(id -un)}"
  _chronopack_local="/tmp/chronopack-${_chronopack_user}-${_chronopack_job_id}"
fi

export CHRONOPACK_LOCAL="${_chronopack_local}"
export TMPDIR="${TMPDIR:-${CHRONOPACK_LOCAL}/tmp}"
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-${CHRONOPACK_LOCAL}/conda-pkgs}"
export MPLCONFIGDIR="${MPLCONFIGDIR:-${TMPDIR}/matplotlib}"
export CHRONOPACK_ENV_PREFIX="${CHRONOPACK_ENV_PREFIX:-${CHRONOPACK_LOCAL}/env}"

mkdir -p "${TMPDIR}" "${CONDA_PKGS_DIRS}" "${MPLCONFIGDIR}" "${CHRONOPACK_ENV_PREFIX}" || return 1

if [[ ! -x "${CHRONOPACK_ENV_PREFIX}/bin/python" ]]; then
  echo "Unpacking ChronoPack environment archive to ${CHRONOPACK_ENV_PREFIX}"
  tar -xzf "${_chronopack_archive}" -C "${CHRONOPACK_ENV_PREFIX}" || return 1
fi

if [[ -x "${CHRONOPACK_ENV_PREFIX}/bin/conda-unpack" ]]; then
  "${CHRONOPACK_ENV_PREFIX}/bin/conda-unpack" || return 1
else
  echo "conda-unpack was not found in ${CHRONOPACK_ENV_PREFIX}/bin" >&2
  return 1
fi

source "${CHRONOPACK_ENV_PREFIX}/bin/activate" || return 1

if [[ "${CHRONOPACK_SKIP_EDITABLE:-0}" != "1" ]]; then
  python -m pip install -e "${_chronopack_repo}" || return 1
fi

echo "ChronoPack environment active:"
echo "  repo:   ${_chronopack_repo}"
echo "  prefix: ${CHRONOPACK_ENV_PREFIX}"
echo "  local:  ${CHRONOPACK_LOCAL}"
