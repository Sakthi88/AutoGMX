#!/usr/bin/env bash
set -Eeuo pipefail

# Portable HPC wrapper for the MD pipeline.
# This script does not submit a job by itself; call it from a scheduler script
# after loading the required modules/conda environment.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

JSON_CONFIG=""
START_STEP=""
END_STEP=""

usage() {
  cat <<'EOF'
Usage: bash hpc/run_hpc_job.sh [--json config.json] [--from N] [--to N]

Environment variables commonly set by scheduler templates:
  GMX_BIN        GROMACS binary, e.g. gmx or gmx_mpi
  PYTHON_BIN     Python binary, default python3
  MPI_LAUNCHER   Optional launcher, e.g. "srun" or "mpirun -np 4"
  NTMPI          GROMACS thread-MPI ranks. Use 0 when using external MPI gmx_mpi.
  NTOMP          OpenMP threads per rank.
  USE_GPU        auto, yes, or no
  GPU_ID         Optional GPU id
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --json)
      JSON_CONFIG="$2"
      shift 2
      ;;
    --from)
      START_STEP="$2"
      shift 2
      ;;
    --to)
      END_STEP="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

export PYTHON_BIN="${PYTHON_BIN:-python3}"
export GMX_BIN="${GMX_BIN:-gmx}"
export USE_GPU="${USE_GPU:-auto}"

# Infer scheduler resources when the submit script has not set them explicitly.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  export NTOMP="${NTOMP:-${SLURM_CPUS_PER_TASK:-auto}}"
  if [[ -z "${MPI_LAUNCHER:-}" && "${GMX_BIN}" == *"gmx_mpi"* ]]; then
    export MPI_LAUNCHER="srun"
    export NTMPI="${NTMPI:-0}"
  else
    export NTMPI="${NTMPI:-${SLURM_NTASKS:-1}}"
  fi
elif [[ -n "${PBS_JOBID:-}" ]]; then
  export NTOMP="${NTOMP:-${OMP_NUM_THREADS:-auto}}"
  if [[ -z "${MPI_LAUNCHER:-}" && "${GMX_BIN}" == *"gmx_mpi"* ]]; then
    export MPI_LAUNCHER="mpirun -np ${PBS_NP:-1}"
    export NTMPI="${NTMPI:-0}"
  else
    export NTMPI="${NTMPI:-1}"
  fi
else
  export NTOMP="${NTOMP:-auto}"
  export NTMPI="${NTMPI:-1}"
fi

echo "HPC wrapper summary"
echo "  Pipeline: ${PIPELINE_DIR}"
echo "  Host: $(hostname)"
echo "  Date: $(date '+%F %T')"
echo "  GMX_BIN: ${GMX_BIN}"
echo "  PYTHON_BIN: ${PYTHON_BIN}"
echo "  MPI_LAUNCHER: ${MPI_LAUNCHER:-<none>}"
echo "  NTMPI: ${NTMPI}"
echo "  NTOMP: ${NTOMP}"
echo "  USE_GPU: ${USE_GPU}"

cmd=(bash "${PIPELINE_DIR}/run_pipeline.sh")
if [[ -n "${JSON_CONFIG}" ]]; then
  cmd+=(--json "${JSON_CONFIG}")
fi
if [[ -n "${START_STEP}" ]]; then
  cmd+=(--from "${START_STEP}")
fi
if [[ -n "${END_STEP}" ]]; then
  cmd+=(--to "${END_STEP}")
fi

echo "Command: ${cmd[*]}"
exec "${cmd[@]}"
