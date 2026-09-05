#!/usr/bin/env bash
#PBS -N autogmx
#PBS -l select=1:ncpus=8:mem=32gb
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -o logs/pbs.out

set -Eeuo pipefail

# Submit from the project root:
#   qsub hpc/pbs.submit.sh
#
# Edit the module/conda lines for your cluster before submission.

cd "${PBS_O_WORKDIR}"
mkdir -p logs

# Example environment setup. Replace these with site-specific module names.
# module purge
# module load gromacs/2024.3
# module load ambertools
# module load apbs
# module load pdb2pqr
# source ~/miniconda3/etc/profile.d/conda.sh
# conda activate md_pipeline

export PYTHON_BIN="${PYTHON_BIN:-python3}"
export GMX_BIN="${GMX_BIN:-gmx}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-8}"
export NTOMP="${NTOMP:-${OMP_NUM_THREADS}}"
export NTMPI="${NTMPI:-1}"
export USE_GPU="${USE_GPU:-auto}"

# For external-MPI GROMACS builds, use something like:
#   export GMX_BIN=gmx_mpi
#   export MPI_LAUNCHER="mpirun -np ${PBS_NP:-1}"
#   export NTMPI=0

bash hpc/run_hpc_job.sh --json hpc/cluster.example.json --to 11
