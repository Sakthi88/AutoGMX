#!/usr/bin/env bash
#SBATCH --job-name=autogmx
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err

set -Eeuo pipefail

# Submit from the project root:
#   sbatch hpc/slurm.submit.sh
#
# Edit the module/conda lines for your cluster before submission.

cd "${SLURM_SUBMIT_DIR}"
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
export NTOMP="${SLURM_CPUS_PER_TASK:-8}"
export NTMPI="${NTMPI:-1}"
export USE_GPU="${USE_GPU:-auto}"

# For external-MPI GROMACS builds, use something like:
#   export GMX_BIN=gmx_mpi
#   export MPI_LAUNCHER=srun
#   export NTMPI=0

bash hpc/run_hpc_job.sh --json hpc/cluster.example.json --to 11
