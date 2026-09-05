# HPC Submission Layer

This directory adds scheduler wrappers around the existing AutoGMX pipeline.
It does not replace `run_pipeline.sh`; it only prepares the scheduler
environment and then calls the normal pipeline.

## Files

- `slurm.submit.sh`: Slurm batch template.
- `pbs.submit.sh`: PBS/Torque batch template.
- `cluster.example.json`: HPC-oriented pipeline configuration.
- `run_hpc_job.sh`: Portable wrapper called by scheduler scripts.

Scheduler resource variables such as `GMX_BIN`, `MPI_LAUNCHER`, `NTMPI`,
`NTOMP`, and `USE_GPU` are intentionally kept in the submit scripts instead of
`cluster.example.json`. JSON overrides are loaded after `config.sh`, so values
inside the JSON file would override the scheduler allocation.

## Slurm

Edit module and conda lines in `hpc/slurm.submit.sh`, then submit from the
project root:

```bash
sbatch hpc/slurm.submit.sh
```

For external-MPI GROMACS builds, set:

```bash
export GMX_BIN=gmx_mpi
export MPI_LAUNCHER=srun
export NTMPI=0
```

Use `NTMPI=0` because external-MPI `gmx_mpi` normally should not receive the
thread-MPI `-ntmpi` option.

## PBS/Torque

Edit module and conda lines in `hpc/pbs.submit.sh`, then submit from the
project root:

```bash
qsub hpc/pbs.submit.sh
```

For external-MPI GROMACS builds, set:

```bash
export GMX_BIN=gmx_mpi
export MPI_LAUNCHER="mpirun -np ${PBS_NP:-1}"
export NTMPI=0
```

## Direct Allocation

Inside an interactive allocation, after loading modules:

```bash
bash hpc/run_hpc_job.sh --json hpc/cluster.example.json --to 11
```

Restart from a failed stage:

```bash
bash hpc/run_hpc_job.sh --json hpc/cluster.example.json --from 9 --to 11
```

## Notes

- Do not run production MD on a login node.
- Keep input files in `inputs/` unless paths are changed in the JSON file.
- Use cluster-specific module names; the template module lines are comments.
- For multi-cluster use, copy `cluster.example.json` to one config per cluster,
  for example `cluster.delta.json`, `cluster.gpu.json`, or `cluster.cpu.json`.
