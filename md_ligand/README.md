<<<<<<< HEAD
# Automated GROMACS Protein-Ligand MD Pipeline (MM-PBSA Version)

> **⚠️ IMPORTANT WARNINGS - READ BEFORE RUNNING**
>
> - **DO NOT run long production simulations on laptops or small systems** - Extended MD runs (100+ ns) generate excessive heat, can cause thermal throttling, hardware degradation, or permanent damage. Use this pipeline for **demo/short runs only** (1-10 ns) on personal machines. For production work, use HPC clusters, cloud instances, or dedicated GPU workstations with proper cooling.
> - **CHARMM force field NOT included in default GROMACS installer** - You must manually download and install the latest CHARMM36 force field from [CHARMM-GUI](https://www.charmm-gui.org/) or [MacKerell Lab](https://mackerell.umaryland.edu/charmm_ff.shtml) into the appropriate GROMACS force field directory (e.g., `$GMXDIR/share/gromacs/top/charmm36-feb2026.ff/`) without errors.
> - **This pipeline is built for GROMACS accessibility, NOT for lazy usage** - First learn the MD process manually (pdb2gmx, editconf, genbox, genion, grompp, mdrun, analysis). Understanding each step makes it easier to debug small errors when using this automation tool.
> - **gmxapi / gromacs_py Python package** - If you plan to use the Python API for GROMACS, install `gmxapi` or `gromacs_py` separately: `conda install -c conda-forge gmxapi` or `pip install gmxapi`. This pipeline uses CLI `gmx` commands by default; Python API is optional for custom extensions.

This pipeline automates a standard protein-ligand molecular dynamics workflow in GROMACS with **MM-PBSA binding free energy analysis** as the final step. It covers system preparation through production MD and MM-PBSA calculation (steps 1-11).

## Features

- **Modular 11-step pipeline** (`step1` to `step11`) with restart-safe execution
- **Per-step sentinel files** and logs under `logs/`
- **Automatic GROMACS version detection** and CPU/GPU `mdrun` flag selection
- **Multiple ligand parametrization options**: ACPYPE (GAFF/AMBER), CGenFF (CHARMM)
- **Automatic topology assembly** and `.gro` merging
- **Automatic force-field/water-model/ligand-prep compatibility rectification**
- **Centered production trajectory** output as `results/production/md_center.xtc`
- **Optional `-pbc nojump` preprocessing** before trajectory centering
- **Automatic index generation** for analysis
- **JSON config overrides** for cluster/HPC runs
- **MM-PBSA binding free energy** with per-residue decomposition option
- **PDF summary report** generation

## Pipeline Modes / Actions

| Mode | Steps | Description |
|------|-------|-------------|
| **mmpbsa** | 1-11 | Standard MD + MM-PBSA binding free energy (default) |
| **standard** | 1-10 | Standard MD up to basic analysis |
| **production** | 1-9 | Up to production MD only |

Run a specific mode:
```bash
# Full pipeline with MM-PBSA (steps 1-11)
bash run_pipeline.sh --to 11

# Standard MD only (steps 1-10)
bash run_pipeline.sh --to 10

# Production MD only (steps 1-9)
bash run_pipeline.sh --to 9
```

## Workflow Diagram

```
┌─────────────┐     ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
│  Inputs     │────▶│  Step 1-3   │────▶│  Step 4-6   │────▶│  Step 7-9   │
│ protein.pdb │     │  Protein    │     │  Box,       │     │  Min,       │
│ ligand.mol2 │     │  Ligand     │     │  Solvate,   │     │  Equil,     │
│             │     │  Assemble   │     │  Ions       │     │  Production │
└─────────────┘     └─────────────┘     └─────────────┘     └─────────────┘
                           │                   │                   │
                           ▼                   ▼                   ▼
                    ┌─────────────┐     ┌─────────────┐     ┌─────────────┐
                    │ Force Field │     │ Water Model │     │  Trajectory │
                    │ Selection   │     │ Selection   │     │  (md.xtc)   │
                    └─────────────┘     └─────────────┘     └─────────────┘
                                                                    │
                    ┌─────────────┐     ┌─────────────┐             │
                    │  Step 10    │◀────│  Step 11    │◀────────────┘
                    │  Analysis   │     │  MM-PBSA    │
                    └─────────────┘     └─────────────┘
                           │                   │
                           ▼                   ▼
                    ┌─────────────┐     ┌─────────────┐
                    │ RMSD, RMSF, │     │ ΔG_bind,    │
                    │ H-bonds,    │     │ Decomposition│
                    │ PCA, FEL    │     │ Per-residue │
                    └─────────────┘     └─────────────┘
                           │                   │
                           ▼                   ▼
                    ┌─────────────┐     ┌─────────────┐
                    │ *.xvg,      │     │ mmpbsa_     │
                    │ *.png       │     │ results.csv,│
                    └─────────────┘     │ summary.csv,│
                                        │ decomposition│
                                        │ _plot.png   │
                                        └─────────────┘
```

## Requirements

### System Dependencies
- **GROMACS** ≥ 2021 (2023+ recommended)
- **Python** ≥ 3.10
- **AmberTools** (for ACPYPE: `antechamber`, `parmchk2`, `tleap`)
- **APBS** ≥ 3.0 (for MM-PBSA polar solvation)
- **PDB2PQR** ≥ 3.0 (for MM-PBSA charge assignment)
- **OpenMPI** (optional, for parallel runs)

### Python Packages
See [`requirements.txt`](requirements.txt) for pip-installable packages and [`environment.yml`](environment.yml) for conda environment.

Key packages:
- **MD/Analysis**: `MDAnalysis` ≥ 2.7, `acpype` ≥ 2023.10
- **Scientific**: `numpy` ≥ 1.26, `scipy` ≥ 1.11, `scikit-learn` ≥ 1.4, `networkx` ≥ 3.0
- **Web**: `Flask` ≥ 3.0, `Werkzeug` ≥ 3.0
- **Visualization**: `matplotlib` ≥ 3.8
- **Web scraping**: `beautifulsoup4` ≥ 4.12, `mechanize` ≥ 0.4.0
- **Reports**: `reportlab` ≥ 4.0

### Installation (Conda Recommended)

**Option 1: Using environment.yml (recommended)**
```bash
conda env create -f environment.yml
conda activate md_pipeline
```

**Option 2: Manual conda + pip**
```bash
# Create environment
conda create -n md_pipeline python=3.10 -y
conda activate md_pipeline

# Core MD tools
conda install -c conda-forge gromacs=2024.3 -y
conda install -c conda-forge ambertools=23 -y
conda install -c conda-forge apbs=3.5 pdb2pqr=3.5 -y

# Python packages
pip install -r requirements.txt
```

### Ubuntu/Debian (Alternative)
```bash
# System packages
sudo apt-get update && sudo apt-get install -y \
  gromacs gromacs-top \
  python3-pip python3-venv \
  ambertools \
  apbs pdb2pqr

# Python packages
pip install -r requirements.txt
```

## Configuration

### Main Config (`config.sh` / `config.example.json`)

| Variable | Default | Description |
|----------|---------|-------------|
| `FORCE_FIELD` | `amber99sb-ildn` | GROMACS force field (amber99sb-ildn, charmm36, oplsaa, gromos54a7, etc.) |
| `WATER_MODEL` | `tip3p` | Water model (tip3p, spce, tip4p, tip4pew, tip5p) |
| `LIGAND_PREP_TOOL` | `acpype` | Ligand parametrization: `acpype`, `cgenff` |
| `LIGAND_RESNAME` | `LIG` | Residue name for ligand in topology |
| `BOX_DISTANCE_NM` | `1.0` | Box margin in nm |
| `ION_CONCENTRATION_M` | `0.15` | Salt concentration (M) |
| `TEMPERATURE_K` | `300` | Simulation temperature |
| `PRESSURE_BAR` | `1.0` | Simulation pressure |
| `MD_TIME_NS` | `100` | Production MD length (ns) |
| `USE_GPU` | `auto` | GPU detection: `auto`, `yes`, `no` |

### MM-PBSA Specific Config

| Variable | Default | Description |
|----------|---------|-------------|
| `MMBSA_FRAME_STEP` | `10` | Process every Nth frame for MM-PBSA |
| `MMBSA_DECOMPOSE` | `no` | Per-residue decomposition (`yes`/`no`) |
| `MMBSA_USE_APBS` | `yes` | Use APBS for polar solvation |
| `MMBSA_IONIC_STRENGTH` | `0.15` | Ionic strength for PB (M) |
| `MMBSA_SURFACE_TENSION` | `0.0072` | Surface tension for SA (kcal/mol/Å²) |

## How to Use

### 1. Prepare Inputs
```bash
mkdir -p inputs
cp /path/to/protein.pdb inputs/protein.pdb
cp /path/to/ligand.mol2 inputs/ligand.mol2
```

### 2. Quick Start (MM-PBSA Pipeline)
```bash
cd md_pipeline
bash run_pipeline.sh --to 11
```

### 3. Custom Configuration
```bash
# Using JSON config
bash run_pipeline.sh --json config.example.json

# Override specific parameters
FORCE_FIELD=charmm36 WATER_MODEL=tip3p MD_TIME_NS=200 bash run_pipeline.sh --to 11

# Enable per-residue decomposition
MMBSA_DECOMPOSE=yes bash run_pipeline.sh --to 11
```

### 4. Restart from Specific Step
```bash
# Resume from step 7 (minimization)
bash run_pipeline.sh --from 7 --to 11

# Run only MM-PBSA (step 11) after MD is done
bash run_pipeline.sh --from 11 --to 11
```

### 5. HPC/Cluster Submission
```bash
# With MPI
MPI_LAUNCHER="mpirun -np 4" NTMPI=4 NTOMP=8 bash run_pipeline.sh --to 11

# GPU node
USE_GPU=yes GPU_ID=0 bash run_pipeline.sh --to 11
```

### 6. Run MM-PBSA Separately
```bash
# MM-PBSA only (requires completed production MD)
bash steps/step11_mmpbsa.sh
```

## Output Structure

```
md_pipeline/
├── inputs/                 # Input files (protein.pdb, ligand.mol2)
├── work/                   # Intermediate files per step
│   ├── 01_protein/        # Protein prep outputs
│   ├── 02_ligand/         # Ligand parametrization
│   ├── 03_complex/        # Assembled complex
│   └── ...
├── results/
│   ├── minimization/      # EM outputs
│   ├── equilibration/     # NVT/NPT outputs
│   └── production/        # Production MD (md.xtc, md.tpr, md.gro, md_center.xtc)
├── analysis/              # All analysis outputs
│   ├── *.xvg              # GROMACS analysis plots
│   ├── *.pdb              # Representative structures
│   ├── mmpbsa/            # MM-PBSA results
│   │   ├── mmpbsa_results.csv       # Per-frame energy components
│   │   ├── mmpbsa_summary.csv       # Mean ± Std for each component
│   │   ├── mmpbsa_summary.json      # Machine-readable summary
│   │   ├── decomposition.csv        # Per-residue (if enabled)
│   │   ├── mmpbsa_decomposition.png # Bar chart of energy components
│   │   └── mmpbsa_timeseries.png    # Energy components over time
│   └── md_report.pdf      # Summary PDF report
├── logs/                   # Per-step logs
├── state/                  # Sentinel files for restart
└── README.md              # This file
```

## Key Outputs by Step

| Step | Key Outputs |
|------|-------------|
| 1-3 | `work/03_complex/complex.gro`, `work/03_complex/topol.top` |
| 4-6 | `work/06_ions/ions.gro`, `work/06_ions/topol.top` |
| 7 | `results/minimization/em.gro`, `em.edr`, `em.log` |
| 8 | `results/equilibration/nvt.gro`, `npt.gro` |
| 9 | `results/production/md.xtc`, `md.tpr`, `md.gro`, `md_center.xtc` |
| 10 | `analysis/rmsd.xvg`, `rmsf.xvg`, `hbonds.xvg`, `gyrate.xvg`, `pca_*.xvg`, `fes_*.xvg` |
| 11 | `analysis/mmpbsa/mmpbsa_summary.csv`, `mmpbsa_decomposition.png`, `md_report.pdf` |

## MM-PBSA Output Details

### `mmpbsa_results.csv` (Per-frame)
| Column | Description |
|--------|-------------|
| `frame` | Frame index |
| `time_ps` | Simulation time (ps) |
| `vdw` | van der Waals interaction (kcal/mol) |
| `elec` | Electrostatic interaction (kcal/mol) |
| `polar_solv` | Polar solvation (PB, kcal/mol) |
| `nonpolar_solv` | Non-polar solvation (SA, kcal/mol) |
| `delta_G` | Total binding free energy (kcal/mol) |

### `mmpbsa_summary.csv` (Statistics)
| Component | Mean | Std Dev | Median | Min | Max | SEM |
|-----------|------|---------|--------|-----|-----|-----|
| vdw | ... | ... | ... | ... | ... | ... |
| elec | ... | ... | ... | ... | ... | ... |
| polar_solv | ... | ... | ... | ... | ... | ... |
| nonpolar_solv | ... | ... | ... | ... | ... | ... |
| delta_G | **ΔG_bind** | ... | ... | ... | ... | ... |

### `decomposition.csv` (Per-residue, if `MMBSA_DECOMPOSE=yes`)
| Column | Description |
|--------|-------------|
| `residue` | Residue name + number (e.g., `ALA123`) |
| `resid` | Residue index |
| `resname` | Residue name |
| `vdw_mean` | Mean vdW contribution (kcal/mol) |
| `elec_mean` | Mean electrostatic contribution (kcal/mol) |
| `total_mean` | Total mean contribution (kcal/mol) |

## Troubleshooting

### Common Issues

| Issue | Solution |
|-------|----------|
| `ACPYPE` fails | Run `setup_ubuntu.sh` or ensure AmberTools in PATH |
| `CGenFF` not found | Install CGenFF or use `LIGAND_PREP_TOOL=acpype` |
| `APBS`/`PDB2PQR` missing | `conda install -c conda-forge apbs pdb2pqr` |
| GPU not detected | Set `USE_GPU=yes GPU_ID=0` explicitly |
| Out of memory | Reduce `MD_TIME_NS` or use `NTOMP=1 NTMPI=1` |
| MM-PBSA slow | Increase `MMBSA_FRAME_STEP` (e.g., 20 or 50) |

### Restart After Failure
```bash
# Check which step failed
ls state/*.done

# Resume from failed step
bash run_pipeline.sh --from 9 --to 11
```

### Log Inspection
```bash
# View latest log
tail -f logs/step9_production_*.log

# View MM-PBSA log
cat logs/step11_mmpbsa_*.log
```

## MM-PBSA Methodology

The MM-PBSA calculation follows the standard single-trajectory protocol:

1. **Molecular Mechanics (MM) Energy**: Extracted via `gmx energy` using Receptor-Ligand energy groups from a rerun TPR
2. **Polar Solvation (PB)**: Computed using APBS (Poisson-Boltzmann) on PQR files from PDB2PQR
3. **Non-polar Solvation (SA)**: SASA-based using γ·ΔSASA + b (default γ = 0.0072 kcal/mol/Å²)
4. **Binding Free Energy**: ΔG = ΔE_MM + ΔG_PB + ΔG_SA

For details, see `helpers/mmpbsa_calculation.py`.

## Citation

If you use this pipeline, please cite:
- GROMACS: Abraham et al., *SoftwareX* (2015)
- MM-PBSA: Kumari et al., *J. Chem. Inf. Model.* (2014)

## License

MIT License - See LICENSE file for details.
=======
# AutoGMX
GROMACS Automation
>>>>>>> ecb898619f8418077a6c8a09c5728549ed7b80f2
