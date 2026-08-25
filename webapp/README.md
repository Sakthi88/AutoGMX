# Local MD Web Launcher

This web app is a Linux-native browser interface for the Bash-based GROMACS pipeline.

## Start

```bash
cd md_pipeline/webapp
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python app.py
```

Then open [http://127.0.0.1:5000](http://127.0.0.1:5000).

## What It Supports

- Upload `protein.pdb`
- Upload `ligand.mol2` for ACPYPE preparation
- Or upload prepared ligand topology files:
  - `ligand.itp`
  - `ligand.gro`
  - optional `posre_ligand.itp`
- Launch the pipeline as a background local job
- Review recent jobs and live log tails

## Runtime Requirements

The machine running the web app still needs the underlying CLI tools installed:

- `gmx` or `gmx_mpi`
- `acpype` for ACPYPE mode

On Ubuntu with Anaconda/Conda ACPYPE, repair the bundled AmberTools runtime from the pipeline directory:

```bash
PYTHON_BIN="$(command -v python)" ACPYPE_BIN="$(command -v acpype)" bash setup_ubuntu.sh
```

The script installs HDF5 1.14 and HDF4 into the same Conda environment, confirms `libhdf5_hl.so.310` is present, and verifies the ACPYPE runtime. The runner scopes those Conda libraries to ACPYPE and keeps system GROMACS separate.

The launcher uses the Bash pipeline and process groups provided by Linux.
