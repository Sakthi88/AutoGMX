#!/usr/bin/env python3
"""Run LigParGen to generate OPLS-AA ligand topology.

LigParGen is a web service (http://ligpargen.utmb.edu/) that generates
OPLS-AA parameters for small molecules. This script supports both the
web API and a local installation if available.
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path
from urllib import request, parse

try:
    import requests
except Exception:
    requests = None


LIGPARGEN_URL = "http://ligpargen.utmb.edu/ligpargen"


def read_mol2_charge(mol2_path):
    """Extract formal charge from MOL2 file."""
    charge = 0
    with open(mol2_path, "r") as f:
        in_atom = False
        for line in f:
            if line.startswith("@<TRIPOS>ATOM"):
                in_atom = True
                continue
            if line.startswith("@<TRIPOS>") and in_atom:
                break
            if in_atom and line.strip():
                parts = line.split()
                if len(parts) >= 9:
                    try:
                        charge += float(parts[8])
                    except ValueError:
                        pass
    return round(charge)


def submit_ligpargen_job(mol2_path, name, net_charge):
    """Submit job to LigParGen web service."""
    if requests is None:
        raise RuntimeError("requests library required for LigParGen web submission")

    with open(mol2_path, "r") as f:
        mol2_content = f.read()

    data = {
        "mol2": mol2_content,
        "name": name,
        "charge": str(net_charge),
        "forcefield": "opls",
    }

    response = requests.post(LIGPARGEN_URL, data=data, timeout=60)
    response.raise_for_status()

    # Parse response - LigParGen returns a job ID or direct download
    try:
        result = response.json()
        job_id = result.get("job_id")
    except json.JSONDecodeError:
        # Might return HTML with job ID
        import re
        match = re.search(r'job[_\-]?id["\']?\s*[:=]\s*["\']?(\w+)', response.text, re.I)
        if match:
            job_id = match.group(1)
        else:
            raise RuntimeError(f"Could not parse LigParGen response: {response.text[:500]}")

    return job_id


def poll_ligpargen_result(job_id, max_wait=300):
    """Poll for LigParGen job completion."""
    if requests is None:
        raise RuntimeError("requests library required")

    start = time.time()
    while time.time() - start < max_wait:
        response = requests.get(f"{LIGPARGEN_URL}/status/{job_id}", timeout=30)
        if response.status_code == 200:
            try:
                result = response.json()
                if result.get("status") == "completed":
                    return result.get("download_url")
                elif result.get("status") == "failed":
                    raise RuntimeError(f"LigParGen job failed: {result.get('error')}")
            except json.JSONDecodeError:
                pass
        time.sleep(5)

    raise RuntimeError("LigParGen job timed out")


def download_ligpargen_files(download_url, output_dir):
    """Download and extract LigParGen output files."""
    if requests is None:
        raise RuntimeError("requests library required")

    response = requests.get(download_url, timeout=60)
    response.raise_for_status()

    # LigParGen returns a tar.gz or zip with topology files
    import tarfile
    import zipfile
    import io

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if download_url.endswith(".tar.gz") or download_url.endswith(".tgz"):
        with tarfile.open(fileobj=io.BytesIO(response.content), mode="r:gz") as tar:
            tar.extractall(output_dir)
    elif download_url.endswith(".zip"):
        with zipfile.ZipFile(io.BytesIO(response.content)) as zf:
            zf.extractall(output_dir)
    else:
        # Assume direct file content
        (output_dir / "ligand.itp").write_bytes(response.content)

    return output_dir


def find_output_files(output_dir, name):
    """Find generated topology files."""
    output_dir = Path(output_dir)

    # Common naming patterns from LigParGen
    itp_files = list(output_dir.glob("*.itp"))
    gro_files = list(output_dir.glob("*.gro"))
    top_files = list(output_dir.glob("*.top"))

    ligand_itp = None
    ligand_gro = None

    for f in itp_files:
        if name.lower() in f.name.lower() or "ligand" in f.name.lower():
            ligand_itp = f
            break
    if not ligand_itp and itp_files:
        ligand_itp = itp_files[0]

    for f in gro_files:
        if name.lower() in f.name.lower() or "ligand" in f.name.lower():
            ligand_gro = f
            break
    if not ligand_gro and gro_files:
        ligand_gro = gro_files[0]

    return ligand_itp, ligand_gro


def run_local_ligpargen(mol2_path, name, output_dir, net_charge):
    """Run local LigParGen installation if available."""
    ligpargen_bin = os.environ.get("LIGPARGEN_BIN", "ligpargen")
    import subprocess

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        ligpargen_bin,
        "-i", str(mol2_path),
        "-n", name,
        "-c", str(net_charge),
        "-o", str(output_dir),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=output_dir)
    if result.returncode != 0:
        raise RuntimeError(f"Local LigParGen failed: {result.stderr}")

    return find_output_files(output_dir, name)


def write_ligand_top(top_path, itp_path, molecule_name):
    """Write minimal topology file including the ligand ITP."""
    lines = [
        '#include "oplsaa.ff/forcefield.itp"',
        f'#include "{itp_path.name}"',
        "",
        "[ system ]",
        "Ligand",
        "",
        "[ molecules ]",
        f"{molecule_name} 1",
        "",
    ]
    top_path.write_text("\n".join(lines), encoding="utf-8", newline="\n")


def extract_moleculetype(itp_path):
    """Extract molecule type name from ITP file."""
    with open(itp_path, "r") as f:
        in_atoms = False
        for line in f:
            line = line.strip()
            if line.startswith("[ moleculetype ]"):
                in_atoms = True
                continue
            if in_atoms and line and not line.startswith(";"):
                return line.split()[0]
    raise RuntimeError(f"Could not determine [ moleculetype ] from {itp_path}")


def main():
    parser = argparse.ArgumentParser(description="Run LigParGen for OPLS-AA ligand topology")
    parser.add_argument("--mol2", required=True, help="Input MOL2 file")
    parser.add_argument("--name", required=True, help="Ligand name")
    parser.add_argument("--charge", type=int, default=0, help="Net charge")
    parser.add_argument("--output-itp", required=True, help="Output ITP file")
    parser.add_argument("--output-gro", required=True, help="Output GRO file")
    parser.add_argument("--output-top", required=True, help="Output TOP file")
    parser.add_argument("--work-dir", required=True, help="Working directory")
    parser.add_argument("--use-web", action="store_true", help="Use web service (default: try local first)")

    args = parser.parse_args()

    mol2_path = Path(args.mol2)
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    # Detect charge if not provided
    net_charge = args.charge
    if net_charge == 0:
        net_charge = read_mol2_charge(mol2_path)
        print(f"Detected charge from MOL2: {net_charge}")

    # Try local first, then web
    ligand_itp = None
    ligand_gro = None

    if not args.use_web:
        try:
            ligand_itp, ligand_gro = run_local_ligpargen(mol2_path, args.name, work_dir, net_charge)
        except Exception as e:
            print(f"Local LigParGen failed: {e}, trying web service...")
            args.use_web = True

    if args.use_web:
        if requests is None:
            raise RuntimeError("requests library required for web LigParGen. Install with: pip install requests")
        job_id = submit_ligpargen_job(mol2_path, args.name, net_charge)
        print(f"Submitted LigParGen job: {job_id}")
        download_url = poll_ligpargen_result(job_id)
        print(f"Job completed, downloading from: {download_url}")
        download_ligpargen_files(download_url, work_dir)
        ligand_itp, ligand_gro = find_output_files(work_dir, args.name)

    if not ligand_itp or not ligand_itp.exists():
        raise RuntimeError("LigParGen did not produce expected ITP file")

    if not ligand_gro or not ligand_gro.exists():
        # Generate GRO from MOL2 using Open Babel if needed
        import subprocess
        try:
            subprocess.run(["obabel", "-imol2", str(mol2_path), "-ogro", "-O", str(args.output_gro)], check=True)
            ligand_gro = Path(args.output_gro)
        except Exception:
            raise RuntimeError("LigParGen did not produce GRO file and Open Babel conversion failed")

    # Copy to final output locations
    import shutil
    shutil.copy2(ligand_itp, args.output_itp)
    shutil.copy2(ligand_gro, args.output_gro)

    # Extract molecule name and write topology
    molecule_name = extract_moleculetype(args.output_itp)
    write_ligand_top(args.output_top, args.output_itp, molecule_name)

    print(f"Generated: {args.output_itp}, {args.output_gro}, {args.output_top}")


if __name__ == "__main__":
    main()