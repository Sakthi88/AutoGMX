#!/usr/bin/env python3
"""Run ATB (Automated Topology Builder) to generate ligand topology.

ATB (https://atb.uq.edu.au/) generates GROMOS-compatible topologies
and supports multiple force fields including GROMOS 54A7, 53A6,
AMBER, CHARMM, and OPLS.
"""
import argparse
import json
import os
import sys
import time
from pathlib import Path

try:
    import requests
except Exception:
    requests = None


ATB_API_URL = "https://atb.uq.edu.au/api"


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


def submit_atb_job(mol2_path, name, net_charge, force_field, email=None):
    """Submit job to ATB API."""
    if requests is None:
        raise RuntimeError("requests library required for ATB submission")

    with open(mol2_path, "r") as f:
        mol2_content = f.read()

    # Map force field names to ATB format
    ff_map = {
        "gromos54a7": "GROMOS54A7",
        "gromos53a6": "GROMOS53A6",
        "amber99sb-ildn": "AMBER99SB-ILDN",
        "amber14sb": "AMBER14SB",
        "charmm36": "CHARMM36",
        "oplsaa": "OPLSAA",
    }
    atb_ff = ff_map.get(force_field.lower(), "GROMOS54A7")

    data = {
        "molecule": {
            "name": name,
            "format": "mol2",
            "content": mol2_content,
            "charge": net_charge,
            "forcefield": atb_ff,
        }
    }
    if email:
        data["email"] = email

    headers = {"Content-Type": "application/json"}
    response = requests.post(f"{ATB_API_URL}/molecule", json=data, headers=headers, timeout=60)
    response.raise_for_status()

    result = response.json()
    molecule_id = result.get("molecule_id")
    if not molecule_id:
        raise RuntimeError(f"ATB submission failed: {result}")

    return molecule_id


def poll_atb_result(molecule_id, max_wait=600):
    """Poll for ATB job completion."""
    if requests is None:
        raise RuntimeError("requests library required")

    start = time.time()
    while time.time() - start < max_wait:
        response = requests.get(f"{ATB_API_URL}/molecule/{molecule_id}", timeout=30)
        if response.status_code == 200:
            result = response.json()
            status = result.get("status")
            if status == "completed":
                return result
            elif status == "failed":
                raise RuntimeError(f"ATB job failed: {result.get('error')}")
        elif response.status_code == 404:
            # Might still be processing
            pass
        time.sleep(10)

    raise RuntimeError("ATB job timed out")


def download_atb_topology(molecule_id, output_dir, force_field):
    """Download topology files from ATB."""
    if requests is None:
        raise RuntimeError("requests library required")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ATB provides different format downloads
    ff_suffix = force_field.lower().replace("-", "")
    download_url = f"{ATB_API_URL}/molecule/{molecule_id}/download/gromos/{ff_suffix}"

    response = requests.get(download_url, timeout=60)
    response.raise_for_status()

    # ATB returns a zip file with topology
    import zipfile
    import io

    with zipfile.ZipFile(io.BytesIO(response.content)) as zf:
        zf.extractall(output_dir)

    return output_dir


def find_output_files(output_dir, name, force_field):
    """Find generated topology files."""
    output_dir = Path(output_dir)

    # ATB naming convention
    ff_lower = force_field.lower()
    itp_pattern = f"*{name}*.itp"
    gro_pattern = f"*{name}*.gro"

    itp_files = list(output_dir.glob(itp_pattern))
    gro_files = list(output_dir.glob(gro_pattern))

    # Fallback to any itp/gro
    if not itp_files:
        itp_files = list(output_dir.glob("*.itp"))
    if not gro_files:
        gro_files = list(output_dir.glob("*.gro"))

    ligand_itp = itp_files[0] if itp_files else None
    ligand_gro = gro_files[0] if gro_files else None

    return ligand_itp, ligand_gro


def run_local_atb(mol2_path, name, output_dir, net_charge, force_field):
    """Run local ATB installation if available (atb.py script)."""
    atb_bin = os.environ.get("ATB_BIN", "atb")
    import subprocess

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Local ATB typically uses a Python script
    cmd = [
        "python3", atb_bin,
        "-i", str(mol2_path),
        "-n", name,
        "-c", str(net_charge),
        "-f", force_field,
        "-o", str(output_dir),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, cwd=output_dir)
    if result.returncode != 0:
        raise RuntimeError(f"Local ATB failed: {result.stderr}")

    return find_output_files(output_dir, name, force_field)


def write_ligand_top(top_path, itp_path, molecule_name, force_field):
    """Write minimal topology file including the ligand ITP."""
    ff_include = {
        "gromos54a7": 'gromos54a7.ff/forcefield.itp',
        "gromos53a6": 'gromos53a6.ff/forcefield.itp',
        "amber99sb-ildn": 'amber99sb-ildn.ff/forcefield.itp',
        "amber14sb": 'amber14sb.ff/forcefield.itp',
        "charmm36": 'charmm36.ff/forcefield.itp',
        "oplsaa": 'oplsaa.ff/forcefield.itp',
    }
    ff_file = ff_include.get(force_field.lower(), f'{force_field}.ff/forcefield.itp')

    lines = [
        f'#include "{ff_file}"',
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
        in_moleculetype = False
        for line in f:
            line = line.strip()
            if line.startswith("[ moleculetype ]"):
                in_moleculetype = True
                continue
            if in_moleculetype and line and not line.startswith(";"):
                return line.split()[0]
    raise RuntimeError(f"Could not determine [ moleculetype ] from {itp_path}")


def main():
    parser = argparse.ArgumentParser(description="Run ATB for ligand topology generation")
    parser.add_argument("--mol2", required=True, help="Input MOL2 file")
    parser.add_argument("--name", required=True, help="Ligand name")
    parser.add_argument("--charge", type=int, default=0, help="Net charge")
    parser.add_argument("--force-field", required=True, help="Force field (gromos54a7, amber99sb-ildn, etc.)")
    parser.add_argument("--output-itp", required=True, help="Output ITP file")
    parser.add_argument("--output-gro", required=True, help="Output GRO file")
    parser.add_argument("--output-top", required=True, help="Output TOP file")
    parser.add_argument("--work-dir", required=True, help="Working directory")
    parser.add_argument("--email", default="", help="Email for ATB notifications")
    parser.add_argument("--use-web", action="store_true", help="Use web API (default: try local first)")

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
            ligand_itp, ligand_gro = run_local_atb(mol2_path, args.name, work_dir, net_charge, args.force_field)
        except Exception as e:
            print(f"Local ATB failed: {e}, trying web API...")
            args.use_web = True

    if args.use_web:
        if requests is None:
            raise RuntimeError("requests library required for web ATB. Install with: pip install requests")
        molecule_id = submit_atb_job(mol2_path, args.name, net_charge, args.force_field, args.email or None)
        print(f"Submitted ATB job: {molecule_id}")
        result = poll_atb_result(molecule_id)
        print(f"Job completed")
        download_atb_topology(molecule_id, work_dir, args.force_field)
        ligand_itp, ligand_gro = find_output_files(work_dir, args.name, args.force_field)

    if not ligand_itp or not ligand_itp.exists():
        raise RuntimeError("ATB did not produce expected ITP file")

    if not ligand_gro or not ligand_gro.exists():
        # Generate GRO from MOL2 using Open Babel if needed
        import subprocess
        try:
            subprocess.run(["obabel", "-imol2", str(mol2_path), "-ogro", "-O", str(args.output_gro)], check=True)
            ligand_gro = Path(args.output_gro)
        except Exception:
            raise RuntimeError("ATB did not produce GRO file and Open Babel conversion failed")

    # Copy to final output locations
    import shutil
    shutil.copy2(ligand_itp, args.output_itp)
    shutil.copy2(ligand_gro, args.output_gro)

    # Extract molecule name and write topology
    molecule_name = extract_moleculetype(args.output_itp)
    write_ligand_top(args.output_top, args.output_itp, molecule_name, args.force_field)

    print(f"Generated: {args.output_itp}, {args.output_gro}, {args.output_top}")


if __name__ == "__main__":
    main()