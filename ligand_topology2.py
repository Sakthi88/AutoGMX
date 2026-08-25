#!/usr/bin/env python3
"""
Ligand Topology Generator + Automatic Topology Merger for Gromacs
Supports ACPYPE (GAFF/GAFF2) and CGenFF workflows.

Designed for integration into Gromacs MD automation pipelines.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional, List, Tuple

# ----------------------------------------------------------------------
# Utility
# ----------------------------------------------------------------------
def run_cmd(cmd: list, cwd: Optional[Path] = None, check: bool = True) -> subprocess.CompletedProcess:
    print(f"\n>>> {' '.join(map(str, cmd))}")
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.stdout.strip():
        print(result.stdout)
    if result.stderr.strip():
        print(result.stderr, file=sys.stderr)
    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed (exit {result.returncode})")
    return result

# ----------------------------------------------------------------------
# Ligand preparation
# ----------------------------------------------------------------------
def prepare_ligand(input_file: Path, out_dir: Path, add_h: bool = True) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    mol2 = out_dir / f"{input_file.stem}_prep.mol2"
    cmd = ["obabel", str(input_file), "-O", str(mol2)]
    if add_h:
        cmd.append("-h")
    run_cmd(cmd)
    print(f"[OK] Prepared: {mol2}")
    return mol2

# ----------------------------------------------------------------------
# ACPYPE
# ----------------------------------------------------------------------
def run_acpype(
    mol2: Path,
    out_dir: Path,
    net_charge: int,
    atom_type: str = "gaff2",
    charge_method: str = "bcc",
    basename: str = "LIG",
    multiplicity: int = 1,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "acpype",
        "-i", str(mol2),
        "-n", str(net_charge),
        "-m", str(multiplicity),
        "-a", atom_type,
        "-c", charge_method,
        "-b", basename,
        "-o", "gmx",
        "-d",
    ]
    run_cmd(cmd, cwd=out_dir)

    # ACPYPE creates <basename>.acpype
    acpype_dir = out_dir / f"{basename}.acpype"
    if not acpype_dir.exists():
        acpype_dir = out_dir / f"{mol2.stem}.acpype"
    if not acpype_dir.exists():
        raise FileNotFoundError("ACPYPE output directory not found")

    print(f"[OK] ACPYPE results → {acpype_dir}")
    return acpype_dir

# ----------------------------------------------------------------------
# CGenFF helpers (same as before)
# ----------------------------------------------------------------------
def prepare_cgenff(mol2: Path, out_dir: Path, resname: str = "LIG") -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / f"{resname}.mol2"
    shutil.copy(mol2, target)
    print(f"[CGenFF] Prepared mol2 → {target}")
    print("Upload this file to https://cgenff.umaryland.edu/ and download the .str")
    return target

def convert_cgenff(
    mol2: Path,
    str_file: Path,
    ff_dir: Path,
    resname: str,
    out_dir: Path,
    script_path: Optional[Path] = None,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    if script_path is None:
        for candidate in [Path("cgenff_charmm2gmx.py"), Path.home() / "cgenff_charmm2gmx.py"]:
            if candidate.exists():
                script_path = candidate
                break
        if script_path is None:
            raise FileNotFoundError("cgenff_charmm2gmx.py not found")

    cmd = [sys.executable, str(script_path), resname, str(mol2), str(str_file), str(ff_dir)]
    run_cmd(cmd, cwd=out_dir)
    print(f"[OK] CGenFF conversion finished in {out_dir}")
    return out_dir

# ----------------------------------------------------------------------
# ★ Automatic Topology Merger
# ----------------------------------------------------------------------
def merge_ligand_topology(
    protein_top: Path,
    ligand_itp: Path,
    output_top: Path,
    ligand_name: str = "LIG",
    ligand_count: int = 1,
    posre_itp: Optional[Path] = None,
    forcefield_pattern: str = r'#include\s+".*forcefield\.itp"',
) -> Path:
    """
    Merge ligand topology into a protein topology file.

    Parameters
    ----------
    protein_top : Path
        Original protein topology (from pdb2gmx)
    ligand_itp : Path
        Ligand .itp file (e.g. LIG_GMX.itp)
    output_top : Path
        Where to write the merged Complex.top
    ligand_name : str
        Residue / molecule name used in [ molecules ]
    ligand_count : int
        Number of ligand molecules
    posre_itp : Path, optional
        Position restraint file for the ligand (posre_LIG.itp)
    forcefield_pattern : str
        Regex to locate the force-field include line

    Returns
    -------
    Path to the written complex topology
    """
    if not protein_top.exists():
        raise FileNotFoundError(f"Protein topology not found: {protein_top}")
    if not ligand_itp.exists():
        raise FileNotFoundError(f"Ligand ITP not found: {ligand_itp}")

    # Read original topology
    with open(protein_top, "r") as f:
        lines = f.readlines()

    def merge_ligand_topology(
    protein_top: Path,
    ligand_itp: Path,
    output_top: Path,
    ligand_name: str = "LIG",
    ligand_count: int = 1,
    posre_itp: Optional[Path] = None,
    forcefield_pattern: str = r'#include\s+".*forcefield\.itp"',
) -> Path:
    
    # Merge ligand topology into a protein topology file.
    
    if not protein_top.exists():
        raise FileNotFoundError(f"Protein topology not found: {protein_top}")
    if not ligand_itp.exists():
        raise FileNotFoundError(f"Ligand ITP not found: {ligand_itp}")

    with open(protein_top, "r") as f:
        lines = f.readlines()

    # --------------------------------------------------------------
    # 1. Insert #include for ligand ITP right after force-field include
    # --------------------------------------------------------------
    new_lines = []
    inserted_include = False
    ff_re = re.compile(forcefield_pattern, re.IGNORECASE)

    for line in lines:
        new_lines.append(line)
        if not inserted_include and ff_re.search(line):
            new_lines.append(f'#include "{ligand_itp.name}"\n')
            if posre_itp is not None and posre_itp.exists():
                new_lines.append(f'#include "{posre_itp.name}"\n')
            inserted_include = True

    if not inserted_include:
        print("[WARNING] Could not find forcefield.itp include. Inserting after first #include.")
        new_lines = []
        first_include_done = False
        for line in lines:
            new_lines.append(line)
            if not first_include_done and line.strip().startswith("#include"):
                new_lines.append(f'#include "{ligand_itp.name}"\n')
                if posre_itp is not None and posre_itp.exists():
                    new_lines.append(f'#include "{posre_itp.name}"\n')
                first_include_done = True

    # --------------------------------------------------------------
    # 2. Add ligand to [ molecules ] section
    # --------------------------------------------------------------
    final_lines = []
    in_molecules = False
    molecules_added = False

    for line in new_lines:
        stripped = line.strip()

        if stripped.lower().startswith("[ molecules ]"):
            in_molecules = True
            final_lines.append(line)
            continue

        if in_molecules:
            if stripped and not stripped.startswith(";") and not stripped.startswith("["):
                final_lines.append(line)
            elif stripped.startswith("[") or stripped == "":
                if not molecules_added:
                    final_lines.append(f"{ligand_name:<10} {ligand_count}\n")
                    molecules_added = True
                in_molecules = False
                final_lines.append(line)
            else:
                final_lines.append(line)
        else:
            final_lines.append(line)

    if in_molecules and not molecules_added:
        final_lines.append(f"{ligand_name:<10} {ligand_count}\n")

    # --------------------------------------------------------------
    # 3. Write output
    # --------------------------------------------------------------
    output_top.parent.mkdir(parents=True, exist_ok=True)
    with open(output_top, "w") as f:
        f.writelines(final_lines)

    print(f"[OK] Merged topology written ? {output_top}")
    print(f"     Ligand include : {ligand_itp.name}")
    if posre_itp is not None:
        print(f"     Posre include  : {posre_itp.name}")
    print(f"     Molecules entry: {ligand_name}  {ligand_count}")
    return output_top

def merge_coordinates(
    protein_gro: Path,
    ligand_gro: Path,
    output_gro: Path,
) -> Path:
    """
    Simple coordinate merge: concatenate ATOM/HETATM records.
    (Useful after pdb2gmx + ACPYPE)
    """
    def extract_atoms(gro_path: Path) -> List[str]:
        lines = gro_path.read_text().splitlines()
        # GRO format: title, natoms, then atom lines, then box
        if len(lines) < 3:
            raise ValueError(f"Invalid GRO file: {gro_path}")
        natoms = int(lines[1])
        return lines[2 : 2 + natoms]

    protein_atoms = extract_atoms(protein_gro)
    ligand_atoms = extract_atoms(ligand_gro)

    title = f"Complex: {protein_gro.stem} + {ligand_gro.stem}"
    total_atoms = len(protein_atoms) + len(ligand_atoms)

    # Keep the box from the protein (or ligand if preferred)
    box_line = protein_gro.read_text().splitlines()[-1]

    with open(output_gro, "w") as f:
        f.write(f"{title}\n")
        f.write(f"{total_atoms}\n")
        for line in protein_atoms + ligand_atoms:
            f.write(line + "\n")
        f.write(box_line + "\n")

    print(f"[OK] Merged coordinates → {output_gro}  ({total_atoms} atoms)")
    return output_gro

# ----------------------------------------------------------------------
# High-level convenience function for automation pipelines
# ----------------------------------------------------------------------
def build_complex_topology(
    protein_top: Path,
    protein_gro: Path,
    acpype_dir: Path,
    output_dir: Path,
    ligand_name: str = "LIG",
    ligand_count: int = 1,
) -> Tuple[Path, Path]:
    """
    One-call helper for typical ACPYPE workflow.

    Returns
    -------
    (complex_top, complex_gro)
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Locate ACPYPE files (common naming)
    candidates_itp = list(acpype_dir.glob(f"*{ligand_name}*GMX.itp")) + list(acpype_dir.glob("*_GMX.itp"))
    candidates_gro = list(acpype_dir.glob(f"*{ligand_name}*GMX.gro")) + list(acpype_dir.glob("*_GMX.gro"))
    candidates_posre = list(acpype_dir.glob(f"posre*{ligand_name}*.itp")) + list(acpype_dir.glob("posre_*.itp"))

    if not candidates_itp:
        raise FileNotFoundError(f"No GMX.itp found in {acpype_dir}")
    ligand_itp = candidates_itp[0]

    ligand_gro = candidates_gro[0] if candidates_gro else None
    posre_itp = candidates_posre[0] if candidates_posre else None

    # Copy ITP (and posre) next to the output topology so relative includes work
    dest_itp = output_dir / ligand_itp.name
    shutil.copy(ligand_itp, dest_itp)
    if posre_itp:
        shutil.copy(posre_itp, output_dir / posre_itp.name)

    complex_top = output_dir / "Complex.top"
    merge_ligand_topology(
        protein_top=protein_top,
        ligand_itp=dest_itp,
        output_top=complex_top,
        ligand_name=ligand_name,
        ligand_count=ligand_count,
        posre_itp=output_dir / posre_itp.name if posre_itp else None,
    )

    complex_gro = None
    if ligand_gro and protein_gro.exists():
        complex_gro = output_dir / "Complex.gro"
        merge_coordinates(protein_gro, ligand_gro, complex_gro)

    return complex_top, complex_gro

# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Ligand topology generation + automatic merging for Gromacs",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--input", type=Path, help="Ligand structure (.pdb/.mol2/.sdf)")
    parser.add_argument("-o", "--outdir", type=Path, default=Path("ligand_topo"))
    parser.add_argument("-n", "--charge", type=int, default=0)
    parser.add_argument("--method", choices=["acpype", "cgenff", "both"], default="acpype")
    parser.add_argument("--atom-type", choices=["gaff", "gaff2"], default="gaff2")
    parser.add_argument("--charge-method", choices=["bcc", "gas", "user"], default="bcc")
    parser.add_argument("--basename", default="LIG")
    parser.add_argument("--no-add-h", action="store_true")
    parser.add_argument("--multiplicity", type=int, default=1)

    # CGenFF
    parser.add_argument("--cgenff-str", type=Path)
    parser.add_argument("--charmm-ff", type=Path)
    parser.add_argument("--cgenff-script", type=Path)

    # ★ Merging options
    parser.add_argument("--protein-top", type=Path, help="Protein topology from pdb2gmx")
    parser.add_argument("--protein-gro", type=Path, help="Protein coordinates from pdb2gmx")
    parser.add_argument("--merge", action="store_true",
                        help="Automatically merge ligand into protein topology")
    parser.add_argument("--ligand-count", type=int, default=1)

    args = parser.parse_args()

    if not args.input and not (args.protein_top and args.merge):
        parser.error("Either --input (for topology generation) or --protein-top + --merge is required")

    args.outdir.mkdir(parents=True, exist_ok=True)
    acpype_dir = None

    # ---- Topology generation ----
    if args.input:
        mol2 = prepare_ligand(args.input, args.outdir / "prep", add_h=not args.no_add_h)

        if args.method in ("acpype", "both"):
            acpype_dir = run_acpype(
                mol2=mol2,
                out_dir=args.outdir / "acpype",
                net_charge=args.charge,
                atom_type=args.atom_type,
                charge_method=args.charge_method,
                basename=args.basename,
                multiplicity=args.multiplicity,
            )

        if args.method in ("cgenff", "both"):
            if args.cgenff_str and args.charmm_ff:
                convert_cgenff(
                    mol2=mol2,
                    str_file=args.cgenff_str,
                    ff_dir=args.charmm_ff,
                    resname=args.basename,
                    out_dir=args.outdir / "cgenff",
                    script_path=args.cgenff_script,
                )
            else:
                prepare_cgenff(mol2, args.outdir / "cgenff", resname=args.basename)

    # ---- Automatic merging ----
    if args.merge:
        if not args.protein_top:
            sys.exit("--merge requires --protein-top")
        if acpype_dir is None:
            # Try to locate an existing ACPYPE directory
            possible = list(args.outdir.glob("**/*.acpype"))
            if not possible:
                sys.exit("No ACPYPE directory found. Run topology generation first or provide path.")
            acpype_dir = possible[0]

        complex_top, complex_gro = build_complex_topology(
            protein_top=args.protein_top,
            protein_gro=args.protein_gro if args.protein_gro else Path("dummy.gro"),
            acpype_dir=acpype_dir,
            output_dir=args.outdir / "complex",
            ligand_name=args.basename,
            ligand_count=args.ligand_count,
        )
        print("\n=== Complex ready ===")
        print(f"Topology : {complex_top}")
        if complex_gro:
            print(f"Coords   : {complex_gro}")

    print("\nDone.")

if __name__ == "__main__":
    main()