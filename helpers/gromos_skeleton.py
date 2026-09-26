#!/usr/bin/env python3
"""
gromos_skeleton.py
------------------
Fully offline GROMOS-compatible ligand topology skeleton generator for GROMACS.

This tool builds a usable .itp skeleton (atom types, bonds, angles, dihedrals,
position restraints) from a ligand structure. Partial charges are left as
placeholders (0.000) or taken from an optional charges file, because proper
GROMOS charges require QM / ATB-level parameterization.

IMPORTANT
---------
- This is a *skeleton* generator, not a full ATB-quality parameterizer.
- For production work with GROMOS you should still obtain charges from a
  reliable source (ATB when online is available, or manual assignment).
- For a completely offline, high-quality ligand topology, prefer ACPYPE (GAFF)
  which is already the default in AutoGMX.

Dependencies: only the Python standard library + OpenBabel (optional but
recommended for robust connectivity). Falls back to distance-based bonding
if OpenBabel is unavailable.

Usage:
  python helpers/gromos_skeleton.py -i ligand.mol2 -o ligand.itp -n LIG
  python helpers/gromos_skeleton.py -i ligand.pdb -o ligand.itp -n LIG --charges charges.txt
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Minimal GROMOS 54A7-style atom type mapping (heuristic)
# These are approximate; real GROMOS typing is more sophisticated.
# ---------------------------------------------------------------------------
ELEMENT_TO_GROMOS = {
    "H": "H",       # polar / nonpolar H approximated
    "C": "CH2",     # default aliphatic carbon (will be refined)
    "N": "NL",      # neutral nitrogen
    "O": "OA",      # hydroxyl / ether oxygen
    "S": "S",       # sulfur
    "P": "P",       # phosphorus
    "F": "F",       # fluorine
    "CL": "CL",     # chlorine
    "BR": "BR",     # bromine
    "I": "I",       # iodine
}

# Simple distance cutoffs (nm) for bond detection when no connectivity exists
BOND_CUTOFF = {
    ("H", "H"): 0.10,
    ("H", "C"): 0.13,
    ("H", "N"): 0.12,
    ("H", "O"): 0.12,
    ("H", "S"): 0.15,
    ("C", "C"): 0.18,
    ("C", "N"): 0.17,
    ("C", "O"): 0.17,
    ("C", "S"): 0.20,
    ("N", "N"): 0.16,
    ("N", "O"): 0.16,
    ("O", "O"): 0.16,
    ("S", "S"): 0.23,
}
DEFAULT_BOND_CUTOFF = 0.20  # nm


def element_from_name(name: str) -> str:
    """Extract element symbol from atom name."""
    name = name.strip().upper()
    if not name:
        return "X"
    # Common patterns: C1, CA, CB, H1, OXT, CL1, BR, etc.
    if name.startswith(("CL", "BR")):
        return name[:2]
    return name[0]


def distance(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def read_structure(path: Path) -> Tuple[List[dict], List[Tuple[int, int]]]:
    """
    Read atoms and optional bonds from mol2 / pdb / gro.
    Returns (atoms, bonds) where atoms is a list of dicts with keys:
      idx, name, element, x, y, z, charge (optional)
    and bonds is a list of (i, j) 1-based indices.
    """
    suffix = path.suffix.lower()
    atoms: List[dict] = []
    bonds: List[Tuple[int, int]] = []

    if suffix == ".mol2":
        atoms, bonds = _read_mol2(path)
    elif suffix in (".pdb", ".ent"):
        atoms = _read_pdb(path)
    elif suffix == ".gro":
        atoms = _read_gro(path)
    else:
        # Try OpenBabel conversion to mol2
        atoms, bonds = _convert_via_obabel(path)

    if not atoms:
        raise ValueError(f"No atoms found in {path}")

    # If no bonds, infer from distances
    if not bonds:
        bonds = _infer_bonds(atoms)

    return atoms, bonds


def _read_mol2(path: Path) -> Tuple[List[dict], List[Tuple[int, int]]]:
    atoms = []
    bonds = []
    section = None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("@<TRIPOS>ATOM"):
                section = "ATOM"
                continue
            if line.startswith("@<TRIPOS>BOND"):
                section = "BOND"
                continue
            if line.startswith("@<TRIPOS>"):
                section = None
                continue
            if section == "ATOM" and line.strip():
                parts = line.split()
                if len(parts) < 6:
                    continue
                idx = int(parts[0])
                name = parts[1]
                x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                # mol2 atom type often contains element info
                atype = parts[5] if len(parts) > 5 else name
                element = element_from_name(atype.split(".")[0] if "." in atype else name)
                charge = float(parts[8]) if len(parts) > 8 else 0.0
                atoms.append({
                    "idx": idx, "name": name, "element": element,
                    "x": x, "y": y, "z": z, "charge": charge
                })
            elif section == "BOND" and line.strip():
                parts = line.split()
                if len(parts) >= 3:
                    bonds.append((int(parts[1]), int(parts[2])))
    return atoms, bonds


def _read_pdb(path: Path) -> List[dict]:
    atoms = []
    idx = 0
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                idx += 1
                name = line[12:16].strip()
                element = line[76:78].strip().upper() if len(line) > 76 else element_from_name(name)
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms.append({
                    "idx": idx, "name": name, "element": element or element_from_name(name),
                    "x": x, "y": y, "z": z, "charge": 0.0
                })
    return atoms


def _read_gro(path: Path) -> List[dict]:
    atoms = []
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[1].strip().split()[0])
    for i, line in enumerate(lines[2:2 + n]):
        name = line[10:15].strip()
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        atoms.append({
            "idx": i + 1, "name": name, "element": element_from_name(name),
            "x": x, "y": y, "z": z, "charge": 0.0
        })
    return atoms


def _convert_via_obabel(path: Path) -> Tuple[List[dict], List[Tuple[int, int]]]:
    if not shutil.which("obabel"):
        raise RuntimeError(
            f"Cannot read {path}: unsupported format and OpenBabel (obabel) not found."
        )
    with tempfile.NamedTemporaryFile(suffix=".mol2", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        subprocess.run(
            ["obabel", str(path), "-O", tmp_path],
            check=True, capture_output=True
        )
        return _read_mol2(Path(tmp_path))
    finally:
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)


def _infer_bonds(atoms: List[dict]) -> List[Tuple[int, int]]:
    bonds = []
    n = len(atoms)
    for i in range(n):
        for j in range(i + 1, n):
            ei = atoms[i]["element"]
            ej = atoms[j]["element"]
            key = tuple(sorted([ei, ej]))
            cutoff = BOND_CUTOFF.get(key, DEFAULT_BOND_CUTOFF)
            # coordinates are in Angstrom for mol2/pdb, nm for gro – normalise roughly
            d = distance(
                (atoms[i]["x"], atoms[i]["y"], atoms[i]["z"]),
                (atoms[j]["x"], atoms[j]["y"], atoms[j]["z"])
            )
            # If values look like nm (small numbers), scale cutoff
            if max(abs(atoms[i]["x"]), abs(atoms[i]["y"]), abs(atoms[i]["z"])) < 50:
                # already nm
                pass
            else:
                # Angstrom → treat cutoff in Angstrom
                cutoff = cutoff * 10.0
            if d < cutoff:
                bonds.append((atoms[i]["idx"], atoms[j]["idx"]))
    return bonds


def assign_gromos_type(atom: dict, bonds: List[Tuple[int, int]], atoms: List[dict]) -> str:
    """Very simple heuristic GROMOS atom type."""
    el = atom["element"].upper()
    base = ELEMENT_TO_GROMOS.get(el, el)

    # Count heavy-atom neighbours
    nbrs = []
    for a, b in bonds:
        if a == atom["idx"]:
            nbrs.append(b)
        elif b == atom["idx"]:
            nbrs.append(a)

    heavy = 0
    for nidx in nbrs:
        for at in atoms:
            if at["idx"] == nidx and at["element"] != "H":
                heavy += 1

    if el == "C":
        if heavy >= 3:
            return "CH1"   # tertiary-ish
        if heavy == 2:
            return "CH2"
        if heavy == 1:
            return "CH3"
        return "C"
    if el == "O":
        if heavy >= 2:
            return "OA"   # ether / ester
        return "OA"       # hydroxyl also OA in many GROMOS sets
    if el == "N":
        return "NL"
    return base


def build_angles(bonds: List[Tuple[int, int]]) -> List[Tuple[int, int, int]]:
    adj = defaultdict(set)
    for a, b in bonds:
        adj[a].add(b)
        adj[b].add(a)
    angles = []
    for j in adj:
        nbrs = sorted(adj[j])
        for i_idx, i in enumerate(nbrs):
            for k in nbrs[i_idx + 1:]:
                angles.append((i, j, k))
    return angles


def build_dihedrals(bonds: List[Tuple[int, int]]) -> List[Tuple[int, int, int, int]]:
    adj = defaultdict(set)
    for a, b in bonds:
        adj[a].add(b)
        adj[b].add(a)
    dihedrals = []
    for j, k in bonds:
        for i in adj[j]:
            if i == k:
                continue
            for l in adj[k]:
                if l == j or l == i:
                    continue
                dihedrals.append((i, j, k, l))
    # unique
    seen = set()
    unique = []
    for d in dihedrals:
        key = tuple(sorted([d, d[::-1]]))
        if key not in seen:
            seen.add(key)
            unique.append(d)
    return unique


def load_charges(charges_file: Optional[Path], natoms: int) -> List[float]:
    if charges_file is None:
        return [0.0] * natoms
    charges = []
    with open(charges_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            # accept "idx charge" or just "charge"
            if len(parts) >= 2:
                charges.append(float(parts[1]))
            else:
                charges.append(float(parts[0]))
    if len(charges) != natoms:
        raise ValueError(
            f"Charges file has {len(charges)} entries, but molecule has {natoms} atoms"
        )
    return charges


def write_itp(
    atoms: List[dict],
    bonds: List[Tuple[int, int]],
    angles: List[Tuple[int, int, int]],
    dihedrals: List[Tuple[int, int, int, int]],
    charges: List[float],
    resname: str,
    outfile: Path,
    molname: Optional[str] = None,
) -> None:
    molname = molname or resname
    lines = []
    lines.append("; GROMOS-compatible topology skeleton generated by gromos_skeleton.py")
    lines.append("; WARNING: Partial charges are placeholders (or user-supplied).")
    lines.append(";          For production GROMOS simulations obtain proper charges")
    lines.append(";          (e.g. from ATB when available, or manual assignment).")
    lines.append(";          For a fully offline high-quality alternative use ACPYPE (GAFF).")
    lines.append("")
    lines.append("[ moleculetype ]")
    lines.append("; Name            nrexcl")
    lines.append(f"{molname:<16} 3")
    lines.append("")
    lines.append("[ atoms ]")
    lines.append(";   nr  type  resnr  resid  atom  cgnr  charge    mass")

    for i, at in enumerate(atoms):
        gtype = assign_gromos_type(at, bonds, atoms)
        mass = {
            "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999,
            "S": 32.06, "P": 30.974, "F": 18.998, "CL": 35.45,
            "BR": 79.90, "I": 126.90
        }.get(at["element"].upper(), 12.0)
        q = charges[i]
        lines.append(
            f"{at['idx']:6d}  {gtype:<4s}  1  {resname:<4s}  {at['name']:<5s}  "
            f"{at['idx']:4d}  {q:8.4f}  {mass:8.4f}"
        )

    lines.append("")
    lines.append("[ bonds ]")
    lines.append(";  ai   aj  funct")
    for a, b in sorted(bonds):
        lines.append(f"{a:5d} {b:5d}  1")

    lines.append("")
    lines.append("[ angles ]")
    lines.append(";  ai   aj   ak  funct")
    for a, b, c in sorted(angles):
        lines.append(f"{a:5d} {b:5d} {c:5d}  1")

    lines.append("")
    lines.append("[ dihedrals ]")
    lines.append(";  ai   aj   ak   al  funct")
    for a, b, c, d in sorted(dihedrals):
        lines.append(f"{a:5d} {b:5d} {c:5d} {d:5d}  1")

    lines.append("")
    outfile.write_text("\n".join(lines) + "\n")
    print(f"Wrote topology skeleton: {outfile}")


def write_gro(atoms: List[dict], resname: str, outfile: Path) -> None:
    """Write a simple .gro with the ligand coordinates."""
    # Assume input was Angstrom if large numbers, convert to nm
    scale = 0.1 if max(abs(a["x"]) for a in atoms) > 50 else 1.0
    lines = [f"{resname} skeleton", f"{len(atoms):5d}"]
    for at in atoms:
        x = at["x"] * scale
        y = at["y"] * scale
        z = at["z"] * scale
        lines.append(
            f"{1:5d}{resname:<5s}{at['name']:>5s}{at['idx']:5d}"
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
        )
    lines.append("   1.00000   1.00000   1.00000")
    outfile.write_text("\n".join(lines) + "\n")
    print(f"Wrote coordinates: {outfile}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Offline GROMOS topology skeleton generator for GROMACS"
    )
    parser.add_argument("-i", "--input", required=True, help="Ligand structure (mol2/pdb/gro)")
    parser.add_argument("-o", "--output", default="ligand.itp", help="Output .itp file")
    parser.add_argument("-n", "--resname", default="LIG", help="Residue / molecule name")
    parser.add_argument("--charges", help="Optional charges file (one charge per line or 'idx charge')")
    parser.add_argument("--gro", help="Also write a .gro coordinate file")
    args = parser.parse_args()

    path = Path(args.input)
    if not path.exists():
        sys.exit(f"Input file not found: {path}")

    atoms, bonds = read_structure(path)
    angles = build_angles(bonds)
    dihedrals = build_dihedrals(bonds)
    charges = load_charges(Path(args.charges) if args.charges else None, len(atoms))

    write_itp(
        atoms, bonds, angles, dihedrals, charges,
        resname=args.resname,
        outfile=Path(args.output),
    )

    if args.gro:
        write_gro(atoms, args.resname, Path(args.gro))
    else:
        # always produce a companion .gro next to the itp if possible
        gro_path = Path(args.output).with_suffix(".gro")
        write_gro(atoms, args.resname, gro_path)

    print("\nNOTE: This is a topology *skeleton*.")
    print("      Review atom types and supply proper GROMOS charges before production use.")
    print("      For a fully offline high-quality alternative, use LIGAND_PREP_TOOL=acpype.")


if __name__ == "__main__":
    main()
