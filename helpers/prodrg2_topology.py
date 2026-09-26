#!/usr/bin/env python3
"""
prodrg2_topology.py
-------------------
Standalone extraction of PRODRG2 (GROMOS) ligand topology handling
from CHAPERONg (abeebyekeen/CHAPERONg).

Original functions:
  - s0PrdrgSta()  → prepare protein-ligand complex coordinates
  - s0PrdrgStb()  → prepare/insert ligand topology into topol.top

Expected input structure:
  ./prodrg/
      DRGFIN.GRO
      DRGGMX.ITP

Usage example:
  from prodrg2_topology import prepare_prodrg2_complex, insert_prodrg2_topology

  prepare_prodrg2_complex(
      protein_gro="protein_processed.gro",
      ligname="LIG",
      prodrg_dir="prodrg"
  )

  insert_prodrg2_topology(
      topol_file="topol.top",
      ligname="LIG",
      prodrg_dir="prodrg"
  )
"""

import os
import shutil
import re
from pathlib import Path
from typing import Tuple, Optional


def _extract_molecule_name_from_itp(itp_path: Path) -> str:
    """
    Extract the molecule name (residue name) from a PRODRG2 DRGGMX.ITP file.
    Looks for the first entry after the [ moleculetype ] section.
    """
    moltyp_found = False
    name_header_found = False

    with open(itp_path, "r") as f:
        for line in f:
            line = line.strip()
            if line == "[ moleculetype ]":
                moltyp_found = True
                continue
            if moltyp_found and "Name" in line and "nrexcl" in line:
                name_header_found = True
                continue
            if name_header_found and line and not line.startswith(";"):
                # First non-comment token is the molecule name
                return line.split()[0]

    raise ValueError(f"Could not find molecule name in {itp_path}")


def prepare_prodrg2_complex(
    protein_gro: str,
    ligname: str,
    prodrg_dir: str = "prodrg",
    output_complex_gro: Optional[str] = None,
) -> str:
    """
    Stage A – Prepare protein–ligand complex coordinates (s0PrdrgSta equivalent).

    Parameters
    ----------
    protein_gro : str
        Path to the processed protein .gro file (usually *_processed.gro)
    ligname : str
        Desired 3-letter (or short) residue name for the ligand
    prodrg_dir : str
        Directory containing PRODRG2 output files (default: "prodrg")
    output_complex_gro : str, optional
        Name of the resulting complex .gro file.
        If None, it is auto-generated as {protein_base}-{ligname}.gro

    Returns
    -------
    str
        Path to the generated complex .gro file
    """
    prodrg_path = Path(prodrg_dir)
    drgfin = prodrg_path / "DRGFIN.GRO"
    drgitp = prodrg_path / "DRGGMX.ITP"

    if not drgfin.exists():
        raise FileNotFoundError(f"Missing {drgfin}")
    if not drgitp.exists():
        raise FileNotFoundError(f"Missing {drgitp}")

    # 1. Extract original molecule name from ITP and rename in GRO
    original_lig_id = _extract_molecule_name_from_itp(drgitp)

    # Copy and rename residue name in the GRO file
    ligand_gro = Path(f"{ligname}.gro")
    shutil.copy(drgfin, ligand_gro)

    # Simple string replacement of the residue name
    content = ligand_gro.read_text()
    content = content.replace(original_lig_id, ligname)
    ligand_gro.write_text(content)

    # 2. Count atoms in protein and ligand
    with open(protein_gro) as f:
        lines = f.readlines()
    atom_count_protein = int(lines[1].strip().split()[0])

    with open(ligand_gro) as f:
        lig_lines = f.readlines()
    atom_count_ligand = int(lig_lines[1].strip().split()[0])

    total_atoms = atom_count_protein + atom_count_ligand
    print(f"Protein atoms : {atom_count_protein}")
    print(f"Ligand atoms  : {atom_count_ligand}")
    print(f"Total atoms   : {total_atoms}")

    # 3. Build complex .gro
    if output_complex_gro is None:
        base = Path(protein_gro).stem
        if base.endswith("_processed"):
            base = base[:-10]
        output_complex_gro = f"{base}-{ligname}.gro"

    complex_path = Path(output_complex_gro)

    # Extract ligand atom lines (skip title + atom count + box)
    lig_atom_lines = []
    for line in lig_lines[2:]:
        # Skip empty lines or box vectors (contain multiple dots)
        tokens = line.split()
        if not tokens:
            continue
        if len(tokens) >= 2 and "." in tokens[0] and "." in tokens[1]:
            continue  # box line
        lig_atom_lines.append(line.rstrip("\n"))

    with open(complex_path, "w") as out:
        # Title
        out.write(f"{Path(protein_gro).stem}-{ligname}_complex\n")
        # Atom count
        out.write(f"{total_atoms}\n")

        # Protein atoms (skip title + count)
        for line in lines[2:]:
            tokens = line.split()
            if not tokens:
                continue
            # Detect box line
            if len(tokens) >= 2 and "." in tokens[0] and "." in tokens[1]:
                # Write ligand atoms just before the box
                for lline in lig_atom_lines:
                    out.write(lline + "\n")
                out.write(line)  # original box
                break
            out.write(line)

    print(f"Complex coordinates written to: {complex_path}")
    return str(complex_path)


def insert_prodrg2_topology(
    topol_file: str = "topol.top",
    ligname: str = "LIG",
    prodrg_dir: str = "prodrg",
    backup: bool = True,
) -> None:
    """
    Stage B – Insert ligand topology into topol.top (s0PrdrgStb equivalent).

    Parameters
    ----------
    topol_file : str
        Path to the system topology file (default: topol.top)
    ligname : str
        Residue name that will be used in the topology
    prodrg_dir : str
        Directory containing PRODRG2 output (default: "prodrg")
    backup : bool
        Whether to keep a backup of the original topol.top
    """
    prodrg_path = Path(prodrg_dir)
    drgitp = prodrg_path / "DRGGMX.ITP"

    if not drgitp.exists():
        raise FileNotFoundError(f"Missing {drgitp}")

    # 1. Rename molecule name inside the ITP and copy to working directory
    original_lig_id = _extract_molecule_name_from_itp(drgitp)
    ligand_itp = Path(f"{ligname}.itp")

    content = drgitp.read_text()
    content = content.replace(original_lig_id, ligname)
    ligand_itp.write_text(content)
    print(f"Ligand topology written to: {ligand_itp}")

    # 2. Modify topol.top
    topol = Path(topol_file)
    if not topol.exists():
        raise FileNotFoundError(f"{topol_file} not found")

    if backup:
        backup_path = topol.with_name("topolPreLigParmod.top")
        shutil.copy(topol, backup_path)
        print(f"Original topology backed up to: {backup_path}")

    lines = topol.read_text().splitlines(keepends=True)
    new_lines = []

    # State machine for insertion points
    after_posres = False
    molecules_section = False
    molecules_header_seen = False
    inserted_topology = False
    inserted_molecule = False

    for i, line in enumerate(lines):
        new_lines.append(line)

        # Insert ligand topology after the position-restraint block
        if "; Include Position restraint file" in line:
            after_posres = True
            continue

        if after_posres and not inserted_topology:
            # Look for the end of the #ifdef POSRES block
            if line.strip().startswith("#endif"):
                new_lines.append("\n; Include ligand topology\n")
                new_lines.append(f'#include "{ligand_itp.name}"\n')
                inserted_topology = True
                after_posres = False

        # Insert ligand into [ molecules ] section
        if line.strip() == "[ molecules ]":
            molecules_section = True
            continue

        if molecules_section and not inserted_molecule:
            # After the comment line "; Compound        #mols"
            if line.strip().startswith(";") or line.strip() == "":
                continue
            # First real molecule entry → insert ligand right after it
            new_lines.append(f"{ligname:<20} 1\n")
            inserted_molecule = True
            molecules_section = False

    if not inserted_topology:
        print("WARNING: Could not locate insertion point for ligand topology.")
    if not inserted_molecule:
        print("WARNING: Could not locate [ molecules ] section.")

    topol.write_text("".join(new_lines))
    print(f"Updated topology written to: {topol_file}")


# ----------------------------------------------------------------------
# Convenience wrapper that runs both stages
# ----------------------------------------------------------------------
def run_prodrg2_pipeline(
    protein_gro: str,
    ligname: str,
    prodrg_dir: str = "prodrg",
    topol_file: str = "topol.top",
) -> str:
    """
    Full PRODRG2 pipeline: prepare complex coordinates + insert topology.
    Returns the path of the generated complex .gro file.
    """
    complex_gro = prepare_prodrg2_complex(
        protein_gro=protein_gro,
        ligname=ligname,
        prodrg_dir=prodrg_dir,
    )
    insert_prodrg2_topology(
        topol_file=topol_file,
        ligname=ligname,
        prodrg_dir=prodrg_dir,
    )
    return complex_gro


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="PRODRG2 (GROMOS) ligand topology helper extracted from CHAPERONg"
    )
    parser.add_argument("-p", "--protein", required=True, help="Processed protein .gro file")
    parser.add_argument("-l", "--ligname", required=True, help="Desired ligand residue name")
    parser.add_argument("-d", "--prodrg-dir", default="prodrg", help="PRODRG2 output directory")
    parser.add_argument("-t", "--topol", default="topol.top", help="System topology file")
    parser.add_argument("--stage", choices=["a", "b", "both"], default="both",
                        help="a = coordinates only, b = topology only, both = full pipeline")

    args = parser.parse_args()

    if args.stage in ("a", "both"):
        prepare_prodrg2_complex(args.protein, args.ligname, args.prodrg_dir)

    if args.stage in ("b", "both"):
        insert_prodrg2_topology(args.topol, args.ligname, args.prodrg_dir)
