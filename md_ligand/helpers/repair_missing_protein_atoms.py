#!/usr/bin/env python3
"""Add missing atoms to protein residues without constructing missing loops."""
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Rebuild missing atoms in standard residues using PDBFixer."
    )
    parser.add_argument("input_pdb", type=Path)
    parser.add_argument("output_pdb", type=Path)
    args = parser.parse_args()

    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError as exc:
        raise SystemExit(
            "Protein atom repair requires PDBFixer and OpenMM. Install them in this "
            "environment with: conda install -c conda-forge pdbfixer"
        ) from exc

    if not args.input_pdb.is_file():
        raise SystemExit(f"Protein PDB not found: {args.input_pdb}")

    fixer = PDBFixer(filename=str(args.input_pdb))
    fixer.findMissingResidues()
    # Do not build unresolved sequence segments; only repair atoms in residues
    # already present in the uploaded structure.
    skipped_residue_segments = len(fixer.missingResidues)
    fixer.missingResidues = {}
    fixer.findMissingAtoms()

    missing_atoms = sum(len(atoms) for atoms in fixer.missingAtoms.values())
    missing_terminal_atoms = sum(len(atoms) for atoms in fixer.missingTerminals.values())
    fixer.addMissingAtoms()

    remaining = sum(len(atoms) for atoms in fixer.missingAtoms.values())
    remaining += sum(len(atoms) for atoms in fixer.missingTerminals.values())
    if remaining:
        raise SystemExit(
            f"PDBFixer could not rebuild {remaining} missing atom(s) in {args.input_pdb}. "
            "Inspect nonstandard residues or gaps in the structure."
        )

    args.output_pdb.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = args.output_pdb.with_name(args.output_pdb.name + ".tmp")
    try:
        with temporary_path.open("w", encoding="utf-8", newline="\n") as handle:
            PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
        temporary_path.replace(args.output_pdb)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise

    repaired = missing_atoms + missing_terminal_atoms
    if repaired:
        print(f"Added {repaired} missing atom(s) to {args.output_pdb}.")
    else:
        print(f"No missing atoms found in {args.input_pdb}.")
    if skipped_residue_segments:
        print(
            f"Left {skipped_residue_segments} unresolved residue segment(s) unbuilt; "
            "only missing atoms within existing residues were repaired."
        )


if __name__ == "__main__":
    main()
