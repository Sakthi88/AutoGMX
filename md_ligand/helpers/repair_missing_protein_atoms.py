#!/usr/bin/env python3
"""Add missing atoms to protein residues without constructing missing loops."""
import argparse
from pathlib import Path


def residue_key(residue):
    return (
        residue.chain.id,
        residue.id,
        getattr(residue, "insertionCode", ""),
        residue.name,
    )


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

    planned_atoms = {}
    for residue, atoms in fixer.missingAtoms.items():
        planned_atoms.setdefault(residue_key(residue), set()).update(atom.name for atom in atoms)
    for residue, atoms in fixer.missingTerminals.items():
        planned_atoms.setdefault(residue_key(residue), set()).update(atoms)

    fixer.addMissingAtoms()

    rebuilt_atoms = {}
    for chain in fixer.topology.chains():
        for residue in chain.residues():
            rebuilt_atoms[residue_key(residue)] = {atom.name for atom in residue.atoms()}

    still_missing = []
    for key, atom_names in planned_atoms.items():
        absent = sorted(atom_names - rebuilt_atoms.get(key, set()))
        if absent:
            still_missing.append(f"{key[3]} {key[1]}: {', '.join(absent)}")
    if still_missing:
        raise SystemExit(
            "PDBFixer could not rebuild the following protein atom(s): "
            + "; ".join(still_missing)
            + ". Inspect nonstandard residues or structural gaps."
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

    repaired = sum(len(atom_names) for atom_names in planned_atoms.values())
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
