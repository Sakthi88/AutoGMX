#!/usr/bin/env python3
"""Split a docked PDB complex into protein and ligand coordinate files."""
import sys
from collections import defaultdict
from pathlib import Path


VALID_RECORDS = ("ATOM  ", "HETATM", "TER   ", "END")
WATER_AND_ION_NAMES = {
    "HOH", "WAT", "SOL", "SPC", "TIP", "TIP3", "TIP4", "TIP5",
    "NA", "NA+", "CL", "CL-", "K", "K+", "LI", "LI+", "MG", "MG2",
    "CA", "CA2", "ZN", "FE", "FE2", "FE3", "MN", "CU", "CO", "NI",
    "ACN", "ACT", "BEN", "DMS", "DMSO", "EDO", "EOH", "GOL", "IPA", "MEO",
    "MOH", "MPD", "PEG", "PGE", "PG4", "1PE", "2PE", "P4C", "P6G",
    "RB", "CS", "SR", "BA", "CD", "HG", "AL", "CR", "MO", "W", "V", "Y",
    "AG", "AU", "PT",
}


def residue_name_from_line(line):
    # PDB columns 18-20 contain the three-character residue name.
    # The following blank and chain-ID columns must not be included.
    parts = line.split()
    if len(parts) >= 4 and parts[0] in {"ATOM", "HETATM"}:
        token_name = parts[3].upper()
        fixed_name = line[17:20].strip().upper() if len(line) >= 20 else ""
        # Prefer whitespace fields when a loosely formatted PDB has shifted columns.
        if token_name and token_name != fixed_name:
            return token_name
        if fixed_name:
            return fixed_name
    return ""


def residue_key_from_line(line):
    return (
        residue_name_from_line(line),
        line[21:22].strip() if len(line) > 21 else "",
        line[22:26].strip() if len(line) > 25 else "",
        line[26:27].strip() if len(line) > 26 else "",
    )


def format_residue_key(key):
    name, chain, number, insertion = key
    return f"{name} chain {chain or '-'} residue {number}{insertion}"


def main():
    if len(sys.argv) != 6:
        raise SystemExit(
            "Usage: split_docking_pdb.py input.pdb protein_out.pdb ligand_out.pdb ligand_resname strip_hetatm"
        )

    source_path = Path(sys.argv[1])
    protein_out = Path(sys.argv[2])
    ligand_out = Path(sys.argv[3])
    ligand_names = {name.strip().upper() for name in sys.argv[4].split(",") if name.strip()}
    strip_hetatm = sys.argv[5].strip().lower() == "yes"

    if not ligand_names:
        raise SystemExit("At least one docking ligand residue name is required.")

    protein_lines = []
    ligand_lines = []
    heterogen_groups = defaultdict(list)
    atom_residue_groups = defaultdict(list)
    ligand_atoms = 0
    protein_atoms = 0

    with source_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.startswith(VALID_RECORDS):
                continue

            record = raw[:6]
            resname = residue_name_from_line(raw)
            is_ligand = record in ("ATOM  ", "HETATM") and resname in ligand_names

            if is_ligand:
                ligand_lines.append(raw)
                ligand_atoms += 1
            elif record == "HETATM":
                heterogen_groups[residue_key_from_line(raw)].append(raw)
                if not strip_hetatm:
                    protein_lines.append(raw)
                    protein_atoms += 1
            elif record == "ATOM  ":
                atom_residue_groups[residue_key_from_line(raw)].append(raw)
                protein_lines.append(raw)
                protein_atoms += 1
            elif record == "TER   " and protein_lines:
                protein_lines.append(raw)

    inferred_key = None
    if ligand_atoms == 0:
        candidates = {
            key: lines
            for key, lines in heterogen_groups.items()
            if key[0] not in WATER_AND_ION_NAMES
        }
        if not candidates:
            # Some docking/export tools write ligand coordinates as ATOM records.
            # Consider only nonstandard residue names after ordinary HETATM detection.
            protein_residue_names = {
                "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
                "ASH", "GLH", "CYM", "CYX", "HID", "HIE", "HIP", "LYN", "MSE",
            }
            candidates = {
                key: lines
                for key, lines in atom_residue_groups.items()
                if key[0] not in protein_residue_names | WATER_AND_ION_NAMES
            }
        if len(candidates) == 1:
            inferred_key, ligand_lines = next(iter(candidates.items()))
            ligand_atoms = len(ligand_lines)
            source_kind = "ATOM" if inferred_key in atom_residue_groups else "HETATM"
            sys.stderr.write(
                f"Warning: no atoms matched residue name(s) {', '.join(sorted(ligand_names))}; "
                f"using the only inferred ligand residue: {format_residue_key(inferred_key)} "
                f"({source_kind} records).\n"
            )
            if source_kind == "ATOM":
                protein_lines = [
                    line for line in protein_lines
                    if not (line.startswith("ATOM  ") and residue_key_from_line(line) == inferred_key)
                ]
                protein_atoms -= ligand_atoms
        elif len(candidates) > 1:
            available = ", ".join(format_residue_key(key) for key in sorted(candidates))
            raise SystemExit(
                f"No ligand atoms matched residue name(s) {', '.join(sorted(ligand_names))} in {source_path}. "
                f"Possible ligand residues are: {available}. Choose automatic HETATM detection or enter the exact residue name."
            )

    if ligand_atoms == 0:
        observed_hetatm = sorted({key[0] for key in heterogen_groups})
        observed_atom_nonstandard = sorted({
            key[0] for key in atom_residue_groups
            if key[0] not in {
                "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
                "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
            }
        })
        joined_names = ", ".join(sorted(ligand_names))
        found = f"HETATM residues: {', '.join(observed_hetatm) or 'none'}; nonstandard ATOM residues: {', '.join(observed_atom_nonstandard) or 'none'}"
        raise SystemExit(
            f"No ligand atoms found in {source_path} for residue name(s): {joined_names}. "
            f"Parsed {found}. Enter the ligand residue name or use automatic detection."
        )
    if protein_atoms == 0:
        raise SystemExit(f"No protein atoms found in {source_path}")

    protein_out.parent.mkdir(parents=True, exist_ok=True)
    ligand_out.parent.mkdir(parents=True, exist_ok=True)

    with protein_out.open("w", encoding="utf-8", newline="\n") as handle:
        for line in protein_lines:
            handle.write(line if line.endswith("\n") else f"{line}\n")
        handle.write("END\n")

    with ligand_out.open("w", encoding="utf-8", newline="\n") as handle:
        for line in ligand_lines:
            handle.write(line if line.endswith("\n") else f"{line}\n")
        handle.write("END\n")


if __name__ == "__main__":
    main()
