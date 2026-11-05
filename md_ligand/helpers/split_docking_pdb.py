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
}


def residue_name_from_line(line):
    # PDB columns 18-20 contain the three-character residue name.
    # The following blank and chain-ID columns must not be included.
    if len(line) >= 20:
        candidate = line[17:20].strip().upper()
        if candidate:
            return candidate
    parts = line.split()
    if len(parts) >= 4:
        return parts[3].upper()
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
                protein_lines.append(raw)
                protein_atoms += 1
            elif record == "TER   " and protein_lines:
                protein_lines.append(raw)

    if ligand_atoms == 0:
        candidates = {
            key: lines
            for key, lines in heterogen_groups.items()
            if key[0] not in WATER_AND_ION_NAMES
        }
        if len(candidates) == 1:
            key, ligand_lines = next(iter(candidates.items()))
            ligand_atoms = len(ligand_lines)
            sys.stderr.write(
                f"Warning: no atoms matched residue name(s) {', '.join(sorted(ligand_names))}; "
                f"using the only non-water, non-ion HETATM residue: {format_residue_key(key)}.\n"
            )
            if not strip_hetatm:
                protein_keys = set(candidates)
                protein_lines = [
                    line for line in protein_lines
                    if not (line.startswith("HETATM") and residue_key_from_line(line) in protein_keys)
                ]
                protein_atoms -= ligand_atoms
        elif len(candidates) > 1:
            available = ", ".join(format_residue_key(key) for key in sorted(candidates))
            raise SystemExit(
                f"No ligand atoms matched residue name(s) {', '.join(sorted(ligand_names))} in {source_path}. "
                f"Non-water HETATM residues are: {available}. Set the ligand residue name explicitly."
            )

    if ligand_atoms == 0:
        joined_names = ", ".join(sorted(ligand_names))
        raise SystemExit(
            f"No ligand atoms found in {source_path} for residue name(s): {joined_names}. "
            "Check the PDB residue-name columns and the configured ligand residue name."
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
