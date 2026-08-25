#!/usr/bin/env python3
import sys
from pathlib import Path


VALID_RECORDS = ("ATOM  ", "HETATM", "TER   ", "END")


def residue_name_from_line(line):
    if len(line) >= 21:
        candidate = line[17:21].strip().upper()
        if candidate:
            return candidate
    parts = line.split()
    if len(parts) >= 4:
        return parts[3].upper()
    return ""


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
    fallback_ligand_lines = []
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
                continue

            if record == "HETATM":
                fallback_ligand_lines.append(raw)
                continue

            if record == "ATOM  ":
                protein_lines.append(raw)
                protein_atoms += 1
                continue

            if record == "HETATM" and not strip_hetatm:
                protein_lines.append(raw)
                protein_atoms += 1
                continue

            if record == "TER   " and protein_lines:
                protein_lines.append(raw)

    if ligand_atoms == 0 and fallback_ligand_lines:
        ligand_lines = list(fallback_ligand_lines)
        ligand_atoms = len(ligand_lines)
        sys.stderr.write(
            f"Warning: no atoms matched residue name(s) {', '.join(sorted(ligand_names))}; "
            "falling back to all HETATM records in the docking complex.\n"
        )

    if ligand_atoms == 0:
        joined_names = ", ".join(sorted(ligand_names))
        raise SystemExit(f"No ligand atoms found in {source_path} for residue name(s): {joined_names}")
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
