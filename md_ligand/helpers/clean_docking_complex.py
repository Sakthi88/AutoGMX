#!/usr/bin/env python3
"""Clean a docking PDB before the MD pipeline splits protein and ligand."""
import json
import os
import sys
import tempfile
from collections import defaultdict
from pathlib import Path


WATER_NAMES = {
    "HOH", "WAT", "DOD", "SOL", "H2O", "TIP", "TP3", "TP4", "TP5", "SPC",
}
SOLVENT_NAMES = {
    "ACN", "ACT", "BEN", "DMS", "DMSO", "EDO", "EOH", "GOL", "IPA", "MEO",
    "MOH", "MPD", "PEG", "PGE", "PG4", "P6G", "1PE", "2PE", "P4C",
}
METAL_ELEMENTS = {
    "LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "MN", "FE", "CO",
    "NI", "CU", "ZN", "CD", "HG", "AL", "CR", "MO", "W", "V", "Y", "AG",
    "AU", "PT",
}
METAL_RESNAMES = {
    "LI", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA", "MN", "FE", "CO",
    "NI", "CU", "ZN", "CD", "HG", "AL", "CR", "MO", "W", "V", "Y", "AG",
    "AU", "PT",
}


def residue_key(line):
    return (
        line[17:20].strip().upper(),
        line[21:22].strip(),
        line[22:26].strip(),
        line[26:27].strip(),
    )


def atom_name(line):
    return line[12:16].strip().upper()


def element_name(line):
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    if element:
        return element
    inferred = atom_name(line).lstrip("0123456789")
    return inferred[:1]


def is_hydrogen(line):
    return element_name(line) in {"H", "D"}


def occupancy(line):
    try:
        return float(line[54:60].strip())
    except (ValueError, IndexError):
        return 0.0


def altloc_priority(line):
    altloc = line[16:17].strip().upper() if len(line) > 16 else ""
    if not altloc:
        return 2
    if altloc == "A":
        return 1
    return 0


def read_first_model(path):
    lines = []
    saw_model = False
    first_model_done = False
    active = True
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith("MODEL "):
                if not saw_model:
                    saw_model = True
                    active = True
                else:
                    active = False
                continue
            if raw.startswith("ENDMDL"):
                if active:
                    first_model_done = True
                    active = False
                continue
            if saw_model and first_model_done:
                continue
            if active and raw.startswith(("ATOM  ", "HETATM", "TER   ")):
                lines.append(raw if raw.endswith("\n") else raw + "\n")
    return lines, saw_model


def clean_complex(source, output, ligand_resname):
    ligand_resname = ligand_resname.strip().upper()
    if not ligand_resname:
        raise ValueError("A ligand residue name is required for automatic docking-complex cleanup.")

    records, used_first_model = read_first_model(source)
    atoms = [line for line in records if line.startswith(("ATOM  ", "HETATM"))]
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records were found in {source}.")

    # Identify complete non-ligand HETATM residues that should be removed.
    heterogens = defaultdict(list)
    for line in atoms:
        if line.startswith("HETATM") and residue_key(line)[0] != ligand_resname:
            heterogens[residue_key(line)].append(line)

    remove_residues = {}
    for key, group in heterogens.items():
        name = key[0]
        reason = None
        if name in WATER_NAMES:
            reason = "water"
        elif name in SOLVENT_NAMES:
            reason = "solvent"
        elif any(element_name(line) in METAL_ELEMENTS for line in group) or name in METAL_RESNAMES:
            reason = "metal"
        if reason:
            remove_residues[key] = reason

    counts = {
        "input_atoms": len(atoms),
        "hydrogens_removed": 0,
        "water_atoms_removed": 0,
        "solvent_atoms_removed": 0,
        "metal_atoms_removed": 0,
        "duplicate_or_altloc_atoms_removed": 0,
        "output_atoms": 0,
        "used_first_model": used_first_model,
    }

    filtered = []
    for line in records:
        if not line.startswith(("ATOM  ", "HETATM")):
            filtered.append(line)
            continue
        key = residue_key(line)
        reason = remove_residues.get(key) if line.startswith("HETATM") else None
        if reason:
            counts[f"{reason}_atoms_removed"] += 1
            continue
        if is_hydrogen(line) and not (
            line.startswith("HETATM") and key[0] == ligand_resname
        ):
            counts["hydrogens_removed"] += 1
            continue
        filtered.append(line)

    # For repeated atom names and alternate locations, keep the highest
    # occupancy record, preferring blank altloc and then altloc A on ties.
    candidates = defaultdict(list)
    for index, line in enumerate(filtered):
        if line.startswith(("ATOM  ", "HETATM")):
            candidates[(residue_key(line), atom_name(line))].append((index, line))
    keep_indices = set()
    for group in candidates.values():
        winner = max(
            group,
            key=lambda item: (occupancy(item[1]), altloc_priority(item[1]), -item[0]),
        )
        keep_indices.add(winner[0])
        counts["duplicate_or_altloc_atoms_removed"] += len(group) - 1

    output_lines = []
    for index, line in enumerate(filtered):
        if line.startswith(("ATOM  ", "HETATM")) and index not in keep_indices:
            continue
        if line.startswith("TER   "):
            output_lines.append(line)
        elif line.startswith(("ATOM  ", "HETATM")):
            output_lines.append(line)
    output_lines.append("END\n")
    counts["output_atoms"] = sum(line.startswith(("ATOM  ", "HETATM")) for line in output_lines)

    if not any(line.startswith("ATOM  ") for line in output_lines):
        raise ValueError("Cleanup removed all protein ATOM records.")
    if not any(
        line.startswith(("ATOM  ", "HETATM")) and residue_key(line)[0] == ligand_resname
        for line in output_lines
    ):
        raise ValueError(
            f"The requested ligand residue '{ligand_resname}' is absent after cleanup."
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    fd, temp_name = tempfile.mkstemp(prefix=output.name + ".", suffix=".tmp", dir=output.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as handle:
            handle.writelines(output_lines)
        Path(temp_name).replace(output)
    except Exception:
        Path(temp_name).unlink(missing_ok=True)
        raise
    return counts


def main():
    if len(sys.argv) != 4:
        raise SystemExit(
            "Usage: clean_docking_complex.py input.pdb output.pdb ligand_resname"
        )
    source = Path(sys.argv[1])
    output = Path(sys.argv[2])
    ligand_resname = sys.argv[3]
    if not source.is_file():
        raise SystemExit(f"Docking complex not found: {source}")
    try:
        counts = clean_complex(source, output, ligand_resname)
    except (OSError, ValueError) as exc:
        raise SystemExit(str(exc)) from exc
    print(json.dumps(counts, sort_keys=True))
    print(f"Cleaned docking complex written to {output}")


if __name__ == "__main__":
    main()
