#!/usr/bin/env python3
"""Validate and normalize files produced by cgenff_charmm2gmx."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


HEADER = re.compile(r"^\s*\[\s*([^\]]+)\s*\]", re.IGNORECASE)


def require_nonempty(path: Path, label: str) -> list[str]:
    if not path.is_file():
        raise ValueError(f"Missing CGenFF conversion output ({label}): {path}")
    text = path.read_text(encoding="utf-8", errors="replace")
    if not text.strip():
        raise ValueError(f"CGenFF conversion output is empty ({label}): {path}")
    return text.splitlines()


def normalize_prm(path: Path) -> None:
    lines = require_nonempty(path, "parameters")
    sections: list[tuple[str | None, list[str]]] = []
    current_name: str | None = None
    current_lines: list[str] = []
    for line in lines:
        match = HEADER.match(line)
        if match:
            if current_lines:
                sections.append((current_name, current_lines))
            current_name = match.group(1).strip().lower()
            current_lines = [line]
        else:
            current_lines.append(line)
    if current_lines:
        sections.append((current_name, current_lines))

    section_names = {name for name, _ in sections if name}
    if not section_names.intersection({"bondtypes", "angletypes", "dihedraltypes", "impropertypes"}):
        raise ValueError(f"CGenFF parameter file has no bonded parameter sections: {path}")

    order = ["atomtypes", "bondtypes", "angletypes", "dihedraltypes", "impropertypes", "nonbond_params"]
    ordered: list[str] = []
    for name in (None, *order):
        for section_name, section_lines in sections:
            if section_name == name:
                ordered.extend(section_lines)
    for section_name, section_lines in sections:
        if section_name is not None and section_name not in order:
            ordered.extend(section_lines)
    path.write_text("\n".join(ordered).rstrip() + "\n", encoding="utf-8", newline="\n")


def validate_itp(path: Path) -> None:
    lines = require_nonempty(path, "topology")
    sections = {match.group(1).strip().lower() for line in lines if (match := HEADER.match(line))}
    required = {"moleculetype", "atoms"}
    missing = required - sections
    if missing:
        raise ValueError(f"CGenFF topology is missing section(s) {', '.join(sorted(missing))}: {path}")

    atoms = 0
    in_atoms = False
    for line in lines:
        match = HEADER.match(line)
        if match:
            in_atoms = match.group(1).strip().lower() == "atoms"
            continue
        if in_atoms and line.strip() and not line.lstrip().startswith(";"):
            if line.split()[0].isdigit():
                atoms += 1
    if atoms == 0:
        raise ValueError(f"CGenFF topology has no atom records: {path}")


def validate_coordinates(path: Path) -> None:
    lines = require_nonempty(path, "initial coordinates")
    if not any(line.startswith(("ATOM  ", "HETATM")) for line in lines):
        raise ValueError(f"CGenFF initial-coordinate PDB has no atom records: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("itp")
    parser.add_argument("prm")
    parser.add_argument("initial_pdb")
    args = parser.parse_args()

    validate_itp(Path(args.itp))
    normalize_prm(Path(args.prm))
    validate_coordinates(Path(args.initial_pdb))
    print("Validated CGenFF conversion outputs.")


if __name__ == "__main__":
    main()