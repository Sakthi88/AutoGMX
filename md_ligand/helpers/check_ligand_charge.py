#!/usr/bin/env python3
"""Reject an ACPYPE closed-shell charge that gives an odd electron count."""

from __future__ import annotations

import re
import sys
from pathlib import Path

ATOMIC_NUMBERS = {
    "H": 1, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "P": 15, "S": 16,
    "Cl": 17, "Br": 35, "I": 53, "Si": 14, "Na": 11, "Mg": 12, "Ca": 20,
    "Zn": 30, "Fe": 26,
}
TWO_LETTER_ELEMENTS = {symbol.lower() for symbol in ATOMIC_NUMBERS if len(symbol) == 2}


def atom_element(atom_name: str, atom_type: str) -> str | None:
    type_token = re.sub(r"[^A-Za-z]", "", atom_type.split(".", 1)[0])
    if type_token:
        normalized = type_token.lower()
        if normalized in TWO_LETTER_ELEMENTS:
            return normalized.title()
        return normalized[0].upper()
    name_token = re.sub(r"[^A-Za-z]", "", atom_name)
    if name_token:
        normalized = name_token.lower()
        if normalized[:2] in TWO_LETTER_ELEMENTS:
            return normalized[:2].title()
        return normalized[0].upper()
    return None


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: check_ligand_charge.py ligand.mol2 net_charge")

    try:
        net_charge = int(sys.argv[2])
    except ValueError:
        raise SystemExit("Ligand net charge must be an integer for ACPYPE.")

    lines = Path(sys.argv[1]).read_text(encoding="utf-8", errors="replace").splitlines()
    in_atoms = False
    atomic_number_sum = 0
    unrecognized = []
    for line in lines:
        if line.startswith("@<TRIPOS>ATOM"):
            in_atoms = True
            continue
        if in_atoms and line.startswith("@<TRIPOS>"):
            break
        if not in_atoms or not line.strip():
            continue
        fields = line.split()
        if len(fields) < 6:
            continue
        element = atom_element(fields[1], fields[5])
        if element not in ATOMIC_NUMBERS:
            unrecognized.append(fields[1])
            continue
        atomic_number_sum += ATOMIC_NUMBERS[element]

    if atomic_number_sum == 0 or unrecognized:
        print("Ligand charge preflight skipped because one or more MOL2 atom elements could not be identified.")
        return 0

    electrons = atomic_number_sum - net_charge
    if electrons % 2 == 0:
        print(f"Ligand charge preflight passed: charge {net_charge} gives {electrons} electrons.")
        return 0

    raise SystemExit(
        f"Selected ligand charge {net_charge} gives {electrons} electrons, but ACPYPE runs a closed-shell calculation by default. "
        "Choose the chemically validated formal charge (an adjacent integer charge gives an even electron count), "
        "or parameterize a confirmed radical separately with its correct spin multiplicity."
    )


if __name__ == "__main__":
    raise SystemExit(main())