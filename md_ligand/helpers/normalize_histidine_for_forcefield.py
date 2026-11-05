#!/usr/bin/env python3
"""Translate common histidine PDB names to CHARMM36 residue names."""
import argparse
from pathlib import Path


CHARMM_HISTIDINE_NAMES = {"HID": "HSD", "HIE": "HSE", "HIP": "HSP"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb")
    parser.add_argument("--force-field", required=True)
    args = parser.parse_args()
    if not args.force_field.lower().startswith("charmm"):
        return
    path = Path(args.pdb)
    changed = 0
    output = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 20:
            old = line[17:20].strip().upper()
            new = CHARMM_HISTIDINE_NAMES.get(old)
            if new:
                line = f"{line[:17]}{new:>3}{line[20:]}"
                changed += 1
        output.append(line)
    path.write_text("".join(output), encoding="utf-8", newline="\n")
    if changed:
        print(f"Normalized {changed} CHARMM histidine atom records in {path}.")


if __name__ == "__main__":
    main()
