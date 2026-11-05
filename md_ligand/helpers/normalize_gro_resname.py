#!/usr/bin/env python3
import sys


def main():
    if len(sys.argv) != 4:
        raise SystemExit("Usage: normalize_gro_resname.py input.gro output.gro RES")

    input_path, output_path, residue_name = sys.argv[1], sys.argv[2], sys.argv[3][:5]
    with open(input_path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    natoms = int(lines[1].strip())
    atoms = lines[2:2 + natoms]
    box = lines[2 + natoms]

    updated = []
    for line in atoms:
        resnum = line[:5]
        resname = f"{residue_name:>5s}"
        updated.append(f"{resnum}{resname}{line[10:]}")

    with open(output_path, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(lines[0])
        handle.write(lines[1])
        for line in updated:
            handle.write(line if line.endswith("\n") else f"{line}\n")
        handle.write(box if box.endswith("\n") else f"{box}\n")


if __name__ == "__main__":
    main()
