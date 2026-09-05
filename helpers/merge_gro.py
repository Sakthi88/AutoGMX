#!/usr/bin/env python3
import sys


def read_gro(path):
    with open(path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()
    title = lines[0].rstrip("\n")
    natoms = int(lines[1].strip())
    atoms = [line.rstrip("\n") for line in lines[2:2 + natoms]]
    box = lines[2 + natoms].rstrip("\n")
    return title, atoms, box


def main():
    if len(sys.argv) != 4:
        raise SystemExit("Usage: merge_gro.py protein.gro ligand.gro output.gro")

    protein_title, protein_atoms, protein_box = read_gro(sys.argv[1])
    ligand_title, ligand_atoms, _ = read_gro(sys.argv[2])
    total_atoms = len(protein_atoms) + len(ligand_atoms)

    with open(sys.argv[3], "w", encoding="utf-8", newline="\n") as handle:
        handle.write(f"{protein_title} + {ligand_title}\n")
        handle.write(f"{total_atoms:5d}\n")
        for atom_line in protein_atoms + ligand_atoms:
            handle.write(f"{atom_line}\n")
        handle.write(f"{protein_box}\n")


if __name__ == "__main__":
    main()
