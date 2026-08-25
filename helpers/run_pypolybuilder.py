#!/usr/bin/env python3
import argparse
import shlex
import subprocess
from pathlib import Path


def existing_file(value, label):
    if not value:
        raise RuntimeError(f"{label} is required for pyPolyBuilder.")
    path = Path(value)
    if not path.exists():
        raise RuntimeError(f"{label} not found: {path}")
    return path


def existing_comma_files(value, label):
    if not value:
        raise RuntimeError(f"{label} is required for pyPolyBuilder.")
    paths = [Path(item.strip()) for item in str(value).split(",") if item.strip()]
    if not paths:
        raise RuntimeError(f"{label} is required for pyPolyBuilder.")
    for path in paths:
        if not path.exists():
            raise RuntimeError(f"{label} not found: {path}")
    return ",".join(str(path) for path in paths)


def extract_moleculetype(itp_path):
    found = False
    with open(itp_path, "r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if raw.startswith("[ moleculetype ]"):
                found = True
                continue
            if found and (not line or line.startswith(";")):
                continue
            if found:
                return line.split()[0]
    raise RuntimeError(f"Could not determine [ moleculetype ] from {itp_path}")


def write_ligand_top(top_path, force_field, itp_path, molecule_name):
    lines = []
    if force_field:
        lines.append(f'#include "{force_field}.ff/forcefield.itp"')
    lines.extend(
        [
            f'#include "{itp_path.name}"',
            "",
            "[ system ]",
            "Ligand",
            "",
            "[ molecules ]",
            f"{molecule_name} 1",
            "",
        ]
    )
    top_path.write_text("\n".join(lines), encoding="utf-8", newline="\n")


def main():
    parser = argparse.ArgumentParser(description="Run pyPolyBuilder and normalize ligand outputs.")
    parser.add_argument("--bin", default="pypolybuilder")
    parser.add_argument("--work-dir", required=True)
    parser.add_argument("--mode", choices=("dendrimer", "polymer"), default="dendrimer")
    parser.add_argument("--params", required=True)
    parser.add_argument("--forcefield-path", default="")
    parser.add_argument("--force-field", default="")
    parser.add_argument("--core", default="")
    parser.add_argument("--terminal", default="")
    parser.add_argument("--inter", default="")
    parser.add_argument("--bbs", default="")
    parser.add_argument("--connections", default="")
    parser.add_argument("--ngen", default="0")
    parser.add_argument("--nsteps", default="200")
    parser.add_argument("--ngenga", default="20")
    parser.add_argument("--npop", default="25")
    parser.add_argument("--name", default="ligand")
    parser.add_argument("--output-itp", required=True)
    parser.add_argument("--output-gro", required=True)
    parser.add_argument("--output-top", required=True)
    parser.add_argument("--extra-args", default="")
    args = parser.parse_args()

    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    output_itp = Path(args.output_itp)
    output_gro = Path(args.output_gro)
    output_top = Path(args.output_top)

    params_path = existing_file(args.params, "pyPolyBuilder parameter file")
    cmd = [
        args.bin,
        "--gromacs",
        "--params",
        str(params_path),
        "--output",
        str(output_itp),
        "--gro",
        str(output_gro),
        "--name",
        args.name,
        "--nsteps",
        str(args.nsteps),
        "--ngenga",
        str(args.ngenga),
        "--npop",
        str(args.npop),
    ]

    if args.forcefield_path:
        ff_path = existing_file(Path(args.forcefield_path) / "ffbonded.itp", "pyPolyBuilder force-field bonded file").parent
        existing_file(ff_path / "ffnonbonded.itp", "pyPolyBuilder force-field nonbonded file")
        cmd.extend(["--forcefield", str(ff_path)])

    if args.mode == "dendrimer":
        cmd.append("--dendrimer")
        cmd.extend(["--core", str(existing_file(args.core, "pyPolyBuilder core .itp"))])
        cmd.extend(["--ter", str(existing_file(args.terminal, "pyPolyBuilder terminal .itp"))])
        if args.inter:
            cmd.extend(["--inter", str(existing_file(args.inter, "pyPolyBuilder intermediate .itp"))])
        cmd.extend(["--ngen", str(args.ngen)])
    else:
        cmd.append("--polymer")
        cmd.extend(["--bbs", existing_comma_files(args.bbs, "pyPolyBuilder building-block .itp list")])
        cmd.extend(["--in", str(existing_file(args.connections, "pyPolyBuilder connections file"))])
        cmd.extend(["--ngen", str(args.ngen)])

    if args.extra_args:
        cmd.extend(shlex.split(args.extra_args))

    result = subprocess.run(cmd, cwd=work_dir, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError("pyPolyBuilder topology generation failed.")
    if not output_itp.exists():
        raise RuntimeError(f"pyPolyBuilder did not create expected topology: {output_itp}")
    if not output_gro.exists():
        raise RuntimeError(f"pyPolyBuilder did not create expected coordinates: {output_gro}")

    molecule_name = extract_moleculetype(output_itp)
    write_ligand_top(output_top, args.force_field, output_itp, molecule_name)


if __name__ == "__main__":
    main()
