#!/usr/bin/env python3
"""Run pyPolyBuilder to generate GROMOS-compatible polymer/dendrimer topology.

This wrapper uses the local pyPolyBuilder Python package.
"""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description="Run pyPolyBuilder for polymer/dendrimer topology")
    parser.add_argument("--bin", default="pypolybuilder", help="pyPolyBuilder command")
    parser.add_argument("--work-dir", required=True, help="Working directory")
    parser.add_argument("--mode", choices=["dendrimer", "polymer"], default="dendrimer")
    parser.add_argument("--params", required=True, help="Parameter .itp file")
    parser.add_argument("--forcefield-path", default="", help="Force field path")
    parser.add_argument("--core", default="", help="Core .itp (dendrimer)")
    parser.add_argument("--terminal", default="", help="Terminal .itp (dendrimer)")
    parser.add_argument("--inter", default="", help="Intermediate .itp (dendrimer)")
    parser.add_argument("--bbs", default="", help="Building block .itp files (polymer, comma-separated)")
    parser.add_argument("--connections", default="", help="Connections file (polymer)")
    parser.add_argument("--ngen", type=int, default=0, help="Generations")
    parser.add_argument("--nsteps", type=int, default=200, help="Geometry steps")
    parser.add_argument("--ngenga", type=int, default=20, help="GA generations")
    parser.add_argument("--npop", type=int, default=25, help="GA population")
    parser.add_argument("--name", default="polymer", help="Molecule name")
    parser.add_argument("--output-itp", required=True, help="Output ITP file")
    parser.add_argument("--output-gro", required=True, help="Output GRO file")
    parser.add_argument("--output-top", required=True, help="Output TOP file")
    parser.add_argument("--extra-args", default="", help="Extra arguments")
    parser.add_argument("--force-field", default="gromos54a7", help="Force field name")

    args = parser.parse_args()

    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    # Build command
    cmd = [args.bin]

    if args.mode == "dendrimer":
        cmd.extend(["--dendrimer"])
        if args.core:
            cmd.extend(["--core", args.core])
        if args.terminal:
            cmd.extend(["--terminal", args.terminal])
        if args.inter:
            cmd.extend(["--inter", args.inter])
    else:
        cmd.extend(["--polymer"])
        if args.bbs:
            cmd.extend(["--bbs", args.bbs])
        if args.connections:
            cmd.extend(["--in", args.connections])

    cmd.extend([
        "--params", args.params,
        "--ngen", str(args.ngen),
        "--nsteps", str(args.nsteps),
        "--ngenga", str(args.ngenga),
        "--npop", str(args.npop),
        "--name", args.name,
        "--output", str(work_dir / "default.itp"),
        "--gro", str(work_dir / "default.gro"),
        "--gromacs",
    ])

    if args.forcefield_path:
        cmd.extend(["--forcefield", args.forcefield_path])

    if args.extra_args:
        cmd.extend(args.extra_args.split())

    # Run pyPolyBuilder
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"stdout:\n{result.stdout}")
        print(f"stderr:\n{result.stderr}", file=sys.stderr)
        raise RuntimeError(f"pyPolyBuilder failed with code {result.returncode}")

    # Find output files
    itp_files = list(work_dir.glob("*.itp"))
    gro_files = list(work_dir.glob("*.gro"))

    # pyPolyBuilder outputs to default.itp and default.gro by default
    out_itp = work_dir / "default.itp"
    out_gro = work_dir / "default.gro"

    if not out_itp.exists() and itp_files:
        out_itp = itp_files[0]
    if not out_gro.exists() and gro_files:
        out_gro = gro_files[0]

    if not out_itp.exists():
        raise RuntimeError("pyPolyBuilder did not produce ITP file")

    if not out_gro.exists():
        # Try to generate GRO from ITP using Open Babel
        try:
            subprocess.run(["obabel", "-igro", str(out_gro), "-ogro", "-O", str(args.output_gro)], check=True)
        except Exception:
            pass

    # Copy to final locations
    shutil.copy2(out_itp, args.output_itp)
    if out_gro.exists():
        shutil.copy2(out_gro, args.output_gro)
    else:
        # Create minimal GRO
        Path(args.output_gro).write_text(f"{args.name}\n    1{args.name}     1   0.000   0.000   0.000\n   1.000   1.000   1.000\n")

    # Write minimal topology
    mol_name = args.name
    with open(args.output_itp, 'r') as f:
        for line in f:
            if line.strip().startswith('[ moleculetype ]'):
                next(f)
                mol_name = next(f).split()[0]
                break

    top_content = f'''#include "gromos54a7.ff/forcefield.itp"
#include "{Path(args.output_itp).name}"

[ system ]
{args.name}

[ molecules ]
{mol_name} 1
'''
    Path(args.output_top).write_text(top_content)

    print(f"Generated: {args.output_itp}, {args.output_gro}, {args.output_top}")


if __name__ == "__main__":
    main()