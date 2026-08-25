#!/usr/bin/env python3
import os
import pathlib
import re
import sys


SOLVENT_ION_NAMES = {
    "SOL",
    "WAT",
    "HOH",
    "TIP3",
    "TIP4",
    "TIP5",
    "NA",
    "NA+",
    "CL",
    "CL-",
    "K",
    "K+",
    "CA",
    "CA2+",
    "MG",
    "MG2+",
}


def strip_existing_ligand_includes(text, ligand_include, ligand_posre, ligand_mol_name):
    ligand_names = {
        pathlib.PurePosixPath(ligand_include.replace("\\", "/")).name,
        pathlib.PurePosixPath(ligand_posre.replace("\\", "/")).name,
        "ligand.prm",
        f"{ligand_mol_name}.itp",
        f"posre_{ligand_mol_name}.itp",
    }
    rebuilt = []
    skip_posres_block = False
    for line in text.splitlines(keepends=True):
        match = re.match(r'\s*#include\s+"([^"]+)"', line)
        include_name = ""
        if match:
            include_name = pathlib.PurePosixPath(match.group(1).replace("\\", "/")).name

        if include_name in ligand_names:
            continue

        if line.strip() == "#ifdef POSRES_LIG":
            skip_posres_block = True
            pending = [line]
            continue

        if skip_posres_block:
            pending.append(line)
            if line.strip() == "#endif":
                if not any(
                    re.match(r'\s*#include\s+"([^"]+)"', item)
                    and pathlib.PurePosixPath(re.match(r'\s*#include\s+"([^"]+)"', item).group(1).replace("\\", "/")).name in ligand_names
                    for item in pending
                ):
                    rebuilt.extend(pending)
                skip_posres_block = False
            continue

        rebuilt.append(line)
    return "".join(rebuilt)


def clean_molecules_section(text, ligand_mol_name, ligand_count):
    marker = "[ molecules ]"
    if marker not in text:
        raise SystemExit("Could not locate [ molecules ] section in topol.top")
    head, tail = text.split(marker, 1)
    cleaned = []
    ligand_added = False
    for line in tail.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith(";"):
            cleaned.append(line)
            continue
        name = stripped.split()[0]
        if name == ligand_mol_name:
            if not ligand_added:
                cleaned.append(f"{ligand_mol_name:<18} {ligand_count}")
                ligand_added = True
            continue
        if name.upper() in SOLVENT_ION_NAMES:
            continue
        cleaned.append(line)
    if not ligand_added:
        cleaned.append(f"{ligand_mol_name:<18} {ligand_count}")
    return head + marker + "\n".join(cleaned).rstrip() + "\n"


def find_forcefield_dir(pipeline_dir, force_field):
    """Find forcefield directory in pipeline or GROMACS data dir."""
    # Check pipeline directory first
    local_dir = pathlib.Path(pipeline_dir) / f"{force_field}.ff"
    if local_dir.exists():
        return local_dir

    # Check GROMACS data directory
    gmx_executable = os.environ.get("GMX_BIN", "gmx")
    try:
        import subprocess
        result = subprocess.run(
            [gmx_executable, "--version"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        for line in result.stdout.splitlines():
            if "Data prefix:" in line:
                data_prefix = line.split("Data prefix:", 1)[1].strip()
                gmx_data = pathlib.Path(data_prefix, "share", "gromacs", "top", f"{force_field}.ff")
                if gmx_data.exists():
                    return gmx_data
    except Exception:
        pass
    return None


def fix_forcefield_include_path(text, pipeline_dir, force_field):
    """Replace relative forcefield includes with absolute paths."""
    forcefield_dir = find_forcefield_dir(pipeline_dir, force_field)
    if not forcefield_dir:
        return text

    abs_dir_posix = pathlib.PurePosixPath(forcefield_dir.resolve()).as_posix()

    lines = text.splitlines(keepends=True)
    rebuilt = []
    for line in lines:
        # Match both './forcefield.ff/file' and 'forcefield.ff/file'
        match = re.match(r'^(\s*#include\s+")(\./)?' + re.escape(f"{force_field}.ff") + r'/([^"]+)"', line)
        if match:
            rebuilt.append(f'{match.group(1)}{abs_dir_posix}/{match.group(3)}"\n')
        else:
            rebuilt.append(line)
    return "".join(rebuilt)


def main():
    if len(sys.argv) != 8:
        raise SystemExit(
            "Usage: assemble_topology.py topol.top ligand_include ligand_posre "
            "ligand_mol_name ligand_count pipeline_dir force_field"
        )

    top_path = pathlib.Path(sys.argv[1])
    ligand_include = sys.argv[2]
    ligand_posre = sys.argv[3]
    ligand_mol_name = sys.argv[4]
    ligand_count = sys.argv[5]
    pipeline_dir = sys.argv[6]
    force_field = sys.argv[7]

    text = top_path.read_text(encoding="utf-8")
    text = strip_existing_ligand_includes(text, ligand_include, ligand_posre, ligand_mol_name)

    text = fix_forcefield_include_path(text, pipeline_dir, force_field)

    ligand_prm = pathlib.Path(top_path).with_name("ligand.prm")
    extra_include = '#include "ligand.prm"\n' if ligand_prm.exists() else ""

    insertion = (
        '; Include ligand topology\n'
        f'{extra_include}'
        f'#include "{ligand_include}"\n'
        '\n'
        '#ifdef POSRES_LIG\n'
        f'#include "{ligand_posre}"\n'
        '#endif\n\n'
    )
    lines = text.splitlines(keepends=True)
    rebuilt = []
    inserted = False
    for line in lines:
        rebuilt.append(line)
        if line.startswith('#include "') and 'forcefield.itp' in line and not inserted:
            rebuilt.append(insertion)
            inserted = True
    if not inserted:
        raise SystemExit("Could not locate forcefield include in topol.top")
    text = clean_molecules_section("".join(rebuilt), ligand_mol_name, ligand_count)
    top_path.write_text(text, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    main()
