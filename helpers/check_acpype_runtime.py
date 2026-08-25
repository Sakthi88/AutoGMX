#!/usr/bin/env python3
"""Fail early when ACPYPE's bundled AmberTools binaries cannot start."""

from __future__ import annotations

import pathlib
import os
import re
import shutil
import subprocess
import sys


def acpype_python_candidates(acpype_command: str) -> list[str]:
    candidates = [sys.executable]
    executable = pathlib.Path(shutil.which(acpype_command) or acpype_command)
    if executable.is_file():
        try:
            first_line = executable.read_text(encoding="utf-8", errors="replace").splitlines()[0]
        except (IndexError, OSError):
            first_line = ""
        if first_line.startswith("#!"):
            interpreter = first_line[2:].strip().split()[0]
            if pathlib.Path(interpreter).is_file():
                candidates.append(interpreter)

    unique = []
    for candidate in candidates:
        if candidate not in unique:
            unique.append(candidate)
    return unique


def ambertools_bin_dir(acpype_command: str) -> pathlib.Path | None:
    probe = (
        "import pathlib, acpype; "
        "print(pathlib.Path(acpype.__file__).resolve().parent / 'amber_linux' / 'bin')"
    )
    for python_bin in acpype_python_candidates(acpype_command):
        try:
            result = subprocess.run(
                [python_bin, "-c", probe],
                capture_output=True,
                text=True,
                timeout=10,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired):
            continue
        candidate = pathlib.Path(result.stdout.strip())
        if result.returncode == 0 and candidate.is_dir():
            return candidate
    return None


def ambertools_library_dir(bin_dir: pathlib.Path) -> pathlib.Path:
    """Return ACPYPE's bundled library directory used by its own launcher."""
    return bin_dir.parent / "lib"


def linked_library_errors(executable: pathlib.Path, library_dir: pathlib.Path) -> list[str]:
    """Return unresolved ELF dependencies reported by ldd, when available."""
    ldd = shutil.which("ldd")
    if not ldd:
        return []
    try:
        runtime_env = dict(os.environ)
        runtime_env["LD_LIBRARY_PATH"] = str(library_dir)
        result = subprocess.run(
            [ldd, str(executable)],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
            env=runtime_env,
        )
    except (OSError, subprocess.TimeoutExpired):
        return []

    output = f"{result.stdout}\n{result.stderr}"
    return [
        f"{executable.name}: ldd: {line.strip()}"
        for line in output.splitlines()
        if "not found" in line
    ]


def main() -> int:
    arguments = sys.argv[1:]
    require_detected = "--require-detected" in arguments
    arguments = [argument for argument in arguments if argument != "--require-detected"]
    acpype_command = arguments[0] if arguments else "acpype"
    bin_dir = ambertools_bin_dir(acpype_command)
    if bin_dir is None:
        message = "ACPYPE runtime preflight skipped: its bundled AmberTools directory was not detected."
        if require_detected:
            print(message, file=sys.stderr)
            print("Use the Python environment that installed ACPYPE, or select an installed CGenFF executable.", file=sys.stderr)
            return 1
        print(message)
        return 0

    loader_errors = []
    library_dir = ambertools_library_dir(bin_dir)
    runtime_env = dict(os.environ)
    # ACPYPE itself replaces LD_LIBRARY_PATH with amber_linux/lib before it
    # launches tleap. Test that exact environment instead of the caller's one.
    runtime_env["LD_LIBRARY_PATH"] = str(library_dir)

    for name in ("sqm", "teLeap"):
        executable = bin_dir / name
        if not executable.exists():
            continue
        loader_errors.extend(linked_library_errors(executable, library_dir))
        try:
            result = subprocess.run(
                [str(executable), "-h"],
                capture_output=True,
                text=True,
                timeout=10,
                check=False,
                env=runtime_env,
            )
        except OSError as exc:
            loader_errors.append(f"{name}: {exc}")
            continue
        except subprocess.TimeoutExpired:
            continue
        output = f"{result.stdout}\n{result.stderr}"
        if "error while loading shared libraries" in output:
            loader_errors.append(f"{name}: {output.strip()}")

    if not loader_errors:
        print("ACPYPE runtime preflight passed.")
        return 0

    print("ACPYPE cannot start its bundled AmberTools binaries.", file=sys.stderr)
    for error in loader_errors:
        print(error, file=sys.stderr)

    missing_libraries = sorted(
        {
            match.group(1) or match.group(2)
            for error in loader_errors
            for match in re.finditer(
                r"error while loading shared libraries: ([^:]+):|([^\s]+)\s+=>\s+not found",
                error,
            )
        }
    )
    if missing_libraries:
        print(f"Missing shared libraries: {', '.join(missing_libraries)}", file=sys.stderr)

    pipeline_dir = pathlib.Path(__file__).resolve().parents[1]
    print(
        "Install the matching runtime in the same environment as ACPYPE, then restart the job:\n"
        f"  bash {pipeline_dir}/setup_ubuntu.sh\n"
        "The setup script installs HDF5 1.14 and HDF4 for Conda ACPYPE installs, "
        "or the matching Ubuntu packages when Conda is not in use.",
        file=sys.stderr,
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
