#!/usr/bin/env bash
set -Eeuo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ACPYPE_BIN="${ACPYPE_BIN:-acpype}"
PYTHON_BIN="${PYTHON_BIN:-}"

if [[ -z "${PYTHON_BIN}" ]]; then
  acpype_path="$(command -v "${ACPYPE_BIN}" 2>/dev/null || true)"
  if [[ -n "${acpype_path}" && -f "${acpype_path}" ]]; then
    read -r acpype_shebang < "${acpype_path}" || true
    if [[ "${acpype_shebang}" == '#!'* ]]; then
      acpype_python="${acpype_shebang:2}"
      acpype_python="${acpype_python%% *}"
      [[ -x "${acpype_python}" ]] && PYTHON_BIN="${acpype_python}"
    fi
  fi
fi
PYTHON_BIN="${PYTHON_BIN:-python3}"

if ! command -v "${ACPYPE_BIN}" >/dev/null 2>&1; then
  echo "ACPYPE command '${ACPYPE_BIN}' was not found. Install ACPYPE in the Python environment first." >&2
  exit 1
fi

ACPYPE_AMBER_LIB="$("${PYTHON_BIN}" -c '
import pathlib
import acpype
print(pathlib.Path(acpype.__file__).resolve().parent / "amber_linux" / "lib")
')"

repair_acpype_sonames() {
  [[ -d "${ACPYPE_AMBER_LIB}" ]] || return 0
  local soname target
  for soname in libhdf5.so.310 libhdf5_hl.so.310 libzip.so.5; do
    target="$(find "${ACPYPE_AMBER_LIB}" -maxdepth 1 -type f -name "${soname}.*" -printf '%f\n' | sort -V | tail -n 1)"
    [[ -n "${target}" ]] || continue
    if [[ ! -e "${ACPYPE_AMBER_LIB}/${soname}" ]]; then
      ln -s "${target}" "${ACPYPE_AMBER_LIB}/${soname}"
      echo "Linked ACPYPE runtime library: ${soname} -> ${target}"
    fi
  done
}

# ACPYPE's Python wrapper resets LD_LIBRARY_PATH to amber_linux/lib. Its
# wheel includes versioned libraries but can omit their required SONAME links.
repair_acpype_sonames

PYTHON_PREFIX="$("${PYTHON_BIN}" -c 'import sys; print(sys.prefix)')"
CONDA_BIN="${CONDA_BIN:-}"
for candidate in "${PYTHON_PREFIX}/bin/conda" "$(dirname "${PYTHON_BIN}")/conda"; do
  if [[ -z "${CONDA_BIN}" && -x "${candidate}" ]]; then
    CONDA_BIN="${candidate}"
  fi
done

if [[ -n "${CONDA_BIN}" ]]; then
  echo "Installing ACPYPE HDF runtimes into Conda environment: ${PYTHON_PREFIX}"
  "${CONDA_BIN}" install --yes --prefix "${PYTHON_PREFIX}" -c conda-forge 'hdf5>=1.14,<1.15' hdf4
  if [[ -d "${PYTHON_PREFIX}/lib" ]]; then
    export LD_LIBRARY_PATH="${PYTHON_PREFIX}/lib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
  fi
else
  if ! command -v apt-get >/dev/null 2>&1; then
    echo "No Conda environment was detected and apt-get is unavailable." >&2
    echo "Install HDF5 1.14 (libhdf5_hl.so.310) and HDF4 in ACPYPE's Python environment, then run this script again." >&2
    exit 1
  fi

  if [[ "${EUID}" -eq 0 ]]; then
    SUDO=()
  else
    if ! command -v sudo >/dev/null 2>&1; then
      echo "Run this script as root, or install sudo and run it again." >&2
      exit 1
    fi
    SUDO=(sudo)
  fi

  "${SUDO[@]}" apt-get update
  if ! apt-cache show libhdf5-hl-310 >/dev/null 2>&1; then
    echo "This Ubuntu release does not provide libhdf5-hl-310." >&2
    echo "Install ACPYPE in a Conda environment and rerun this script so it can install HDF5 1.14." >&2
    exit 1
  fi

  packages=(libhdf5-hl-310)
  if apt-cache show libhdf4-0 >/dev/null 2>&1; then
    packages+=(libhdf4-0)
  fi
  "${SUDO[@]}" apt-get install -y "${packages[@]}"
fi

"${PYTHON_BIN}" "${PIPELINE_DIR}/helpers/check_acpype_runtime.py" "${ACPYPE_BIN}"
echo "Ubuntu ACPYPE runtime setup completed."
