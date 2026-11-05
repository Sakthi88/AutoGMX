#!/usr/bin/env bash
set -Eeuo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PIPELINE_DIR}/config.sh"
source "${PIPELINE_DIR}/lib/common.sh"

STEP_NAME="step0_docking_cleanup"
if ! step_should_run "${STEP_NAME}"; then
  log "Skipping ${STEP_NAME}; already completed."
  exit 0
fi

init_step_logging "${STEP_NAME}"
prepare_directories

if [[ -z "${DOCKING_COMPLEX_PDB}" ]]; then
  log "No docking complex supplied; cleanup step has nothing to process."
  mark_step_done "${STEP_NAME}"
  exit 0
fi

require_file "${DOCKING_COMPLEX_PDB}"
[[ -n "${DOCKING_LIGAND_RESNAME}" ]] || die "DOCKING_LIGAND_RESNAME is required for docking-complex cleanup."

CLEANUP_DIR="${WORK_DIR}/00_cleanup"
CLEANED_COMPLEX="${CLEANUP_DIR}/docking_complex.cleaned.pdb"
mkdir -p "${CLEANUP_DIR}"

run_cmd "${PYTHON_BIN}" "${PIPELINE_DIR}/helpers/clean_docking_complex.py" \
  "${DOCKING_COMPLEX_PDB}" "${CLEANED_COMPLEX}" "${DOCKING_LIGAND_RESNAME}"

CLEANUP_REPORT="${CLEANED_COMPLEX}.report.json"
CLEANED_LIGAND_RESNAME="$("${PYTHON_BIN}" -c 'import json,sys; print(json.load(open(sys.argv[1], encoding="utf-8"))["ligand_resname_used"])' "${CLEANUP_REPORT}")"
log "Using ligand residue name ${CLEANED_LIGAND_RESNAME} after cleanup."

CLEANUP_ENV="${STATE_DIR}/docking_cleanup.env"
mkdir -p "${STATE_DIR}"
printf 'export DOCKING_COMPLEX_PDB=%q\nexport DOCKING_LIGAND_RESNAME=%q\n' "${CLEANED_COMPLEX}" "${CLEANED_LIGAND_RESNAME}" > "${CLEANUP_ENV}.tmp"
mv "${CLEANUP_ENV}.tmp" "${CLEANUP_ENV}"

mark_step_done "${STEP_NAME}"
log "Downstream steps will use ${CLEANED_COMPLEX}"
