#!/usr/bin/env bash
set -Eeuo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PIPELINE_DIR}/config.sh"
source "${PIPELINE_DIR}/lib/common.sh"

STEP_NAME="step3_assemble_complex"
if ! step_should_run "${STEP_NAME}"; then
  log "Skipping ${STEP_NAME}; already completed."
  exit 0
fi

init_step_logging "${STEP_NAME}"
prepare_directories

PROTEIN_WORK_DIR="${WORK_DIR}/01_protein"
LIGAND_WORK_DIR="${WORK_DIR}/02_ligand"
COMPLEX_WORK_DIR="${WORK_DIR}/03_complex"
mkdir -p "${COMPLEX_WORK_DIR}"

require_file "${PROTEIN_WORK_DIR}/protein_processed.gro"
require_file "${PROTEIN_WORK_DIR}/topol.top"
require_file "${LIGAND_WORK_DIR}/ligand.gro"
require_file "${LIGAND_WORK_DIR}/ligand.itp"
require_file "${LIGAND_WORK_DIR}/posre_ligand.itp"
require_file "${LIGAND_WORK_DIR}/ligand_metadata.env"

source "${LIGAND_WORK_DIR}/ligand_metadata.env"

cp "${PROTEIN_WORK_DIR}/topol.top" "${COMPLEX_WORK_DIR}/topol.top"
cp "${LIGAND_WORK_DIR}/ligand.itp" "${COMPLEX_WORK_DIR}/ligand.itp"
cp "${LIGAND_WORK_DIR}/posre_ligand.itp" "${COMPLEX_WORK_DIR}/posre_ligand.itp"
if [[ -f "${LIGAND_WORK_DIR}/ligand.prm" ]]; then
  cp "${LIGAND_WORK_DIR}/ligand.prm" "${COMPLEX_WORK_DIR}/ligand.prm"
fi
cp "${PROTEIN_WORK_DIR}/posre.itp" "${COMPLEX_WORK_DIR}/posre.itp"

run_cmd "${PYTHON_BIN}" "${PIPELINE_DIR}/helpers/merge_gro.py" \
  "${PROTEIN_WORK_DIR}/protein_processed.gro" \
  "${LIGAND_WORK_DIR}/ligand.gro" \
  "${COMPLEX_WORK_DIR}/complex.gro"

run_cmd "${PYTHON_BIN}" "${PIPELINE_DIR}/helpers/assemble_topology.py" \
  "${COMPLEX_WORK_DIR}/topol.top" \
  "ligand.itp" \
  "posre_ligand.itp" \
  "${LIGAND_MOLECULE_NAME}" \
  "1" \
  "${PIPELINE_DIR}" \
  "${FORCE_FIELD}"

mark_step_done "${STEP_NAME}"
log "Finished ${STEP_NAME}"
