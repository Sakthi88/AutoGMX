#!/usr/bin/env bash
set -Eeuo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PIPELINE_DIR}/config.sh"
source "${PIPELINE_DIR}/lib/common.sh"

STEP_NAME="step5_solvate"
if ! step_should_run "${STEP_NAME}"; then
  log "Skipping ${STEP_NAME}; already completed."
  exit 0
fi

init_step_logging "${STEP_NAME}"
prepare_directories

BOX_WORK_DIR="${WORK_DIR}/04_box"
SOLV_WORK_DIR="${WORK_DIR}/05_solvate"
mkdir -p "${SOLV_WORK_DIR}"

require_file "${BOX_WORK_DIR}/boxed.gro"
require_file "${BOX_WORK_DIR}/topol.top"

cp "${BOX_WORK_DIR}/topol.top" "${SOLV_WORK_DIR}/topol.top"
cp "${BOX_WORK_DIR}/ligand.itp" "${SOLV_WORK_DIR}/ligand.itp"
cp "${BOX_WORK_DIR}/posre.itp" "${SOLV_WORK_DIR}/posre.itp"
cp "${BOX_WORK_DIR}/posre_ligand.itp" "${SOLV_WORK_DIR}/posre_ligand.itp"
if [[ -f "${BOX_WORK_DIR}/ligand.prm" ]]; then
  cp "${BOX_WORK_DIR}/ligand.prm" "${SOLV_WORK_DIR}/ligand.prm"
fi

SOLVENT_BOX="$(solvent_box_for_water_model)"

run_cmd gmx_cmd solvate \
  -cp "${BOX_WORK_DIR}/boxed.gro" \
  -cs "${SOLVENT_BOX}" \
  -o "${SOLV_WORK_DIR}/solvated.gro" \
  -p "${SOLV_WORK_DIR}/topol.top"

# Validate ligand parameters before the more expensive ion, minimization, and MD steps.
# This structure is intentionally charged at this point; step 6 adds ions to
# neutralize it. Permit only that expected grompp warning during the preflight.
render_mdp "${IONS_MDP}" "${SOLV_WORK_DIR}/topology_preflight.mdp"
if ! gmx_cmd grompp \
  -f "${SOLV_WORK_DIR}/topology_preflight.mdp" \
  -c "${SOLV_WORK_DIR}/solvated.gro" \
  -p "${SOLV_WORK_DIR}/topol.top" \
  -o "${SOLV_WORK_DIR}/topology_preflight.tpr" \
  -maxwarn 1; then
  die "Topology preflight failed before ion insertion. Inspect the GROMACS error above for a missing atom type, topology include, or incompatible ligand force field."
fi

mark_step_done "${STEP_NAME}"
log "Finished ${STEP_NAME}"
