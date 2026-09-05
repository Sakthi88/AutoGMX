#!/usr/bin/env bash
set -Eeuo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PIPELINE_DIR}/config.sh"
source "${PIPELINE_DIR}/lib/common.sh"

STEP_NAME="step6_add_ions"
if ! step_should_run "${STEP_NAME}"; then
  log "Skipping ${STEP_NAME}; already completed."
  exit 0
fi

init_step_logging "${STEP_NAME}"
prepare_directories

SOLV_WORK_DIR="${WORK_DIR}/05_solvate"
IONS_WORK_DIR="${WORK_DIR}/06_ions"
mkdir -p "${IONS_WORK_DIR}"

require_file "${SOLV_WORK_DIR}/solvated.gro"
require_file "${SOLV_WORK_DIR}/topol.top"

cp "${SOLV_WORK_DIR}/topol.top" "${IONS_WORK_DIR}/topol.top"
# Copy every auxiliary *.itp include topol.top may reference: ligand.itp,
# posre_ligand.itp, and the protein's own includes (posre.itp for a
# single-chain protein, or posre_Protein_chain_A.itp/
# topol_Protein_chain_A.itp/... for multi-chain proteins). topol.top's
# #include list already matches whichever files earlier steps produced.
shopt -s nullglob
itp_files=("${SOLV_WORK_DIR}"/*.itp)
shopt -u nullglob
if [[ ${#itp_files[@]} -eq 0 ]]; then
  die "No *.itp include files found in ${SOLV_WORK_DIR}"
fi
for itp_file in "${itp_files[@]}"; do
  cp "${itp_file}" "${IONS_WORK_DIR}/$(basename "${itp_file}")"
done
if [[ -f "${SOLV_WORK_DIR}/ligand.prm" ]]; then
  cp "${SOLV_WORK_DIR}/ligand.prm" "${IONS_WORK_DIR}/ligand.prm"
fi

render_mdp "${IONS_MDP}" "${IONS_WORK_DIR}/ions.mdp"

# The .tpr must be generated from the charged system so genion can calculate
# and neutralize that charge. Permit the single expected PME net-charge warning.
run_cmd gmx_cmd grompp \
  -f "${IONS_WORK_DIR}/ions.mdp" \
  -c "${SOLV_WORK_DIR}/solvated.gro" \
  -p "${IONS_WORK_DIR}/topol.top" \
  -o "${IONS_WORK_DIR}/ions.tpr" \
  -maxwarn 1

printf "SOL\n" | gmx_cmd genion \
  -s "${IONS_WORK_DIR}/ions.tpr" \
  -o "${IONS_WORK_DIR}/solv_ions.gro" \
  -p "${IONS_WORK_DIR}/topol.top" \
  -pname NA \
  -nname CL \
  -neutral \
  -conc "${ION_CONCENTRATION_M}"

mark_step_done "${STEP_NAME}"
log "Finished ${STEP_NAME}"
