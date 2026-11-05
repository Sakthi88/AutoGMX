#!/usr/bin/env bash
set -Eeuo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PIPELINE_DIR}/config.sh"
source "${PIPELINE_DIR}/lib/common.sh"

STEP_NAME="step7_minimize"
if ! step_should_run "${STEP_NAME}"; then
  log "Skipping ${STEP_NAME}; already completed."
  exit 0
fi

init_step_logging "${STEP_NAME}"
prepare_directories

IONS_WORK_DIR="${WORK_DIR}/06_ions"
EM_WORK_DIR="${RESULTS_DIR}/minimization"
mkdir -p "${EM_WORK_DIR}"

require_file "${IONS_WORK_DIR}/solv_ions.gro"
require_file "${IONS_WORK_DIR}/topol.top"

cp "${IONS_WORK_DIR}/topol.top" "${EM_WORK_DIR}/topol.top"
# Copy every auxiliary *.itp include topol.top may reference: ligand.itp,
# posre_ligand.itp, and the protein's own includes (posre.itp for a
# single-chain protein, or posre_Protein_chain_A.itp/
# topol_Protein_chain_A.itp/... for multi-chain proteins). topol.top's
# #include list already matches whichever files earlier steps produced.
shopt -s nullglob
itp_files=("${IONS_WORK_DIR}"/*.itp)
shopt -u nullglob
if [[ ${#itp_files[@]} -eq 0 ]]; then
  die "No *.itp include files found in ${IONS_WORK_DIR}"
fi
for itp_file in "${itp_files[@]}"; do
  cp "${itp_file}" "${EM_WORK_DIR}/$(basename "${itp_file}")"
done
if [[ -f "${IONS_WORK_DIR}/ligand.prm" ]]; then
  cp "${IONS_WORK_DIR}/ligand.prm" "${EM_WORK_DIR}/ligand.prm"
fi

render_mdp "${EM_MDP}" "${EM_WORK_DIR}/em.mdp"

run_cmd gmx_cmd grompp \
  -f "${EM_WORK_DIR}/em.mdp" \
  -c "${IONS_WORK_DIR}/solv_ions.gro" \
  -p "${EM_WORK_DIR}/topol.top" \
  -o "${EM_WORK_DIR}/em.tpr" \
  -maxwarn "${GMX_MAXWARN}"

run_cmd run_mdrun_with_fallback "${EM_WORK_DIR}/em"
make_protein_ligand_index "${EM_WORK_DIR}/em.gro" "${EM_WORK_DIR}/indexlig.ndx"

mark_step_done "${STEP_NAME}"
log "Finished ${STEP_NAME}"
