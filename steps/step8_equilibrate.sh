#!/usr/bin/env bash
set -Eeuo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PIPELINE_DIR}/config.sh"
source "${PIPELINE_DIR}/lib/common.sh"

STEP_NAME="step8_equilibrate"
if ! step_should_run "${STEP_NAME}"; then
  log "Skipping ${STEP_NAME}; already completed."
  exit 0
fi

init_step_logging "${STEP_NAME}"
prepare_directories

EM_WORK_DIR="${RESULTS_DIR}/minimization"
EQUIL_WORK_DIR="${RESULTS_DIR}/equilibration"
mkdir -p "${EQUIL_WORK_DIR}/nvt" "${EQUIL_WORK_DIR}/npt"

require_file "${EM_WORK_DIR}/em.gro"
require_file "${EM_WORK_DIR}/topol.top"

cp "${EM_WORK_DIR}/topol.top" "${EQUIL_WORK_DIR}/topol.top"
# Copy every auxiliary *.itp include topol.top may reference: ligand.itp,
# posre_ligand.itp, and the protein's own includes (posre.itp for a
# single-chain protein, or posre_Protein_chain_A.itp/
# topol_Protein_chain_A.itp/... for multi-chain proteins). topol.top's
# #include list already matches whichever files earlier steps produced.
shopt -s nullglob
itp_files=("${EM_WORK_DIR}"/*.itp)
shopt -u nullglob
if [[ ${#itp_files[@]} -eq 0 ]]; then
  die "No *.itp include files found in ${EM_WORK_DIR}"
fi
for itp_file in "${itp_files[@]}"; do
  cp "${itp_file}" "${EQUIL_WORK_DIR}/$(basename "${itp_file}")"
done
if [[ -f "${EM_WORK_DIR}/ligand.prm" ]]; then
  cp "${EM_WORK_DIR}/ligand.prm" "${EQUIL_WORK_DIR}/ligand.prm"
fi
if [[ -f "${EM_WORK_DIR}/indexlig.ndx" ]]; then
  cp "${EM_WORK_DIR}/indexlig.ndx" "${EQUIL_WORK_DIR}/indexlig.ndx"
fi

render_mdp "${NVT_MDP}" "${EQUIL_WORK_DIR}/nvt/nvt.mdp"
run_cmd gmx_cmd grompp \
  -f "${EQUIL_WORK_DIR}/nvt/nvt.mdp" \
  -c "${EM_WORK_DIR}/em.gro" \
  -r "${EM_WORK_DIR}/em.gro" \
  -p "${EQUIL_WORK_DIR}/topol.top" \
  -o "${EQUIL_WORK_DIR}/nvt/nvt.tpr" \
  -maxwarn "${GMX_MAXWARN}"

run_cmd run_mdrun_with_fallback "${EQUIL_WORK_DIR}/nvt/nvt"

render_mdp "${NPT_MDP}" "${EQUIL_WORK_DIR}/npt/npt.mdp"
run_cmd gmx_cmd grompp \
  -f "${EQUIL_WORK_DIR}/npt/npt.mdp" \
  -c "${EQUIL_WORK_DIR}/nvt/nvt.gro" \
  -r "${EQUIL_WORK_DIR}/nvt/nvt.gro" \
  -t "${EQUIL_WORK_DIR}/nvt/nvt.cpt" \
  -p "${EQUIL_WORK_DIR}/topol.top" \
  -o "${EQUIL_WORK_DIR}/npt/npt.tpr" \
  -maxwarn "${GMX_MAXWARN}"

run_cmd run_mdrun_with_fallback "${EQUIL_WORK_DIR}/npt/npt"

mark_step_done "${STEP_NAME}"
log "Finished ${STEP_NAME}"
