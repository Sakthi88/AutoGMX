#!/usr/bin/env bash
set -Eeuo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PIPELINE_DIR}/config.sh"
source "${PIPELINE_DIR}/lib/common.sh"

STEP_NAME="step9_production"
if ! step_should_run "${STEP_NAME}"; then
  log "Skipping ${STEP_NAME}; already completed."
  exit 0
fi

init_step_logging "${STEP_NAME}"
prepare_directories

EQUIL_WORK_DIR="${RESULTS_DIR}/equilibration"
if [[ "${ENSEMBLE_REPLICAS}" != "1" ]]; then
  [[ "${ENSEMBLE_REPLICAS}" =~ ^[2-9][0-9]*$ ]] || die "ENSEMBLE_REPLICAS must be a positive integer"
  log "Running ${ENSEMBLE_REPLICAS} deterministic production replicas."
  bash "${PIPELINE_DIR}/steps/step9_ensemble.sh"
  exit 0
fi
PROD_WORK_DIR="${RESULTS_DIR}/production"
mkdir -p "${PROD_WORK_DIR}"

require_file "${EQUIL_WORK_DIR}/npt/npt.gro"
require_file "${EQUIL_WORK_DIR}/npt/npt.cpt"
require_file "${EQUIL_WORK_DIR}/topol.top"

cp "${EQUIL_WORK_DIR}/topol.top" "${PROD_WORK_DIR}/topol.top"
# Copy every auxiliary *.itp include topol.top may reference: ligand.itp,
# posre_ligand.itp, and the protein's own includes (posre.itp for a
# single-chain protein, or posre_Protein_chain_A.itp/
# topol_Protein_chain_A.itp/... for multi-chain proteins). topol.top's
# #include list already matches whichever files earlier steps produced.
shopt -s nullglob
itp_files=("${EQUIL_WORK_DIR}"/*.itp)
shopt -u nullglob
if [[ ${#itp_files[@]} -eq 0 ]]; then
  die "No *.itp include files found in ${EQUIL_WORK_DIR}"
fi
for itp_file in "${itp_files[@]}"; do
  cp "${itp_file}" "${PROD_WORK_DIR}/$(basename "${itp_file}")"
done
if [[ -f "${EQUIL_WORK_DIR}/ligand.prm" ]]; then
  cp "${EQUIL_WORK_DIR}/ligand.prm" "${PROD_WORK_DIR}/ligand.prm"
fi
if [[ -f "${EQUIL_WORK_DIR}/indexlig.ndx" ]]; then
  cp "${EQUIL_WORK_DIR}/indexlig.ndx" "${PROD_WORK_DIR}/indexlig.ndx"
fi

render_mdp "${MD_MDP}" "${PROD_WORK_DIR}/md.mdp"

run_cmd gmx_cmd grompp \
  -f "${PROD_WORK_DIR}/md.mdp" \
  -c "${EQUIL_WORK_DIR}/npt/npt.gro" \
  -t "${EQUIL_WORK_DIR}/npt/npt.cpt" \
  -p "${PROD_WORK_DIR}/topol.top" \
  -o "${PROD_WORK_DIR}/md.tpr" \
  -maxwarn "${GMX_MAXWARN}"

run_cmd run_mdrun_with_fallback "${PROD_WORK_DIR}/md"

mark_step_done "${STEP_NAME}"
log "Finished ${STEP_NAME}"
