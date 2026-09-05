#!/usr/bin/env bash
set -Eeuo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PIPELINE_DIR}/config.sh"
source "${PIPELINE_DIR}/lib/common.sh"

STEP_NAME="step4_define_box"
if ! step_should_run "${STEP_NAME}"; then
  log "Skipping ${STEP_NAME}; already completed."
  exit 0
fi

init_step_logging "${STEP_NAME}"
prepare_directories

COMPLEX_WORK_DIR="${WORK_DIR}/03_complex"
BOX_WORK_DIR="${WORK_DIR}/04_box"
mkdir -p "${BOX_WORK_DIR}"

require_file "${COMPLEX_WORK_DIR}/complex.gro"

cp "${COMPLEX_WORK_DIR}/topol.top" "${BOX_WORK_DIR}/topol.top"
# Copy every auxiliary *.itp include topol.top may reference: ligand.itp,
# posre_ligand.itp, and the protein's own includes (posre.itp for a
# single-chain protein, or posre_Protein_chain_A.itp/
# topol_Protein_chain_A.itp/... for multi-chain proteins). topol.top's
# #include list already matches whichever files step3 actually produced.
shopt -s nullglob
itp_files=("${COMPLEX_WORK_DIR}"/*.itp)
shopt -u nullglob
if [[ ${#itp_files[@]} -eq 0 ]]; then
  die "No *.itp include files found in ${COMPLEX_WORK_DIR}"
fi
for itp_file in "${itp_files[@]}"; do
  cp "${itp_file}" "${BOX_WORK_DIR}/$(basename "${itp_file}")"
done
if [[ -f "${COMPLEX_WORK_DIR}/ligand.prm" ]]; then
  cp "${COMPLEX_WORK_DIR}/ligand.prm" "${BOX_WORK_DIR}/ligand.prm"
fi

run_cmd gmx_cmd editconf \
  -f "${COMPLEX_WORK_DIR}/complex.gro" \
  -o "${BOX_WORK_DIR}/boxed.gro" \
  -bt "${BOX_TYPE}" \
  -d "${BOX_DISTANCE_NM}" \
  -c

mark_step_done "${STEP_NAME}"
log "Finished ${STEP_NAME}"
