#!/usr/bin/env bash
set -Eeuo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PIPELINE_DIR}/config.sh"
source "${PIPELINE_DIR}/lib/common.sh"

STEP_NAME="step1_protein_prep"
if ! step_should_run "${STEP_NAME}"; then
  log "Skipping ${STEP_NAME}; already completed."
  exit 0
fi

init_step_logging "${STEP_NAME}"
prepare_directories
auto_rectify_config
auto_select_ligand_generator
require_command "${GMX_BIN}"

PROTEIN_WORK_DIR="${WORK_DIR}/01_protein"
DOCKING_WORK_DIR="${WORK_DIR}/00_docking"
mkdir -p "${PROTEIN_WORK_DIR}"

CLEAN_PDB="${PROTEIN_WORK_DIR}/protein.cleaned.pdb"
PROCESSED_GRO="${PROTEIN_WORK_DIR}/protein_processed.gro"
PROTEIN_TOP="${PROTEIN_WORK_DIR}/topol.top"
PROTEIN_POSRE="${PROTEIN_WORK_DIR}/posre.itp"
PROTEIN_PDB_OUT="${PROTEIN_WORK_DIR}/protein_processed.pdb"

if [[ -n "${DOCKING_COMPLEX_PDB}" ]]; then
  require_file "${DOCKING_COMPLEX_PDB}"
  [[ -n "${DOCKING_LIGAND_RESNAME}" ]] || die "DOCKING_LIGAND_RESNAME is required when DOCKING_COMPLEX_PDB is provided"
  mkdir -p "${DOCKING_WORK_DIR}"
  run_cmd "${PYTHON_BIN}" "${PIPELINE_DIR}/helpers/split_docking_pdb.py" \
    "${DOCKING_COMPLEX_PDB}" "${CLEAN_PDB}" "${DOCKING_WORK_DIR}/${LIGAND_NAME}.pdb" "${DOCKING_LIGAND_RESNAME}" "${STRIP_HETATM}"
elif [[ "${STRIP_HETATM}" == "yes" ]]; then
  require_file "${PROTEIN_PDB}"
  awk '/^(ATOM  |TER   |END   )/ {print}' "${PROTEIN_PDB}" > "${CLEAN_PDB}"
else
  require_file "${PROTEIN_PDB}"
  cp "${PROTEIN_PDB}" "${CLEAN_PDB}"
fi

run_cmd "${PYTHON_BIN}" "${PIPELINE_DIR}/helpers/normalize_histidine_for_forcefield.py" \
  "${CLEAN_PDB}" --force-field "${FORCE_FIELD}"

log "CMD: ${GMX_BIN} pdb2gmx -f ${CLEAN_PDB} -o ${PROCESSED_GRO} -p ${PROTEIN_TOP} -i ${PROTEIN_POSRE} -ff ${FORCE_FIELD} -water ${WATER_MODEL} -ignh"
PDB2GMX_ARGS=(-f "${CLEAN_PDB}" -o "${PROCESSED_GRO}" -p "${PROTEIN_TOP}" -i "${PROTEIN_POSRE}" -ff "${FORCE_FIELD}" -water "${WATER_MODEL}" -ignh)
if [[ "${FORCE_FIELD,,}" == charmm* ]]; then
  # The bundled CHARMM port's automatic 1MET terminal block is incomplete.
  log "CHARMM termini: standard NH3+ / COO-"
  PDB2GMX_ARGS+=(-ter)
  if ! for ((chain=0; chain<100; chain++)); do printf '1\n0\n'; done | gmx_cmd pdb2gmx "${PDB2GMX_ARGS[@]}"; then
    die "pdb2gmx failed for ${FORCE_FIELD}/${WATER_MODEL}. Check the preceding GROMACS fatal error for the specific residue, termini, or force-field cause; this is not necessarily a water-model error."
  fi
elif ! gmx_cmd pdb2gmx "${PDB2GMX_ARGS[@]}"; then
  die "pdb2gmx failed for ${FORCE_FIELD}/${WATER_MODEL}. Check the preceding GROMACS fatal error for the specific residue, termini, or force-field cause; this is not necessarily a water-model error."
fi

run_cmd gmx_cmd editconf -f "${PROCESSED_GRO}" -o "${PROTEIN_PDB_OUT}"

mark_step_done "${STEP_NAME}"
log "Finished ${STEP_NAME}"
