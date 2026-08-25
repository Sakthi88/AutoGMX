#!/usr/bin/env bash
set -Eeuo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PIPELINE_DIR}/config.sh"
source "${PIPELINE_DIR}/lib/common.sh"

STEP_NAME="step11_mmpbsa"
if ! step_should_run "${STEP_NAME}"; then
  log "Skipping ${STEP_NAME}; already completed."
  exit 0
fi

init_step_logging "${STEP_NAME}"
prepare_directories
load_runtime_overrides

MMBSA_DIR="${ANALYSIS_DIR}/mmpbsa"
mkdir -p "${MMBSA_DIR}"

TPR_FILE="${RESULTS_DIR}/production/md.tpr"
XTC_FILE="${RESULTS_DIR}/production/md.xtc"
GRO_FILE="${RESULTS_DIR}/production/md.gro"

require_file "${TPR_FILE}"
require_file "${XTC_FILE}"

# ── Determine ligand selection ────────────────────────────────────────────
# Use index group from step10 analysis (created by make_protein_ligand_index)
INDEX_FILE="${ANALYSIS_DIR}/index.ndx"
LIGAND_SEL="group ${ANALYSIS_LIGAND_GROUP:-Ligand}"
RECEPTOR_SEL="protein"

log "MM-PBSA configuration:"
log "  TPR: ${TPR_FILE}"
log "  XTC: ${XTC_FILE}"
log "  Receptor: ${RECEPTOR_SEL}"
log "  Ligand: ${LIGAND_SEL}"
log "  Output: ${MMBSA_DIR}"

# ── Step A: Run MM-PBSA calculation ───────────────────────────────────────
log "CMD: ${PYTHON_BIN} ${PIPELINE_DIR}/helpers/mmpbsa_calculation.py"

MMBSA_ARGS=(
  -t "${XTC_FILE}"
  -s "${TPR_FILE}"
  -o "${MMBSA_DIR}"
  -n "${INDEX_FILE}"
  --gmx-bin "${GMX_BIN}"
  --receptor-selection "${RECEPTOR_SEL}"
  --ligand-selection "${LIGAND_SEL}"
  --frame-step "${MMBSA_FRAME_STEP:-10}"
  --temperature "${TEMPERATURE_K:-300}"
  --ionic-strength "${MMBSA_IONIC_STRENGTH:-0.15}"
  --surface-tension "${MMBSA_SURFACE_TENSION:-0.0072}"
)

if [[ "${MMBSA_DECOMPOSE:-no}" == "yes" ]]; then
  MMBSA_ARGS+=(--decompose)
fi

if [[ "${MMBSA_USE_APBS:-yes}" == "no" ]]; then
  MMBSA_ARGS+=(--no-apbs)
fi

"${PYTHON_BIN}" "${PIPELINE_DIR}/helpers/mmpbsa_calculation.py" "${MMBSA_ARGS[@]}"

# ── Step B: Generate PDF report ───────────────────────────────────────────
log "CMD: ${PYTHON_BIN} ${PIPELINE_DIR}/helpers/generate_report.py"

PDF_PATH="${RESULTS_DIR}/md_report.pdf"

"${PYTHON_BIN}" "${PIPELINE_DIR}/helpers/generate_report.py" \
  -a "${ANALYSIS_DIR}" \
  -o "${PDF_PATH}" \
  -n "${PROJECT_NAME}"

log "PDF report: ${PDF_PATH}"

mark_step_done "${STEP_NAME}"
log "Finished ${STEP_NAME}"
