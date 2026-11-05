#!/usr/bin/env bash
set -Eeuo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PIPELINE_DIR}/config.sh"
source "${PIPELINE_DIR}/lib/common.sh"
EQUIL_WORK_DIR="${RESULTS_DIR}/equilibration"; PROD_ROOT="${RESULTS_DIR}/production"
require_file "${EQUIL_WORK_DIR}/npt/npt.gro"; require_file "${EQUIL_WORK_DIR}/topol.top"; mkdir -p "${PROD_ROOT}"
for ((replica=1; replica<=ENSEMBLE_REPLICAS; replica++)); do
  printf -v replica_dir '%s/replica_%03d' "${PROD_ROOT}" "${replica}"; seed=$((ENSEMBLE_BASE_SEED + replica - 1)); mkdir -p "${replica_dir}"
  cp "${EQUIL_WORK_DIR}/topol.top" "${replica_dir}/topol.top"; shopt -s nullglob
  for source_file in "${EQUIL_WORK_DIR}"/*.itp; do cp "${source_file}" "${replica_dir}/$(basename "${source_file}")"; done
  shopt -u nullglob; [[ -f "${EQUIL_WORK_DIR}/ligand.prm" ]] && cp "${EQUIL_WORK_DIR}/ligand.prm" "${replica_dir}/ligand.prm"; [[ -f "${EQUIL_WORK_DIR}/indexlig.ndx" ]] && cp "${EQUIL_WORK_DIR}/indexlig.ndx" "${replica_dir}/indexlig.ndx"
  if [[ -f "${replica_dir}/md.cpt" && "${CHECKPOINT_POLICY}" != "off" ]]; then validate_checkpoint "${replica_dir}/md.cpt"; require_file "${replica_dir}/md.tpr"; else
    MD_CONTINUATION=no MD_GEN_VEL=yes MD_GEN_SEED="${seed}"; export MD_CONTINUATION MD_GEN_VEL MD_GEN_SEED; render_mdp "${MD_MDP}" "${replica_dir}/md.mdp"
    run_cmd gmx_cmd grompp -f "${replica_dir}/md.mdp" -c "${EQUIL_WORK_DIR}/npt/npt.gro" -p "${replica_dir}/topol.top" -o "${replica_dir}/md.tpr" -maxwarn "${GMX_MAXWARN}"
  fi
  run_cmd run_mdrun_with_fallback "${replica_dir}/md"
done
printf '{"replicas": %s, "base_seed": %s, "gromacs_version": "%s", "mode": "independent production velocities from common NPT structure"}\n' "${ENSEMBLE_REPLICAS}" "${ENSEMBLE_BASE_SEED}" "$(detect_gromacs_version)" > "${PROD_ROOT}/ensemble_manifest.json"
mark_step_done step9_production; log "Finished production ensemble; manifest: ${PROD_ROOT}/ensemble_manifest.json"
