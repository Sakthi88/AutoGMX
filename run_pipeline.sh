#!/usr/bin/env bash
set -Eeuo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${PIPELINE_DIR}/config.sh"

JSON_CONFIG=""
START_STEP=1
END_STEP=11

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      source "$2"
      shift 2
      ;;
    --json)
      JSON_CONFIG="$2"
      shift 2
      ;;
    --from)
      START_STEP="$2"
      shift 2
      ;;
    --to)
      END_STEP="$2"
      shift 2
      ;;
    -h|--help)
      cat <<'EOF'
Usage: bash run_pipeline.sh [--config custom.sh] [--json overrides.json] [--from N] [--to N]
EOF
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

if [[ -n "${JSON_CONFIG}" ]]; then
  eval "$("${PYTHON_BIN}" "${PIPELINE_DIR}/helpers/json_to_env.py" "${JSON_CONFIG}")"
fi

source "${PIPELINE_DIR}/lib/common.sh"
prepare_directories
auto_rectify_config
auto_select_ligand_generator
print_runtime_summary

steps=(
  "1:${PIPELINE_DIR}/steps/step1_protein_prep.sh"
  "2:${PIPELINE_DIR}/steps/step2_ligand_prep.sh"
  "3:${PIPELINE_DIR}/steps/step3_assemble_complex.sh"
  "4:${PIPELINE_DIR}/steps/step4_define_box.sh"
  "5:${PIPELINE_DIR}/steps/step5_solvate.sh"
  "6:${PIPELINE_DIR}/steps/step6_add_ions.sh"
  "7:${PIPELINE_DIR}/steps/step7_minimize.sh"
  "8:${PIPELINE_DIR}/steps/step8_equilibrate.sh"
  "9:${PIPELINE_DIR}/steps/step9_production.sh"
  "10:${PIPELINE_DIR}/steps/step10_analysis.sh"
  "11:${PIPELINE_DIR}/steps/step11_mmpbsa.sh"
)

for entry in "${steps[@]}"; do
  step_number="${entry%%:*}"
  step_script="${entry#*:}"
  if (( step_number < START_STEP || step_number > END_STEP )); then
    continue
  fi
  bash "${step_script}"
done

log "Pipeline completed successfully."
