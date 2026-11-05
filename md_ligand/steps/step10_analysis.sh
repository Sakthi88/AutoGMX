#!/usr/bin/env bash
set -Eeuo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PIPELINE_DIR}/config.sh"
source "${PIPELINE_DIR}/lib/common.sh"

STEP_NAME="step10_analysis"
if ! step_should_run "${STEP_NAME}"; then
  log "Skipping ${STEP_NAME}; already completed."
  exit 0
fi

init_step_logging "${STEP_NAME}"
prepare_directories

PROD_WORK_DIR="$(production_primary_dir)"
INDEX_FILE="${ANALYSIS_DIR}/index.ndx"
CENTERED_TRAJ="${PROD_WORK_DIR}/md_center.xtc"
NOJUMP_TRAJ="${PROD_WORK_DIR}/md_nojump.xtc"
TRAJ_INPUT="${PROD_WORK_DIR}/md.xtc"

require_file "${PROD_WORK_DIR}/md.tpr"
require_file "${PROD_WORK_DIR}/md.xtc"
require_file "${PROD_WORK_DIR}/md.gro"

if [[ -f "${PROD_WORK_DIR}/indexlig.ndx" ]]; then
  log "Using imported index file: ${PROD_WORK_DIR}/indexlig.ndx"
  cp "${PROD_WORK_DIR}/indexlig.ndx" "${INDEX_FILE}"
else
  make_protein_ligand_index "${PROD_WORK_DIR}/md.gro" "${INDEX_FILE}"
fi

case "${MD_CENTER_NOJUMP}" in
  yes)
    printf "%s\n" "${MD_CENTER_OUTPUT_GROUP}" | \
      gmx_cmd trjconv -s "${PROD_WORK_DIR}/md.tpr" -f "${PROD_WORK_DIR}/md.xtc" -n "${INDEX_FILE}" -o "${NOJUMP_TRAJ}" -pbc nojump
    TRAJ_INPUT="${NOJUMP_TRAJ}"
    ;;
  no)
    ;;
  *)
    die "MD_CENTER_NOJUMP must be either yes or no"
    ;;
esac

printf "%s\n%s\n" "${MD_CENTER_GROUP}" "${MD_CENTER_OUTPUT_GROUP}" | \
  gmx_cmd trjconv -s "${PROD_WORK_DIR}/md.tpr" -f "${TRAJ_INPUT}" -n "${INDEX_FILE}" -o "${CENTERED_TRAJ}" -center -pbc mol -ur compact

run_cmd gmx_cmd editconf -f "${PROD_WORK_DIR}/md.gro" -o "${ANALYSIS_DIR}/md_final.pdb"

printf "%s\n%s\n" "${ANALYSIS_FIT_GROUP}" "${ANALYSIS_RMS_GROUP}" | \
  gmx_cmd rms -s "${PROD_WORK_DIR}/md.tpr" -f "${CENTERED_TRAJ}" -n "${INDEX_FILE}" -o "${ANALYSIS_DIR}/rmsd_protein.xvg"

printf "Protein_Ligand\nProtein_Ligand\n" | \
  gmx_cmd rms -s "${PROD_WORK_DIR}/md.tpr" -f "${CENTERED_TRAJ}" -n "${INDEX_FILE}" -o "${ANALYSIS_DIR}/rmsd_protein_ligand.xvg"

printf "%s\n%s\n" "${ANALYSIS_FIT_GROUP}" "${ANALYSIS_LIGAND_GROUP}" | \
  gmx_cmd rms -s "${PROD_WORK_DIR}/md.tpr" -f "${CENTERED_TRAJ}" -n "${INDEX_FILE}" -o "${ANALYSIS_DIR}/rmsd_ligand.xvg"

printf "%s\n" "${ANALYSIS_RMS_GROUP}" | \
  gmx_cmd rmsf -s "${PROD_WORK_DIR}/md.tpr" -f "${CENTERED_TRAJ}" -n "${INDEX_FILE}" -o "${ANALYSIS_DIR}/rmsf_protein.xvg" -res

printf "%s\n" "${ANALYSIS_RMS_GROUP}" | \
  gmx_cmd gyrate -s "${PROD_WORK_DIR}/md.tpr" -f "${CENTERED_TRAJ}" -n "${INDEX_FILE}" -o "${ANALYSIS_DIR}/gyration_protein.xvg"

HBOND_XVG="${ANALYSIS_DIR}/hbonds_protein_ligand.xvg"
HBOND_LOG="${ANALYSIS_DIR}/hbond_attempt.log"

# gmx hbond (GROMACS >=2024) only treats N and O as donor/acceptor elements
# by default. Ligands whose polar groups are S/F/Cl/etc. (or whose element
# metadata isn't set the way the tool expects) get a hard
# "no donors AND has no acceptors" abort even though real H-bonds may exist.
# Try progressively more permissive methods before giving up.
run_hbond_attempt() {
  local label="$1"; shift
  log "Trying Protein-Ligand hbond analysis (${label})..."
  printf "%s\n%s\n" "${ANALYSIS_RMS_GROUP}" "${ANALYSIS_LIGAND_GROUP}" | \
    gmx_cmd hbond -s "${PROD_WORK_DIR}/md.tpr" -f "${CENTERED_TRAJ}" -n "${INDEX_FILE}" \
      -num "${HBOND_XVG}" "$@" >"${HBOND_LOG}" 2>&1
}

if run_hbond_attempt "default N/O elements"; then
  :
elif grep -qi "no donors\|no acceptors" "${HBOND_LOG}"; then
  log "Default N/O element set found nothing on the Ligand selection; retrying with a wider element set (N,O,S,F,Cl,Br,I)"
  if run_hbond_attempt "wide element set" -de "N O S F Cl Br I" -ae "N O S F Cl Br I"; then
    :
  elif grep -qi "no donors\|no acceptors" "${HBOND_LOG}"; then
    log "Wider element set still found nothing; falling back to gmx hbond-legacy (name-based donor/acceptor rules)"
    if printf "%s\n%s\n" "${ANALYSIS_RMS_GROUP}" "${ANALYSIS_LIGAND_GROUP}" | \
        gmx_cmd hbond-legacy -s "${PROD_WORK_DIR}/md.tpr" -f "${CENTERED_TRAJ}" -n "${INDEX_FILE}" \
          -num "${HBOND_XVG}" >"${HBOND_LOG}" 2>&1; then
      :
    else
      log "WARNING: Protein-Ligand hbond analysis failed with all methods (ligand likely has no true polar heteroatoms); skipping hbond plot. Details: ${HBOND_LOG}"
    fi
  else
    log "WARNING: Protein-Ligand hbond analysis failed; skipping hbond plot. Details: ${HBOND_LOG}"
  fi
else
  log "WARNING: Protein-Ligand hbond analysis failed; skipping hbond plot. Details: ${HBOND_LOG}"
fi

printf "Potential\nTemperature\nPressure\nDensity\n0\n" | \
  gmx_cmd energy -f "${PROD_WORK_DIR}/md.edr" -o "${ANALYSIS_DIR}/energies.xvg"

printf "Potential\n0\n" | \
  gmx_cmd energy -f "${PROD_WORK_DIR}/md.edr" -o "${ANALYSIS_DIR}/potential.xvg"

printf "Temperature\n0\n" | \
  gmx_cmd energy -f "${PROD_WORK_DIR}/md.edr" -o "${ANALYSIS_DIR}/temperature.xvg"

printf "Pressure\n0\n" | \
  gmx_cmd energy -f "${PROD_WORK_DIR}/md.edr" -o "${ANALYSIS_DIR}/pressure.xvg"

printf "Density\n0\n" | \
  gmx_cmd energy -f "${PROD_WORK_DIR}/md.edr" -o "${ANALYSIS_DIR}/density.xvg"

if command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
  if "${PYTHON_BIN}" -c "import matplotlib" >/dev/null 2>&1; then
    plot_optional() {
      local input_xvg="$1"
      local title="$2"
      local output_png="$3"
      if [[ ! -s "${input_xvg}" ]]; then
        log "Skipping plot; missing XVG data: ${input_xvg}"
        return 0
      fi
      log "CMD: ${PYTHON_BIN} ${PIPELINE_DIR}/helpers/plot_analysis.py ${input_xvg} ${title} ${output_png}"
      if ! "${PYTHON_BIN}" "${PIPELINE_DIR}/helpers/plot_analysis.py" "${input_xvg}" "${title}" "${output_png}"; then
        log "Plot generation failed for ${input_xvg}; continuing because MD production outputs are complete."
      fi
    }
    plot_optional "${ANALYSIS_DIR}/rmsd_protein.xvg" "Protein RMSD" "${ANALYSIS_DIR}/rmsd_protein.png"
    plot_optional "${ANALYSIS_DIR}/rmsd_ligand.xvg" "Ligand RMSD" "${ANALYSIS_DIR}/rmsd_ligand.png"
    if [[ -s "${ANALYSIS_DIR}/rmsd_protein.xvg" && -s "${ANALYSIS_DIR}/rmsd_ligand.xvg" ]]; then
      log "CMD: ${PYTHON_BIN} ${PIPELINE_DIR}/helpers/plot_analysis.py --rmsd-comparison ${ANALYSIS_DIR}/rmsd_protein.xvg ${ANALYSIS_DIR}/rmsd_ligand.xvg ${ANALYSIS_DIR}/rmsd_protein_ligand.png"
      if ! "${PYTHON_BIN}" "${PIPELINE_DIR}/helpers/plot_analysis.py" --rmsd-comparison \
        "${ANALYSIS_DIR}/rmsd_protein.xvg" "${ANALYSIS_DIR}/rmsd_ligand.xvg" \
        "${ANALYSIS_DIR}/rmsd_protein_ligand.png"; then
        log "Combined Protein/Ligand RMSD plot generation failed; continuing because MD production outputs are complete."
      fi
    else
      log "Skipping combined Protein/Ligand RMSD plot; protein or ligand RMSD XVG data is missing."
    fi
    plot_optional "${ANALYSIS_DIR}/rmsf_protein.xvg" "Protein RMSF" "${ANALYSIS_DIR}/rmsf_protein.png"
    plot_optional "${ANALYSIS_DIR}/gyration_protein.xvg" "Protein Radius of Gyration" "${ANALYSIS_DIR}/gyration_protein.png"
    plot_optional "${ANALYSIS_DIR}/hbonds_protein_ligand.xvg" "Protein-Ligand Hydrogen Bonds" "${ANALYSIS_DIR}/hbonds_protein_ligand.png"
    plot_optional "${ANALYSIS_DIR}/potential.xvg" "Potential Energy" "${ANALYSIS_DIR}/potential.png"
    plot_optional "${ANALYSIS_DIR}/temperature.xvg" "Temperature" "${ANALYSIS_DIR}/temperature.png"
    plot_optional "${ANALYSIS_DIR}/pressure.xvg" "Pressure" "${ANALYSIS_DIR}/pressure.png"
    plot_optional "${ANALYSIS_DIR}/density.xvg" "Density" "${ANALYSIS_DIR}/density.png"
  else
    log "matplotlib not available; skipping PNG plot generation."
  fi
fi

if [[ "${RUN_FREE_ENERGY}" == "yes" ]]; then
  FREE_ENERGY_DIR="${ANALYSIS_DIR}/free_energy"
  mkdir -p "${FREE_ENERGY_DIR}"
  log "CMD: ${PYTHON_BIN} ${PIPELINE_DIR}/helpers/pca_free_energy.py --top ${PROD_WORK_DIR}/md.tpr --traj ${CENTERED_TRAJ} --outdir ${FREE_ENERGY_DIR}"
  if ! "${PYTHON_BIN}" "${PIPELINE_DIR}/helpers/pca_free_energy.py" \
    --top "${PROD_WORK_DIR}/md.tpr" \
    --fallback-top "${PROD_WORK_DIR}/md.gro" \
    --traj "${CENTERED_TRAJ}" \
    --outdir "${FREE_ENERGY_DIR}" \
    --selection "${FREE_ENERGY_SELECTION}" \
    --temperature "${FREE_ENERGY_TEMPERATURE_K}" \
    --components "${FREE_ENERGY_COMPONENTS}" \
    --bins "${FREE_ENERGY_BINS}" \
    --step "${FREE_ENERGY_FRAME_STEP}"; then
    log "PCA free-energy analysis was skipped or incomplete; see ${FREE_ENERGY_DIR}/analysis_summary.txt"
  fi
  mark_step_done "free_energy"
fi

if [[ "${RUN_UMBRELLA}" == "yes" ]]; then
  UMBRELLA_DIR="${ANALYSIS_DIR}/umbrella_sampling"
  mkdir -p "${UMBRELLA_DIR}"
  {
    printf 'GROMACS WHAM Umbrella Sampling Analysis\n'
    printf 'Status: waiting for umbrella-window inputs\n'
    printf 'Generated: %s\n' "$(date '+%F %T')"
    printf 'Source checked: %s\n\n' "${PROD_WORK_DIR}"
    printf 'This linked task cannot run from ordinary production MD files alone.\n'
    printf 'Upload matching umbrella-window .tpr files plus pullx*.xvg or pullf*.xvg files in the Umbrella tab.\n'
    printf 'The app will then create tpr-files.dat and pullx-files.dat or pullf-files.dat and run gmx wham.\n'
    printf 'Reference: https://manual.gromacs.org/current/onlinehelp/gmx-wham.html\n'
  } > "${UMBRELLA_DIR}/analysis_summary.txt"
fi

mark_step_done "${STEP_NAME}"
log "Finished ${STEP_NAME}"
