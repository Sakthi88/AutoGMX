#!/usr/bin/env bash
set -Eeuo pipefail

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "${PIPELINE_DIR}/config.sh"
source "${PIPELINE_DIR}/lib/common.sh"

STEP_NAME="step2_ligand_prep"
if ! step_should_run "${STEP_NAME}"; then
  log "Skipping ${STEP_NAME}; already completed."
  exit 0
fi

init_step_logging "${STEP_NAME}"
prepare_directories
auto_rectify_config
ensure_forcefield_compatibility

LIGAND_WORK_DIR="${WORK_DIR}/02_ligand"
DOCKING_WORK_DIR="${WORK_DIR}/00_docking"
mkdir -p "${LIGAND_WORK_DIR}"

LIGAND_MOL2="${LIGAND_WORK_DIR}/${LIGAND_NAME}.mol2"
LIGAND_PDB="${LIGAND_WORK_DIR}/${LIGAND_NAME}.pdb"
LIGAND_ITP="${LIGAND_WORK_DIR}/ligand.itp"
LIGAND_GRO="${LIGAND_WORK_DIR}/ligand.gro"
LIGAND_POSRE="${LIGAND_WORK_DIR}/posre_ligand.itp"
LIGAND_PRM="${LIGAND_WORK_DIR}/ligand.prm"
LIGAND_META="${LIGAND_WORK_DIR}/ligand_metadata.env"

LIGAND_SOURCE="${LIGAND_INPUT}"
LIGAND_SOURCE_FORMAT="${LIGAND_INPUT_FORMAT}"
if [[ -n "${DOCKING_COMPLEX_PDB}" ]]; then
  LIGAND_SOURCE="${DOCKING_WORK_DIR}/${LIGAND_NAME}.pdb"
  LIGAND_SOURCE_FORMAT="pdb"
fi

generate_ligand_gro_from_source() {
  local source_path="$1"
  local source_format="$2"
  local output_gro="$3"
  local coords_pdb="${LIGAND_WORK_DIR}/${LIGAND_NAME}_coords.pdb"

  case "${source_format}" in
    pdb)
      cp "${source_path}" "${coords_pdb}"
      ;;
    *)
      require_command "${OBABEL_BIN}"
      run_cmd "${OBABEL_BIN}" "${source_path}" -O "${coords_pdb}"
      ;;
  esac

  require_command "${GMX_BIN}"
  run_cmd gmx_cmd editconf -f "${coords_pdb}" -o "${output_gro}"
}

write_ligand_posre() {
  local itp_file="$1"
  local output_file="$2"
  awk '
    BEGIN {print "[ position_restraints ]"; print "; ai funct fcx fcy fcz"}
    /^\[ atoms \]/ {in_atoms=1; next}
    /^\[/ && $0 !~ /^\[ atoms \]/ {in_atoms=0}
    in_atoms && $1 ~ /^[0-9]+$/ {printf "%6d %6d %6d %6d %6d\n", $1, 1, 1000, 1000, 1000}
  ' "${itp_file}" > "${output_file}"
}

if [[ -n "${PREPARED_LIGAND_ITP}" ]]; then
  require_file "${PREPARED_LIGAND_ITP}"
  cp "${PREPARED_LIGAND_ITP}" "${LIGAND_ITP}"
  if [[ -n "${PREPARED_LIGAND_GRO}" ]]; then
    require_file "${PREPARED_LIGAND_GRO}"
    cp "${PREPARED_LIGAND_GRO}" "${LIGAND_GRO}"
  else
    require_file "${LIGAND_SOURCE}"
    generate_ligand_gro_from_source "${LIGAND_SOURCE}" "${LIGAND_SOURCE_FORMAT}" "${LIGAND_GRO}"
  fi
  if [[ -n "${PREPARED_LIGAND_POSRE}" ]]; then
    require_file "${PREPARED_LIGAND_POSRE}"
    cp "${PREPARED_LIGAND_POSRE}" "${LIGAND_POSRE}"
  else
    awk '
      BEGIN {print "[ position_restraints ]"; print "; ai funct fcx fcy fcz"}
      /^\[ atoms \]/ {in_atoms=1; next}
      /^\[/ && $0 !~ /^\[ atoms \]/ {in_atoms=0}
      in_atoms && $1 ~ /^[0-9]+$/ {printf "%6d %6d %6d %6d %6d\n", $1, 1, 1000, 1000, 1000}
    ' "${LIGAND_ITP}" > "${LIGAND_POSRE}"
  fi
  if [[ -n "${PREPARED_LIGAND_PRM}" ]]; then
    require_file "${PREPARED_LIGAND_PRM}"
    cp "${PREPARED_LIGAND_PRM}" "${LIGAND_PRM}"
  fi

  run_cmd "${PYTHON_BIN}" "${PIPELINE_DIR}/helpers/normalize_gro_resname.py" \
    "${LIGAND_GRO}" "${LIGAND_GRO}" "${LIGAND_RESNAME}"

  LIGAND_MOLECULE_NAME_DETECTED="$(extract_moleculetype_name "${LIGAND_ITP}")"
  [[ -n "${LIGAND_MOLECULE_NAME_DETECTED}" ]] || die "Failed to determine ligand [ moleculetype ] from ${LIGAND_ITP}"
  printf 'LIGAND_MOLECULE_NAME=%q\n' "${LIGAND_MOLECULE_NAME_DETECTED}" > "${LIGAND_META}"
  mark_step_done "${STEP_NAME}"
  log "Finished ${STEP_NAME} using prepared ligand topology files"
  exit 0
fi

require_file "${LIGAND_SOURCE}"

if [[ "${LIGAND_PREP_TOOL}" == "acpype" || "${LIGAND_PREP_TOOL}" == "cgenff" ]]; then
  case "${LIGAND_SOURCE_FORMAT}" in
    mol2)
      cp "${LIGAND_SOURCE}" "${LIGAND_MOL2}"
      ;;
    pdb)
      cp "${LIGAND_SOURCE}" "${LIGAND_PDB}"
      require_command "${OBABEL_BIN}"
      run_cmd "${OBABEL_BIN}" "${LIGAND_PDB}" -O "${LIGAND_MOL2}"
      ;;
    *)
      require_command "${OBABEL_BIN}"
      run_cmd "${OBABEL_BIN}" "${LIGAND_SOURCE}" -O "${LIGAND_MOL2}"
      ;;
  esac
fi

case "${LIGAND_PREP_TOOL}" in
  acpype)
    require_command "${ACPYPE_BIN}"
    configure_acpype_runtime_environment
    run_acpype_cmd "${PYTHON_BIN}" "${PIPELINE_DIR}/helpers/check_acpype_runtime.py" "${ACPYPE_BIN}"
    run_cmd "${PYTHON_BIN}" "${PIPELINE_DIR}/helpers/check_ligand_charge.py" \
      "${LIGAND_MOL2}" "${LIGAND_NET_CHARGE}"
    pushd "${LIGAND_WORK_DIR}" >/dev/null
    run_acpype_cmd "${ACPYPE_BIN}" -i "${LIGAND_MOL2}" -b "${LIGAND_NAME}" -a "${ACPYPE_ATOMTYPE}" -n "${LIGAND_NET_CHARGE}"
    ACPYPE_DIR="${LIGAND_NAME}.acpype"
    cp "${ACPYPE_DIR}/${LIGAND_NAME}_GMX.itp" "${LIGAND_ITP}"
    cp "${ACPYPE_DIR}/${LIGAND_NAME}_GMX.gro" "${LIGAND_GRO}"
    if [[ -f "${ACPYPE_DIR}/posre_${LIGAND_NAME}.itp" ]]; then
      cp "${ACPYPE_DIR}/posre_${LIGAND_NAME}.itp" "${LIGAND_POSRE}"
    else
      awk '
        BEGIN {print "[ position_restraints ]"; print "; ai funct fcx fcy fcz"}
        /^\[ atoms \]/ {in_atoms=1; next}
        /^\[/ && $0 !~ /^\[ atoms \]/ {in_atoms=0}
        in_atoms && $1 ~ /^[0-9]+$/ {printf "%6d %6d %6d %6d %6d\n", $1, 1, 1000, 1000, 1000}
      ' "${LIGAND_ITP}" > "${LIGAND_POSRE}"
    fi
    popd >/dev/null
    ;;
  cgenff)
    require_command "${CGENFF_BIN}"
    require_command "${CGENFF_PYTHON}"
    require_command "${CGENFF_CONVERTER}"
    pushd "${LIGAND_WORK_DIR}" >/dev/null
    run_cmd "${CGENFF_BIN}" "${LIGAND_MOL2}" > "${LIGAND_RESNAME}.str"
    run_cmd "${CGENFF_PYTHON}" "${CGENFF_CONVERTER}" "${LIGAND_RESNAME}" "${LIGAND_MOL2}" "${LIGAND_RESNAME}.str" "${CGENFF_FORCEFIELD_DIR}"
    CGENFF_BASE="$(printf '%s' "${LIGAND_RESNAME}" | tr '[:upper:]' '[:lower:]')"
    cp "${CGENFF_BASE}.itp" "${LIGAND_ITP}"
    if [[ -f "${CGENFF_BASE}.prm" ]]; then
      cp "${CGENFF_BASE}.prm" "${LIGAND_WORK_DIR}/ligand.prm"
    fi
    run_cmd gmx_cmd editconf -f "${CGENFF_BASE}_ini.pdb" -o "${LIGAND_GRO}"
    awk '
      BEGIN {print "[ position_restraints ]"; print "; ai funct fcx fcy fcz"}
      /^\[ atoms \]/ {in_atoms=1; next}
      /^\[/ && $0 !~ /^\[ atoms \]/ {in_atoms=0}
      in_atoms && $1 ~ /^[0-9]+$/ {printf "%6d %6d %6d %6d %6d\n", $1, 1, 1000, 1000, 1000}
    ' "${LIGAND_ITP}" > "${LIGAND_POSRE}"
    popd >/dev/null
    ;;
  pypolybuilder)
    require_command "${PYPOLYBUILDER_BIN}"
    run_cmd "${PYTHON_BIN}" "${PIPELINE_DIR}/helpers/run_pypolybuilder.py" \
      --bin "${PYPOLYBUILDER_BIN}" \
      --work-dir "${LIGAND_WORK_DIR}" \
      --mode "${PYPOLYBUILDER_MODE}" \
      --params "${PYPOLYBUILDER_PARAMS}" \
      --forcefield-path "${PYPOLYBUILDER_FORCEFIELD_PATH}" \
      --force-field "${FORCE_FIELD}" \
      --core "${PYPOLYBUILDER_CORE}" \
      --terminal "${PYPOLYBUILDER_TERMINAL}" \
      --inter "${PYPOLYBUILDER_INTERMEDIATE}" \
      --bbs "${PYPOLYBUILDER_BBS}" \
      --connections "${PYPOLYBUILDER_CONNECTIONS}" \
      --ngen "${PYPOLYBUILDER_NGEN}" \
      --nsteps "${PYPOLYBUILDER_NSTEPS}" \
      --ngenga "${PYPOLYBUILDER_NGENGA}" \
      --npop "${PYPOLYBUILDER_NPOP}" \
      --name "${LIGAND_NAME}" \
      --output-itp "${LIGAND_ITP}" \
      --output-gro "${LIGAND_GRO}" \
      --output-top "${LIGAND_WORK_DIR}/ligand.top" \
      --extra-args "${PYPOLYBUILDER_EXTRA_ARGS}"
    awk '
      BEGIN {print "[ position_restraints ]"; print "; ai funct fcx fcy fcz"}
      /^\[ atoms \]/ {in_atoms=1; next}
      /^\[/ && $0 !~ /^\[ atoms \]/ {in_atoms=0}
      in_atoms && $1 ~ /^[0-9]+$/ {printf "%6d %6d %6d %6d %6d\n", $1, 1, 1000, 1000, 1000}
    ' "${LIGAND_ITP}" > "${LIGAND_POSRE}"
    ;;
  *)
    die "Unsupported ligand preparation tool: ${LIGAND_PREP_TOOL}"
    ;;
esac

run_cmd "${PYTHON_BIN}" "${PIPELINE_DIR}/helpers/normalize_gro_resname.py" \
  "${LIGAND_GRO}" "${LIGAND_GRO}" "${LIGAND_RESNAME}"

LIGAND_MOLECULE_NAME_DETECTED="$(extract_moleculetype_name "${LIGAND_ITP}")"
[[ -n "${LIGAND_MOLECULE_NAME_DETECTED}" ]] || die "Failed to determine ligand [ moleculetype ] from ${LIGAND_ITP}"
printf 'LIGAND_MOLECULE_NAME=%q\n' "${LIGAND_MOLECULE_NAME_DETECTED}" > "${LIGAND_META}"

mark_step_done "${STEP_NAME}"
log "Finished ${STEP_NAME}"
