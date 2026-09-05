#!/usr/bin/env python3
import os
import pathlib
import shlex
import sys


COMPATIBLE_SETS = {
    "cgenff": [
        ("charmm36-feb2026_cgenff-5.0", "tip3p"),
        ("charmm36-jul2022", "tip3p"),
        ("charmm36", "tip3p"),
        ("charmm36-mar2019", "tip3p"),
        ("charmm27", "tip3p"),
        ("oplsaa", "tip4p"),
    ],
    "acpype": [
        ("amber99sb-ildn", "tip3p"),
        ("amber99sb-ildn", "spce"),
        ("amber14sb", "tip3p"),
        ("amber99sb", "tip3p"),
        ("amber99", "tip3p"),
        ("amber96", "tip3p"),
        ("amber94", "tip3p"),
        ("amber03", "tip3p"),
        ("amberGS", "tip3p"),
        ("charmm36-feb2026_cgenff-5.0", "tip3p"),
        ("charmm36-jul2022", "tip3p"),
        ("charmm36-mar2019", "tip3p"),
        ("charmm36", "tip3p"),
        ("charmm27", "tip3p"),
        ("oplsaa", "tip4p"),
        ("oplsaa", "spce"),
    ],
    "atb": [
        ("gromos54a7", "spc"),
        ("gromos54a7", "spce"),
        ("gromos53a6", "spc"),
        ("gromos53a5", "spc"),
        ("gromos45a3", "spc"),
        ("gromos43a2", "spc"),
        ("gromos43a1", "spc"),
        ("amber99sb-ildn", "tip3p"),
        ("amber14sb", "tip3p"),
        ("charmm36", "tip3p"),
        ("oplsaa", "tip4p"),
    ],
}

WATER_BY_FAMILY = {
    "amber": {"tip3p", "tip3p_original", "spce", "spc", "tip4p", "tip4pew", "tip5p"},
    "charmm": {"tip3p", "tip3p_original", "spce", "spc", "tip4p", "tip4pew", "tip5p"},
    "gromos": {"spc", "spce"},
    "oplsaa": {"tip4p", "tip3p", "spce", "spc", "tip4pew"},
}

# Use one documented water model for each force-field family.  This avoids
# creating jobs that silently drift between water models while pdb2gmx retries
# force fields, and ensures solvation uses the same model that pdb2gmx chose.
CANONICAL_WATER_BY_FAMILY = {
    "amber": "tip3p",
    "charmm": "tip3p",
    "gromos": "spc",
    "oplsaa": "tip4p",
}


def quote(value):
    return shlex.quote(str(value))


def family(force_field):
    lowered = force_field.lower()
    for item in WATER_BY_FAMILY:
        if lowered.startswith(item):
            return item
    return ""


def forcefield_exists(pipeline_dir, force_field):
    if not force_field:
        return False
    local_dir = pathlib.Path(pipeline_dir, f"{force_field}.ff")
    if local_dir.exists():
        return True
    gmxlib = os.environ.get("GMXLIB")
    if gmxlib and pathlib.Path(gmxlib, f"{force_field}.ff").exists():
        return True
    # Check GROMACS default data directory
    gmx_executable = os.environ.get("GMX_BIN", "gmx")
    try:
        import subprocess
        result = subprocess.run(
            [gmx_executable, "--version"],
            capture_output=True,
            text=True,
            timeout=5,
        )
        for line in result.stdout.splitlines():
            if "Data prefix:" in line:
                data_prefix = line.split("Data prefix:", 1)[1].strip()
                gmx_data = pathlib.Path(data_prefix, "share", "gromacs", "top", f"{force_field}.ff")
                if gmx_data.exists():
                    return True
    except Exception:
        pass
    return False


def first_available(tool, pipeline_dir):
    for force_field, water_model in COMPATIBLE_SETS[tool]:
        if forcefield_exists(pipeline_dir, force_field):
            return force_field, water_model
    return COMPATIBLE_SETS[tool][0]


def is_suitable_force_field(tool, force_field):
    ff_family = family(force_field)
    if tool == "cgenff":
        return ff_family == "charmm"
    if tool == "acpype":
        # ACPYPE is the app's available ligand generator for the AMBER,
        # CHARMM, and OPLS choices listed in COMPATIBLE_SETS.
        return ff_family in {"amber", "charmm", "oplsaa"}
    if tool == "pypolybuilder":
        return ff_family == "gromos"
    if tool == "atb":
        # ATB supports GROMOS, AMBER, CHARMM, and OPLS
        return ff_family in {"gromos", "amber", "charmm", "oplsaa"}
    return False


def canonical_water_model(force_field):
    return CANONICAL_WATER_BY_FAMILY.get(family(force_field))


def main():
    pipeline_dir = pathlib.Path(os.environ.get("PIPELINE_DIR", ".")).resolve()
    tool = os.environ.get("LIGAND_PREP_TOOL", "acpype").lower()
    force_field = os.environ.get("FORCE_FIELD", "amber99sb-ildn")
    water_model = os.environ.get("WATER_MODEL", "tip3p").lower()
    auto_rectify = os.environ.get("AUTO_RECTIFY", "yes").lower()
    out_path = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else None

    prepared_topology = tool == "prepared_topology"
    updates = {}
    notes = []

    # A supplied topology already contains force-field-specific ligand parameters.
    # Do not silently replace its force field based on a topology-generator heuristic.
    if not prepared_topology:
        if tool not in COMPATIBLE_SETS:
            updates["LIGAND_PREP_TOOL"] = "acpype"
            tool = "acpype"
            notes.append("Unsupported ligand-prep tool was replaced with acpype.")

        if auto_rectify in {"yes", "true", "1", "auto"}:
            # The selected force field is authoritative: report an invalid
            # choice instead of silently substituting the Amber default.
            if not forcefield_exists(pipeline_dir, force_field):
                raise SystemExit(
                    f"Selected force field '{force_field}' is not installed for this GROMACS setup. "
                    "Install its .ff directory or choose an installed option."
                )
            if not is_suitable_force_field(tool, force_field):
                raise SystemExit(
                    f"Selected force field '{force_field}' is not supported by ligand preparation tool '{tool}'. "
                    "Choose a compatible ligand tool or force field; the selection was not changed."
                )

            supported_water_models = WATER_BY_FAMILY.get(family(force_field), set())
            if water_model not in supported_water_models:
                allowed = ", ".join(sorted(supported_water_models)) or "none"
                raise SystemExit(
                    f"Selected water model '{water_model}' is not supported by force field '{force_field}'. "
                    f"Choose one of: {allowed}. The selection was not changed."
                )

        # Keep the resolved settings for later pipeline steps, which start in
        # fresh shells and otherwise reload the original job configuration.
        # Persisting all three also keeps an ACPYPE fallback internally
        # consistent from protein preparation through production MD.
        updates.setdefault("FORCE_FIELD", force_field)
        updates.setdefault("WATER_MODEL", water_model)
        updates["LIGAND_PREP_TOOL"] = tool

    if updates.get("FORCE_FIELD", force_field).lower().startswith("charmm"):
        ff_dir = pathlib.Path(pipeline_dir, f"{updates.get('FORCE_FIELD', force_field)}.ff")
        updates["CGENFF_FORCEFIELD_DIR"] = str(ff_dir)

    lines = ["# Generated by helpers/auto_rectify_config.py"]
    for key, value in updates.items():
        lines.append(f"export {key}={quote(value)}")
    for note in notes:
        lines.append(f"# {note}")
    text = "\n".join(lines) + "\n"

    if out_path:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text, encoding="utf-8", newline="\n")
    print(text, end="")


if __name__ == "__main__":
    main()
