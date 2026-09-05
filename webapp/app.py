#!/usr/bin/env python3
import json
import os
import re
import signal
import shutil
import subprocess
import sys
import uuid
import zipfile
from datetime import datetime
from pathlib import Path

from flask import Flask, abort, flash, redirect, render_template, request, send_file, url_for
from werkzeug.exceptions import RequestEntityTooLarge


APP_DIR = Path(__file__).resolve().parent
PIPELINE_DIR = APP_DIR.parent
JOBS_DIR = PIPELINE_DIR / "web_jobs"
JOBS_DIR.mkdir(parents=True, exist_ok=True)

app = Flask(__name__)
app.config["SECRET_KEY"] = os.environ.get("MD_WEBAPP_SECRET", "local-md-demo-secret")
app.config["MAX_CONTENT_LENGTH"] = int(os.environ.get("MD_WEBAPP_MAX_UPLOAD_MB", "5120")) * 1024 * 1024

STEP_LABELS = [
    ("step1_protein_prep", "Protein preparation"),
    ("step2_ligand_prep", "Ligand preparation"),
    ("step3_assemble_complex", "Complex assembly"),
    ("step4_define_box", "Box definition"),
    ("step5_solvate", "Solvation"),
    ("step6_add_ions", "Ion addition"),
    ("step7_minimize", "Energy minimization"),
    ("step8_equilibrate", "Equilibration"),
    ("step9_production", "Production MD"),
    ("step10_analysis", "Analysis"),
    ("step11_mmpbsa", "MM-PBSA & Report"),
]

LINKED_ANALYSIS_LABELS = [
    ("free_energy", "Free Energy"),
]

STATUS_LABELS = {
    "awaiting_topology": "Awaiting Topology",
    "queued": "Queued",
    "running": "Running",
    "failed": "Failed",
    "completed": "Completed",
    "finished_unknown": "Finished (unknown result)",
    "analysis_ready": "Analysis Ready",
}

AUTO_REFRESH_SECONDS = 60


@app.errorhandler(RequestEntityTooLarge)
def handle_request_entity_too_large(_error):
    max_mb = app.config["MAX_CONTENT_LENGTH"] // (1024 * 1024)
    flash(f"Upload is too large for the current limit ({max_mb} MB). Set MD_WEBAPP_MAX_UPLOAD_MB higher if needed.", "error")
    return redirect(url_for("index")), 413

FORCE_FIELD_OPTIONS = [
    ("charmm36-feb2026_cgenff-5.0", "CHARMM36 Feb 2026 (CGenFF 5.0)"),
    ("charmm36-jul2022", "CHARMM36 Jul 2022"),
    ("charmm36-mar2019", "CHARMM36 Mar 2019"),
    ("charmm36", "CHARMM36"),
    ("charmm27", "CHARMM27"),
    ("amber14sb", "AMBER14SB"),
    ("amber99sb-ildn", "AMBER99SB-ILDN"),
    ("amber99sb", "AMBER99SB"),
    ("amber99", "AMBER99"),
    ("amber96", "AMBER96"),
    ("amber94", "AMBER94"),
    ("amber03", "AMBER03"),
    ("amberGS", "AMBERGS"),
    ("oplsaa", "OPLS-AA/L"),
    ("gromos54a7", "GROMOS54A7"),
    ("gromos53a6", "GROMOS53A6"),
    ("gromos53a5", "GROMOS53A5"),
    ("gromos45a3", "GROMOS45A3"),
    ("gromos43a2", "GROMOS43A2"),
    ("gromos43a1", "GROMOS43A1"),
]

WATER_MODEL_OPTIONS = [
    ("tip3p", "TIP3P"),
    ("tip3p_original", "TIP3P Original"),
    ("spce", "SPC/E"),
    ("spc", "SPC"),
    ("tip4p", "TIP4P"),
    ("tip4pew", "TIP4P-Ew"),
    ("tip5p", "TIP5P"),
]


def now():
    return datetime.now()


def list_jobs():
    jobs = []
    for job_dir in sorted(JOBS_DIR.iterdir(), reverse=True):
        meta_path = job_dir / "job.json"
        if meta_path.exists():
            jobs.append(read_job(job_dir.name))
    return jobs


def parse_timestamp(value):
    if not value:
        return None
    try:
        return datetime.fromisoformat(value)
    except ValueError:
        return None


def format_runtime(delta_seconds):
    if delta_seconds is None:
        return "n/a"
    total_seconds = max(int(delta_seconds), 0)
    days, remainder = divmod(total_seconds, 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, seconds = divmod(remainder, 60)
    clock = f"{hours:02}:{minutes:02}:{seconds:02}"
    if days:
        return f"{days}d {clock}"
    return clock


def format_step_progress(current_step, total_steps):
    if current_step is None:
        return f"{total_steps}/{total_steps}"
    return f"{current_step}/{total_steps}"


def load_job_meta(job_dir):
    meta_path = job_dir / "job.json"
    if not meta_path.exists():
        raise FileNotFoundError(job_dir.name)
    with meta_path.open("r", encoding="utf-8-sig") as handle:
        return json.load(handle)


def save_job_meta(job_dir, meta):
    with (job_dir / "job.json").open("w", encoding="utf-8") as handle:
        json.dump(meta, handle, indent=2)


def infer_stage(job_dir, status):
    if status == "awaiting_topology":
        return {
            "completed_steps": 0,
            "current_step": 1,
            "current_stage": "Waiting for ligand topology upload",
            "progress_label": "Prep",
            "progress_percent": 5,
        }
    if status == "analysis_ready":
        return {
            "completed_steps": 0,
            "current_step": None,
            "current_stage": "Files saved for specialist analysis",
            "progress_label": "Analysis",
            "progress_percent": 100,
        }

    completed_steps = 0
    current_step = None
    current_stage = "Waiting to start"

    for index, (step_name, label) in enumerate(STEP_LABELS, start=1):
        if (job_dir / "state" / f"{step_name}.done").exists():
            completed_steps = index
            continue
        current_step = index
        if status == "failed":
            current_stage = f"Failed during {label}"
        elif status == "queued":
            current_stage = f"Waiting for {label}"
        elif status == "finished_unknown":
            current_stage = f"Stopped before {label}"
        else:
            current_stage = label
        break
    else:
        current_stage = "Analysis complete"

    return {
        "completed_steps": completed_steps,
        "current_step": current_step,
        "current_stage": current_stage,
        "progress_label": format_step_progress(current_step, len(STEP_LABELS)),
        "progress_percent": int((completed_steps / len(STEP_LABELS)) * 100),
    }


def collect_linked_analysis_status(job_dir, meta):
    analysis_plan = meta.get("analysis_plan", {})
    linked = []
    result_enabled = analysis_plan.get("result_analysis", "yes") == "yes"
    result_done = (job_dir / "state" / "step10_analysis.done").exists()
    linked.append({
        "key": "result_analysis",
        "label": "Result Analysis",
        "enabled": result_enabled,
        "done": result_done,
    })
    for key, label in LINKED_ANALYSIS_LABELS:
        enabled = analysis_plan.get(key, "no") == "yes"
        linked.append({
            "key": key,
            "label": label,
            "enabled": enabled,
            "done": (job_dir / "state" / f"{key}.done").exists(),
        })
    return linked


def infer_runtime(job_dir, meta, status):
    started_at = parse_timestamp(meta.get("created_at"))
    if started_at is None:
        return {}

    if status in {"running", "awaiting_topology"}:
        ended_at = now()
    elif status == "completed":
        ended_at = datetime.fromtimestamp((job_dir / "state" / "step10_analysis.done").stat().st_mtime)
    elif (job_dir / "run.failed").exists():
        ended_at = datetime.fromtimestamp((job_dir / "run.failed").stat().st_mtime)
    elif (job_dir / "run.log").exists():
        ended_at = datetime.fromtimestamp((job_dir / "run.log").stat().st_mtime)
    elif (job_dir / "prep.log").exists():
        ended_at = datetime.fromtimestamp((job_dir / "prep.log").stat().st_mtime)
    else:
        ended_at = started_at

    seconds = max((ended_at - started_at).total_seconds(), 0)
    return {
        "started_at_display": started_at.strftime("%Y-%m-%d %H:%M:%S"),
        "runtime_seconds": int(seconds),
        "runtime_display": format_runtime(seconds),
    }


def normalize_positive_int(value, field_name, required=True):
    text = (value or "").strip()
    if not text:
        if required:
            raise ValueError(f"{field_name} is required.")
        return ""
    try:
        number = int(text)
    except ValueError as exc:
        raise ValueError(f"{field_name} must be an integer.") from exc
    if number < 1:
        raise ValueError(f"{field_name} must be at least 1.")
    return str(number)


def normalize_nonnegative_int(value, field_name, required=True):
    text = (value or "").strip()
    if not text:
        if required:
            raise ValueError(f"{field_name} is required.")
        return ""
    try:
        number = int(text)
    except ValueError as exc:
        raise ValueError(f"{field_name} must be an integer.") from exc
    if number < 0:
        raise ValueError(f"{field_name} must be 0 or greater.")
    return str(number)


def infer_status(job_dir, meta):
    if (job_dir / "state" / "step10_analysis.done").exists():
        return "completed"
    if (job_dir / "run.failed").exists():
        return "failed"

    pid = meta.get("pid")
    if pid:
        try:
            os.kill(pid, 0)
        except OSError:
            return "finished_unknown"
        return "running"

    status_override = meta.get("status_override")
    if status_override:
        return status_override
    return "queued"


def tail_file(path, max_lines=30):
    if not path.exists():
        return ""
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    return "".join(lines[-max_lines:])


def resolve_log_tail(job_dir, status):
    if status == "awaiting_topology":
        return tail_file(job_dir / "prep.log")
    return tail_file(job_dir / "run.log")


def read_runtime_overrides(job_dir):
    override_path = job_dir / "state" / "auto_rectification.env"
    values = {}
    if not override_path.exists():
        return values
    with override_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line.startswith("export ") or "=" not in line:
                continue
            key, value = line[len("export "):].split("=", 1)
            values[key] = value.strip("'\"")
    return values


def collect_analysis_plots(job_id):
    analysis_dir = JOBS_DIR / job_id / "analysis"
    if not analysis_dir.exists():
        return []
    return [
        {"name": path.stem.replace("_", " ").title(), "filename": path.name}
        for path in sorted(analysis_dir.glob("*.png"))
    ]


def collect_linked_analysis_files(job_id):
    analysis_dir = JOBS_DIR / job_id / "analysis"
    groups = []
    for key, label in LINKED_ANALYSIS_LABELS:
        group_dir = analysis_dir / key
        if not group_dir.exists():
            continue
        files = [
            {"name": path.name, "path": str(path)}
            for path in sorted(group_dir.iterdir())
            if path.is_file()
        ]
        if files:
            groups.append({"key": key, "label": label, "files": files})
    return groups


def save_upload(file_storage, destination):
    destination.parent.mkdir(parents=True, exist_ok=True)
    file_storage.save(destination)


def localize_topology_includes(root_dir):
    local_includes = {"ligand.itp", "posre.itp", "posre_ligand.itp", "ligand.prm"}
    for top_path in Path(root_dir).rglob("topol.top"):
        text = top_path.read_text(encoding="utf-8", errors="replace")
        changed = False

        def replacement(match):
            nonlocal changed
            include_path = match.group(1)
            include_name = include_path.replace("\\", "/").rsplit("/", 1)[-1]
            if include_name in local_includes and (top_path.parent / include_name).exists() and include_path != include_name:
                changed = True
                return f'#include "{include_name}"'
            return match.group(0)

        updated = re.sub(r'#include\s+"([^"]+)"', replacement, text)
        if changed:
            top_path.write_text(updated, encoding="utf-8", newline="\n")


def write_posre_from_itp(itp_path, posre_path):
    in_atoms = False
    with open(itp_path, "r", encoding="utf-8", errors="replace") as src, open(posre_path, "w", encoding="utf-8", newline="\n") as dst:
        dst.write("[ position_restraints ]\n")
        dst.write("; ai funct fcx fcy fcz\n")
        for raw in src:
            stripped = raw.strip()
            if stripped.startswith("[ atoms ]"):
                in_atoms = True
                continue
            if in_atoms and stripped.startswith("[") and not stripped.startswith("[ atoms ]"):
                break
            if in_atoms and stripped and not stripped.startswith(";"):
                parts = stripped.split()
                if parts[0].isdigit():
                    dst.write(f"{int(parts[0]):6d} {1:6d} {1000:6d} {1000:6d} {1000:6d}\n")


def count_itp_atoms(itp_path):
    in_atoms = False
    count = 0
    with open(itp_path, "r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            stripped = raw.strip()
            if stripped.startswith("[ atoms ]"):
                in_atoms = True
                continue
            if in_atoms and stripped.startswith("[") and not stripped.startswith("[ atoms ]"):
                break
            if in_atoms and stripped and not stripped.startswith(";"):
                parts = stripped.split()
                if parts and parts[0].isdigit():
                    count += 1
    if count == 0:
        raise RuntimeError(f"Could not count atoms in {itp_path}")
    return count


def sanitize_ligand_posre(ligand_itp, posre_path):
    posre_path = Path(posre_path)
    if not posre_path.exists():
        write_posre_from_itp(ligand_itp, posre_path)
        return
    max_atom = count_itp_atoms(ligand_itp)
    kept = []
    kept_restraints = 0
    with posre_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            stripped = raw.strip()
            if not stripped or stripped.startswith(";") or stripped.startswith("["):
                kept.append(raw)
                continue
            parts = stripped.split()
            if parts and parts[0].isdigit() and int(parts[0]) <= max_atom:
                kept.append(raw)
                kept_restraints += 1
    if kept_restraints == 0:
        write_posre_from_itp(ligand_itp, posre_path)
    else:
        posre_path.write_text("".join(kept), encoding="utf-8", newline="\n")


def extract_moleculetype(itp_path):
    found = False
    with open(itp_path, "r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if raw.startswith("[ moleculetype ]"):
                found = True
                continue
            if found and (not line or line.startswith(";")):
                continue
            if found:
                return line.split()[0]
    raise RuntimeError(f"Could not determine [ moleculetype ] from {itp_path}")


def extract_topology_forcefield(top_path):
    text = Path(top_path).read_text(encoding="utf-8", errors="replace")
    match = re.search(r'#include\s+["<]([^">]+\.ff)/forcefield\.itp[">]', text)
    if not match:
        return ""
    return Path(match.group(1).replace("\\", "/")).name.removesuffix(".ff")


def python_has_analysis_packages(python_path):
    if not python_path:
        return False
    try:
        result = subprocess.run(
            [str(python_path), "-c", "import MDAnalysis, matplotlib, sklearn, scipy"],
            cwd=PIPELINE_DIR,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    except OSError:
        return False
    return result.returncode == 0


def resolve_python_bin(preferred=""):
    candidates = [
        preferred,
        os.environ.get("PYTHON_BIN", ""),
        sys.executable,
        shutil.which("python"),
        shutil.which("py"),
        r"C:\Python314\python.exe",
    ]
    seen = set()
    for candidate in candidates:
        if not candidate or candidate in seen:
            continue
        seen.add(candidate)
        if python_has_analysis_packages(candidate):
            return str(candidate)
    for candidate in candidates:
        if candidate and candidate not in {"python3"}:
            return str(candidate)
    return sys.executable

def normalize_gmx_bin(value):
    text = (value or "gmx").strip() or "gmx"
    parts = text.split()
    if len(parts) > 1 and parts[0] in {"gmx", "gmx_mpi"}:
        return parts[0]
    return text


def normalize_ntomp(value):
    text = (value or "auto").strip().lower()
    if text in {"", "auto", "0"}:
        return "auto"
    return normalize_positive_int(text, "CPU threads")


def launch_pipeline(job_dir, config_path, metadata):
    run_log = (job_dir / "run.log").open("w", encoding="utf-8")
    process = subprocess.Popen(
        [
            resolve_python_bin(),
            str(APP_DIR / "launch_job.py"),
            str(PIPELINE_DIR),
            str(config_path),
            str(job_dir / "run.failed"),
        ],
        cwd=PIPELINE_DIR,
        stdout=run_log,
        stderr=subprocess.STDOUT,
        start_new_session=True,
        env=os.environ.copy(),
    )
    run_log.close()
    metadata["pid"] = process.pid
    save_job_meta(job_dir, metadata)
    return process


def new_job_dir():
    job_id = datetime.now().strftime("%Y%m%d-%H%M%S") + "-" + uuid.uuid4().hex[:8]
    job_dir = JOBS_DIR / job_id
    (job_dir / "inputs").mkdir(parents=True, exist_ok=True)
    return job_id, job_dir


def command_exists(command_name):
    return bool(command_name and shutil.which(command_name))


def split_docking_complex(job_dir, docking_path):
    inputs_dir = job_dir / "inputs"
    prep_log = job_dir / "prep.log"
    protein_out = inputs_dir / "protein.pdb"
    ligand_pdb = inputs_dir / "ligand.pdb"
    ligand_mol2 = inputs_dir / "ligand.mol2"

    split_cmd = [
        resolve_python_bin(),
        str(PIPELINE_DIR / "helpers" / "split_docking_pdb.py"),
        str(docking_path),
        str(protein_out),
        str(ligand_pdb),
        "LIG",
        "yes",
    ]
    split_result = subprocess.run(
        split_cmd,
        cwd=PIPELINE_DIR,
        text=True,
        capture_output=True,
        check=False,
    )
    with prep_log.open("a", encoding="utf-8") as handle:
        handle.write(f"[{now().strftime('%Y-%m-%d %H:%M:%S')}] CMD: {' '.join(split_cmd)}\n")
        if split_result.stdout:
            handle.write(split_result.stdout)
        if split_result.stderr:
            handle.write(split_result.stderr)
    if split_result.returncode != 0:
        raise RuntimeError(split_result.stderr.strip() or split_result.stdout.strip() or "Docking complex split failed.")

    if not command_exists("obabel"):
        raise RuntimeError("Open Babel (`obabel`) is required to generate ligand MOL2 from the docked PDB.")

    obabel_cmd = ["obabel", str(ligand_pdb), "-O", str(ligand_mol2)]
    obabel_result = subprocess.run(
        obabel_cmd,
        cwd=PIPELINE_DIR,
        text=True,
        capture_output=True,
        check=False,
    )
    with prep_log.open("a", encoding="utf-8") as handle:
        handle.write(f"[{now().strftime('%Y-%m-%d %H:%M:%S')}] CMD: {' '.join(obabel_cmd)}\n")
        if obabel_result.stdout:
            handle.write(obabel_result.stdout)
        if obabel_result.stderr:
            handle.write(obabel_result.stderr)
    if obabel_result.returncode != 0:
        raise RuntimeError(obabel_result.stderr.strip() or obabel_result.stdout.strip() or "Open Babel ligand conversion failed.")

    return {
        "protein_pdb": str(protein_out),
        "ligand_pdb": str(ligand_pdb),
        "ligand_mol2": str(ligand_mol2),
        "docking_complex": str(docking_path),
    }


def build_override_config(job_dir, form):
    inputs_dir = job_dir / "inputs"
    config_path = job_dir / "web_config.json"

    docking_path = inputs_dir / "docking_complex.pdb"
    ligand_structure = inputs_dir / "ligand.mol2"
    ligand_itp = inputs_dir / "ligand.itp"
    ligand_gro = inputs_dir / "ligand.gro"
    ligand_posre = inputs_dir / "posre_ligand.itp"
    ligand_prm = inputs_dir / "ligand.prm"

    payload = {
        "PROJECT_NAME": form["project_name"],
        "INPUT_DIR": str(inputs_dir),
        "WORK_DIR": str(job_dir / "work"),
        "RESULTS_DIR": str(job_dir / "results"),
        "LOG_DIR": str(job_dir / "logs"),
        "STATE_DIR": str(job_dir / "state"),
        "ANALYSIS_DIR": str(job_dir / "analysis"),
        "TMP_DIR": str(job_dir / "tmp"),
        "FORCE_FIELD": form["force_field"],
        "WATER_MODEL": form["water_model"],
        "LIGAND_PREP_TOOL": form["ligand_prep_tool"],
        "AUTO_RECTIFY": form.get("auto_rectify", "yes"),
        "LIGAND_RESNAME": form["ligand_resname"],
        "LIGAND_NAME": form["ligand_name"],
        "LIGAND_NET_CHARGE": form["ligand_charge"],
        "NVT_STEPS": form["nvt_steps"],
        "NPT_STEPS": form["npt_steps"],
        "MD_TIME_NS": form["md_time_ns"],
        "NTOMP": form["ntomp"],
        "USE_GPU": form["use_gpu"],
        "STRIP_HETATM": form.get("strip_hetatm", "yes"),
        "MD_CENTER_NOJUMP": form.get("md_center_nojump", "yes"),
        "MD_CENTER_GROUP": form.get("md_center_group", "Protein_Ligand"),
        "MD_CENTER_OUTPUT_GROUP": form.get("md_center_output_group", "System"),
        "ANALYSIS_FIT_GROUP": form.get("analysis_fit_group", "Backbone"),
        "ANALYSIS_RMS_GROUP": form.get("analysis_rms_group", "Protein"),
        "ANALYSIS_LIGAND_GROUP": form.get("analysis_ligand_group", "Ligand"),
        "RUN_RESULT_ANALYSIS": form.get("run_result_analysis", "yes"),
        "RUN_FREE_ENERGY": form.get("run_free_energy", "no"),
        "RUN_MMBSA": form.get("run_mmpbsa", "yes"),
        "FREE_ENERGY_SELECTION": form.get("free_energy_selection", "backbone"),
        "FREE_ENERGY_TEMPERATURE_K": form.get("free_energy_temperature_k") or form.get("temperature_k", "300"),
        "FREE_ENERGY_COMPONENTS": form.get("free_energy_components", "3"),
        "FREE_ENERGY_BINS": form.get("free_energy_bins", "180"),
        "FREE_ENERGY_FRAME_STEP": form.get("free_energy_frame_step", "1"),
        "GMX_BIN": normalize_gmx_bin(form.get("gmx_bin", "gmx")),
        "PYTHON_BIN": resolve_python_bin(form.get("python_bin", "")),
        "DOCKING_COMPLEX_PDB": str(docking_path),
        "DOCKING_LIGAND_RESNAME": form["ligand_resname"],
        "LIGAND_INPUT": str(ligand_structure),
        "LIGAND_INPUT_FORMAT": "mol2",
    }

    if form["ligand_prep_tool"] == "prepared_topology":
        payload["PREPARED_LIGAND_ITP"] = str(ligand_itp)
        if ligand_gro.exists():
            payload["PREPARED_LIGAND_GRO"] = str(ligand_gro)
        if ligand_posre.exists():
            payload["PREPARED_LIGAND_POSRE"] = str(ligand_posre)
        if ligand_prm.exists():
            payload["PREPARED_LIGAND_PRM"] = str(ligand_prm)

    config_path.write_text(json.dumps(payload, indent=2), encoding="utf-8", newline="\n")
    return config_path


def save_optional_upload(file_storage, destination):
    if file_storage and file_storage.filename:
        save_upload(file_storage, destination)
        return str(destination)
    return ""


def build_analysis_plan(form):
    return {
        "result_analysis": form.get("run_result_analysis", "yes"),
        "free_energy": form.get("run_free_energy", "no"),
    }


def write_linked_analysis_summary(job_dir, analysis_key, title, source_dir):
    output_dir = job_dir / "analysis" / analysis_key
    output_dir.mkdir(parents=True, exist_ok=True)
    lines = [
        title,
        f"Generated: {datetime.now().isoformat(timespec='seconds')}",
        f"Source: {source_dir}",
        "",
        "Available files:",
    ]
    files = [path for path in sorted(Path(source_dir).iterdir()) if path.is_file()] if Path(source_dir).exists() else []
    if files:
        lines.extend(f"- {path.name} ({path.stat().st_size} bytes)" for path in files)
    else:
        lines.append("- No files were available in the source directory.")
    (output_dir / "analysis_summary.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    write_done_stamps(job_dir, [analysis_key])


def natural_sort_key(path):
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", Path(path).name)]


def safe_extract_zip(zip_path, target_dir):
    target_dir = Path(target_dir).resolve()
    with zipfile.ZipFile(zip_path) as archive:
        for member in archive.infolist():
            destination = (target_dir / member.filename).resolve()
            if target_dir != destination and target_dir not in destination.parents:
                raise ValueError(f"Unsafe ZIP path: {member.filename}")
        archive.extractall(target_dir)


def unpack_analysis_archives(analysis_dir):
    for zip_path in list(Path(analysis_dir).glob("*.zip")):
        target_dir = Path(analysis_dir) / zip_path.stem
        target_dir.mkdir(parents=True, exist_ok=True)
        safe_extract_zip(zip_path, target_dir)


def write_path_list(paths, destination):
    Path(destination).write_text(
        "\n".join(str(Path(path).resolve()) for path in paths) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def read_xvg_xy(path):
    x_vals = []
    y_vals = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip() or raw[0] in {"@", "#"}:
                continue
            parts = raw.split()
            if len(parts) < 2:
                continue
            try:
                x_vals.append(float(parts[0]))
                y_vals.append(float(parts[1]))
            except ValueError:
                continue
    return x_vals, y_vals



def run_optional_plot(input_xvg, title, output_png):
    if not Path(input_xvg).exists():
        return False
    try:
        subprocess.run(
            [resolve_python_bin(), str(PIPELINE_DIR / "helpers" / "plot_analysis.py"), str(input_xvg), title, str(output_png)],
            cwd=PIPELINE_DIR,
            check=True,
            text=True,
            capture_output=True,
        )
        return True
    except Exception:
        return False


def find_duplicate_names(paths):
    seen = set()
    duplicates = []
    for path in paths:
        name = Path(path).name.lower()
        if name in seen:
            duplicates.append(Path(path).name)
        seen.add(name)
    return duplicates


    tpr_files = sorted(Path(analysis_dir).rglob("*.tpr"), key=natural_sort_key)
    pullf_files = sorted([path for path in Path(analysis_dir).rglob("*.xvg") if "pullf" in path.name.lower()], key=natural_sort_key)
    pullx_files = sorted([path for path in Path(analysis_dir).rglob("*.xvg") if "pullx" in path.name.lower()], key=natural_sort_key)
    selected_pull_files = pullf_files if pullf_files else pullx_files
    selected_pull_flag = "-if" if pullf_files else "-ix"
    selected_pull_list_name = "pullf-files.dat" if pullf_files else "pullx-files.dat"

    with summary_path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write(f"TPR windows found: {len(tpr_files)}\n")
        handle.write(f"pullf files found: {len(pullf_files)}\n")
        handle.write(f"pullx files found: {len(pullx_files)}\n")
        handle.write("Selected pull input: " + selected_pull_list_name + "\n")

    duplicate_tprs = find_duplicate_names(tpr_files)
    duplicate_pulls = find_duplicate_names(selected_pull_files)
    if duplicate_tprs or duplicate_pulls:
        with summary_path.open("a", encoding="utf-8", newline="\n") as handle:
            handle.write("\nStatus: waiting for corrected files\n")
            handle.write("The PDF notes WHAM files must have unique names.\n")
            if duplicate_tprs:
                handle.write(f"Duplicate TPR names: {', '.join(duplicate_tprs)}\n")
            if duplicate_pulls:
                handle.write(f"Duplicate pull names: {', '.join(duplicate_pulls)}\n")
        return

    if len(tpr_files) != len(selected_pull_files):
        with summary_path.open("a", encoding="utf-8", newline="\n") as handle:
            handle.write("\nStatus: waiting for matching window files\n")
            handle.write("The number of .tpr files must match the number of selected pull files.\n")
            handle.write("List/order must correspond window-by-window as shown in the PDF.\n")
        return

    tpr_list = output_dir / "tpr-files.dat"
    pull_list = output_dir / selected_pull_list_name
    write_path_list(tpr_files, tpr_list)
    write_path_list(selected_pull_files, pull_list)

    profile_xvg = output_dir / "profile.xvg"
    histo_xvg = output_dir / "histo.xvg"
    gmx_bin = normalize_gmx_bin(form.get("gmx_bin", "gmx"))
    command = [
        gmx_bin,
        "wham",
        "-it",
        str(tpr_list),
        selected_pull_flag,
        str(pull_list),
        "-o",
        str(profile_xvg),
        "-hist",
        str(histo_xvg),
        "-unit",
        "kCal",
    ]

    with summary_path.open("a", encoding="utf-8", newline="\n") as handle:
        handle.write("\nGenerated input lists:\n")
        handle.write(f"- {tpr_list.name}\n")
        handle.write(f"- {pull_list.name}\n")
        handle.write("\nWHAM command:\n")
        handle.write(" ".join(command) + "\n")

    if not shutil.which(gmx_bin) and not Path(gmx_bin).exists():
        with summary_path.open("a", encoding="utf-8", newline="\n") as handle:
            handle.write("\nStatus: waiting for GROMACS\n")
            handle.write("GROMACS command was not found. Run the WHAM command manually after correcting the GROMACS path.\n")
        return

    result = subprocess.run(command, cwd=output_dir, text=True, capture_output=True, check=False)
    (output_dir / "wham.stdout.log").write_text(result.stdout or "", encoding="utf-8", newline="\n")
    (output_dir / "wham.stderr.log").write_text(result.stderr or "", encoding="utf-8", newline="\n")
    if result.returncode != 0:
        with summary_path.open("a", encoding="utf-8", newline="\n") as handle:
            handle.write(f"\nStatus: failed\nWHAM exited with code {result.returncode}. See wham.stdout.log and wham.stderr.log.\n")
        return


def read_job(job_id):
    job_dir = JOBS_DIR / job_id
    meta = load_job_meta(job_dir)
    overrides = read_runtime_overrides(job_dir)
    if "FORCE_FIELD" in overrides:
        meta["force_field"] = overrides["FORCE_FIELD"]
    if "WATER_MODEL" in overrides:
        meta["water_model"] = overrides["WATER_MODEL"]
    if "LIGAND_PREP_TOOL" in overrides:
        meta["ligand_prep_tool"] = overrides["LIGAND_PREP_TOOL"]
    meta["auto_rectification"] = "Applied" if overrides else meta.get("auto_rectification", "Enabled")
    meta["id"] = job_id
    meta["status"] = infer_status(job_dir, meta)
    meta["status_label"] = STATUS_LABELS.get(meta["status"], meta["status"].replace("_", " ").title())
    meta.update(infer_stage(job_dir, meta["status"]))
    meta.update(infer_runtime(job_dir, meta, meta["status"]))
    meta["run_log_tail"] = resolve_log_tail(job_dir, meta["status"])
    meta["analysis_plots"] = collect_analysis_plots(job_id)
    meta["linked_analyses"] = collect_linked_analysis_status(job_dir, meta)
    meta["linked_analysis_files"] = collect_linked_analysis_files(job_id)
    meta["auto_refresh_seconds"] = AUTO_REFRESH_SECONDS if meta["status"] == "running" else None
    meta["has_report"] = (job_dir / "results" / "md_report.pdf").exists()
    return meta


@app.route("/", methods=["GET"])
def index():
    return render_template(
        "index.html",
        jobs=list_jobs(),
        pipeline_dir=str(PIPELINE_DIR),
        force_field_options=FORCE_FIELD_OPTIONS,
        water_model_options=WATER_MODEL_OPTIONS,
    )


@app.route("/jobs/<job_id>/analysis/<filename>")
def analysis_file(job_id, filename):
    path = JOBS_DIR / job_id / "analysis" / filename
    if not path.exists() or path.suffix.lower() != ".png":
        abort(404)
    return send_file(path)


@app.route("/jobs/prepare", methods=["POST"])
def prepare_job():
    docking_complex = request.files.get("docking_file")
    project_name = request.form.get("project_name", "protein_ligand_job").strip() or "protein_ligand_job"

    if not docking_complex or not docking_complex.filename:
        flash("Docked PDB upload is required.", "error")
        return redirect(url_for("index"))

    job_id, job_dir = new_job_dir()
    inputs_dir = job_dir / "inputs"
    docking_path = inputs_dir / "docking_complex.pdb"
    save_upload(docking_complex, docking_path)

    try:
        generated_files = split_docking_complex(job_dir, docking_path)
    except RuntimeError as exc:
        flash(str(exc), "error")
        return redirect(url_for("index"))

    metadata = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "project_name": project_name,
        "input_mode": "docking_complex",
        "status_override": "awaiting_topology",
        "generated_files": generated_files,
    }
    save_job_meta(job_dir, metadata)
    flash(f"Docked complex prepared for job {job_id}. Upload the ligand topology to continue.", "success")
    return redirect(url_for("topology_upload", job_id=job_id))


def collect_run_form(default_project):
    return {
        "project_name": request.form.get("project_name", default_project).strip() or default_project,
        "force_field": request.form.get("force_field", "amber99sb-ildn"),
        "water_model": request.form.get("water_model", "tip3p"),
        "ligand_resname": request.form.get("ligand_resname", "LIG").upper()[:5],
        "ligand_name": request.form.get("ligand_name", "ligand"),
        "ligand_charge": request.form.get("ligand_charge", "0"),
        "ligand_prep_tool": request.form.get("ligand_prep_tool", "acpype"),
        "nvt_steps": request.form.get("nvt_steps", ""),
        "npt_steps": request.form.get("npt_steps", ""),
        "md_time_ns": request.form.get("md_time_ns", "100"),
        "ntomp": request.form.get("ntomp", "auto"),
        "use_gpu": request.form.get("use_gpu", "auto"),
        "strip_hetatm": "yes" if request.form.get("strip_hetatm", "on") else "no",
        "md_center_nojump": "yes" if request.form.get("md_center_nojump", "on") else "no",
        "gmx_bin": request.form.get("gmx_bin", r"gmx"),
        "python_bin": request.form.get("python_bin", ""),
        "auto_rectify": "yes" if request.form.get("auto_rectify", "on") else "no",
        "run_result_analysis": "yes" if request.form.get("run_result_analysis") else "no",
        "run_free_energy": "yes" if request.form.get("run_free_energy") else "no",
        "analysis_fit_group": request.form.get("analysis_fit_group", "Backbone"),
        "analysis_rms_group": request.form.get("analysis_rms_group", "Protein"),
        "analysis_ligand_group": request.form.get("analysis_ligand_group", "Ligand"),
        "md_center_group": request.form.get("md_center_group", "Protein_Ligand"),
        "md_center_output_group": request.form.get("md_center_output_group", "System"),
        "free_energy_selection": request.form.get("free_energy_selection", "backbone"),
        "free_energy_temperature_k": request.form.get("free_energy_temperature_k", "300"),
        "free_energy_components": request.form.get("free_energy_components", "3"),
        "free_energy_bins": request.form.get("free_energy_bins", "180"),
        "free_energy_frame_step": request.form.get("free_energy_frame_step", "1"),
    }


@app.route("/jobs/full-auto", methods=["POST"])
def full_auto_job():
    docking_complex = request.files.get("docking_file")
    if not docking_complex or not docking_complex.filename:
        flash("Docked PDB upload is required for full automation.", "error")
        return redirect(url_for("index"))

    form = collect_run_form("protein_ligand_auto")
    if form["ligand_prep_tool"] not in {"acpype", "cgenff"}:
        flash("Unsupported ligand topology generation option.", "error")
        return redirect(url_for("index"))

    try:
        form["md_time_ns"] = normalize_positive_int(form["md_time_ns"], "Production time")
        form["ntomp"] = normalize_ntomp(form["ntomp"])
    except ValueError as exc:
        flash(str(exc), "error")
        return redirect(url_for("index"))
    for field, label in (("nvt_steps", "NVT nsteps"), ("npt_steps", "NPT nsteps")):
        try:
            form[field] = normalize_positive_int(form[field], label, required=False)
        except ValueError as exc:
            flash(str(exc), "error")
            return redirect(url_for("index"))

    job_id, job_dir = new_job_dir()
    docking_path = job_dir / "inputs" / "docking_complex.pdb"
    save_upload(docking_complex, docking_path)
    try:
        generated_files = split_docking_complex(job_dir, docking_path)
    except RuntimeError as exc:
        flash(str(exc), "error")
        return redirect(url_for("index"))

    config_path = build_override_config(job_dir, form)
    metadata = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "project_name": form["project_name"],
        "force_field": form["force_field"],
        "water_model": form["water_model"],
        "ligand_prep_tool": form["ligand_prep_tool"],
        "auto_rectification": "Enabled" if form["auto_rectify"] == "yes" else "Disabled",
        "analysis_plan": build_analysis_plan(form),
        "input_mode": "full_auto_docking_complex",
        "generated_files": generated_files,
    }
    launch_pipeline(job_dir, config_path, metadata)
    flash(f"Full automatic MD job {job_id} started.", "success")
    return redirect(url_for("job_detail", job_id=job_id))


def write_done_stamps(job_dir, step_names):
    state_dir = job_dir / "state"
    state_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    for step_name in step_names:
        (state_dir / f"{step_name}.done").write_text(stamp + "\n", encoding="utf-8")


@app.route("/jobs/old-method", methods=["POST"])
def old_method_job():
    required_uploads = {
        "protein_gro": "protein_processed.gro",
        "protein_top": "topol.top",
        "ligand_gro": "ligand.gro",
        "ligand_itp": "ligand.itp",
    }
    missing = [name for name in required_uploads if not request.files.get(name) or not request.files[name].filename]
    if missing:
        flash("Old method requires protein GRO, protein topology, ligand GRO, and ligand ITP.", "error")
        return redirect(url_for("index"))

    form = collect_run_form("protein_ligand_old_method")
    form["ligand_prep_tool"] = "prepared_topology"
    job_id, job_dir = new_job_dir()
    work_dir = job_dir / "work"
    protein_dir = work_dir / "01_protein"
    ligand_dir = work_dir / "02_ligand"
    complex_dir = work_dir / "03_complex"
    for path in (protein_dir, ligand_dir, complex_dir):
        path.mkdir(parents=True, exist_ok=True)

    save_upload(request.files["protein_gro"], protein_dir / "protein_processed.gro")
    save_upload(request.files["protein_top"], protein_dir / "topol.top")
    detected_forcefield = extract_topology_forcefield(protein_dir / "topol.top")
    if detected_forcefield and detected_forcefield != form["force_field"]:
        form["force_field"] = detected_forcefield
        flash(f"Using force field {detected_forcefield} referenced by the uploaded protein topology.", "success")
    save_upload(request.files["ligand_gro"], ligand_dir / "ligand.gro")
    save_upload(request.files["ligand_itp"], ligand_dir / "ligand.itp")
    ligand_posre = request.files.get("ligand_posre")
    if ligand_posre and ligand_posre.filename:
        save_upload(ligand_posre, ligand_dir / "posre_ligand.itp")
    else:
        write_posre_from_itp(ligand_dir / "ligand.itp", ligand_dir / "posre_ligand.itp")
    sanitize_ligand_posre(ligand_dir / "ligand.itp", ligand_dir / "posre_ligand.itp")
    ligand_prm = request.files.get("ligand_prm")
    if ligand_prm and ligand_prm.filename:
        save_upload(ligand_prm, ligand_dir / "ligand.prm")
    protein_posre = request.files.get("protein_posre")
    if protein_posre and protein_posre.filename:
        save_upload(protein_posre, protein_dir / "posre.itp")
    else:
        (protein_dir / "posre.itp").write_text("; optional protein position restraints were not uploaded\n", encoding="utf-8")

    shutil.copy2(protein_dir / "topol.top", complex_dir / "topol.top")
    shutil.copy2(protein_dir / "posre.itp", complex_dir / "posre.itp")
    shutil.copy2(ligand_dir / "ligand.itp", complex_dir / "ligand.itp")
    shutil.copy2(ligand_dir / "posre_ligand.itp", complex_dir / "posre_ligand.itp")
    copy_prm = ligand_dir / "ligand.prm"
    if copy_prm.exists():
        shutil.copy2(copy_prm, complex_dir / "ligand.prm")
    subprocess.run(
        [
            resolve_python_bin(),
            str(PIPELINE_DIR / "helpers" / "merge_gro.py"),
            str(protein_dir / "protein_processed.gro"),
            str(ligand_dir / "ligand.gro"),
            str(complex_dir / "complex.gro"),
        ],
        check=True,
        cwd=PIPELINE_DIR,
    )
    molecule_name = extract_moleculetype(ligand_dir / "ligand.itp")
    (ligand_dir / "ligand_metadata.json").write_text(json.dumps({"LIGAND_MOLECULE_NAME": molecule_name}), encoding="utf-8")
    subprocess.run(
        [
            resolve_python_bin(),
            str(PIPELINE_DIR / "helpers" / "assemble_topology.py"),
            str(complex_dir / "topol.top"),
            "ligand.itp",
            "posre_ligand.itp",
            molecule_name,
            "1",
            str(PIPELINE_DIR),
            form["force_field"],
        ],
        check=True,
        cwd=complex_dir,
    )
    write_done_stamps(job_dir, ["step1_protein_prep", "step2_ligand_prep", "step3_assemble_complex"])
    config_path = build_override_config(job_dir, form)
    metadata = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "project_name": form["project_name"],
        "force_field": form["force_field"],
        "water_model": form["water_model"],
        "ligand_prep_tool": "prepared_topology",
        "auto_rectification": "Enabled" if form["auto_rectify"] == "yes" else "Disabled",
        "analysis_plan": build_analysis_plan(form),
        "input_mode": "old_method_gro_topology",
        "generated_files": {},
    }
    launch_pipeline(job_dir, config_path, metadata)
    flash(f"Old-method continuation job {job_id} started from assembled complex.", "success")
    return redirect(url_for("job_detail", job_id=job_id))


STAGE_START_STEPS = {
    "equilibration": "step8_equilibrate",
    "nvt": "step8_equilibrate",
    "npt": "step8_equilibrate",
    "md": "step9_production",
    "analysis": "step10_analysis",
}


@app.route("/jobs/resume-stage", methods=["POST"])
def resume_stage_job():
    stage_type = request.form.get("stage_type", "analysis")
    previous_job_id = request.form.get("previous_job_id", "").strip()
    restart_zip = request.files.get("restart_zip")
    if stage_type not in STAGE_START_STEPS:
        flash("Unsupported continuation stage.", "error")
        return redirect(url_for("index"))
    if not previous_job_id and (not restart_zip or not restart_zip.filename):
        flash("Select a previous run or upload a restart ZIP for this stage.", "error")
        return redirect(url_for("index"))

    job_id, job_dir = new_job_dir()
    if previous_job_id:
        previous_dir = JOBS_DIR / previous_job_id
        if not previous_dir.exists():
            flash("Previous run was not found.", "error")
            return redirect(url_for("index"))
        shutil.rmtree(job_dir)
        shutil.copytree(previous_dir, job_dir, ignore=shutil.ignore_patterns("run.log", "run.failed", "job.json"))
        localize_topology_includes(job_dir)
    else:
        zip_path = job_dir / "inputs" / "restart.zip"
        save_upload(restart_zip, zip_path)
        with zipfile.ZipFile(zip_path) as archive:
            archive.extractall(job_dir)
        localize_topology_includes(job_dir)

    meta = load_job_meta(job_dir) if (job_dir / "job.json").exists() else {}
    form = collect_run_form(meta.get("project_name", f"{stage_type}_restart"))
    if not request.form.get("force_field") and meta.get("force_field"):
        form["force_field"] = meta["force_field"]
    if not request.form.get("water_model") and meta.get("water_model"):
        form["water_model"] = meta["water_model"]
    if not request.form.get("ligand_prep_tool") and meta.get("ligand_prep_tool"):
        form["ligand_prep_tool"] = meta["ligand_prep_tool"]
    config_path = build_override_config(job_dir, form)

    start_step = STAGE_START_STEPS[stage_type]
    remove = False
    for step_name, _label in STEP_LABELS:
        if step_name == start_step:
            remove = True
        if remove:
            stamp = job_dir / "state" / f"{step_name}.done"
            if stamp.exists():
                stamp.unlink()
    if stage_type in {"md", "analysis"}:
        for linked_key, _label in LINKED_ANALYSIS_LABELS:
            stamp = job_dir / "state" / f"{linked_key}.done"
            if stamp.exists():
                stamp.unlink()

    metadata = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "project_name": form["project_name"],
        "force_field": form["force_field"],
        "water_model": form["water_model"],
        "ligand_prep_tool": meta.get("ligand_prep_tool", form["ligand_prep_tool"]),
        "auto_rectification": "Enabled" if form["auto_rectify"] == "yes" else "Disabled",
        "analysis_plan": build_analysis_plan(form),
        "input_mode": f"resume_{stage_type}",
        "generated_files": meta.get("generated_files", {}),
    }
    launch_pipeline(job_dir, config_path, metadata)
    flash(f"{stage_type.upper()} continuation job {job_id} started.", "success")
    return redirect(url_for("job_detail", job_id=job_id))



@app.route("/jobs/special-analysis", methods=["POST"])
def special_analysis_job():
    analysis_type = request.form.get("analysis_type", "free_energy")
    project_name = request.form.get("project_name", f"{analysis_type}_analysis").strip() or f"{analysis_type}_analysis"
    analysis_plan = {
        "result_analysis": "yes" if analysis_type == "result_analysis" else "no",
        "free_energy": "yes" if request.form.get("run_free_energy") or analysis_type == "free_energy_3d" else "no",
    }
    job_id, job_dir = new_job_dir()
    analysis_dir = job_dir / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    for upload in request.files.getlist("analysis_files"):
        if upload and upload.filename:
            save_upload(upload, analysis_dir / Path(upload.filename).name)
    if analysis_type == "result_analysis":
        write_done_stamps(job_dir, ["step10_analysis"])
        run_requested_linked_analyses(job_dir, analysis_plan, analysis_dir)
    elif analysis_type == "free_energy_3d":
        run_requested_linked_analyses(job_dir, analysis_plan, analysis_dir)
    metadata = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "project_name": project_name,
        "input_mode": analysis_type,
        "status_override": "analysis_ready",
        "analysis_plan": analysis_plan,
        "generated_files": {"analysis_dir": str(analysis_dir)},
    }
    save_job_meta(job_dir, metadata)
    flash(f"{analysis_type.replace('_', ' ').title()} files saved in job {job_id}.", "success")
    return redirect(url_for("job_detail", job_id=job_id))


@app.route("/jobs/mmpbsa", methods=["POST"])
def mmpbsa_job():
    """Launch MM-PBSA binding free energy calculation."""
    project_name = request.form.get("project_name", "mmpbsa_analysis").strip() or "mmpbsa_analysis"
    previous_job_id = request.form.get("previous_job_id", "")

    if not previous_job_id:
        flash("Please select a previous run with production MD data.", "error")
        return redirect(url_for("index"))

    prev_job_dir = JOBS_DIR / previous_job_id
    if not prev_job_dir.exists():
        flash(f"Job {previous_job_id} not found.", "error")
        return redirect(url_for("index"))

    # Verify production files exist
    tpr_file = prev_job_dir / "results" / "production" / "md.tpr"
    xtc_file = prev_job_dir / "results" / "production" / "md.xtc"
    if not tpr_file.exists() or not xtc_file.exists():
        flash("Production MD files (md.tpr, md.xtc) not found. Run production MD first.", "error")
        return redirect(url_for("index"))

    job_id, job_dir = new_job_dir()
    results_dir = job_dir / "results" / "production"
    results_dir.mkdir(parents=True, exist_ok=True)
    analysis_dir = job_dir / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    # Symlink production files
    (results_dir / "md.tpr").symlink_to(tpr_file.resolve())
    (results_dir / "md.xtc").symlink_to(xtc_file.resolve())

    # Copy topology files (.itp, .top, .ff) from previous job's production dir
    # Required for MM-PBSA grompp energy group setup
    prev_prod_dir = prev_job_dir / "results" / "production"
    if prev_prod_dir.exists():
        for item in prev_prod_dir.iterdir():
            if item.is_file() and item.suffix in (".itp", ".top", ".prm"):
                shutil.copy2(item, results_dir / item.name)
            elif item.is_dir() and item.suffix == ".ff":
                shutil.copytree(item, results_dir / item.name, dirs_exist_ok=True)

    # Copy web_config.json from previous job
    prev_config = prev_job_dir / "web_config.json"
    if prev_config.exists():
        shutil.copy2(prev_config, job_dir / "web_config.json")

    # Collect MM-PBSA parameters
    frame_step = request.form.get("frame_step", "10")
    ionic_strength = request.form.get("ionic_strength", "0.15")
    surface_tension = request.form.get("surface_tension", "0.0072")
    decompose = "yes" if request.form.get("decompose") else "no"
    use_apbs = "yes" if request.form.get("use_apbs", "on") else "no"

    # Build override config
    config_data = {}
    if prev_config.exists():
        try:
            config_data = json.loads(prev_config.read_text())
        except Exception:
            pass

    config_data.update({
        "RUN_MMBSA": "yes",
        "MMBSA_FRAME_STEP": frame_step,
        "MMBSA_IONIC_STRENGTH": ionic_strength,
        "MMBSA_SURFACE_TENSION": surface_tension,
        "MMBSA_DECOMPOSE": decompose,
        "MMBSA_USE_APBS": use_apbs,
        "ANALYSIS_DIR": str(analysis_dir),
    })

    config_path = job_dir / "web_config.json"
    config_path.write_text(json.dumps(config_data, indent=2), encoding="utf-8")

    # Copy index.ndx and analysis files from previous job for MM-PBSA and report generation
    prev_analysis_dir = prev_job_dir / "analysis"
    if prev_analysis_dir.exists():
        # Copy index.ndx (required for receptor/ligand selection in MM-PBSA)
        prev_index = prev_analysis_dir / "index.ndx"
        if prev_index.exists():
            shutil.copy2(prev_index, analysis_dir / "index.ndx")
        # Copy all analysis plots and files for PDF report
        for item in prev_analysis_dir.iterdir():
            if item.is_file():
                shutil.copy2(item, analysis_dir / item.name)
            elif item.is_dir():
                shutil.copytree(item, analysis_dir / item.name, dirs_exist_ok=True)

    # Launch MM-PBSA step
    write_done_stamps(job_dir, ["step1_protein_prep", "step2_ligand_prep",
                                "step3_assemble_complex", "step4_define_box",
                                "step5_solvate", "step6_add_ions",
                                "step7_minimize", "step8_equilibrate",
                                "step9_production", "step10_analysis"])

    metadata = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "project_name": project_name,
        "input_mode": "mmpbsa",
        "status_override": "running",
        "force_field": config_data.get("FORCE_FIELD", "amber99sb-ildn"),
        "water_model": config_data.get("WATER_MODEL", "tip3p"),
        "previous_job": previous_job_id,
    }
    save_job_meta(job_dir, metadata)

    # Run MM-PBSA as background process
    run_log = (job_dir / "run.log").open("w", encoding="utf-8")
    subprocess.Popen(
        [
            resolve_python_bin(),
            str(APP_DIR / "launch_job.py"),
            str(PIPELINE_DIR),
            str(config_path),
            str(job_dir / "run.failed"),
        ],
        cwd=PIPELINE_DIR,
        stdout=run_log,
        stderr=subprocess.STDOUT,
        start_new_session=True,
        env=os.environ.copy(),
    )

    flash(f"MM-PBSA calculation started for job {job_id} (based on {previous_job_id}).", "success")
    return redirect(url_for("job_detail", job_id=job_id))


@app.route("/jobs/<job_id>/download-report")
def download_report(job_id):
    """Download the generated PDF report."""
    job_dir = JOBS_DIR / job_id
    report_path = job_dir / "results" / "md_report.pdf"
    if not report_path.exists():
        flash("PDF report not yet generated.", "error")
        return redirect(url_for("job_detail", job_id=job_id))
    return send_file(str(report_path), as_attachment=True,
                     download_name=f"{job_id}_md_report.pdf")


@app.route("/jobs/<job_id>/topology", methods=["GET", "POST"])
def topology_upload(job_id):
    job_dir = JOBS_DIR / job_id
    try:
        meta = load_job_meta(job_dir)
    except FileNotFoundError:
        abort(404)

    if request.method == "GET":
        return render_template(
            "topology.html",
            job=read_job(job_id),
            force_field_options=FORCE_FIELD_OPTIONS,
            water_model_options=WATER_MODEL_OPTIONS,
        )

    ligand_prep_tool = request.form.get("ligand_prep_tool", "prepared_topology")
    ligand_itp = request.files.get("ligand_itp")
    ligand_gro = request.files.get("ligand_gro")
    ligand_posre = request.files.get("ligand_posre")
    ligand_prm = request.files.get("ligand_prm")

    if ligand_prep_tool == "prepared_topology" and (not ligand_itp or not ligand_itp.filename):
        flash("Ligand topology (.itp) is required to continue.", "error")
        return redirect(url_for("topology_upload", job_id=job_id))

    try:
        md_time_ns = normalize_positive_int(request.form.get("md_time_ns"), "Production time (ns)")
        ntomp = normalize_ntomp(request.form.get("ntomp", "auto"))
        nvt_steps = normalize_positive_int(request.form.get("nvt_steps"), "NVT nsteps", required=False)
        npt_steps = normalize_positive_int(request.form.get("npt_steps"), "NPT nsteps", required=False)
    except ValueError as exc:
        flash(str(exc), "error")
        return redirect(url_for("topology_upload", job_id=job_id))

    inputs_dir = job_dir / "inputs"
    if ligand_prep_tool == "prepared_topology":
        save_upload(ligand_itp, inputs_dir / "ligand.itp")
        if ligand_gro and ligand_gro.filename:
            save_upload(ligand_gro, inputs_dir / "ligand.gro")
        if ligand_posre and ligand_posre.filename:
            save_upload(ligand_posre, inputs_dir / "posre_ligand.itp")
        if ligand_prm and ligand_prm.filename:
            save_upload(ligand_prm, inputs_dir / "ligand.prm")
    elif ligand_prep_tool not in {"acpype", "cgenff", "prepared_topology"}:
        flash("Unsupported ligand topology generation option.", "error")
        return redirect(url_for("topology_upload", job_id=job_id))

    form = {
        "project_name": request.form.get("project_name", meta.get("project_name", f"md_job_{job_id}")).strip() or meta.get("project_name", f"md_job_{job_id}"),
        "force_field": request.form.get("force_field", "amber99sb-ildn"),
        "water_model": request.form.get("water_model", "tip3p"),
        "ligand_resname": request.form.get("ligand_resname", "LIG").upper()[:5],
        "ligand_name": request.form.get("ligand_name", "ligand"),
        "ligand_charge": request.form.get("ligand_charge", "0"),
        "ligand_prep_tool": ligand_prep_tool,
        "nvt_steps": nvt_steps,
        "npt_steps": npt_steps,
        "md_time_ns": md_time_ns,
        "ntomp": ntomp,
        "use_gpu": request.form.get("use_gpu", "auto"),
        "strip_hetatm": "yes" if request.form.get("strip_hetatm") else "no",
        "md_center_nojump": "yes" if request.form.get("md_center_nojump") else "no",
        "gmx_bin": request.form.get("gmx_bin", "gmx"),
        "auto_rectify": "yes" if request.form.get("auto_rectify") else "no",
        "run_result_analysis": "yes" if request.form.get("run_result_analysis") else "no",
        "run_free_energy": "yes" if request.form.get("run_free_energy") else "no",
        "free_energy_selection": request.form.get("free_energy_selection", "backbone"),
        "free_energy_temperature_k": request.form.get("free_energy_temperature_k", "300"),
        "free_energy_components": request.form.get("free_energy_components", "3"),
        "free_energy_bins": request.form.get("free_energy_bins", "180"),
        "free_energy_frame_step": request.form.get("free_energy_frame_step", "1"),
    }
    config_path = build_override_config(job_dir, form)

    run_log = (job_dir / "run.log").open("w", encoding="utf-8")
    process = subprocess.Popen(
        [
            resolve_python_bin(),
            str(APP_DIR / "launch_job.py"),
            str(PIPELINE_DIR),
            str(config_path),
            str(job_dir / "run.failed"),
        ],
        cwd=PIPELINE_DIR,
        stdout=run_log,
        stderr=subprocess.STDOUT,
        start_new_session=True,
        env=os.environ.copy(),
    )
    run_log.close()

    metadata = {
        "created_at": meta.get("created_at", datetime.now().isoformat(timespec="seconds")),
        "project_name": form["project_name"],
        "force_field": form["force_field"],
        "water_model": form["water_model"],
        "ligand_prep_tool": form["ligand_prep_tool"],
        "auto_rectification": "Enabled" if form["auto_rectify"] == "yes" else "Disabled",
        "analysis_plan": build_analysis_plan(form),
        "input_mode": "docking_complex",
        "pid": process.pid,
        "generated_files": meta.get("generated_files", {}),
    }
    save_job_meta(job_dir, metadata)

    flash(f"Job {job_id} started.", "success")
    return redirect(url_for("job_detail", job_id=job_id))


@app.route("/jobs/<job_id>", methods=["GET"])
def job_detail(job_id):
    try:
        job = read_job(job_id)
    except FileNotFoundError:
        abort(404)
    return render_template("job.html", job=job)


@app.route("/jobs/<job_id>/stop", methods=["POST"])
def stop_job(job_id):
    try:
        job = read_job(job_id)
    except FileNotFoundError:
        abort(404)
    pid = job.get("pid")
    if pid and job["status"] == "running":
        os.killpg(pid, signal.SIGTERM)
        flash(f"Stop signal sent to job {job_id}.", "success")
    else:
        flash("Job is not currently running.", "error")
    return redirect(url_for("job_detail", job_id=job_id))


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=True)
