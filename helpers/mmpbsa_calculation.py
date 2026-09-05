#!/usr/bin/env python3
"""MM-PBSA binding free-energy calculation (GROMACS MM + APBS PB + SASA).

Drop-in replacement for helpers/mmpbsa_calculation.py.

Fixes vs. the previous helper
-----------------------------
* Do not write PDB of the solvated system for GROMACS. Write GRO instead so
  MDAnalysis does not spam stderr (altLocs / chainIDs / occupancies) and so
  grompp sees native coordinates with a box.
* PB/SA use the dry complex (receptor + ligand) only — never waters/ions.
* PQR is written from TPR charges + Bondi radii (ligands work; PDB2PQR is
  optional).
* APBS input is valid mg-auto with calcenergy; polar energy is parsed from
  APBS output (no broken DX interpolation).
* gmx energy XVG legends are parsed as `@ s0 legend "..."`.
* Receptor-Ligand *and* Ligand-Receptor terms, including 1-4.
* Automation-friendly: MDAnalysis warnings suppressed, MMPBSA_STATUS line,
  non-zero exit only on hard failure.

Outputs (in --output-dir)
    mmpbsa_results.csv
    mmpbsa_summary.csv
    mmpbsa_summary.json
    mmpbsa_decomposition.png / mmpbsa_timeseries.png  (if matplotlib)
    decomposition.csv  (if --decompose)
"""
from __future__ import annotations

import argparse
import csv
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)
os.environ.setdefault("MDA_SILENT", "1")

try:
    import MDAnalysis as mda
    logging.getLogger("MDAnalysis").setLevel(logging.ERROR)
except ImportError as exc:
    raise SystemExit(f"MDAnalysis is required: {exc}") from exc


NM_TO_ANG = 10.0
KJ_TO_KCAL = 1.0 / 4.184
COULOMB_KCAL = 332.0636

log = logging.getLogger("mmpbsa")


# ── process helpers ──────────────────────────────────────────────────────────
def _run(cmd, cwd=None, input_text=None, check=True, log_name="cmd"):
    result = subprocess.run(
        cmd,
        cwd=cwd,
        capture_output=True,
        text=True,
        input=input_text,
        check=False,
        env={**os.environ, "GMX_MAXBACKUP": "-1", "GMX_NO_TERM": "1"},
    )
    if check and result.returncode != 0:
        err = (result.stderr or result.stdout or "")[-2500:]
        raise RuntimeError(
            f"{log_name} failed ({result.returncode}): {' '.join(str(c) for c in cmd)}\n{err}"
        )
    return result


def _gmx_version(gmx_bin):
    r = _run([gmx_bin, "--version"], check=False, log_name="gmx --version")
    m = re.search(r"VERSION\s+(\S+)", (r.stdout or "") + (r.stderr or ""))
    return m.group(1) if m else "unknown"


def _which(name):
    return shutil.which(name) is not None


# ── index / topology discovery ───────────────────────────────────────────────
def _parse_ndx(path):
    groups = {}
    current = None
    with open(path, "r") as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("[") and line.endswith("]"):
                current = line[1:-1].strip()
                groups[current] = []
            elif current is not None:
                groups[current].extend(int(x) - 1 for x in line.split() if x.lstrip("-").isdigit())
    return groups


def _write_ndx(path, groups):
    lines = []
    for name, indices in groups.items():
        lines.append(f"[ {name} ]")
        row = []
        for i, idx in enumerate(indices):
            row.append(str(int(idx) + 1))
            if len(row) == 15:
                lines.append(" ".join(row))
                row = []
        if row:
            lines.append(" ".join(row))
    path.write_text("\n".join(lines) + "\n")


def _find_file(start: Path, names, extra_globs=None, depth=3):
    dirs = [start]
    cur = start
    for _ in range(depth):
        cur = cur.parent
        dirs.append(cur)
    for d in dirs:
        for n in names:
            p = d / n
            if p.exists():
                return p
        if extra_globs:
            for g in extra_globs:
                hits = sorted(d.glob(g))
                if hits:
                    return hits[0]
    return None


def _stage_topology(tpr_or_top: Path, temp_dir: Path):
    """Copy .top/.itp/.ff next to a working directory so grompp can #include them."""
    src_dir = tpr_or_top.parent
    for pattern in ("*.itp", "*.top", "*.rtp", "*.dat"):
        for f in src_dir.glob(pattern):
            dest = temp_dir / f.name
            if not dest.exists():
                shutil.copy2(f, dest)
    for ff in src_dir.glob("*.ff"):
        if ff.is_dir() and not (temp_dir / ff.name).exists():
            shutil.copytree(ff, temp_dir / ff.name)

    if tpr_or_top.suffix == ".top":
        dest = temp_dir / tpr_or_top.name
        if tpr_or_top.resolve() != dest.resolve():
            shutil.copy2(tpr_or_top, dest)
        return dest

    found = _find_file(src_dir, ["topol.top", "system.top", "complex.top"], extra_globs=["*.top"])
    if found:
        dest = temp_dir / found.name
        shutil.copy2(found, dest)
        # also copy includes from the .top's directory if different
        if found.parent.resolve() != src_dir.resolve():
            for f in found.parent.glob("*.itp"):
                shutil.copy2(f, temp_dir / f.name)
            for ff in found.parent.glob("*.ff"):
                if ff.is_dir() and not (temp_dir / ff.name).exists():
                    shutil.copytree(ff, temp_dir / ff.name)
        return dest
    return None


# ── GROMACS MM ───────────────────────────────────────────────────────────────
def _build_energy_mdp(mdp_template_path, out_mdp_path):
    content = mdp_template_path.read_text() if mdp_template_path else ""
    drop_prefixes = (
        "energygrps",
        "continuation",
        "nsteps",
        "nstxout",
        "nstvout",
        "nstfout",
        "nstenergy",
        "nstlog",
        "nstxout-compressed",
        "nstcalcenergy",
        # Anything below references index-group NAMES from the original
        # simulation (e.g. "Protein", "Non-Protein", "SOL"). Our rerun only
        # defines System/Receptor/Ligand/Complex in index.ndx, so leaving
        # these in makes grompp fail with
        # "Group X referenced in the .mdp file was not found".
        "tc-grps",
        "tcoupl",
        "tau-t",
        "ref-t",
        "nsttcouple",
        "nh-chain-length",
        "pcoupl",
        "pcoupltype",
        "tau-p",
        "ref-p",
        "compressibility",
        "refcoord-scaling",
        "acc-grps",
        "accelerate",
        "freezegrps",
        "freezedim",
        "energygrp-excl",
        "energygrp-flags",
        "deform",
        "wall-atomtype",
        "wall-density",
        "qmmm-grps",
        "gen-vel",
        "gen-temp",
        "gen-seed",
    )
    filtered = []
    skip_cont = False
    for line in content.splitlines():
        stripped = line.strip()
        if skip_cont:
            skip_cont = stripped.endswith("\\")
            continue
        # GROMACS treats "_" and "-" as interchangeable in mdp parameter
        # names (tau_t == tau-t, gen_vel == gen-vel, ...). Normalize before
        # matching drop_prefixes so an underscore-spelled template doesn't
        # slip past the filter and collide with the hyphen-spelled
        # replacements added below ("doubly defined" grompp error).
        normalized = stripped.replace("_", "-")
        if any(normalized.startswith(p) for p in drop_prefixes):
            skip_cont = stripped.endswith("\\")
            continue
        filtered.append(line)

    additions = [
        "nsteps = 0",
        "nstenergy = 1",
        "nstcalcenergy = 1",
        "nstlog = 1",
        "nstxout = 0",
        "nstvout = 0",
        "nstfout = 0",
        "nstxout-compressed = 0",
        "energygrps = Receptor Ligand",
        "continuation = yes",
        # Group-free equivalents: nsteps=0 rerun does no integration, so
        # coupling algorithms are irrelevant, but grompp still parses these
        # keywords and needs *some* valid value.
        "tcoupl = no",
        "pcoupl = no",
        "gen-vel = no",
        "tc-grps = System",
        "tau-t = 0.1",
        "ref-t = 300",
    ]
    body = "\n".join(filtered).rstrip()
    out_mdp_path.write_text(body + ("\n" if body else "") + "\n".join(additions) + "\n")


def _discover_energy_terms(gmx_bin, edr_path, tpr_path):
    result = _run(
        [gmx_bin, "energy", "-f", str(edr_path), "-s", str(tpr_path)],
        input_text="\n",
        check=False,
        log_name="gmx energy (discover)",
    )
    text = (result.stdout or "") + "\n" + (result.stderr or "")
    available = {}
    for line in text.splitlines():
        for idx, name in re.findall(r"(\d+)\s+([A-Za-z][\w().\-:/]*)", line):
            n = int(idx)
            if n <= 0:
                continue
            if name in ("kJ/mol", "K", "bar", "nm", "nm^3", "kg/m^3", "End", "Select"):
                continue
            available[name] = n
    return available


def _pick_pair_terms(available):
    """Return {kind: term_name} for receptor-ligand pair interactions."""
    picked = {}
    kinds = [
        ("coul", ("Coul-SR", "Coulomb-SR", "Coulomb-(SR)", "Coul-SR:")),
        ("lj", ("LJ-SR", "LJ-(SR)", "LJ-SR:")),
        ("coul14", ("Coul-14", "Coulomb-14")),
        ("lj14", ("LJ-14",)),
    ]
    for kind, needles in kinds:
        for name in available:
            low = name
            if not any(n in low for n in needles):
                continue
            has_rec = "Receptor" in name
            has_lig = "Ligand" in name
            if has_rec and has_lig:
                picked[kind] = name
                break
    return picked


def _parse_xvg(path):
    legend = []
    rows = []
    with open(path, "r") as fh:
        for line in fh:
            if line.startswith("@"):
                m = re.search(r"@\s*s(\d+)\s+legend\s+\"(.+)\"", line)
                if m:
                    legend.append((int(m.group(1)), m.group(2)))
                continue
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    return legend, rows


def _extract_interaction_energy(gmx_bin, edr_path, tpr_path):
    """Return (vdw_kj, elec_kj) from Receptor-Ligand pair terms."""
    edr_path = Path(edr_path)
    tpr_path = Path(tpr_path)
    if not edr_path.exists():
        raise RuntimeError(f"EDR not written: {edr_path}")

    available = _discover_energy_terms(gmx_bin, edr_path, tpr_path)
    picked = _pick_pair_terms(available)
    if "coul" not in picked and "lj" not in picked:
        raise RuntimeError(
            "No Receptor-Ligand pair terms in EDR. "
            f"Available: {list(available.keys())[:40]}"
        )

    order = [k for k in ("coul", "lj", "coul14", "lj14") if k in picked]
    indices = [str(available[picked[k]]) for k in order] + ["0"]
    energy_xvg = edr_path.parent / "energy.xvg"
    _run(
        [gmx_bin, "energy", "-f", str(edr_path), "-s", str(tpr_path), "-o", str(energy_xvg)],
        cwd=str(edr_path.parent),
        input_text="\n".join(indices) + "\n",
        check=False,
        log_name="gmx energy (extract)",
    )
    if not energy_xvg.exists():
        raise RuntimeError("gmx energy did not write energy.xvg")

    legend, rows = _parse_xvg(energy_xvg)
    if not rows:
        raise RuntimeError("energy.xvg contains no data rows")

    # Map legend series index -> kind. Series s0 is the first energy column
    # (column 1 of the numeric row; column 0 is time).
    name_by_series = {s: name for s, name in legend}
    col_for_kind = {}
    for series, name in name_by_series.items():
        for kind, term in picked.items():
            if term in name or name in term:
                col_for_kind[kind] = series + 1  # +1 skip time
    if not col_for_kind:
        # Fall back to selection order: first energy col = order[0], ...
        for i, kind in enumerate(order):
            col_for_kind[kind] = i + 1

    def mean_col(kind):
        c = col_for_kind.get(kind)
        if c is None:
            return 0.0
        vals = [r[c] for r in rows if len(r) > c]
        return float(np.mean(vals)) if vals else 0.0

    vdw = mean_col("lj") + mean_col("lj14")
    elec = mean_col("coul") + mean_col("coul14")
    return vdw, elec


# ── PQR / APBS ───────────────────────────────────────────────────────────────
_BONDI = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
    "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98, "NA": 2.27, "MG": 1.73,
    "ZN": 1.39, "FE": 1.80, "K": 2.75, "CA": 2.31, "MN": 1.73, "CU": 1.40,
}


def _element_of(atom):
    el = ""
    try:
        el = (atom.element or "").upper()
    except Exception:
        el = ""
    if el and el in _BONDI:
        return el
    name = re.sub(r"[0-9].*$", "", (getattr(atom, "name", "") or "").upper())
    typ = re.sub(r"[0-9].*$", "", (getattr(atom, "type", "") or "").upper())
    for cand in (name, typ):
        for n in sorted(_BONDI, key=len, reverse=True):
            if cand == n or cand.startswith(n):
                return n
    return "C"


def _vdw_radius(atom):
    return _BONDI.get(_element_of(atom), 1.70)


def _atom_charge(atom):
    try:
        return float(atom.charge)
    except Exception:
        return 0.0


def _format_pqr(serial, atom, chain, x, y, z):
    name = (atom.name or "X")[:4]
    if len(name) < 4:
        aname = f" {name:<3s}"
    else:
        aname = name
    resname = (atom.resname or "UNK")[:4]
    resid = int(getattr(atom, "resid", 1)) % 10000
    q = _atom_charge(atom)
    r = _vdw_radius(atom)
    return (
        f"ATOM  {serial:5d} {aname} {resname:>3s} {chain} {resid:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f} {q:7.4f} {r:6.4f}"
    )


def _write_pqr(ag, path, chain="A"):
    lines = []
    for i, atom in enumerate(ag, start=1):
        x, y, z = atom.position
        lines.append(_format_pqr(i, atom, chain, x, y, z))
    lines.append("TER")
    lines.append("END")
    Path(path).write_text("\n".join(lines) + "\n")


def _prepare_pdb_attrs(ag, chain):
    n = len(ag)
    if n == 0:
        return
    try:
        ag.chainIDs = np.array([chain] * n)
    except Exception:
        pass


def _write_pdb_quiet(ag, path, chain="A"):
    _prepare_pdb_attrs(ag, chain)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ag.write(str(path))


def _nice_dime(n):
    cands = sorted({2**a * 3**b * 5**c + 1 for a in range(8) for b in range(6) for c in range(4)})
    n = max(33, int(n))
    for d in cands:
        if d >= n:
            return d
    return 129


def _grid_from_coords(coords_ang, spacing=0.5, fadd=10.0):
    mins = coords_ang.min(axis=0)
    maxs = coords_ang.max(axis=0)
    center = 0.5 * (mins + maxs)
    flen = np.maximum((maxs - mins) + 2.0 * fadd, 20.0)
    clen = flen * 1.7
    dime = [_nice_dime(flen[i] / spacing + 1) for i in range(3)]
    return dime, clen, flen, center


def _write_apbs_in(pqr_path, in_path, dime, cglen, fglen, center, temp, ionic, pdie, sdie):
    dx, dy, dz = dime
    cx, cy, cz = cglen
    fx, fy, fz = fglen
    ox, oy, oz = center
    pqr_name = Path(pqr_path).name
    body = textwrap.dedent(f"""\
        read
            mol pqr {pqr_name}
        end
        elec name polar
            mg-auto
            dime {dx} {dy} {dz}
            cglen {cx:.3f} {cy:.3f} {cz:.3f}
            fglen {fx:.3f} {fy:.3f} {fz:.3f}
            cgcent {ox:.3f} {oy:.3f} {oz:.3f}
            fgcent {ox:.3f} {oy:.3f} {oz:.3f}
            mol 1
            lpbe
            bcfl sdh
            pdie {pdie}
            sdie {sdie}
            srfm smol
            chgm spl2
            srad 1.4
            swin 0.3
            sdens 10.0
            temp {temp}
            ion charge 1 conc {ionic} radius 2.0
            ion charge -1 conc {ionic} radius 2.0
            calcenergy total
            calcforce no
        end
        print elecEnergy polar end
        quit
    """)
    Path(in_path).write_text(body)


def _parse_apbs_energy(stdout, stderr):
    text = (stdout or "") + "\n" + (stderr or "")
    patterns = [
        r"Global net ELEC energy\s*=\s*([-+0-9.eE]+)\s*(kJ/mol|kJ\/mol|kcal/mol)?",
        r"Total electrostatic energy\s*=\s*([-+0-9.eE]+)\s*(kJ/mol|kJ\/mol|kcal/mol)?",
        r"elecEnergy\s+([-+0-9.eE]+)\s*(kJ/mol|kJ\/mol|kcal/mol)?",
    ]
    for pat in patterns:
        matches = re.findall(pat, text)
        if matches:
            val, unit = matches[-1]
            energy = float(val)
            unit = (unit or "kJ/mol").lower()
            if "kcal" in unit:
                return energy
            return energy * KJ_TO_KCAL
    raise RuntimeError("Could not parse APBS electrostatic energy from output.")


def _run_apbs_energy(pqr_path, work_dir, dime, cglen, fglen, center, temp, ionic, pdie, sdie):
    pqr_path = Path(pqr_path)
    in_path = Path(work_dir) / (pqr_path.stem + ".in")
    _write_apbs_in(pqr_path, in_path, dime, cglen, fglen, center, temp, ionic, pdie, sdie)
    result = _run(["apbs", str(in_path.name)], cwd=str(work_dir), check=True, log_name="apbs")
    return _parse_apbs_energy(result.stdout, result.stderr)


# ── SASA ─────────────────────────────────────────────────────────────────────
def _fibonacci_sphere(n_points):
    indices = np.arange(n_points, dtype=float) + 0.5
    phi = np.arccos(1.0 - 2.0 * indices / n_points)
    theta = np.pi * (1.0 + 5**0.5) * indices
    x = np.cos(theta) * np.sin(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(phi)
    return np.column_stack([x, y, z])


def _sasa_shrake(positions, radii, probe=1.4, n_points=92):
    if len(positions) == 0:
        return 0.0
    r = np.asarray(radii, dtype=float) + probe
    try:
        from scipy.spatial import cKDTree
        tree = cKDTree(positions)
        query = lambda i, rad: tree.query_ball_point(positions[i], rad)
    except Exception:
        def query(i, rad):
            d = np.linalg.norm(positions - positions[i], axis=1)
            return np.where(d <= rad)[0].tolist()

    sphere = _fibonacci_sphere(n_points)
    max_r = float(r.max())
    total = 0.0
    for i in range(len(positions)):
        nbs = query(i, r[i] + max_r)
        pts = positions[i] + sphere * r[i]
        buried = np.zeros(n_points, dtype=bool)
        for j in nbs:
            if j == i:
                continue
            d = np.linalg.norm(pts - positions[j], axis=1)
            buried |= d < r[j]
        total += (1.0 - buried.mean()) * 4.0 * np.pi * r[i] ** 2
    return float(total)


def _compute_sa(ag, probe=1.4):
    if len(ag) == 0:
        return 0.0
    radii = np.array([_vdw_radius(atom) for atom in ag], dtype=float)
    return _sasa_shrake(ag.positions, radii, probe=probe)


# ── MM fallback from charges ─────────────────────────────────────────────────
def _mm_coulomb_fallback(rec_ag, lig_ag, box_ang):
    rec_pos = rec_ag.positions
    lig_pos = lig_ag.positions
    rec_q = np.array([_atom_charge(a) for a in rec_ag])
    lig_q = np.array([_atom_charge(a) for a in lig_ag])
    box = np.asarray(box_ang[:3], dtype=float)
    box[box == 0] = 1e9
    elec = 0.0
    vdw = 0.0
    rec_r = np.array([_vdw_radius(a) for a in rec_ag])
    lig_r = np.array([_vdw_radius(a) for a in lig_ag])
    for i in range(len(lig_pos)):
        diffs = rec_pos - lig_pos[i]
        diffs -= box * np.round(diffs / box)
        dist = np.sqrt((diffs ** 2).sum(axis=1))
        dist = np.clip(dist, 0.12, None)
        elec += COULOMB_KCAL * np.sum(rec_q * lig_q[i] / dist)
        sigma = 0.5 * (rec_r + lig_r[i])
        eps = 0.08
        sr6 = (sigma / dist) ** 6
        vdw += np.sum(4.0 * eps * (sr6 ** 2 - sr6))
    return float(vdw), float(elec)


# ── plotting ─────────────────────────────────────────────────────────────────
def _generate_plots(results, summary, out_dir):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        log.info("matplotlib not available; skipping plots")
        return

    components = ["vdw", "elec", "polar_solv", "nonpolar_solv", "delta_G"]
    labels = ["vdW", "Electrostatics", "Polar solv.", "Nonpolar solv.", "Total ΔG"]
    colors = ["#9aa890", "#8a9aa8", "#c17b6a", "#b8a58a", "#c5d0b8"]

    fig, ax = plt.subplots(figsize=(10, 5.5))
    means = [summary[c]["mean"] for c in components]
    stds = [summary[c]["std"] for c in components]
    bars = ax.bar(labels, means, yerr=stds, capsize=4, color=colors, edgecolor="#0e100e", linewidth=0.4)
    ax.set_ylabel("Energy (kcal/mol)")
    ax.set_title("MM-PBSA binding free energy")
    ax.axhline(0, color="#2c322c", linewidth=0.8)
    ax.grid(axis="y", alpha=0.25)
    for bar, val in zip(bars, means):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.4 * np.sign(bar.get_height() or 1),
            f"{val:.2f}",
            ha="center",
            va="bottom" if val >= 0 else "top",
            fontsize=9,
        )
    fig.tight_layout()
    fig.savefig(str(out_dir / "mmpbsa_decomposition.png"), dpi=180)
    plt.close(fig)

    times = [r["time_ps"] for r in results]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    series = [
        ("vdw", "van der Waals", "#9aa890", axes[0, 0]),
        ("elec", "Electrostatics", "#8a9aa8", axes[0, 1]),
        ("polar_solv", "Polar solvation", "#c17b6a", axes[1, 0]),
        ("delta_G", "Total ΔG", "#c5d0b8", axes[1, 1]),
    ]
    for key, title, color, ax in series:
        vals = [r[key] for r in results]
        ax.plot(times, vals, color=color, linewidth=1.0)
        ax.set_title(title)
        ax.set_ylabel("kcal/mol")
        ax.grid(alpha=0.25)
        mean_val = float(np.mean(vals))
        ax.axhline(mean_val, color=color, linestyle="--", alpha=0.7, label=f"mean {mean_val:.2f}")
        ax.legend(fontsize=8)
    axes[1, 0].set_xlabel("Time (ps)")
    axes[1, 1].set_xlabel("Time (ps)")
    fig.suptitle("MM-PBSA energy components")
    fig.tight_layout()
    fig.savefig(str(out_dir / "mmpbsa_timeseries.png"), dpi=180)
    plt.close(fig)


def _per_residue_decomposition(u, rec_ag, lig_ag, frame_indices, out_dir):
    residues = rec_ag.residues
    decomp = []
    max_frames = min(20, len(frame_indices))
    for res in residues:
        vdw_list, elec_list = [], []
        res_atoms = res.atoms
        for fidx in frame_indices[:max_frames]:
            u.trajectory[fidx]
            box = u.dimensions[:3]
            vdw, elec = _mm_coulomb_fallback(res_atoms, lig_ag, box)
            vdw_list.append(vdw)
            elec_list.append(elec)
        vdw_m = float(np.mean(vdw_list)) if vdw_list else 0.0
        elec_m = float(np.mean(elec_list)) if elec_list else 0.0
        decomp.append({
            "residue": f"{res.resname}{res.resid}",
            "resid": int(res.resid),
            "resname": str(res.resname),
            "vdw_mean": vdw_m,
            "elec_mean": elec_m,
            "total_mean": vdw_m + elec_m,
        })
    csv_path = out_dir / "decomposition.csv"
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["residue", "resid", "resname", "vdw_mean", "elec_mean", "total_mean"]
        )
        writer.writeheader()
        writer.writerows(decomp)
    log.info("Per-residue decomposition: %s", csv_path)
    ordered = sorted(decomp, key=lambda x: x["total_mean"])
    (out_dir / "decomposition_summary.json").write_text(json.dumps({
        "top_binding": ordered[:10],
        "top_repulsive": ordered[-10:][::-1],
    }, indent=2))


# ── main calculation ─────────────────────────────────────────────────────────
def calculate_mmpbsa(args):
    gmx_bin = args.gmx_bin
    log.info("GROMACS: %s (%s)", gmx_bin, _gmx_version(gmx_bin))
    log.info("Trajectory: %s", args.trajectory)
    log.info("Topology:   %s", args.topology)

    tpr_path = Path(args.topology)
    traj_path = Path(args.trajectory)

    ndx_path_in = args.index_file
    if not ndx_path_in:
        found = _find_file(tpr_path.parent, ["index.ndx", "md.ndx", "prod.ndx"], extra_globs=["*.ndx"])
        if found:
            ndx_path_in = str(found)
            log.info("Auto-detected index file: %s", found)

    ndx_groups = {}
    if ndx_path_in and Path(ndx_path_in).exists():
        ndx_groups = _parse_ndx(ndx_path_in)
        log.info("Index groups: %s", ", ".join(ndx_groups.keys()))

    u = mda.Universe(str(tpr_path), str(traj_path))
    n_frames = len(u.trajectory)
    log.info("Frames in trajectory: %d", n_frames)

    def resolve_selection(sel, role):
        if not sel:
            return "protein" if role == "receptor" else "resname LIG UNK MOL INH LIGAND"
        if sel.lower().startswith("group "):
            name = sel.split(None, 1)[1].strip()
            # case-insensitive group lookup
            match = next((g for g in ndx_groups if g.lower() == name.lower()), None)
            if match is None:
                log.warning("Group '%s' not in index. Available: %s", name, list(ndx_groups.keys()))
                if role == "ligand":
                    return "not protein and not resname SOL HOH WAT NA CL K NA+ CL- ION NAW CLW"
                return "protein"
            indices = ndx_groups[match]
            if not indices:
                raise SystemExit(f"Index group '{match}' is empty.")
            return "index " + " ".join(str(i) for i in indices)
        return sel

    sel_receptor = resolve_selection(args.receptor_selection or "protein", "receptor")
    sel_ligand = resolve_selection(args.ligand_selection or "resname LIG", "ligand")

    receptor_atoms = u.select_atoms(sel_receptor)
    ligand_atoms = u.select_atoms(sel_ligand)
    log.info("Receptor: %d atoms | Ligand: %d atoms", len(receptor_atoms), len(ligand_atoms))
    if len(receptor_atoms) == 0 or len(ligand_atoms) == 0:
        raise SystemExit(
            "One or more selections returned zero atoms. "
            f"receptor='{sel_receptor}' ligand='{sel_ligand}'. "
            "Pass -n index.ndx and --ligand-selection 'group Ligand'."
        )

    complex_atoms = receptor_atoms + ligand_atoms

    frame_step = max(1, args.frame_step)
    frame_indices = list(range(0, n_frames, frame_step))
    if args.max_frames and len(frame_indices) > args.max_frames:
        frame_indices = frame_indices[: args.max_frames]
    log.info("Processing %d frames (step=%d)", len(frame_indices), frame_step)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(tempfile.mkdtemp(prefix="mmpbsa_"))
    keep_temp = args.keep_temp

    energy_tpr = None
    top_for_grompp = None
    energy_mdp = None
    ndx_path = temp_dir / "index.ndx"

    try:
        _write_ndx(ndx_path, {
            "System": list(range(u.atoms.n_atoms)),
            "Receptor": [a.index for a in receptor_atoms],
            "Ligand": [a.index for a in ligand_atoms],
            "Complex": [a.index for a in complex_atoms],
        })

        base_mdp = Path(args.mdp_file) if args.mdp_file else None
        if not base_mdp or not base_mdp.exists():
            found_mdp = _find_file(tpr_path.parent, ["md.mdp", "prod.mdp", "mdout.mdp"], extra_globs=["*.mdp"])
            base_mdp = found_mdp
        if not base_mdp or not Path(base_mdp).exists():
            base_mdp = temp_dir / "minimal.mdp"
            base_mdp.write_text(textwrap.dedent("""\
                integrator = md
                dt = 0.002
                cutoff-scheme = Verlet
                nstlist = 10
                rlist = 1.2
                coulombtype = PME
                rcoulomb = 1.2
                vdwtype = Cut-off
                rvdw = 1.2
                pbc = xyz
                constraints = h-bonds
                constraint_algorithm = lincs
            """))
        energy_mdp = temp_dir / "energy.mdp"
        _build_energy_mdp(Path(base_mdp), energy_mdp)

        top_for_grompp = _stage_topology(tpr_path, temp_dir)
        if top_for_grompp:
            log.info("Using topology for grompp: %s", top_for_grompp.name)
        else:
            log.warning(
                "No .top file found next to the TPR. MM will use a Coulomb fallback "
                "from topology charges (approximate LJ)."
            )

        # Frame-0 GRO of the FULL system (atom order matches TPR/topology).
        u.trajectory[frame_indices[0]]
        gro0 = temp_dir / "frame0.gro"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u.atoms.write(str(gro0))

        if top_for_grompp is not None:
            energy_tpr = temp_dir / "energy.tpr"
            gp = _run(
                [
                    gmx_bin, "grompp",
                    "-f", str(energy_mdp),
                    "-c", str(gro0),
                    "-p", str(top_for_grompp),
                    "-n", str(ndx_path),
                    "-o", str(energy_tpr),
                    "-maxwarn", "10",
                ],
                cwd=str(temp_dir),
                check=False,
                log_name="grompp",
            )
            if not energy_tpr.exists():
                log.warning("grompp failed; MM fallback will be used.\n%s", (gp.stderr or "")[-1500:])
                energy_tpr = None
            else:
                log.info("Energy-groups TPR ready")

        use_apbs = bool(args.use_apbs) and _which("apbs")
        if args.use_apbs and not _which("apbs"):
            log.warning("apbs not on PATH; polar solvation will be 0. Install APBS or pass --no-apbs.")
            use_apbs = False

        results = []
        for i, fidx in enumerate(frame_indices):
            u.trajectory[fidx]
            frame = int(u.trajectory.frame)
            time_ps = float(u.trajectory.time)
            log.info("Frame %d/%d  (frame=%d  t=%.1f ps)", i + 1, len(frame_indices), frame, time_ps)

            frame_work = temp_dir / f"frame_{fidx:06d}"
            frame_work.mkdir(exist_ok=True)

            # Full-system GRO for MM rerun
            gro_full = frame_work / "full.gro"
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                u.atoms.write(str(gro_full))

            vdw = elec = 0.0
            mm_ok = False
            if energy_tpr is not None:
                rerun_edr = frame_work / "frame.edr"
                rerun_log = frame_work / "frame.log"
                md = _run(
                    [
                        gmx_bin, "mdrun",
                        "-s", str(energy_tpr),
                        "-rerun", str(gro_full),
                        "-e", str(rerun_edr),
                        "-g", str(rerun_log),
                        "-o", str(frame_work / "frame.trr"),
                        "-cpo", str(frame_work / "state.cpt"),
                        "-ntomp", "1",
                        "-nb", "cpu",
                    ],
                    cwd=str(frame_work),
                    check=False,
                    log_name="mdrun -rerun",
                )
                try:
                    vdw_kj, elec_kj = _extract_interaction_energy(gmx_bin, rerun_edr, energy_tpr)
                    vdw = vdw_kj * KJ_TO_KCAL
                    elec = elec_kj * KJ_TO_KCAL
                    mm_ok = True
                except Exception as exc:
                    log.warning("MM extraction failed: %s", exc)
                    if md.stderr:
                        log.debug("%s", md.stderr[-800:])
            if not mm_ok:
                vdw, elec = _mm_coulomb_fallback(
                    receptor_atoms, ligand_atoms, u.dimensions[:3]
                )
                log.info("    MM fallback (charges): vdW=%+.2f  Elec=%+.2f kcal/mol", vdw, elec)
            else:
                log.info("    MM: vdW=%+.2f  Elec=%+.2f kcal/mol", vdw, elec)

            polar_solv = 0.0
            if use_apbs:
                try:
                    pqr_c = frame_work / "complex.pqr"
                    pqr_r = frame_work / "receptor.pqr"
                    pqr_l = frame_work / "ligand.pqr"
                    _write_pqr(complex_atoms, pqr_c, chain="C")
                    _write_pqr(receptor_atoms, pqr_r, chain="R")
                    _write_pqr(ligand_atoms, pqr_l, chain="L")
                    dime, cglen, fglen, center = _grid_from_coords(complex_atoms.positions)
                    kwargs = dict(
                        dime=dime, cglen=cglen, fglen=fglen, center=center,
                        temp=args.temperature, ionic=args.ionic_strength,
                        pdie=args.pdie, sdie=args.sdie,
                    )
                    e_c = _run_apbs_energy(pqr_c, frame_work, **kwargs)
                    e_r = _run_apbs_energy(pqr_r, frame_work, **kwargs)
                    e_l = _run_apbs_energy(pqr_l, frame_work, **kwargs)
                    polar_solv = e_c - e_r - e_l
                    log.info("    PB: %+.2f kcal/mol", polar_solv)
                except Exception as exc:
                    log.warning("    PB failed: %s", exc)

            sa_c = _compute_sa(complex_atoms)
            sa_r = _compute_sa(receptor_atoms)
            sa_l = _compute_sa(ligand_atoms)
            dsa = sa_c - sa_r - sa_l
            nonpolar_solv = args.surface_tension * dsa
            log.info("    SA: ΔSA=%.1f Å²  ΔG_SA=%+.2f kcal/mol", dsa, nonpolar_solv)

            delta_g = vdw + elec + polar_solv + nonpolar_solv
            log.info("    ΔG_total = %+.2f kcal/mol", delta_g)

            results.append({
                "frame": frame,
                "time_ps": time_ps,
                "vdw": vdw,
                "elec": elec,
                "polar_solv": polar_solv,
                "nonpolar_solv": nonpolar_solv,
                "delta_G": delta_g,
                "SA_complex": sa_c,
                "SA_receptor": sa_r,
                "SA_ligand": sa_l,
            })

            # keep disk light
            for leftover in frame_work.glob("*.trr"):
                leftover.unlink(missing_ok=True)
            for leftover in frame_work.glob("*.cpt"):
                leftover.unlink(missing_ok=True)

        if not results:
            raise SystemExit("No frames were processed.")

        csv_path = out_dir / "mmpbsa_results.csv"
        with open(csv_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(results[0].keys()))
            writer.writeheader()
            writer.writerows(results)
        log.info("Per-frame results: %s", csv_path)

        components = ["vdw", "elec", "polar_solv", "nonpolar_solv", "delta_G"]
        summary = {}
        for comp in components:
            values = np.array([r[comp] for r in results], dtype=float)
            summary[comp] = {
                "mean": float(np.mean(values)),
                "std": float(np.std(values)),
                "median": float(np.median(values)),
                "min": float(np.min(values)),
                "max": float(np.max(values)),
                "sem": float(np.std(values) / np.sqrt(len(values))),
            }

        with open(out_dir / "mmpbsa_summary.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["Component", "Mean (kcal/mol)", "Std Dev", "Median", "Min", "Max", "SEM"])
            for comp in components:
                s = summary[comp]
                w.writerow([comp, f"{s['mean']:.4f}", f"{s['std']:.4f}", f"{s['median']:.4f}",
                            f"{s['min']:.4f}", f"{s['max']:.4f}", f"{s['sem']:.4f}"])

        summary_json = {
            "n_frames": len(results),
            "frame_step": frame_step,
            "temperature": args.temperature,
            "pH": args.pH,
            "ionic_strength": args.ionic_strength,
            "surface_tension": args.surface_tension,
            "energy_method": "GROMACS energy groups" if energy_tpr is not None else "charge-based Coulomb fallback",
            "pb_method": "APBS Poisson-Boltzmann (calcenergy)" if use_apbs else "skipped",
            "selections": {"receptor": sel_receptor[:200], "ligand": sel_ligand[:200]},
            "components": summary,
        }
        (out_dir / "mmpbsa_summary.json").write_text(json.dumps(summary_json, indent=2))
        log.info("JSON summary: %s", out_dir / "mmpbsa_summary.json")

        if args.decompose:
            log.info("Running per-residue decomposition...")
            _per_residue_decomposition(u, receptor_atoms, ligand_atoms, frame_indices, out_dir)

        _generate_plots(results, summary, out_dir)
        log.info("MM-PBSA calculation complete.")
        print("MMPBSA_STATUS=success")
        print(f"MMPBSA_DELTA_G={summary['delta_G']['mean']:.4f}")
        return summary_json

    finally:
        if keep_temp:
            log.info("Temp dir kept: %s", temp_dir)
        else:
            shutil.rmtree(temp_dir, ignore_errors=True)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="MM-PBSA binding free energy calculation")
    p.add_argument("-t", "--trajectory", required=True, help="GROMACS trajectory (.xtc)")
    p.add_argument("-s", "--topology", required=True, help="GROMACS topology (.tpr or .gro)")
    p.add_argument("-o", "--output-dir", required=True, help="Output directory")
    p.add_argument("-n", "--index-file", default=None, help="GROMACS index file (.ndx)")
    p.add_argument("-f", "--mdp-file", default=None, help="MDP template file")
    p.add_argument("--gmx-bin", default="gmx", help="GROMACS binary")
    p.add_argument("--receptor-selection", default=None, help="MDAnalysis selection or 'group NAME'")
    p.add_argument("--ligand-selection", default=None, help="MDAnalysis selection or 'group NAME'")
    p.add_argument("--frame-step", type=int, default=10, help="Process every N-th frame")
    p.add_argument("--max-frames", type=int, default=None, help="Maximum frames")
    p.add_argument("--temperature", type=float, default=300.0, help="Temperature (K)")
    p.add_argument("--pH", type=float, default=7.0, help="pH (kept for compatibility)")
    p.add_argument("--ionic-strength", type=float, default=0.15, help="Ionic strength (M)")
    p.add_argument("--surface-tension", type=float, default=0.0072,
                   help="Surface tension (kcal/mol/Å²)")
    p.add_argument("--pdie", type=float, default=2.0, help="Solute dielectric for APBS")
    p.add_argument("--sdie", type=float, default=78.54, help="Solvent dielectric for APBS")
    p.add_argument("--use-apbs", action="store_true", default=True, help="Use APBS for polar solvation")
    p.add_argument("--no-apbs", action="store_false", dest="use_apbs", help="Skip APBS")
    p.add_argument("--decompose", action="store_true", default=False, help="Per-residue decomposition")
    p.add_argument("--keep-temp", action="store_true", help="Keep per-frame working directory")
    return p.parse_args(argv)


def main(argv=None):
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stdout,
    )
    args = parse_args(argv)
    try:
        calculate_mmpbsa(args)
    except SystemExit:
        print("MMPBSA_STATUS=failed")
        raise
    except Exception as exc:
        log.error("%s", exc)
        print("MMPBSA_STATUS=failed")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
