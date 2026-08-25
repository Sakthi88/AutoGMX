#!/usr/bin/env python3
"""MM-PBSA binding free-energy calculation using GROMACS energy extraction.

Uses GROMACS tools for accurate MM energies (real force-field parameters,
proper PBC, correct exclusions) and APBS for Poisson-Boltzmann solvation.

Workflow per frame:
  1. Write frame to PDB
  2. Run PDB2PQR for charge/radius assignment
  3. Run APBS for polar solvation energy
  4. Compute SASA for non-polar solvation
  5. Extract MM interaction energy via GROMACS gmx energy

Output:
    mmpbsa_results.csv   -- per-frame component breakdown
    mmpbsa_summary.csv   -- averaged results with standard deviations
    mmpbsa_summary.json  -- machine-readable summary
    decomposition.csv    -- per-residue energy decomposition (if enabled)
"""
import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

import numpy as np

try:
    import MDAnalysis as mda
except ImportError as exc:
    raise SystemExit(f"MDAnalysis is required: {exc}")


# ── Constants ──────────────────────────────────────────────────────────────────
NM_TO_ANG = 10.0
KJ_TO_KCAL = 1.0 / 4.184
BOLTZMANN_KCAL = 0.001987204  # kT in kcal/mol at 300K


# ── GROMACS helpers ────────────────────────────────────────────────────────────
def _run(cmd, cwd=None, input_text=None, check=True):
    result = subprocess.run(
        cmd, cwd=cwd, capture_output=True, text=True,
        input=input_text, check=False,
    )
    if check and result.returncode != 0:
        raise RuntimeError(
            f"Command failed ({result.returncode}): {' '.join(cmd)}\n"
            f"stderr: {result.stderr[:2000]}"
        )
    return result


def _gmx_version(gmx_bin):
    r = _run([gmx_bin, "--version"], check=False)
    m = re.search(r"VERSION\s+(\S+)", r.stdout + r.stderr)
    return m.group(1) if m else "unknown"


def _build_energy_mdp(mdp_template_path, out_mdp_path):
    """Create an MDP file with energy groups for energy extraction."""
    content = mdp_template_path.read_text()

    # Remove any existing energygrps or continuation lines to avoid duplicates
    lines = content.splitlines()
    filtered = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith(("energygrps", "continuation")):
            continue
        filtered.append(line)

    additions = [
        "energygrps = Receptor Ligand",
        "continuation = yes",
    ]
    content = "\n".join(filtered).rstrip() + "\n" + "\n".join(additions) + "\n"
    out_mdp_path.write_text(content)


def _extract_interaction_energy(gmx_bin, edr_path, tpr_path):
    """Extract Receptor-Ligand interaction energy using gmx energy.

    Uses gmx energy with -s (TPR) to discover pair-specific energy terms,
    then extracts them. The TPR must have energygrps = Receptor Ligand.
    Returns (vdw_kj, elec_kj, total_kj) in kJ/mol.
    """
    edr_path = Path(edr_path)
    tpr_path = Path(tpr_path)
    energy_xvg = edr_path.parent / "energy.xvg"

    # Step 1: Discover available energy terms
    # Run gmx energy with -s to get the interactive term list
    discover_result = _run(
        [gmx_bin, "energy", "-f", str(edr_path), "-s", str(tpr_path)],
        input_text="\n",  # empty selection to get term list then exit
        check=False,
    )

    # Parse the term list from combined output
    # GROMACS prints the term list to its output (stdout or stderr varies)
    output_text = discover_result.stdout + "\n" + discover_result.stderr

    available = {}
    for line in output_text.split("\n"):
        line = line.strip()
        # Match: "  7  LJ-(SR)" or "  7  Coul-SR:Receptor-Ligand"
        # GROMACS 2025.2 format: "  N  Term-Name"
        m = re.match(r"^\s*(\d+)\s+(\S+(?:\s+\S+)?)\s*$", line)
        if m:
            idx = int(m.group(1))
            name = m.group(2).strip()
            if idx > 0 and name:
                available[name] = idx
        # Also handle tabular format: "  7  LJ-(SR)     8  Angle"
        for pair in re.findall(r"(\d+)\s+([\w().\-:/]+)", line):
            idx = int(pair[0])
            name = pair[1].strip()
            if idx > 0 and name and name not in ("kJ/mol", "K", "bar", "nm",
                                                   "nm^3", "kg/m^3"):
                available[name] = idx

    # Find Receptor-Ligand interaction terms
    coul_key = None
    lj_key = None
    for name, idx in available.items():
        if "Coul-SR" in name and "Receptor" in name and "Ligand" in name:
            coul_key = name
        if "LJ-SR" in name and "Receptor" in name and "Ligand" in name:
            lj_key = name

    if not coul_key and not lj_key:
        # Broader matching
        for name, idx in available.items():
            if "Coul" in name and "Receptor" in name and "Ligand" in name:
                coul_key = name
            if "LJ" in name and "Receptor" in name and "Ligand" in name:
                lj_key = name

    if not coul_key and not lj_key:
        raise RuntimeError(
            f"Could not find Receptor-Ligand interaction terms. "
            f"Available terms: {list(available.keys())[:30]}"
        )

    # Build selection: term indices followed by 0 to exit
    sel_indices = []
    if coul_key:
        sel_indices.append(str(available[coul_key]))
    if lj_key:
        sel_indices.append(str(available[lj_key]))
    sel_indices.append("0")

    # Step 2: Extract the selected energies
    _run(
        [gmx_bin, "energy", "-f", str(edr_path), "-s", str(tpr_path),
         "-o", str(energy_xvg)],
        cwd=str(edr_path.parent),
        input_text="\n".join(sel_indices) + "\n",
        check=False,
    )

    if not energy_xvg.exists():
        return 0.0, 0.0, 0.0

    # Step 3: Parse XVG - use legend to match columns to terms
    legend = []
    with open(energy_xvg, "r") as f:
        for line in f:
            if line.startswith("@"):
                if "legend" in line:
                    m = re.search(r'@"(\d+)"\s+"(.+)"', line)
                    if m:
                        legend.append((int(m.group(1)), m.group(2)))
            elif line.startswith("#") or not line.strip():
                continue
            else:
                break

    # Map legend columns to our terms
    coul_col = -1
    lj_col = -1
    for col_idx, name in legend:
        if "Coul-SR" in name and "Receptor" in name and "Ligand" in name:
            coul_col = col_idx
        if "LJ-SR" in name and "Receptor" in name and "Ligand" in name:
            lj_col = col_idx

    # Parse data rows
    coul_vals = []
    lj_vals = []
    with open(energy_xvg, "r") as f:
        for line in f:
            if line.startswith(("@", "#")) or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    vals = [float(x) for x in parts[1:]]  # skip time column
                    if coul_col > 0 and (coul_col - 1) < len(vals):
                        coul_vals.append(vals[coul_col - 1])
                    if lj_col > 0 and (lj_col - 1) < len(vals):
                        lj_vals.append(vals[lj_col - 1])
                    # Fallback for 2-column data: assume LJ, Coul order
                    if coul_col < 0 and lj_col < 0:
                        if len(vals) >= 2:
                            lj_vals.append(vals[0])
                            coul_vals.append(vals[1])
                        elif len(vals) == 1:
                            lj_vals.append(vals[0])
                except ValueError:
                    continue

    vdw = float(np.mean(lj_vals)) if lj_vals else 0.0
    elec = float(np.mean(coul_vals)) if coul_vals else 0.0
    return vdw, elec, vdw + elec


# ── PDB2PQR / APBS helpers ───────────────────────────────────────────────────
def _run_pdb2pqr(pdb_path, pqr_path, pH=7.0):
    """Assign charges and radii with PDB2PQR."""
    _run([
        "pdb2pqr", "--ff", "AMBER",
        "--with-ph", str(pH),
        str(pdb_path), str(pqr_path),
    ])


def _run_apbs(pqr_path, output_prefix, temp=300.0, ionic=0.15, pdie=1.0, sdie=78.5):
    """Run APBS to solve the Poisson-Boltzmann equation."""
    pqr = Path(pqr_path)
    pqr_content = pqr.read_text()

    coords = []
    for line in pqr_content.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            coords.append([x, y, z])

    coords = np.array(coords)
    coords_nm = coords / NM_TO_ANG
    lo = coords_nm.min(axis=0) - 0.3
    hi = coords_nm.max(axis=0) + 0.3

    # Correct APBS DX write syntax: "write dx <name> <filename>"
    apbs_input = textwrap.dedent(f"""\
        read mol pqr {pqr.name}
        elec name comp
            mg-auto
            molecule
            npb Solver
            bcfl mdh
            pdie {pdie}
            sdie {sdie}
            srfm smol
            chgm spl2
            ion 1 {ionic} 2.0
            ion -1 {ionic} 2.0
            temp {temp}
            write dx comp {output_prefix}
        end
        quit
    """)

    apbs_in = pqr.with_suffix(".in")
    apbs_in.write_text(apbs_input)

    # Run APBS from the PQR's directory so DX is written there
    _run(["apbs", str(apbs_in)], cwd=str(pqr.parent))

    # APBS writes to <output_prefix>.dx in the cwd
    dx_path = pqr.parent / f"{output_prefix}.dx"
    return dx_path


def _compute_elec_energy_pqr(pqr_path, dx_path, temp=300.0):
    """Compute electrostatic energy from PQR charges and DX potential grid.

    APBS DX potential is in kT/e units. We convert to kcal/mol using:
        E = q * phi * (RT/F) where RT/F ≈ 0.596 kcal/mol per e at 300K
        Or equivalently: E_kcal = q * phi_kT * kT_kcal
    """
    kT_kcal = BOLTZMANN_KCAL * temp  # kT in kcal/mol

    pqr_content = Path(pqr_path).read_text()
    coords = []
    charges = []
    for line in pqr_content.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            q = float(line[54:66])
            # Convert Å to nm for grid interpolation
            coords.append([x / NM_TO_ANG, y / NM_TO_ANG, z / NM_TO_ANG])
            charges.append(q)

    coords = np.array(coords)
    charges = np.array(charges)

    dx_content = Path(dx_path).read_text()
    header_lines = []
    data_values = []
    for line in dx_content.splitlines():
        stripped = line.strip()
        if stripped.startswith(("object", "component", "type", "grid", "origin",
                               "delta", "label")):
            header_lines.append(stripped)
        else:
            try:
                vals = [float(v) for v in stripped.split()]
                data_values.extend(vals)
            except ValueError:
                continue

    if not data_values:
        return 0.0

    nx = ny = nz = 0
    origin = np.zeros(3)
    delta = np.zeros(3)
    for line in header_lines:
        parts = line.split()
        if "object 1" in line and "contents" in line:
            nx, ny, nz = int(parts[-3]), int(parts[-2]), int(parts[-1])
        elif "origin" in line:
            origin = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
        elif "delta" in line and len(parts) >= 5:
            idx = int(parts[1])
            delta[idx] = float(parts[2])

    if nx * ny * nz == 0 or len(data_values) < nx * ny * nz:
        return 0.0

    # DX files use column-major (Fortran) order: x changes fastest
    grid = np.array(data_values[:nx * ny * nz]).reshape((nz, ny, nx), order="F")

    energy = 0.0
    for i in range(len(charges)):
        x, y, z = coords[i] - origin
        ix = x / delta[0] if delta[0] > 0 else 0
        iy = y / delta[1] if delta[1] > 0 else 0
        iz = z / delta[2] if delta[2] > 0 else 0

        ix0 = int(np.clip(np.floor(ix), 0, nx - 2))
        iy0 = int(np.clip(np.floor(iy), 0, ny - 2))
        iz0 = int(np.clip(np.floor(iz), 0, nz - 2))

        fx = ix - ix0
        fy = iy - iy0
        fz = iz - iz0

        # Trilinear interpolation
        pot = (
            grid[iz0, iy0, ix0] * (1 - fx) * (1 - fy) * (1 - fz) +
            grid[iz0, iy0, ix0 + 1] * fx * (1 - fy) * (1 - fz) +
            grid[iz0, iy0 + 1, ix0] * (1 - fx) * fy * (1 - fz) +
            grid[iz0, iy0 + 1, ix0 + 1] * fx * fy * (1 - fz) +
            grid[iz0 + 1, iy0, ix0] * (1 - fx) * (1 - fy) * fz +
            grid[iz0 + 1, iy0, ix0 + 1] * fx * (1 - fy) * fz +
            grid[iz0 + 1, iy0 + 1, ix0] * (1 - fx) * fy * fz +
            grid[iz0 + 1, iy0 + 1, ix0 + 1] * fx * fy * fz
        )
        # pot is in kT/e, charges in e, result in kT -> convert to kcal/mol
        energy += charges[i] * pot * kT_kcal

    return energy


# ── SASA ──────────────────────────────────────────────────────────────────────
def _compute_sa(universe, selection, probe=1.4):
    """Compute SASA using MDAnalysis geometry-based approach."""
    ag = universe.select_atoms(selection)
    if len(ag) == 0:
        return 0.0

    # Use MDAnalysis' built-in SASA if available
    try:
        from MDAnalysis.analysis import sa
        sasa = sa.SASA(ag, probe_radius=probe)
        sasa.run()
        return float(sasa.results.area)
    except (ImportError, AttributeError, TypeError):
        pass

    # Fallback: Shrake-Rupley with KDTree neighbor search
    from scipy.spatial import cKDTree
    positions = ag.positions
    n_atoms = len(ag)
    
    # Use element-based radii
    radii = np.array([_vdw_radius(ag[i].type) for i in range(n_atoms)])
    probe_r = probe
    
    # Generate points on sphere for each atom
    n_points = 960
    sphere_pts = _fibonacci_sphere(n_points)
    
    # Build KDTree for neighbor search
    tree = cKDTree(positions)
    
    total_area = 0.0
    for i in range(n_atoms):
        r = radii[i] + probe_r
        atom_pts = positions[i] + sphere_pts * r
        
        # Find neighbors within 2*max_radius
        max_r = np.max(radii) + probe_r
        neighbors = tree.query_ball_point(positions[i], 2 * max_r)
        
        # Check each point for burial
        accessible = 0
        for pt in atom_pts:
            buried = False
            for j in neighbors:
                if j == i:
                    continue
                d = np.linalg.norm(pt - positions[j])
                if d < radii[j] + probe_r:
                    buried = True
                    break
            if not buried:
                accessible += 1
        
        exposed_fraction = accessible / n_points
        total_area += exposed_fraction * 4.0 * np.pi * r ** 2
    
    return float(total_area)


def _vdw_radius(atom_type):
    """Get van der Waals radius for atom type (in Å)."""
    radii = {
        'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80,
        'P': 1.80, 'F': 1.47, 'CL': 1.75, 'BR': 1.85, 'I': 1.98,
    }
    at = atom_type.upper()
    for elem, rad in radii.items():
        if at.startswith(elem):
            return rad
    return 1.70  # default carbon


def _fibonacci_sphere(n_points):
    """Generate points on unit sphere using Fibonacci lattice."""
    indices = np.arange(n_points, dtype=float) + 0.5
    phi = np.arccos(1 - 2 * indices / n_points)
    theta = np.pi * (1 + 5**0.5) * indices
    x = np.cos(theta) * np.sin(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(phi)
    return np.column_stack([x, y, z])
    return 1.70  # default carbon


def _fibonacci_sphere(n_points):
    """Generate points on unit sphere using Fibonacci lattice."""
    indices = np.arange(n_points, dtype=float) + 0.5
    phi = np.arccos(1 - 2 * indices / n_points)
    theta = np.pi * (1 + 5**0.5) * indices
    x = np.cos(theta) * np.sin(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(phi)
    return np.column_stack([x, y, z])


# ── Main MM-PBSA calculation ──────────────────────────────────────────────────
def calculate_mmpbsa(args):
    gmx_bin = args.gmx_bin
    print(f"GROMACS: {gmx_bin} ({_gmx_version(gmx_bin)})")
    print(f"Trajectory: {args.trajectory}")
    print(f"Topology: {args.topology}")

    # Load index file if provided for group selections
    ndx_groups = {}
    if args.index_file and Path(args.index_file).exists():
        # Parse the index file manually
        with open(args.index_file, 'r') as f:
            current_group = None
            for line in f:
                line = line.strip()
                if line.startswith('[') and line.endswith(']'):
                    current_group = line[1:-1].strip()
                    ndx_groups[current_group] = []
                elif current_group and line:
                    ndx_groups[current_group].extend([int(x) - 1 for x in line.split()])  # 0-based
        print(f"Loaded index groups: {list(ndx_groups.keys())}")

    u = mda.Universe(args.topology, args.trajectory)
    n_frames = len(u.trajectory)
    print(f"Trajectory contains {n_frames} frames")

    sel_receptor = args.receptor_selection or "protein"
    sel_ligand = args.ligand_selection or "resname LIG"

    # Handle "group X" selections using loaded index
    def resolve_selection(sel):
        if sel.startswith("group "):
            group_name = sel[6:].strip()
            if group_name in ndx_groups:
                indices = ndx_groups[group_name]
                # Create a selection string from indices
                return f"index {' '.join(str(i) for i in indices)}"
            else:
                print(f"Warning: Group '{group_name}' not found in index file, available: {list(ndx_groups.keys())}")
                return "resname UNK"  # fallback
        return sel

    sel_receptor = resolve_selection(sel_receptor)
    sel_ligand = resolve_selection(sel_ligand)

    n_receptor = u.select_atoms(sel_receptor).n_atoms
    n_ligand = u.select_atoms(sel_ligand).n_atoms
    print(f"Receptor: {n_receptor} atoms | Ligand: {n_ligand} atoms")

    if n_receptor == 0 or n_ligand == 0:
        raise SystemExit("One or more selections returned zero atoms.")

    # Frame sampling
    frame_step = max(1, args.frame_step)
    frame_indices = list(range(0, n_frames, frame_step))
    if args.max_frames and len(frame_indices) > args.max_frames:
        frame_indices = frame_indices[:args.max_frames]
    print(f"Processing {len(frame_indices)} frames (every {frame_step} frame(s))")

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(tempfile.mkdtemp(prefix="mmpbsa_"))

    try:
        # ── Step 1: Create index with energy groups ───────────────────────
        print("\n--- Setting up energy groups ---")
        receptor_atoms = u.select_atoms(sel_receptor)
        ligand_atoms = u.select_atoms(sel_ligand)

        # Write index file with proper line wrapping (15 atoms per line)
        ndx_path = temp_dir / "index.ndx"
        ndx_lines = []

        def _write_ndx_group(name, atom_indices):
            ndx_lines.append(f"[ {name} ]")
            line = ""
            for i, idx in enumerate(atom_indices):
                line += f" {idx + 1}"  # GROMACS is 1-based
                if (i + 1) % 15 == 0:
                    ndx_lines.append(line.strip())
                    line = ""
            if line.strip():
                ndx_lines.append(line.strip())

        _write_ndx_group("System", range(u.atoms.n_atoms))
        _write_ndx_group("Receptor", (a.index for a in receptor_atoms))
        _write_ndx_group("Ligand", (a.index for a in ligand_atoms))
        ndx_path.write_text("\n".join(ndx_lines) + "\n")

        print(f"Index groups: Receptor ({n_receptor} atoms), Ligand ({n_ligand} atoms)")

        # ── Step 2: Create MDP with energy groups ─────────────────────────
        base_mdp = Path(args.mdp_file) if args.mdp_file else None
        if not base_mdp or not base_mdp.exists():
            base_mdp = temp_dir / "minimal.mdp"
            base_mdp.write_text(textwrap.dedent("""\
                integrator = md
                nsteps = 0
                dt = 0.002
                nstxout = 0
                nstvout = 0
                nstfout = 0
                nstlog = 1
                nstenergy = 1
                nstxout-compressed = 0
                continuation = yes
                constraint_algorithm = lincs
                constraints = h-bonds
                cutoff-scheme = Verlet
                nstlist = 1
                rlist = 1.2
                coulombtype = PME
                rcoulomb = 1.2
                vdwtype = Cut-off
                rvdw = 1.2
                pbc = xyz
            """))

        energy_mdp = temp_dir / "energy.mdp"
        _build_energy_mdp(base_mdp, energy_mdp)

        # ── Step 3: Find .top file for grompp (TPR cannot be used with -p) ───────────
        top_path = None
        if args.topology.endswith(".tpr"):
            top_dir = Path(args.topology).parent
            # Copy all topology files from production directory for grompp
            for itp_file in top_dir.glob("*.itp"):
                shutil.copy2(itp_file, temp_dir / itp_file.name)
            for top_file in top_dir.glob("*.top"):
                shutil.copy2(top_file, temp_dir / top_file.name)
            # Also copy local forcefield directory if exists
            for ff_dir in top_dir.glob("*.ff"):
                if ff_dir.is_dir():
                    shutil.copytree(ff_dir, temp_dir / ff_dir.name, dirs_exist_ok=True)

            # Look for existing .top file
            for candidate in [top_dir / "topol.top", top_dir / "system.top"]:
                if candidate.exists():
                    top_path = candidate
                    break
            if not top_path:
                # Extract topology from TPR using gmx dump
                top_path = temp_dir / "extracted.top"
                _run([
                    gmx_bin, "dump", "-s", str(args.topology), "-o", str(top_path)
                ], check=False)
                if not top_path.exists() or top_path.stat().st_size == 0:
                    raise RuntimeError("Could not extract topology from TPR")
        else:
            top_path = Path(args.topology)

        # Pre-generate a TPR with energy groups for reuse
        energy_tpr = temp_dir / "energy.tpr"
        _run([
            gmx_bin, "grompp",
            "-f", str(energy_mdp),
            "-c", str(args.topology),  # TPR contains coordinates
            "-p", str(top_path),       # Use .top file
            "-n", str(ndx_path),
            "-o", str(energy_tpr),
            "-maxwarn", "5",
        ], check=False)

        base_tpr = energy_tpr if energy_tpr.exists() else None

        # ── Step 4: For each frame, extract energies ──────────────────────
        results = []
        pdb_complex = temp_dir / "complex.pdb"
        pdb_receptor = temp_dir / "receptor.pdb"
        pdb_ligand = temp_dir / "ligand.pdb"

        for idx, fidx in enumerate(frame_indices):
            u.trajectory[fidx]
            frame = u.trajectory.frame
            time_ps = u.trajectory.time
            print(f"\n  Frame {idx + 1}/{len(frame_indices)} (frame={frame}, t={time_ps:.1f} ps)")

            # Write frame PDB
            u.select_atoms("all").write(str(pdb_complex))
            receptor_atoms.write(str(pdb_receptor))
            ligand_atoms.write(str(pdb_ligand))

            # ── MM energy via GROMACS ─────────────────────────────────────
            frame_work = temp_dir / f"frame_{fidx:06d}"
            frame_work.mkdir(exist_ok=True)
            shutil.copy2(pdb_complex, frame_work / "complex.pdb")
            shutil.copy2(ndx_path, frame_work / "index.ndx")

            # grompp to create frame TPR with energy groups
            tpr_frame = frame_work / "frame.tpr"
            # Use topology file from temp_dir (copied with all includes)
            frame_top = temp_dir / top_path.name
            _run([
                gmx_bin, "grompp",
                "-f", str(energy_mdp),
                "-c", str(frame_work / "complex.pdb"),
                "-p", str(frame_top),
                "-n", str(frame_work / "index.ndx"),
                "-o", str(tpr_frame),
                "-maxwarn", "5",
            ], cwd=str(frame_work), check=False)

            if not tpr_frame.exists():
                print("    Warning: grompp failed, skipping MM energy for this frame")
                vdw_kj, elec_kj = 0.0, 0.0
            else:
                # mdrun -rerun to compute energies (no -n flag; index is in TPR)
                rerun_edr = frame_work / "frame.edr"
                _run([
                    gmx_bin, "mdrun",
                    "-s", str(tpr_frame),
                    "-rerun", str(frame_work / "complex.pdb"),
                    "-e", str(rerun_edr),
                    "-o", str(frame_work / "frame.xtc"),
                ], cwd=str(frame_work), check=False)

                # Extract interaction energies
                vdw_kj, elec_kj = 0.0, 0.0
                try:
                    vdw_kj, elec_kj, _ = _extract_interaction_energy(
                        gmx_bin, rerun_edr, tpr_frame
                    )
                    print(f"    MM: vdW={vdw_kj * KJ_TO_KCAL:+.2f}  Elec={elec_kj * KJ_TO_KCAL:+.2f} kcal/mol")
                except Exception as exc:
                    print(f"    Warning: MM energy extraction failed: {exc}")

            # Convert to kcal/mol
            vdw = vdw_kj * KJ_TO_KCAL
            elec = elec_kj * KJ_TO_KCAL

            # ── Polar solvation (PB) via APBS ────────────────────────────
            polar_solv = 0.0
            if args.use_apbs:
                try:
                    pqr_complex = frame_work / "complex.pqr"
                    pqr_receptor = frame_work / "receptor.pqr"
                    pqr_ligand = frame_work / "ligand.pqr"

                    _run_pdb2pqr(pdb_complex, pqr_complex, pH=args.pH)
                    _run_pdb2pqr(pdb_receptor, pqr_receptor, pH=args.pH)
                    _run_pdb2pqr(pdb_ligand, pqr_ligand, pH=args.pH)

                    dx_c = _run_apbs(pqr_complex, "complex",
                                     temp=args.temperature, ionic=args.ionic_strength)
                    dx_r = _run_apbs(pqr_receptor, "receptor",
                                     temp=args.temperature, ionic=args.ionic_strength)
                    dx_l = _run_apbs(pqr_ligand, "ligand",
                                     temp=args.temperature, ionic=args.ionic_strength)

                    e_c = _compute_elec_energy_pqr(pqr_complex, dx_c, temp=args.temperature)
                    e_r = _compute_elec_energy_pqr(pqr_receptor, dx_r, temp=args.temperature)
                    e_l = _compute_elec_energy_pqr(pqr_ligand, dx_l, temp=args.temperature)
                    polar_solv = e_c - e_r - e_l
                    print(f"    PB: {polar_solv:+.2f} kcal/mol")
                except Exception as exc:
                    print(f"    Warning: PB calculation failed: {exc}")

            # ── Non-polar solvation (SA) ──────────────────────────────────
            sa_c = _compute_sa(u, "all")
            sa_r = _compute_sa(u, sel_receptor)
            sa_l = _compute_sa(u, sel_ligand)
            gamma = args.surface_tension
            nonpolar_solv = gamma * (sa_c - sa_r - sa_l)
            print(f"    SA: ΔSA={sa_c - sa_r - sa_l:.1f} Å²  ΔG_SA={nonpolar_solv:+.2f} kcal/mol")

            # ── Total ─────────────────────────────────────────────────────
            delta_g = vdw + elec + polar_solv + nonpolar_solv
            print(f"    ΔG_total = {delta_g:+.2f} kcal/mol")

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

    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)

    if not results:
        raise SystemExit("No frames were processed.")

    # ── Write per-frame CSV ───────────────────────────────────────────────
    csv_path = out_dir / "mmpbsa_results.csv"
    fieldnames = list(results[0].keys())
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    print(f"\nPer-frame results: {csv_path}")

    # ── Summary statistics ────────────────────────────────────────────────
    components = ["vdw", "elec", "polar_solv", "nonpolar_solv", "delta_G"]
    summary = {}
    for comp in components:
        values = np.array([r[comp] for r in results])
        summary[comp] = {
            "mean": float(np.mean(values)),
            "std": float(np.std(values)),
            "median": float(np.median(values)),
            "min": float(np.min(values)),
            "max": float(np.max(values)),
            "sem": float(np.std(values) / np.sqrt(len(values))),
        }

    summary_csv = out_dir / "mmpbsa_summary.csv"
    with open(summary_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["Component", "Mean (kcal/mol)", "Std Dev", "Median", "Min", "Max", "SEM"])
        for comp in components:
            s = summary[comp]
            writer.writerow([comp, f"{s['mean']:.4f}", f"{s['std']:.4f}",
                             f"{s['median']:.4f}", f"{s['min']:.4f}",
                             f"{s['max']:.4f}", f"{s['sem']:.4f}"])
    print(f"Summary: {summary_csv}")

    summary_json = {
        "n_frames": len(results),
        "frame_step": frame_step,
        "temperature": args.temperature,
        "pH": args.pH,
        "ionic_strength": args.ionic_strength,
        "surface_tension": args.surface_tension,
        "energy_method": "GROMACS (force-field parameters from TPR)",
        "pb_method": "APBS Poisson-Boltzmann" if args.use_apbs else "skipped",
        "selections": {
            "receptor": sel_receptor,
            "ligand": sel_ligand,
        },
        "components": summary,
    }
    json_path = out_dir / "mmpbsa_summary.json"
    json_path.write_text(json.dumps(summary_json, indent=2))
    print(f"JSON summary: {json_path}")

    # ── Per-residue decomposition ─────────────────────────────────────────
    if args.decompose:
        print("\nRunning per-residue decomposition...")
        _per_residue_decomposition(u, sel_receptor, sel_ligand, frame_indices, out_dir, args)

    # ── Plots ─────────────────────────────────────────────────────────────
    _generate_plots(results, summary, out_dir)

    print("\nMM-PBSA calculation complete.")
    return summary_json


def _per_residue_decomposition(u, sel_receptor, sel_ligand, frame_indices, out_dir, args):
    """Per-residue decomposition using pair-wise distances and charges.

    Note: This uses generic LJ parameters since force-field-specific
    parameters are not easily extracted per-atom from the trajectory.
    The electrostatic component uses actual charges from the topology.
    """
    rec_atoms = u.select_atoms(sel_receptor)
    residues = rec_atoms.residues
    decomp_results = []

    # Pre-select ligand atoms once
    lig_atoms = u.select_atoms(sel_ligand)

    for res in residues:
        res_name = f"{res.resname}{res.resid}"
        res_atoms = res.atoms

        vdw_list = []
        elec_list = []

        # Use configurable frame limit
        max_decomp_frames = min(20, len(frame_indices))
        for fidx in frame_indices[:max_decomp_frames]:
            u.trajectory[fidx]
            rec_pos = res_atoms.positions / NM_TO_ANG  # Å to nm
            lig_pos = lig_atoms.positions / NM_TO_ANG

            # PBC-aware minimum image distances
            box = u.trajectory.dimensions[:3] / NM_TO_ANG  # box in nm
            energy = 0.0
            for i in range(len(lig_pos)):
                diffs = rec_pos - lig_pos[i]
                # Minimum image convention
                diffs = diffs - box * np.round(diffs / box)
                dists = np.sqrt(np.sum(diffs ** 2, axis=1))
                mask = dists < 1.2  # 1.2 nm cutoff
                if not np.any(mask):
                    continue
                close_dists = dists[mask]
                # Generic LJ params (not force-field specific, but avoids crash)
                sigma = 0.35
                epsilon = 0.1
                sr6 = (sigma / close_dists) ** 6
                energy += 4.0 * epsilon * np.sum(sr6 ** 2 - sr6)
            vdw_list.append(energy)

            # Electrostatics using actual charges if available
            try:
                rec_q = res_atoms.charges
                lig_q = lig_atoms.charges
                if len(rec_q) > 0 and len(lig_q) > 0:
                    e_elec = 0.0
                    for i in range(len(lig_pos)):
                        diffs = rec_pos - lig_pos[i]
                        diffs = diffs - box * np.round(diffs / box)
                        dists = np.sqrt(np.sum(diffs ** 2, axis=1))
                        mask = dists > 0.01
                        if np.any(mask):
                            e_elec += 332.0636 * np.sum(rec_q[mask] * lig_q[i] / dists[mask])
                    elec_list.append(e_elec)
                else:
                    elec_list.append(0.0)
            except Exception:
                elec_list.append(0.0)

        # Safe mean calculation
        vdw_mean = float(np.mean(vdw_list)) if vdw_list else 0.0
        elec_mean = float(np.mean(elec_list)) if elec_list else 0.0
        if vdw_list and elec_list and len(vdw_list) == len(elec_list):
            total_mean = float(np.mean(np.array(vdw_list) + np.array(elec_list)))
        else:
            total_mean = vdw_mean + elec_mean

        decomp_results.append({
            "residue": res_name,
            "resid": int(res.resid),
            "resname": res.resname,
            "vdw_mean": vdw_mean,
            "elec_mean": elec_mean,
            "total_mean": total_mean,
        })

    csv_path = out_dir / "decomposition.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["residue", "resid", "resname",
                                                "vdw_mean", "elec_mean", "total_mean"])
        writer.writeheader()
        writer.writerows(decomp_results)
    print(f"Per-residue decomposition: {csv_path}")

    decomp_results.sort(key=lambda x: x["total_mean"])
    summary_path = out_dir / "decomposition_summary.json"
    summary_path.write_text(json.dumps({
        "top_binding": decomp_results[:10],
        "top_repulsive": decomp_results[-10:][::-1],
    }, indent=2))

    # Plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(14, 6))
        totals = [d["total_mean"] for d in decomp_results]
        colors = ["#16a34a" if t < 0 else "#dc2626" for t in totals]
        ax.bar(range(len(totals)), totals, color=colors, width=1.0, edgecolor="none")
        res_names = [d["residue"] for d in decomp_results]
        step = max(1, len(res_names) // 20)
        ax.set_xticks(range(0, len(res_names), step))
        ax.set_xticklabels([res_names[i] for i in range(0, len(res_names), step)],
                           rotation=45, ha="right", fontsize=8)
        ax.set_ylabel("ΔG Contribution (kcal/mol)")
        ax.set_title("Per-Residue Energy Decomposition")
        ax.axhline(y=0, color="black", linewidth=0.5)
        ax.grid(axis="y", alpha=0.3)
        fig.tight_layout()
        fig.savefig(str(out_dir / "decomposition_plot.png"), dpi=200)
        plt.close(fig)
    except Exception as exc:
        print(f"Warning: Could not generate decomposition plot: {exc}")


def _generate_plots(results, summary, out_dir):
    """Generate analysis plots."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("Warning: matplotlib not available, skipping plots.")
        return

    components = ["vdw", "elec", "polar_solv", "nonpolar_solv", "delta_G"]
    labels = ["vdW", "Electrostatics", "Polar Solv.", "Nonpolar Solv.", "Total ΔG"]
    colors = ["#2563eb", "#dc2626", "#16a34a", "#ca8a04", "#7c3aed"]

    # Bar chart
    fig, ax = plt.subplots(figsize=(10, 6))
    means = [summary[c]["mean"] for c in components]
    stds = [summary[c]["std"] for c in components]
    bars = ax.bar(labels, means, yerr=stds, capsize=5, color=colors,
                  edgecolor="black", linewidth=0.5)
    ax.set_ylabel("Energy (kcal/mol)", fontsize=12)
    ax.set_title("MM-PBSA Binding Free Energy Decomposition", fontsize=14)
    ax.axhline(y=0, color="black", linewidth=0.5)
    ax.grid(axis="y", alpha=0.3)
    for bar, val in zip(bars, means):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.5 * np.sign(bar.get_height()),
                f"{val:.2f}", ha="center",
                va="bottom" if val >= 0 else "top", fontsize=10)
    fig.tight_layout()
    fig.savefig(str(out_dir / "mmpbsa_decomposition.png"), dpi=200)
    plt.close(fig)
    print(f"Decomposition plot: {out_dir / 'mmpbsa_decomposition.png'}")

    # Time series
    times = [r["time_ps"] for r in results]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    plot_data = [
        ("vdw", "van der Waals", "#2563eb", axes[0, 0]),
        ("elec", "Electrostatics", "#dc2626", axes[0, 1]),
        ("polar_solv", "Polar Solvation", "#16a34a", axes[1, 0]),
        ("delta_G", "Total ΔG Binding", "#7c3aed", axes[1, 1]),
    ]
    for key, title, color, ax in plot_data:
        vals = [r[key] for r in results]
        ax.plot(times, vals, color=color, linewidth=0.8, alpha=0.7)
        ax.set_title(title, fontsize=11)
        ax.set_ylabel("kcal/mol")
        ax.grid(alpha=0.3)
        mean_val = np.mean(vals)
        ax.axhline(y=mean_val, color=color, linestyle="--", alpha=0.5,
                   label=f"mean={mean_val:.2f}")
        ax.legend(fontsize=9)
    axes[1, 0].set_xlabel("Time (ps)")
    axes[1, 1].set_xlabel("Time (ps)")
    fig.suptitle("MM-PBSA Energy Components Over Time", fontsize=13)
    fig.tight_layout()
    fig.savefig(str(out_dir / "mmpbsa_timeseries.png"), dpi=200)
    plt.close(fig)
    print(f"Time series plot: {out_dir / 'mmpbsa_timeseries.png'}")


# ── CLI ────────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="MM-PBSA binding free energy calculation")
    p.add_argument("-t", "--trajectory", required=True, help="GROMACS trajectory (.xtc)")
    p.add_argument("-s", "--topology", required=True, help="GROMACS topology (.tpr or .gro)")
    p.add_argument("-o", "--output-dir", required=True, help="Output directory")
    p.add_argument("-n", "--index-file", default=None, help="GROMACS index file (.ndx)")
    p.add_argument("-f", "--mdp-file", default=None, help="MDP template file")
    p.add_argument("--gmx-bin", default="gmx", help="GROMACS binary")
    p.add_argument("--receptor-selection", default=None, help="MDAnalysis selection for receptor")
    p.add_argument("--ligand-selection", default=None, help="MDAnalysis selection for ligand")
    p.add_argument("--frame-step", type=int, default=10, help="Process every N-th frame")
    p.add_argument("--max-frames", type=int, default=None, help="Maximum frames")
    p.add_argument("--temperature", type=float, default=300.0, help="Temperature (K)")
    p.add_argument("--pH", type=float, default=7.0, help="pH for protonation")
    p.add_argument("--ionic-strength", type=float, default=0.15, help="Ionic strength (M)")
    p.add_argument("--surface-tension", type=float, default=0.0072,
                   help="Surface tension (kcal/mol/Å²)")
    p.add_argument("--use-apbs", action="store_true", default=True,
                   help="Use APBS for polar solvation")
    p.add_argument("--no-apbs", action="store_false", dest="use_apbs",
                   help="Skip APBS")
    p.add_argument("--decompose", action="store_true", default=False,
                   help="Per-residue decomposition")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    calculate_mmpbsa(args)
