#!/usr/bin/env python3
import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np


R_KJ_MOL_K = 0.0083144621


def require_modules():
    missing = []
    modules = {}
    for name in ("MDAnalysis", "matplotlib", "sklearn", "scipy"):
        try:
            modules[name] = __import__(name)
        except ImportError:
            missing.append(name)
    if missing:
        raise RuntimeError(
            "Missing Python package(s): "
            + ", ".join(missing)
            + f". Python used: {sys.executable}. Install with: {sys.executable} -m pip install MDAnalysis matplotlib scikit-learn scipy"
        )
    return modules


def load_aligned_matrix(topology, trajectory, selection, start, stop, step, fallback_topology=None):
    import MDAnalysis as mda
    from MDAnalysis.analysis import align

    topology_used = Path(topology)
    try:
        universe = mda.Universe(str(topology_used), str(trajectory))
    except Exception:
        if not fallback_topology:
            raise
        topology_used = Path(fallback_topology)
        universe = mda.Universe(str(topology_used), str(trajectory))
    atoms = universe.select_atoms(selection)
    if atoms.n_atoms == 0:
        raise RuntimeError(f"No atoms matched selection: {selection}")

    reference = universe.copy()
    reference.trajectory[start]
    aligner = align.AlignTraj(universe, reference, select=selection, in_memory=True)
    aligner.run(start=start, stop=stop, step=step)

    coords = []
    for ts in universe.trajectory[start:stop:step]:
        coords.append(atoms.positions.astype(np.float64).ravel())
    matrix = np.asarray(coords)
    if matrix.shape[0] < 3:
        raise RuntimeError("Need at least 3 trajectory frames for PCA/free-energy analysis.")
    return matrix, atoms.n_atoms, len(universe.trajectory), topology_used


def run_pca(matrix, n_components):
    from sklearn.decomposition import PCA

    max_components = min(n_components, matrix.shape[0], matrix.shape[1])
    pca = PCA(n_components=max_components)
    scores = pca.fit_transform(matrix)
    return pca, scores


def kde_grid(scores, pcx, pcy, bins):
    from scipy.stats import gaussian_kde

    x = scores[:, pcx]
    y = scores[:, pcy]
    if np.allclose(x.min(), x.max()) or np.allclose(y.min(), y.max()):
        raise RuntimeError("PC projections have near-zero range; cannot build a 2D KDE landscape.")

    pad_x = 0.08 * (x.max() - x.min())
    pad_y = 0.08 * (y.max() - y.min())
    xi = np.linspace(x.min() - pad_x, x.max() + pad_x, bins)
    yi = np.linspace(y.min() - pad_y, y.max() + pad_y, bins)
    xg, yg = np.meshgrid(xi, yi)
    kde = gaussian_kde(np.vstack([x, y]))
    density = kde(np.vstack([xg.ravel(), yg.ravel()])).reshape(xg.shape)
    return xg, yg, density


def free_energy_from_density(density, temperature):
    k_t = R_KJ_MOL_K * temperature
    probability = density / np.sum(density)
    with np.errstate(divide="ignore", invalid="ignore"):
        energy = -k_t * np.log(probability)
    finite = np.isfinite(energy)
    if not np.any(finite):
        raise RuntimeError("Could not compute finite free energies from the density grid.")
    energy = energy - np.nanmin(energy[finite])
    energy[~finite] = np.nanmax(energy[finite])
    return energy


def one_dimensional_free_energy(values, temperature, bins):
    hist, edges = np.histogram(values, bins=bins, density=True)
    centers = 0.5 * (edges[1:] + edges[:-1])
    k_t = R_KJ_MOL_K * temperature
    with np.errstate(divide="ignore", invalid="ignore"):
        energy = -k_t * np.log(hist)
    finite = np.isfinite(energy)
    if np.any(finite):
        energy = energy - np.nanmin(energy[finite])
        energy[~finite] = np.nanmax(energy[finite])
    else:
        energy[:] = 0.0
    return centers, energy


def save_plots(outdir, scores, pca, xg, yg, density, free_energy, temperature, bins):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm

    frames = np.arange(scores.shape[0])

    fig, ax = plt.subplots(figsize=(6, 4))
    ratios = pca.explained_variance_ratio_ * 100.0
    ax.bar(np.arange(1, len(ratios) + 1), ratios, alpha=0.75)
    ax.plot(np.arange(1, len(ratios) + 1), np.cumsum(ratios), marker="o", color="black")
    ax.set_xlabel("Principal component")
    ax.set_ylabel("Variance explained (%)")
    ax.set_title("PCA variance")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "pca_scree_variance.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5))
    sc = ax.scatter(scores[:, 0], scores[:, 1], c=frames, cmap="viridis", s=10, linewidths=0)
    fig.colorbar(sc, ax=ax, label="Frame")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PC1 vs PC2 projection")
    fig.tight_layout()
    fig.savefig(outdir / "pc1_pc2_scatter.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5))
    cf = ax.contourf(xg, yg, density, levels=24, cmap="jet")
    ax.contour(xg, yg, density, levels=12, colors="black", linewidths=0.5)
    fig.colorbar(cf, ax=ax, label="Density")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("PC1/PC2 density")
    fig.tight_layout()
    fig.savefig(outdir / "pc1_pc2_kde_density.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5))
    fe = ax.contourf(xg, yg, free_energy, levels=24, cmap="jet")
    ax.contour(xg, yg, free_energy, levels=12, colors="white", linewidths=0.5)
    fig.colorbar(fe, ax=ax, label="Delta G (kJ/mol)")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(f"Free energy landscape at {temperature:g} K")
    fig.tight_layout()
    fig.savefig(outdir / "free_energy_landscape_2d.png", dpi=300)
    plt.close(fig)

    if scores.shape[1] >= 3:
        fig = plt.figure(figsize=(7, 6))
        ax = fig.add_subplot(111, projection="3d")
        p = ax.scatter(scores[:, 0], scores[:, 1], scores[:, 2], c=frames, cmap="viridis", s=10)
        fig.colorbar(p, ax=ax, shrink=0.65, label="Frame")
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.set_zlabel("PC3")
        ax.set_title("3D PCA projection")
        fig.tight_layout()
        fig.savefig(outdir / "pc1_pc2_pc3_scatter_3d.png", dpi=300)
        plt.close(fig)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection="3d")
    surf = ax.plot_surface(xg, yg, free_energy, cmap=cm.rainbow, linewidth=0, antialiased=True)
    fig.colorbar(surf, ax=ax, shrink=0.65, label="Delta G (kJ/mol)")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_zlabel("Delta G (kJ/mol)")
    ax.set_title("3D free energy surface")
    fig.tight_layout()
    fig.savefig(outdir / "free_energy_surface_3d.png", dpi=300)
    plt.close(fig)

    centers, energy_1d = one_dimensional_free_energy(scores[:, 0], temperature, bins)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(centers, energy_1d, lw=1.5)
    ax.set_xlabel("PC1")
    ax.set_ylabel("Delta G (kJ/mol)")
    ax.set_title("1D free energy along PC1")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "free_energy_pc1_1d.png", dpi=300)
    plt.close(fig)

    np.savetxt(outdir / "free_energy_pc1_1d.csv", np.column_stack([centers, energy_1d]), delimiter=",", header="PC1,DeltaG_kJ_mol", comments="")


def write_failure(outdir, message):
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "analysis_summary.txt").write_text(
        "PCA Free Energy Analysis\n"
        "Status: skipped\n\n"
        f"{message}\n",
        encoding="utf-8",
    )


def main():
    parser = argparse.ArgumentParser(description="PCA-based free-energy landscape analysis for GROMACS trajectories")
    parser.add_argument("--top", required=True, help="Topology file, usually md.tpr")
    parser.add_argument("--fallback-top", default=None, help="Fallback topology, usually md.gro, if the TPR parser fails")
    parser.add_argument("--traj", required=True, help="Trajectory file, preferably centered md_center.xtc")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--selection", default="backbone", help="MDAnalysis selection string")
    parser.add_argument("--temperature", type=float, default=300.0, help="Temperature in K")
    parser.add_argument("--components", type=int, default=3, help="Number of PCs")
    parser.add_argument("--bins", type=int, default=180, help="Grid/histogram bins")
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--stop", type=int, default=None)
    parser.add_argument("--step", type=int, default=1)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    try:
        require_modules()
        outdir.mkdir(parents=True, exist_ok=True)
        matrix, n_atoms, total_frames, topology_used = load_aligned_matrix(
            Path(args.top),
            Path(args.traj),
            args.selection,
            args.start,
            args.stop,
            args.step,
            fallback_topology=args.fallback_top,
        )
        pca, scores = run_pca(matrix, args.components)
        if scores.shape[1] < 2:
            raise RuntimeError("Need at least two principal components for a PC1/PC2 free-energy landscape.")
        xg, yg, density = kde_grid(scores, 0, 1, args.bins)
        free_energy = free_energy_from_density(density, args.temperature)

        np.save(outdir / "pc_projections.npy", scores)
        np.save(outdir / "pca_components.npy", pca.components_)
        np.save(outdir / "pca_mean.npy", pca.mean_)
        np.save(outdir / "explained_variance.npy", pca.explained_variance_)
        np.save(outdir / "explained_variance_ratio.npy", pca.explained_variance_ratio_)
        np.save(outdir / "free_energy_grid.npy", free_energy)
        np.savetxt(outdir / "pc_scores.csv", scores, delimiter=",", header=",".join(f"PC{i + 1}" for i in range(scores.shape[1])), comments="")
        np.savetxt(outdir / "explained_variance_ratio.csv", pca.explained_variance_ratio_, delimiter=",")
        np.savetxt(outdir / "free_energy_grid.csv", free_energy, delimiter=",")

        save_plots(outdir, scores, pca, xg, yg, density, free_energy, args.temperature, args.bins)

        summary = {
            "status": "completed",
            "topology": str(topology_used),
            "requested_topology": str(Path(args.top)),
            "trajectory": str(Path(args.traj)),
            "selection": args.selection,
            "temperature_K": args.temperature,
            "total_trajectory_frames": total_frames,
            "frames_analyzed": int(scores.shape[0]),
            "selected_atoms": int(n_atoms),
            "components": int(scores.shape[1]),
            "explained_variance_ratio": [float(x) for x in pca.explained_variance_ratio_],
            "minimum_delta_g_kj_mol": float(np.nanmin(free_energy)),
            "maximum_delta_g_kj_mol": float(np.nanmax(free_energy)),
        }
        (outdir / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        (outdir / "analysis_summary.txt").write_text(
            "PCA Free Energy Analysis\n"
            "Status: completed\n"
            f"Topology: {topology_used}\n"
            f"Requested topology: {args.top}\n"
            f"Trajectory: {args.traj}\n"
            f"Selection: {args.selection}\n"
            f"Temperature: {args.temperature:g} K\n"
            f"Frames analyzed: {scores.shape[0]} of {total_frames}\n"
            f"Selected atoms: {n_atoms}\n"
            f"Explained variance ratio: {', '.join(f'{x:.4f}' for x in pca.explained_variance_ratio_)}\n",
            encoding="utf-8",
        )
    except Exception as exc:
        write_failure(outdir, str(exc))
        print(f"PCA free-energy analysis skipped: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
