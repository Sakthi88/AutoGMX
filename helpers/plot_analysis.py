#!/usr/bin/env python3
import pathlib
import sys

try:
    import matplotlib.pyplot as plt
except Exception as exc:
    raise SystemExit(f"matplotlib is required for plotting: {exc}")


def read_xvg(path):
    x_vals = []
    y_vals = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line[0] in {"@", "#"}:
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    x_vals.append(float(parts[0]))
                    y_vals.append(float(parts[1]))
                except ValueError:
                    continue
    return x_vals, y_vals


def plot_rmsd_comparison(protein_path, ligand_path, output_path):
    protein_x, protein_y = read_xvg(protein_path)
    ligand_x, ligand_y = read_xvg(ligand_path)
    if not protein_x or not protein_y:
        raise SystemExit(f"No plottable protein RMSD data found in {protein_path}")
    if not ligand_x or not ligand_y:
        raise SystemExit(f"No plottable ligand RMSD data found in {ligand_path}")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(protein_x, protein_y, color="#2563eb", linewidth=1.6, label="Protein RMSD")
    ax.plot(ligand_x, ligand_y, color="#ea580c", linewidth=1.6, label="Ligand RMSD")
    ax.set_title("Protein and Ligand RMSD")
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("RMSD (nm)")
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def main():
    if len(sys.argv) == 5 and sys.argv[1] == "--rmsd-comparison":
        plot_rmsd_comparison(
            pathlib.Path(sys.argv[2]),
            pathlib.Path(sys.argv[3]),
            pathlib.Path(sys.argv[4]),
        )
        return

    if len(sys.argv) != 4:
        raise SystemExit(
            "Usage: plot_analysis.py input.xvg title output.png\n"
            "   or: plot_analysis.py --rmsd-comparison protein.xvg ligand.xvg output.png"
        )
    input_path = pathlib.Path(sys.argv[1])
    title = sys.argv[2]
    output_path = pathlib.Path(sys.argv[3])

    x_vals, y_vals = read_xvg(input_path)
    if not x_vals or not y_vals:
        raise SystemExit(f"No plottable numeric data found in {input_path}")
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x_vals, y_vals, linewidth=1.5)
    ax.set_title(title)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
