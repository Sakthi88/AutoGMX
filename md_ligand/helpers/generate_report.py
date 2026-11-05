#!/usr/bin/env python3
"""Generate a comprehensive PDF report from MD simulation results.

Collects analysis plots, MM-PBSA results, simulation parameters,
and energy data into a single publication-quality PDF.
"""
import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_RIGHT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch, mm
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    Image,
    NextPageTemplate,
    PageBreak,
    PageTemplate,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


# ── Styles ─────────────────────────────────────────────────────────────────────
def _build_styles():
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(
        "ReportTitle", parent=styles["Title"],
        fontSize=22, spaceAfter=6, textColor=colors.HexColor("#1e293b"),
    ))
    styles.add(ParagraphStyle(
        "SectionHead", parent=styles["Heading1"],
        fontSize=16, spaceBefore=18, spaceAfter=8,
        textColor=colors.HexColor("#1e40af"),
        borderWidth=0, borderPadding=0,
    ))
    styles.add(ParagraphStyle(
        "SubSectionHead", parent=styles["Heading2"],
        fontSize=13, spaceBefore=12, spaceAfter=6,
        textColor=colors.HexColor("#334155"),
    ))
    styles.add(ParagraphStyle(
        "BodyText2", parent=styles["BodyText"],
        fontSize=10, leading=14, spaceAfter=6,
    ))
    styles.add(ParagraphStyle(
        "SmallText", parent=styles["BodyText"],
        fontSize=8, leading=10, textColor=colors.HexColor("#64748b"),
    ))
    styles.add(ParagraphStyle(
        "TableHeader", parent=styles["Normal"],
        fontSize=9, textColor=colors.white, alignment=TA_CENTER,
        fontName="Helvetica-Bold",
    ))
    styles.add(ParagraphStyle(
        "TableCell", parent=styles["Normal"],
        fontSize=9, alignment=TA_CENTER,
    ))
    return styles


# ── Table helpers ──────────────────────────────────────────────────────────────
def _param_table(params, styles):
    """Create a styled key-value table."""
    data = [[
        Paragraph("<b>Parameter</b>", styles["TableHeader"]),
        Paragraph("<b>Value</b>", styles["TableHeader"]),
    ]]
    for key, val in params:
        data.append([
            Paragraph(str(key), styles["TableCell"]),
            Paragraph(str(val), styles["TableCell"]),
        ])
    t = Table(data, colWidths=[180, 280])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1e40af")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#cbd5e1")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f1f5f9")]),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ]))
    return t


def _mmpbsa_table(summary, styles):
    """Create MM-PBSA results table."""
    components = summary.get("components", {})
    data = [[
        Paragraph("<b>Component</b>", styles["TableHeader"]),
        Paragraph("<b>Mean (kcal/mol)</b>", styles["TableHeader"]),
        Paragraph("<b>Std Dev</b>", styles["TableHeader"]),
        Paragraph("<b>SEM</b>", styles["TableHeader"]),
    ]]
    labels = {
        "vdw": "van der Waals",
        "elec": "Electrostatics",
        "polar_solv": "Polar Solvation",
        "nonpolar_solv": "Nonpolar Solvation",
        "delta_G": "Total ΔG Binding",
    }
    for key in ["vdw", "elec", "polar_solv", "nonpolar_solv", "delta_G"]:
        if key in components:
            c = components[key]
            data.append([
                Paragraph(labels.get(key, key), styles["TableCell"]),
                Paragraph(f"{c['mean']:.4f}", styles["TableCell"]),
                Paragraph(f"{c['std']:.4f}", styles["TableCell"]),
                Paragraph(f"{c['sem']:.4f}", styles["TableCell"]),
            ])
    t = Table(data, colWidths=[140, 120, 100, 100])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1e40af")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#cbd5e1")),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f1f5f9")]),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        # Highlight the total ΔG row
        ("BACKGROUND", (0, -1), (-1, -1), colors.HexColor("#ede9fe")),
        ("FONTNAME", (0, -1), (-1, -1), "Helvetica-Bold"),
    ]))
    return t


# ── Page numbering ─────────────────────────────────────────────────────────────
def _add_page_number(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 8)
    canvas.setFillColor(colors.HexColor("#64748b"))
    page_num = canvas.getPageNumber()
    text = f"Page {page_num}"
    canvas.drawCentredString(A4[0] / 2, 15 * mm, text)
    # Header line
    canvas.setStrokeColor(colors.HexColor("#e2e8f0"))
    canvas.setLineWidth(0.5)
    canvas.line(20 * mm, A4[1] - 18 * mm, A4[0] - 20 * mm, A4[1] - 18 * mm)
    canvas.restoreState()


# ── Report builder ─────────────────────────────────────────────────────────────
def build_report(args):
    """Build the PDF report from available analysis outputs."""
    analysis_dir = Path(args.analysis_dir)
    output_path = Path(args.output_pdf)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    styles = _build_styles()

    doc = SimpleDocTemplate(
        str(output_path),
        pagesize=A4,
        topMargin=25 * mm,
        bottomMargin=25 * mm,
        leftMargin=20 * mm,
        rightMargin=20 * mm,
    )

    story = []

    # ── Title page ────────────────────────────────────────────────────────
    story.append(Spacer(1, 60))
    story.append(Paragraph("Molecular Dynamics Simulation Report", styles["ReportTitle"]))
    story.append(Spacer(1, 12))

    project_name = args.project_name or "N/A"
    story.append(Paragraph(f"Project: <b>{project_name}</b>", styles["BodyText2"]))
    story.append(Paragraph(f"Generated: <b>{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</b>",
                           styles["BodyText2"]))
    story.append(Spacer(1, 20))

    # ── Section 1: Simulation Parameters ──────────────────────────────────
    story.append(Paragraph("1. Simulation Parameters", styles["SectionHead"]))

    # Read web_config.json if available
    config_data = {}
    config_path = analysis_dir.parent / "web_config.json"
    if config_path.exists():
        try:
            config_data = json.loads(config_path.read_text())
        except Exception:
            pass

    params = [
        ("Force Field", config_data.get("FORCE_FIELD", "N/A")),
        ("Water Model", config_data.get("WATER_MODEL", "N/A")),
        ("Ligand Prep Tool", config_data.get("LIGAND_PREP_TOOL", "N/A")),
        ("Ligand Residue", config_data.get("LIGAND_RESNAME", "N/A")),
        ("Production Time (ns)", config_data.get("MD_TIME_NS", "N/A")),
        ("Temperature (K)", "300"),
        ("Pressure (bar)", "1.0"),
        ("Timestep (ps)", "0.002"),
        ("Electrostatics", "PME"),
        ("VDW Method", "Force-switch"),
        ("Cutoff (nm)", "1.2"),
        ("TCoupl", "V-rescale"),
        ("PCoupl", "Parrinello-Rahman"),
        ("Constraints", "H-bonds (LINCS)"),
        ("GPU Acceleration", config_data.get("USE_GPU", "auto")),
        ("CPU Threads", config_data.get("NTOMP", "auto")),
    ]
    story.append(_param_table(params, styles))
    story.append(Spacer(1, 12))

    manifest_path = analysis_dir.parent / "results" / "production" / "ensemble_manifest.json"
    provenance = [
        ("Checkpoint policy", config_data.get("CHECKPOINT_POLICY", "validate")),
        ("Checkpoint interval (min)", config_data.get("CHECKPOINT_INTERVAL_MIN", 15)),
        ("Production replicas", config_data.get("ENSEMBLE_REPLICAS", 1)),
        ("Base random seed", config_data.get("ENSEMBLE_BASE_SEED", "N/A")),
    ]
    if manifest_path.exists():
        try:
            manifest = json.loads(manifest_path.read_text())
            provenance.extend([("Recorded replicas", manifest.get("replicas", "N/A")), ("GROMACS version", manifest.get("gromacs_version", "unknown")), ("Ensemble method", manifest.get("mode", "N/A"))])
        except (OSError, json.JSONDecodeError):
            provenance.append(("Ensemble manifest", "Unreadable"))
    story.append(Paragraph("2. Reproducibility and Restart Record", styles["SectionHead"]))
    story.append(Paragraph("Deterministic replica seeds, periodic checkpoints, and checkpoint validation make this run restartable and auditable.", styles["BodyText2"]))
    story.append(_param_table(provenance, styles))
    story.append(Spacer(1, 12))

    # ── Section 2: Analysis Plots ─────────────────────────────────────────
    story.append(Paragraph("2. Structural Analysis", styles["SectionHead"]))

    plot_configs = [
        ("rmsd_protein.png", "RMSD - Protein"),
        ("rmsd_ligand.png", "RMSD - Ligand"),
        ("rmsd_protein_ligand.png", "RMSD Comparison (Protein vs Ligand)"),
        ("rmsf_protein.png", "RMSF - Protein"),
        ("gyration_protein.png", "Radius of Gyration"),
        ("hbonds_protein_ligand.png", "Hydrogen Bonds (Protein-Ligand)"),
    ]

    for filename, title in plot_configs:
        plot_path = analysis_dir / filename
        if plot_path.exists():
            story.append(Paragraph(f"2.{plot_configs.index((filename, title)) + 1} {title}",
                                   styles["SubSectionHead"]))
            try:
                img = Image(str(plot_path), width=460, height=290)
                img.hAlign = "CENTER"
                story.append(img)
            except Exception:
                story.append(Paragraph(f"[Plot file exists but could not be embedded: {filename}]",
                                       styles["SmallText"]))
            story.append(Spacer(1, 8))

    # ── Section 3: Energetics ─────────────────────────────────────────────
    story.append(Paragraph("3. Energetics", styles["SectionHead"]))

    energy_plots = [
        ("potential.png", "Potential Energy"),
        ("temperature.png", "Temperature"),
        ("pressure.png", "Pressure"),
        ("density.png", "Density"),
    ]

    for filename, title in energy_plots:
        plot_path = analysis_dir / filename
        if plot_path.exists():
            story.append(Paragraph(f"3.{energy_plots.index((filename, title)) + 1} {title}",
                                   styles["SubSectionHead"]))
            try:
                img = Image(str(plot_path), width=460, height=290)
                img.hAlign = "CENTER"
                story.append(img)
            except Exception:
                story.append(Paragraph(f"[Plot file exists but could not be embedded: {filename}]",
                                       styles["SmallText"]))
            story.append(Spacer(1, 8))

    # ── Section 4: MM-PBSA ────────────────────────────────────────────────
    mmpbsa_json = analysis_dir / "mmpbsa" / "mmpbsa_summary.json"
    if mmpbsa_json.exists():
        story.append(PageBreak())
        story.append(Paragraph("4. MM-PBSA Binding Free Energy", styles["SectionHead"]))

        try:
            mmpbsa_data = json.loads(mmpbsa_json.read_text())
        except Exception:
            mmpbsa_data = {}

        # MM-PBSA parameters
        mmpbsa_params = [
            ("Energy Method", mmpbsa_data.get("energy_method", "N/A")),
            ("PB Method", mmpbsa_data.get("pb_method", "N/A")),
            ("Frames Analyzed", mmpbsa_data.get("n_frames", "N/A")),
            ("Frame Step", mmpbsa_data.get("frame_step", "N/A")),
            ("Temperature (K)", mmpbsa_data.get("temperature", "N/A")),
            ("pH", mmpbsa_data.get("pH", "N/A")),
            ("Ionic Strength (M)", mmpbsa_data.get("ionic_strength", "N/A")),
            ("Surface Tension", mmpbsa_data.get("surface_tension", "N/A")),
            ("Receptor Selection", mmpbsa_data.get("selections", {}).get("receptor", "protein")),
            ("Ligand Selection", mmpbsa_data.get("selections", {}).get("ligand", "resname LIG")),
        ]
        story.append(Paragraph("4.1 MM-PBSA Parameters", styles["SubSectionHead"]))
        story.append(_param_table(mmpbsa_params, styles))
        story.append(Spacer(1, 12))

        # Results table
        story.append(Paragraph("4.2 Binding Free Energy Decomposition", styles["SubSectionHead"]))
        story.append(_mmpbsa_table(mmpbsa_data, styles))
        story.append(Spacer(1, 12))

        # Plots
        mmpbsa_plots = [
            ("mmpbsa_decomposition.png", "4.3 Energy Component Bar Chart"),
            ("mmpbsa_timeseries.png", "4.4 Energy Components Over Time"),
        ]
        for filename, title in mmpbsa_plots:
            plot_path = analysis_dir / "mmpbsa" / filename
            if plot_path.exists():
                story.append(Paragraph(title, styles["SubSectionHead"]))
                try:
                    img = Image(str(plot_path), width=460, height=290)
                    img.hAlign = "CENTER"
                    story.append(img)
                except Exception:
                    story.append(Paragraph(f"[Plot: {filename}]", styles["SmallText"]))
                story.append(Spacer(1, 8))

        # Per-residue decomposition
        decomp_plot = analysis_dir / "mmpbsa" / "decomposition_plot.png"
        if decomp_plot.exists():
            story.append(Paragraph("4.5 Per-Residue Energy Decomposition",
                                   styles["SubSectionHead"]))
            try:
                img = Image(str(decomp_plot), width=460, height=250)
                img.hAlign = "CENTER"
                story.append(img)
            except Exception:
                story.append(Paragraph("[Decomposition plot]", styles["SmallText"]))
            story.append(Spacer(1, 8))

            # Top contributing residues table
            decomp_summary = analysis_dir / "mmpbsa" / "decomposition_summary.json"
            if decomp_summary.exists():
                try:
                    decomp_data = json.loads(decomp_summary.read_text())
                    story.append(Paragraph("4.6 Top Contributing Residues",
                                           styles["SubSectionHead"]))
                    decomp_table_data = [[
                        Paragraph("<b>Residue</b>", styles["TableHeader"]),
                        Paragraph("<b>vdW (kcal/mol)</b>", styles["TableHeader"]),
                        Paragraph("<b>Elec (kcal/mol)</b>", styles["TableHeader"]),
                        Paragraph("<b>Total (kcal/mol)</b>", styles["TableHeader"]),
                    ]]
                    for res in decomp_data.get("top_binding", [])[:10]:
                        decomp_table_data.append([
                            Paragraph(res.get("residue", ""), styles["TableCell"]),
                            Paragraph(f"{res.get('vdw_mean', 0):.4f}", styles["TableCell"]),
                            Paragraph(f"{res.get('elec_mean', 0):.4f}", styles["TableCell"]),
                            Paragraph(f"{res.get('total_mean', 0):.4f}", styles["TableCell"]),
                        ])
                    dt = Table(decomp_table_data, colWidths=[100, 120, 120, 120])
                    dt.setStyle(TableStyle([
                        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#1e40af")),
                        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                        ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                        ("GRID", (0, 0), (-1, -1), 0.5, colors.HexColor("#cbd5e1")),
                        ("ROWBACKGROUNDS", (0, 1), (-1, -1),
                         [colors.white, colors.HexColor("#f1f5f9")]),
                        ("TOPPADDING", (0, 0), (-1, -1), 3),
                        ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
                    ]))
                    story.append(dt)
                except Exception:
                    pass

    # ── Section 5: Free Energy Landscape ──────────────────────────────────
    fe_dir = analysis_dir / "free_energy"
    fe_plots = list(fe_dir.glob("*.png")) if fe_dir.exists() else []
    if fe_plots:
        story.append(PageBreak())
        story.append(Paragraph("5. Free Energy Landscape", styles["SectionHead"]))
        fe_summary = fe_dir / "analysis_summary.json"
        if fe_summary.exists():
            try:
                fe_data = json.loads(fe_summary.read_text())
                fe_params = [
                    ("Selection", fe_data.get("selection", "N/A")),
                    ("Temperature (K)", fe_data.get("temperature", "N/A")),
                    ("PCA Components", fe_data.get("n_components", "N/A")),
                    ("Grid Bins", fe_data.get("bins", "N/A")),
                ]
                story.append(_param_table(fe_params, styles))
                story.append(Spacer(1, 8))
            except Exception:
                pass

        for i, plot_file in enumerate(sorted(fe_plots)):
            if plot_file.suffix == ".png":
                try:
                    img = Image(str(plot_file), width=440, height=300)
                    img.hAlign = "CENTER"
                    story.append(img)
                    story.append(Spacer(1, 8))
                except Exception:
                    pass

    # ── Section 6: Conclusion ─────────────────────────────────────────────
    story.append(PageBreak())
    story.append(Paragraph("6. Summary", styles["SectionHead"]))
    story.append(Paragraph(
        "This report summarizes the molecular dynamics simulation results generated by the "
        "automated MD pipeline. All analyses were performed using GROMACS tools and custom "
        "Python scripts. The MM-PBSA binding free energy was calculated using APBS for "
        "Poisson-Boltzmann electrostatics and MDAnalysis for trajectory processing.",
        styles["BodyText2"],
    ))

    # Extract key result
    if mmpbsa_json.exists():
        try:
            mmpbsa_data = json.loads(mmpbsa_json.read_text())
            total_dg = mmpbsa_data.get("components", {}).get("delta_G", {})
            if total_dg:
                story.append(Spacer(1, 12))
                story.append(Paragraph(
                    f"<b>Estimated Binding Free Energy (ΔG): "
                    f"{total_dg['mean']:.2f} ± {total_dg['std']:.2f} kcal/mol</b>",
                    styles["BodyText2"],
                ))
        except Exception:
            pass

    story.append(Spacer(1, 20))
    story.append(Paragraph(
        "For questions about this report or the analysis methodology, please refer to the "
        "pipeline documentation or contact the project administrator.",
        styles["SmallText"],
    ))

    # ── Build PDF ─────────────────────────────────────────────────────────
    doc.build(story, onFirstPage=_add_page_number, onLaterPages=_add_page_number)
    print(f"PDF report generated: {output_path}")
    return str(output_path)


# ── CLI ────────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="Generate MD simulation PDF report")
    p.add_argument("-a", "--analysis-dir", required=True, help="Analysis output directory")
    p.add_argument("-o", "--output-pdf", required=True, help="Output PDF path")
    p.add_argument("-n", "--project-name", default=None, help="Project name for title")
    return p.parse_args()


if __name__ == "__main__":
    args = parse_args()
    build_report(args)
