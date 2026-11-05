#!/usr/bin/env python3
"""Unified ligand topology generator with backend abstraction.

This module provides a single entry point for generating ligand topologies
using various force fields with proper backend separation.

Architecture:
- Each force field family has dedicated backends
- Backends handle parameterization (not just file conversion)
- Mandatory validation: stereochemistry, charges, gmx grompp
- Provenance tracking for reproducibility
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path

# Backend registry - maps tool names to backend classes
BACKEND_REGISTRY = {}


class TopologyBackend(ABC):
    """Abstract base class for topology generation backends."""
    
    name: str = None
    force_field_family: str = None
    supported_force_fields: list = []
    
    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if cls.name:
            BACKEND_REGISTRY[cls.name] = cls
    
    @abstractmethod
    def check_dependencies(self) -> list:
        """Return list of missing dependencies."""
        pass
    
    @abstractmethod
    def parameterize(self, smiles: str, name: str, work_dir: Path, 
                     force_field: str, net_charge: int) -> tuple[Path, Path]:
        """Generate parameterized ITP and GRO files.
        
        Returns:
            Tuple of (itp_path, gro_path)
        """
        pass
    
    def validate_output(self, itp_path: Path, gro_path: Path, 
                        formal_charge: int) -> tuple[bool, str]:
        """Validate generated topology. Override for backend-specific validation."""
        return True, "Basic validation passed"
    
    def get_force_field_include(self, force_field: str) -> str:
        """Get the GROMACS force field include path."""
        return f"{force_field}.ff/forcefield.itp"


class AmberBackend(TopologyBackend):
    """AMBER GAFF/GAFF2 via ACPYPE."""
    name = "acpype"
    force_field_family = "AMBER"
    supported_force_fields = ["gaff", "gaff2"]
    
    def check_dependencies(self) -> list:
        missing = []
        if not shutil.which("acpype"):
            missing.append("acpype (conda install -c conda-forge acpype)")
        if not shutil.which("obabel"):
            missing.append("openbabel (conda install -c conda-forge openbabel)")
        return missing
    
    def parameterize(self, smiles: str, name: str, work_dir: Path, 
                     force_field: str, net_charge: int) -> tuple[Path, Path]:
        pipedir = Path(__file__).parent.parent
        atom_type = force_field if force_field in ["gaff", "gaff2"] else "gaff2"
        
        # Write SMILES to MOL2 for ACPYPE
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        mol2_path = work_dir / f"{name}.mol2"
        writer = Chem.Mol2Writer(str(mol2_path))
        writer.write(mol)
        writer.close()
        
        cmd = [
            sys.executable, str(pipedir / "helpers" / "run_acpype.py"),
            "-i", str(mol2_path),
            "-b", name,
            "-a", atom_type,
            "-n", str(net_charge),
        ]
        
        result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"ACPYPE failed: {result.stderr}")
        
        acpype_dir = work_dir / f"{name}.acpype"
        itp_path = acpype_dir / f"{name}_GMX.itp"
        gro_path = acpype_dir / f"{name}_GMX.gro"
        
        if not itp_path.exists() or not gro_path.exists():
            raise RuntimeError("ACPYPE did not produce expected output files")
        
        return itp_path, gro_path


class CharmmBackend(TopologyBackend):
    """CHARMM/CGenFF via cgenff_charmm2gmx."""
    name = "cgenff"
    force_field_family = "CHARMM"
    supported_force_fields = ["charmm36", "cgenff"]
    
    def check_dependencies(self) -> list:
        missing = []
        pipedir = Path(__file__).parent.parent
        if not shutil.which(os.environ.get("CGENFF_BIN", "cgenff")):
            missing.append("cgenff binary (from ParamChem)")
        if not shutil.which(os.environ.get("CGENFF_PYTHON", "python")):
            missing.append("python for cgenff converter")
        if not Path(os.environ.get("CGENFF_CONVERTER", str(pipedir / "cgenff_charmm2gmx_py3_nx3.py"))).exists():
            missing.append("cgenff_charmm2gmx converter script")
        return missing
    
    def parameterize(self, smiles: str, name: str, work_dir: Path, 
                     force_field: str, net_charge: int) -> tuple[Path, Path]:
        pipedir = Path(__file__).parent.parent
        
        # Write SMILES to MOL2
        from rdkit import Chem
        from rdkit.Chem import AllChem
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        
        mol2_path = work_dir / f"{name}.mol2"
        writer = Chem.Mol2Writer(str(mol2_path))
        writer.write(mol)
        writer.close()
        
        cgenff_bin = os.environ.get("CGENFF_BIN", "cgenff")
        cgenff_python = os.environ.get("CGENFF_PYTHON", "python")
        cgenff_converter = os.environ.get("CGENFF_CONVERTER", str(pipedir / "cgenff_charmm2gmx_py3_nx3.py"))
        cgenff_ff_dir = os.environ.get("CGENFF_FORCEFIELD_DIR", "")
        
        # Run cgenff
        str_file = work_dir / f"{name}.str"
        result = subprocess.run([cgenff_bin, str(mol2_path)], cwd=work_dir, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"CGenFF failed: {result.stderr}")
        str_file.write_text(result.stdout)
        
        # Run converter
        cmd = [cgenff_python, cgenff_converter, name, str(mol2_path), str(str_file), cgenff_ff_dir]
        result = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"CGenFF converter failed: {result.stderr}")
        
        base_name = name.lower()
        itp_path = work_dir / f"{base_name}.itp"
        gro_path = work_dir / f"{base_name}_ini.pdb"
        
        if not itp_path.exists():
            raise RuntimeError("CGenFF did not produce ITP file")
        
        # Convert PDB to GRO
        final_gro = work_dir / f"{name}.gro"
        subprocess.run(["gmx", "editconf", "-f", str(gro_path), "-o", str(final_gro)], check=True, capture_output=True)
        
        return itp_path, final_gro


class GromosBackend(TopologyBackend):
    """GROMOS via GroLigFF/ATB (true parameterization)."""
    name = "gromos"
    force_field_family = "GROMOS"
    supported_force_fields = ["gromos54a7", "gromos2016h66"]
    
    def check_dependencies(self) -> list:
        missing = []
        # GroLigFF is the primary backend
        if not shutil.which("groligff"):
            missing.append("GroLigFF (groligff binary from GROMOS) - required for novel ligand parameterization")
        return missing
    
    def parameterize(self, smiles: str, name: str, work_dir: Path, 
                     force_field: str, net_charge: int) -> tuple[Path, Path]:
        # Import and use the GroLigFF parameterization from run_gromos
        import run_gromos
        return run_gromos.run_groligff_parameterization(smiles, name, work_dir, force_field, net_charge)
    
    def get_force_field_include(self, force_field: str) -> str:
        ff_map = {"gromos54a7": "gromos54a7", "gromos2016h66": "gromos2016h66"}
        return f"{ff_map.get(force_field, force_field)}.ff/forcefield.itp"


class OplsBackend(TopologyBackend):
    """OPLS-AA via LigParGen (true OPLS-AA parameterization)."""
    name = "oplsaa"
    force_field_family = "OPLS-AA"
    supported_force_fields = ["oplsaa"]
    
    def check_dependencies(self) -> list:
        missing = []
        try:
            import requests
        except ImportError:
            missing.append("requests (pip install requests)")
        if not shutil.which("obabel"):
            missing.append("openbabel (conda install -c conda-forge openbabel)")
        return missing
    
    def parameterize(self, smiles: str, name: str, work_dir: Path, 
                     force_field: str, net_charge: int) -> tuple[Path, Path]:
        import run_oplsaa
        return run_oplsaa.run_ligpargen_parameterization(smiles, name, work_dir, net_charge)


class OpenFFBackend(TopologyBackend):
    """OpenFF Sage (NOT OPLS-AA - separate force field family)."""
    name = "openff"
    force_field_family = "OpenFF"
    supported_force_fields = ["openff-2.0.0", "openff-2.1.0"]
    
    def check_dependencies(self) -> list:
        missing = []
        try:
            from openff.toolkit import Molecule
            from openff.toolkit.typing.engines.smirnoff import ForceField
            from openff.interchange import Interchange
        except ImportError:
            missing.append("openff-toolkit, openff-interchange (conda install -c conda-forge openff-toolkit openff-interchange)")
        return missing
    
    def parameterize(self, smiles: str, name: str, work_dir: Path, 
                     force_field: str, net_charge: int) -> tuple[Path, Path]:
        import run_oplsaa
        return run_oplsaa.run_openff_parameterization(smiles, name, work_dir, force_field, net_charge)
    
    def get_force_field_include(self, force_field: str) -> str:
        return "openff.ff/forcefield.itp"


class LigParGenBackend(TopologyBackend):
    """LigParGen web service for OPLS-AA (alias for oplsaa backend)."""
    name = "ligpargen"
    force_field_family = "OPLS-AA"
    supported_force_fields = ["oplsaa"]
    
    def check_dependencies(self) -> list:
        missing = []
        try:
            import requests
        except ImportError:
            missing.append("requests (pip install requests)")
        if not shutil.which("obabel"):
            missing.append("openbabel (conda install -c conda-forge openbabel)")
        return missing
    
    def parameterize(self, smiles: str, name: str, work_dir: Path, 
                     force_field: str, net_charge: int) -> tuple[Path, Path]:
        import run_oplsaa
        return run_oplsaa.run_ligpargen_parameterization(smiles, name, work_dir, net_charge)


def extract_partial_charges_from_itp(itp_path: Path) -> list:
    """Extract partial charges from ITP [ atoms ] section."""
    charges = []
    in_atoms = False
    with open(itp_path, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("[ atoms ]"):
                in_atoms = True
                continue
            if in_atoms and stripped.startswith("["):
                break
            if in_atoms and stripped and not stripped.startswith(";"):
                parts = stripped.split()
                if len(parts) >= 7 and parts[0].isdigit():
                    try:
                        charges.append(float(parts[6]))
                    except (ValueError, IndexError):
                        pass
    return charges


def validate_charges(partial_charges: list, formal_charge: int, tolerance: float = 0.01) -> tuple[bool, str]:
    """Validate that partial charges sum to formal charge."""
    total = sum(partial_charges)
    diff = abs(total - formal_charge)
    if diff > tolerance:
        return False, f"Partial charges sum to {total:.4f}, expected formal charge {formal_charge} (diff={diff:.4f})"
    return True, f"Charge validated: sum={total:.4f}, formal={formal_charge}"


def extract_moleculetype(itp_path: Path) -> str:
    """Extract molecule type name from ITP file."""
    with open(itp_path, "r") as f:
        in_moleculetype = False
        for line in f:
            stripped = line.strip()
            if stripped.startswith("[ moleculetype ]"):
                in_moleculetype = True
                continue
            if in_moleculetype and stripped and not stripped.startswith(";"):
                return stripped.split()[0]
    raise RuntimeError(f"Could not determine [ moleculetype ] from {itp_path}")


def write_ligand_top(top_path: Path, itp_path: Path, molecule_name: str, force_field: str, backend: TopologyBackend):
    """Write minimal topology file including the ligand ITP."""
    ff_include = backend.get_force_field_include(force_field)
    
    lines = [
        f'#include "{ff_include}"',
        f'#include "{itp_path.name}"',
        "",
        "[ system ]",
        f"{molecule_name} in vacuum",
        "",
        "[ molecules ]",
        f"{molecule_name} 1",
        "",
    ]
    top_path.write_text("\n".join(lines), encoding="utf-8", newline="\n")


def count_itp_atoms(itp_path: Path) -> int:
    """Count atoms in ITP [ atoms ] section."""
    count = 0
    in_atoms = False
    with open(itp_path, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("[ atoms ]"):
                in_atoms = True
                continue
            if in_atoms and stripped.startswith("["):
                break
            if in_atoms and stripped and not stripped.startswith(";"):
                parts = stripped.split()
                if parts and parts[0].isdigit():
                    count += 1
    return count


def count_gro_atoms(gro_path: Path) -> int:
    """Count atoms in GRO file."""
    with open(gro_path, "r") as f:
        lines = f.readlines()
    if len(lines) >= 2:
        try:
            return int(lines[1].strip())
        except ValueError:
            pass
    return 0


def validate_stereochemistry(mol: "Chem.Mol", allow_undefined: bool = False) -> list:
    """Validate stereochemistry, return list of warnings/errors."""
    issues = []
    from rdkit import Chem
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    for atom_idx, stereo in chiral_centers:
        if stereo == "?":
            msg = f"Undefined stereochemistry at atom {atom_idx}"
            if allow_undefined:
                issues.append(("warning", msg))
            else:
                issues.append(("error", msg))
    
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            stereo = bond.GetStereo()
            if stereo == Chem.BondStereo.STEREONONE:
                a1 = bond.GetBeginAtom()
                a2 = bond.GetEndAtom()
                if a1.GetDegree() > 1 and a2.GetDegree() > 1:
                    msg = f"Undefined double bond stereochemistry at bond {bond.GetIdx()}"
                    if allow_undefined:
                        issues.append(("warning", msg))
                    else:
                        issues.append(("error", msg))
    
    return issues


def convert_to_smiles(input_path: Path) -> str:
    """Convert MOL2/PDB/SDF file to SMILES using RDKit."""
    from rdkit import Chem
    mol = None
    suffix = input_path.suffix.lower()
    if suffix == ".mol2":
        mol = Chem.MolFromMol2File(str(input_path), removeHs=False)
    elif suffix == ".pdb":
        mol = Chem.MolFromPDBFile(str(input_path), removeHs=False)
    elif suffix in [".sdf", ".mol"]:
        mol = Chem.MolFromMolFile(str(input_path), removeHs=False)
    
    if mol is None:
        raise RuntimeError(f"Could not parse {input_path} with RDKit")
    
    return Chem.MolToSmiles(mol, canonical=True)


def validate_topology_with_grompp(top_path: Path, gro_path: Path, work_dir: Path, gmx_bin: str = "gmx") -> tuple[bool, str]:
    """Validate topology with gmx grompp."""
    mdp_path = work_dir / "validate.mdp"
    mdp_path.write_text("""
integrator  = steep
nsteps      = 0
emtol       = 1000
emstep      = 0.01
nstlist     = 1
cutoff-scheme = Verlet
ns_type     = grid
coulombtype = PME
rcoulomb    = 1.0
rvdw        = 1.0
pbc         = xyz
""", encoding="utf-8")
    
    tpr_path = work_dir / "validate.tpr"
    cmd = [gmx_bin, "grompp", "-f", str(mdp_path), "-c", str(gro_path), "-p", str(top_path), "-o", str(tpr_path), "-maxwarn", "10"]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=work_dir)
    
    if result.returncode == 0:
        return True, "gmx grompp validation PASSED"
    else:
        return False, f"gmx grompp validation FAILED:\n{result.stderr}\n{result.stdout}"


def write_provenance(output_dir: Path, metadata: dict):
    """Write parameterization provenance file."""
    prov_path = output_dir / f"{metadata['ligand']}.parameterization.json"
    metadata["timestamp"] = datetime.now().isoformat()
    prov_path.write_text(json.dumps(metadata, indent=2))


def list_available_backends():
    """Print available backends and their force fields."""
    print("Available backends:")
    for name, backend_class in BACKEND_REGISTRY.items():
        print(f"  {name}: {backend_class.force_field_family}")
        for ff in backend_class.supported_force_fields:
            print(f"    - {ff}")
    print()


def main():
    parser = argparse.ArgumentParser(
        description="Generate ligand topology for various force fields using proper backends",
        add_help=False
    )
    parser.add_argument("-h", "--help", action="help", help="Show this help message and exit")
    parser.add_argument("--list-backends", action="store_true", help="List available backends and exit")
    parser.add_argument("--backend", choices=list(BACKEND_REGISTRY.keys()),
                        help="Parameterization backend")
    parser.add_argument("--input", help="Input file (MOL2, PDB, SDF, or SMILES string)")
    parser.add_argument("--name", help="Ligand name (residue name, max 5 chars)")
    parser.add_argument("--force-field", help="Force field variant (backend-specific)")
    parser.add_argument("--charge", type=int, default=0, help="Net formal charge of ligand")
    parser.add_argument("--output-itp", help="Output ITP file")
    parser.add_argument("--output-gro", help="Output GRO file")
    parser.add_argument("--output-top", help="Output TOP file")
    parser.add_argument("--work-dir", help="Working directory")
    parser.add_argument("--smiles", action="store_true", help="Input is SMILES string")
    parser.add_argument("--allow-undefined-stereo", action="store_true", 
                        help="Allow undefined stereochemistry (not recommended)")
    parser.add_argument("--skip-validation", action="store_true", help="Skip gmx grompp validation")
    parser.add_argument("--gmx-bin", default="gmx", help="GROMACS binary")
    
    args = parser.parse_args()
    
    if args.list_backends:
        list_available_backends()
        return
    
    # Validate required arguments
    required = ["backend", "input", "name", "force_field", "output_itp", "output_gro", "output_top", "work_dir"]
    missing = [arg for arg in required if getattr(args, arg.replace("-", "_")) is None]
    if missing:
        parser.error(f"the following arguments are required: {', '.join('--' + m.replace('_', '-') for m in missing)}")
    
    # Get backend instance
    backend_class = BACKEND_REGISTRY[args.backend]
    backend = backend_class()
    
    # Validate force field is supported by backend
    if args.force_field not in backend.supported_force_fields:
        print(f"ERROR: Backend '{args.backend}' does not support force field '{args.force_field}'", file=sys.stderr)
        print(f"Supported: {backend.supported_force_fields}", file=sys.stderr)
        sys.exit(1)
    
    # Check dependencies
    missing = backend.check_dependencies()
    if missing:
        print(f"Missing dependencies for {args.backend}: {', '.join(missing)}", file=sys.stderr)
        sys.exit(1)
    
    input_path = Path(args.input)
    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    
    output_itp = Path(args.output_itp)
    output_gro = Path(args.output_gro)
    output_top = Path(args.output_top)
    
    # Get SMILES
    if args.smiles:
        smiles = args.input
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise RuntimeError(f"Invalid SMILES: {smiles}")
    else:
        if not input_path.exists():
            print(f"Input file not found: {input_path}", file=sys.stderr)
            sys.exit(1)
        print(f"Converting {input_path} to SMILES...")
        smiles = convert_to_smiles(input_path)
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles)
        print(f"SMILES: {smiles}")
    
    # Validate stereochemistry
    stereo_issues = validate_stereochemistry(mol, allow_undefined=args.allow_undefined_stereo)
    errors = [msg for level, msg in stereo_issues if level == "error"]
    warnings = [msg for level, msg in stereo_issues if level == "warning"]
    
    if errors:
        print("STEREOCHEMISTRY ERRORS:", file=sys.stderr)
        for err in errors:
            print(f"  - {err}", file=sys.stderr)
        print("Use --allow-undefined-stereo to override (not recommended)", file=sys.stderr)
        sys.exit(1)
    
    for warn in warnings:
        print(f"WARNING: {warn}")
    
    print(f"Generating {backend.force_field_family} {args.force_field} topology for {args.name} using {args.backend}...")
    
    try:
        # Parameterize
        final_itp, final_gro = backend.parameterize(
            smiles, args.name, work_dir, args.force_field, args.charge
        )
        
        # Validate atom count match
        itp_atoms = count_itp_atoms(final_itp)
        gro_atoms = count_gro_atoms(final_gro)
        if itp_atoms != gro_atoms:
            raise RuntimeError(f"Atom count mismatch: ITP has {itp_atoms}, GRO has {gro_atoms}")
        
        # Validate charges
        partial_charges = extract_partial_charges_from_itp(final_itp)
        if partial_charges:
            charge_valid, charge_msg = validate_charges(partial_charges, args.charge)
            print(f"Charge validation: {charge_msg}")
            if not charge_valid:
                print(f"ERROR: {charge_msg}", file=sys.stderr)
                sys.exit(1)
        else:
            print("WARNING: Could not extract partial charges from ITP")
        
        # Backend-specific validation
        valid, msg = backend.validate_output(final_itp, final_gro, args.charge)
        if not valid:
            print(f"Backend validation failed: {msg}", file=sys.stderr)
            sys.exit(1)
        print(f"Backend validation: {msg}")
        
        # Copy to final locations
        shutil.copy2(final_itp, output_itp)
        shutil.copy2(final_gro, output_gro)
        
        # Extract molecule name and write topology
        molecule_name = extract_moleculetype(output_itp)
        write_ligand_top(output_top, output_itp, molecule_name, args.force_field, backend)
        
        # Validate with gmx grompp
        if not args.skip_validation:
            print("Validating topology with gmx grompp...")
            valid, msg = validate_topology_with_grompp(output_top, output_gro, work_dir, args.gmx_bin)
            print(msg)
            if not valid:
                sys.exit(1)
        
        # Write provenance
        provenance = {
            "ligand": args.name,
            "input": str(input_path) if not args.smiles else "SMILES",
            "smiles": smiles,
            "force_field": args.force_field,
            "force_field_family": backend.force_field_family,
            "backend": args.backend,
            "formal_charge": args.charge,
            "partial_charges_sum": sum(partial_charges) if partial_charges else None,
            "atom_count": itp_atoms,
            "parameterization_status": "validated" if not args.skip_validation else "not_validated",
            "stereochemistry_warnings": warnings,
            "gromacs_version": subprocess.run([args.gmx_bin, "--version"], capture_output=True, text=True).stdout.split()[2] if shutil.which(args.gmx_bin) else "unknown",
            "rdkit_version": Chem.__version__,
        }
        write_provenance(work_dir, provenance)
        
        print(f"Generated: {output_itp}, {output_gro}, {output_top}")
        print(f"Provenance: {work_dir / f'{args.name}.parameterization.json'}")
        
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()