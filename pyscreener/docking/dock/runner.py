import os
from pathlib import Path
import re
import subprocess as sp
import sys
from typing import Mapping, Optional, Tuple, Union

from openbabel import pybel
from rdkit.Chem import AllChem as Chem
import ray

from pyscreener.exceptions import MisconfiguredDirectoryError, MissingEnvironmentVariableError
from pyscreener.utils import calc_score
from pyscreener.docking import CalculationData, DockingRunner, Result
from pyscreener.docking.dock import utils
from pyscreener.docking.dock.metadata import DOCKMetadata

try:
    DOCK6 = Path(os.environ["DOCK6"])
except KeyError:
    raise MissingEnvironmentVariableError(
        "DOCK6 environment variable not set! "
        "See https://github.com/coleygroup/pyscreener#specifying-an-environment-variable "
        "for more information."
    )

VDW_DEFN_FILE = DOCK6 / "parameters" / "vdw_AMBER_parm99.defn"
FLEX_DEFN_FILE = DOCK6 / "parameters" / "flex.defn"
FLEX_DRIVE_FILE = DOCK6 / "parameters" / "flex_drive.tbl"
DOCK = DOCK6 / "bin" / "dock6"

for f in (VDW_DEFN_FILE, FLEX_DEFN_FILE, FLEX_DRIVE_FILE, DOCK):
    if not f.exists():
        raise MisconfiguredDirectoryError(
            f'$DOCK6 directory not configured properly! DOCK6 path is set as "{DOCK6}", but there '
            f'is no "{f.name}" located under the "{f.parents}" directory. '
            "See https://github.com/coleygroup/pyscreener#specifying-an-environment-variable for more information."
        )


class DOCKRunner(DockingRunner):
    @staticmethod
    def prepare(data: CalculationData) -> CalculationData:
        data = DOCKRunner.prepare_receptor(data)
        data = DOCKRunner.prepare_ligand(data)

        return data

    @staticmethod
    def prepare_receptor(data: CalculationData) -> CalculationData:
        """Prepare the files necessary to dock ligands against the input receptor using this
        Screener's parameters

        Parameter
        ---------
        receptor : str
            the filepath of a file containing a receptor. Must be in a format that is readable by
            Chimera

        Returns
        -------
        rec_sph : str
            the filepath of the file containing the selected spheres
        grid_prefix : str
            the prefix of all prepared grid files.
        None
            if receptor preparation fails at any point
        """
        rec_mol2 = utils.prepare_mol2(data.receptor, data.in_path)
        rec_pdb = utils.prepare_pdb(data.receptor, data.in_path)
        
        if rec_mol2 is None or rec_pdb is None:
            return data

        rec_dms = utils.prepare_dms(rec_pdb, data.metadata.probe_radius, data.in_path)
        if rec_dms is None:
            return data

        rec_sph = utils.prepare_sph(
            rec_dms,
            data.metadata.steric_clash_dist,
            data.metadata.min_radius,
            data.metadata.max_radius,
            data.in_path,
        )
        if rec_sph is None:
            return data

        rec_sph = utils.select_spheres(
            rec_sph,
            data.metadata.sphere_mode,
            data.center,
            data.size,
            data.metadata.docked_ligand_file,
            data.metadata.buffer,
            data.in_path,
        )

        rec_box = utils.prepare_box(
            rec_sph,
            data.center,
            data.size,
            data.metadata.enclose_spheres,
            data.metadata.buffer,
            data.in_path,
        )
        if rec_box is None:
            return data

        grid_stem = utils.prepare_grid(rec_mol2, rec_box, data.in_path)
        if grid_stem is None:
            return data

        data.metadata.prepared_receptor = rec_sph, grid_stem
        return data

    @staticmethod
    def prepare_and_run(data: CalculationData) -> CalculationData:
        DOCKRunner.prepare_ligand(data)
        DOCKRunner.run(data)

        return data

    @staticmethod
    def prepare_ligand(data: CalculationData) -> CalculationData:
        if data.smi is not None:
            DOCKRunner.prepare_from_smi(data)
        else:
            DOCKRunner.prepare_from_file(data)

        return data

    @staticmethod
    def prepare_from_smi(data: CalculationData) -> CalculationData:
        """Prepare an input ligand file from the ligand's SMILES string

        Parameters
        ----------
        data: CalculationData

        Returns
        -------
        CalculationData
        """
        mol2 = Path(data.in_path) / f"{data.name}.mol2"

        mol = Chem.AddHs(Chem.MolFromSmiles(data.smi))
        Chem.EmbedMolecule(mol)
        Chem.MMFFOptimizeMolecule(mol)

        try:
            mol = pybel.readstring("mol", Chem.MolToMolBlock(mol))
            # mol = pybel.readstring(format="smi", string=data.smi)
            # mol.make3D()
            # mol.addh()
            mol.calccharges(model="gasteiger")
        except Exception:
            pass

        mol.write(format="mol2", filename=str(mol2), overwrite=True, opt={"h": None})

        data.metadata.prepared_ligand = mol2
        return data

    @staticmethod
    def prepare_from_file(data: CalculationData) -> Optional[Tuple]:
        """Convert a single ligand to the appropriate input format with specified geometry"""
        fmt = Path(data.input_file).suffix.strip(".")
        mols = list(pybel.readfile(fmt, data.input_file))
        mol = mols[0]

        mol2 = Path(data.in_path) / f"{mol.title or data.name}.mol2"
        data.smi = mol.write()
        try:
            mol.addh()
            mol.calccharges(model="gasteiger")
        except Exception:
            pass

        mol.write(format="mol2", filename=mol2.absolute().as_posix(), overwrite=True, opt={"h": None})
        data.metadata.prepared_ligand = mol2

        return data

    @staticmethod
    def run(data: CalculationData) -> Optional[float]:
        """Dock this ligand into the ensemble of receptors

        Returns
        -------
        score : Optional[float]
            the ligand's docking score. None if DOCKing failed.
        """
        p_ligand = Path(data.metadata.prepared_ligand)
        ligand_name = p_ligand.stem
        sph_file, grid_prefix = data.metadata.prepared_receptor

        name = f"{Path(sph_file).stem}_{ligand_name}"

        infile, outfile_prefix = DOCKRunner.prepare_input_file(
            p_ligand, sph_file, grid_prefix, name, data.in_path, data.out_path
        )

        logfile = Path(outfile_prefix).parent / f"{name}.log"
        argv = [str(DOCK), "-i", str(infile), "-o", str(logfile)]

        ret = sp.run(argv, stdout=sp.PIPE, stderr=sp.PIPE)
        try:
            ret.check_returncode()
        except sp.SubprocessError:
            print(f"ERROR: docking failed. argv: {argv}", file=sys.stderr)
            print(f'Message: {ret.stderr.decode("utf-8")}', file=sys.stderr)

        scores = DOCKRunner.parse_logfile(logfile)
        score = None if scores is None else calc_score(scores, data.score_mode, data.k)

        data.result = Result(
            data.smi, name, re.sub("[:,.]", "", ray.state.current_node_id()), score
        )

        return scores

    @staticmethod
    def validate_metadata(metadata: DOCKMetadata):
        return

    @staticmethod
    def parse_logfile(outfile: Union[str, Path]) -> Optional[float]:
        """parse a DOCK log file for the scores of the conformations

        Parameters
        ----------
        outfile : Union[str, Path]
            the filepath of a scored log file generated by DOCK6

        Returns
        -------
        Optional[float]
            the DOCKing score of the ligand. None if no scores were parsed or the log file was
            unparseable
        """
        try:
            with open(outfile) as fid:
                score_lines = [line for line in fid if "Grid_Score" in line]
        except OSError:
            pass

        scores = []
        for line in score_lines:
            try:
                scores.append(float(line.split()[1]))
            except ValueError:
                continue

        return scores or None

    @staticmethod
    def prepare_input_file(
        ligand_file: Union[str, Path],
        sph_file: str,
        grid_prefix: str,
        name: Optional[str] = None,
        in_path: Union[str, Path] = ".",
        out_path: Union[str, Path] = ".",
        **kwargs,
    ) -> Tuple[Path, Path]:
        """Prepare an input file with which to run DOCK

        Parameters
        ----------
        ligand_file : str
            the MOL2 file corresponding to the ligand that will be docked
        sph_file : str
            the SPH file containing the DOCK spheres of the receptor
        grid_prefix : str
            the prefix of the prepared grid files (as was passed to the grid program)
        name : Optional[str], default=None
            the name to use for the input file and output file
        in_path : Union[str, os.PathLike], default="inputs"
            the path under which to write the input files both the input file and output
        out_path : Union[str, os.PathLike], default="outputs"
            the path under which to write the output files
        **kwargs
            keyword options DOCKing parameters

        Returns
        -------
        infile: Path
            the filepath of the input file
        outfile_prefix: Path
            the prefix of the outfile name. DOCK will automatically name outfiles
            as <outfile_prefix>_scored.mol2
        """
        name = name or f"{Path(sph_file).stem}_{Path(ligand_file).stem}"

        infile = in_path / f"{name}.in"
        outfile_prefix = out_path / name

        with open(infile, "w") as fid:
            # fid.write("conformer_search_type flex\n")
            # fid.write("write_fragment_libraries no\n")
            # fid.write("user_specified_anchor no\n")
            # fid.write("limit_max_anchors no\n")
            # fid.write("min_anchor_size 5\n")

            # fid.write("pruning_use_clustering yes\n")
            # fid.write("pruning_max_orients 100\n")
            # fid.write("pruning_clustering_cutoff 100\n")
            # fid.write("pruning_conformer_score_cutoff 100.0\n")
            # fid.write("pruning_conformer_score_scaling_factor 1.0\n")

            # fid.write("use_clash_overlap no\n")
            # fid.write("write_growth_tree no\n")
            # fid.write("use_internal_energy yes\n")
            # fid.write("internal_energy_rep_exp 12\n")
            # fid.write("internal_energy_cutoff 100.0\n")

            # fid.write(f"ligand_atom_file {ligand_file}\n")
            # fid.write("limit_max_ligands no\n")
            # fid.write("skip_molecule no\n")
            # fid.write("read_mol_solvation no\n")
            # fid.write("calculate_rmsd no\n")
            # fid.write("use_rmsd_reference_mol no\n")
            # fid.write("use_database_filter no\n")
            # fid.write("orient_ligand yes\n")
            # fid.write("automated_matching yes\n")
            # fid.write(f"receptor_site_file {sph_file}\n")
            # fid.write("max_orientations 1000\n")
            # fid.write("critical_points no\n")
            # fid.write("chemical_matching no\n")
            # fid.write("use_ligand_spheres no\n")
            # fid.write("bump_filter no\n")

            # fid.write("score_molecules yes\n")
            # fid.write("contact_score_primary no\n")
            # fid.write("contact_score_secondary no\n")
            # fid.write("grid_score_primary yes\n")
            # fid.write("grid_score_secondary no\n")
            # fid.write("grid_score_rep_rad_scale 1\n")
            # fid.write("grid_score_vdw_scale 1\n")
            # fid.write("grid_score_es_scale 1\n")
            # fid.write(f"grid_score_grid_prefix {grid_prefix}\n")
            # fid.write("multigrid_score_secondary no\n")
            # fid.write("dock3.5_score_secondary no\n")
            # fid.write("continuous_score_secondary no\n")
            # fid.write("footprint_similarity_score_secondary no\n")
            # fid.write("pharmacophore_score_secondary no\n")
            # fid.write("descriptor_score_secondary no\n")
            # fid.write("gbsa_zou_score_secondary no\n")
            # fid.write("gbsa_hawkins_score_secondary no\n")
            # fid.write("SASA_score_secondary no\n")
            # fid.write("amber_score_secondary no\n")

            # fid.write("minimize_ligand yes\n")
            # fid.write("minimize_anchor yes\n")
            # fid.write("minimize_flexible_growth yes\n")
            # fid.write("use_advanced_simplex_parameters no\n")

            # fid.write("simplex_max_cycles 1\n")
            # fid.write("simplex_score_converge 0.1\n")
            # fid.write("simplex_cycle_converge 1.0\n")
            # fid.write("simplex_trans_step 1.0\n")
            # fid.write("simplex_rot_step 0.1\n")
            # fid.write("simplex_tors_step 10.0\n")
            # fid.write("simplex_anchor_max_iterations 500\n")
            # fid.write("simplex_grow_max_iterations 500\n")
            # fid.write("simplex_grow_tors_premin_iterations 0\n")
            # fid.write("simplex_random_seed 0\n")
            # fid.write("simplex_restraint_min no\n")

            # fid.write("atom_model all\n")
            # fid.write(f"vdw_defn_file {VDW_DEFN_FILE}\n")
            # fid.write(f"flex_defn_file {FLEX_DEFN_FILE}\n")
            # fid.write(f"flex_drive_file {FLEX_DRIVE_FILE}\n")

            # fid.write(f"ligand_outfile_prefix {outfile_prefix}\n")
            # fid.write("write_orientations no\n")
            # fid.write("num_scored_conformers 5\n")
            # fid.write("write_conformations no\n")
            # fid.write("rank_ligands no\n")

            fid.write(DOCKRunner.infile_line(kwargs, "conformer_search_type", "flex"))
            fid.write(DOCKRunner.infile_line(kwargs, "write_fragment_libraries", "no"))

            fid.write(DOCKRunner.infile_line(kwargs, "user_specified_anchor", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "limit_max_anchors", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "min_anchor_size", "5"))
            fid.write(DOCKRunner.infile_line(kwargs, "pruning_use_clustering", "yes"))
            fid.write(DOCKRunner.infile_line(kwargs, "pruning_max_orients", "100"))
            fid.write(DOCKRunner.infile_line(kwargs, "pruning_clustering_cutoff", "100"))
            fid.write(DOCKRunner.infile_line(kwargs, "pruning_conformer_score_cutoff", "100"))
            fid.write(DOCKRunner.infile_line(kwargs, "pruning_conformer_score_scaling_factor", "1"))
            fid.write(DOCKRunner.infile_line(kwargs, "use_clash_overlap", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "write_growth_tree", "no"))

            fid.write(DOCKRunner.infile_line(kwargs, "use_internal_energy", "yes"))
            fid.write(DOCKRunner.infile_line(kwargs, "internal_energy_rep_exp", "12"))
            fid.write(DOCKRunner.infile_line(kwargs, "internal_energy_cutoff", "100"))

            fid.write(f"ligand_atom_file {ligand_file}\n")
            fid.write(DOCKRunner.infile_line(kwargs, "limit_max_ligands", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "skip_molecule", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "read_mol_solvation", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "calculate_rmsd", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "use_rmsd_reference_mol", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "use_database_filter", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "orient_ligand", "yes"))
            fid.write(DOCKRunner.infile_line(kwargs, "automated_matching", "yes"))
            fid.write(f"receptor_site_file {sph_file}\n")
            fid.write(DOCKRunner.infile_line(kwargs, "max_orientations", "1000"))
            fid.write(DOCKRunner.infile_line(kwargs, "critical_points", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "chemical_matching", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "use_ligand_spheres", "no"))

            fid.write(DOCKRunner.infile_line(kwargs, "bump_filter", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "score_molecules", "yes"))
            fid.write(DOCKRunner.infile_line(kwargs, "contact_score_primary", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "contact_score_secondary", "no"))

            fid.write(DOCKRunner.infile_line(kwargs, "grid_score_primary", "yes"))
            fid.write(DOCKRunner.infile_line(kwargs, "grid_score_secondary", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "grid_score_rep_rad_scale", "1"))
            fid.write(DOCKRunner.infile_line(kwargs, "grid_score_vdw_scale", "1"))
            fid.write(DOCKRunner.infile_line(kwargs, "grid_score_es_scale", "1"))
            fid.write(f"grid_score_grid_prefix {grid_prefix}\n")
            
            fid.write(DOCKRunner.infile_line(kwargs, "multigrid_score_secondary", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "5_score_secondary", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "continuous_score_secondary", "no"))
            fid.write(DOCKRunner.infile_line(kwargs,"footprint_similarity_score_secondary", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "pharmacophore_score_secondary", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "descriptor_score_secondary", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "gbsa_zou_score_secondary", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "gbsa_hawkins_score_secondary", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "SASA_score_secondary", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "amber_score_secondary", "no"))

            fid.write(DOCKRunner.infile_line(kwargs, "minimize_ligand", "yes"))
            fid.write(DOCKRunner.infile_line(kwargs, "minimize_anchor", "yes"))
            fid.write(DOCKRunner.infile_line(kwargs, "minimize_flexible_growth", "yes"))
            fid.write(DOCKRunner.infile_line(kwargs, "use_advanced_simplex_parameters", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "simplex_max_cycles", "1"))
            fid.write(DOCKRunner.infile_line(kwargs, "simplex_score_converge", "0.1"))
            fid.write(DOCKRunner.infile_line(kwargs, "simplex_cycle_converge", "1.0"))
            fid.write(DOCKRunner.infile_line(kwargs, "simplex_trans_step", "1.0"))
            fid.write(DOCKRunner.infile_line(kwargs, "simplex_rot_step", "0.1"))
            fid.write(DOCKRunner.infile_line(kwargs, "simplex_tors_step", "10"))
            fid.write(DOCKRunner.infile_line(kwargs, "simplex_anchor_max_iterations", "500"))
            fid.write(DOCKRunner.infile_line(kwargs, "simplex_grow_max_iterations", "500"))
            fid.write(DOCKRunner.infile_line(kwargs, "simplex_grow_tors_premin_iterations", "0"))
            fid.write(DOCKRunner.infile_line(kwargs, "simplex_random_seed", "0"))
            fid.write(DOCKRunner.infile_line(kwargs, "simplex_restraint_min", "no"))

            fid.write(DOCKRunner.infile_line(kwargs, "atom_model", "all"))
            fid.write(f"vdw_defn_file {VDW_DEFN_FILE}\n")
            fid.write(f"flex_defn_file {FLEX_DEFN_FILE}\n")
            fid.write(f"flex_drive_file {FLEX_DRIVE_FILE}\n")
        
            fid.write(f"ligand_outfile_prefix {outfile_prefix}\n")
            fid.write(DOCKRunner.infile_line(kwargs, "write_orientations", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "num_scored_conformers", "5"))
            fid.write(DOCKRunner.infile_line(kwargs, "write_conformations", "no"))
            fid.write(DOCKRunner.infile_line(kwargs, "rank_ligands", "no"))

        return infile, outfile_prefix

    def infile_line(options: Mapping, param: str, default: str) -> str:
        """generate a line in the infile for the parameter and its default value. If the parameter
        is present in the options dictionary, substitute that value."""
        return f"{param} {options.get(param, default)}\n"
