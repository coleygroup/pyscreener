import os
from pathlib import Path
import re
import subprocess as sp
from typing import List, Mapping, Optional, Tuple, Union
import warnings

from openbabel import pybel
from rdkit.Chem import AllChem as Chem
import ray

from pyscreener.exceptions import (
    MisconfiguredDirectoryError,
    MissingEnvironmentVariableError,
    ReceptorPreparationError,
)
from pyscreener.utils import reduce_scores
from pyscreener.warnings import ChargeWarning, ConformerWarning, SimulationFailureWarning
from pyscreener.docking import Simulation, DockingRunner, Result
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
            f'is no "{f.name}" located under the "{f.parents[0].name}" directory. '
            "See https://github.com/coleygroup/pyscreener#specifying-an-environment-variable for more information."
        )


class DOCKRunner(DockingRunner):
    @classmethod
    def is_multithreaded(cls) -> bool:
        return False

    @staticmethod
    def prepare(sim: Simulation) -> Simulation:
        _ = DOCKRunner.prepare_receptor(sim)
        _ = DOCKRunner.prepare_ligand(sim)

        return sim

    @staticmethod
    def prepare_receptor(sim: Simulation) -> Simulation:
        """Prepare the files necessary to dock ligands against the input receptor using this
        Screener's parameters

        Parameter
        ---------
        simulation : Simulation
            the simulation for which to prepare the DOCK6 receptor inputs

        Returns
        -------
        Simulation
        """
        in_path = sim.in_path / "receptors"
        in_path.mkdir(parents=True, exist_ok=True)

        try:
            rec_mol2 = utils.prepare_mol2(sim.receptor, in_path)
            rec_pdb = utils.prepare_pdb(sim.receptor, in_path)
            rec_dms = utils.prepare_dms(rec_pdb, sim.metadata.probe_radius, in_path)
            rec_sph = utils.prepare_sph(
                rec_dms,
                sim.metadata.steric_clash_dist,
                sim.metadata.min_radius,
                sim.metadata.max_radius,
                in_path,
            )
            rec_sph = utils.select_spheres(
                rec_sph,
                sim.metadata.sphere_mode,
                sim.center,
                sim.size,
                sim.metadata.docked_ligand_file,
                sim.metadata.buffer,
                in_path,
            )
            rec_box = utils.prepare_box(
                rec_sph,
                sim.center,
                sim.size,
                sim.metadata.enclose_spheres,
                sim.metadata.buffer,
                in_path,
            )
            grid_stem = utils.prepare_grid(rec_mol2, rec_box, in_path, sim.metadata.grid_params)
        except ReceptorPreparationError:
            raise
            # return None  # should think about whether to handle or just raise

        sim.metadata.prepared_receptor = rec_sph, grid_stem
        return sim

    @staticmethod
    def prepare_and_run(sim: Simulation) -> Optional[Result]:
        if not DOCKRunner.prepare_ligand(sim):
            return None

        _ = DOCKRunner.run(sim)

        return sim.result

    @staticmethod
    def prepare_ligand(sim: Simulation) -> bool:
        if sim.smi is not None:
            return DOCKRunner.prepare_from_smi(sim)

        return DOCKRunner.prepare_from_file(sim)

    @staticmethod
    def prepare_from_smi(sim: Simulation) -> bool:
        """Prepare an input ligand file from the ligand's SMILES string

        Parameters
        ----------
        sim: CalculationData

        Returns
        -------
        CalculationData
        """
        mol2 = Path(sim.in_path) / f"{sim.name}.mol2"

        mol = Chem.MolFromSmiles(sim.smi)
        if mol is None:
            return False

        mol = Chem.AddHs(mol)

        try:
            Chem.EmbedMolecule(mol)
            Chem.MMFFOptimizeMolecule(mol)
        except ValueError:
            warnings.warn("Could not generate 3D conformer of molecule!", ConformerWarning)

        try:
            mol = pybel.readstring("mol", Chem.MolToMolBlock(mol))
        except IOError:
            return False

        try:
            mol.calccharges(model="gasteiger")
        except Exception:
            warnings.warn("Could not calculate charges for ligand!", ChargeWarning)

        mol.write(format="mol2", filename=str(mol2), overwrite=True, opt={"h": None})
        sim.metadata.prepared_ligand = mol2

        return True

    @staticmethod
    def prepare_from_file(sim: Simulation) -> bool:
        """Convert a single ligand to the appropriate input format with specified geometry"""
        fmt = Path(sim.input_file).suffix.strip(".")

        try:
            mol = next(pybel.readfile(fmt, sim.input_file))
        except StopIteration:
            return False

        p_mol2 = Path(sim.in_path) / f"{mol.title or sim.name}.mol2"

        sim.smi = mol.write()
        mol.addh()
        try:
            mol.calccharges(model="gasteiger")
        except Exception:
            warnings.warn("Could not calculate charges for ligand!", ChargeWarning)

        mol.write(format="mol2", filename=str(p_mol2), overwrite=True, opt={"h": None})
        sim.metadata.prepared_ligand = p_mol2

        return True

    @staticmethod
    def run(sim: Simulation) -> Optional[List[float]]:
        """Dock this ligand into the ensemble of receptors

        Returns
        -------
        score : Optional[float]
            the ligand's docking score. None if DOCKing failed.
        """
        if sim.metadata.prepared_receptor is None or sim.metadata.prepared_ligand is None:
            return None

        p_ligand = Path(sim.metadata.prepared_ligand)
        ligand_name = p_ligand.stem
        sph_file, grid_prefix = sim.metadata.prepared_receptor

        name = f"{Path(sph_file).stem}_{ligand_name}"

        infile, outfile_prefix = DOCKRunner.prepare_input_file(
            p_ligand,
            sph_file,
            grid_prefix,
            name,
            sim.in_path,
            sim.out_path,
            sim.metadata.dock_params,
        )

        logfile = Path(outfile_prefix).parent / f"{name}.log"
        argv = [str(DOCK), "-i", str(infile), "-o", str(logfile)]

        ret = sp.run(argv, stdout=sp.PIPE, stderr=sp.PIPE)
        try:
            ret.check_returncode()
        except sp.SubprocessError:
            warnings.warn(f'Message: {ret.stderr.decode("utf-8")}', SimulationFailureWarning)

        scores = DOCKRunner.parse_logfile(logfile)
        score = None if scores is None else reduce_scores(scores, sim.reduction, k=sim.k)

        sim.result = Result(sim.smi, name, re.sub("[:,.]", "", ray.state.current_node_id()), score)

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
        params: Optional[Mapping] = None,
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
        DEFAULT_PARAMS = {
            "conformer_search_type": "flex",
            "write_fragment_libraries": "no",
            "user_specified_anchor": "no",
            "limit_max_anchors": "no",
            "min_anchor_size": 5,
            "pruning_use_clustering": "yes",
            "pruning_max_orients": 100,
            "pruning_clustering_cutoff": 100,
            "pruning_conformer_score_cutoff": 100,
            "pruning_conformer_score_scaling_factor": 1,
            "use_clash_overlap": "no",
            "write_growth_tree": "no",
            "use_internal_energy": "yes",
            "internal_energy_rep_exp": 12,
            "internal_energy_cutoff": 100,
            "limit_max_ligands": "no",
            "skip_molecule": "no",
            "read_mol_solvation": "no",
            "calculate_rmsd": "no",
            "use_rmsd_reference_mol": "no",
            "use_database_filter": "no",
            "orient_ligand": "yes",
            "automated_matching": "yes",
            "max_orientations": 1000,
            "critical_points": "no",
            "chemical_matching": "no",
            "use_ligand_spheres": "no",
            "bump_filter": "no",
            "score_molecules": "yes",
            "contact_score_primary": "no",
            "contact_score_secondary": "no",
            "grid_score_primary": "yes",
            "grid_score_secondary": "no",
            "grid_score_rep_rad_scale": 1,
            "grid_score_vdw_scale": 1,
            "grid_score_es_scale": 1,
            "multigrid_score_secondary": "no",
            "5_score_secondary": "no",
            "continuous_score_secondary": "no",
            "footprint_similarity_score_secondary": "no",
            "pharmacophore_score_secondary": "no",
            "descriptor_score_secondary": "no",
            "gbsa_zou_score_secondary": "no",
            "gbsa_hawkins_score_secondary": "no",
            "SASA_score_secondary": "no",
            "amber_score_secondary": "no",
            "minimize_ligand": "yes",
            "minimize_anchor": "yes",
            "minimize_flexible_growth": "yes",
            "use_advanced_simplex_parameters": "no",
            "simplex_max_cycles": 1,
            "simplex_score_converge": 0.1,
            "simplex_cycle_converge": 1.0,
            "simplex_trans_step": 1.0,
            "simplex_rot_step": 0.1,
            "simplex_tors_step": 10,
            "simplex_anchor_max_iterations": 500,
            "simplex_grow_max_iterations": 500,
            "simplex_grow_tors_premin_iterations": 0,
            "simplex_random_seed": 0,
            "simplex_restraint_min": "no",
            "atom_model": "all",
            "write_orientations": "no",
            "num_scored_conformers": 5,
            "write_conformations": "no",
            "rank_ligands": "no",
        }
        if params is None:
            params = DEFAULT_PARAMS
        else:
            params = {**params, **DEFAULT_PARAMS}
        name = name or f"{Path(sph_file).stem}_{Path(ligand_file).stem}"

        infile = in_path / f"{name}.in"
        outfile_prefix = out_path / name

        with open(infile, "w") as fid:
            [print(f"{k} {v}", file=fid) for k, v in params.items()]
            print(f"vdw_defn_file {VDW_DEFN_FILE}", file=fid)
            print(f"flex_defn_file {FLEX_DEFN_FILE}", file=fid)
            print(f"flex_drive_file {FLEX_DRIVE_FILE}", file=fid)
            print(f"ligand_atom_file {ligand_file}", file=fid)
            print(f"receptor_site_file {sph_file}", file=fid)
            print(f"grid_score_grid_prefix {grid_prefix}", file=fid)
            print(f"ligand_outfile_prefix {outfile_prefix}", file=fid)

        return infile, outfile_prefix
