import os
from pathlib import Path
import re
import subprocess as sp
import sys
from typing import Optional, Sequence, Tuple, Union

from openbabel import pybel
import ray

from pyscreener.exceptions import (
    MisconfiguredDirectoryError,
    MissingEnvironmentVariableError,
)
from pyscreener.utils import calc_score
from pyscreener.docking import CalculationData, DockingRunner, Result
from pyscreener.docking.dock import utils

try:
    DOCK6 = Path(os.environ["DOCK6"])
except KeyError:
    raise MissingEnvironmentVariableError(
        "DOCK6 environment variable not set! "
        "See https://github.com/coleygroup/pyscreener#specifying-an-environment-variable for more information"
    )

VDW_DEFN_FILE = DOCK6 / "parameters" / "vdw_AMBER_parm99.defn"
FLEX_DEFN_FILE = DOCK6 / "parameters" / "flex.defn"
FLEX_DRIVE_FILE = DOCK6 / "parameters" / "flex_drive.tbl"
DOCK = DOCK6 / "bin" / "dock6"

for f in (VDW_DEFN_FILE, FLEX_DEFN_FILE, FLEX_DRIVE_FILE, DOCK):
    if not f.exists():
        raise MisconfiguredDirectoryError(
            "DOCK6 directory not configured properly! "
            f'DOCK6 path is set as "{DOCK6}", but there is no '
            f'"{f.name}" located in the "{f.parents[0].name}" subdirectory '
            "under the DOCK6 path. See https://github.com/coleygroup/pyscreener#specifying-an-environment-variable for more information"
        )


class DOCKRunner(DockingRunner):
    @staticmethod
    def prepare(data: CalculationData) -> CalculationData:
        data = DOCKRunner.prepare_receptor(data)
        # TODO(degraff): fix this to accept input ligand files
        data = DOCKRunner.prepare_from_smi(data)

        return data

    @staticmethod
    def prepare(data: CalculationData) -> CalculationData:
        data = DOCKRunner.prepare_receptor(data)
        data = DOCKRunner.prepare_ligand(data)

        return data

    @staticmethod
    def prepare_receptor(data: CalculationData) -> CalculationData:
        """Prepare the files necessary to dock ligands against the input
        receptor using this Screener's parameters

        Parameter
        ---------
        receptor : str
            the filepath of a file containing a receptor. Must be in a format
            that is readable by Chimera

        Returns
        -------
        rec_sph : str
            the filepath of the file containing the selected spheres
        grid_prefix : str
            the prefix of all prepared grid files.
        None
            if receptor preparation fails at any point
        """
        # receptor_pdbqt = Path(data.receptor).with_suffix('.pdbqt')
        # receptor_pdbqt = Path(data.in_path) / receptor_pdbqt.name
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
            data.metadata.in_path,
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

        try:
            mol = pybel.readstring(format="smi", string=data.smi)
            mol.make3D()
            mol.addh()
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

        mol.write(format="mol2", filename=mol2, overwrite=True, opt={"h": None})
        data.metadata.prepared_ligand = mol2

        return data
        # name = name or Path(filename).stem

        # ret = sp.run(['obabel', filename, '-osmi'], stdout=sp.PIPE, check=True)
        # lines = ret.stdout.decode('utf-8').splitlines()
        # smis = [line.split()[0] for line in lines]

        # if not use_3d:
        #     ligands = [
        #         DOCKRunner.prepare_from_smi(smi, f'{name}_{i}', path)
        #         for i, smi in enumerate(smis)
        #     ]
        #     return [lig for lig in ligands if lig]

        # path = Path(path)
        # if not path.is_dir():
        #     path.mkdir()

        # mol2 = f'{path}/{name}_.mol2'
        # argv = ['obabel', filename, '-omol2', '-O', mol2, '-m',
        #         '-h', '--partialcharge', 'gasteiger']

        # ret = sp.run(argv, check=False, stderr=sp.PIPE)
        # try:
        #     ret.check_returncode()
        # except sp.SubprocessError:
        #     return None

        # stderr = ret.stderr.decode('utf-8')
        # for line in stderr.splitlines():
        #     if 'converted' not in line:
        #         continue
        #     n_mols = int(line.split()[0])

        # mol2s = [f'{path}/{name}_{i}.mol2' for i in range(1, n_mols)]

        # return list(zip(smis, mol2s))

    @staticmethod
    def run(data: CalculationData) -> Optional[Sequence[float]]:
        """Dock this ligand into the ensemble of receptors

        Parameters
        ----------
        ligand : Tuple[str, str]
            a tuple containing the ligand's SMILES string and its prepared
            .mol2 file that will be docked against each receptor
        receptors : List[Tuple[str, str]]
            a list of tuples containing the sphere file and grid file prefix
            corresponding to each receptor in the ensemble.
        in_path : Union[str, os.PathLike] (Default = 'inputs')
            the path under which to write the input files
        out_path : Union[str, os.PathLike] (Default = 'outputs')
            the path under which to write the output files
        repeats : int (Default = 1)
            the number of times each docking run should be repeated

        Returns
        -------
        ensemble_rowss : List[List[Dict]]
            an MxO list of dictionaries where each dictionary is a record of an
            individual docking run and:
            - M is the number of receptors each ligand is docked against
            - O is the number of times each docking run is repeated.
            Each dictionary contains the following keys:
            - smiles: the ligand's SMILES string
            - name: the name of the ligand
            - in: the filename of the input ligand file
            - out: the filename of the output docked ligand file
            - log: the filename of the output log file
            - score: the ligand's docking score
        """
        p_ligand = Path(data.metadata.prepared_ligand)
        ligand_name = p_ligand.stem
        sph_file, grid_prefix = data.metadata.prepared_receptor

        name = f"{Path(sph_file).stem}_{ligand_name}"

        infile, outfile_prefix = DOCKRunner.prepare_input_file(
            p_ligand, sph_file, grid_prefix, name, data.in_path, data.out_path
        )

        # out = Path(f'{outfile_prefix}_scored.mol2')
        log = Path(outfile_prefix).parent / f"{name}.out"
        argv = [str(DOCK), "-i", infile, "-o", log]

        ret = sp.run(argv, stdout=sp.PIPE, stderr=sp.PIPE)
        try:
            ret.check_returncode()
        except sp.SubprocessError:
            print(f"ERROR: docking failed. argv: {argv}", file=sys.stderr)
            print(f'Message: {ret.stderr.decode("utf-8")}', file=sys.stderr)

        scores = DOCKRunner.parse_log_file(log)
        if scores is None:
            score = None
        else:
            score = calc_score(scores, data.score_mode, data.k)

        data.result = Result(
            data.smi, name, re.sub("[:,.]", "", ray.state.current_node_id()), score
        )

        return scores

    def parse_outfile(outfile: Union[str, Path]) -> Optional[float]:
        """parse a DOCK out file for the scores of the conformations

        Parameters
        ----------
        outfile : Union[str, PathLike]
            the filename of a scored outfile file generated by DOCK6 or a
            PathLike object pointing to that file

        Returns
        -------
        Optional[List[float]]
            the scores of the docked conformations in the ordering of the
            out file. None if no scores were parsed or the log file was
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
                scores.append(float(line.split()[2]))
            except ValueError:
                continue
        # scores = [float(line.split()[2]) for line in score_lines]

        return scores or None

    def prepare_input_file(
        ligand_file: Union[str, Path],
        sph_file: str,
        grid_prefix: str,
        name: Optional[str] = None,
        in_path: Union[str, Path] = ".",
        out_path: Union[str, Path] = ".",
    ) -> Tuple[str, str]:
        """Prepare the input file with which to run DOCK

        Parameters
        ----------
        ligand_file : str
            the MOL2 file corresponding to the ligand that will be docked
        sph_file : str
            the .sph file containing the DOCK spheres of the receptor
        grid_prefix : str
            the prefix of the prepared grid files (as was passed to
            the grid program)
        name : Optional[str] (Default = None)
            the name to use for the input file and output file
        in_path : Union[str, os.PathLike] (Default = 'inputs')
            the path under which to write the input files
            both the input file and output
        out_path : Union[str, os.PathLike] (Default = 'outputs')
            the path under which to write the output files

        Returns
        -------
        infile: str
            the name of the input file
        outfile_prefix: str
            the prefix of the outfile name. DOCK will automatically name outfiles
            as <outfile_prefix>_scored.mol2
        """
        name = name or f"{Path(sph_file).stem}_{Path(ligand_file).stem}"

        infile = in_path / f"{name}.in"
        outfile_prefix = out_path / name

        with open(infile, "w") as fid:
            fid.write("conformer_search_type flex\n")
            fid.write("write_fragment_libraries no\n")
            fid.write("user_specified_anchor no\n")
            fid.write("limit_max_anchors no\n")
            fid.write("min_anchor_size 5\n")

            fid.write("pruning_use_clustering yes\n")
            fid.write("pruning_max_orients 100\n")
            fid.write("pruning_clustering_cutoff 100\n")
            fid.write("pruning_conformer_score_cutoff 100.0\n")
            fid.write("pruning_conformer_score_scaling_factor 1.0\n")

            fid.write("use_clash_overlap no\n")
            fid.write("write_growth_tree no\n")
            fid.write("use_internal_energy yes\n")
            fid.write("internal_energy_rep_exp 12\n")
            fid.write("internal_energy_cutoff 100.0\n")

            fid.write(f"ligand_atom_file {ligand_file}\n")
            fid.write("limit_max_ligands no\n")
            fid.write("skip_molecule no\n")
            fid.write("read_mol_solvation no\n")
            fid.write("calculate_rmsd no\n")
            fid.write("use_rmsd_reference_mol no\n")
            fid.write("use_database_filter no\n")
            fid.write("orient_ligand yes\n")
            fid.write("automated_matching yes\n")
            fid.write(f"receptor_site_file {sph_file}\n")
            fid.write("max_orientations 1000\n")
            fid.write("critical_points no\n")
            fid.write("chemical_matching no\n")
            fid.write("use_ligand_spheres no\n")
            fid.write("bump_filter no\n")
            fid.write("score_molecules yes\n")

            fid.write("contact_score_primary no\n")
            fid.write("contact_score_secondary no\n")

            fid.write("grid_score_primary yes\n")
            fid.write("grid_score_secondary no\n")
            fid.write("grid_score_rep_rad_scale 1\n")
            fid.write("grid_score_vdw_scale 1\n")
            fid.write("grid_score_es_scale 1\n")
            fid.write(f"grid_score_grid_prefix {grid_prefix}\n")

            fid.write("multigrid_score_secondary no\n")
            fid.write("dock3.5_score_secondary no\n")
            fid.write("continuous_score_secondary no\n")
            fid.write("footprint_similarity_score_secondary no\n")
            fid.write("pharmacophore_score_secondary no\n")
            fid.write("descriptor_score_secondary no\n")
            fid.write("gbsa_zou_score_secondary no\n")
            fid.write("gbsa_hawkins_score_secondary no\n")
            fid.write("SASA_score_secondary no\n")
            fid.write("amber_score_secondary no\n")

            fid.write("minimize_ligand yes\n")
            fid.write("minimize_anchor yes\n")
            fid.write("minimize_flexible_growth yes\n")
            fid.write("use_advanced_simplex_parameters no\n")

            fid.write("simplex_max_cycles 1\n")
            fid.write("simplex_score_converge 0.1\n")
            fid.write("simplex_cycle_converge 1.0\n")
            fid.write("simplex_trans_step 1.0\n")
            fid.write("simplex_rot_step 0.1\n")
            fid.write("simplex_tors_step 10.0\n")
            fid.write("simplex_anchor_max_iterations 500\n")
            fid.write("simplex_grow_max_iterations 500\n")
            fid.write("simplex_grow_tors_premin_iterations 0\n")
            fid.write("simplex_random_seed 0\n")
            fid.write("simplex_restraint_min no\n")

            fid.write("atom_model all\n")
            fid.write(f"vdw_defn_file {VDW_DEFN_FILE}\n")
            fid.write(f"flex_defn_file {FLEX_DEFN_FILE}\n")
            fid.write(f"flex_drive_file {FLEX_DRIVE_FILE}\n")

            fid.write(f"ligand_outfile_prefix {outfile_prefix}\n")
            fid.write("write_orientations no\n")
            fid.write("num_scored_conformers 5\n")
            fid.write("write_conformations no\n")
            fid.write("rank_ligands no\n")

        return infile, outfile_prefix
