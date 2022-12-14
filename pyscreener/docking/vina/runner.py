from itertools import takewhile
from pathlib import Path
import re
import shutil
import subprocess as sp
from typing import List, Optional, Tuple, Union
import warnings

from openbabel import pybel
import numpy as np
from rdkit.Chem import AllChem as Chem
import ray

from pyscreener import utils
from pyscreener.exceptions import MissingExecutableError, ReceptorPreparationError
from pyscreener.warnings import ChargeWarning, ConformerWarning, SimulationFailureWarning
from pyscreener.docking.sim import Simulation
from pyscreener.docking.runner import DockingRunner
from pyscreener.docking.result import Result
from pyscreener.docking.vina.metadata import VinaMetadata
from pyscreener.docking.vina.utils import Software

if shutil.which("prepare_receptor") is None:
    raise MissingExecutableError(
        'Could not find "prepare_receptor" on PATH! '
        "See https://github.com/coleygroup/pyscreener#adding-an-executable-to-your-path for more information."
    )


class VinaRunner(DockingRunner):
    @classmethod
    def is_multithreaded(cls) -> bool:
        return True

    @staticmethod
    def prepare(sim: Simulation) -> Simulation:
        _ = VinaRunner.prepare_receptor(sim)
        _ = VinaRunner.prepare_ligand(sim)

        return sim

    @staticmethod
    def prepare_receptor(sim: Simulation) -> Simulation:
        """Prepare a receptor PDBQT file from its input file

        Parameters
        ----------
        receptor : str
            the filename of a file containing a receptor

        Returns
        -------
        receptor_pdbqt : Optional[str]
            the filepath of the resulting PDBQT file. None if preparation failed
        """
        name = Path(sim.receptor).with_suffix(".pdbqt").name
        receptor_pdbqt = Path(sim.in_path) / name

        argv = ["prepare_receptor", "-r", sim.receptor, "-o", receptor_pdbqt]
        try:
            ret = sp.run(argv, stderr=sp.PIPE)
            ret.check_returncode()
        except sp.SubprocessError:
            raise ReceptorPreparationError(
                f'failed to prepare "{sim.receptor}". Message: {ret.stderr.decode("utf-8")}'
            )

        sim.metadata.prepared_receptor = receptor_pdbqt

        return sim

    @staticmethod
    def prepare_and_run(sim: Simulation) -> Optional[Result]:
        if not VinaRunner.prepare_ligand(sim):
            return None

        _ = VinaRunner.run(sim)

        return sim.result

    @staticmethod
    def prepare_ligand(sim: Simulation) -> bool:
        if sim.smi is not None:
            return VinaRunner.prepare_from_smi(sim)

        return VinaRunner.prepare_from_file(sim)

    @staticmethod
    def prepare_from_smi(sim: Simulation) -> bool:
        """Prepare the ligand PDQBT file from its SMILES string

        If successful, sets the `prepared_ligand` attribute of the metadata  and return True.
        Otherwise, do nothing and return False.

        Parameters
        ----------
        sim : Simulation

        Returns
        -------
        bool
            whether the ligand preparation succeeded
        """
        pdbqt = Path(sim.in_path) / f"{sim.name}.pdbqt"

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

        mol.write("pdbqt", str(pdbqt), overwrite=True, opt={"h": None})
        sim.metadata.prepared_ligand = pdbqt

        return True

    @staticmethod
    def prepare_from_file(sim: Simulation) -> bool:
        """Prepare the ligand PDQBT file from its input chemical file, retaining the input
        geometry.

        Sets the `prepared_ligand` attribute of the metadata. If ligand preparation fails for any
        reason, set the value to `None` and return.

        Parameters
        ----------
        sim : Simulation

        Returns
        -------
        bool
            whether the ligand preparation succeeded
        """
        fmt = Path(sim.input_file).suffix.strip(".")

        mol = next(pybel.readfile(fmt, sim.input_file))

        pdbqt = Path(sim.in_path) / f"{mol.title or sim.name}.pdbqt"
        sim.smi = mol.write()
        try:
            mol.addh()
            mol.calccharges(model="gasteiger")
        except Exception:
            warnings.warn("Could not calculate charges for ligand!", ChargeWarning)

        mol.write(format="pdbqt", filename=str(pdbqt), overwrite=True, opt={"h": None})
        sim.metadata.prepared_ligand = pdbqt

        return True

    @staticmethod
    def run(sim: Simulation) -> Optional[List[float]]:
        """run the given ligand using the specified vina-type docking program and parameters

        If the simulation is not possible due to the prepared inputs not being set beforehand, do
        nothing and return None. Otherwise, if the simulation itself fails, set the `result`
        attribute of the simulation with a score of `None` **and** return `None`.

        Returns
        -------
        scores : Optional[List[float]]
            the conformer scores parsed from the log file. None if no scores were parseable from
            the logfile due to simulation failure.
        """
        if sim.metadata.prepared_receptor is None or sim.metadata.prepared_ligand is None:
            return None

        p_ligand = Path(sim.metadata.prepared_ligand)
        ligand_name = p_ligand.stem

        name = f"{Path(sim.receptor).stem}_{ligand_name}"

        argv, _, log = VinaRunner.build_argv(
            ligand=sim.metadata.prepared_ligand,
            receptor=sim.metadata.prepared_receptor,
            software=sim.metadata.software,
            center=sim.center,
            size=sim.size,
            ncpu=sim.ncpu,
            exhaustiveness=sim.metadata.exhaustiveness,
            num_modes=sim.metadata.num_modes,
            energy_range=sim.metadata.energy_range,
            name=name,
            path=Path(sim.out_path),
            extra=sim.metadata.extra,
        )

        ret = sp.run(argv, stdout=sp.PIPE, stderr=sp.PIPE)
        try:
            ret.check_returncode()
        except sp.SubprocessError:
            warnings.warn(f'Message: {ret.stderr.decode("utf-8")}', SimulationFailureWarning)

        scores = VinaRunner.parse_logfile(log)
        if scores is None:
            score = None
        else:
            score = utils.reduce_scores(np.array(scores), sim.reduction, k=sim.k)

        sim.result = Result(sim.smi, name, re.sub("[:,.]", "", ray.state.current_node_id()), score)

        return scores

    @staticmethod
    def build_argv(
        ligand: str,
        receptor: str,
        software: Software,
        center: Tuple[float, float, float],
        size: Tuple[float, float, float] = (10, 10, 10),
        ncpu: int = 1,
        exhaustiveness: int = 8,
        num_modes: int = 9,
        energy_range: float = 3.0,
        name: Optional[str] = None,
        path: Path = Path("."),
        extra: Optional[List[str]] = None,
    ) -> Tuple[List[str], Path, Path]:
        """Builds the argument vector to run a vina-type docking program

        Parameters
        ----------
        ligand : str
            the filename of the input ligand PDBQT file
        receptor : str
            the filename of the input receptor PDBQT file
        software : Software
            the docking program to run
        center : Tuple[float, float, float]
            the x-, y-, and z-coordinates of the center of the search box
        size : Tuple[float, float, float], default=(10, 10, 10)
            the  x-, y-, and z-radii, respectively, of the search box
        ncpu : int, default=1
            the number of cores to allocate to the docking program
        exhaustiveness: int
            the exhaustiveness of the global search. Larger values are more exhaustive
        num_modes: int
            the number of output modes
        energy_range: float
            the maximum energy difference (in kcal/mol) between the best and worst output binding
            modes
        extra : Optional[List[str]]
            additional command line arguments that will be passed to the docking calculation
        name : string, default=<receptor>_<ligand>)
            the base name to use for both the log and out files
        path : Path, default=Path('.')
            the path under which both the log and out files should be written
        extra : Optional[List[str]], default=None
            additional command line arguments to pass to each run

        Returns
        -------
        argv : List[str]
            the argument vector with which to run an instance of a vina-type
            docking program
        out : Path
            the filepath of the out file which the docking program will write to
        log : Path
            the filepath of the log file which the docking program will write to
        """
        name = name or (Path(receptor).stem + "_" + Path(ligand).stem)
        extra = extra or []

        out = path / f"{software.value}_{name}_out.pdbqt"
        log = path / f"{software.value}_{name}.log"

        argv = [
            software.value,
            f"--receptor={receptor}",
            f"--ligand={ligand}",
            f"--center_x={center[0]}",
            f"--center_y={center[1]}",
            f"--center_z={center[2]}",
            f"--size_x={size[0]}",
            f"--size_y={size[1]}",
            f"--size_z={size[2]}",
            f"--cpu={ncpu}",
            f"--out={out}",
            f"--log={log}",
            f"--exhaustiveness={exhaustiveness}",
            f"--num_modes={num_modes}",
            f"--energy_range={energy_range}",
            *extra,
        ]

        return argv, out, log

    @staticmethod
    def parse_logfile(logfile: Union[str, Path]) -> Optional[List[float]]:
        """parse a Vina-type log file for the scores of the binding modes

        Parameters
        ----------
        logfile : Union[str, Path]
            the path to a Vina-type log file

        Returns
        -------
        Optional[List[float]]
            the scores of the docked binding modes in the ordering of the log file. None if no
            scores were parsed or the log file was unparsable
        """
        TABLE_BORDER = "-----+------------+----------+----------"
        try:
            with open(logfile) as fid:
                for line in fid:
                    if TABLE_BORDER in line:
                        break

                score_lines = list(takewhile(lambda line: "Writing" not in line, fid))
        except OSError:
            return None

        scores = []
        for line in score_lines:
            try:
                scores.append(float(line.split()[1]))
            except ValueError:
                continue

        return scores or None

    @staticmethod
    def parse_outfile(outfile: Union[str, Path]) -> Optional[List[float]]:
        """parse a Vina-type output file for the scores of the binding modes

        Paramaters
        ----------
        outfile : Union[str, Path]
            the filepath a vina-type output file

        Returns
        -------
        Optional[List[float]]
            the scores of the binding in the ordering of the output file. None if no scores were
            parsed or the log file was unparseable
        """
        try:
            with open(outfile) as fid:
                score_lines = [line for line in fid.readlines() if "REMARK VINA RESULT" in line]
        except OSError:
            return None

        scores = []
        for line in score_lines:
            try:
                scores.append(float(line.split()[3]))
            except ValueError:
                continue

        return scores or None

    @staticmethod
    def validate_metadata(metadata: VinaMetadata):
        if shutil.which(metadata.software.value) is None:
            raise MissingExecutableError(
                f'Could not find "{metadata.software.value}" on PATH! '
                "See https://github.com/coleygroup/pyscreener/tree/refactor#adding-an-executable-to-your-path for more information."
            )
