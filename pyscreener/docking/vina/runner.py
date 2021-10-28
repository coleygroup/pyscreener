from itertools import takewhile
from pathlib import Path
import re
import shutil
import subprocess as sp
import sys
from typing import List, Optional, Tuple, Union

from openbabel import pybel
from rdkit.Chem import AllChem as Chem
import ray

from pyscreener import utils
from pyscreener.exceptions import MissingExecutableError
from pyscreener.docking.data import CalculationData
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
    @staticmethod
    def prepare(data: CalculationData) -> CalculationData:
        data = VinaRunner.prepare_receptor(data)
        data = VinaRunner.prepare_ligand(data)

        return data

    @staticmethod
    def prepare_receptor(data: CalculationData) -> CalculationData:
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
        name = Path(data.receptor).with_suffix(".pdbqt").name
        receptor_pdbqt = Path(data.in_path) / name

        argv = ["prepare_receptor", "-r", data.receptor, "-o", receptor_pdbqt]
        try:
            ret = sp.run(argv, stderr=sp.PIPE)
            ret.check_returncode()
        except sp.SubprocessError:
            print(f'ERROR: failed to convert "{data.receptor}"', file=sys.stderr)
            print(ret.stderr.decode("utf-8"))
            return None

        data.metadata.prepared_receptor = receptor_pdbqt
        return data

    @staticmethod
    def prepare_and_run(data: CalculationData) -> CalculationData:
        VinaRunner.prepare_ligand(data)
        VinaRunner.run(data)

        return data

    @staticmethod
    def prepare_ligand(data: CalculationData) -> CalculationData:
        if data.smi is not None:
            VinaRunner.prepare_from_smi(data)
        else:
            VinaRunner.prepare_from_file(data)

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
        pdbqt = Path(data.in_path) / f"{data.name}.pdbqt"

        mol = Chem.AddHs(Chem.MolFromSmiles(data.smi))
        Chem.EmbedMolecule(mol)
        Chem.MMFFOptimizeMolecule(mol)

        try:
            mol = pybel.readstring("mol", Chem.MolToMolBlock(mol))
            # mol.make3D()
            # mol.addh()
            mol.calccharges(model="gasteiger")
        except Exception:
            pass

        mol.write("pdbqt", str(pdbqt), overwrite=True, opt={"h": None})
        data.metadata.prepared_ligand = pdbqt

        return data

    @staticmethod
    def prepare_from_file(data: CalculationData) -> CalculationData:
        """Convert the molecule contained in the arbitrary chemical file to a PDBQT file with the
        same geometry

        Parameters
        ----------
        data: CalculationData

        Returns
        -------
        CalculationData
        """
        fmt = Path(data.input_file).suffix.strip(".")
        mols = list(pybel.readfile(fmt, data.input_file))
        mol = mols[0]

        pdbqt = Path(data.in_path) / f"{mol.title or data.name}.pdbqt"
        data.smi = mol.write()
        try:
            mol.addh()
            mol.calccharges(model="gasteiger")
        except Exception:
            pass

        mol.write(format="pdbqt", filename=str(pdbqt), overwrite=True, opt={"h": None})
        data.metadata.prepared_ligand = pdbqt

        return data

    @staticmethod
    def run(data: CalculationData) -> Optional[List[float]]:
        """Dock the given ligand using the specified vina-type docking program
        and parameters

        Returns
        -------
        scores : Optional[List[float]]
            the conformer scores parsed from the log file
        """
        p_ligand = Path(data.metadata.prepared_ligand)
        ligand_name = p_ligand.stem

        name = f"{Path(data.receptor).stem}_{ligand_name}"

        argv, out, log = VinaRunner.build_argv(
            ligand=data.metadata.prepared_ligand,
            receptor=data.metadata.prepared_receptor,
            software=data.metadata.software,
            center=data.center,
            size=data.size,
            ncpu=data.ncpu,
            exhaustiveness=data.metadata.exhaustiveness,
            num_modes=data.metadata.num_modes,
            energy_range=data.metadata.energy_range,
            name=name,
            path=Path(data.out_path),
            extra=data.metadata.extra,
        )

        ret = sp.run(argv, stdout=sp.PIPE, stderr=sp.PIPE)
        try:
            ret.check_returncode()
        except sp.SubprocessError:
            print(f'ERROR: docking failed. Message: {ret.stderr.decode("utf-8")}', file=sys.stderr)

        scores = VinaRunner.parse_logfile(log)
        if scores is None:
            score = None
        else:
            score = utils.calc_score(scores, data.score_mode, data.k)

        data.result = Result(
            data.smi, name, re.sub("[:,.]", "", ray.state.current_node_id()), score
        )

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
        energy_range: float = 3.,
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
            scores were parsed or the log file was unparseable
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
                score_lines = [line for line in fid.readlines() if 'REMARK VINA RESULT' in line]
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