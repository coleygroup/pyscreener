from itertools import takewhile
from pathlib import Path
import re
import shutil
import subprocess as sp
import sys
from typing import List, Optional, Tuple

from openbabel import pybel
import ray

from pyscreener import utils
from pyscreener.docking import CalculationData, DockingRunner, Result
from pyscreener.docking.vina.metadata import VinaMetadata
from pyscreener.exceptions import MissingExecutableError
from pyscreener.docking.vina.utils import Software

if shutil.which("prepare_receptor") is None:
    raise MissingExecutableError(
        'Could not find "prepare_receptor" on PATH! '
        "See https://github.com/coleygroup/pyscreener/tree/refactor#adding-an-executable-to-your-path for more information."
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
        receptor_pdbqt = Path(data.receptor).with_suffix(".pdbqt")
        receptor_pdbqt = Path(data.in_path) / receptor_pdbqt.name

        argv = ["prepare_receptor", "-r", data.receptor, "-o", receptor_pdbqt]
        try:
            sp.run(argv, stderr=sp.PIPE, check=True)
        except sp.SubprocessError:
            print(f'ERROR: failed to convert "{data.receptor}"', file=sys.stderr)
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

        mol = pybel.readstring(format="smi", string=data.smi)
        try:
            mol.make3D()
            mol.addh()
            mol.calccharges(model="gasteiger")
        except Exception:
            pass

        mol.write(format="pdbqt", filename=str(pdbqt), overwrite=True, opt={"h": None})
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

        mol.write(format="pdbqt", filename=pdbqt, overwrite=True, opt={"h": None})
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

        argv, _, log = VinaRunner.build_argv(
            ligand=data.metadata.prepared_ligand,
            receptor=data.metadata.prepared_receptor,
            software=data.metadata.software,
            center=data.center,
            size=data.size,
            ncpu=data.ncpu,
            extra=data.metadata.extra,
            name=name,
            path=Path(data.out_path),
        )

        ret = sp.run(argv, stdout=sp.PIPE, stderr=sp.PIPE)
        try:
            ret.check_returncode()
        except sp.SubprocessError:
            print(f"ERROR: docking failed. argv: {argv}", file=sys.stderr)
            print(f'Message: {ret.stderr.decode("utf-8")}', file=sys.stderr)

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
        name: Optional[str] = None,
        path: Path = Path("."),
        extra: Optional[List[str]] = None,
    ) -> Tuple[List[str], str, str]:
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
        out : str
            the filepath of the out file which the docking program will write to
        log : str
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
            *extra,
        ]

        return argv, out, log

    @staticmethod
    def parse_logfile(logfile: str) -> Optional[List[float]]:
        """parse a Vina-type log file for the scores of the conformations

        Parameters
        ----------
        logfile : str
            the path to a Vina-type log file

        Returns
        -------
        Optional[List[float]]
            the scores of the docked conformations in the ordering of the
            log file. None if no scores were parsed or the log file was
            unparseable
        """
        TABLE_BORDER = "-----+------------+----------+----------"
        try:
            with open(logfile) as fid:
                for line in fid:
                    if TABLE_BORDER in line:
                        break

                score_lines = list(takewhile(lambda line: "Writing" not in line, fid))
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
    def validate_metadata(metadata: VinaMetadata):
        if shutil.which(metadata.software.value) is None:
            raise MissingExecutableError(
                f'Could not find "{metadata.software.value}" on PATH! '
                "See https://github.com/coleygroup/pyscreener/tree/refactor#adding-an-executable-to-your-path for more information."
            )