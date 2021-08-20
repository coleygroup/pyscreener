from dataclasses import dataclass, field
from itertools import takewhile
from math import ceil, log10
from pathlib import Path
import re
import shlex
import subprocess as sp
import sys
from typing import Dict, List, Mapping, Optional, Tuple, Union

# from openbabel import pybel
import ray

from pyscreener import utils
from pyscreener.exceptions import NotSimulatedError
from pyscreener.docking.simulation import DockingSimulation

@dataclass(repr=False, eq=False)
class VinaCalculationData:
    ligand: str
    receptor: str
    software: str
    center: Tuple[float, float, float]
    size: Tuple[float, float, float] = (10., 10., 10.)
    ncpu: int = 1
    extra: Optional[str] = None
    name: str = None
    in_path: Union[str, Path] = '.'
    out_path: Union[str, Path] = '.'
    prepared_ligand: Optional[Union[str, Path]] = None,
    prepared_receptor: Optional[Union[str, Path]] = None

    def __post_init__(self):
        # if self.software not in ('vina', 'smina', 'psovina', 'qvina'):
        #     raise ValueError(f'Invalid docking software: "{self.software}"')
        self.extra = shlex.split(self.extra) if self.extra else []
        self.in_path = Path(self.in_path)
        self.out_path = Path(self.out_path)

class VinaCalculation(DockingSimulation):
    def __init__(
        self, ligand: Union[str, Path], receptor: Union[str, Path],
        software: str, center: Tuple[float, float, float],
        size: Tuple[int, int, int] = (10, 10, 10), ncpu: int = 1,
        extra: Optional[str] = None, name: Optional[str] = None,
        in_path: Union[str, Path] = './inputs',
        out_path: Union[str, Path] = './outputs',
        prepared_ligand: Optional[Union[str, Path]] = None,
        prepared_receptor: Optional[Union[str, Path]] = None
    ):
        if software not in ('vina', 'smina', 'psovina', 'qvina'):
            raise ValueError(f'Invalid docking software: "{software}"')

        self.software = software
        self.center = center
        self.size = size
        self.ncpu = ncpu
        self.extra = shlex.split(extra) if extra else []

        self.name = name
        self.in_path = Path(in_path)
        self.out_path = Path(out_path)

        self.__prepared_ligand = Path(prepared_ligand) if prepared_ligand else None
        self.__prepared_receptor = Path(prepared_receptor) if prepared_receptor else None

        self.__result = None
        
        super().__init__('.')
    
    @property
    def receptor(self) -> Path:
        return self.__prepared_receptor

    @property
    def smi(self) -> str:
        return self.smi

    @property
    def ligand(self) -> Path:
        return self.__prepared_ligand

    def score(self) -> Optional[float]:
        try:
            return self.result['score']
        except TypeError:
            raise NotSimulatedError(
                'Simulation has not been run!'
            )

    def result(self) -> Dict:
        return self.__result

    def prepare(self):
        self.__prepared_receptor = self.prepare_receptor(self.receptor)
        self.__prepared_ligand = self.prepare_ligand(self.ligand)

    def run(self) -> float:
        argv, _, log = self.build_argv()

        ret = sp.run(argv, stdout=sp.PIPE, stderr=sp.PIPE)
        try:
            ret.check_returncode()
        except sp.SubprocessError:
            print(f'ERROR: docking failed. argv: {argv}', file=sys.stderr)
            print(f'Message: {ret.stderr.decode("utf-8")}', file=sys.stderr)

        self.__result = {
            'smiles': self.smi,
            'name': self.ligand.stem,
            'node_id': re.sub('[:,.]', '', ray.state.current_node_id()),
            'score': self.parse_log_file(log)
        }

        return self.score

    def dock_ligand(
        self, ligand: Tuple[str, str], software: str,
        receptors: List[str],
        center: Tuple[float, float, float],
        size: Tuple[int, int, int] = (10, 10, 10), ncpu: int = 1, 
        path: str = '.', extra: Optional[List[str]] = None,
        repeats: int = 1, score_mode: str = 'best', k: int = 1
    ) -> List[List[Dict]]:
        """Dock the given ligand using the specified vina-type docking program 
        and parameters into the ensemble of receptors repeatedly
        
        Parameters
        ----------
        software : str
            the docking program to run
        ligand : Ligand
            a tuple containing a ligand's SMILES string and associated docking
            input file
        receptors : List[str]
            the filesnames of PDBQT files corresponding to 
            various receptor poses
        center : Tuple[float, float, float]
            the x-, y-, and z-coordinates, respectively, of the search box 
            center
        size : Tuple[int, int, int], default=(10, 10, 10)
            the x, y, and z-radii, respectively, of the search box
        path : string, default='.'
            the path under which both the log and out files should be written to
        ncpu : int, default=1
            the number of cores to allocate to the docking program
        repeats : int, default=1
            the number of times to repeat a docking run
        score_mode : str
            the mode used to calculate a score for an individual docking run 
            given multiple output scored conformations
        k : int, default=1
            the number of top scores to average if using "top-k"score mode
        Returns
        -------
        ensemble_rowss : List[List[Dict]]
            an MxO list of dictionaries where each dictionary is a record of an 
            individual docking run and:

            * M is the number of receptors each ligand is docked against
            * O is the number of times each docking run is repeated.

            Each dictionary contains the following keys:

            * smiles: the ligand's SMILES string
            * name: the name of the ligand
            * in: the filename of the input ligand file
            * out: the filename of the output docked ligand file
            * log: the filename of the output log file
            * score: the ligand's docking score
        """
        raise NotImplementedError

    def prepare_receptor(self, receptor: str) -> Optional[str]:
        """Prepare a receptor PDBQT file from its input file

        Parameters
        ----------
        receptor : str
            the filename of a file containing a receptor

        Returns
        -------
        receptor_pdbqt : Optional[str]
            the filepath of the resulting PDBQT file.
            None if preparation failed
        """
        receptor_pdbqt = Path(receptor).with_suffix('.pdbqt').name
        receptor_pdbqt = str(self.receptors_dir / receptor_pdbqt)

        argv = ['prepare_receptor', '-r', receptor, '-o', receptor_pdbqt]
        try:
            sp.run(argv, stderr=sp.PIPE, check=True)
        except sp.SubprocessError:
            print(f'ERROR: failed to convert "{receptor}"', file=sys.stderr)
            return None

        return receptor_pdbqt
    
    @staticmethod
    def prepare_from_smi(
        smi: str, name: str = 'ligand', path: str = '.'
    ) -> Tuple[str, str]:
        """Prepare an input ligand file from the ligand's SMILES string

        Parameters
        ----------
        smi : str
            the SMILES string of the ligand
        name : Optional[str], default='ligand'
            the name of the ligand.
        path : str, default='.'
            the path under which the output PDBQT file should be written

        Returns
        -------
        smi : str
            the ligand's SMILES string
        pdbqt : str
            the filepath of the coresponding input file
        """
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)

        pdbqt = str(path / f'{name}.pdbqt')

        mol = pybel.readstring(format='smi', string=smi)
        
        mol.make3D()
        mol.addh()

        try:
            mol.calccharges(model='gasteiger')
        except Exception:
            pass
        mol.write(format='pdbqt', filename=pdbqt,
                overwrite=True, opt={'h': None})

        return smi, pdbqt

    @staticmethod
    def prepare_from_file(
        filename: str, path: str = '.', offset: int = 0
    ) -> List[Tuple[str, str]]:
        """Convert the molecules contained in the arbitrary chemical file to
        indidvidual PDBQT files with their same geometry

        Parameters
        ----------
        filename : str
            the name of the file containing the molecules
        path : str, default='.'
            the path under which the output PDBQT files should be written
        offset : int, default=0
            the numbering offset for ligand naming

        Returns
        -------
        List[Tuple[str, str]]
            a list of tuples of the SMILES string and the filepath of
            the corresponding PDBQT file
        """
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)

        fmt = Path(filename).suffix.strip('.')
        mols = list(pybel.readfile(fmt, filename))
        width = ceil(log10(len(mols) + offset)) + 1

        smis = []
        pdbqts = []
        for i, mol in enumerate(mols):
            if mol.title:
                name = mol.title
            else:
                name = f'ligand_{i+offset:0{width}}'
            
            smi = mol.write()
            pdbqt = str(path / f'{name}.pdbqt')
            mol.addh()
            mol.make3D()
            mol.calccharges(model='gasteiger')
            mol.write(format='pdbqt', filename=pdbqt,
                      overwrite=True, opt={'h': None})

            smis.append(smi)
            pdbqts.append(pdbqt)

        return list(zip(smis, pdbqts))
    
    def build_argv(
        self, name: Optional[str] = None,
        path: str = '.',
    ) -> Tuple[List[str], str, str]:
        """Builds the argument vector to run a vina-type docking simulation

        Parameters
        ----------
        name : string, default=<receptor>_<ligand>)
            the base name to use for both the log and out files
        path : string, default='.'
            the path under which both the log and out files should be written

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
        name = self.name or (
            self.software + '_' + self.receptor.stem+ '_' +self.ligand.stem
        )

        out = path / f'{name}_out.pdbqt'
        log = path / f'{name}.log'
        
        argv = [
            self.software,
            f'--receptor={self.receptor}', f'--ligand={self.ligand}',
            f'--center_x={self.center[0]}',
            f'--center_y={self.center[1]}',
            f'--center_z={self.center[2]}',
            f'--size_x={self.size[0]}',
            f'--size_y={self.size[1]}',
            f'--size_z={self.size[2]}',
            f'--cpu={self.ncpu}', f'--out={out}', f'--log={log}', *self.extra
        ]

        return argv, out, log

    def parse_log_file(self, log_file: str) -> Optional[float]:
        """parse a Vina-type log file to calculate the overall ligand score
        from a single docking run

        Parameters
        ----------
        log_file : str
            the path of a Vina-type log file

        Returns
        -------
        Optional[float]
            the score parsed from the log file. None if no score was parsed
        """
        TABLE_BORDER = '-----+------------+----------+----------'
        try:
            with open(log_file) as fid:
                for line in fid:
                    if TABLE_BORDER in line:
                        break

                score_lines = takewhile(
                    lambda line: 'Writing' not in line, fid
                )
                scores = [float(line.split()[1]) for line in score_lines]

            if len(scores) == 0:
                score = None
            else:
                score = utils.calc_score(scores, self.score_mode, self.k)
        except OSError:
            score = None
        
        return score