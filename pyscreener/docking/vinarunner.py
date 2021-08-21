from itertools import takewhile
from math import ceil, log10
from pathlib import Path
import re
import subprocess as sp
import sys
from typing import Dict, List, Optional, Tuple, Union

from openbabel import pybel
import ray

from pyscreener import utils

class VinaRunner(object):
    @staticmethod
    def prepare_receptor(
        receptor: str, path: Union[str, Path], name: Optional[str] = None
    ) -> Optional[str]:
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
        name = name or Path(receptor).stem
        receptor_pdbqt = (Path(path) / name).with_suffix('.pdbqt')
        # receptor_pdbqt = name or Path(receptor).with_suffix('.pdbqt').name
        # receptor_pdbqt = str(self.receptors_dir / receptor_pdbqt)

        argv = ['prepare_receptor', '-r', receptor, '-o', receptor_pdbqt]
        try:
            sp.run(argv, stderr=sp.PIPE, check=True)
        except sp.SubprocessError:
            print(f'ERROR: failed to convert "{receptor}"', file=sys.stderr)
            return None

        return receptor_pdbqt

    @staticmethod
    def prepare_from_smi(
        smi: str, path: str = '.', name: str = 'ligand', 
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
        **kwargs
            additional and unused keyword arguments

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
        filename: str, path: str = '.',  offset: int = 0
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
        **kwargs
            additional and unused keyword arguments

        Returns
        -------
        List[Tuple[str, str]]
            a list of tuples of the SMILES string and the filepath of
            the corresponding PDBQT file
        """
        path = Path(path)
        if not path.is_dir():
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
            # mol.make3D()
            mol.calccharges(model='gasteiger')
            mol.write(format='pdbqt', filename=pdbqt,
                      overwrite=True, opt={'h': None})

            smis.append(smi)
            pdbqts.append(pdbqt)

        return list(zip(smis, pdbqts))
    
    @staticmethod
    def dock_ligand(ligand: Tuple[str, str], receptor: str,
                    software: str,
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
        if repeats <= 0:
            raise ValueError(f'Repeats must be greater than 0! ({repeats})')

        path = Path(path)
        # path.mkdir(parents=True, exist_ok=True)

        smi, pdbqt = ligand

        p_pdbqt = Path(pdbqt)
        ligand_name = p_pdbqt.stem

        name = f'{Path(receptor).stem}_{ligand_name}'

        argv, _, log = VinaRunner.build_argv(
            ligand=pdbqt, receptor=receptor, software=software,
            name=name, center=center, size=size, ncpu=ncpu, 
            extra=extra, path=path
        )

        ret = sp.run(argv, stdout=sp.PIPE, stderr=sp.PIPE)
        try:
            ret.check_returncode()
        except sp.SubprocessError:
            print(
                f'ERROR: docking failed. argv: {argv}', file=sys.stderr
            )
            print(
                f'Message: {ret.stderr.decode("utf-8")}', file=sys.stderr
            )

        return {
            'smiles': smi,
            'name': ligand_name,
            'node_id': re.sub('[:,.]', '', ray.state.current_node_id()),
            'score': VinaRunner.parse_log_file(log, score_mode, k)
        }

    @staticmethod
    def build_argv(
        ligand: str, receptor: str, software: str, 
        center: Tuple[float, float, float],
        size: Tuple[int, int, int] = (10, 10, 10),
        ncpu: int = 1, name: Optional[str] = None,
        path: str = '.', extra: Optional[List[str]] = None
    ) -> Tuple[List[str], str, str]:
        """Builds the argument vector to run a vina-type docking program

        Parameters
        ----------
        software : str
            the name of the docking program to run
        receptor : str
            the filename of the input receptor file
        ligand : str
            the filename of the input ligand file
        center : Tuple[float, float, float]
            the coordinates (x,y,z) of the center of the vina search box
        size : Tuple[int, int, int], default=(10, 10, 10)
            the size of the vina search box in angstroms for the x, y, and z-
            dimensions, respectively
        ncpu : int, default=1
            the number of cores to allocate to the docking program
        name : string, default=<receptor>_<ligand>)
            the base name to use for both the log and out files
        path : string, default='.'
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
        if software not in ('vina', 'smina', 'psovina', 'qvina'):
            raise ValueError(f'Invalid docking program: "{software}"')

        name = name or (Path(receptor).stem+'_'+Path(ligand).stem)
        extra = extra or []

        out = path / f'{software}_{name}_out.pdbqt'
        log = path / f'{software}_{name}.log'
        
        argv = [
            software, f'--receptor={receptor}', f'--ligand={ligand}',
            f'--center_x={center[0]}',
            f'--center_y={center[1]}',
            f'--center_z={center[2]}',
            f'--size_x={size[0]}', f'--size_y={size[1]}', f'--size_z={size[2]}',
            f'--cpu={ncpu}', f'--out={out}', f'--log={log}', *extra
        ]

        return argv, out, log

    @staticmethod
    def parse_log_file(
        log_file: str, score_mode: str = 'best', k: int = 1
    ) -> Optional[float]:
        """parse a Vina-type log file to calculate the overall ligand score
        from a single docking run

        Parameters
        ----------
        log_file : str
            the path to a Vina-type log file
        score_mode : str, default='best'
            the method by which to calculate a score from multiple scored
            conformations. Choices include: 'best', 'average', and 'boltzmann'.
            See Screener.calc_score for further explanation of these choices.
        k : int, default=1
            the number of top scores to average if using "top-k" score mode
        Returns
        -------
        Optional[float]
            the score parsed from the log file. None if no score was parsed
        """
        # vina-type log files have scoring information between this 
        # table border and the line: "Writing output ... done."
        TABLE_BORDER = '-----+------------+----------+----------'
        try:
            with open(log_file) as fid:
                for line in fid:
                    if TABLE_BORDER in line:
                        break

                score_lines = takewhile(
                    lambda line: 'Writing' not in line, fid)
                scores = [float(line.split()[1])
                            for line in score_lines]

            if len(scores) == 0:
                score = None
            else:
                score = utils.calc_score(scores, score_mode, k)
        except OSError:
            score = None
        
        return score