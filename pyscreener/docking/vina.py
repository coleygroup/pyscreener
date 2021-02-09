from functools import partial
from itertools import chain, takewhile
import os
from pathlib import Path
import re
import subprocess as sp
import sys
import timeit
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from openbabel import pybel
import ray
from tqdm import tqdm

from pyscreener.preprocessing import autobox
from pyscreener.docking import Screener

class Vina(Screener):
    """A wrapper around Vina-type docking software to perform docking runs via
    python function calls

    Parameters
    ----------
    receptors : Optional[List[str]], default=None
        the filepath(s) of receptors to prepare for docking
    pdbids : Optiona[List[str]], default=None
        a list of PDB IDs corresponding to receptors to prepare for DOCKing.
    center : Optional[Tuple[float, float, float]], default=None
        the center of the docking box (if known)
    size : Tuple[float, float, float], default=(10., 10., 10.)
        the x-, y-, and z-radii of the docking box
    ncpu : int, default=1
        the number of cpu cores over which to parallelize the docking run
    extra : Optional[List[str]], default=None
        additional command line arguments that will be passed to each docking
        run
    docked_ligand_file : Optional[str], default=None
        the filepath of a PDB file containing the coordinates of a docked ligand
    buffer : float, default=10.
        the amount of buffer space to be added around the docked ligand when
        selecting spheres and when constructing the docking box if 
        enclose_spheres is True
    score_mode : str, default='best'
        the mode used to calculate a score for an individual docking run given
        multiple output scored conformations
    repeats : int, default=1
        the number of times each docking run should be repeated
    receptor_score_mode : str
        the mode used to calculate an overall score for a single receptor
        given repeated docking runs against that receptor
    ensemble_score_mode : str
        the mode used to calculate an overall score for an ensemble of receptors
        given multiple receptors in an ensemble
    distributed : bool
        True if the computation will parallelized over a distributed setup.
        False if the computation will parallelized over a local setup
    num_workers : int
        the number of worker processes to initialize when
        distributing computation
    path : os.PathLike
        the path under which input and output folders will be placed
    verbose : int
        the level of output this Screener should output
    
    Attributes
    ----------
    software : str
    center : Optional[Tuple[float, float, float]]
    size : Tuple[float, float, float]
    ncpu : int
    extra : Optional[List[str]]
    """
    def __init__(self, software: str, receptors: Optional[List[str]] = None,
                 pdbids: Optional[List[str]] = None,
                 center: Optional[Tuple[float, float, float]] = None,
                 size: Tuple = (10., 10., 10.), ncpu: int = 1,
                 extra: Optional[List[str]] = None,
                 docked_ligand_file: Optional[str] = None, buffer: float = 10.,
                 score_mode: str = 'best', repeats: int = 1,
                 receptor_score_mode: str = 'best', 
                 ensemble_score_mode: str = 'best',
                 distributed: bool = False, num_workers: int = -1,
                 path: str = '.', verbose: int = 0, **kwargs):
        if software not in ('vina', 'qvina', 'smina', 'psovina'):
            raise ValueError(f'Unrecognized docking software: "{software}"')
        if center is None and docked_ligand_file is None:
            raise ValueError(
                'Args "center" and "docked_ligand_file" were None!')
        if center is None:
            print('Autoboxing ...', end=' ', flush=True)
            center, size = autobox.docked_ligand(docked_ligand_file, buffer)
            print('Done!')
            s_center = f'({center[0]:0.1f}, {center[1]:0.1f}, {center[2]:0.1f})'
            s_size = f'({size[0]:0.1f}, {size[1]:0.1f}, {size[2]:0.1f})'
            print(f'Autoboxed ligand from "{docked_ligand_file}" with', 
                  f'center={s_center} and size={s_size}', flush=True) 

        self.software = software
        self.center = center
        self.size = size
        self.ncpu = ncpu
        self.extra = extra or []

        super().__init__(receptors=receptors, pdbids=pdbids,
                         repeats=repeats, score_mode=score_mode,
                         receptor_score_mode=receptor_score_mode,
                         ensemble_score_mode=ensemble_score_mode,
                         distributed=distributed,
                         num_workers=num_workers, ncpu=ncpu,
                         path=path, verbose=verbose, **kwargs)

    def __call__(self, *args, **kwargs):
        return self.dock(*args, **kwargs)

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
        receptor_pdbqt = str(Path(receptor).with_suffix('.pdbqt'))
        argv = ['prepare_receptor', '-r', receptor, '-o', receptor_pdbqt]
        try:
            sp.run(argv, stderr=sp.PIPE, check=True)
        except sp.SubprocessError:
            print(f'ERROR: failed to convert {receptor}', file=sys.stderr)
            return None

        return receptor_pdbqt

    @staticmethod
    def prepare_from_smi(smi: str, name: str = 'ligand',
                         path: str = '.', **kwargs) -> Tuple[str, str]:
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
        if not path.is_dir():
            path.mkdir(parents=True, exist_ok=True)
        pdbqt = str(path / f'{name}.pdbqt')

        mol = pybel.readstring(format='smi', string=smi)
        mol.addh()
        mol.make3D()
        mol.calccharges(model='gasteiger')
        mol.write(format='pdbqt', filename=pdbqt,
                  overwrite=True, opt={'h': None})

        return smi, pdbqt
    
    def prepare_and_dock(
        self, smis: Sequence[str], names: Sequence[str]
    ) -> List[List[List[Dict]]]:
        
        @ray.remote(num_cpus=self.ncpu)
        def prepare_and_dock_(smi, name):
            ligand = Vina.prepare_from_smi(smi, name, self.tmp_in)

            return Vina.dock_ligand(
                ligand, software=self.software, receptors=self.receptors,
                center=self.center, size=self.size, ncpu=self.ncpu,
                extra=self.extra, path=self.tmp_out,
                repeats=self.repeats, score_mode=self.score_mode
            )

        refs = list(map(prepare_and_dock_.remote, smis, names))
        ligs_recs_reps = [
            ray.get(r) for r in tqdm(refs, desc='Docking ligands')
        ]
        # with self.Pool(self.distributed, self.num_workers, self.ncpu) as pool:
        #     ligs_recs_reps = pool.map(dock_ligand, ligands, 
        #                                 chunksize=2)
        #     ligs_recs_reps = list(tqdm(
        #         ligs_recs_reps, total=len(ligands),
        #         desc='Docking', unit='ligand')
        #     )

        return ligs_recs_reps

    @staticmethod
    def prepare_from_file(filepath: str, use_3d: bool = False,
                          name: Optional[str] = None, path: str = '.', 
                          **kwargs) -> Tuple:
        """Convert a single ligand to the appropriate input format

        Parameters
        ----------
        filename : str
            the name of the file containing the ligand
        use_3d : bool, default=False
            whether to use the 3D information in the input file (if possible)
        name : Optional[str], default=None
            the name of the ligand. If None, use the stem of the input file
        path : str, default='.'
            the path under which the output .pdbqt file should be written
        **kwargs
            additional and unused keyword arguments

        Returns
        -------
        List[Tuple]
            a list of tuples of the SMILES string and the filepath of
            the corresponding input file
        """
        name = name or Path(filepath).stem

        ret = sp.run(['obabel', filepath, '-osmi'], stdout=sp.PIPE, check=True)
        lines = ret.stdout.decode('utf-8').splitlines()
        smis = [line.split()[0] for line in lines]

        if not use_3d:
            ligands = [
                Vina.prepare_from_smi(smi, f'{name}_{i}', path) 
                for i, smi in enumerate(smis)
            ]
            return [lig for lig in ligands if lig]
        
        path = Path(path)
        if not path.is_dir():
            path.mkdir(parents=True, exist_ok=True)

        pdbqt = f'{path}/{name}_.pdbqt'
        argv = ['obabel', filepath, '-opdbqt', '-O', pdbqt, '-m']
        ret = sp.run(argv, check=False, stderr=sp.PIPE)
        
        try:
            ret.check_returncode()
        except sp.SubprocessError:
            return None

        n_mols = 0
        stderr = ret.stderr.decode('utf-8')
        for line in stderr.splitlines():
            if 'converted' not in line:
                continue
            n_mols = int(line.split()[0])

        # have to think about some molecules failing and
        # how that affects numbering
        pdbqts = [f'{path}/{name}_{i}.pdbqt' for i in range(1, n_mols)]

        return list(zip(smis, pdbqts))

    def run_docking(self, ligands: Sequence[Tuple[str, str]]
                   ) -> List[List[List[Dict]]]:
        # dock_ligand = partial(
        #     Vina.dock_ligand,
        #     software=self.software, receptors=self.receptors,
        #     center=self.center, size=self.size, ncpu=self.ncpu,
        #     extra=self.extra, path=self.out_path,
        #     repeats=self.repeats, score_mode=self.score_mode
        # )

        @ray.remote(num_cpus=self.ncpu)
        def dock_ligand_(ligand):
            return Vina.dock_ligand(
                ligand, software=self.software, receptors=self.receptors,
                center=self.center, size=self.size, ncpu=self.ncpu,
                extra=self.extra, path=self.out_path,
                repeats=self.repeats, score_mode=self.score_mode
            )

        refs = list(map(dock_ligand_.remote, ligands))
        ligs_recs_reps = [
            ray.get(r) for r in tqdm(refs, desc='Docking ligands')
        ]
        # with self.Pool(self.distributed, self.num_workers, self.ncpu) as pool:
        #     ligs_recs_reps = pool.map(dock_ligand, ligands, 
        #                                 chunksize=2)
        #     ligs_recs_reps = list(tqdm(
        #         ligs_recs_reps, total=len(ligands),
        #         desc='Docking', unit='ligand')
        #     )

        return ligs_recs_reps

    @staticmethod
    def dock_ligand(ligand: Tuple[str, str], software: str,
                    receptors: List[str],
                    center: Tuple[float, float, float],
                    size: Tuple[int, int, int] = (10, 10, 10), ncpu: int = 1, 
                    path: str = '.', extra: Optional[List[str]] = None,
                    repeats: int = 1, score_mode: str = 'best'
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

        smi, pdbqt = ligand

        p_pdbqt = Path(pdbqt)
        ligand_name = p_pdbqt.stem

        ensemble_rowss = []
        for receptor in receptors:
            repeat_rows = []
            for repeat in range(repeats):
                name = f'{Path(receptor).stem}_{ligand_name}_{repeat}'

                argv, p_out, p_log = Vina.build_argv(
                    software=software, receptor=receptor, ligand=pdbqt, 
                    name=name, center=center, size=size, ncpu=ncpu, 
                    extra=extra, path=path
                )

                ret = sp.run(argv, stdout=sp.PIPE, stderr=sp.PIPE)
                try:
                    ret.check_returncode()
                except sp.SubprocessError:
                    print(f'ERROR: docking failed. argv: {argv}',
                          file=sys.stderr)
                    print(f'Message: {ret.stderr.decode("utf-8")}',
                          file=sys.stderr)

                repeat_rows.append({
                    'smiles': smi,
                    'name': ligand_name,
                    'node_id': re.sub('[:,.]', '', ray.state.current_node_id()),
                    # 'in': p_pdbqt,
                    # 'out': p_out,
                    # 'log': p_log,
                    'score': Vina.parse_log_file(p_log, score_mode)
                })

            ensemble_rowss.append(repeat_rows)

        return ensemble_rowss

    @staticmethod
    def build_argv(software: str, receptor: str, ligand: str,
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
        if software not in {'vina', 'smina', 'psovina', 'qvina'}:
            raise ValueError(f'Invalid docking program: "{software}"')

        path = Path(path)
        if not path.is_dir():
            path.mkdir(parents=True, exist_ok=True)

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
    def parse_log_file(log_file: str,
                       score_mode: str = 'best') -> Optional[float]:
        """parse a Vina-type log file to calculate the overall ligand score
        from a single docking run

        Parameters
        ----------
        log_file : str
            the path to a Vina-type log file
        score_mode : str (Default = 'best')
            the method by which to calculate a score from multiple scored
            conformations. Choices include: 'best', 'average', and 'boltzmann'.
            See Screener.calc_score for further explanation of these choices.
        
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
                score = Screener.calc_score(scores, score_mode)
        except OSError:
            score = None
        
        return score

    @staticmethod
    def parse_ligand_results(ligand_results: List[List[Dict]],
                             score_mode: str = 'best'):
        for receptor_results in ligand_results:
            for repeat_result in receptor_results:
                score = Vina.parse_log_file(repeat_result['log'], score_mode)
                repeat_result['score'] = score

                # repeat_result['score'] = score
                # p_in = repeat_result['in']
                # repeat_result['in'] = Path(p_in.parent.name) / p_in.name
                # p_out = repeat_result['out']
                # repeat_result['out'] = Path(p_out.parent.name) / p_out.name
                # p_log = repeat_result['log']
                # repeat_result['log'] = Path(p_log.parent.name) / p_log.name
        return ligand_results