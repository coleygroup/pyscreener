from abc import ABC, abstractmethod
import csv
import datetime
import glob
from itertools import chain
from math import ceil, exp, log10
import os
from pathlib import Path
import re
import shutil
import tarfile
import tempfile
import timeit
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union

from openbabel import pybel
import ray
from rdkit import Chem

from pyscreener.preprocessing import pdbfix
from pyscreener import utils

pybel.ob.obErrorLog.SetOutputLevel(0)

class Screener(ABC):
    """A Screener conducts virtual screens against an ensemble of receptors.

    Classes that implement the Screener interface are responsible for
    defining the following methods:

    * prepare_receptor
    * prepare_from_smi
    * prepare_from_file
    * prepare_and_dock

    NOTE: This is an abstract base class and cannot be instantiated.

    Parameters
    ----------
    receptors : List[str]
        the filepath(s) of receptors to prepare for docking
    pdbids : List[str]
        a list of PDB IDs corresponding to receptors to prepare for DOCKing.
    repeats : int, default=1
        the number of times each docking run will be repeated
    score_mode : str, default='best'
        the mode used to calculate a score for an individual docking run given
        multiple output scored conformations
    receptor_score_mode : str, default='best'
        the mode used to calculate an overall score for a single receptor
        given repeated docking runs against that receptor
    ensemble_score_mode : str, default='best'
        the mode used to calculate an overall score for an ensemble of receptors
        given multiple receptors in an ensemble
    distributed : bool, default=False
        True if the computation will parallelized over a distributed setup.
        False if the computation will parallelized over a local setup
    num_workers : int, default=-1
        the number of worker processes to initialize when
        distributing computation
    ncpu : int, default=1
        the number of cores allocated to each worker process
    path : os.PathLike, default='.'
        the directory under which input and output folders will be placed
    tmp_dir : os.PathLike, default=tempfile.gettempdir()
        the temp directory under which input and output will be placed
    verbose : int
        the level of output this Screener should output

    Attributes
    ----------
    repeats : int
    score_mode : str
    receptor_score_mode : str
    ensemble_score_mode : str
    distributed : bool
    num_workers : int, default=-1
    ncpu : int, default=1
    verbose : int, default=0
    **kwargs
        additional and unused keyword arguments
    """
    def __init__(self, receptors: Optional[Sequence[str]] = None,
                 pdbids: Optional[Sequence[str]] = None,
                 repeats: int = 1, score_mode: str = 'best',
                 receptor_score_mode: str = 'best', 
                 ensemble_score_mode: str = 'best',
                 ncpu: int = 1,
                 path: str = '.', tmp_dir = tempfile.gettempdir(),
                 verbose: int = 0, **kwargs):
        self.path = Path(path)
        self.tmp_dir = tmp_dir

        receptors = receptors or []
        if pdbids:
            receptors.extend([
                pdbfix.pdbfix(pdbid=pdbid, path=self.receptors_dir)
                for pdbid in pdbids
            ])
        if len(receptors) == 0:
            raise ValueError('No receptors or PDBIDs provided!')

        self.receptors = receptors
        self.repeats = repeats
        self.score_mode = score_mode
        self.receptor_score_mode = receptor_score_mode
        self.ensemble_score_mode = ensemble_score_mode
        
        self.ncpu = ncpu

        self.verbose = verbose

        self.num_docked_ligands = 0
        
        if not ray.is_initialized():
            try:
                print('No Ray cluster started! attempting to autoconnect...')
                ray.init(address='auto')
            except ConnectionError:
                print('No Ray cluster detected! Starting local cluster...')
                ray.init()

    def __len__(self) -> int:
        """The number of ligands this screener has simulated"""
        return self.num_docked_ligands

    def __call__(self, *args, **kwargs) -> Dict[str, Optional[float]]:
        return self.dock(*args, **kwargs)
    
    @property
    def path(self) -> Path:
        return self.__path
    
    @path.setter
    def path(self, path: os.PathLike):
        path = Path(path)
        receptors_dir = path / 'receptors'
        if not receptors_dir.is_dir():
            receptors_dir.mkdir(parents=True, exist_ok=True)
        
        self.__path = path
        self.__receptors_dir = receptors_dir

    @property
    def receptors_dir(self) -> os.PathLike:
        return self.__receptors_dir

    @property
    def tmp_dir(self) -> Path:
        """the Screener's temp directory"""
        return self.__tmp_dir
        
    @tmp_dir.setter
    def tmp_dir(self, path: str):
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        path = Path(path) / 'pyscreener' / f'session_{timestamp}'
        if not path.is_dir():
            path.mkdir(parents=True)
        self.__tmp_dir = path
        self.tmp_in = path / 'inputs'
        self.tmp_out = path / 'outputs'

    @property
    def tmp_in(self) -> Path:
        return self.__tmp_in

    @tmp_in.setter
    def tmp_in(self, path: str):
        path = Path(path)
        if not path.is_dir():
            path.mkdir(parents=True)
        self.__tmp_in = path
    
    @property
    def tmp_out(self) -> Path:
        return self.__tmp_out

    @tmp_out.setter
    def tmp_out(self, path: str):
        path = Path(path)
        if not path.is_dir():
            path.mkdir(parents=True)
        self.__tmp_out = path

    def collect_files(self, out_path: Optional[str]=None):
        """Collect all the files from the local disks of the respective nodes

        For I/O purposes, input and output files for each simulation are 
        created on the local disk of each node. If these files are desired at 
        the end, they must be copied over from the node's local file system to 
        the final destination.
        
        This is achieved by creating a gzipped tar file of the temp directory 
        (the one that contains all of the input and output files for 
        simulations conducted on that node) and moving these tar files under 
        the desired path. Each tar file is named according the node ID from 
        which it originates.

        This function should ideally only be called once during the lifetime
        of a Screener because it is slow and early calls will yield nothing
        over a single, final call.

        Parameters
        ----------
        out_path : Optional[str], defualt=None
            the path under which the tar files should be collected to. If None,
            use self.path
        """
        out_path = out_path or self.path
        out_path = Path(out_path)
        if not out_path.is_dir():
            out_path.mkdir(parents=True)

        refs = []
        for node in ray.nodes():    # run on all nodes
            address = node["NodeManagerAddress"]
            @ray.remote(resources={f'node:{address}': 0.1})
            def zip_and_move_tmp():
                output_id = re.sub('[:,.]', '', ray.state.current_node_id())
                tmp_tar = (self.tmp_dir / output_id).with_suffix('.tar.gz')

                with tarfile.open(tmp_tar, 'w:gz') as tar:
                    tar.add(self.tmp_dir, arcname=output_id)

                shutil.copy(str(tmp_tar), str(out_path))

            refs.append(zip_and_move_tmp.remote())
        ray.wait(refs)
        
    @property
    def receptors(self) -> List[str]:
        return self.__receptors

    @receptors.setter
    def receptors(self, receptors):
        receptors = [self.prepare_receptor(rec) for rec in receptors]
        receptors = [rec for rec in receptors if rec is not None]
        if len(receptors) == 0:
            raise RuntimeError('Preparation failed for all receptors!')
        # self.__receptors = receptors

        refs = []
        for node in ray.nodes():    # run on all nodes
            address = node["NodeManagerAddress"]

            @ray.remote(resources={f'node:{address}': 0.1})
            def copy_receptors():
                self.tmp_dir.mkdir(parents=True, exist_ok=True)
                return self.copy_receptors(receptors)
            # def copy_receptors():
            #     

            #     copied_receptors = []
            #     for rec in receptors:
            #         if isinstance(rec, str):
            #             copied_receptor = shutil.copy(rec, str(self.tmp_dir))
            #         else:
            #             rec = list(chain(*(glob.glob(f'{r}*') for r in rec)))
            #             copied_receptor = tuple([
            #                 shutil.copy(r, str(self.tmp_dir)) for r in rec
            #             ])
            #         copied_receptors.append(copied_receptor)

            #     return copied_receptors
                # return [
                #     shutil.copy(rec, str(self.tmp_dir)) for rec in receptors
                # ]

            refs.append(copy_receptors.remote())
        ray.wait(refs)

        # print(ray.get(refs[0]))
        self.__receptors = ray.get(refs[0])

        # def copy_receptors():
        #     return [shutil.copy(rec, str(self.tmp_dir)) for rec in receptors]
        # rs = utils.run_on_all_nodes(copy_receptors)
        # print(rs[0])
        

    @abstractmethod
    def prepare_and_dock(
        self, smis: Sequence[str], names: Sequence[str]
    ) -> List[List[List[Dict]]]:
        """Run the docking simulations for the input molecules
        
        Parameters
        ----------
        smis : Sequence[str]
            the SMILES strings of the molecules to dock
        names: Sequence[str]
            a parallel list containing the ID or name of each molecule
        
        Returns
        -------
        List[List[List[Dict]]]
            an NxMxO list of dictionaries where each individual dictionary is a 
            record of an individual docking run and
            
            * N is the number of molecules that were docked
            * M is the number of receptors in the ensemble
            * O is the number of times each docking run was repeated
        """

    @abstractmethod
    def prepare_receptor(self, *args, **kwargs) -> Optional[str]:
        """Prepare a receptor input file for the docking software and return
        the corresponding filepath or None if preparation failed"""

    @abstractmethod
    def copy_receptors(self, receptors) -> List[str]:
        """Copy the prepared receptors to self.tmp_dir"""

    @staticmethod
    @abstractmethod
    def prepare_from_smi(*args, **kwargs) -> Tuple[str, str]:
        """Prepare a ligand input file from a SMILES string and return both
        the SMILES string and the corresponding filepath"""

    @staticmethod
    @abstractmethod
    def prepare_from_file(*args, **kwargs) -> List[Tuple[str, str]]:
        """Prepare the corresponding ligand input files from an
        arbitrary input file type"""

    def dock(self, *smis_or_files: Iterable,
             full_results: bool = False,
             **kwargs) -> Dict[str, Optional[float]]:
        """dock the ligands contained in sources

        NOTE: the star operator, *, in the function signature.
            If intending to pass multiple filepaths as an iterable, first 
            unpack the iterable in the function call by prepending a *.
            If passing multiple SMILES strings, either option is acceptable,
            but it is marginally more efficient to NOT unpack the iterable.

        Parameters
        ----------
        smis_or_files: Iterable
            an iterable of ligand sources, where each ligand source may be
            one of the following:
            
            * an arbitrary format file containing molecules,
            * a list of SMILES strings
            * a single SMILES string
        full_results : bool, default=False
            whether the full dataframe of every single run should be returned
        **kwargs
            keyword arguments to pass to the appropriate prepare_from_*
            function(s)

        Returns
        -------
        d_smi_score : Dict[str, Optional[float]]
            a dictionary mapping SMILES string to the best score among the
            corresponding ligands. (None if all corresponding ligands
            failed failed to dock)
        records : List[Dict]
            a list of dictionaries containing the record of every single
            docking run performed. Each dictionary contains the following keys:
            
            * 'smiles': the ligand's SMILES string
            * 'name': the name of the ligand
            * 'node_id': the ID of the node that ran the docking simulations and
                    the ID of the TAR file containing all of the
                    corresponding files if collect_files() is called
            * 'score': the ligand's docking score
        """
        recordsss = self.dock_ensemble(*smis_or_files, **kwargs)

        smis_scores = []
        for ligand_results in recordsss:
            smi = ligand_results[0][0]['smiles']
            score = self.calc_ligand_score(
                ligand_results, self.receptor_score_mode,
                self.ensemble_score_mode
            )
            smis_scores.append((smi, score))

        d_smi_score = {}
        for smi_score in smis_scores:
            smi, score = smi_score
            if smi not in d_smi_score:
                d_smi_score[smi] = score
            elif score is None:
                continue
            else:
                curr_score = d_smi_score[smi]
                if curr_score is None:
                    d_smi_score[smi] = score
                else:
                    d_smi_score[smi] = min(d_smi_score[smi], score)

        if full_results:
            return d_smi_score, list(chain(*list(chain(*recordsss))))

        return d_smi_score

    def dock_ensemble(self, *smis_or_files: Iterable,
                      **kwargs) -> List[List[List[Dict]]]:
        """Run the docking program with the ligands contained in *smis_or_files

        NOTE: the star operator, *, in the function signature
            If intending to pass multiple filepaths as an iterable, first 
            unpack the iterable in the function call by prepending a *

        Parameters
        ----------
        smis_or_files: Iterable
            an unpacked iterable of ligand sources, where each ligand 
            source may be one of the following:

            * an arbitrary chemical supply file (e.g., CSV, SMI, SDF, etc.)
            * a list of SMILES strings
            * a single SMILES string
        **kwargs
            keyword arguments to pass to the appropriate prepare_from_*
            function(s)

        Returns
        -------
        recordsss : List[List[List[Dict]]]
            an NxMxO list of dictionaries where each dictionary is a record of 
            an individual docking run and:

            * N is the number of total ligands that were docked
            * M is the number of receptors each ligand was docked against
            * O is the number of times each docking run was repeated.
            
            Each dictionary contains the following keys:

            * 'smiles': the ligand's SMILES string (not canonicalized)
            * 'name': the name of the ligand
            * 'node_id': the ID of the node that ran the docking simulations and
                    the ID of the TAR file containing all of the
                    corresponding files if collect_files() is called
            * 'score': the ligand's docking score

        """
        begin = timeit.default_timer()

        smis = []
        names = []
        for smi_or_file in smis_or_files:
            smis_, names_ = self.get_smis(smi_or_file, **kwargs)
            smis.extend(smis_), names.extend(names_)

        width = ceil(log10(len(smis) + len(self))) + 1
        names = [
            name or f'ligand_{i + len(self):0{width}}'
            for i, name in enumerate(names)
        ]

        assert len(smis) == len(names)

        recordsss = self.prepare_and_dock(smis, names)

        self.num_docked_ligands += len(recordsss)

        total = timeit.default_timer() - begin

        mins, secs = divmod(int(total), 60)
        hrs, mins = divmod(mins, 60)
        if self.verbose > 0 and len(recordsss) > 0:
            print(f'  Time to dock {len(recordsss)} ligands:',
                  f'{hrs:d}h {mins:d}m {secs:d}s ' +
                  f'({total/len(recordsss):0.3f} s/ligand)', flush=True)

        return recordsss

    def get_smis(self, source: Union[str, Sequence[str]],
                 **kwargs) -> Tuple[List[str], List[str]]:
        """Get the SMILES strings and names from an input source file

        Parameters
        ----------
        source : Union[str, Sequence[str]]
            the file to parse, SMILES string, or list of SMILES strings
        **kwargs
            keyword arguments to pass to the appropriate
            get_smis_from_* method

        Returns
        -------
        smis : List[str]
            the SMILES strings
        names : List[Optional[[str]]
            the names of the molecules

        Raises
        ------
        ValueError
            if the file does not exit
        """
        if not isinstance(source, str):
            smis = [smi for smi in source]
            names = [None for _ in smis]

        p_source = Path(source)
        if not p_source.exists():
            smis, names = [source], [None]
        elif p_source.suffix == '.csv':
            smis, names = self.get_smis_from_csv(source, **kwargs)
        elif p_source.suffix in ('.smi', '.sdf'):
            smis, names = self.get_smis_from_supply(source, **kwargs)
        else:
            smis, names = self.get_smis_from_file(source, **kwargs)
            
        return smis, names

    @staticmethod
    def get_smis_from_csv(csv_filename: str, title_line: bool = True,
                          smiles_col: int = 0, name_col: Optional[int] = None,
                          **kwargs) -> Tuple[List[str], List[Optional[str]]]:
        """Prepare the input files corresponding to the SMILES strings
        contained in a CSV file

        Parameters
        ----------
        csv_filename : str
            the filename of the CSV file containing the ligands to convert
        title_line : bool, default=True
            does the CSV file contain a title line?
        smiles_col : int, default=0
            the column containing the SMILES strings
        name_col : Optional[int], default=None
            the column containing the molecule name
        **kwargs
            additional and unused keyword arguments

        Returns
        -------
        smis : List[str]
            the SMILES strings parsed from the CSV
        names : List[Optional[str]]
            the names parsed from the CSV. None's if there were none.
        """
        with open(csv_filename) as fid:
            reader = csv.reader(fid)
            if title_line:
                next(reader)

            if name_col is None:
                smis = [row[smiles_col] for row in reader]
                names = [None for _ in smis]
            else:
                smis, names = zip(*[
                    (row[smiles_col], row[name_col]) for row in reader
                ])
        
        return smis, names

    @staticmethod
    def get_smis_from_supply(supply: str,
                             id_prop_name: Optional[str] = None,
                             **kwargs) -> Tuple[List[str], List[Optional[str]]]:
        """Prepare the input files corresponding to the molecules contained in 
        a molecular supply file

        Parameters
        ----------
        supply : str
            the filename of the SDF or SMI file containing
            the ligands to convert
        id_prop_name : Optional[str]
            the name of the property containing the ID, if one exists
            (e.g., "CatalogID", "Chemspace_ID", "Name", etc...)
        **kwargs
            additional and unused keyword arguments

        Returns
        -------
        smis : List[str]
            the SMILES strings parsed from the supply file
        names : List[Optional[str]]
            the names parsed from the supply file. None's if there were none.
        """
        p_supply = Path(supply)
        if p_supply.suffix == '.sdf':
            mols = Chem.SDMolSupplier(supply)
        elif p_supply.suffix == '.smi':
            mols = Chem.SmilesMolSupplier(supply)
        else:
            raise ValueError(
                f'input file: "{supply}" does not have .sdf or .smi extension')

        smis = []
        names = []

        if id_prop_name:
            for mol in mols:
                if mol is None:
                    continue

                smis.append(Chem.MolToSmiles(mol))
                names.append(mol.GetProp(id_prop_name))
        else:
            for mol in mols:
                if mol is None:
                    continue

                smis.append(Chem.MolToSmiles(mol))
            names = [None for _ in smis]
        return smis, names

    @staticmethod
    def get_smis_from_file(filename: str,
                           **kwargs) -> Tuple[List[str], List[Optional[str]]]:
        """get the SMILES strings and names from an arbitrary chemical 
        file format

        Parameters
        ----------
        filename : str
            the file to parse
        **kwargs
            additional and unused keyword arguments

        Returns
        -------
        smis : List[str]
            the SMILES strings
        names : List[Optional[str]]
            the names. None for each molecule without a name
        """
        p = Path(filename)
        fmt = p.suffix.strip('.')

        smis = []
        names = []
        for mol in pybel.readfile(fmt, filename):
            smis.append(mol.write())

            if mol.title:
                names.append(mol.title)
            else:
                names.append(None)

        return smis, names

    @staticmethod
    def calc_ligand_score(ligand_results: List[List[Dict]],
                          receptor_score_mode: str = 'best',
                          ensemble_score_mode: str = 'best') -> Optional[float]:
        """Calculate the overall score of a ligand given all of its docking
        runs

        Parameters
        ----------
        ligand_results : List[List[Dict]]
            an MxO list of list of dictionaries where each individual 
            dictionary is a record of an individual docking run and
            
            * M is the number of receptors the ligand was docked against
            * O is the number of times each docking run was repeated
        receptor_score_mode : str, default='best'
            the mode used to calculate the overall score for a given receptor
            pose with multiple, repeated runs
        ensemble_score_mode : str, default='best'
            the mode used to calculate the overall score for a given ensemble
            of receptors

        Returns
        -------
        ensemble_score : Optional[float]
            the overall score of a ligand's ensemble docking. None if no such
            score was calculable
        
        See also
        --------
        calc_score
            for documentation on possible values for receptor_score_mode
            and ensemble_score_mode
        """
        receptor_scores = []
        for receptor in ligand_results:
            successful_rep_scores = [
                repeat['score']
                for repeat in receptor if repeat['score'] is not None
            ]
            if successful_rep_scores:
                receptor_scores.append(Screener.calc_score(
                    successful_rep_scores, receptor_score_mode
                ))

        receptor_scores = [score for score in receptor_scores]
        if receptor_scores:
            ensemble_score = Screener.calc_score(
                receptor_scores, ensemble_score_mode)
        else:
            ensemble_score = None
        
        return ensemble_score
    
    @staticmethod
    def calc_score(scores: Sequence[float], score_mode: str = 'best') -> float:
        """Calculate an overall score from a sequence of scores

        Parameters
        ----------
        scores : Sequence[float]
        score_mode : str, default='best'
            the method used to calculate the overall score. Choices include:

            * 'best': return the top score
            * 'avg': return the average of the scores
            * 'boltzmann': return the boltzmann average of the scores

        Returns
        -------
        score : float
        """
        scores = sorted(scores)
        if score_mode in ('best', 'top'):
            score = scores[0]
        elif score_mode in ('avg', 'mean'):
            score = sum(score for score in scores) / len(scores)
        elif score_mode == 'boltzmann':
            Z = sum(exp(-score) for score in scores)
            score = sum(score * exp(-score) / Z for score in scores)
        else:
            score = scores[0]
            
        return score
