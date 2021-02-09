from abc import ABC, abstractmethod
from concurrent.futures import Executor, ProcessPoolExecutor as Pool
import csv
import datetime
from itertools import chain
from math import ceil, exp, log10
import os
from pathlib import Path
import re
import shutil
import tarfile
import tempfile
import timeit
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Type

from openbabel import pybel
import ray
from rdkit import Chem
from tqdm import tqdm

from pyscreener.preprocessing import pdbfix

if not ray.is_initialized():
    try:
        ray.init(address='auto')
    except ConnectionError:
        ray.init()

class Screener(ABC):
    """A Screener conducts virtual screens against an ensemble of receptors.

    Classes that implement the Screener interface are responsible for
    defining the following methods:

    * prepare_receptor
    * prepare_from_smi
    * prepare_from_file
    * run_docking
    * parse_ligand_results

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

    Parameters
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
                 distributed: bool = False,
                 num_workers: int = -1, ncpu: int = 1,
                 path: str = '.', tmp_dir = tempfile.gettempdir(),
                 verbose: int = 0, **kwargs):
        self.path = path
        self.tmp_dir = tmp_dir

        receptors = receptors or []
        if pdbids:
            receptors.extend((
                pdbfix.pdbfix(pdbid=pdbid, path=self.in_path)
                for pdbid in pdbids
            ))
        if len(receptors) == 0:
            raise ValueError('No receptors or PDBids provided!')

        self.receptors = receptors
        self.repeats = repeats
        self.score_mode = score_mode
        self.receptor_score_mode = receptor_score_mode
        self.ensemble_score_mode = ensemble_score_mode
        
        self.distributed = distributed
        self.num_workers = num_workers
        self.ncpu = ncpu

        self.verbose = verbose

        self.num_docked_ligands = 0
        
    def __len__(self) -> int:
        """The number of ligands this screener has simulated"""
        return self.num_docked_ligands

    def __call__(self, *args, **kwargs) -> Dict[str, Optional[float]]:
        return self.dock(*args, **kwargs)
    
    @property
    def path(self) -> os.PathLike:
        """the Screener's parent directory"""
        return self.__path
        
    @path.setter
    def path(self, path: str):
        """set both input and output directories"""
        path = Path(path)
        if not path.is_dir():
            path.mkdir(parents=True)
        self.__path = path
        self.in_path = f'{path}/inputs'
        self.out_path = f'{path}/outputs'

    @property
    def in_path(self) -> os.PathLike:
        return self.__in_path

    @in_path.setter
    def in_path(self, path: str):
        path = Path(path)
        if not path.is_dir():
            path.mkdir(parents=True)
        self.__in_path = path

    @property
    def out_path(self) -> os.PathLike:
        return self.__out_path

    @out_path.setter
    def out_path(self, path: str):
        path = Path(path)
        if not path.is_dir():
            path.mkdir(parents=True)
        self.__out_path = path

    @property
    def tmp_dir(self) -> os.PathLike:
        """the Screener's temp directory"""
        return self.__tmp_dir
        
    @tmp_dir.setter
    def tmp_dir(self, path: str):
        """set both the temp input and output directories"""
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        path = Path(path) / 'pyscreener' / f'session_{timestamp}'
        if not path.is_dir():
            path.mkdir(parents=True)
        self.__tmp_dir = path
        self.tmp_in = path / 'inputs'
        self.tmp_out = path / 'outputs'

    @property
    def tmp_in(self) -> os.PathLike:
        return self.__tmp_in

    @tmp_in.setter
    def tmp_in(self, path: str):
        path = Path(path)
        if not path.is_dir():
            path.mkdir(parents=True)
        self.__tmp_in = path
    
    @property
    def tmp_out(self) -> os.PathLike:
        return self.__tmp_out

    @tmp_out.setter
    def tmp_out(self, path: str):
        path = Path(path)
        if not path.is_dir():
            path.mkdir(parents=True)
        self.__tmp_out = path

    def dock(self, *smis_or_files: Iterable,
             full_results: bool = False,
             **kwargs) -> Dict[str, Optional[float]]:
        """dock the ligands contained in sources

        NOTE: the star operator, *, in the function signature.
            If intending to pass multiple filepaths as an iterable, first 
            unpack the iterable in the function call by prepending a *.
            If passing multiple SMILES strings, either option is acceptable,
            but it is much more efficient to NOT unpack the iterable.

        Parameters
        ----------
        smis_or_files: Iterable
            an iterable of ligand sources, where each ligand source may be
            one of the following:
            - a ligand supply file,
            - a list of SMILES strings
            - a single SMILES string
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
            - smiles: the ligand's SMILES string
            - name: the name of the ligand
            - in: the filename of the input ligand file
            - out: the filename of the output docked ligand file
            - log: the filename of the output log file
            - score: the ligand's docking score
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

            * a ligand supply file
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

            * N is the number of total ligands that will be docked
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
        begin = timeit.default_timer()

        # ligands = self.prepare_ligands(*smis_or_files, **kwargs)
        # recordsss = self.run_docking(ligands)

        smis = []
        names = []
        for smi_or_file in smis_or_files:
            smis_, names_ = self.get_smis(
                smi_or_file, len(smis) + len(self), **kwargs
            )
            smis.extend(smis_), names.extend(names_)

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

    def collect_files(self, out_path: str):
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
        """
        out_path = Path(out_path)
        if not out_path.is_dir():
            out_path.mkdir(parents=True)

        refs = []
        for node in ray.nodes():    # run on all nodes
            address = node["NodeManagerAddress"]

            @ray.remote(resources={f'node:{address}': 0.1})
            def copy_tmp_dir():
                output_id = re.sub('[:,.]', '', ray.state.current_node_id())
                tmp_tar = (self.tmp_dir / output_id).with_suffix('.tar.gz')

                with tarfile.open(tmp_tar, 'w:gz') as tar:
                    tar.add(self.tmp_dir, arcname=output_id)

                shutil.move(str(tmp_tar), str(out_path))

            refs.append(copy_tmp_dir.remote())
        ray.wait(refs)

    @abstractmethod
    def prepare_and_dock(
        self, smis: Sequence[str], names: Sequence[str]
    ) -> List[List[List[Dict]]]:
        pass

    @abstractmethod
    def run_docking(self, ligands: Sequence[Tuple[str, str]]
                   ) -> List[List[List[Dict]]]:
        """Run the docking simulations for the input ligands
        
        Parameters
        ----------
        ligands : Sequence[Tuple[str, str]]
            a sequence of tuples containing a ligand's SMILES string and the 
            filepath of the corresponding input file
        
        Returns
        -------
        List[List[List[Dict]]]
            an NxMxO list of dictionaries where each individual dictionary is a 
            record of an individual docking run and
            
            * N is the number of ligands contained in the ligand sources
            * M is the number of receptors in the ensemble against which each \
                ligand should be docked
            * O is the number of times each docking run should be repeated
            
            NOTE: the records contain a 'score' that is None for each entry
            as the log/out files must first be parsed to obtain the value
        """

    @staticmethod
    @abstractmethod
    def parse_ligand_results(recs_reps: List[List[Dict]],
                             score_mode: str = 'best') -> List[List[Dict]]:
        """Parse the results of the docking simulations for a single ligand
        
        Parameters
        ----------
        recs_reps : List[List[Dict]]
            an MxO list of list of dictionaries where each individual 
            dictionary is a record of an individual docking run and
            
            * M is the number of receptors in the ensemble against which each ligand should be docked
            * O is the number of times each docking run should be repeated
        

        Returns
        -------
        recs_reps : List[List[Dict]]
            the same List as the input argument, but with the 
            'score' key of record updated to reflect the desired 
            score parsed from each docking run
        """

    @property
    def receptors(self):
        return self.__receptors

    @receptors.setter
    def receptors(self, receptors):
        receptors = [self.prepare_receptor(receptor) for receptor in receptors]
        receptors = [receptor for receptor in receptors if receptor is not None]
        if len(receptors) == 0:
            raise RuntimeError('Preparation failed for all receptors!')
        self.__receptors = receptors
    
    @abstractmethod
    def prepare_receptor(self, *args, **kwargs) -> Optional[str]:
        """Prepare a receptor input file for the docking software and return
        the corresponding filepath"""

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

    def prepare_ligands(self, *sources, path: Optional[str] = None,
                        **kwargs) -> List[Tuple[str, str]]:
        path = path or self.in_path
        return list(chain(*(
            self._prepare_ligands(source, i+len(self), path, **kwargs)
            for i, source in enumerate(sources)
        )))

    def _prepare_ligands(self, source, i: int, path: Optional[str] = None,
                         **kwargs) -> List[Tuple[str, str]]:
        if isinstance(source, str):
            p_source = Path(source)

            if not p_source.exists():
                return [self.prepare_from_smi(source, f'ligand_{i}', path)]

            if p_source.suffix == '.csv':
                return self.prepare_from_csv(source, **kwargs)
            if p_source.suffix == '.smi':
                return self.prepare_from_supply(source, **kwargs)
            if p_source.suffix == '.sdf':
                if kwargs['use_3d']:
                    return self.prepare_from_file(source, path=path,
                                                  **kwargs)
                else:
                    return self.prepare_from_supply(source, **kwargs)
            
            return self.prepare_from_file(source, path=path, **kwargs)

        if isinstance(source, Sequence):
            return self.prepare_from_smis(source, **kwargs)
        
        raise TypeError('Arg "source" must be of type str or ', 
                        f'Sequence[str]. Got: {type(source)}')

    def prepare_from_smis(self, smis: Sequence[str],
                          names: Optional[Sequence[str]] = None, 
                          start: int = 0, nconvert: Optional[int] = None,
                          **kwargs) -> List[Tuple]:
        """Convert the list of SMILES strings to their corresponding input files

        Parameters
        ----------
        smis : Sequence[str]
            a sequence of SMILES strings
        names : Optional[Sequence[str]], default=None
            a parallel sequence of names for each ligand
        start : int, default=0
            the index at which to start ligand preparation
        nconvert : Optional[int], default=None
            the number of ligands to convert. If None, convert all ligands
        **kwargs
            additional and unused keyword arguments

        Returns
        -------
        ligands : List[Tuple]
            a list of tuples containing a ligand's SMILES string and the 
            filepath of the corresponding input file
        """
        begin = timeit.default_timer()
        
        stop = min(len(smis), start+nconvert) if nconvert else len(smis)

        if names is None:
            width = ceil(log10(len(smis))) + 1
            names = (f'ligand_{i:0{width}}' for i in range(start, stop))
        else:
            # could theoretically handle empty strings
            names = names[start:stop]
        smis = smis[start:stop]
        paths = (self.in_path for _ in range(len(smis)))

        @ray.remote
        def prepare_from_smi_(smi, name, path):
            return self.prepare_from_smi(smi, name, path)
        
        refs = list(map(prepare_from_smi_.remote, smis, names, paths))
        ligands = [ray.get(r) for r in tqdm(refs, desc='Preparing ligands')]
        
        total = timeit.default_timer() - begin
        if self.verbose > 1 and len(ligands) > 0:
            m, s = divmod(int(total), 60)
            h, m = divmod(m, 60)
            print(f'    Time to prepare {len(ligands)} ligands: ',
                    f'{h}h {m}m {s}s ({total/len(ligands):0.4f} s/ligand)', 
                    flush=True)
            
        return ligands

    # @staticmethod
    # def pmap(f, *args, ncpu=1, chunksize=4):
    #     with Pool(max_workers=ncpu) as client:
    #         xs = [x for x in client.map(f, *args, chunksize=chunksize) if x]
    #     return xs

    def get_smis(self, source: str, offset: int = 0,
                 **kwargs) -> Tuple[List[str], List[str]]:
        """
        get the SMILES strings and names from an input source file

        Parameters
        ----------
        source : str
            the file to parse
        offset : int, default=0
            the offset for numbering ligands
        **kwargs
            keyword arguments to pass to the appropriate
            get_smis_from_* method

        Returns
        -------
        smis : List[str]
            the SMILES strings
        names : List[str]
            the names of the molecules

        Raises
        ------
        ValueError
            if the file does not exit
        """
        source = Path(source)
        if not source.exists():
            return [str(source)], [f'ligand_{offset}']

        if source.suffix == '.csv':
            smis, names = self.get_smis_from_csv(source, **kwargs)
        elif source.suffix in ('.smi', '.sdf'):
            smis, names = self.get_smis_from_supply(source, **kwargs)
        else:
            smis, names = self.get_smis_from_file(source, offset, **kwargs)

        if names is None:
            names = [f'ligand_{i+offset}' for i in range(len(smis))]
        
        return smis, names

    @staticmethod
    def get_smis_from_csv(csv_filename: str, title_line: bool = True,
                          smiles_col: int = 0, name_col: Optional[int] = None,
                          **kwargs) -> Tuple[List[str], Optional[List[str]]]:
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
        names : Optional[List[str]]
            the names parsed from the CSV. None if there were none.
        """
        with open(csv_filename) as fid:
            reader = csv.reader(fid)
            if title_line:
                next(reader)

            if name_col is None:
                smis = [row[smiles_col] for row in reader]
                names = None
            else:
                smis, names = zip(*[
                    (row[smiles_col], row[name_col]) for row in reader
                ])
        
        return smis, names

    @staticmethod
    def get_smis_from_supply(supply: str,
                             id_prop_name: Optional[str] = None,
                             **kwargs) -> List[Tuple]:
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
        names : Optional[List[str]]
            the names parsed from the supply file. None if there were none.
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
        names = None

        if id_prop_name:
            names = []
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

        return smis, names

    @staticmethod
    def get_smis_from_file(filename: str, offset: int = 0,
                           **kwargs) -> Tuple[List[str], List[str]]:
        """get the SMILES strings and names from an arbitrary chemical 
        file format

        Parameters
        ----------
        filename : str
            the file to parse
        offset : int, default=0
        **kwargs
            additional and unused keyword arguments

        Returns
        -------
        smis : List[str]
            the SMILES strings
        names : List[str]
            the names
        """
        p = Path(filename)
        fmt = p.suffix.strip('.')

        smis = []
        names = []
        for i, mol in enumerate(pybel.readfile(fmt, filename)):
            smis.append(mol.write())

            if mol.title:
                names.append(mol.title)
            else:
                names.append(f'ligand_{i+offset}')

        return smis, names

    def prepare_from_csv(self, csv_filename: str, title_line: bool = True,
                         smiles_col: int = 0, name_col: Optional[int] = None,
                         start: int = 0, nconvert: Optional[int] = None,
                         **kwargs) -> List[Tuple]:
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
        start : int, default=0
            the index at which to start conversion
        nconvert : Optional[int], default=None
            the number of ligands to convert. If None, convert all molecules
        **kwargs
            additional and unused keyword arguments

        Returns
        -------
        ligands : List[Tuple]
            a list of tuples containing a ligand's SMILES string and the 
            filepath of the corresponding input file. Files are named 
            <compound_id>.<suffix> if compound_id property exists in the 
            original supply file. Otherwise, they are named:
                
                lig0.<suffix>, lig1.<suffix>, ...
        """
        with open(csv_filename) as fid:
            reader = csv.reader(fid)
            if title_line:
                next(reader)

            if name_col is None:
                smis = [row[smiles_col] for row in reader]
                names = None
            else:
                smis, names = zip(*[
                    (row[smiles_col], row[name_col]) for row in reader
                ])
                # smis, names = zip(*smis_names)
        
        return self.prepare_from_smis(smis, names=names,
                                      start=start, nconvert=nconvert)

    def prepare_from_supply(self, supply: str,
                            id_prop_name: Optional[str] = None,
                            start: int = 0, nconvert: Optional[int] = None,  
                            **kwargs) -> List[Tuple]:
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
        start : int, default=0
            the index at which to start ligand conversion
        nconvert : Optional[int], default=None
            the number of ligands to convert. If None, convert all molecules
        **kwargs
            additional and unused keyword arguments

        Returns
        -------
        ligands : List[Tuple[str, str]]
            a list of tuples containing a ligand's SMILES string and the 
            filepath of the corresponding input file. Files are named 
            <compound_id>.<suffix> if compound_id property exists in the 
            original supply file. Otherwise, they are named:
                
                lig0.<suffix>, lig1.<suffix>, ...
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
        names = None

        if id_prop_name:
            names = []
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

        return self.prepare_from_smis(smis, names=names,
                                      start=start, nconvert=nconvert)

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
            the method used to calculate the overall score

            Choices:

            * 'best' - return the top score
            * 'avg' - return the average of the scores
            * 'boltzmann' - return the boltzmann average of the scores

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
    
    @staticmethod
    def Pool(distributed: bool = False, num_workers: int = -1, ncpu: int = 1,
             all_cores: bool = False) -> Type[Executor]:
        """build a process pool to parallelize computation over

        Parameters
        ----------
        distributed : bool, default=False
            whether to return a distributed or a local process pool
        num_workers : int, default=-1
            if distributed is True, then this argument is ignored. If False,
            then it should be equal to the total number of worker processes
            desired. Using a value of -1 will spawn as many worker processes
            as cores available on this machine.

            NOTE: this is usually not a good idea and it's much better to
            specify the number of processes explicitly.
        ncpu : int, default=1
            if distributed is True, then this argument should be the number of 
            cores allocated to each worker. if False, then this should be the
            number of cores that is desired to be allocated to each worker.

            NOTE: this is an implicit argument because Screener.dock() will   
            make subprocess calls to progams that themselves can utilize 
            multiple cores. It will not actually assign <ncpu> cores to 
            each worker process.
        all_cores : bool (Default = False)
            whether to initialize as many processes as cores available
            (= num_workers * ncpu).
        
        Returns
        -------
        Executor
            the initialized process pool
        
        Notes
        ------
        in some cases, as shown in the examples below, the values specified for
        num_workers and ncpu will be inconsequential. Regardless, it is good
        practice for this function to always be called the same way, with only
        all_cores changing, depending on the context in which the initialized Executor will be used

        **Ex. 1**

        *Given:* a single machine with 16 cores, screening using vina-type
        docking software (via the docking.Vina class)
        
        the function should be called with distributed=False, all_cores=False, 
        and both num_workers and ncpu should be specified such that the product 
        of the two is equal to 16.
        Choices: (1, 16), (2, 8), (4, 4), (8, 2), and (16, 1). You will often have to determine the optimal values empirically.

        **Ex. 2**

        *Given:* a cluster of machines where you've requested resources for 8
        tasks with 2 cores each. The software was then initialized with
        8 separate MPI processes and screening using vina-type docking
        software is to be performed.

        the function should be called with distributed=True and all_cores=False
        (neither num_workers or ncpu needs to be specified)

        **Ex. 3**

        *Given:* a single machine with 16 cores, and pure python code is to be
        executed in parallel

        the function should be called with distributed=False, all_cores=True,
        and both num_workers and ncpu should be specified such that the product 
        of the two is equal to 16.
        Choices: see Ex. 1
        """
        if distributed:
            from mpi4py import MPI
            from mpi4py.futures import MPIPoolExecutor as Pool

            num_workers = MPI.COMM_WORLD.size
        else:
            if num_workers == -1:
                try:
                    num_workers = len(os.sched_getaffinity(0))
                except AttributeError:
                    num_workers = os.cpu_count()
                ncpu = 1

        if all_cores and not distributed:
            num_workers *= ncpu

        return Pool(max_workers=num_workers)
