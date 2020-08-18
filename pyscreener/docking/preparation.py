"""This module contains template functions for preparing the inputs necessary
to run virtual screeens in a number of docking programs"""

import csv
from itertools import chain
from math import ceil, log10
import os
from pathlib import Path
import subprocess as sp
import sys
import timeit
from typing import Callable, Iterable, List, Optional, Sequence, Tuple, Union

from rdkit import Chem
from tqdm import tqdm

from .utils import Ligand, OBABEL
from ..utils import Input

def prepare_receptors(receptors: Iterable[str], 
                      prepare_receptor: Callable[[str], str],
                      **kwargs) -> List[str]:   
    """Prepare a receptor input file from its supply file

    Parameter
    ---------
    receptors : Iterable[str]
        the filenames of files containing various poses of the receptor
    prepare_receptor : Callable[[str], str]
        a function that prepares a single receptor file
    **kwargs
        keyword arguments passed to prepare_receptor

    Returns
    -------
    receptor_inputs : List[str]
        the filenames of the resulting input receptor files
    """
    receptor_inputs = []
    for receptor in receptors:
        receptor_input = prepare_receptor(receptor, **kwargs)
        if receptor_input:
            receptor_inputs.append(receptor_input)

    if len(receptor_inputs) == 0:
        raise RuntimeError('Preparation failed for each receptor!')

    return receptor_inputs

def prepare_ligands(ligands: Iterable, *args, **kwargs):
    return list(chain(*(
        prepare_ligand(ligand, *args, **kwargs) for ligand in ligands
    )))

def prepare_ligand(ligand: Union[str, Sequence[str]],
                   prepare_from_smi: Callable[..., Ligand],
                   prepare_from_file: Callable[..., Ligand],
                   use_3d: bool = False, **kwargs) -> List[Ligand]:
    """Prepare the input for a single ligand file

    Parameters
    ----------
    ligand : Union[str, Sequence[str]]
        a single SMILES string, the filepath of a file containing one or many 
        ligands, or a Sequence of SMILES strings each corresponding to one
        ligand
    prepare_from_smi : Callable[..., Ligand]
        a function to prepare a ligand input from a single SMILES string
    prepare_from_file : Callable[..., Ligand]
        a function to prepare a ligand input from a molecule file
    use_3d : bool (Default = False)
        whether to use the 3D coordinates contained in a ligand file
    **kwargs
        keyword arguments for the appropriate function

    Returns
    -------
    List[Ligand]
        the prepared ligands

    Raises
    ------
    TypeError
        if ligand is not a string representing a SMILES string, a string
        representing a filepath, or a Sequence of SMILES strings
    """
    if isinstance(ligand, str):
        p_ligand = Path(ligand)

        if not p_ligand.exists():
            return [prepare_from_smi(ligand, **kwargs)]

        if p_ligand.suffix == '.csv':
            return prepare_from_csv(ligand, prepare_from_smi, **kwargs)
        if p_ligand.suffix == '.smi':
            return prepare_from_supply(ligand, prepare_from_smi, **kwargs)
        if p_ligand.suffix == '.sdf':
            if use_3d:
                return prepare_from_file(ligand, **kwargs)
            else:
                return prepare_from_supply(ligand, prepare_from_smi, **kwargs)
        
        return prepare_from_file(ligand, **kwargs)

    elif isinstance(ligand, Sequence):
        return prepare_from_smis(ligand, prepare_from_smi, **kwargs)
    
    raise TypeError('Arg "ligand" must be of type str or Sequence[str].',
                    f'Got: {type(ligand)}')

def prepare_from_smis(smis: Sequence[str],
                      prepare_from_smi: Callable[..., Ligand],
                      names: Optional[Sequence[str]] = None, 
                      start: int = 0, nconvert: Optional[int] = None,
                      path: str = '.', ncpu: int = 1, distributed: bool = True, 
                      verbose: int = 0, **kwargs) -> List[Ligand]:
    """Convert the list of SMILES strings to their corresponding PDBQT files

    Parameters
    ----------
    smis : Sequence[str]
        a sequence of SMILES strings
    prepare_from_smi: Callable[..., Tuple[str, str]]
        a function that prepares an input ligand file from a SMILES string
    names : Optional[Sequence[str]] (Default = None)
        a parallel sequence of names for each ligand
    path : str (Default = '.')
        the path under all converted files should be organized. By default, 
        uses the current directory
    start : int (Default = 0)
        the index at which to start ligand preparation
    nconvert : Optional[int] (Default = None)
        the number of ligands to convert. If None, convert all ligands
    ncpu : int (Default = 1)
        the number of cores available to each worker
    verbose : int, (Default = 0)
        whether or not to print performance data
    **kwargs
        additional and unused keyword arguments

    Returns
    -------
    ligands : List[Ligand]
        a list of tuples containing a ligand's SMILES string and the filepath
        of the corresponding .pdbqt file, which are named:
            lig0.pdbqt, lig1.pdbqt, ...
    """
    begin = timeit.default_timer()
    
    if distributed:
        from mpi4py import MPI
        from mpi4py.futures import MPIPoolExecutor as Pool

        n_workers = MPI.COMM_WORLD.size * ncpu
    else:
        from concurrent.futures import ProcessPoolExecutor as Pool
        try:
            n_workers = len(os.sched_getaffinity(0))
        except AttributeError:
            n_workers = os.cpu_count()

    path = Path(path)
    if not path.is_dir():
        path.mkdir(parents=True)
    
    stop = min(len(smis), start+nconvert) if nconvert else len(smis)

    if names is None:
        width = ceil(log10(len(smis))) + 1
        names = (f'ligand_{i:0{width}}' for i in range(start, stop))
    else:
        # could theoretically handle empty strings, but I think that's the
        # caller's responsiblity
        names = names[start:stop]

    smis = smis[start:stop]

    CHUNKS_PER_WORKER = 16
    batch_size = ceil(len(smis) / (CHUNKS_PER_WORKER*n_workers))

    with Pool(max_workers=n_workers) as client:
        paths = (path for _ in range(len(smis)))
        ligands = [
            ligand for ligand in tqdm(
                client.map(prepare_from_smi, smis, names, paths, 
                           chunksize=batch_size), total=len(smis),
                desc='Preparing ligands', unit='ligand', smoothing=0.
            ) if ligand
        ]
    
    total = timeit.default_timer() - begin
    if verbose > 1:
        print_summary(total, len(ligands))
        
    return ligands

def prepare_from_csv(csv_filename: str,
                     prepare_from_smi: Callable[..., Ligand],
                     title_line: bool = True,
                     smiles_col: int = 0, name_col: Optional[int] = None,
                     start: int = 0, nconvert: Optional[int] = None,
                     njobs: int = 1, path: str = '.', verbose: int = 0,
                     **kwargs) -> List[Ligand]:
    """Prepare the input files corresponding to the molecules contained in a 
    CSV file

    Parameters
    ----------
    csv_filename : str
        the filename of the CSV file containing the ligands to convert
    prepare_from_smi: Callable[..., Tuple[str, str]]
        a function that prepares an input ligand file from a SMILES string
    title_line : bool (Default = True)
        does the CSV file contain a title line?
    smiles_col : int (Default = 0)
        the column containing the SMILES strings
    name_col : Optional[int] (Default = None)
        the column containing the molecule name
    path : str (Default = '.')
        the path under which the directory containing all converted files
        should be organized. By default, uses the current directory
    dir_name : Optional[str] (Default = None)
        the name of the directory to which all resulting files should be
        written. If None, default to using '<lig_supply>_pdbqt/'
    start : int (Default = 0)
        the index at which to start conversion
    nconvert : Optional[int] (Default = None)
        the number of ligands to convert. If None, convert all molecules
    njobs : int (Default = 1)
        the number of jobs to distribute input preparation over
    verbose : int (Default = 0)
        the level of output to print
    **kwargs
        additional and unused keyword arguments

    Returns
    -------
    ligands : List[Ligand]
        a list of tuples containing a ligand's SMILES string and the filepath
        of the corresponding .pdbqt file.
        PDBQT files are named <compound_id>.pdbqt if compound_id property
        exists in the original supply file. Otherwise, they are named:
            lig0.pdbqt, lig1.pdbqt, ...
    """
    begin = timeit.default_timer()

    with open(csv_filename) as fid:
        reader = csv.reader(fid)
        if title_line:
            next(reader)

        if name_col is None:
            smis = [row[smiles_col] for row in reader]
            names = None
        else:
            smis_names = [(row[smiles_col], row[name_col]) for row in reader]
            smis, names = zip(*smis_names)
    
    ligands = prepare_from_smis(smis, prepare_from_smi,
                                names=names, path=path,
                                start=start, nconvert=nconvert, njobs=njobs)

    total = timeit.default_timer() - begin
    if verbose > 1:
        print_summary(total, len(ligands))

    return ligands

def prepare_from_supply(supply: str,
                        prepare_from_smi: Callable[..., Ligand],
                        id_prop_name: Optional[str] = None,
                        start: int = 0, nconvert: Optional[int] = None,  
                        njobs: int = 1, path: str = '.', verbose: int = 0,
                        **kwargs) -> List[Ligand]:
    """Prepare the input files corresponding to the molecules contained in a 
    molecule supply file

    Parameters
    ----------
    supply : str
        the filename of the SDF or SMI file containing the ligands to convert
    id_prop_name : Optional[str]
        the name of the property containing the ID, if one exists
            e.g., "CatalogID", "Chemspace_ID", "Name", etc...
    path : str (Default = '.')
        the path under which the directory containing all prepared files
        should be organized. By default, uses the current directory
    dir_name : Optional[str] (Default = None)
        the name of the directory to which all resulting files should be
        written. If None, default to using '<lig_supply>_pdbqt/'
    start : int (Default = 0)
        the index at which to start ligand conversion
    nconvert : Optional[int] (Default = None)
        the number of ligands to convert. If None, convert all molecules
    njobs : int (Default = 1)
        the number of jobs to distribute input preparation over
    verbose : int (Default = 0)
        the level of output to print
    **kwargs
        additional and unused keyword arguments

    Returns
    -------
    ligands : List[Ligand]
        a list of tuples containing a ligand's SMILES string and the filepath
        of the corresponding .pdbqt file.
        PDBQT files are named <compound_id>.pdbqt if compound_id property
        exists in the original supply file. Otherwise, they are named:
            lig0.pdbqt, lig1.pdbqt, ...
    """
    begin = timeit.default_timer()

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

    ligands = prepare_from_smis(smis, prepare_from_smi, 
                                names=names, path=path,
                                start=start, nconvert=nconvert, n_jobs=njobs)

    total = timeit.default_timer() - begin
    if verbose > 1:
        print_summary(total, len(ligands))

    return ligands

def print_summary(total, n_ligands):
    m, s = divmod(int(total), 60)
    h, m = divmod(m, 60)
    if n_ligands > 0:
        print(f'    Time to prepare {n_ligands} ligands: {h}h {m}m {s}s',
              f'({total/n_ligands:0.4f} s/ligand)', flush=True)
