"""This module contains functions for preparing input files of vina-type
docking software"""

import csv
from math import ceil, log10
import os
from pathlib import Path
import subprocess as sp
import sys
import timeit
from typing import List, Optional, Sequence, Tuple

from rdkit import Chem
from tqdm import tqdm

OBABEL = sp.run('which obabel', shell=True, encoding='utf-8', 
                stdout=sp.PIPE, check=True).stdout.strip()

def prepare_receptor(receptor: str) -> str:
    """Prepare a receptor PDBQT file from its input file

    Parameter
    ---------
    receptor : str
        the filename of the file containing the receptor to convert

    Returns
    -------
    receptor_pdbqt : str
        the filename of the resulting pdbqt file
    """
    receptor_pdbqt = str(Path(receptor).with_suffix('.pdbqt'))
    args = [OBABEL, receptor, '-O', receptor_pdbqt,
            '-xh', '-xr', '--partialcharge', 'gasteiger']
    sp.run(args, stderr=sp.PIPE, check=True)

    return receptor_pdbqt

def prepare_ligands(ligands, *args, **kwargs) -> List[Tuple[str, str]]:
    if isinstance(ligands, str):
        p_ligand = Path(ligands)

        if not p_ligand.exists():
            return [prepare_from_smi(ligands, *args, **kwargs)]
        if p_ligand.suffix == '.csv':
            return prepare_from_csv(ligands, *args, **kwargs)
        if p_ligand.suffix in {'.sdf', '.smi'}:
            return prepare_from_supply(ligands, *args, **kwargs)
        
        return [prepare_from_file(ligands, *args, **kwargs)]

    elif isinstance(ligands, Sequence):
        return prepare_from_smis(ligands, **kwargs)
    
    raise TypeError('argument "ligand" must be of type str or Sequence[str]!')

def run_obabel(smi: str, pdbqt: str):
    argv = [OBABEL, f'-:{smi}', '-O', pdbqt,
            '-xh', '--gen3d', '--partialcharge', 'gasteiger']
    return sp.run(argv, check=False, stderr=sp.PIPE)

def prepare_from_smi(smi: str, name: str = 'ligand',
                     path: str = '.', **kwargs) -> Tuple[str, str]:
    path = Path(path)
    if not path.is_dir():
        path.mkdir()
    
    pdbqt = path / f'{name}.pdbqt'

    run_obabel(smi, pdbqt)

    return smi, pdbqt

def prepare_from_file(ligand: str, name: Optional[str] = None,
                      path: str = '.', **kwargs) -> Tuple[str, str]:
    """Convert a single ligand to PDBQT format

    Parameters
    ----------
    ligand : str
        either the filename containing the ligand or a SMILES string
        representing the ligand
    name : Optional[str] (Default = None)
        the name of the ligand
        If None, use the stem of the input file
    is_file : bool
        is ligand a filename or a SMILES string?
    path : str (Default = '.')
        the path under which the output .pdbqt file should be written
    name : Optional[string]
        the stem of the file to write this ligand to (i.e., no extension)
    **kwargs
        additional and unused keyword arguments

    Returns
    -------
    smi : str
        the SMILES string corresponding to the converted ligand
    pdbqt : string
        the filename of the output .pdbqt file
    """
    path = Path(path)
    if not path.is_dir():
        path.mkdir()

    name = name or Path(ligand).stem
    pdbqt = str(path / f'{name}.pdbqt')

    res = sp.run([OBABEL, ligand, '-osmi'], stdout=sp.PIPE, check=True)
    smi = res.stdout.decode(sys.getdefaultencoding()).split('\t')[0]

    sp.run([OBABEL, ligand, '-O', pdbqt], stderr=sp.PIPE, check=True)

    return smi, pdbqt

def prepare_from_smis(smis: Sequence[str], 
                      names: Optional[Sequence[str]] = None, 
                      start: int = 0, nconvert: Optional[int] = None,
                      path: str = '.', ncpu: int = 1, distributed: bool = True, 
                      verbose: int = 0, **kwargs) -> List[Tuple[str, str]]:
    """Convert the list of SMILES strings to their corresponding PDBQT files

    Parameters
    ----------
    smis : Sequence[str]
        a sequence of SMILES strings
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
    ligands : List[Tuple[str, str]]
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
        pdbqts = [f'{path}/{name}.pdbqt' for name in names]
        ligands = list(zip(smis, pdbqts))
        results = client.map(run_obabel, smis, pdbqts, chunksize=batch_size)

        ligands = [
            ligand for res, ligand in tqdm(
                zip(results, ligands), total=len(ligands),
                desc='Converting ligands', unit='ligand', smoothing=0.
            ) if res.returncode == 0
        ]
    
    total = timeit.default_timer() - begin
    if verbose > 1:
        print_summary(total, len(ligands))
        
    return ligands

def prepare_from_csv(csv_filename: str, title_line: bool = True,
                     smiles_col: int = 0, name_col: Optional[int] = None,
                     start: int = 0, nconvert: Optional[int] = None,
                     njobs: int = 1, path: str = '.', verbose: int = 0,
                     **kwargs) -> List[Tuple[str, str]]:
    """Prepare the PDBQT files corresponding to the molecules contained in a 
    CSV file

    Parameters
    ----------
    csv_filename : str
        the filename of the CSV file containing the ligands to convert
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
    ligands : List[Tuple[str, str]]
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
    
    ligands = prepare_from_smis(smis, names=names, path=path,
                                start=start, nconvert=nconvert, njobs=njobs)

    total = timeit.default_timer() - begin
    if verbose > 1:
        print_summary(total, len(ligands))

    return ligands

def prepare_from_supply(supply: str, id_prop_name: Optional[str] = None,
                        start: int = 0, nconvert: Optional[int] = None,  
                        njobs: int = 1, path: str = '.', verbose: int = 0,
                        **kwargs) -> List[Tuple[str, str]]:
    """Prepare the PDBQT files corresponding to the molecules contained in a 
    chemical supply file

    Parameters
    ----------
    lig_supply : str
        the filename of the SDF or SMI file containing the ligands to convert
    id_prop_name : Optional[str]
        the name of the property containing the ID, if one exists
            e.g., "CatalogID", "Chemspace_ID", "Name", etc...
    path : str (Default = '.')
        the path under which the directory containing all converted files
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
    ligands : List[Tuple[str, str]]
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
    names = []

    for mol in mols:
        if mol is None:
            continue

        smis.append(Chem.MolToSmiles(mol))
        if id_prop_name:
            names.append(mol.GetProp(id_prop_name))

    ligands = prepare_from_smis(smis, names=names, start=start, 
                                nconvert=nconvert, n_jobs=njobs, path=path)

    total = timeit.default_timer() - begin
    if verbose > 1:
        print_summary(total, len(ligands))

    return ligands

def print_summary(total, n_ligands):
    m, s = divmod(int(total), 60)
    h, m = divmod(m, 60)
    if n_ligands > 0:
        print(f'    Time to prepare {n_ligands} ligands: {h}h {m}m {s}s')
        print(f'    Average time: {total/n_ligands:0.4f} s/ligand', flush=True)
