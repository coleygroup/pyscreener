import csv
from math import ceil, log10
import os
from pathlib import Path
import subprocess as sp
import sys
import timeit
from typing import Callable, List, Optional, Sequence, Tuple

from rdkit import Chem
from tqdm import tqdm

from . import autodock, ucsfdock, rdock

OBABEL = sp.run('which obabel', shell=True, encoding='utf-8', 
                stdout=sp.PIPE, check=True).stdout.strip()

def prepare_receptor(docker: str, receptors, **kwargs):
    if docker in {'vina', 'smina', 'psovina', 'qvina'}:
        return autodock.prepare_receptors(receptors)
    if docker == 'dock':
        return ucsfdock.prepare_receptors(receptors, **kwargs)
    if docker == 'rdock':
        return rdock.prepare_receptors(receptors, **kwargs)

    raise ValueError(f'{docker} is not a recognized docking program')

# def prepare_ligands(docker: str, ligands, **kwargs):
#     if docker in {'vina', 'smina', 'psovina', 'qvina'}:
#         return autodock.prepare_ligands(ligands, **kwargs)
#     if docker == 'dock':
#         return ucsfdock.prepare_ligands(ligands, **kwargs)
#     if docker == 'rdock':
#         return rdock.prepare_ligands(ligands, **kwargs)

#     raise ValueError(f'{docker} is not a recognized docking program')

def prepare_ligands(docker, ligands, **kwargs) -> List[Tuple[str, str]]:
    prepare_from_smi = {
        'vina': autodock.prepare_from_smi,
        'smina': autodock.prepare_from_smi,
        'qvina': autodock.prepare_from_smi,
        'psovina': autodock.prepare_from_smi,
        'ucsfdock': ucsfdock.prepare_from_smi
    }.get(docker)
    
    if prepare_from_smi is None:
        raise ValueError(f'Unrecognized docking program: {docker}')
    
    if isinstance(ligands, str):
        p_ligand = Path(ligands)

        if not p_ligand.exists():
            return [prepare_from_smi(ligands, **kwargs)]
        if p_ligand.suffix == '.csv':
            return prepare_from_csv(ligands, prepare_from_smi, **kwargs)
        if p_ligand.suffix in {'.sdf', '.smi'}:
            return prepare_from_supply(ligands, prepare_from_smi, **kwargs)
        
        return [prepare_from_file(ligands, prepare_from_smi, **kwargs)]

    elif isinstance(ligands, Sequence):
        return prepare_from_smis(ligands, prepare_from_smi, **kwargs)
    
    raise TypeError('argument "ligand" must be of type str or Sequence[str]!')

def prepare_from_file(filename: str, prepare_from_smi: Callable,
                      name: Optional[str] = None, path: str = '.', 
                      **kwargs) -> Tuple[str, str]:
    """Convert a single ligand to the appropriate input format

    Parameters
    ----------
    filename : str
        the name of the file containing the ligand
    prepare_from_smi: Callable[..., Tuple[str, str]]
        a function that prepares an input ligand file from a SMILES string
    name : Optional[str] (Default = None)
        the name of the ligand. If None, use the stem of the input file
    path : str (Default = '.')
        the path under which the output .pdbqt file should be written
    **kwargs
        additional and unused keyword arguments

    Returns
    -------
    smi : str
        the SMILES string corresponding to the converted ligand
    input_filename : string
        the filename of the prepared input file
    """
    name = name or Path(filename).stem

    res = sp.run([OBABEL, filename, '-osmi'], stdout=sp.PIPE, check=True)
    smi = res.stdout.decode(sys.getdefaultencoding()).split('\t')[0]

    return prepare_from_smi(smi, name, path)

def prepare_from_smis(smis: Sequence[str], prepare_from_smi: Callable,
                      names: Optional[Sequence[str]] = None, 
                      start: int = 0, nconvert: Optional[int] = None,
                      path: str = '.', ncpu: int = 1, distributed: bool = True, 
                      verbose: int = 0, **kwargs) -> List[Tuple[str, str]]:
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
        paths = (path for _ in range(len(smis)))
        # pdbqts = [f'{path}/{name}.pdbqt' for name in names]
        # ligands = list(zip(smis, pdbqts))
        ligands = [
            ligand for ligand in client.map(
                prepare_from_smi, smis, names, paths, chunksize=batch_size
            ) if ligand
        ]
        # ligands = [ligand for ligand in ligands if ligand]
        # ligands = [
        #     ligand for res, ligand in tqdm(
        #         zip(results, ligands), total=len(ligands),
        #         desc='Converting ligands', unit='ligand', smoothing=0.
        #     ) if res.returncode == 0
        # ]
    
    total = timeit.default_timer() - begin
    if verbose > 1:
        print_summary(total, len(ligands))
        
    return ligands

def prepare_from_csv(csv_filename: str, prepare_from_smi: Callable,
                     title_line: bool = True,
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
    
    ligands = prepare_from_smis(smis, prepare_from_smi,
                                names=names, path=path,
                                start=start, nconvert=nconvert, njobs=njobs)

    total = timeit.default_timer() - begin
    if verbose > 1:
        print_summary(total, len(ligands))

    return ligands

def prepare_from_supply(supply: str, prepare_from_smi: Callable,
                        id_prop_name: Optional[str] = None,
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
