from pathlib import Path
import subprocess as sp
import sys
from typing import Iterable, Optional, Tuple

from ..utils import OBABEL, Ligand
from ..preparation import prepare_receptors, prepare_ligands
from ...utils import Input

def prepare_inputs(docker: str, receptors: Iterable[str], ligands: Iterable,
                   center: Tuple, size: Tuple[int, int, int] = (20, 20, 20), 
                   ncpu: int = 1, path: str = '.',
                   **kwargs) -> Input:
    # 1 preparation of ligand and receptor
    receptors = prepare_receptors(receptors, prepare_receptor)
    ligands = prepare_ligands(ligands, prepare_from_smi,
                              prepare_from_file, **kwargs)
    # 2 generating receptors surfaces and spheres
    # 3 generating boxes and grids
    # 4 generating ensemble input files
    
    # ligands is type List[Tuple[str, List[str]]] where the innermost
    # List[str] is the ensemble_infiles
    return {'ligands': ligands}

def prepare_receptor(receptor: str):
    """Prepare a receptor mol2 file from its input file

    Parameter
    ---------
    receptor : str
        the filename of a file containing the receptor

    Returns
    -------
    receptor_mol2 : str
        the filename of the resulting MOL2 file
    """
    receptor_mol2 = str(Path(receptor).with_suffix('.mol2'))
    args = [OBABEL, receptor, '-omol2', '-O', receptor_mol2,
            '-h', '--partialcharge', 'gasteiger']
    
    try:
        sp.run(args, stderr=sp.PIPE, check=True)
    except sp.SubprocessError:
        print(f'ERROR: failed to convert {receptor}, skipping...')

    return receptor_mol2

def prepare_from_smi(smi: str, name: str = 'ligand',
                     path: str = '.', **kwargs) -> Optional[Ligand]:
    """Prepare an input ligand file from the ligand's SMILES string

    Parameters
    ----------
    smi : str
        the SMILES string of the ligand
    name : Optional[str] (Default = None)
        the name of the ligand.
    path : str (Default = '.')
        the path under which the output PDBQT file should be written
    **kwargs
        additional and unused keyword arguments

    Returns
    -------
    Optional[Ligand]
        a tuple of the SMILES string and the corresponding prepared input file.
        None if preparation failed for any reason
    """
    path = Path(path)
    if not path.is_dir():
        path.mkdir()
    
    mol2 = str(path / f'{name}.mol2')

    argv = [OBABEL, f'-:{smi}', '-omol2', '-O', mol2,
            '-h', '--gen3d', '--partialcharge', 'gasteiger']
    ret = sp.run(argv, check=False, stderr=sp.PIPE)

    try:
        ret.check_returncode()
        return smi, mol2
    except sp.SubprocessError:
        return None

def prepare_from_file(filename: str, use_3d: bool = False,
                      name: Optional[str] = None, path: str = '.', 
                      **kwargs) -> Ligand:
    """Convert a single ligand to the appropriate input format

    Parameters
    ----------
    filename : str
        the name of the file containing the ligand
    use_3d : bool (Default = False)
        whether to use the 3D information in the input file (if possible)
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
    List[Ligand]
        a tuple of the SMILES string the prepared input file corresponding
        to the molecule contained in filename
    """
    name = name or Path(filename).stem

    ret = sp.run([OBABEL, filename, '-osmi'], stdout=sp.PIPE, check=True)
    lines = ret.stdout.decode('utf-8').splitlines()
    smis = [line.split()[0] for line in lines]

    if not use_3d:
        ligands = [prepare_from_smi(smi, f'{name}_{i}', path) 
                   for i, smi in enumerate(smis)]
        return [lig for lig in ligands if lig]
    
    path = Path(path)
    if not path.is_dir():
        path.mkdir()

    mol2 = f'{path}/{name}_.mol2'
    argv = [OBABEL, filename, '-omol2', '-O', mol2, '-m']
    ret = sp.run(argv, check=False, stderr=sp.PIPE)
    
    try:
        ret.check_returncode()
    except sp.SubprocessError:
        return None

    stderr = ret.stderr.decode('utf-8')
    for line in stderr.splitlines():
        if 'converted' not in line:
            continue
        n_mols = int(line.split()[0])

    mol2s = [f'{path}/{name}_{i}.mol2' for i in range(1, n_mols)]

    return list(zip(smis, mol2s))
