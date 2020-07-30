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

from pyscreener.docking.utils import Ligand, OBABEL

def prepare_receptor(receptor: str, **kwargs) -> Optional[str]:
    """Prepare a receptor PDBQT file from its input file

    Parameter
    ---------
    receptor : str
        the filename of a file containing a receptor

    Returns
    -------
    receptor_pdbqt : Optional[str]
        the filenames of the resulting PDBQT files. None if preparation failed
    """
    receptor_pdbqt = str(Path(receptor).with_suffix('.pdbqt'))
    args = [OBABEL, receptor, '-O', receptor_pdbqt,
            '-xh', '-xr', '--partialcharge', 'gasteiger']
    try:
        sp.run(args, stderr=sp.PIPE, check=True)
    except sp.SubprocessError:
        print(f'ERROR: failed to convert {receptor}', file=sys.stderr)
        return None

    return receptor_pdbqt

def prepare_from_smi(smi: str, name: str = 'ligand',
                     path: str = '.', **kwargs) -> Optional[Ligand]:
    path = Path(path)
    if not path.is_dir():
        path.mkdir()
    
    pdbqt = str(path / f'{name}.pdbqt')

    argv = [OBABEL, f'-:{smi}', '-O', pdbqt,
            '-xh', '--gen3d', '--partialcharge', 'gasteiger']
    ret = sp.run(argv, check=False, stderr=sp.PIPE)

    try:
        ret.check_returncode()
        return smi, pdbqt
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
        return [prepare_from_smi(smi, name, path) for smi in smis]
    
    path = Path(path)
    if not path.is_dir():
        path.mkdir()

    pdbqt = f'{path}/{name}_.pdbqt'
    argv = [OBABEL, filename, '-opdbqt', '-O', pdbqt]
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

    pdbqts = [f'{path}/{name}_{i}.pdbqt' for i in range(1, n_mols)]

    return list(zip(smis, pdbqts))