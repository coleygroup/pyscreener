from pathlib import Path
import subprocess as sp
import sys

from pyscreener.docking.utils import OBABEL

def prepare_receptors(receptors, **kwargs):
    """Prepare a receptor mol2 file from its input file

    Parameter
    ---------
    receptors : List[str]
        the filenames of files containing various poses of the receptor

    Returns
    -------
    receptor_pdbqts : List[str]
        the filenames of the resulting PDBQT files
    """
    receptor_mol2s = []
    for receptor in receptors:
        receptor_mol2 = str(Path(receptor).with_suffix('.mol2'))
        args = [OBABEL, receptor, '-omol2', '-O', receptor_mol2, '-h']
        
        try:
            sp.run(args, stderr=sp.PIPE, check=True)
            receptor_mol2s.append(receptor_mol2)
        except sp.SubprocessError:
            print(f'ERROR: failed to convert {receptor}, skipping...')

    if len(receptor_mol2s) == 0:
        raise RuntimeError('Preparation failed for each receptor!')
    
    return receptor_mol2s

def prepare_ligands(ligands, **kwargs):
    pass
