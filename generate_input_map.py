from collections import defaultdict
import csv
from pathlib import Path
from typing import Dict, Optional, Tuple
import sys

from rdkit import Chem

def parse_pdbqt(filepath: Path, smiles_col: Optional[int] = None
                validate: bool = False) -> Optional[str]
    """Parse a pdbqt file to extract the SMILES string

    Parameters
    ----------
    filepath : Path
        the filepath of a pdbqt file from which to extract the SMILES string
    smiles_col : Optional[int] (Default = None)
        the column containing the SMILES string. If None, return the first
        column containing a valid SMILES string
    validate : bool (Default = False)
        whether the SMILES string should be validated

    Returns
    -------
    smi : Optional[str]
        the SMILES string of the ligand contained in the file. None if
        no SMILES string was found or it was invalid
    """
    with open(filepath) as fid:
        for line in fid:
            if 'SMILES' not in line:
                continue

            tokens = line.split()

            if smiles_col:
                smi = line[smiles_col]
                if validate:
                    if Chem.MolFromSmiles(smi):
                        return smi
                    else:
                        return None
                return smi

            # if smiles_col is unspecified, scan through line for first
            # valid SMILES string and return it
            for token in tokens:
                smi = Chem.MolFromSmiles(token)
                if smi:
                    return smi

    return None

def collect_pdbqts(dir: Path, smiles_col: Optional[int] = None,
                   validate: bool = False) -> Dict[str, List[str]]:
    """Recursively crawl a directory collecting

    Parameters
    ----------
    dir : Path
        the directory to crawl
    smiles_col : Optional[int] (Default = None)
        the column in the PDBQT files containing the SMILES string
    validate : bool (Default = False)
        whether the SMILES strings should be validated

    Returns
    -------
    d_smi_pdbqts : Dict[str, List[str]]
        a dictionary mapping from SMILES string to all input files
        that correspond to it.
    """
    d_smi_pdbqts = defaultdict(list)
    for child in dir.iterdir():
        if child.is_dir():
            d_smi_pdbqts.update(collect_pdbqts(child))
        elif child.suffix == '.pdbqt':
            smi = parse_pdbqt(child, smiles_col, validate)
            if smi:
                d_smi_pdbqts[smi].append(str(child))

    return d_smi_pdbqts

def main():
    if len(sys.argv) < 2:
        raise ValueError('no root directory specified')
    elif len(sys.argv) == 2:
        d_smi_pdbqts = collect_pdbqts(sys.argv[1])
    elif len(sys.argv) == 3:
        root_dir = sys.argv[1]
        smiles_col = sys.argv[2]
        d_smi_pdbqts = collect_pdbqts(root_dir, smiles_col)
    else:
        root_dir = sys.argv[1]
        smiles_col = sys.argv[2]
        validate = bool(sys.argv[3])
        d_smi_pdbqts = collect_pdbqts(root_dir, smiles_col, validate)

    input_map_filename = 'smi_pdbqts.csv'
    with open(input_map_filename, 'w') as fid:
        writer = csv.writer(fid)
        writer.writerow(['smiles', 'pdbqts'])
        for smi, pdbqts in d_smi_pdbqts.items():
            writer.writerow([smi, *pdbqts])

    print(f'Inputs map has been written to {input_map_filename}')
if __name__ == '__main__':
    main()
