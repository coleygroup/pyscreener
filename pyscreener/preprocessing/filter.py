import csv
from itertools import zip_longest
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from rdkit import Chem
from rdkit.Chem import QED
from tqdm import tqdm

def filter_ligands(ligands: str,
                   **kwargs) -> Tuple[List[str], Optional[List[str]]]:
    if isinstance(ligands, str):
        p_ligand = Path(ligands)

        if p_ligand.suffix == '.csv':
            return filter_csv(ligands, **kwargs)
        if p_ligand.suffix in {'.sdf', '.smi'}:
            return filter_supply(ligands, **kwargs)
        
        return [ligands], None

    elif isinstance(ligands, Sequence):
        return filter_smis(ligands, **kwargs)
    
    raise TypeError('argument "ligand" must be of type str or Sequence[str]!')

def filter_smis(smis: List[str], names: Optional[List[str]] = None,
                **kwargs) -> Tuple[List[str], Optional[List[str]]]:
    mols = [Chem.MolFromSmiles(smi)
            for smi in tqdm(smis, desc='Reading in mols', unit='mol')]

    return filter_mols(mols, names, **kwargs)

def filter_csv(csvfile: str, title_line: bool = True,
               smiles_col: int = 0, name_col: Optional[int] = None,
               **kwargs) -> Tuple[List[str], Optional[List[str]]]:
    with open(csvfile) as fid:
        reader = csv.reader(fid)
        if title_line:
            next(reader)

        reader = tqdm(reader, desc='Reading in mols', unit='mol')
        if name_col is None:
            mols = [Chem.MolFromSmiles(row[smiles_col]) for row in reader]
            names = None
        else:
            mols_names = [(Chem.MolFromSmiles(row[smiles_col]), row[name_col]) 
                           for row in reader]
            mols, names = zip(*mols_names)

    return filter_mols(mols, names, **kwargs)


def filter_supply(supplyfile: str, id_prop_name: Optional[str],
                  **kwargs) -> Tuple[List[str], Optional[List[str]]]:
    p_supply = Path(supplyfile)
    if p_supply.suffix == '.sdf':
        supply = Chem.SDMolSupplier(supplyfile)
    elif p_supply.suffix == '.smi':
        supply = Chem.SmilesMolSupplier(supplyfile)
    else:
        raise ValueError(
            f'input file "{supplyfile}" does not have .sdf or .smi extension')

    supply = tqdm(supply, desc='Reading in mols', unit='mol')
    mols = []
    names = None

    if id_prop_name:
        names = []
        for mol in supply:
            if mol is None:
                continue

            mols.append(Chem.MolToSmiles(mol))
            names.append(mol.GetProp(id_prop_name))
    else:
        for mol in supply:
            if mol is None:
                continue

            mols.append(Chem.MolToSmiles(mol))

    return filter_mols(mols, names, **kwargs)

def filter_mols(mols: List[Chem.Mol], names: Optional[List[str]] = None,
                max_atoms: int = 1000, max_weight: float = 10000.,
                max_logP: float = 10.,
                **kwargs) -> Tuple[List[str], Optional[List[str]]]:
    names = names or []
    
    smis_filtered = []
    names_filtered = []

    for mol, name in tqdm(zip_longest(mols, names), total=len(mols),
                          desc='Filtering mols', unit='mol'):
        if mol.GetNumHeavyAtoms() > max_atoms:
            continue

        props = QED.properties(mol)
        if props.MW > max_weight:
            continue
        if props.ALOGP > max_logP:
            continue

        smis_filtered.append(Chem.MolToSmiles(mol))
        names_filtered.append(name)

    return smis_filtered, names_filtered