"""This module contains functions for generating the feature matrix of a set of
molecules located either in a sequence of SMILES strings or in a file"""
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Set, Tuple

import h5py
import numpy as np
import ray
from rdkit.Chem import AllChem as Chem
from tqdm import tqdm

from pyscreener.utils import chunks


@ray.remote
def smis_to_fps(
    smis: Iterable[str], radius: int = 2, length: int = 2048
) -> List[Optional[np.ndarray]]:
    fps = []
    for smi in smis:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            fps.append(None)
        else:
            fps.append(
                np.array(
                    Chem.GetMorganFingerprintAsBitVect(mol, radius, nBits=length, useChirality=True)
                )
            )

    return fps


def gen_fps_h5(
    smis: Sequence[str], path: str = ".", name: str = "fps", radius: int = 2, length: int = 2048
) -> Tuple[str, Set[int]]:
    """Generate an hdf5 file containing the feature matrix of the list of
    SMILES strings

    Parameters
    ----------
    smis : Sequence[str]
        the SMILES strings from which to generate fingerprints
    path : str, default='.'
        the path under which the hdf5 file should be written
    name : str, default ='fps'
        the name of the output hdf5 file
    radius : int, default=2
        the radius of the fingerprints
    length : int, default=2048
        the length of the fingerprints
    **kwargs
        additional and unused keyword arguments

    Returns
    -------
    fps_h5 : str
        the filename of the hdf5 file containing the feature matrix of the
        representations generated from the input SMILES strings
    invalid_idxs : Set[int]
        the set of indexes in the iterable containing invalid SMILES strings
    """
    if not ray.is_initialized():
        ray.init()

    CHUNKSIZE = 1024
    fps_h5 = Path(path) / f"{name}.h5"

    with h5py.File(str(fps_h5), "w") as h5f:
        fps_dset = h5f.create_dataset(
            "fps", (len(smis), length), "int8", chunks=(CHUNKSIZE, length)
        )

        invalid_idxs = set()
        i = 0
        offset = 0

        ray.put(smis)
        refs = [
            smis_to_fps.remote(smis_chunk, radius, length) for smis_chunk in chunks(smis, CHUNKSIZE)
        ]

        for ref in tqdm(refs, desc="Calculating fingerprints", unit="chunk"):
            fps_chunk = ray.get(ref)
            for fp in fps_chunk:
                while fp is None:
                    invalid_idxs.add(i + offset)
                    offset += 1
                    fp = next(fps_chunk)

                fps_dset[i] = fp
                i += 1

        if len(invalid_idxs) != 0:
            fps_dset.resize(len(smis) - len(invalid_idxs), axis=0)

    return fps_h5, invalid_idxs
