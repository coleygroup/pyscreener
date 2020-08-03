import csv
from functools import partial
import gzip
import os
from pathlib import Path
# import sys
# import timeit
from typing import Iterable, Optional, Set, Tuple

import h5py
import numpy as np
from rdkit.Chem import AllChem as Chem
from tqdm import tqdm

def smi_to_fp(smi: str, radius: int = 2,
              length: int = 2048) -> Optional[np.ndarray]:
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None

    return np.array(Chem.GetMorganFingerprintAsBitVect(
        mol, radius, nBits=length, useChirality=True))

def gen_fps_h5(smis: Iterable[str], n_mols: Optional[int] = None,
               path: str = '.', name: str = 'fps',
               radius: int = 2, length: int = 2048,
               distributed: bool = True,
               n_workers: int = -1) -> Tuple[str, Set[int]]:
    """Generate an hdf5 file containing the feature matrix of the list of
    SMILES strings

    Parameters
    ----------
    smis : Iterable[str]
        the SMILES strings from which to generate fingerprints
    n_mols : Optional[int]
        the length of the iterable
    path : str (Default = '.')
        the path under which the H5 file should be written
    name : str (Default = 'fps')
        the name of the output H5 file
    radius : int (Default = 2)
        the radius of the fingerprints
    length : int (Default = 2048)
        the length of the fingerprints
    distributed : bool (Default = True)
        whether to parallelize fingerprint calculation over a distributed
        computing setup
    n_workers : int (Default = -1)
        how many jobs to parellize file parsing over.
        A value of -1 defaults to using all cores

    Returns
    -------
    fps_h5 : str
        the filename of the hdf5 file containing the feature matrix of the
        representations generated from the input SMILES strings
    invalid_idxs : Set[int]
        the set of indexes in the iterable containing invalid SMILES strings
    """
    if distributed:
        from mpi4py import MPI
        from mpi4py.futures import MPIPoolExecutor as Pool

        n_workers = MPI.COMM_WORLD.size
    else:
        from concurrent.futures import ProcessPoolExecutor as Pool
        if n_workers == -1:
            try:
                n_workers = len(os.sched_getaffinity(0))
            except AttributeError:
                n_workers = os.cpu_count()

    fps_h5 = f'{path}/{name}.h5'
    compression = None

    with Pool(max_workers=n_workers) as pool, h5py.File(fps_h5, 'w') as h5f:
        CHUNKSIZE = 1024
        fps_dset = h5f.create_dataset(
            'fps', (n_mols, length), compression=compression,
            chunks=(CHUNKSIZE, length), dtype='int8'
        )
        # semaphore = mp.Semaphore((njobs+2) * chunksize)
        # def gen_rows(reader: csv.reader, sem: mp.Semaphore):
        #     for row in reader:
        #         sem.acquire()
        #         yield row
        # rows = gen_rows(reader, semaphore)
        smi_to_fp_partial = partial(smi_to_fp, radius=radius, length=length)

        invalid_idxs = set()
        offset = 0

        fps = pool.map(smi_to_fp_partial, smis, chunksize=CHUNKSIZE)
        for i, fp in tqdm(enumerate(fps), total=n_mols,
                          desc='Calculating fingerprints', unit='fp'):
            while fp is None:
                invalid_idxs.add(i+offset)
                offset += 1
                fp = next(fps)
            
            fps_dset[i] = fp
            # semaphore.release()
        # original dataset size included potentially invalid SMILES
        n_mols_valid = n_mols - len(invalid_idxs)
        if n_mols_valid != n_mols:
            fps_dset.resize(n_mols_valid, axis=0)

    return fps_h5, invalid_idxs

def gen_fps_h5_from_file(filepath: str, delimiter: str = ',',
                         title_line: bool = True, smiles_col: int = 0,
                         path: Optional[str] = None, name: Optional[str] = None,
                         **kwargs) -> Tuple[str, Set[int]]:
    """Parses a CSV file containing SMILES strings to generate an hdf5 file 
    containing the feature matrix of the library.

    Parameters
    ----------
    filepath : str
        the filepath of a (compressed) CSV file containing the SMILES strings
        for which to generate fingerprints
    delimiter : str (Default = ',')
        the column separator for each row
    title_line : bool (Default = True)
        does the file contain a title line?
    smiles_col : int (Default = 0)
        the index of the column containing the SMILES string of the molecule
    path : Optional[str] (Default = None)
        the path under which the H5 file should be written
    name : Optional[str] (Default = None)
        the name of the output H5 file
    **kwargs
        keyword arguments to gen_fps_h5_from_smis

    Returns
    -------
    fps_h5 : str
        the filename of an hdf5 file containing the feature matrix of the
        representations generated from the molecules in the input file.
        The row ordering corresponds to the ordering of smis
    invalid_rows : Set[int]
        the set of rows in filepath containing invalid SMILES strings
    
    See also
    --------
    gen_fps_h5_from_smis
    """
    if os.stat(filepath).st_size == 0:
        raise ValueError(f'"{filepath} is empty!"')

    path = path or '.'
    name = name or Path(filepath).stem

    if Path(filepath).suffix == '.gz':
        open_ = partial(gzip.open, mode='rt')
    else:
        open_ = open

    with open_(filepath) as fid:
        reader = csv.reader(fid, delimiter=delimiter)
        n_mols = sum(1 for _ in reader); fid.seek(0)
        if title_line:
            next(reader)
            n_mols -= 1

        smis = (row[smiles_col] for row in reader)
        fps_h5, invalid_rows = gen_fps_h5(
            smis, n_mols, path, name, **kwargs)

    return fps_h5, invalid_rows

# def main():
#     filepath = sys.argv[1]
#     fps_h5 = parse_smiles_par(filepath, njobs=sys.argv[2])
#     print(fps_h5)

# if __name__ == '__main__':
#     main()
