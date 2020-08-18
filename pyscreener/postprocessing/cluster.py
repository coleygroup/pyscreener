"""This module contains functions for clutering a set of molecules"""

import csv
from itertools import chain
import os
from pathlib import Path
from random import sample
import sys
import timeit
from typing import Dict, Iterable, List, Optional

import h5py
import numpy as np
from scipy import sparse
from sklearn.cluster import MiniBatchKMeans

from . import fingerprints

def cluster(d_smi_score: Dict[str, Optional[float]],
            path: str = '.', **kwargs) -> List[Dict]:
    d_smi_cid = cluster_smis(d_smi_score.keys(), path=path, **kwargs)

    smi_score_clusters = {}
    for smi, score in d_smi_score.items():
        cid = d_smi_cid[smi]
        smi_score_clusters[cid][smi] = score

    return list(smi_score_clusters.values())

def cluster_smis(smis: Iterable[str], n_cluster: int = 100, 
                 path: str = '.', name: str = 'fps', 
                 **kwargs) -> Dict[str, int]:
    """Cluster the SMILES strings

    Parameters
    ----------
    smis : Iterable[str]
        the SMILES strings to cluster
    n_cluster : int (Default = 100)
        the number of clusters to generate
    path : str (Default = '.')
        the path under which to write the fingerprint file
    name : str (Default = '.')
        the name of the output fingerprint file
    **kwargs
        keyword arguments to fingerprints.gen_fps_h5

    Returns
    -------
    d_smi_cid : Dict[str, int]
        a mapping from SMILES string to cluster ID
    
    See also
    --------
    fingerprints.gen_fps_h5
    """
    fps_h5, invalid_idxs = fingerprints.gen_fps_h5(
        smis, path=path, name=name, **kwargs)

    smis = [smi for i, smi in enumerate(smis) if i not in invalid_idxs]
    cids = cluster_fps_h5(fps_h5, n_cluster)
    
    return dict(zip(smis, cids))

def cluster_file(filepath: str, delimiter: str = ',',
                 title_line: bool = True, smiles_col: int = 0,
                 distributed: bool = True, n_workers: int = -1, 
                 n_cluster: int = 100, path: Optional[str] = None):
    """
    Cluster the molecules

    Parameters
    ----------
    filepath : str
        the filepath of the .smi file to cluster
    smiles_col : int (Default = 0)
        the column containing the SMILES string in each row
    title_line : bool (Default = True)
        does the file contain a title line?
    n_workers : int (Default = -1)
        the number of jobs to parallelize fingerprint preparation over. 
        A value of -1 uses all available cores
    n_cluster : int (Default = 100)
        the number of clusters to group the molecules into
    path : Optional[str]
        the path under which all clustering files should be written.
    """
    start = timeit.default_timer()

    fps_h5, invalid_rows = fingerprints.gen_fps_h5_from_file(
        filepath, delimiter=delimiter, smiles_col=smiles_col,
        title_line=title_line, n_workers=n_workers, path=path)

    cluster_ids = cluster_fps_h5(fps_h5, n_cluster=n_cluster)

    if path is None:
        path = Path(filepath)
        path = path.with_name(f'{path.stem}_clusters')

    with open(filepath) as fid_in, open(path, 'w') as fid_out:
        reader = csv.reader(fid_in)
        if title_line:
            next(reader)

        offset = 0
        fid_out.write('smiles,cluster\n')
        for i, row in enumerate(reader):
            while (i+offset) in invalid_rows:
                offset += 1
                row = next(reader)
            smi = row[smiles_col]
            fid_out.write(f'{smi},{cluster_ids[i]}\n')

    elapsed = timeit.default_timer() - start
    print(f'Total time to cluster the {len(cluster_ids)} mols in "{filepath}"',
          f'over: {elapsed:0.3f}s')

def cluster_fps_h5(fps_h5: str, n_cluster: int = 100) -> List[int]:
    """Cluster the feature matrix of fingerprints in fps_h5

    Parameters
    ----------
    fps : str
        the filepath of an h5py file containing the NxM matrix of
        molecular fingerprints, where N is the number of molecules and
        M is the length of the fingerprint (feature representation)
    ncluster : int (Default = 100)
        the number of clusters to form with the given fingerprints (if the
        input method requires this parameter)

    Returns
    -------
    cids : List[int]
        the cluster id corresponding to a given fingerprint
    """
    begin = timeit.default_timer()

    BATCH_SIZE = 1000
    ITER = 1000

    clusterer = MiniBatchKMeans(n_clusters=n_cluster, batch_size=BATCH_SIZE)

    with h5py.File(fps_h5, 'r') as h5f:
        fps = h5f['fps']
        chunk_size = fps.chunks[0]

        for _ in range(ITER):
            rand_idxs = sorted(sample(range(len(fps)), BATCH_SIZE))
            batch_fps = fps[rand_idxs]
            clusterer.partial_fit(batch_fps)

        cidss = [clusterer.predict(fps[i:i+chunk_size])
                 for i in range(0, len(fps), chunk_size)]

    elapsed = timeit.default_timer() - begin

    print(f'Clustering took: {elapsed:0.3f}s')

    return list(chain(*cidss))

# def write_clusters(path: str, d_cluster_smidxs: Dict[int, List[int]],
#                    smis: h5py.Dataset) -> str:
#     """Write each cluster to a separate file under the given path"""
#     p = Path(path)
#     if not p.is_dir():
#         p.mkdir(parents=True)

#     for _, smidxs in d_cluster_smidxs.items():
#         write_cluster(str(p / f'cluster_{id:003d}.smi'), smidxs, smis)

# def write_cluster(filepath: str, smidxs: List[int],
#                   smis: h5py.Dataset) -> None:
#     """Write each string in smis to filepath"""
#     with open(filepath, 'w', newline='') as f:
#         f.write('smiles\n')
#         for idx in smidxs:
#             smi = smis[idx]
#             f.write(f'{smi}\n')

# def cluster_fps(fps: np.ndarray,
#                 n_cluster: int = 100, method: str = 'minibatch',
#                 n_workers: Optional[int] = None) -> np.ndarray:
#     """
#     Cluster the molecular fingerprints, fps, by a given method

#     Parameters
#     ----------
#     fps : np.ndarray
#         the feature matrix of the molecules to cluster
#     n_cluster : int (Default = 100)
#         the number of clusters to form with the given fingerprints (if the
#         input method requires this parameter)
#     method : str (Default = 'kmeans')
#         the clusering method to use.
#         Choices include:
#         - 'kmeans': k-means clustering
#         - 'minibatch': mini-batch k-means clustering: 
#         - 'optics': OPTICS clustering
#     n_workers : Optional[int]
#         the number of jobs to parallelize clustering over, if possible

#     Returns
#     -------
#     cluster_ids : np.ndarray
#         the cluster id corresponding to a given fingerprint
#     """
#     begin = timeit.default_timer()

#     fps = sparse.vstack(fps, format='csr')

#     if method == 'kmeans':
#         clusterer = cluster.KMeans(n_clusters=n_cluster, n_jobs=n_workers)
#     elif method == 'minibatch':
#         clusterer = cluster.MiniBatchKMeans(
#             n_clusters=n_cluster, n_init=10, batch_size=100, init_size=1000)
#     elif method == 'optics':
#         clusterer = cluster.OPTICS(
#             min_samples=0.01, metric='jaccard', n_jobs=n_workers)
#     else:
#         raise ValueError(f'Unrecognized clustering method: "{method}"')

#     cluster_ids = clusterer.fit_predict(fps)

#     elapsed = timeit.default_timer() - begin
#     print(f'Clustering and predictions took: {elapsed:0.3f}s')

#     return cluster_ids

# if __name__ == '__main__':
#     cluster_smiles_file(sys.argv[1], njobs=MAX_CPU)
