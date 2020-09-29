"""This module contains the singular function dock, which is a template
function for performing docking runs over a distributed system"""
from collections import defaultdict
import datetime
from itertools import chain
from math import ceil
import os
from os import PathLike
import timeit
from typing import Dict, List, Tuple

from pyscreener.utils import calc_score

def dock(docker: str, inputs: Dict,
         path: str = './docking_results_'+str(datetime.date.today()),
         score_mode: str = 'best', repeats: int = 1,
         repeat_score_mode: str = 'best', ensemble_score_mode: str = 'best',
         distributed: bool = False, ncpu: int = 1, num_workers: int = -1, 
         verbose: int = 0, **kwargs) -> Tuple[Dict[str, float], List[Dict]]:
    """Run the specified docking program with the given inputs

    Parameters
    ----------
    docker : str
        the docking program that will be used
    inputs : Input
        a dictionary docking program inputs
    path : string (Default = './docking_results_YYYY-MM-DD/')
        the path under which docking outputs should be organized
    score_mode : str (Default = 'best')
        the method by which to calculate the score of a docking run
    repeats : int (Default = 1)
        the number of times a docking run should be repeated
    repeat_score_mode : str (Default = 'best')
        the method used to calculate the overall docking score of a molecule
        for repeated runs
    ensemble_score_mode : str (Default = 'best')
        the method by which to calculate the overall docking score of a
        molecule in an ensemble of docking runs (multiple structures)
    distributed: bool (Default = True)
        whether the work should be performed over a distributed system. 
        Setting this value to false will result in all simulations being
        performed on the local machine
    ncpu : int (Default = 1)
        the number of cores available to each worker. Only used if distributed
        is False.
    num_workers : int (Default = -1)
        the number of workers to distribute docking simulations over. A value
        of -1 will use all available cores. This argument only needs to be set
        if distributed is False.
    verbose : int (Default = 0)
        whether or not to print performance data
    **kwargs
        additional and unused keyword arguments

    Returns
    -------
    d_smi_score : Dict[str, float]
        a dictionary mapping SMILES string to the best score among the
        corresponding ligands. Does not include ligands for which all
        inputs failed to dock.
    rows : List[Dict]
        A dataframe containing the following information for each 
        docking run:
            smiles  - the ligand's SMILES string
            name    - the name of the ligand
            in      - the filepath of the input ligand file
            out     - the filepath of the output docked ligand
            log     - the filepath of the log file
            score   - the ligand's docking score
    """
    begin = timeit.default_timer()

    if distributed:
        from mpi4py import MPI
        from mpi4py.futures import MPIPoolExecutor as Pool

        num_workers = MPI.COMM_WORLD.size
    else:
        from concurrent.futures import ProcessPoolExecutor as Pool
        if num_workers == -1:
            try:
                num_workers = len(os.sched_getaffinity(0)) // ncpu
            except AttributeError:
                num_workers = os.cpu_count()

    # BATCHES_PER_PROCESS = 32
    # batch_size = ceil(size / (BATCHES_PER_PROCESS*n_workers))
    CHUNKSIZE = 32

    if docker in {'vina', 'smina', 'psovina', 'qvina'}:
        from pyscreener.docking.vina import dock_inputs
    elif docker == 'dock':
        from pyscreener.docking.ucsfdock import dock_inputs
    else:
        raise ValueError(f'Unrecognized docking program: "{docker}"')

    with Pool(max_workers=num_workers) as client:
        rowsss = dock_inputs(
            docker=docker, **inputs, path=path, score_mode=score_mode, 
            repeats=repeats, chunksize=CHUNKSIZE, client=client
        )
    
    # record only one score for a given SMILES string across all docking runs
    d_smi_score = defaultdict(lambda: float('inf'))
    for ligand_rowss in rowsss:
        if len(ligand_rowss) == 0:
            continue

        receptor_scores = [calc_score(
            [row['score'] for row in receptor_rows], repeat_score_mode
        ) for receptor_rows in ligand_rowss]
        ensemble_score = calc_score(receptor_scores, ensemble_score_mode)

        smi = ligand_rowss[0][0]['smiles']
        curr_score = d_smi_score[smi]
        d_smi_score[smi] = min(curr_score, ensemble_score)
    d_smi_score = dict(d_smi_score)

    rows = list(chain(*list(chain(*rowsss))))

    total = timeit.default_timer() - begin
    m, s = divmod(int(total), 60)
    h, m = divmod(m, 60)
    if verbose > 0 and len(rowsss) > 0:
        print(f'  Time to dock {len(rowsss)} ligands:',
              f'{h:d}h {m:d}m {s:d}s ({total/len(rowsss):0.3f} s/ligand)', 
              flush=True)

    return d_smi_score, rows
