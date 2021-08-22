from typing import Any, Callable, List, Mapping, Optional, Sequence

import numpy as np
import ray

def run_on_all_nodes(f: Callable[[], Any]) -> List:
    """Run a function on all nodes in the ray cluster

    Parameters
    ----------
    f : Callable[[], Any]
        the function to run

    Returns
    -------
    List
        a list of the function's result on each node
    """        
    refs = []
    for node in ray.nodes():
        address = node["NodeManagerAddress"]
        g = ray.remote(resources={f'node:{address}': 0.1})(f)
        refs.append(g.remote())

    return ray.get(refs)

def calc_ligand_score(
    resultss: Sequence[Sequence[Mapping]],
    receptor_score_mode: str = 'best',
    ensemble_score_mode: str = 'best',
    k: int = 1,
) -> Optional[float]:
    """Calculate the overall score of a ligand given all of its simulations

    Parameters
    ----------
    ligand_results : Sequence[Sequence[Mapping]]
        an MxN list of list of mappings where each individual mapping is the
        result of an individual simulation and
        
        * M is the number of receptors the ligand was docked against
        * N is the number of times each docking run was repeated
    receptor_score_mode : str, default='best'
        the mode used to calculate the overall score for a given receptor
        pose with multiple, repeated runs
    ensemble_score_mode : str, default='best'
        the mode used to calculate the overall score for a given ensemble
        of receptors
    k : int, default=1
        the number of scores to consider, if averaging the top-k
    Returns
    -------
    ensemble_score : Optional[float]
        the overall score of a ligand's ensemble docking. None if no such
        score was calculable
    
    See also
    --------
    calc_score
        for documentation on possible values for receptor_score_mode
        and ensemble_score_mode
    """
    receptor_scores = []
    for results in resultss:
        rep_scores = [
            repeat['score']
            for repeat in results if repeat['score'] is not None
        ]
        if len(rep_scores) > 0:
            receptor_scores.append(calc_score(
                rep_scores, receptor_score_mode, k
            ))

    if len(receptor_scores) > 0:
        ensemble_score = calc_score(
            receptor_scores, ensemble_score_mode, k
        )
    else:
        ensemble_score = None
    
    return ensemble_score

def calc_score(
    scores: Sequence[float], score_mode: str = 'best', k: int = 1
) -> float:
    """Calculate an overall score from a sequence of scores

    Parameters
    ----------
    scores : Sequence[float]
    score_mode : str, default='best'
        the method used to calculate the overall score. Choices include:

        * 'best', 'top': return the top score
        * 'avg', 'mean': return the average of the scores
        * 'top-k': return the average of the top-k scores
        * 'boltzmann': return the boltzmann average of the scores
    k : int, default=1
        the number of top scores to average, if using the top-k average

    Returns
    -------
    float
    """
    Y = np.array(scores)
    if score_mode in ('best', 'top'):
        return Y.min()

    elif score_mode in ('avg', 'mean'):
        return Y.mean()

    elif score_mode == 'boltzmann':
        Y_e = np.exp(-Y)
        Z = Y_e / Y_e.sum()
        return (Y * Z).sum()

    elif score_mode == 'top-k':
        return Y.sort()[:k].mean()
        
    return Y.min()