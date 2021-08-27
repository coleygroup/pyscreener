from enum import Enum, auto
from typing import Any, Callable, List, Mapping, Optional, Sequence

import numpy as np
import ray

class ScoreMode(Enum):
    """The method by which to calculate a score from multiple possible scores.
    Used when calculating an overall docking score from multiple conformations,
    multiple repeated runs, or docking against an ensemble of receptors."""
    AVG = auto()
    BEST = auto()
    BOLTZMANN = auto()
    TOP_K_AVG = auto()

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
    scores: Sequence[float],
    score_mode: ScoreMode = ScoreMode.BEST, k: int = 1
) -> float:
    """Calculate an overall score from a sequence of scores

    Parameters
    ----------
    scores : Sequence[float]
    score_mode : ScoreMode, default=ScoreMode.BEST
        the method used to calculate the overall score. See ScoreMode for
        choices
    k : int, default=1
        the number of top scores to average, if using ScoreMode.TOP_K_AVG

    Returns
    -------
    float
    """
    Y = np.array(scores)

    if score_mode == ScoreMode.BEST:
        return Y.min()
    elif score_mode == ScoreMode.AVG:
        return Y.mean()
    elif score_mode == ScoreMode.BOLTZMANN:
        Y_e = np.exp(-Y)
        Z = Y_e / Y_e.sum()
        return (Y * Z).sum()
    elif score_mode == ScoreMode.TOP_K_AVG:
        return Y.sort()[:k].mean()
        
    return Y.min()