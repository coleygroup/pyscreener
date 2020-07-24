"""Utility functions"""

from math import exp
from typing import Sequence, TypeVar

def calc_score(scores: Sequence[float], score_mode: str = 'best') -> float:
    """
    Calculate an overall score from a sequence of scores

    Parameters
    ----------
    scores : Sequence[float]
    score_mode : str (Default = 'best')
        the method used to calculate the overall score
        Choices:
            'best' - return the top docking score
            'avg' - return the average of all the docking scores
            'boltzmann' - return the boltzmann average of all the docking scores

    Returns
    -------
    score : float
    """
    if score_mode not in {'best', 'avg', 'boltzmann'}:
        raise ValueError(f'Unrecognized score mode: "{score_mode}"')

    if score_mode == 'best':
        score = scores[0]
    elif score_mode == 'avg':
        score = sum(score for score in scores) / len(scores)
    elif score_mode == 'boltzmann':
        Z = sum(exp(-score) for score in scores)
        score = sum(score * exp(-score) / Z for score in scores)
    
    return score

def time_stats(total_time, n):
    m, s = divmod(int(total_time), 60)
    h, m = divmod(m, 60)
    avg = total_time / n
    
    return h, m, s, avg