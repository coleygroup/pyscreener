from enum import auto
import functools
from typing import Callable, Optional

import numpy as np
import ray

from pyscreener.utils import AutoName, ScoreMode


class ScreenType(AutoName):
    DOCK = auto()
    VINA = auto()


def run_on_all_nodes(func: Callable) -> Callable:
    """Run a function on all nodes in the ray cluster"""

    @functools.wraps(func)
    def wrapper_run_on_all_nodes(*args, **kwargs):
        refs = []
        for node in ray.nodes():
            address = node["NodeManagerAddress"]
            g = ray.remote(resources={f"node:{address}": 0.1})(func)
            refs.append(g.remote(*args, **kwargs))

        ray.wait(refs)
        return ray.get(refs[-1])

    return wrapper_run_on_all_nodes


def reduce_scores(
    S: np.ndarray,
    repeat_score_mode: ScoreMode = ScoreMode.BEST,
    ensemble_score_mode: ScoreMode = ScoreMode.BEST,
    k: int = 1,
) -> Optional[float]:
    """Calculate the overall score of each ligand given all of its simulations

    Parameters
    ----------
    S : np.ndarray
        an `n x r x t` array of docking scores, where n is the number of ligands that were docked,
        r is the number of receptors each ligand was docked against, and t is the number of repeated
        docking attempts against each receptor, and each value is the docking score calculated for
        the given run
    repeat_score_mode : ScoreMode, default=ScoreMode.BEST,
        the mode used to calculate the overall score for from repeated runs
    ensemble_score_mode : ScoreMode, default=ScoreMode.BEST,
        the mode used to calculate the overall score for a given ensemble of receptors
    k : int, default=1
        the number of scores to consider, if averaging the top-k

    Returns
    -------
    S : np.ndarray
        an array of shape `n`, containing the reduced docking score for each ligand
    """
    if repeat_score_mode == ScoreMode.BEST:
        S = np.nanmin(S, axis=2)
    elif repeat_score_mode == ScoreMode.AVG:
        S = np.nanmean(S, axis=2)
    elif repeat_score_mode == ScoreMode.BOLTZMANN:
        S_e = np.exp(-S)
        Z = S_e / np.nansum(S_e, axis=2)[:, :, None]
        S = np.nansum((S * Z), axis=2)
    elif repeat_score_mode == ScoreMode.TOP_K:
        S = np.nanmean(np.sort(S, axis=2)[:, :k], axis=2)

    if ensemble_score_mode == ScoreMode.BEST:
        S = np.nanmin(S, axis=1)
    elif ensemble_score_mode == ScoreMode.AVG:
        S = np.nanmean(S, axis=1)
    elif ensemble_score_mode == ScoreMode.BOLTZMANN:
        S_e = np.exp(-S)
        Z = S_e / np.nansum(S_e, axis=1)[:, None]
        S = np.nansum((S * Z), axis=1)
    elif ensemble_score_mode == ScoreMode.TOP_K:
        S = np.nanmean(np.sort(S, axis=1)[:, :, :k], axis=1)

    return S
