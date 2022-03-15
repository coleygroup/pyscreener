__all__ = [
    "AutoName",
    "Reduction",
    "FileFormat",
    "chunks",
    "calc_score",
    "reduce_scores",
    "run_on_all_nodes",
]

from enum import Enum, auto
import functools
from itertools import islice
from typing import Callable, Iterable, Iterator, List, Optional, Sequence
import warnings

import numpy as np
import ray


class AutoName(Enum):
    def _generate_next_value_(name, start, count, last_values):
        return name

    @classmethod
    def from_str(cls, s):
        return cls[s.replace("-", "_").upper()]


class Reduction(AutoName):
    """The method by which to reduce multiple scores into a single score.
    Used when calculating an overall docking score from multiple conformations,
    multiple repeated runs, or docking against an ensemble of receptors."""

    AVG = auto()
    BEST = auto()
    BOLTZMANN = auto()
    TOP_K = auto()


class FileFormat(AutoName):
    """The format of a molecular suppy file. FILE represents the format of all molecular supply
    files with no explicit support (i.e., CSV, SDF, and SMI.)"""

    CSV = auto()
    FILE = auto()
    SDF = auto()
    SMI = auto()


def chunks(it: Iterable, size: int) -> Iterator[List]:
    """chunk an iterable into chunks of given size, with the last chunk being potentially smaller"""
    it = iter(it)
    return iter(lambda: list(islice(it, size)), [])


def calc_score(
    scores: Sequence[float], score_mode: Reduction = Reduction.BEST, k: int = 1
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

    if score_mode == Reduction.BEST:
        return Y.min()
    elif score_mode == Reduction.AVG:
        return np.nanmean(Y)
    elif score_mode == Reduction.BOLTZMANN:
        Y_e = np.exp(-Y)
        Z = Y_e / np.nansum(Y_e)
        return np.nansum(Y * Z)
    elif score_mode == Reduction.TOP_K:
        return np.nanmean(Y.sort()[:k])

    raise ValueError(f"Invalid ScoreMode! got: {score_mode}")


def reduce_scores(
    S: np.ndarray, reduction: Reduction = Reduction.BEST, k: int = 1
) -> Optional[float]:
    """Calculate the overall score of each ligand given all of its simulations

    Parameters
    ----------
    S : np.ndarray
        an `n x r` array of docking scores, where n is the number of ligands that were docked,
        r is the number of receptors each ligand was docked against, and t is the number of repeated
        docking attempts against each receptor, and each value is the docking score calculated for
        the given run
    ensemble_score_mode : ScoreMode, default=ScoreMode.BEST,
        the mode used to calculate the overall score for a given ensemble of receptors
    k : int, default=1
        the number of scores to consider, if averaging the top-k

    Returns
    -------
    S : np.ndarray
        an array of shape `n` containing the reduced docking score for each ligand
    """
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", r"All-NaN (slice|axis) encountered")

        if reduction == Reduction.BEST:
            S = np.nanmin(S, axis=1)
        elif reduction == Reduction.AVG:
            S = np.nanmean(S, axis=1)
        elif reduction == Reduction.BOLTZMANN:
            S_e = np.exp(-S)
            Z = S_e / np.nansum(S_e, axis=1)[:, None]
            S = np.nansum((S * Z), axis=1)
        elif reduction == Reduction.TOP_K:
            S = np.nanmean(np.sort(S, axis=1)[:, :, :k], axis=1)

    return S


def run_on_all_nodes(func: Callable) -> Callable:
    """A decorator to run a function on all nodes in a ray cluster.

    Ex:
    >>> @run_on_all_nodes
    >>> def f():
    ...     print("hello world")

    will calling the function `f()` will print "hello world" from each node in the ray cluster
    """

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
