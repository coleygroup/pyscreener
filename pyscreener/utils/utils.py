__all__ = [
    "AutoName",
    "Reduction",
    "FileFormat",
    "chunks",
    "reduce_scores",
    "reduce_scores",
    "run_on_all_nodes",
]

from enum import Enum, auto
import functools
from itertools import islice
from typing import Callable, Iterable, Iterator, List, Optional

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


def reduce_scores(
    S: np.ndarray, reduction: Reduction = Reduction.BEST, axis: int = -1, k: int = 1
) -> Optional[float]:
    """Calculate the overall score of each ligand given all its scores against multiple receptors

    Parameters
    ----------
    S : np.ndarray
        an array of docking scores
    reduction : Reduction, default=Reduction.BEST,
        the mode used to calculate the overall score for a given ensemble of receptors
    axis : int, default=-1
        the axis along which to reduce
    k : int, default=1
        the number of scores to consider, if using a TOP_K reduction

    Returns
    -------
    S : np.ndarray
        an array of shape `n` containing the reduced docking score for each ligand

    Raises
    ------
    ValueError
        if an invalid `reduction` was passed
    """
    if np.isnan(S).all():
        return S.sum(axis)

    if reduction == Reduction.BEST:
        return np.nanmin(S, axis)
    elif reduction == Reduction.AVG:
        return np.nanmean(S, axis)
    elif reduction == Reduction.BOLTZMANN:
        S_e = np.exp(-S)
        Z = S_e / np.nansum(S_e, axis, keepdims=True)
        return np.nansum((S * Z), axis)
    elif reduction == Reduction.TOP_K:
        return np.nanmean(np.sort(S, axis)[..., :k], axis)

    raise ValueError(f"Invalid reduction specified! got: {reduction}")


def run_on_all_nodes(func: Callable) -> Callable:
    """A decorator to run a function on all nodes in a ray cluster.

    Ex:
    >>> @run_on_all_nodes
    >>> def f():
    ...     print("hello world")

    will call the function `f()` will print "hello world" from each node in the ray cluster
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
