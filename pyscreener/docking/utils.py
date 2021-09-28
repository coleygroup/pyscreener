from dataclasses import asdict
from enum import auto
import functools
import shutil
from typing import Callable, Optional

import numpy as np
import ray

from pyscreener.exceptions import MissingExecutableError
from pyscreener.utils import AutoName, ScoreMode
from pyscreener.docking.metadata import CalculationMetadata


class ScreenType(AutoName):
    DOCK = auto()
    VINA = auto()


def build_metadata(screen_type: ScreenType, **kwargs) -> CalculationMetadata:
    if screen_type == ScreenType.DOCK:
        from pyscreener.docking import dock

        d_md = asdict(dock.DOCKMetadata())
        d_md.update((k, kwargs[k]) for k in d_md.keys() & kwargs.keys())

        return dock.DOCKMetadata(**d_md)

    elif screen_type == ScreenType.VINA:
        from pyscreener.docking import vina

        d_md = asdict(vina.VinaMetadata())
        d_md.update((k, kwargs[k]) for k in d_md.keys() & kwargs.keys())

        return vina.VinaMetadata(**d_md)

    raise ValueError(f"Invalid screen type specified! got: {screen_type}.")


def valiate_metadata(screen_type: ScreenType, metadata: CalculationMetadata):
    if screen_type == ScreenType.DOCK:
        return
    else:
        if shutil.which(metadata.software.value) is None:
            raise MissingExecutableError(
                f'Could not find "{metadata.software.value}" on PATH! '
                "See https://github.com/coleygroup/pyscreener/tree/refactor#adding-an-executable-to-your-path for more information."
            )


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
