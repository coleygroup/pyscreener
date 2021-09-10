from enum import auto, Enum
import functools
from pyscreener.docking.data import CalculationData
from pyscreener.utils import ScoreMode, calc_score
from typing import Callable, Optional, Sequence

import ray

# class ScreenType(Enum):
#     def _generate_next_value_(name, start, count, last_values):
#         return name.upper()

#     VINA = auto()
#     DOCK = auto()

def run_on_all_nodes(func: Callable) -> Callable:
    """Run a function on all nodes in the ray cluster"""
    @functools.wraps(func)
    def wrapper_run_on_all_nodes(*args, **kwargs):
        refs = []
        for node in ray.nodes():
            address = node["NodeManagerAddress"]
            g = ray.remote(resources={f'node:{address}': 0.1})(func)
            refs.append(g.remote(*args, **kwargs))

        ray.wait(refs)
        return ray.get(refs[-1])
    
    return wrapper_run_on_all_nodes

def calc_ligand_score(
    receptors_repeats: Sequence[Sequence[CalculationData]],
    receptor_score_mode: ScoreMode = ScoreMode.BEST,
    ensemble_score_mode: ScoreMode = ScoreMode.BEST,
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
    receptor_score_mode : ScoreMode, default=ScoreMode.BEST,
        the mode used to calculate the overall score for a given receptor
        pose with multiple, repeated runs
    ensemble_score_mode : ScoreMode, default=ScoreMode.BEST,
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
    for repeats in receptors_repeats:
        rep_scores = [
            rep.result.score
            for rep in repeats if rep.result.score is not None
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