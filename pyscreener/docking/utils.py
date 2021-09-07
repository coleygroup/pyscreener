from enum import auto, Enum
import functools
from typing import Callable

import ray

class ScreenType(Enum):
    def _generate_next_value_(name, start, count, last_values):
        return name.upper()

    VINA = auto()
    DOCK = auto()

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