from typing import Any, Callable, List

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