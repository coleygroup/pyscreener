from typing import Dict

from pyscreener.docking import dock

def prepare(docker, **kwargs) -> Dict:
    """Prepare all of the inputs for the specified docking program"""
    if docker in {'vina', 'smina', 'psovina', 'qvina'}:
        from . import vina
        return vina.prepare_inputs(docker=docker, **kwargs)

    if docker == 'dock':
        from . import ucsfdock
        return ucsfdock.prepare_inputs(**kwargs)

    raise ValueError(f'Unrecognized docking program: "{docker}"')