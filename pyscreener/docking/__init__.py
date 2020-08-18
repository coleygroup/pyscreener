from typing import Dict

from . import vina, ucsfdock
from .docking import dock

def prepare(docker, **kwargs) -> Dict:
    """Prepare all of the inputs for the specified docking program"""
    if docker in {'vina', 'smina', 'psovina', 'qvina'}:
        return vina.prepare_inputs(docker=docker, **kwargs)

    if docker == 'dock':
        raise NotImplementedError

    if docker == 'rdock':
        raise NotImplementedError

    raise ValueError(f'Unrecognized docking program: "{docker}"')