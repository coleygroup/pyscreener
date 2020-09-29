from typing import Dict

from pyscreener.docking.docking import dock

def prepare(docker, **kwargs) -> Dict:
    """Prepare all of the inputs for the specified docking program"""
    if docker in {'vina', 'smina', 'psovina', 'qvina'}:
        from pyscreener.docking import vina
        return vina.prepare_inputs(**kwargs)

    if docker == 'dock':
        from pyscreener.docking import ucsfdock
        return ucsfdock.prepare_inputs(**kwargs)

    raise ValueError(f'Unrecognized docking program: "{docker}"')