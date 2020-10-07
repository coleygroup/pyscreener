from typing import Dict

from pyscreener.docking.base import Screener
from pyscreener.docking.vina import Vina
from pyscreener.docking.ucsfdock import DOCK

# from pyscreener.docking.docking import dock
# from pyscreener.docking.preparation import prepare

def screener(software, **kwargs):
    if software in ('vina', 'qvina', 'smina', 'psovina'):
        return Vina(software=software, **kwargs)
    
    if software in ('dock', 'ucsfdock', 'DOCK'):
        return DOCK(**kwargs)
    
    raise ValueError(f'Unrecognized docking software: "{software}"')