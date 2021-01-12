from typing import Dict

from pyscreener.docking.base import Screener

def screener(software, **kwargs):
    if software in ('vina', 'qvina', 'smina', 'psovina'):
        from pyscreener.docking.vina import Vina
        return Vina(software=software, **kwargs)
    
    if software in ('dock', 'ucsfdock', 'DOCK'):
        from pyscreener.docking.dock import DOCK
        return DOCK(**kwargs)
    
    raise ValueError(f'Unrecognized docking software: "{software}"')