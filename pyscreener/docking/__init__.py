from .data import CalculationData
from .metadata import CalculationMetadata
from .result import Result
from .runner import DockingRunner
from .screen import VirtualScreen

def screener(software, **kwargs):
    if software in ('vina', 'qvina', 'smina', 'psovina'):
        from pyscreener.docking.vina import Vina
        return Vina(software=software, **kwargs)
    
    if software in ('dock', 'ucsfdock', 'DOCK'):
        from pyscreener.docking.dock import DOCK
        return DOCK(**kwargs)
    
    raise ValueError(f'Unrecognized docking software: "{software}"')