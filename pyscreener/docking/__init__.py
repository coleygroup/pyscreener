from .data import CalculationData
from .screen import DockingVirtualScreen

def screener(software, **kwargs):
    software = software.lower()

    if software.lower() in ('vina', 'qvina', 'smina', 'psovina'):
        from pyscreener.docking.vina import VinaRunner
        runner = VinaRunner
    elif software.lower() in ('dock', 'ucsfdock'):
        from pyscreener.docking.dock import DOCKRunner
        runner = DOCKRunner
    else:
        raise ValueError(f'Unrecognized docking software: "{software}"')

    return DockingVirtualScreen(runner, **kwargs)