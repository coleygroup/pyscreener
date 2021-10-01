from .data import CalculationData
from .result import Result
from .runner import DockingRunner
from .screen import DockingVirtualScreen
from .utils import ScreenType

def screener(software, *args, **kwargs):
    software = software.lower()

    if software.lower() in ('vina', 'qvina', 'smina', 'psovina'):
        from pyscreener.docking.vina import VinaRunner
        runner = VinaRunner
        screen_type = ScreenType.VINA
    elif software.lower() in ('dock', 'ucsfdock'):
        from pyscreener.docking.dock import DOCKRunner
        runner = DOCKRunner
        screen_type = ScreenType.DOCK
    else:
        raise ValueError(f'Unrecognized docking software: "{software}"')

    return DockingVirtualScreen(runner, *args, **kwargs)