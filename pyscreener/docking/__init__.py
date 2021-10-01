from dataclasses import asdict

from .data import CalculationData
from .metadata import CalculationMetadata
from .result import Result
from .runner import DockingRunner
from .screen import DockingVirtualScreen
from .utils import ScreenType

def screen_type(software) -> ScreenType:
    if software.lower() in ('vina', 'qvina', 'smina', 'psovina'):
        return ScreenType.VINA
    elif software.lower() in ('dock', 'dock6', 'ucsfdock'):
        return ScreenType.DOCK
    else:
        raise ValueError(f'Unrecognized docking software: "{software}"')

def build_metadata(software: str, **kwargs) -> CalculationMetadata:
    if software.lower() in ('vina', 'qvina', 'smina', 'psovina'):
        from pyscreener.docking.vina.metadata import VinaMetadata

        d_md = asdict(VinaMetadata())
        d_md.update((k, kwargs[k]) for k in d_md.keys() & kwargs.keys())

        return VinaMetadata(**d_md)

    if software.lower() in ('dock', 'ucsfdock'):
        from pyscreener.docking.dock.metadata import DOCKMetadata

        d_md = asdict(DOCKMetadata())
        d_md.update((k, kwargs[k]) for k in d_md.keys() & kwargs.keys())

        return DOCKMetadata(**d_md)

    raise ValueError(f'Unrecognized docking software: "{software}"')

def virtual_screen(software: str, *args, **kwargs):
    if software.lower() in ('vina', 'qvina', 'smina', 'psovina'):
        from pyscreener.docking.vina import VinaRunner

        runner = VinaRunner
    elif software.lower() in ('dock', 'ucsfdock'):
        from pyscreener.docking.dock import DOCKRunner

        runner = DOCKRunner
    else:
        raise ValueError(f'Unrecognized docking software: "{software}"')

    return DockingVirtualScreen(runner, *args, **kwargs)