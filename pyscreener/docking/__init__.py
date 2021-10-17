from dataclasses import asdict
from typing import Dict, Optional

from pyscreener.exceptions import (
    MisconfiguredDirectoryError,
    MissingEnvironmentVariableError,
    MissingExecutableError,
    UnsupportedSoftwareError,
)

from .data import CalculationData
from .metadata import CalculationMetadata
from .result import Result
from .runner import DockingRunner
from .screen import DockingVirtualScreen
from .utils import ScreenType



def build_metadata(software: str, metadata: Optional[Dict] = None) -> CalculationMetadata:
    metadata = metadata or {}

    if software.lower() in ("vina", "qvina", "smina", "psovina"):
        from pyscreener.docking.vina.metadata import VinaMetadata

        d_md = asdict(VinaMetadata())
        d_md.update((k, metadata[k]) for k in d_md.keys() & metadata.keys())

        return VinaMetadata(**d_md)

    if software.lower() in ("dock", "dock6", "ucsfdock"):
        from pyscreener.docking.dock.metadata import DOCKMetadata

        d_md = asdict(DOCKMetadata())
        d_md.update((k, metadata[k]) for k in d_md.keys() & metadata.keys())

        return DOCKMetadata(**d_md)

    raise ValueError(f'Unrecognized docking software: "{software}"')


def get_runner(software: str) -> DockingRunner:
    if software.lower() in ("vina", "qvina", "smina", "psovina"):
        from pyscreener.docking.vina import VinaRunner
        return VinaRunner

    if software.lower() in ("dock", "dock6", "ucsfdock"):
        from pyscreener.docking.dock import DOCKRunner
        return DOCKRunner

    raise ValueError(f'Unrecognized docking software: "{software}"')


def check_env(software, metadata: Optional[Dict] = None):
    print("Checking environment for input screen ...", end=" ")
    try:
        runner = get_runner(software)
        metadata = build_metadata(software, metadata)
        runner.validate_metadata(metadata)
    except (
        MissingExecutableError,
        MisconfiguredDirectoryError,
        MissingEnvironmentVariableError,
    ):
        print("FAIL")
        print("Environment not set up properly! See the exception for more details", flush=True)
        raise

    print("PASS")
    print("Environment is properly set up!")


def virtual_screen(software: str, *args, **kwargs) -> DockingVirtualScreen:
    return DockingVirtualScreen(get_runner(software), *args, **kwargs)
