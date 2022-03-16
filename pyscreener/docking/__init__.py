from dataclasses import asdict
from typing import Dict, Optional

from colorama import init, Fore, Style

from pyscreener.exceptions import (
    MisconfiguredDirectoryError,
    MissingEnvironmentVariableError,
    MissingExecutableError,
    UnsupportedSoftwareError,
)
from .sim import Simulation
from .metadata import SimulationMetadata
from .result import Result
from .runner import DockingRunner
from .screen import DockingVirtualScreen
from .utils import ScreenType

init(autoreset=True)


def build_metadata(software: str, metadata: Optional[Dict] = None) -> SimulationMetadata:
    metadata = metadata or {}

    if software.lower() in ("vina", "qvina", "smina", "psovina"):
        from pyscreener.docking.vina import VinaMetadata

        d_md = asdict(VinaMetadata())
        d_md.update((k, metadata[k]) for k in d_md.keys() & metadata.keys())

        return VinaMetadata(**d_md)

    if software.lower() in ("dock", "dock6", "ucsfdock"):
        from pyscreener.docking.dock.metadata import DOCKMetadata

        d_md = asdict(DOCKMetadata())
        d_md.update((k, metadata[k]) for k in d_md.keys() & metadata.keys())

        return DOCKMetadata(**d_md)

    raise UnsupportedSoftwareError(f'Unrecognized screen type: "{software}"')


def get_runner(software: str) -> DockingRunner:
    if software.lower() in ("vina", "qvina", "smina", "psovina"):
        from pyscreener.docking.vina import VinaRunner

        return VinaRunner

    if software.lower() in ("dock", "dock6", "ucsfdock"):
        from pyscreener.docking.dock import DOCKRunner

        return DOCKRunner

    raise UnsupportedSoftwareError(f'Unrecognized screen type: "{software}"')


def check_env(software, metadata: Optional[Dict] = None):
    print(f'Checking environment and metadata for "{software}" screen')
    try:
        valid_env = False
        print("  Checking PATH and environment variables ...", end=" ")
        runner = get_runner(software)
        print(Style.BRIGHT + Fore.GREEN + "PASS")
        valid_env = True
        print("  Validating metadata ... ", end=" ")
        metadata = build_metadata(software, metadata)
        runner.validate_metadata(metadata)
        print(Style.BRIGHT + Fore.GREEN + "PASS")
    except (
        MisconfiguredDirectoryError,
        MissingEnvironmentVariableError,
        MissingExecutableError,
        UnsupportedSoftwareError,
    ):
        print(Style.BRIGHT + Fore.RED + "FAIL")
        if not valid_env:
            print("Environment not set up properly!", end=" ")
        else:
            print("Invalid metadata supplied!", end=" ")
        print("See the exception for more details", flush=True)
        raise

    print("Environment is properly set up!")


def virtual_screen(software: str, *args, **kwargs) -> DockingVirtualScreen:
    return DockingVirtualScreen(get_runner(software), *args, **kwargs)
