"""pyscreener
pythonic interface to virtual screening software
"""
from modulefinder import Module
from . import args
from .docking import build_metadata, check_env, virtual_screen
from .supply import LigandSupply

try:
    from ._version import version

    __version__ = version
except ModuleNotFoundError:
    __version__ = "1.1.1"
