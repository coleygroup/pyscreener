"""pyscreener
pythonic interface to virtual screening software
"""
from importlib.metadata import PackageNotFoundError, version

from . import args
from .docking import build_metadata, check_env, virtual_screen
from .supply import LigandSupply

try:
    __version__ = version("pyscreener")
except PackageNotFoundError:
    __version__ = "1.2.0"
