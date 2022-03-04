"""pyscreener
pythonic interface to virtual screening software
"""
from . import args
from .docking import build_metadata, check_env, virtual_screen
from .supply import LigandSupply

try:
    from . import _version

    __version__ = _version.version
except ModuleNotFoundError:
    __version__ = "1.1.1"
