"""pyscreener
pythonic interface to virtual screening software
"""

from . import args
from .docking import build_metadata, check_env, virtual_screen
from .supply import LigandSupply

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions