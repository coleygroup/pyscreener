"""pyscreener
pythonic interface to virtual screening software
"""

from pyscreener._version import __version__
from . import args
from .docking import build_metadata, check_env, virtual_screen
from .supply import LigandSupply

# Handle versioneer
from pyscreener._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions