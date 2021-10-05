"""pyscreener
pythonic interface to virtual screening software
"""

from pyscreener._version import __version__
from .docking import build_metadata, virtual_screen
from .postprocessing import postprocess
from .supply import LigandSupply

# Handle versioneer
from pyscreener._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions