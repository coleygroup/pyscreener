"""
pyscreener
pythonic interface to virtual screening software
"""
from typing import Dict, List, Tuple

from pyscreener._version import __version__
# from pyscreener.preprocessing import preprocess
# from pyscreener.postprocessing import postprocess
from . import docking

# Handle versioneer
from pyscreener._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

def build_screener(mode, **kwargs) -> Tuple[Dict, List]:
    if mode == 'docking':
        from pyscreener import docking
        return docking.screener(**kwargs)
    if mode == 'md':
        from pyscreener import md
        raise NotImplementedError
    if mode == 'dft':
        from pyscreener import dft
        raise NotImplementedError
    
    raise ValueError(f'Unrecognized screening mode: "{mode}"')
