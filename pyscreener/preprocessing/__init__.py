from typing import List

from .autobox import autobox
from .filter import filter_ligands

def preprocess(preprocessing_options: List[str], **kwargs):
    if 'none' in preprocessing_options:
        return kwargs

    if 'autobox' in preprocessing_options:
        kwargs['center'], kwargs['size'] = autobox(
            pdbfile=kwargs['receptors'][0], **kwargs)

    if 'pdbfix' in preprocessing_options:
        from .pdbfix import pdbfix
        kwargs['receptors'] = [pdbfix(pdbfile=kwargs['receptors'][0], **kwargs)]
    
    if 'filter' in preprocessing_options:
        kwargs['ligands'], kwargs['names'] = filter_ligands(**kwargs)
    
    if 'tautomers' in preprocessing_options:
        pass

    return kwargs