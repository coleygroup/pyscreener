from typing import List

def preprocess(preprocessing_options: List[str], **kwargs):
    if 'none' in preprocessing_options:
        return kwargs

    if 'autobox' in preprocessing_options:
        from .autobox import autobox
        kwargs['center'], kwargs['size'] = autobox(
            pdbfile=kwargs['receptors'][0], **kwargs)

    if 'pdbfix' in preprocessing_options:
        from .pdbfix import pdbfix
        kwargs['receptors'] = [pdbfix(pdbfile=kwargs['receptors'][0], **kwargs)]
    
    if 'tautomers' in preprocessing_options:
        from .tautomers import tautomers
        kwargs['ligands'] = tautomers(**kwargs)
        
    if 'filter' in preprocessing_options:
        from .filter import filter_ligands
        kwargs['ligands'], kwargs['names'] = filter_ligands(**kwargs)

    return kwargs