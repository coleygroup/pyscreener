from typing import TypeVar

from . import docking, md, dft

Input = TypeVar('Input')

def prepare(mode, **kwargs) -> Input:
    if mode == 'docking':
        receptor = docking.prepare_receptor(**kwargs)
        ligands = docking.prepare_ligands(**kwargs)

        inputs = (receptor, ligands)
    elif mode == 'md':
        raise NotImplementedError

        # receptor = md.prepare_receptor(**kwargs)
        # ligands = md.prepare_ligands(**kwargs)

        # inputs = receptor, ligands
    elif mode == 'dft':
        raise NotImplementedError
    else:
        raise ValueError(f'Unrecognized mode: "{mode}"')

    return inputs

# def prepare_receptor(mode, **kwargs):
#     if mode == 'docking':
#         return docking.prepare_receptor(**kwargs)
#     if mode == 'md':
#         return md.prepare_receptor(**kwargs)
#     if mode == 'dft':
#         return None
    
#     raise ValueError(f'Unrecognized mode: "{mode}"')

# def prepare_ligands(mode, **kwargs):
#     if mode == 'docking':
#         return docking.prepare_ligands(**kwargs)
#     if mode == 'md':
#         return md.prepare_ligands(**kwargs)
#     if mode == 'dft':
#         return dft.prepare_ligands(**kwargs)
    
#     raise ValueError(f'Unrecognized mode: "{mode}"')