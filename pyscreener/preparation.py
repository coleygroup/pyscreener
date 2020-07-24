from . import docking, md, dft

def prepare_receptor(mode, **kwargs):
    if mode == 'docking':
        return docking.prepare_receptor(**kwargs)
    if mode == 'md':
        return md.prepare_receptor(**kwargs)
    if mode == 'dft':
        return None
    
    raise ValueError(f'Unrecognized mode: "{mode}"')

def prepare_ligands(mode, **kwargs):
    if mode == 'docking':
        return docking.prepare_ligands(**kwargs)
    if mode == 'md':
        return md.prepare_ligands(**kwargs)
    if mode == 'dft':
        return dft.prepare_ligands(**kwargs)
    
    raise ValueError(f'Unrecognized mode: "{mode}"')