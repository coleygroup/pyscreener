from . import docking, md

def prepare_receptor(mode, **kwargs):
    if mode == 'docking':
        return docking.prepare_receptor(**kwargs)
    if mode == 'md':
        return md.prepare_receptor(**kwargs)
    
    raise ValueError(f'"{mode}" is not a recognized mode!')

def prepare_ligands(mode, **kwargs):
    if mode == 'docking':
        return docking.prepare_ligands(**kwargs)
    if mode == 'md':
        return md.prepare_ligands(**kwargs)
    
    raise ValueError(f'"{mode}" is not a recognized mode!')