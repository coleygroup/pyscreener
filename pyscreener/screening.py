from . import docking
from . import md

def screen(mode, **kwargs):
    if mode == 'docking':
        return docking.dock_ligands(**kwargs)
    if mode == 'md':
        return md.run_simulations(**kwargs)
    
    raise ValueError