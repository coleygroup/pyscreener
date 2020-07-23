from . import autodock
from . import ucsfdock
from . import rdock

def prepare_receptor(docker: str, receptor, **kwargs):
    if docker in {'vina', 'smina', 'psovina', 'qvina'}:
        return autodock.prepare_receptor(receptor)
    if docker == 'dock':
        return ucsfdock.prepare_receptor(receptor, **kwargs)
    if docker == 'rdock':
        return rdock.prepare_receptor(receptor, **kwargs)

    raise ValueError(f'{docker} is not a recognized docking program')

def prepare_ligands(docker: str, ligands, **kwargs):
    if docker in {'vina', 'smina', 'psovina', 'qvina'}:
        return autodock.prepare_ligands(ligands, **kwargs)
    if docker == 'dock':
        return ucsfdock.prepare_ligands(ligands, **kwargs)
    if docker == 'rdock':
        return rdock.prepare_ligands(ligands, **kwargs)

    raise ValueError(f'{docker} is not a recognized docking program')