from typing import Optional

from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

def pdbfix(pdbfile, pdbid: Optional[str] = None, pH: float = 7.0) -> str:
    if pdbid:
        fixer = PDBFixer(pdbid=pdbid)
    else:
        fixer = PDBFixer(filename=pdbfile)

    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH)

    PDBFile.writeFile(fixer.topology, fixer.positions, open(pdbfile, 'w'))
    
    return pdbfile