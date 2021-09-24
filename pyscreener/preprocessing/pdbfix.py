"""This module contains the function pdbfix, which is used to fix input PDB
files or retrieve them based on their PDB ID from the PDB"""
from pathlib import Path
from typing import Optional

from pdbfixer import PDBFixer
from openmm.app import PDBFile


def pdbfix(
    receptor: Optional[str] = None,
    pdbid: Optional[str] = None,
    pH: float = 7.0,
    path: str = ".",
) -> str:
    if pdbid is not None:
        fixer = PDBFixer(pdbid=pdbid)
    else:
        fixer = PDBFixer(filename=receptor)

    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH)

    if receptor:
        outfile = receptor
    else:
        outfile = Path(path) / f"{pdbid}.pdb"

    PDBFile.writeFile(fixer.topology, fixer.positions, open(outfile, "w"))

    return outfile


def get_pdb(pdbid: str, pH: float = 7.0, path: str = ".") -> str:
    return pdbfix(pdbid=pdbid, pH=pH, path=path)
