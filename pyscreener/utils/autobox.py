"""This module contains functions for ligand autoboxing in docking simulations"""
from enum import Enum
from itertools import takewhile
from pyscreener.exceptions import BadPDBFileError
from typing import Iterable, List, Optional, Tuple

import numpy as np


class PDBRecord(Enum):
    NAME = slice(0, 6)
    ATOM = slice(12, 16)
    RES_SEQ = slice(22, 26)
    X_COORD = slice(30, 38)
    Y_COORD = slice(38, 46)
    Z_COORD = slice(46, 54)


def autobox(
    receptors: Optional[List[str]] = None,
    residues: Optional[List[int]] = None,
    docked_ligand_file: Optional[str] = None,
    buffer: int = 10,
) -> Tuple[Tuple, Tuple]:
    if residues:
        center, size = residues(receptors[0], residues, buffer)
        print("Autoboxing from residues with", end=" ")
    else:
        # allow user to only specify one receptor file
        docked_ligand_file = docked_ligand_file or receptors[0]
        center, size = docked_ligand(docked_ligand_file, buffer)
        print("Autoboxing from docked ligand with", end=" ")

    s_center = f"({center[0]:0.1f}, {center[1]:0.1f}, {center[2]:0.1f})"
    s_size = f"({size[0]:0.1f}, {size[1]:0.1f}, {size[2]:0.1f})"
    print(f"center={s_center} and size={s_size}")

    return center, size


def residues(pdbfile: str, residues: List[int], buffer: float = 0.0) -> Tuple[Tuple, Tuple]:
    """Generate a ligand autobox from a list of protein residues

    The ligand autobox is the minimum bounding box of the alpha carbons of the
    listed residues

    Parameters
    ----------
    pdbfile : str
        a PDB-format file containing the protein of interest
    residues: List[int]
        the residue numbers corresponding to the residues from which to calculate the autobox
    buffer : int, default=0.
        the amount of buffer to add around the ligand autobox, in Angstroms

    Returns
    -------
    center: Tuple[float, float, float]
        the x-, y-, and z-coordinates of the ligand autobox center
    size: Tuple[float, float, float]
        the x-, y-, and z-radii of the ligand autobox
    """
    lines = extract_residues_lines(pdbfile, residues)
    coords = [parse_coordinates(line) for line in lines]

    if len(coords) == 0:
        raise ValueError(
            f'PDB file "{pdbfile} does not contain any of the residue numbers: {residues}!"'
        )

    return minimum_bounding_box(np.array(coords), buffer)


def extract_residues_lines(pdb_filepath, residues: Iterable[int]):
    """Extract the CA lines for the numbered residues in the input PDB file"""
    residues = set(residues)
    if len(residues) == 0:
        raise ValueError("No residues were specified!")

    lines = []
    with open(pdb_filepath) as fid:
        for line in fid:
            if (
                line[PDBRecord.NAME.value].strip() == "ATOM"
                and line[PDBRecord.ATOM.value].strip() == "CA"
                and int(line[PDBRecord.RES_SEQ.value].strip()) in residues
            ):
                lines.append(line)

    return lines


def docked_ligand(docked_ligand_file: str, buffer: int = 10) -> Tuple[Tuple, Tuple]:
    """Generate a ligand autobox from a PDB file containing a docked ligand

    The ligand autobox is the minimum bounding box of the docked ligand with
    an equal buffer in each dimension. The PDB file should be either just the
    protein with a single docked ligand or a PDB file containing just the
    docked ligand.

    Parameters
    ----------
    docked_ligand_file : str
        a PDB-format file containing the coordinates of a docked ligand
    buffer : int, default=10.
        the amount of buffer to add around the ligand autobox, in Angstroms

    Returns
    -------
    center: Tuple[float, float, float]
        the x-, y-, and z-coordinates of the ligand autobox center
    size: Tuple[float, float, float]
        the x-, y-, and z-radii of the ligand autobox
    """
    hetatm_lines = extract_hetatm_lines(docked_ligand_file)
    coords = [parse_coordinates(line) for line in hetatm_lines]

    if len(coords) == 0:
        raise BadPDBFileError(f'No HETATM coordinates could be parsed from "{docked_ligand_file}"')

    return minimum_bounding_box(np.array(coords), buffer)


def extract_hetatm_lines(pdb_filepath: str):
    """extract the lines of the first HETAM in the PDB file"""
    lines = []
    with open(pdb_filepath) as fid:
        for line in fid:
            if "HETATM" == line[PDBRecord.NAME.value].strip():
                lines.append(line)
                break
        lines.extend(takewhile(lambda line: "HETATM" == line[PDBRecord.NAME.value].strip(), fid))

    return lines


def parse_coordinates(line: str) -> Tuple[float, float, float]:
    """Parse the x-, y-, and z-coordinates from an ATOM/HETATM record in a PDB file"""
    try:
        x = float(line[PDBRecord.X_COORD.value])
        y = float(line[PDBRecord.Y_COORD.value])
        z = float(line[PDBRecord.Z_COORD.value])
    except ValueError:
        raise BadPDBFileError(f'could not parse line: "{line}"')

    return x, y, z


def minimum_bounding_box(X: np.ndarray, buffer: float = 10.0) -> Tuple[Tuple, Tuple]:
    """Calculate the minimum bounding box for a list of coordinates

    Parameters
    ----------
    X : np.ndarray
        an `n x 3` array of x-, y-, z-coordinates, respectively, of points from which to construct
        a minimum bounding box
    buffer : float, default=10.
        the amount of buffer to add to the minimum bounding box

    Returns
    -------
    center: Tuple[float, float, float]
        the x-, y-, and z-coordinates of the center of the minimum bounding box
    radii: Tuple[float, float, float]
        the x-, y-, and z-radii of the minimum bounding box
    """
    center = (X.max(axis=0) + X.min(axis=0)) / 2
    radii = X.max(axis=0) - center + buffer

    return tuple(center), tuple(radii)
