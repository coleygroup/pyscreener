"""This module contains functions for ligand autoboxing in docking
simulations"""
from enum import Enum
from itertools import chain, takewhile
from typing import List, Optional, Tuple

import numpy as np

class PDBRecord(Enum):
    ATOM_TYPE = slice(1, 5)
    ATOM_NAME = slice(13, 17)
    RES_NUMBER = slice(23, 27)
    X_COORD = slice(31, 39)
    Y_COORD = slice(39, 47)
    Z_COORD = slice(47, 55)

def autobox(
    receptors: Optional[List[str]] = None,
    residues: Optional[List[int]] = None,
    docked_ligand_file: Optional[str] = None,
    buffer: int = 10
) -> Tuple[Tuple, Tuple]:
    if residues:
        center, size = residues(receptors[0], residues)
        print('Autoboxing from residues with', end=' ')
    else:
        # allow user to only specify one receptor file
        docked_ligand_file = docked_ligand_file or receptors[0]
        center, size = docked_ligand(docked_ligand_file, buffer)
        print('Autoboxing from docked ligand with', end=' ')

    s_center = f'({center[0]:0.1f}, {center[1]:0.1f}, {center[2]:0.1f})'
    s_size = f'({size[0]:0.1f}, {size[1]:0.1f}, {size[2]:0.1f})'
    print(f'center={s_center} and size={s_size}')

    return center, size

def residues(pdbfile: str, residues: List[int]) -> Tuple[Tuple, Tuple]:
    """Generate a ligand autobox from a list of protein residues

    The ligand autobox is the minimum bounding box of the alpha carbons of the
    listed residues

    Parameters
    ----------
    pdbfile : str
        a PDB-format file containing the protein of interest
    residues: List[int]
        the residue number corresponding to the residues from which to
        calculate the autobox

    Returns
    -------
    center: Tuple[float, float, float]
        the x-, y-, and z-coordinates of the ligand autobox center
    size: Tuple[float, float, float]
        the x-, y-, and z-radii of the ligand autobox
    """
    # ATOM_RECORD_COLUMNS = slice(1, 5)
    # ATOM_NAME_COLUMNS = slice(13, 17)
    # RES_NUMBER_COLUMNS = slice(23, 27)

    residues = set(residues)
    residue_coords = []

    with open(pdbfile) as fid:
        for line in fid:
            if 'ATOM' == line[PDBRecord.ATOM_TYPE.value]:
                break
        lines = chain([line], fid)
        for line in lines:
            if 'ATOM' != line[PDBRecord.ATOM_TYPE.value]:
                break

            if (
                line[PDBRecord.ATOM_NAME.value] == 'CA'
                and line[PDBRecord.RES_NUMBER.value] in residues
            ):
                residue_coords.append(parse_coordinates(line))
    
    return minimum_bounding_box(residue_coords)

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
    buffer : int (Default = 10)
        the buffer to add around the ligand autobox, in Angstroms

    Returns
    -------
    center: Tuple[float, float, float]
        the x-, y-, and z-coordinates of the ligand autobox center
    size: Tuple[float, float, float]
        the x-, y-, and z-radii of the ligand autobox
    """
    with open(docked_ligand_file) as fid:
        for line in fid:
            if 'HETATM' == line[PDBRecord.ATOM_TYPE.value]:
                break
        fid = chain([line], fid)
        coords = [
            parse_coordinates(line)
            for line in takewhile(lambda line: 'HETATM' == line[PDBRecord.ATOM_TYPE.value], fid)
        ]

    return minimum_bounding_box(np.array(coords), buffer)

def parse_coordinates(line: str) -> Tuple[float, float, float]:
    """Parse the x-, y-, and z-coordinates from an ATOM record in a PDB file"""
    # X_COORD_COLUMNS = slice(31, 39)
    # Y_COORD_COLUMNS = slice(39, 47)
    # Z_COORD_COLUMNS = slice(47, 55)

    x = float(line[PDBRecord.X_COORD.value])
    y = float(line[PDBRecord.Y_COORD.value])
    z = float(line[PDBRecord.Y_COORD.value])

    return x, y, z

def minimum_bounding_box(X: np.ndarray, buffer: float = 10.) -> Tuple[Tuple, Tuple]:
    """Calculate the minimum bounding box for a list of coordinates

    Parameters
    ----------
    X : np.ndarray
        an `n x 3` array of x-, y-, z-coordinates, respectively, of points from which to construct
        a minimum bounding box
    buffer : float (Default = 10.)
        the amount of buffer to add to the minimum bounding box

    Returns
    -------
    center: Tuple[float, float, float]
        the x-, y-, and z-coordinates of the center of the minimum bounding box
    radii: Tuple[float, float, float]
        the x-, y-, and z-radii of the minimum bounding box
    """
    # xs, ys, zs = zip(*X)
    # min_x, max_x = min(xs), max(xs)
    # min_y, max_y = min(ys), max(ys)
    # min_z, max_z = min(zs), max(zs)

    # center_x = (max_x + min_x) / 2
    # center_y = (max_y + min_y) / 2
    # center_z = (max_z + min_z) / 2

    # size_x = max_x - center_x + buffer
    # size_y = max_y - center_y + buffer
    # size_z = max_z - center_z + buffer

    # center = center_x, center_y, center_z
    # size = size_x, size_y, size_z

    # return center, size

    center = (X.max(axis=0) - X.min(axis=0)) / 2
    radii = X.max(0) - center + buffer

    return tuple(center), tuple(radii)