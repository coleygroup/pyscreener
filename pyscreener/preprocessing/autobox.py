from itertools import takewhile
from typing import Tuple

def autobox(pdbfile: str, buffer: int = 10) -> Tuple[Tuple, Tuple]:
    """
    Generate a ligand autobox from a PDB file containing a docked ligand

    The ligand autobox is the minimum bounding box of the docked ligand with
    an equal buffer in each dimension

    Parameters
    ----------
    pdbfile : str
        a PDB format file containing the coordinates of a docked ligand
    buffer : int (Default = 10)
        the buffer to add around the ligand autobox, in Angstroms

    Returns
    -------
    center: Tuple[float, float, float]
        the x-, y-, and z-coordinates of the ligand autobox center
    size: Tuple[float, float, float]
        the x-, y-, and z-radii of the ligand autobox
    """
    with open(pdbfile) as fid:
        for line in fid:
            if 'HETATM' in line:
                break

        ligand_atom_coords = [
            list(map(float, line.split()[5:8]))
            for line in takewhile(lambda line: 'HETATM' in line, fid)
        ]

    min_x, min_y, min_z = float('inf'), float('inf'), float('inf')
    max_x, max_y, max_z = float('-inf'), float('-inf'), float('-inf')
    for x, y, z in ligand_atom_coords:
        min_x = min(x, min_x)
        max_x = max(x, max_x)
        
        min_y = min(y, min_y)
        max_y = max(y, max_y)

        min_z = min(z, min_z)
        max_z = max(z, max_z)
    
    c_x = (max_x + min_x) / 2
    c_y = (max_y + min_y) / 2
    c_z = (max_z + min_z) / 2

    s_x = int(max_x - c_x) + buffer
    s_y = int(max_y - c_y) + buffer
    s_z = int(max_z - c_z) + buffer

    center = c_x, c_y, c_z
    size = s_x, s_y, s_z

    return center, size
