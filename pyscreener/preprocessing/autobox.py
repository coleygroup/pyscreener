from itertools import takewhile
from typing import List, Optional, Tuple

def extract_xyz(line: str) -> Tuple[float, float, float]:
    return tuple(map(float, line.split()[5:8]))

def minimum_bounding_box(coords: List[Tuple[float, float, float]], 
                         buffer: int = 0) -> Tuple[Tuple, Tuple]:
    """
    Calculate the minimum bounding box for a list of coordinates

    Parameters
    ----------
    coords : List[Tuple[float, float, float]]
        a list of x-, y-, and z-coordinates
    buffer : int (Default = 0)
        the amount of buffer to add to the minimum bounding box

    Returns
    -------
    center: Tuple[float, float, float]
        the x-, y-, and z-coordinates of the center of the minimum bounding box
    size: Tuple[float, float, float]
        the x-, y-, and z-radii of the minimum bounding box
    """
    min_x, min_y, min_z = float('inf'), float('inf'), float('inf')
    max_x, max_y, max_z = float('-inf'), float('-inf'), float('-inf')
    for x, y, z in coords:
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

def autobox(pdbfile: str, residues: Optional[List[int]],
            buffer: int = 10) -> Tuple[Tuple, Tuple]:
    if residues:
        center, size = from_residues(pdbfile, residues)
    else:
        center, size = from_docked_ligand(pdbfile, buffer)

    return center, size

def from_residues(pdbfile: str, residues: List[int]) -> Tuple[Tuple, Tuple]:
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
    residues = set(residues)
    with open(pdbfile) as fid:
        residue_coords = []
        for line in fid:    # advance to the atom information lines
            if 'ATOM' in line:
                break

        for line in fid:
            fields = line.split()
            record, _, atom_type, _, _, res_num, x, y, z, _, _, _ = fields
            if 'ATOM' != record:
                break

            if res_num in residues and atom_type == 'CA':
                residue_coords.append((float(x), float(y), float(z)))
    
    return minimum_bounding_box(residue_coords)

def from_docked_ligand(pdbfile: str, buffer: int = 10) -> Tuple[Tuple, Tuple]:
    """
    Generate a ligand autobox from a PDB file containing a docked ligand

    The ligand autobox is the minimum bounding box of the docked ligand with
    an equal buffer in each dimension

    Parameters
    ----------
    pdbfile : str
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
    with open(pdbfile) as fid:
        for line in fid:
            if 'HETATM' in line:
                break

        ligand_atom_coords = [
            extract_xyz(line)
            for line in takewhile(lambda line: 'HETATM' in line, fid)
        ]

    return minimum_bounding_box(ligand_atom_coords, buffer)
