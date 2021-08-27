try:
    from importlib import resources
except ModuleNotFoundError:
    import importlib_resources as resources
from itertools import takewhile
import os
from pathlib import Path
import shutil
import subprocess as sp
import sys
from typing import Dict, Iterable, Optional, Tuple

from pyscreener.exceptions import (
    MissingEnvironmentVariableError, MissingFileError
)
from pyscreener.docking.dock.data import SphereMode

with resources.path('pyscreener.docking.dock', '.') as p_module:
    PREP_REC = p_module / 'scripts' / 'prep_rec.py'
    WRITE_DMS = p_module / 'scripts' / 'write_dms.py'

try:
    DOCK6 = Path(os.environ['DOCK6'])
except KeyError:
    raise MissingEnvironmentVariableError(
        'DOCK6 environment variable not set! '
        'See https://github.com/coleygroup/pyscreener#specifying-an-environment-variable for more information'
    )

SPHGEN = DOCK6 / 'bin' / 'sphgen_cpp'
SPHERE_SELECTOR = DOCK6 / 'bin' / 'sphere_selector'
SHOWBOX = DOCK6 / 'bin' / 'showbox'
GRID = DOCK6 / 'bin' / 'grid'
VDW_DEFN_FILE = DOCK6 / 'parameters' / 'vdw_AMBER_parm99.defn'

for f in (SPHGEN, SPHERE_SELECTOR, SHOWBOX, GRID, VDW_DEFN_FILE):
    if not f.exists():
        raise MissingFileError(
            f'File: "{f}" is missing! DOCK parent folder not configured '
            'properly. See https://github.com/coleygroup/pyscreener#specifying-an-environment-variable for more information'
        )

def prepare_receptor(
    receptor: str, probe_radius: float = 1.4,
    steric_clash_dist: float = 0.0,
    min_radius: float = 1.4, max_radius: float = 4.0,
    sphere_mode: SphereMode = SphereMode.LARGEST,
    center: Optional[Tuple[float, float, float]] = None,
    size: Tuple[float, float, float] = (10., 10., 10.), 
    docked_ligand_file: Optional[str] = None,
    buffer: float = 10., enclose_spheres: bool = True,
    path: str = '.'
) -> Optional[Tuple[str, str]]:
    """Prepare the DOCK input files corresponding to the given receptor

    Parameters
    ----------
    receptor : str
        the filepath of a file containing a receptor. Must be in a file that
        is readable by Chimera
    center : Tuple[float, float, float]
        the x-, y-, and z-coordinates of the center of the docking box
    size : Tuple[float, float, float] (Default = (20, 20, 20))
        the x-, y-, and z-radii of the docking box
    docked_ligand_file : Optional[str] (Default = None)
        the filepath of a file containing the coordinates of a docked ligand
    use_largest : bool (Default = False)
        whether to use the largest cluster of spheres when selecting spheres
    buffer : float (Default = 10.)
        the amount of buffer space to be added around the docked ligand when
        selecting spheres and when constructing the docking box if 
        enclose_spheres is True
    enclose_spheres : bool (Default = True)
        whether to calculate the docking box by enclosing the selected spheres
        or to use an input center and radii

    Returns
    -------
    sph_grid : Optional[Tuple[str, str]]
        A tuple of strings with the first entry being the filepath of the file 
        containing the selected spheres and the second being entry the prefix 
        of all prepared grid files. None if receptor preparation fails at any 
        point
    """
    rec_mol2 = prepare_mol2(receptor, path)
    rec_pdb = prepare_pdb(receptor, path)
    if rec_mol2 is None or rec_pdb is None:
        return None

    rec_dms = prepare_dms(rec_pdb, probe_radius, path)
    if rec_dms is None:
        return None

    rec_sph = prepare_sph(
        rec_dms, steric_clash_dist, min_radius, max_radius, path
    )
    if rec_sph is None:
        return None

    rec_sph = select_spheres(
        rec_sph, sphere_mode, center, size,
        docked_ligand_file, buffer, path
    )
    rec_box = prepare_box(rec_sph, center, size, enclose_spheres, buffer, path)
    if rec_box is None:
        return None

    grid_prefix = prepare_grid(rec_mol2, rec_box, path)
    if grid_prefix is None:
        return None

    return rec_sph, grid_prefix

def prepare_mol2(receptor: str, path: str = '.') -> Optional[str]:
    """Prepare a receptor mol2 file from its input file

    Parameter
    ---------
    receptor : str
        the filename of a file containing the receptor
    path : str, default='.'
        the path under which to place the receptor

    Returns
    -------
    receptor_mol2 : Optional[str]
        the filename of the prepared mol2 file
    """
    p_rec_mol2 = Path(path) / f'{Path(receptor).stem}_withH.mol2'
    argv = [
        'chimera', '--nogui', '--script', f'{PREP_REC} {receptor} {p_rec_mol2}'
    ]

    ret = sp.run(argv, stdout=sp.PIPE, stderr=sp.PIPE)
    try:
        ret.check_returncode()
    except sp.SubprocessError:
        print(f'ERROR: failed to convert receptor: "{receptor}"')
        if ret.stderr:
            print(f'Message: {ret.stderr.decode("utf-8")}', file=sys.stderr)
        return None

    return str(p_rec_mol2)

def prepare_pdb(receptor: str, path: str = '.') -> Optional[str]:
    """Prepare a receptor PDB file for usage in DOCK runs

    Parameter
    ---------
    receptor : str
        the filename of a file containing the receptor

    Returns
    -------
    receptor_mol2 : Optional[str]
        the filename of the prepared pdb file
    """
    rec_pdb = str(Path(path) / f'DOCK_{Path(receptor).stem}.pdb')
    args = ['obabel', receptor, '-opdb', '-O', rec_pdb]

    ret = sp.run(args, stderr=sp.PIPE)
    try:
        ret.check_returncode()
    except sp.SubprocessError:
        print(f'ERROR: failed to convert receptor: "{receptor}"')
        if ret.stderr:
            print(f'Message: {ret.stderr.decode("utf-8")}', file=sys.stderr)
        return None
    
    return rec_pdb
    
def prepare_dms(rec_pdb: str, probe_radius: float = 1.4,
                path: str = '.') -> Optional[str]:
    p_rec_dms = Path(path) / f'{Path(rec_pdb).stem}.dms'
    argv = [
        'chimera', '--nogui', '--script',
        f'{WRITE_DMS} {rec_pdb} {probe_radius} {str(p_rec_dms)}'
    ]

    ret = sp.run(argv, stdout=sp.PIPE)
    try:
        ret.check_returncode()
    except sp.SubprocessError:
        print(f'ERROR: failed to generate surface from "{rec_pdb}"',
              file=sys.stderr)
        if ret.stderr:
            print(f'Message: {ret.stderr.decode("utf-8")}', file=sys.stderr)
        # return None

    return str(p_rec_dms)

def prepare_sph(rec_dms: str, steric_clash_dist: float = 0.0,
                min_radius: float = 1.4, max_radius: float = 4.0,
                path: str = '.') -> Optional[str]:
    sph_file = str(Path(path) / f'{Path(rec_dms).stem}.sph')
    argv = [
        SPHGEN, '-i', rec_dms, '-o', sph_file, 
        '-s', 'R', 'd', 'X', '-l', f'{steric_clash_dist:0.1f}',
        'm', str(min_radius), '-x', f'{max_radius:0.1f}'
    ]

    ret = sp.run(argv, stdout=sp.PIPE)
    try:
        ret.check_returncode()
    except sp.SubprocessError:
        print(f'ERROR: failed to generate spheres for "{rec_dms}"',
              file=sys.stderr)
        if ret.stderr:
            print(f'Message: {ret.stderr.decode("utf-8")}', file=sys.stderr)
        return None
    
    return sph_file

def select_spheres(
    sph_file: str, sphere_mode: SphereMode = SphereMode.LARGEST,
    center: Optional[Tuple[float, float, float]] = None,
    size: Optional[Tuple[float, float, float]] = None,
    docked_ligand_file: Optional[str] = None,
    buffer: float = 10.0, path: str = '.'
) -> str:
    """Select a subset of spheres from the sphere cluster file and write these
    spheres to a new file

    Parameters
    ----------
    sph_file : str
        the sphere cluster file to select from
    sphere_mode : SphereMode, default=SphereMode.LARGEST
        the method by which to select spheres
    center : Optional[Tuple[float, float, float]], default=None
        the x-, y-, and z-coordinates of the box center to use
        if using SphereMode.BOX
    size : Optional[Tuple[float, float, float]], default=None
        the x-, y-, and z-radii of the box to use if using SphereMode.BOX
    docked_ligand_file : Optional[str], default=None
        a MOL2 file containing the coordinates of a docked/bound ligand if
        using SphereMode.LIGAND
    buffer : float, default=10.0
        a buffer to add around the coordinates of the docked ligand
    path : str, default='.'
        the path under which to write the selected sphere file

    Returns
    -------
    selected_sph : str
        the filepath of the resulting file containing the selected spheres
    """
    selected_sph = str(
        Path(path) / f'{Path(sph_file).stem}_selected_spheres.sph'
    )

    if sphere_mode == sphere_mode.LIGAND:
        argv = [SPHERE_SELECTOR, sph_file, docked_ligand_file, buffer]
        sp.run(argv, check=True)
        shutil.move('selected_spheres.sph', selected_sph)
        return selected_sph

    with open(sph_file, 'r') as fid_in, open(selected_sph, 'w') as fid_out:
        if sphere_mode == sphere_mode.LARGEST:
            clusterline = lambda line: 'cluster' not in line

            for _ in takewhile(clusterline, fid_in):
                continue
            lines = [line for line in takewhile(clusterline, fid_in)]
            fid_out.write(f'DOCK spheres largest cluster\n')

        elif sphere_mode == sphere_mode.BOX:
            lines = [
                line for line in fid_in if inside_box(line, center, size)
            ]
            fid_out.write(f'DOCK spheres within radii {size} of {center}\n')

        fid_out.write(
            f'cluster     1 number of spheres in cluster {len(lines)}\n'
        )
        fid_out.writelines(lines)

    return selected_sph

def inside_box(line: str, center: Tuple[float, float, float],
               size: Tuple[float, float, float]) -> bool:
    """Are the coordinates of sphere corresponding to a line in a sphere
    cluster file are inside the given box?"""
    try:
        xyz = [float(token) for token in line.split()[1:4]]
    except ValueError:
        return False

    x, y, z = xyz
    bx, by, bz = center
    ra, rb, rc = size

    return all([
        bx - ra <= x <= bx + ra,
        by - rb <= y <= by + rb,
        bz - rc <= z <= bz + rc,
    ])
    
def prepare_box(
    sph_file: str, center: Tuple[float, float, float],
    size: Tuple[float, float, float], enclose_spheres: bool = True,
    buffer: float = 10.0, path: str = '.'
) -> Optional[str]:
    shutil.copyfile(sph_file, 'tmp_spheres.sph')
    box_file = str(Path(path) / f'{Path(sph_file).stem}_box.pdb')

    if enclose_spheres:
        showbox_input = f'Y\n{buffer}\ntmp_spheres.sph\n1\n'
    else:
        x, y, z = center
        a, b, c = size
        showbox_input = f'N\nU\n{x} {y} {z}\n{a} {b} {c}\n'

    showbox_input += 'tmp_box.pdb\n'

    # with open('box.in', 'w') as fid:
    #     if enclose_spheres:
    #         fid.write('Y\n')
    #         fid.write(f'{buffer}\n')
    #         fid.write(f'{sph_file}\n')
    #         fid.write(f'1\n')
    #     else:
    #         fid.write('N\n')
    #         fid.write('U\n')
    #         fid.write(f'[{center[0]} {center[1]} {center[2]}]\n')
    #         fid.write(f'[{size[0]} {size[1]} {size[2]}]\n')
    #     fid.write(f'{box_file}\n')
    
    ret = sp.run(
        [SHOWBOX], input=showbox_input,
        universal_newlines=True, stdout=sp.PIPE
    )
    try:
        ret.check_returncode()
    except sp.SubprocessError:
        print(f'ERROR: failed to generate box corresponding to "{sph_file}"',
              file=sys.stderr)
        if ret.stderr:
            print(f'Message: {ret.stderr.decode("utf-8")}', file=sys.stderr)
        return None

    os.unlink('tmp_spheres.sph')
    shutil.move('tmp_box.pdb', box_file)
    return box_file

def prepare_grid(rec_mol2: str, box_file: str,
                 path: str = '.') -> Optional[str]:
    p_grid_stem = Path(path) / f'{Path(rec_mol2).stem}_grid'

    shutil.copy(box_file, 'tmp_box.pdb')
    with open('grid.in', 'w') as fid:
        fid.write('compute_grids yes\n')
        fid.write('grid_spacing 0.4\n')
        fid.write('output_molecule no\n')
        fid.write('contact_score no\n')
        fid.write('energy_score yes\n')
        fid.write('energy_cutoff_distance 9999\n')
        fid.write('atom_model a\n')
        fid.write('attractive_exponent 6\n')
        fid.write('repulsive_exponent 12\n')
        fid.write('distance_dielectric yes\n')
        fid.write('dielectric_factor 4.0\n')
        fid.write('bump_filter yes\n')
        fid.write('bump_overlap 0.75\n')
        fid.write(f'receptor_file {rec_mol2}\n')
        fid.write('box_file tmp_box.pdb\n')
        fid.write(f'vdw_definition_file {VDW_DEFN_FILE}\n')
        fid.write(f'score_grid_prefix  {p_grid_stem}\n')
    
    ret = sp.run([GRID, '-i', 'grid.in', '-o', 'gridinfo.out'], stdout=sp.PIPE)
    try:
        ret.check_returncode()
    except sp.SubprocessError:
        print(f'ERROR: failed to generate grid from {rec_mol2}',
              file=sys.stderr)
        if ret.stderr:
            print(f'Message: {ret.stderr.decode("utf-8")}', file=sys.stderr)
        return None

    os.unlink('tmp_box.pdb')
    return str(p_grid_stem)