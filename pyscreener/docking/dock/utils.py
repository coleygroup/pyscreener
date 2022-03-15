from enum import auto
from itertools import takewhile
from importlib import resources
import os
from pathlib import Path
import shutil
import subprocess as sp
from typing import Mapping, Optional, Tuple

from pyscreener.utils import AutoName
from pyscreener.exceptions import (
    MissingEnvironmentVariableError,
    MisconfiguredDirectoryError,
    ReceptorPreparationError,
)
from pyscreener.docking.dock.exceptions import (
    BoxGenerationError,
    GridGenerationError,
    SphereGenerationError,
    SurfaceGenerationError,
)

with resources.path("pyscreener.docking.dock", ".") as p_module:
    PREP_REC = p_module / "scripts" / "prep_rec.py2"
    WRITE_DMS = p_module / "scripts" / "write_dms.py2"

try:
    DOCK6 = Path(os.environ["DOCK6"])
except KeyError:
    raise MissingEnvironmentVariableError(
        "DOCK6 environment variable not set! "
        "See https://github.com/coleygroup/pyscreener#specifying-an-environment-variable for more information"
    )

SPHGEN = DOCK6 / "bin" / "sphgen_cpp"
SPHERE_SELECTOR = DOCK6 / "bin" / "sphere_selector"
SHOWBOX = DOCK6 / "bin" / "showbox"
GRID = DOCK6 / "bin" / "grid"
VDW_DEFN_FILE = DOCK6 / "parameters" / "vdw_AMBER_parm99.defn"

for f in (SPHGEN, SPHERE_SELECTOR, SHOWBOX, GRID, VDW_DEFN_FILE):
    if not f.exists():
        raise MisconfiguredDirectoryError(
            "DOCK6 directory not configured properly! "
            f'DOCK6 path is set as "{DOCK6}", but there is no '
            f'"{f.name}" located in the "{f.parents[0].name}" subdirectory '
            "under the DOCK6 path. See https://github.com/coleygroup/pyscreener#specifying-an-environment-variable for more information"
        )


class SphereMode(AutoName):
    BOX = auto()
    LARGEST = auto()
    LIGAND = auto()


def prepare_mol2(receptor: str, path: str = ".") -> str:
    """Prepare a receptor mol2 file from its input file

    Parameter
    ---------
    receptor : str
        the filename of a file containing the receptor
    path : str, default="."
        the path under which to place the receptor

    Returns
    -------
    str
        the filepath of the prepared mol2 file
    """
    p_rec_mol2 = Path(path) / f"{Path(receptor).stem}_withH.mol2"
    argv = ["chimera", "--nogui", "--script", f"{PREP_REC} {receptor} {p_rec_mol2}"]

    ret = sp.run(argv, stdout=sp.PIPE, stderr=sp.PIPE)
    try:
        ret.check_returncode()
    except sp.SubprocessError:
        raise ReceptorPreparationError(
            f'failed to prepare MOL2 "{receptor}".' + f'Message: {ret.stderr.decode("utf-8")}'
            if ret.stderr
            else ""
        )

    return str(p_rec_mol2)


def prepare_pdb(receptor: str, path: str = ".") -> str:
    """Prepare a receptor PDB file for usage in DOCK runs

    Parameter
    ---------
    receptor : str
        the filename of a file containing the receptor

    Returns
    -------
    str
        the filepath of the prepared pdb file
    """
    rec_pdb = str(Path(path) / f"DOCK_{Path(receptor).stem}.pdb")
    args = ["obabel", receptor, "-opdb", "-O", rec_pdb]

    ret = sp.run(args, stderr=sp.PIPE)
    try:
        ret.check_returncode()
    except sp.SubprocessError:
        raise ReceptorPreparationError(
            f'failed to prepare PDB for "{receptor}".'
            + +(f'Message: {ret.stderr.decode("utf-8")}' if ret.stderr else "")
        )

    return rec_pdb


def prepare_dms(rec_pdb: str, probe_radius: float = 1.4, path: str = ".") -> str:
    p_rec_dms = Path(path) / f"{Path(rec_pdb).stem}.dms"
    argv = [
        "chimera",
        "--nogui",
        "--script",
        f"{WRITE_DMS} {rec_pdb} {probe_radius} {str(p_rec_dms)}",
    ]

    ret = sp.run(argv, stdout=sp.PIPE)
    try:
        ret.check_returncode()
    except sp.SubprocessError:
        raise SurfaceGenerationError(
            f'failed to generate surface for "{rec_pdb}".'
            + (f'Message: {ret.stderr.decode("utf-8")}' if ret.stderr else "")
        )

    return str(p_rec_dms)


def prepare_sph(
    rec_dms: str,
    steric_clash_dist: float = 0.0,
    min_radius: float = 1.4,
    max_radius: float = 4.0,
    path: str = ".",
) -> Optional[str]:
    sph_file = str(Path(path) / f"{Path(rec_dms).stem}.sph")
    argv = [
        SPHGEN,
        "-i",
        rec_dms,
        "-o",
        sph_file,
        "-s",
        "R",
        "d",
        "X",
        "-l",
        f"{steric_clash_dist:0.1f}",
        "m",
        str(min_radius),
        "-x",
        f"{max_radius:0.1f}",
    ]

    ret = sp.run(argv, stdout=sp.PIPE)
    try:
        ret.check_returncode()
    except sp.SubprocessError:
        raise SphereGenerationError(
            f'failed to generate spheres for "{rec_dms}".'
            + (f'Message: {ret.stderr.decode("utf-8")}' if ret.stderr else "")
        )

    return sph_file


def select_spheres(
    sph_file: str,
    sphere_mode: SphereMode = SphereMode.BOX,
    center: Optional[Tuple[float, float, float]] = None,
    size: Optional[Tuple[float, float, float]] = None,
    docked_ligand_file: Optional[str] = None,
    buffer: float = 10.0,
    path: str = ".",
) -> str:
    """Select a subset of spheres from the sphere cluster file and write these
    spheres to a new file

    Parameters
    ----------
    sph_file : str
        the sphere cluster file to select from
    sphere_mode : SphereMode, default=SphereMode.BOX
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
    selected_sph = str(Path(path) / f"{Path(sph_file).stem}_selected_spheres.sph")

    if sphere_mode == sphere_mode.LIGAND:
        argv = [SPHERE_SELECTOR, sph_file, docked_ligand_file, buffer]
        sp.run(argv, check=True)
        shutil.move("selected_spheres.sph", selected_sph)
        return selected_sph

    with open(sph_file, "r") as fid_in, open(selected_sph, "w") as fid_out:
        if sphere_mode == sphere_mode.LARGEST:
            clusterline = lambda line: "cluster" not in line  # noqa: E731

            for _ in takewhile(clusterline, fid_in):
                continue
            lines = [line for line in takewhile(clusterline, fid_in)]
            fid_out.write("DOCK spheres largest cluster\n")

        elif sphere_mode == sphere_mode.BOX:
            lines = [line for line in fid_in if inside_box(line, center, size)]
            fid_out.write(f"DOCK spheres within radii {size} of {center}\n")

        fid_out.write(f"cluster     1 number of spheres in cluster {len(lines)}\n")
        fid_out.writelines(lines)

    return selected_sph


def inside_box(
    line: str, center: Tuple[float, float, float], size: Tuple[float, float, float]
) -> bool:
    """Are the coordinates of sphere corresponding to a line in a sphere
    cluster file are inside the given box?"""
    try:
        xyz = [float(token) for token in line.split()[1:4]]
    except ValueError:
        return False

    x, y, z = xyz
    bx, by, bz = center
    ra, rb, rc = size

    return all([bx - ra <= x <= bx + ra, by - rb <= y <= by + rb, bz - rc <= z <= bz + rc])


def prepare_box(
    sph_file: str,
    center: Tuple[float, float, float],
    size: Tuple[float, float, float],
    enclose_spheres: bool = True,
    buffer: float = 10.0,
    path: str = ".",
) -> Optional[str]:
    # long file paths causes issues with showbox, so use temp files with shorter names
    TMP_SPH_FILE = "tmp_spheres.sph"
    TMP_BOX_FILE = "tmp_box.pdb"

    shutil.copyfile(sph_file, TMP_SPH_FILE)
    box_file = str(Path(path) / f"{Path(sph_file).stem}_box.pdb")

    if enclose_spheres:
        showbox_input = f"Y\n{buffer}\n{TMP_SPH_FILE}\n1\n"
    else:
        x, y, z = center
        a, b, c = size
        showbox_input = f"N\nU\n{x} {y} {z}\n{a} {b} {c}\n"

    showbox_input += f"{TMP_BOX_FILE}\n"

    ret = sp.run([SHOWBOX], input=showbox_input, universal_newlines=True, stdout=sp.PIPE)
    try:
        ret.check_returncode()
    except sp.SubprocessError:
        raise BoxGenerationError(
            f'failed to generate generate box corresponding to "{sph_file}".'
            + (f'Message: {ret.stderr.decode("utf-8")}' if ret.stderr else "")
        )

    os.unlink(TMP_SPH_FILE)
    shutil.move(TMP_BOX_FILE, box_file)
    return box_file


def prepare_grid(
    rec_mol2: str, box_file: str, path: str = ".", params: Optional[Mapping] = None
) -> str:
    # grid has has problems with long file lengths, so use temp file to shorten
    TMP_BOX_FILE = "tmp_box.pdb"
    shutil.copy(box_file, TMP_BOX_FILE)

    DEFAULT_PARAMS = {
        "compute_grids": "yes",
        "grid_spacing": 0.4,
        "output_molecule": "no",
        "contact_score": "no",
        "energy_score": "yes",
        "energy_cutoff_distance": 9999,
        "atom_model": "a",
        "attractive_exponent": 6,
        "repulsive_exponent": 12,
        "distance_dielectric": "yes",
        "dielectric_factor": 4.0,
        "bump_filter": "yes",
        "bump_overlap": 0.75,
    }
    if params is None:
        params = DEFAULT_PARAMS
    else:
        params = {**params, **DEFAULT_PARAMS}
    p_grid_stem = Path(path) / f"{Path(rec_mol2).stem}_grid"

    with open("grid.in", "w") as fid:
        [print(f"{k} {v}", file=fid) for k, v in params.items()]
        print(f"receptor_file {rec_mol2}", file=fid)
        print("box_file tmp_box.pdb", file=fid)
        print(f"vdw_definition_file {VDW_DEFN_FILE}", file=fid)
        print(f"score_grid_prefix  {p_grid_stem}", file=fid)

    ret = sp.run([GRID, "-i", "grid.in", "-o", "gridinfo.out"], stdout=sp.PIPE)
    try:
        ret.check_returncode()
    except sp.SubprocessError:
        raise GridGenerationError(
            f'failed to generate generate grid from "{rec_mol2}".'
            + (f'Message: {ret.stderr.decode("utf-8")}' if ret.stderr else "")
        )

    os.unlink(TMP_BOX_FILE)
    return str(p_grid_stem)
