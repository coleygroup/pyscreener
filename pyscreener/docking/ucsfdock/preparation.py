from itertools import takewhile
from pathlib import Path
import subprocess as sp
import sys
from typing import Dict, Iterable, Optional, Tuple

from ..utils import OBABEL, Ligand
from ..preparation import prepare_receptors, prepare_ligands
from ...utils import Input

VDW_DEFN_FILE = 'dock6/parameters/vdw_AMBER_parm99.defn'
FLEX_DEFN_FILE = 'dock6/parameters/flex.defn'
FLEX_DRIVE_FILE = 'dock6/parameters/flex_drive.tbl'

def prepare_inputs(docker: str, receptors: Iterable[str], ligands: Iterable,
                   center: Tuple, size: Tuple[int, int, int] = (20, 20, 20), 
                   ncpu: int = 1, path: str = '.',
                   **kwargs) -> Dict:
    # 1 preparation of ligands and receptors
    receptors = prepare_receptors(receptors, prepare_receptor)
    ligands = prepare_ligands(ligands, prepare_from_smi,
                              prepare_from_file, **kwargs)
    # 2 generating receptors surfaces and spheres
    # 3 generating boxes and grids
    # 4 generating ensemble input files
    
    # ligands is type List[Tuple[str, List[str]]] where the innermost
    # List[str] is the ensemble_infiles
    return {'ligands': ligands}

def prepare_receptor(receptor: str):
    """Prepare a receptor mol2 file from its input file

    Parameter
    ---------
    receptor : str
        the filename of a file containing the receptor

    Returns
    -------
    receptor_mol2 : str
        the filename of the resulting MOL2 file
    """
    p_rec = Path(receptor)
    rec_withH = str(p_rec.with_name(f'{p_rec.stem}_withH.mol2'))
    rec_noH = str(p_rec.with_name(f'{p_rec.stem}_noH.pdb'))

    # receptor_file = str(Path(receptor).with_suffix('.mol2'))
    args_withH = [OBABEL, receptor, '-omol2', '-O', rec_withH,
                  '-h', '--partialcharge', 'gasteiger']
    args_noH = [OBABEL, receptor, '-opdb', '-O', rec_noH]
    try:
        sp.run(args_withH, stderr=sp.PIPE, check=True)
        sp.run(args_noH, stderr=sp.PIPE, check=True)
    except sp.SubprocessError:
        print(f'ERROR: failed to convert {receptor}, skipping...')

    return rec_withH, rec_noH

def prepare_from_smi(smi: str, name: str = 'ligand',
                     path: str = '.', **kwargs) -> Optional[Ligand]:
    """Prepare an input ligand file from the ligand's SMILES string

    Parameters
    ----------
    smi : str
        the SMILES string of the ligand
    name : Optional[str] (Default = None)
        the name of the ligand.
    path : str (Default = '.')
        the path under which the output PDBQT file should be written
    **kwargs
        additional and unused keyword arguments

    Returns
    -------
    Optional[Ligand]
        a tuple of the SMILES string and the corresponding prepared input file.
        None if preparation failed for any reason
    """
    path = Path(path)
    if not path.is_dir():
        path.mkdir()
    
    mol2 = str(path / f'{name}.mol2')

    argv = [OBABEL, f'-:{smi}', '-omol2', '-O', mol2,
            '-h', '--gen3d', '--partialcharge', 'gasteiger']
    ret = sp.run(argv, check=False, stderr=sp.PIPE)

    try:
        ret.check_returncode()
        return smi, mol2
    except sp.SubprocessError:
        return None

def prepare_from_file(filename: str, use_3d: bool = False,
                      name: Optional[str] = None, path: str = '.', 
                      **kwargs) -> Ligand:
    """Convert a single ligand to the appropriate input format

    Parameters
    ----------
    filename : str
        the name of the file containing the ligand
    use_3d : bool (Default = False)
        whether to use the 3D information in the input file (if possible)
    prepare_from_smi: Callable[..., Tuple[str, str]]
        a function that prepares an input ligand file from a SMILES string
    name : Optional[str] (Default = None)
        the name of the ligand. If None, use the stem of the input file
    path : str (Default = '.')
        the path under which the output .pdbqt file should be written
    **kwargs
        additional and unused keyword arguments

    Returns
    -------
    List[Ligand]
        a tuple of the SMILES string the prepared input file corresponding
        to the molecule contained in filename
    """
    name = name or Path(filename).stem

    ret = sp.run([OBABEL, filename, '-osmi'], stdout=sp.PIPE, check=True)
    lines = ret.stdout.decode('utf-8').splitlines()
    smis = [line.split()[0] for line in lines]

    if not use_3d:
        ligands = [prepare_from_smi(smi, f'{name}_{i}', path) 
                   for i, smi in enumerate(smis)]
        return [lig for lig in ligands if lig]
    
    path = Path(path)
    if not path.is_dir():
        path.mkdir()

    mol2 = f'{path}/{name}_.mol2'
    argv = [OBABEL, filename, '-omol2', '-O', mol2, '-m']
    ret = sp.run(argv, check=False, stderr=sp.PIPE)
    
    try:
        ret.check_returncode()
    except sp.SubprocessError:
        return None

    stderr = ret.stderr.decode('utf-8')
    for line in stderr.splitlines():
        if 'converted' not in line:
            continue
        n_mols = int(line.split()[0])

    mol2s = [f'{path}/{name}_{i}.mol2' for i in range(1, n_mols)]

    return list(zip(smis, mol2s))

def prepare_dms(rec_noH_pdb: str, probe_radius: float = 1.4):
    rec_dms = str(Path(rec_noH).with_suffix('.dms'))
    argv = ['dms', rec_noH_pdb, '-n', '-w', probe_radius, '-v', '-o', rec_dms]
    sp.run(argv, check=True)

    return rec_dms

# def prepare_sph(dms_file: str, steric_clash_dist: float = 0.0,
#                 max_radius: float = 4.0, min_radius: float = 1.4,
#                 center: Tuple[float, float, float],
#                 size: Tuple[float, float, float],
#                 docked_ligand_file: Optional[str] = None,
#                 use_largest: bool = False, buffer: float = 10.0) -> str:
#     pass

def gen_sph(dms_file: str, steric_clash_dist: float = 0.0,
            max_radius: float = 4.0, min_radius: float = 1.4) -> str:
    sph_file = str(Path(dms_file).with_suffix('.sph'))

    with open('INSPH', 'w') as fid:
        fid.write(f'{dms_file}\n')
        fid.write('R\n')
        fid.write('X\n')
        fid.write(f'{steric_clash_dist:0.1f}\n')
        fid.write(f'{max_radius:0.1f}\n')
        fid.write(f'{min_radius:0.1f}\n')
        fid.write(f'{sph_file}\n')
    
    argv = ['sphgen', '-i', 'INSPH', '-o', 'OUTSPH']
    sp.run(argv, check=True)
    
    # try:
    #     ret.check_returncode()
    # except sp.SubprocessError:
    #     return None

    return sph_file

def select_spheres(sph_file: str, 
                   center: Tuple[float, float, float],
                   size: Tuple[float, float, float],
                   docked_ligand_file: Optional[str] = None,
                   use_largest: bool = False, buffer: float = 10.0) -> str:
    if docked_ligand_file:
        argv = ['sphere_selector', sph_file, docked_ligand_file, buffer]
        sp.run(argv, check=True)
        # sphere_selector always outputs this filename
        return 'selected_spheres.sph'

    p_sph = Path(sph_file)
    selected_sph_file = str(p_sph.with_name(f'{p_sph.stem}_selected'))
    
    def parse_xyz(line):
        tokens = line.split(' ')
        x, y, z = tokens[1:4]
        return float(x), float(y), float(z)

    def inside_docking_box(x, y, z):
        for i, coord in enumerate([x, y, z]):
            if not center[i]-size[i] <= coord <= center[i]+size[i]:
                return False
        return True

    with open(sph_file, 'r') as fid_in, open(selected_sph_file, 'w') as fid_out:
        fid_out.write('DOCK receptor spheres\n')
        fid_out.write('Pyscreener selected spheres\n')
        
        if use_largest:
            # consume all intro lines
            for line in takewhile(lambda line: 'cluster' not in line, fid_in):
                continue

            # write all lines in the first (largest) cluster
            for line in takewhile(lambda line: 'cluster' not in line, fid_in):
                fid_out.write(line)
        else:
            for line in fid_in:
                try:
                    x, y, z = parse_xyz(line)
                except ValueError:
                    continue    # skip unparsable lines

                if inside_docking_box(x, y, z):
                    fid_out.write(line)

    return selected_sph_file

def prepare_box(sph_file: str, enclose_spheres: bool, buffer: float,
                center: Tuple[float, float, float],
                size: Tuple[float, float, float]) -> str:
    p_sph = Path(sph_file)
    p_box = p_sph.with_name(f'{p_sph.stem}_box.pdb')
    box_file = str(p_box)

    with open('box.in', 'r') as fid:
        if enclose_spheres:
            fid.write('Y\n')
            fid.write(f'{buffer}\n')
            fid.write(f'{sph_file}\n')
            fid.write(f'1\n')
        else:
            fid.write('N\n')
            fid.write('U\n')
            fid.write(f'[{center[0]} {center[1]} {center[2]}]\n')
            fid.write(f'[{size[0]} {size[1]} {size[2]}]\n')
        fid.write(f'{box_file}\n')
    
    sp.run('showbox < box.in', shell=True, check=True)

    return box_file

def prepare_grid(rec_withH: str, box_file: str):
    p_rec = Path(rec_withH)
    grid_prefix = f'{p_rec.stem}_grid'  # probably need to include path prefix

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
        fid.write(f'receptor_file {receptor_file}\n')
        fid.write(f'box_file {box_file}\n')
        fid.write(f'vdw_definition_file {VDW_DEFN_FILE}\n')
        fid.write(f'score_grid_prefix  {grid_prefix}\n')

    sp.run(['grid', '-i', 'grid.in', '-o', 'gridinfo.out'], check=True)

    return grid_prefix

def prepare_input_file(ligand_file: str, sph_file, grid_prefix) -> str:
    input_file = f'{Path(ligand_file).stem}_{Path(sph_file).stem}.in'
    input_file = Path(path) / input_file

    with open(input_file, 'w') as fid:
        fid.write('conformer_search_type flex\n')
        fid.write('write_fragment_libraries no\n')
        fid.write('user_specified_anchor no\n')
        fid.write('limit_max_anchors no\n')
        fid.write('min_anchor_size 5\n')

        fid.write('pruning_use_clustering yes\n')
        fid.write('pruning_max_orients 100\n')
        fid.write('pruning_clustering_cutoff 100\n')
        fid.write('pruning_conformer_score_cutoff 100.0\n')
        fid.write('pruning_conformer_score_scaling_factor 1.0\n')

        fid.write('use_clash_overlap no\n')
        fid.write('write_growth_tree no\n')
        fid.write('use_internal_energy yes\n')
        fid.write('internal_energy_rep_exp 12\n')
        fid.write('internal_energy_cutoff 100.0\n')

        fid.write(f'ligand_atom_file {ligand_file}\n')
        fid.write('limit_max_ligands no\n')
        fid.write('skip_molecule no\n')
        fid.write('read_mol_solvation no\n')
        fid.write('calculate_rmsd no\n')
        fid.write('use_database_filter no\n')
        fid.write('orient_ligand yes\n')
        fid.write('automated_matching yes\n')
        fid.write(f'receptor_site_file {sph_file}\n')
        fid.write('max_orientations 1000\n')
        fid.write('critical_points no\n')
        fid.write('chemical_matching no\n')
        fid.write('use_ligand_spheres no\n')
        fid.write('bump_filter no\n')
        fid.write('score_molecules yes\n')

        fid.write('contact_score_primary no\n')
        fid.write('contact_score_secondary no\n')

        fid.write('grid_score_primary yes\n')
        fid.write('grid_score_secondary no\n')
        fid.write('grid_score_rep_rad_scale 1\n')
        fid.write('grid_score_vdw_scale 1\n')
        fid.write('grid_score_es_scale 1\n')
        fid.write(f'grid_score_grid_prefix {grid_prefix}\n')

        fid.write('multigrid_score_secondary no\n')
        fid.write('dock3.5_score_secondary no\n')
        fid.write('continuous_score_secondary no\n')
        fid.write('footprint_similarity_score_secondary no\n')
        fid.write('pharmacophore_score_secondary no\n')
        fid.write('descriptor_score_secondary no\n')
        fid.write('gbsa_zou_score_secondary no\n')
        fid.write('gbsa_hawkins_score_secondary no\n')
        fid.write('SASA_score_secondary no\n')
        fid.write('amber_score_secondary no\n')

        fid.write('minimize_ligand yes\n')
        fid.write('minimize_anchor yes\n')
        fid.write('minimize_flexible_growth yes\n')
        fid.write('use_advanced_simplex_parameters no\n')

        fid.write('simplex_max_cycles 1\n')
        fid.write('simplex_score_converge 0.1\n')
        fid.write('simplex_cycle_converge 1.0\n')
        fid.write('simplex_trans_step 1.0\n')
        fid.write('simplex_rot_step 0.1\n')
        fid.write('simplex_tors_step 10.0\n')
        fid.write('simplex_anchor_max_iterations 500\n')
        fid.write('simplex_grow_max_iterations 500\n')
        fid.write('simplex_grow_tors_premin_iterations 0\n')
        fid.write('simplex_random_seed 0\n')
        fid.write('simplex_restraint_min no\n')

        fid.write('atom_model all\n')
        fid.write(f'vdw_defn_file {VDW_DEFN_FILE}\n')
        fid.write(f'flex_defn_file {FLEX_DEFN_FILE}\n')
        fid.write(f'flex_drive_file {FLEX_DRIVE_FILE}\n')

        fid.write(f'ligand_outfile_prefix {outfile}\n')
        fid.write('write_orientations no\n')
        fid.write('num_scored_conformers 1\n')
        fid.write('rank_ligands no\n')
    
    return str(input_file)
