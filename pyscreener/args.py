import os
from pathlib import Path
import shlex

from configargparse import ArgumentParser, ArgumentTypeError, Namespace

def positive_int(arg: str):
    value = int(arg)
    if value <= 0:
        raise ArgumentTypeError('Value must be greater than 0!')
    return value

def add_general_args(parser: ArgumentParser):
    parser.add_argument('--config', is_config_file=True,
                        help='filepath of a configuration file to use')
    parser.add_argument('--name',
                        help='the base name of the outputs')
    parser.add_argument('--mode', default='docking',
                        choices=['docking', 'md', 'dft'],
                        help='the mode in which to run pyscreener')

    parser.add_argument('--distributed', action='store_true', default=True,
                        help='whether to parallelize computation using a distributed setup')
    parser.add_argument('--nworkers', type=int, metavar='N_WORKERS', default=-1,
                        help='the number of workers to use. (Only used when distributed=False.)')
    
    parser.add_argument('--root', default='.',
                        help='the root directory under which to organize all program outputs')
    parser.add_argument('--tmp', default=os.environ.get('TMP', '.'),
                        help='path to your system\'s TMP directory. If using a distributed setup, this directory must be shared between all machines.')
    parser.add_argument('-v', '--verbose', action='count', default=0,
                        help='the level of output this program should print')

def add_preprocessing_args(parser: ArgumentParser):
    parser.add_argument('--preprocessing-options', nargs='+',
                        default='none',
                        choices=['pdbfix', 'autobox', 'tautomers', 'desalt', 
                                 'filter'],
                        help='the preprocessing options to apply')

def add_preparation_args(parser: ArgumentParser):
    parser.add_argument('--no-title-line', default=False, action='store_true',
                        help='whether there is no title line in the ligands csv file')
    parser.add_argument('--start', type=int, default=0,
                        help='the index at which start conversion')
    parser.add_argument('--nconvert', type=int,
                        help='the number of molecules to convert')
    
def add_screening_args(parser: ArgumentParser):
    ### DOCKING ARGS ###
    parser.add_argument('--docker', default='vina',
                        choices=['vina', 'smina', 'qvina', 'psovina'],
                        help='the name of the docking program to use')
    parser.add_argument('-r', '--receptors', required=True, nargs='+',
                        help='the filenames of the receptors')
    parser.add_argument('-l', '--ligands', required=True, nargs='+',
                        help='the filenames containing the ligands to dock')
    parser.add_argument('--use-3d', action='store_true', default='False',
                        help='how to treat the preparation of ligands from files containing three-dimensional information. If False, use only the 2D graph of each molecule in the SDF file when preparing inputs. Faster, but will result in the loss of conformational/tautomeric information.If True, se the 3D information contained in the file when preparing an input. Slower, but will preserve conformational/tautomeric information.')
    parser.add_argument('-c', '--center', required=True,
                        type=float, nargs=3,
                        metavar=('CENTER_X', 'CENTER_Y', 'CENTER_Z'),
                        help='the x-, y-, and z-coordinates of the center of the docking box')
    parser.add_argument('-s', '--size', required=True,
                        type=int, nargs=3,
                        metavar=('SIZE_X', 'SIZE_Y', 'SIZE_Z'),
                        help='the x-, y-, and z-dimensions of the docking box')
    parser.add_argument('--ncpu', type=int, default=1, metavar='N_CPU',
                        help='the number of cores available to each worker process')
    parser.add_argument('--extra', type=shlex.split,
                        help='extra command line arguments to pass to screening software. E.g., "--exhaustiveness=16"')

    ### SCORING ARGS ###    
    parser.add_argument('--score-mode', default='best',
                        choices={'best', 'avg', 'boltzmann'},
                        help='The method used to calculate the score of a single docking run on a single receptor')
    parser.add_argument('--repeat', type=positive_int, default=1,
                        help='the number of times to repeat each screening run')
    parser.add_argument('--repeat-score-mode', default='best',
                        choices={'best', 'avg', 'boltzmann'},
                        help='The method used to calculate the overall score from multiple docking runs on the same receptor')
    parser.add_argument('--ensemble-score-mode', default='best',
                        choices={'best', 'avg', 'boltzmann'},
                        help='The method used to calculate the overall score from an ensemble of docking runs')

def add_postprocessing_args(parser: ArgumentParser):
    parser.add_argument('--postprocessing-options', nargs='+',
                        default='none',
                        choices=['cluster', 'distribution'],
                        help='the postprocessing options to apply')

def gen_args() -> Namespace:
    parser = ArgumentParser(
        description='Automate virtual screening of compound libraries.')

    add_general_args(parser)
    add_preprocessing_args(parser)
    add_preparation_args(parser)
    add_screening_args(parser)
    add_postprocessing_args(parser)

    args = parser.parse_args()

    if args.name is None:
        args.name = f'{Path(args.receptor).stem}_{Path(args.ligands).stem}'
        
    args.title_line = not args.no_title_line
    del args.no_title_line

    return args
