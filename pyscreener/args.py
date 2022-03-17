from datetime import datetime
import json
from typing import Optional

from configargparse import ArgumentParser, ArgumentTypeError, Namespace

__version__ = "1.2.0"


def gen_args(argv: Optional[str] = None) -> Namespace:
    parser = ArgumentParser(description="Automate virtual screening of compound libraries.")

    add_general_args(parser)
    add_preprocessing_args(parser)
    add_supply_args(parser)
    add_screen_args(parser)
    add_postprocessing_args(parser)

    args = parser.parse_args(argv)
    args.title_line = not args.no_title_line
    del args.no_title_line

    args.metadata_template["buffer"] = args.buffer
    args.metadata_template["docked_ligand_file"] = args.docked_ligand_file

    return args


def add_general_args(parser: ArgumentParser):
    parser.add_argument(
        "--config", is_config_file=True, help="filepath of a configuration file to use"
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--smoke-test",
        action="store_true",
        help="whether to perform a smoke test by checking if the environment is set up properly",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default=f'pyscreener_{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}',
        help="the path of the output directory",
    )
    parser.add_argument(
        "--no-sort",
        action="store_true",
        default=False,
        help="do not sort the output scores CSV file by score",
    )
    parser.add_argument(
        "--collect-all",
        action="store_true",
        default=False,
        help="whether all prepared input files and generated output files should be collected to the final output directory. By default, these files are all stored in a node-local temporary directory that is inaccessible after program completion.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="the level of output this program should print",
    )


def add_preprocessing_args(parser: ArgumentParser):
    parser.add_argument(
        "--preprocessing-options",
        nargs="+",
        choices=("pdbfix", "filter"),
        help="the preprocessing options to apply",
    )
    parser.add_argument(
        "--pH",
        type=float,
        default=7.0,
        help="the pH for which to calculate protonation state for protein and ligand residues",
    )


def add_supply_args(parser: ArgumentParser):
    parser.add_argument("-s", "--smis", nargs="+", help="the SMILES strings of the ligands to dock")
    parser.add_argument(
        "-i", "--input-files", nargs="+", help="the filenames containing ligands to dock"
    )
    parser.add_argument(
        "--input-filetypes",
        nargs="+",
        help="the filetype of each input ligand. If unspecified, will attempt to determine the filetype for each file.",
    )
    parser.add_argument(
        "--no-title-line",
        default=False,
        action="store_true",
        help="whether there is no title line in the ligands CSV file",
    )
    parser.add_argument(
        "--smiles-col",
        type=int,
        default=0,
        help="the column containing the SMILES strings in the CSV file.",
    )
    parser.add_argument(
        "--name-col",
        type=int,
        default=1,
        help="UNUSED the column containing the molecule names/IDs in the CSV file. Molecules will be labeled as ligand_<i> otherwise.",
    )
    parser.add_argument(
        "--id-property",
        help='UNUSED the name of the property containing the molecule names/IDs in a SMI or SDF file (e.g., "CatalogID", "Chemspace_ID", "Name", etc.). Molecules will be labeled as ligand_<i> otherwise.',
    )
    parser.add_argument(
        "--use-3d",
        action="store_true",
        help="whether to use the input 3D geometry of each molecule. Note that, in principle, initial geometry of input molecules to flexible docking simulations is statisically insignificant. This option is useful for presevering tautomeric information about input molecules.",
    )
    parser.add_argument(
        "--optimize",
        action="store_true",
        help="whether the geometry of each molecule should be optimized using the RDKit MMFF94 forcefield first. Note that, in principle, initial geometry of input molecules to flexible docking simulations is statisically insignificant.",
    )


def add_screen_args(parser: ArgumentParser):
    parser.add_argument(
        "--screen-type",
        choices=("dock", "dock6", "ucsfdock", "vina", "qvina", "smina", "psovina"),
        required=True,
        help="the type of docking screen to perform",
    )
    parser.add_argument("--receptors", nargs="+", help="the filenames of the receptors")
    parser.add_argument(
        "--center",
        type=float,
        nargs=3,
        metavar=("CENTER_X", "CENTER_Y", "CENTER_Z"),
        help="the x-, y-, and z-coordinates of the center of the docking box",
    )
    parser.add_argument(
        "--size",
        type=float,
        nargs=3,
        metavar=("SIZE_X", "SIZE_Y", "SIZE_Z"),
        help="the x-, y-, and z-radii of the docking box",
    )
    parser.add_argument("--metadata-template", type=json.loads, default={})
    parser.add_argument(
        "--pdbids", nargs="+", help="the PDB IDs of the crystal structures to dock against"
    )
    parser.add_argument(
        "--docked-ligand-file",
        help="the filepath of a PDB file containing the docked pose of a ligand from which to automatically construct a docking box",
    )
    parser.add_argument(
        "--buffer",
        default=10.0,
        type=float,
        help="the amount of buffer space to add around the docked ligand when calculating the docking box",
    )
    parser.add_argument("-nc", "--ncpu", default=1, type=int)
    parser.add_argument("--base-name", default="ligand")
    parser.add_argument(
        "--score-mode",
        default="best",
        choices=("best", "avg", "boltzmann", "top-k"),
        help="The method used to calculate the score of a single docking run on a single receptor",
    )
    parser.add_argument(
        "--repeat-score-mode",
        default="best",
        choices=("best", "avg", "boltzmann", "top-k"),
        help="The method used to calculate the overall score from repeated docking runs",
    )
    parser.add_argument(
        "--ensemble-score-mode",
        default="best",
        choices=("best", "avg", "boltzmann", "top-k"),
        help="The method used to calculate the overall score from an ensemble of docking runs",
    )
    parser.add_argument(
        "--repeats", default=1, type=int, help="the number of times to repeat each docking run"
    )
    parser.add_argument(
        "-k",
        default=1,
        type=int,
        help="the number of top scores to average if using a top-k score mode",
    )


def add_postprocessing_args(parser: ArgumentParser):
    parser.add_argument(
        "--postprocessing-options",
        nargs="+",
        default="none",
        choices=["visualize"],
        help="the postprocessing options to apply",
    )
    # parser.add_argument(
    #     "--n-cluster", type=int, default=10, help="the number of clusters to form"
    # )
    parser.add_argument(
        "--hist-mode",
        choices=["image", "text"],
        help='the type of histogram to generate. "image" makes a histogram that is output as a PNG file and "text" generates a histogram using terminal output.',
    )


def positive_int(arg: str):
    val = int(arg)
    if val <= 0:
        raise ArgumentTypeError(f"Value must be greater than 0! got: {arg}")

    return val
