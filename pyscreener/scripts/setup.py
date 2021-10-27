import argparse
from contextlib import contextmanager
import os
from pathlib import Path
import shutil
import subprocess as sp
import sys
import tarfile

import colorama
import requests

colorama.init(autoreset=True)


@contextmanager
def cd(newdir):
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)


def setup():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(
        title="software",
        dest="software",
        description="software to set up for pyscreener",
    )
    vina_parser = subparsers.add_parser("vina", aliases=["qvina", "smina", "psovina"])
    vina_parser.add_argument(
        "source_file",
        help="""the source tarball or binary of the selected docking software.
Vina: the Vina 1.1.2 source tarball for your system from https://vina.scripps.edu/downloads/.
QVina2: the source binary from https://github.com/QVina/qvina/blob/master/bin/qvina2.1.
Smina: the source binary from https://sourceforge.net/projects/smina/files/smina.static/download.
PSOVina: the source binary from https://sourceforge.net/projects/psovina/files/binaries/psovina-2.0/Ubuntu-15.04_psovina2/download""",
    )
    vina_parser.add_argument(
        "adfr_tarball",
        help="the ADFR tarball for your system from https://ccsb.scripps.edu/adfr/downloads/",
    )

    dock_parser = subparsers.add_parser("dock", aliases=["dock6", "ucsfdock"])
    dock_parser.add_argument("source_tarball", help="the dock6 source tarball")
    dock_parser.add_argument("chimera_installer", help="the dock6 source tarball")

    args = parser.parse_args()

    if args.software.lower() in ("vina", "qvina", "smina", "psovina"):
        setup_vina(args.software.lower(), args.source_file, args.adfr_tarball)
    else:
        raise NotImplementedError("DOCK6 automatic setup is not supported currently!")


def setup_vina(software, source_file, adfr_tarball):
    path_dir = get_dir_on_path()

    print("Extracting software from tarball and creating symlink ...")
    p_symlink_exe = path_dir / software
    if p_symlink_exe.is_symlink():
        print(f"{software} already exists on path! Skipping this step ...")
    else:
        if software == "vina":
            p_exe = extract_vina_exe(source_file)
            p_symlink_exe.symlink_to(p_exe)
        elif software in {"smina", "qvina", "psovina"}:
            p_exe = Path(source_file)
            p_symlink_exe.symlink_to(p_exe)
        else:
            raise RuntimeError

    p_prep_rec_symlink = path_dir / "prepare_receptor"
    if p_prep_rec_symlink.exists():
        print(f"prepare_receptor already exists on path! Skipping this step ...")
    else:
        p_prep_rec_symlink.symlink_to(extract_adfr_tarball(adfr_tarball))

    print(
        'Finished! To check if everything worked, run "pyscreener check vina METADATA" next!'
    )


def extract_vina_exe(source_tarball):
    p_source = Path(source_tarball)
    with tarfile.open(p_source, "r:*") as tar:
        for m in tar.getmembers():
            if Path(m.name).stem == "vina":
                tar.extract(m, p_source.parent)
                return p_source.parent / Path(m.name)

    return RuntimeError(f"No vina exectuable could be found in tarball {p_source}!")


def extract_adfr_tarball(adfr_tarball):
    DEFAULT = str(Path.home() / "ADFR")
    install_destination = (
        input(f"Destination under which to install ADFR suite (default = {DEFAULT}): ")
        or DEFAULT
    )

    p_adfr = Path(adfr_tarball)
    with tarfile.open(p_adfr, "r:*") as tar:
        install_script = [m for m in tar.getmembers() if "install.sh" in m.name][0].name
        tar.extractall(p_adfr.parent)

    p_install_script = p_adfr.parent / install_script
    with cd(p_install_script.parent):
        sp.run(["./install.sh", "-d", install_destination, "-c", "0"])

    p_prep_rec_py = Path(install_destination) / "bin" / "prepare_receptor"
    if not p_prep_rec_py.exists():
        raise RuntimeError(
            f'Could not find "prepare_receptor" at {p_prep_rec_py}. Setup failed!'
        )

    return p_prep_rec_py


def get_dir_on_path():
    paths = os.environ["PATH"].split(os.pathsep)
    print("Found following directories on PATH:")
    for i, p in enumerate(paths):
        if Path.home() in Path(p).parents:
            print(f"({i}) {p} " + colorama.Fore.RED + colorama.Style.BRIGHT + "*")
        else:
            print(f"({i}) {p}")

    idx = int(
        input(
            "Select index of directory under which to create symlinks "
            + "(preferred directories marked by "
            + colorama.Fore.RED
            + colorama.Style.BRIGHT
            + "*"
            + colorama.Style.RESET_ALL
            + "): "
        )
    )

    return Path(paths[idx])


def setup_dock6(source_tarball, chimera_installer):
    if "DOCK6" in os.environ:
        print(
            "DOCK6 environment variable detected! Checking for proper configuration ..."
        )
        dock6_dir = Path(os.environ["DOCK6"])
        if dock6_dir_configured_properly(dock6_dir):
            print("DOCK6 directory is configured properly!")
        else:
            print(
                f"DOCK6 directory ({dock6_dir}) is not configured properly! ",
                "Reinstalling from provided tarball ...",
                end=" ",
            )
            install_dock6(source_tarball)
            print("Done!")
    else:
        print("Installing dock6 from source tarball ...", end=" ")
        dock6_dir = install_dock6(source_tarball)
        print("Done!")

    p_spghen_symlink = dock6_dir / "bin" / "sphgen_cpp"
    if p_spghen_symlink.exists():
        print("sphgen_cpp already exists under dock6/bin! Skipping this step ...")
    else:
        print("Installing sphgen_cpp ...", end=" ")
        p_spghen_symlink.symlink_to(install_sphgen_cpp())
        print("Done!")

    if shutil.which("chimera") is None:
        print("Installing chimera ...", end=" ")
        install_chimera()
        print("Done!")
    else:
        print("chimera is already installed! Skipping this step ...")

    print(
        'Finished! To check if everything worked, run "pyscreener check dock6 METADATA" next!'
    )


def dock6_dir_configured_properly(dock6_dir):
    pass


def install_dock6(source_tarball) -> Path:
    pass


def install_sphgen_cpp():
    SPHGEN_URL = (
        "http://dock.compbio.ucsf.edu/Contributed_Code/code/sphgen_cpp.1.2.tar.gz"
    )
    try:
        resp = requests.get(SPHGEN_URL)

        DEFAULT = Path.home()
        download_dir = Path(
            input(f"Download destination for sphgen_cpp (default = {DEFAULT}): ")
            or DEFAULT
        )
        p_sphgen_tarball = download_dir / "sphgen_cpp.1.2.tar.gz"
        with open(p_sphgen_tarball, "wb") as fid:
            [fid.write(chunk) for chunk in resp.iter_content(20000)]
    except:
        print("error: could not download and install sphgen_cpp! aborting setup!")
        raise

    with tarfile.open(p_sphgen_tarball, "r:*") as tar:
        makefile = [m for m in tar.getmembers() if "Makefile" in m.name][0].name
        tar.extractall(p_sphgen_tarball.parent)

    p_makefile = p_sphgen_tarball.parent / makefile
    with cd(p_makefile.parent):
        sp.run(["make"], check=True)
        p_sphgen = Path("spghen").absolute()

    if not p_sphgen.exists():
        raise RuntimeError("error: could not install sphgen_cpp! aborting setup!")

    return p_sphgen


def install_chimera():
    pass


if __name__ == "__main__":
    setup()
