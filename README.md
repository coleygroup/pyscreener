[//]: # (Badges)
[![codecov](https://codecov.io/gh/coleygroup/pyscreener/branch/master/graph/badge.svg)](https://codecov.io/gh/coleygroup/pyscreener/branch/master)
![CI](https://github.com/coleygroup/pyscreener/workflows/CI/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/pyscreener/badge/?version=latest)](https://pyscreener.readthedocs.io/en/latest/?badge=latest)

# pyscreener
# A pythonic interface to high-throughput virtual screening software

## Overview
This repository contains the source of pyscreener, both a library and software for conducting HTVS via python calls

## Table of Contents
- [Overview](#overview)
- [Table of Contents](#table-of-contents)
- [Requirements](#requirements)
- [Installation](#installation)
- [Setup](#setup)
- [Running pyscreener](#running-pyscreener-as-a-software)
- [Using pyscreener](#using-pyscreener-as-a-library)

__Note: both of the above steps must be satisifed before using `pyscreener`__

To avoid having to do this every time you start a new shell, you can add whatever commands you typed to your respective shell's startup file (e.g., .bash_profile for a bash shell) (you can also add them to the non-login shell startup file, but it's not good a idea to edit your PATH in these files)

## Installation

### General requirements
- python>=3.6 and the packages in `environment.yml`, `requirements.txt` and [`pdbfixer`](https://github.com/openmm/pdbfixer)
- all corresponding software downloaded and located on your PATH or under the corresponding environment variable PATH (see [external software](#external-software) for more details.)
- (if using DOCK6-based HTVS) the location of the DOCK6 parent directory in your environment variables

### environment setup with conda

0. (if necessary) [install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
1. clone this repository: `git clone git@github.com:coleygroup/pyscreener.git`
1. `cd pyscreener`
1. `conda env create -f environment.yml`
1. `pip install .`
1. follow the corresponding directions below for the intended software

Before running `pyscreener`, be sure to first activate the environment: `conda activate pyscreener` (or whatever you've named your environment)

### external software
* vina-type software
  1. install [ADFR Suite](https://ccsb.scripps.edu/adfr/downloads/) for receptor preparation and follow the directions to add the resulting `bin` folder to your path (you should see a command at the end of the installation process)
  1. install any of the following docking software: [vina](http://vina.scripps.edu/), [qvina2](https://qvina.github.io/), [smina](https://sourceforge.net/projects/smina/), [psovina](https://cbbio.online/software/psovina/index.html) and ensure the desired software executable is in a folder that is located on your path

* [DOCK6](http://dock.compbio.ucsf.edu/)
  0. obtain a license for [DOCK6](http://dock.compbio.ucsf.edu/Online_Licensing/dock_license_application.html)
  1. install DOCK6 from the download link and follow the [installation directions](http://dock.compbio.ucsf.edu/DOCK_6/dock6_manual.htm#Installation)
  1. after ensuring the installation was installed properly, specify the DOCK6 environment variable as the path of the DOCK6 parent directory as detailed [below](#specifying-an-environment-variable). This is the directory that was unzipped from the tarball and is usually named `dock6`. It is the folder that contains the `bin`, `install`, etc. subdirectories.)
  1. install [sphgen_cpp](http://dock.compbio.ucsf.edu/Contributed_Code/sphgen_cpp.htm). On linux systems, this can be done:
    1. `wget http://dock.compbio.ucsf.edu/Contributed_Code/code/sphgen_cpp.1.2.tar.gz`
    1. `tar -xzvf sphgen_cpp.1.2.tar.gz`
    1. `cd sphgen_cpp.1.2`
    1. `make`
  1. place the sphgen_cpp executable (it should be `sphgen_cpp`) inside the `bin` subdirectory of the DOCK6 parent directory. If you've configured the environment variable already, you can run: `mv sphgen_cpp $DOCK6/bin`
  1. install [chimera](https://www.cgl.ucsf.edu/chimera/) and place the file on your PATH as detailed [below](#adding-an-executable-to-your-path)

#### adding an executable to your PATH
<!-- `pyscreener` is not a virtual screening software in itself. Rather, it is a wrapper around common VS software to enable a simple and common interface between them without having to learn the ins and outs of the preparation and simulation pipeline for each different software. With that in mind, it is up to the user to install the appropriate virtual screening software and place them on their PATH. -->

To add a file to your PATH, you have two options:
1. append the directory containing the file to your PATH: `export PATH=$PATH:<DIR>`, where `<DIR>` is the directory containing the file in question. As your PATH must be configured each time run pyscreener, this command should also be placed inside your `~/.bashrc` or `~/.bash_profile` (if using a bash shell) to avoid needing to run the command every time you log in. _Note_: if using a non-bash shell, the specific file will be different.
1. copy the software to a directory that is already on your path. Typically this will be either `~/bin` or `~/.local/bin`. To see what directories are currently on your path, type `echo $PATH`. There will typically be a lot of directories on your path, and it is best to avoid creating files in any directory above your home directory (`$HOME` on most *nix-based systems)

#### specifying an environment variable
To set the `DOCK6` environment variable, run the following command: `export DOCK6=<path/to/dock6>`, where `<path/to/dock6>` is the **full** path of the DOCK6 parent directory mentioned above. As this this environment variable must always be set before running pyscreener, the command should be placed inside your `~/.bashrc` or `~/.bash_profile` (if using a bash shell) to avoid needing to run the command every time you log in. _Note_: if using a non-bash shell, the specific file will be different.

## Setup
pyscreener uses [`ray`](https://docs.ray.io/en/master/index.html) as its parallel backend. If you plan to parallelize the software only across your local machine, you need only run `ray start --head` __before__ running pyscreener. However, if you wish to either (a.) limit the number of cores pyscreener will be run over or (b.) run it over a distributed setup (e.g., an HPC with many distinct nodes), you must manually start a ray cluster.

#### Limiting the number of cores
To do this, simply type `ray start --head --num-cpus <N>` before starting pyscreener (where `N` is the total number of cores you wish to allow pyscreener to utilize). Not performing this step will give pyscreener access to all of the cores on your local machine, potentially slowing down other applications.

#### Distributing across many nodes
While the precise instructions for this will vary with HPC cluster architecture, the general idea is to establish a ray cluster between the nodes allocated to your job. We have provided a sample SLURM submission script ([run_pyscreener_distributed_example.batch](run_pyscreener_distributed_example.batch)) to achieve this, but you may have to alter some commands depending on your system. For more information on this see [here](https://docs.ray.io/en/master/cluster/index.html)

pyscreener writes a lot of intermediate input and output files (due to the inherent specifications of the underlying docking software.) Given that the primary endpoint of pyscreener is a list of ligands and associated scores (rather than the specific binding poses,) these files are written to each node's temporary directory (determined by `tempfile.gettempdir()`) and discarded at the end. If you wish to collect these files, pass the `--collect-all` flag in the program arguments or run the `collect_files()` method of your `Screener` object when your screen is complete. *Note*: the `collect_files()` method is **slow** due to the need to send possibly a **bunch** of files over the network. This method should only be run **once** over the lifetime of a `Screener` object, as several intermediate calls will yield the same result as a single, final call.

## Running pyscreener as a software
__!!please read the entire section before running pyscreener!!__

pyscreener was designed to have a minimal interface under the principal that a high-throughput virtual screen is intended to be a broad strokes technique to gauge ligand favorability. With that in mind, all one really needs to get going are the following:
- the PDB id of your receptor of interest or a PDB format file of the specific structure
- a file containing the ligands you would like to dock, in SDF, SMI, or CSV format
- the coordinates of your docking box (center + size), a PDB format file containing the coordinates of a previously bound ligand, or a numbered list of residues from which to construct the docking box (e.g., [42, 64, 117, 169, 191])

There are a variety of other options you can specify as well (including how to score a ligand given that multiple scored conformations are output, how many times to repeatedly dock a given ligand, etc.) To see all of these options and what they do, use the following command: `python run.py --help`

All of these options may be specified on the command line, but they may also be placed in a configuration file that accepts YAML, INI, and `argparse` syntaxes. Example configuration files are located in [test_configs](test_configs). Assuming everything is working and installed properly, you can run any of these files via the following command: `python run.py --config test_configs/<config>`

## Using pyscreener as a library
The object model of pyscreener relies on four classes:
* [`CalculationData`](pyscreener.docking.data.py): a simple object containing the broadstrokes specifications of a docking calculation common to all types of docking calculations (e.g., Vina, DOCK6, etc.): the SMILES string, the target receptor, the center/size of a docking box, the metadata, and the result.
* [`CalculationMetadata`](pyscreener.docking.metadata.py): a non-descript object that contains software-specific fields. For example, a Vina-type calculation requires a `software` parameter, whereas a DOCK6 calculation requires a number of different parameters for receptor preparation. Most importantly, the metadata will always contain two fields of abstract type: `prepared_ligand` and `prepared_receptor`. You can inspect the respective metadata classes to learn the concrete type, but the intent is that a user will not be inspecting the metadata object
* [`DockingRunner`](pyscreener.docking.runner.py): a static object that takes defines an interface to prepare and run docking calculations. Each calculation type defines its own `DockingRunner` implementation.
* [`DockingVirtualScreen`](pyscreener.docking.screen.py): an object that organizes a virtual screen. At a high level, a virtual is a series of docking calculations with some template set of parameters performed for a collection of molecules and distributed over some set of computational resources. A `DockingVirtualScreen` takes as arguments a `DockingRunner`, a list of receptors (for possible ensemble docking) and a set of template values for a `CalculationData` template. It defines a `__call__(*args)` method that takes an unzipped list of SMILES strings, builds the `CalculationData` objects for each molecule, and submits these objects for preparation and calculation to various resources in the ray cluster (see [setup](#setup)).


To perform docking calls inside your python code using `pyscreener`, you must first initialize a `Screener` object using either of the derived classes: [`Vina`](pyscreener/docking/vina.py) or [`DOCK`](pyscreener/docking/dock.py).
    
`Vina` is the `Screener` class for performing docking simulations using any software derived from AutoDock Vina and accepts the `software` keyword argument to its initializer. Currently, the list of supported Vina-type software is as follows: AutoDock Vina, Smina, QVina2, and PSOVina.
    
`DOCK` is the `Screener` class for performing DOCKing using the DOCK software from UCSF. The input preparation pipeline for this software is a little more involved, so we encourage readers to look at the file to see what these additional parameters are.

For example, the following code snippet will dock benzene (SMILES string c1ccccc1) against the D4 dopamine receptor (PDB ID: 5WIU) using the site of a previously docked ligand.
```python
>>> import ray
>>> ray.init()
[...]
>>> from pyscreener import docking
>>> vina_screener = docking.screener(software='vina', receptors=['testing_inputs/5WIU.pdb'], docked_ligand_file='testing_inputs/5WIU_with_ligand.pdb', buffer=10., path='testing_outputs', ncpu=4)
Autoboxing ... Done!
Autoboxed ligand from "testing_inputs/5WIU_with_ligand.pdb" with center=(-18.2, 14.4, -16.1) and size=(15.4, 13.9, 14.5)
>>> results = vina_screener('c1ccccc1')
>>> results
{'c1ccccc1': -4.4}
>>> results = vina_screener('testing_inputs/test_ligands.csv')
>>> results
{...}
```
    
A few notes from the above example:
- the input PDB file must be *clean* prior to use. You can alternatively pass in a PDB ID (e.g., receptors=['5WIU']) but you must know the coordinates of the docking box for the corresponding PDB file. This usually means downloading the PDB file and manually inspecting it for more reliable results, but it's there if you want it.
- you can manually input the center and size of your docking box, but this must be manually determined before runtime. e.g.
    ```python
    screener = docking.Vina(software='vina', receptors=['testing_inputs/5WIU.pdb'], center=(-18.2, 14.4, -16.1), size=(15.4, 13.9, 14.5), path='testing_outputs', ncpu=4)
    ```
- the prepared input/output files are stored in $TMPDIR by default. You can manually specify this via the `tmp_dir` argument during the `Vina` intialization. If you want these files at the end of execution, call the function `Screener.collect_all()`. This will collect all the input and output folders and move them under the directory specified by the `path` argument.
- If you don't want any files from `pyscreener` at all (only the score dictionary return value), don't set the `path` argument value.
- ray handles task distribution in the backend of the library. You don't need to manually start it if you're just going to call `ray.init()` like we did above. This was only done to highlight that you can initialize ray according to your own needs (i.e., distributed setup).
- you can call the `screener` object on (1) a SMILES string, (2) a csv/SDF/SMI file containing ligands, (3) a list of smiles strings, or (4) any combination of the above (e.g., `screener(ligand1, ligand_source_file, ligands_list)`). It is much more efficient to handle **one large** set of ligands than many small sets (i.e., `screener(one_long_list)` vs `screener(smiles1, smiles2, smiles3, ..., smilesN)`)
    
## Copyright
Copyright (c) 2020, david graff

## Acknowledgements 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
