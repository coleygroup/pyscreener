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
- [Running pyscreener](#running-pyscreener-as-a-software)
- [Using pyscreener](#using-pyscreener-as-a-library)

## Requirements
- python>=3.6
- all virtual screening software and openbabel located on your PATH
- (if using DOCK6-based HTVS) the location of the DOCK6 parent directory in your environment variables

### adding an executable to your PATH
`pyscreener` is not a virtual screening software in itself. Rather, it is a wrapper around common VS software to enable a simple and common interface between them without having to learn the ins and outs of the preparation+simulation pipeline for each different software. With that in mind, it is up to the user to install the appropriate virtual screening software and place them on their PATH.

After you have downloaded the appropriate software relevant to preparing inputs for and actually running your VS software of choice, you have two options:
1. append the directory containing these software to your bin via the command
    `export PATH=$PATH:<dir>`,
    where <dir> is the directory containing the software in question.
2. alternatively, and perhaps easier, is to copy the software to a directory that is already on your path. Typically this will be either `~/bin` or `~/.local/bin`. To see what directories are currently on your path, type `echo $PATH`. There will typically be a lot of directories on your path, and it is best to avoid creating files in any directory above your home directory ($HOME on most *nix-based systems)

### specifying an environment variable
Due to some wonkiness with getting DOCK6 to work with `pyscreener`, it requires that the DOCK6 environment variable be set with the location of the DOCK6 parent folder (the folder that is unpacked after downloading the original zip file and contains both the `bin` and `parameters` subdirectories.) To set the environment variable, enter the following command: `export DOCK6=<path/to/dock6>`. (_note_: this environment vairable must always be set before running pyscreener, so it's probably best to place this inside your `.bashrc` or `.bash_profile`)

__Note: both of the above steps must be satisifed before using `pyscreener`__

To avoid having to do this every time you start a new shell, you can add whatever commands you typed to your respective shell's startup file (e.g., .bash_profile for a bash shell) (you can also add them to the non-login shell startup file, but it's not good a idea to edit your PATH in these files)

## Installation
The first step in installing pyscreener is to clone this repository: `git clone <this_repo>`

### virtual environment setup
The easiest way to install all dependencies is to use conda along with the supplied [environment.yml](environment.yml) file, but you may also install them manually, if desired. All libraries listed in that file are __required__ before using `pyscreener`

0. (if necessary) [install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
1. `cd /path/to/pyscreener`
1. `conda env create -f environment.yml`

Before running `pyscreener`, be sure to first activate the environment: `conda activate pyscreener` (or whatever you've named your environment)

### external software
* vina-type software
  1. install [ADFR Suite](https://ccsb.scripps.edu/adfr/downloads/) for receptor preparation and follow the directions to add the resulting `bin` folder to your path (you should see a command at the end of the installation process)
  1. install any of the following docking software: [vina](http://vina.scripps.edu/), [qvina2](https://qvina.github.io/), [smina](https://sourceforge.net/projects/smina/), [psovina](https://cbbio.online/software/psovina/index.html) and ensure the desired software executable is in a folder that is located on your path
* [DOCK6](http://dock.compbio.ucsf.edu/)
  1. install [DOCK6](http://dock.compbio.ucsf.edu/) and specify the DOCK6 environment variable as the path of the parent folder (the one containing `bin`, `install`, etc.) as detailed [above](###specifying-an-environment-variable)
  1. install [sphgen_cpp](http://dock.compbio.ucsf.edu/Contributed_Code/sphgen_cpp.htm) and place the executable inside the `bin` subdirectory of the DOCK6 parent directory
  1. install [chimera](https://www.cgl.ucsf.edu/chimera/) and place the executable on your PATH (either by moving the exectuable to a folder already on your PATH or adding the folder containing the exectuable to your PATH)


## Running pyscreener as a software
__!!please read the entire section before running pyscreener!!__

pyscreener was designed to have a minimal interface under the principal that a high-throughput virtual screen is intended to be a broad strokes technique to gauge ligand favorability. With that in mind, all one really needs to get going are the following:
- the PDB id of your receptor of interest or a PDB format file of the specific structure
- a file containing the ligands you would like to dock, in SDF, SMI, or CSV format
- the coordinates of your docking box (center + size), a PDB format file containing the coordinates of a previously bound ligand, or a numbered list of residues from which to construct the docking box (e.g., [42, 64, 117, 169, 191])

There are a variety of other options you can specify as well (including how to score a ligand given that multiple scored conformations are output, how many times to repeatedly dock a given ligand, etc.) To see all of these options and what they do, use the following command: `python run.py --help`

All of these options may be specified on the command line, but they may also be placed in a configuration file that accepts YAML, INI, and `argparse` syntaxes. Example configuration files are located in [test_configs](test_configs). Assuming everything is working and installed properly, you can run any of these files via the following command: `python run.py --config test_configs/<config>`

### Parallel setup
pyscreener uses [`ray`](https://docs.ray.io/en/master/index.html) as its parallel backend. If you plan to parallelize the software only across your local machine, nothing additional needs to be done before starting the program. However, if you wish to either (a.) limit the number of cores pyscreener will be run over or (b.) run it over a distributed setup (e.g., an HPC with many distinct nodes), you must manually start a ray cluster.

#### Limiting the number of cores
To do this, simply type `ray start --head --num-cpus <N>` before starting pyscreener (where N is the total number of cores you wish to allow pyscreener to utilize). Not performing this step will give pyscreener access to all of the cores on your local machine, potentially slowing down other applications.

#### Distributing across many nodes
While the precise instructions for this will vary with HPC cluster architecture, the general idea is to establish a ray cluster between the nodes allocated to your job. We have provided a sample SLURM submission script ([run_pyscreener_distributed_example.batch](run_pyscreener_distributed_example.batch)) to achieve this, but you may have to alter some commands depending on your system. For more information on this see [here](https://docs.ray.io/en/master/cluster/index.html)

pyscreener writes a lot of intermediate input and output files. Given that the primary endpoint of pyscreener is a list of ligands and associated scores (rather than the specific binding poses,) these files are written to your system's temporary directory (determined by `tempfile.gettempdir()`). If you are running pyscreener in a distributed setup, either check for yourself or contact your system administrator to see where this directory points to. pyscreener requires that this be a _cluster-global_ directory __rather__ than a _node-local_ directory (i.e., the directory must be visible to all nodes on the cluster rather than only the local node.) If it is the latter, you must specify the proper directory via the `--tmp` or `--tmpdir` argument. It is typically best to avoid pointing this directory to your home directory both for storage and efficiency reasons.

## Using pyscreener as a library
To perform docking calls inside your python code using `pyscreener`, you must first initialize a `Screener` object using either of the derived classes: [`Vina`](pyscreener/docking/vina.py) or [`DOCK`](pyscreener/docking/dock.py).
    
`Vina` is the `Screener` class for performing docking simulations using any software derived from AutoDock Vina and accepts the `software` keyword argument to its initializer. Currently, the list of supported Vina-type software is as follows: AutoDock Vina, Smina, QVina2, and PSOVina.
    
`DOCK` is the `Screener` class for performing DOCKing using the DOCK software from UCSF. The input preparation pipeline for this software is a little more involved, so we encourage readers to look at the file to see what these additional parameters are.

For example, the following code snippet will dock benzene (SMILES string c1ccccc1) against the D4 dopamine receptor (PDB ID: 5WIU) using the site of a previously docked ligand.
```python
>>> import ray
>>> ray.init()
[...]
>>> from pyscreener import docking
>>> screener = docking.Vina(software='vina', receptors=['testing_inputs/5WIU.pdb'], docked_ligand_file='testing_inputs/5WIU_with_ligand.pdb', buffer=10., path='testing_outputs', ncpu=4)
Autoboxing ... Done!
Autoboxed ligand from "testing_inputs/5WIU_with_ligand.pdb" with center=(-18.2, 14.4, -16.1) and size=(15.4, 13.9, 14.5)
>>> results = screener('c1ccccc1')
>>> results
{'c1ccccc1': -4.4}
>>> results = screener('testing_inputs/ligands.csv')
>>> results
{...}
```
    
A few notes from the above example:
- the input PDB file must be *clean* prior to use. You can alternatively pass in a PDB ID (e.g., receptors=['5WIU']) but you must know the coordinates of the docking box for the corresponding PDB file. This usually means downloading the PDB file and manually inspecting it for more reliable results, but it's there if you want it.
- you can manually input the center and size of your docking box, but this must be manually determined before runtime.
- the prepared input/output files are stored in $TMPDIR by default. You can manually specify this via the `tmp_dir` argument during the `Vina` intialization. If you want these files at the end of execution, call the function `Screener.collect_all()`. This will collect all the input and output folders and move them under the directory specified by the `path` argument.
- If you don't want any files from `pyscreener` at all (only the score dictionary return value), don't set the `path` argument value.
- ray handles task distribution in the backend of the library. You don't need to manually start it if you're just going to call `ray.init()` like we did above. This was only done to highlight that you can initialize ray according to your own needs (i.e., distributed setup).
- you can call the `screener` object on (1) a SMILES string, (2) a csv/SDF/SMI file containing ligands, (3) a list of smiles strings, or (4) any combination of the above (e.g., `screener(ligand1, ligand_source_file, ligands_list)`). It is much more efficient to handle **one large** set of ligands than many small sets (i.e., `screener(one_long_list)` vs `screener(smiles1, smiles2, smiles3, ..., smilesN)`)
    
## Copyright
Copyright (c) 2020, david graff

## Acknowledgements 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.
