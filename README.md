pyscreener
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/coleygroup/pyscreener/workflows/CI/badge.svg)](https://github.com/coleygroup/pyscreener/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/coleygroup/pyscreener/branch/master/graph/badge.svg)](https://codecov.io/gh/coleygroup/pyscreener/branch/master)

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
Due to some wonkiness with getting DOCK6 to work with `pyscreener`, it requires that the DOCK6 environment variable be set with the location of the DOCK6 parent folder (the folder that is unpacked after downloading the original zip file and contains both the bin/ and parameters/ subfolders.) To set the environment variable, enter the following command: `export DOCK6=<path/to/dock6>`.


__Note: both of the above steps must be satisifed before using `pyscreener`__

To avoid having to do this every time you start a new shell, you can add whatever commands you typed to your respective shell's startup file (e.g., .bash_profile for a bash shell) (you can also add them to the non-login shell startup file, but it's not good a idea to recursively edit your PATH in these files)

## Installation
The first step in installing pyscreener is to clone this repository: `git clone <this_repo>`

The easiest way to install all dependencies is to use conda along with the supplied [environment.yml](environment.yml) file, but you may also install them manually, if desired. All libraries listed in that file are __required__ before using `pyscreener`

### virtual environment setup via conda 
0. (if necessary) [install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/)
1. `cd /path/to/pyscreener`
1. `conda env create -f environment.yml`

Before running `pyscreener`, be sure to first activate the environment: `conda activate pyscreener`

### external software
* vina-type software
  - [ADFR Suite](https://ccsb.scripps.edu/adfr/downloads/) for receptor preparation
  - [vina](http://vina.scripps.edu/)
  - [qvina2](https://qvina.github.io/)
  - [smina](https://sourceforge.net/projects/smina/)
  - [psovina](https://cbbio.online/software/psovina/index.html)
* [DOCK6](http://dock.compbio.ucsf.edu/)

## Running pyscreener as a software
pyscreener was designed to have a minimal interface under the principal that a high-throughput virtual screen is intended to be a broad strokes technique to gauge ligand favorability. With that in mind, all one really needs to get going are the following:
- the PDB id of your receptor of interest or a PDB format file of the specific structure
- a file containing the ligands you would like to dock, in SDF, SMI, or CSV format
- the coordinates of your docking box (center + size), a PDB format file containing the coordinates of a previously bound ligand, or a numbered list of residues from which to construct the docking box (e.g., [42, 64, 117, 169, 191])

There are a variety of other options you can specify as well (including how to score a ligand given that multiple scored conformations are output, how many times to repeatedly dock a given ligand, etc.) To see all of these options and what they do, use the following command: `python run.py --help`

All of these options may be specified on the command line, but they may also be placed in a configuration file that accepts YAML, INI, and `argparse` syntaxes. Example configuration files are located in [test_configs](test_configs). Assuming everything is working and installed properly, you can run any of these files via the following command: `python run.py --config test_configs/<config>`

## Using pyscreener as a library
At the core of the pyscreener software is the `pyscreener` library that enables the running of docking software from input preparation all the way to output file parsing. The workhorse class is the [`Screener`](pyscreener/docking/screener.py) ABC, which handles all of this for a user. To actually initialize a screener object, either of the derived classes: [`Vina`](pyscreener/docking/vina.py) or [`DOCK`](pyscreener/docking/dock.py). `Vina` is the `Screener` class for performing docking simulations using any software derived from AutoDock Vina and accepts the `software` keyword argument to its initializer. Currently, the list of supported Vina-type software is as follows: AutoDock Vina, Smina, QVina2, and PSOVina. `DOCK` is the `Screener` class for performing DOCKing using the DOCK software from UCSF. The input preparation pipeline for this software is a little more involved, so we encourage readers to look at the file to see what these additional parameters are.

## Copyright
Copyright (c) 2020, david graff

## Acknowledgements 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.5.