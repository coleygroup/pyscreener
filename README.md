# pyscreener
# A pythonic interface to high-throughput virtual screening software

## Overview
This repository contains the source of pyscreener, both a library and software for conducting HTVS via python calls

## Table of Contents
- [Overview](#overview)
- [Table of Contents](#table-of-contents)
- [Requirements](#requirements)
- [Installation](#installation)
- [Object Model](#object-model)
- [Running MolPAL](#running-molpal)
  * [Required Settings](#required-settings)
  * [Optional Settings](#optional-settings)
- [Future Directions](#future-directions)
- [Reproducing Experimental Results](#reproducing-experimental-results)

## Requirements
- Python (>= 3.6)
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
* [openbabel](https://github.com/openbabel/openbabel) _Note_: `pyscreener` currently relies on subprocess calls to the `obabel` executable. However, this will be replaced by internal calls to the `pybel` library in a future release. This shouldn't affect your installation if you used the `environment.yml` file with conda.
* vina-type software
  - [ADFR Suite](https://ccsb.scripps.edu/adfr/downloads/) for receptor preparation
  - [vina](http://vina.scripps.edu/)
  - [qvina2](https://qvina.github.io/)
  - [smina](https://sourceforge.net/projects/smina/)
  - [psovina](https://cbbio.online/software/psovina/index.html)
* [DOCK6](http://dock.compbio.ucsf.edu/)
