Getting Started
===============

This page details how to get started with pyscreener. 

Installation
------------
The first step in installing pyscreener is to clone this repository: ``git clone <this_repo>``

The easiest way to install all dependencies is to use conda along with the supplied [environment.yml](environment.yml) file, but you may also install them manually, if desired. All libraries listed in that file are **require** before using `pyscreener`

virtual environment setup via conda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#. (if necessary) `install conda <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_
#. ``cd /path/to/pyscreener``
#. ``conda env create -f environment.yml``

Before running ``pyscreener``, be sure to first activate the environment: `conda activate pyscreener`

external software
^^^^^^^^^^^^^^^^^
* vina-type software

  #. install `ADFR Suite <https://ccsb.scripps.edu/adfr/downloads/>`_ for receptor preparation
  #. install any of the following docking software: `vina <http://vina.scripps.edu/), [qvina2](https://qvina.github.io/>`_, `smina <https://sourceforge.net/projects/smina/>`_, `psovina <https://cbbio.online/software/psovina/index.html>`_ and ensure the desired software is located on your path

* DOCK6

  #. install `DOCK6 <http://dock.compbio.ucsf.edu/>`_ and specify the DOCK6 environment variable as the path of the parent folder (the one containing ``bin``, ``install``, etc.) as detailed
  #. install `sphgen_cpp <http://dock.compbio.ucsf.edu/Contributed_Code/sphgen_cpp.html>`_ and place the executable inside the ``bin`` subdirectory of the DOCK6 parent directory
  #. install `chimera <https://www.cgl.ucsf.edu/chimera/>`_ and place the executable on your PATH (either by moving the exectuable to a folder already on your PATH or adding the folder containing the exectuable to your PATH)