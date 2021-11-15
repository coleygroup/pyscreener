---
title: 'pyscreener: A Python Wrapper for Computational Docking Software'
tags:
  - Python
  - drug discovery
  - distributed computing
  - molecular docking
authors:
  - name: David E. Graff
    orcid: 0000-0003-1250-3329
    affiliation: "1, 2"
  - name: Connor W. Coley^[corresponding author]
    orcid: 0000-0002-8271-8723
    affiliation: "2, 3"
affiliations:
  - name: Department of Chemistry and Chemical Biology, Harvard University
    index: 1
  - name: Department of Chemical Engineering, Massachusetts Institute of Technology
    index: 2
  - name: Department of Electrical Engineering and Computer Science, Massachusetts Institute of Technology
    index: 3
date: XX November 2021
bibliography: refs.bib
---

# Summary
`pyscreener` is a Python library that seeks to alleviate the challenges of large-scale structure-based design using computational docking. It provides a simple and uniform interface that is agnostic to the backend docking engine with which to calculate the docking score of a given molecule in a specified active site. Additionally, `pyscreener` features first-class support for task distribution, allowing users to seamlessly scale their code from a local, multi-core setup to a large, heterogeneous resource allocation.

# Statement of Need
Computational docking is an important technique in structure-based drug design that enables the rapid approximation of binding affinity for a candidate ligand in a matter of CPU seconds. With the growing popularity of ultra-large ligand libraries, docking is increasingly used to sift through hundreds of millions of compounds to try to identify novel and potent binders for a variety of protein targets \cite{gorgulla_open-source_2020,lyu_ultra-large_2019}. There are many choices of docking software, and certain software are better suited towards specific protein-ligand contexts (e.g., flexible protein side chains or sugar-like ligand molecules). Switching between these software is often not trivial as the input preparation, simulation, and output parsing pipelines differ between each software.

In addition, many of these programs exist only as command-line applications and lack Python bindings.
This presents an additional challenge for their integration into molecular optimization workflows, such as reinforcement learning or genetic algorithms.
Molecular optimization objectives have largely been limited to benchmark tasks, such as penalized logP, QED, JNK3 or GSK3$\beta$ inhibitor classification \cite{li_multi-objective_2018}, and others contained in the GuacaMol library \cite{brown_guacamol_2019}. These benchmarks are useful for comparing molecular design techniques, but they are not representative of true drug discovery tasks in terms of complexity; computational docking is at least one step in the right direction.

While many molecular optimization techniques propose new molecules in the form of SMILES strings \cite{elton_deep_2019}, most docking programs accept input in the form of molecular supply files with predefined 3D geometry (e.g., Mol2 or PDBQT format). Using the docking score of a molecule as a design objective thus requires an ad hoc implementation for which no standardized approach exists.
The `vina` library \cite{eberhardt_autodock_2021} is currently the only library capable of performing molecular docking within Python code, but it is limited to docking molecules using solely AutoDock Vina as the backend docking engine.
Moreover, the object model of the `vina` library accepts input ligands only as PDBQT files or strings and still does not address the need to quickly calculate the docking score of a molecule from its SMILES string.

In our work on the MolPAL software \cite{graff_accelerating_2021}, we required a library that is able to accept molecular inputs as SMILES strings and output their corresponding docking scores for a given receptor and docking box. Our use-case also called for docking large batches of molecules across large and distributed hardware setups. Lastly, we desired that our library be flexible with respect to the underlying docking engine, allowing us to use a variety of backend docking software (e.g., Vina \cite{trott_autodock_2010}, Smina \cite{koes_lessons_2013}, QVina \cite{alhossary_fast_2015}, or DOCK6 \cite{allen_dock_2015}) with minimal changes to client code. To that end, we developed `pyscreener`, a Python library that is flexible with respect to both molecular input format and docking engine that transparently handles the distribution of docking simulations across large resource allocations.
# Implementation and Performance

# Examples

To illustrate `pyscreener`, we consider docking benezene (SMILES string `"c1ccccc1"`) against 5WIU with a docking box centered at (-18.2, 14.4, -16.1) with x-, y-, and z-radii (15.4, 13.9, 14.5). We may perform this docking using AutoDock Vina over 6 CPU cores via `pyscreener` like so:
```python
>>> import pyscreener as ps
>>> metadata = ps.build_metadata("vina")
>>> virtual_screen = ps.virtual_screen(
...   "vina", receptors=["5WIU.pdb"],
...   center=(-18.2, 14.4, -16.1),
...   size=(15.4, 13.9, 14.5),
...   metadata_template=metadata,
...   ncpu=6
... )
>>> scores = virtual_screen("c1ccccc1")
>>> scores
array([-4.4])
```

Alternatively, we may dock many molecules by passing a `List` of SMILES strings to the `DockingVirtualScreen`:
```python
>>> smis = [
...   "c1ccccc1",
...   "O=C(Cc1ccccc1)NC1C(=O)N2C1SC(C2C(=O)O)(C)C",
...   "C=CCN1CCC23C4C(=O)CCC2(C1CC5=C3C(=C(C=C5)O)O4)O"
... ]
>>> scores = virtual_screen(smis)
>>> scores.shape
(3,)
```

By default, AutoDock Vina docks molecules using an `--exhaustiveness` value of 8, but we may specify a higher number in the `metadata`:
```python
>>> metadata = ps.build_metadata("vina", dict(exhaustivness=32))
```
We may also utilize other docking engines in the AutoDock Vina family by specifying the `software` for Vina-type metadata. Here, we use the accelerated optimization routine of QVina for faster docking. Note that we also support `software` values of `"smina"` and `"psovina"` in addition to `"vina"` and `"qvina"`.
```python
>>> metadata = ps.build_metadata("vina", dict(software="qvina"))
```

It is also possible to dock molecules using DOCK6 in `pyscreener`. To do this, we must first construct DOCK6 metadata and specify that we are creating a DOCK6 virtual screen (note that DOCK6 is not multithreaded and thus does not benefit from being assigned multiple CPU cores per task):
```python
>>> metadata = ps.build_metadata("dock")
>>> virtual_screen = ps.virtual_screen(
...   "dock",
...   receptors=["5WIU.pdb"],
...   center=(-18.2, 14.4, -16.1),
...   size=(15.4, 13.9, 14.5),
...   metadata_template=metadata
... )
>>> scores = virtual_screen("c1ccccc1")
>>> scores
array([-12.35])
```

# Acknowledgements

The authors thank Keir Adams and Wenhao Gao for providing feedback on the preparation of this paper and the `pyscreener` code. The computations in this paper were run on the FASRC Cannon cluster supported by the FAS Division of Science Research Computing Group at Harvard University. The authors also acknowledge the MIT SuperCloud and Lincoln Laboratory Supercomputing Center for providing HPC and consultation resources that have contributed to the research results reported within this paper. This work was funded by the MIT-IBM Watson AI Lab.