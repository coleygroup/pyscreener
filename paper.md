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
    affiliation: 2
affiliations:
  - name: Department of Chemistry and Chemical Biology, Harvard University
    index: 1
  - name: Department of Chemical Engineering, Massachusetts Institute of Technology
    index: 2
date: XX November 2021
bibliography: refs.bib
---

# Summary


# Statement of Need

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