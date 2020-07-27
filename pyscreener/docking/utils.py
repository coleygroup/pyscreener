import subprocess as sp
from typing import Tuple

# a Ligand is a tuple of a molecule's SMILES string 
# and its corresponding input file
Ligand = Tuple[str, str]

OBABEL = sp.run('which obabel', shell=True, encoding='utf-8', 
                stdout=sp.PIPE, check=True).stdout.strip()