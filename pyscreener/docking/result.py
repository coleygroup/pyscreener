from dataclasses import dataclass
from typing import Optional, TypeVar

S = TypeVar('S')
T = TypeVar('T')

@dataclass
class Result:
    smiles: str
    name: str
    node_id: str
    score: Optional[float]