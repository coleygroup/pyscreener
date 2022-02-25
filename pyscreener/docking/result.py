from dataclasses import dataclass
from typing import Optional


@dataclass
class Result:
    smiles: str
    name: str
    node_id: str
    score: Optional[float]
