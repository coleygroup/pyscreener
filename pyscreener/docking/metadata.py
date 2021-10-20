from abc import ABC
from typing import Optional, TypeVar

S = TypeVar('S')
T = TypeVar('T')

class CalculationMetadata(ABC):
    prepared_ligand: Optional[T] = None
    prepared_receptor: Optional[S] = None