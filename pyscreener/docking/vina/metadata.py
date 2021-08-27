from dataclasses import dataclass
from enum import auto, Enum
from pathlib import Path
import shlex
from typing import Optional, Tuple

from pyscreener.docking.metadata import CalculationMetadata

class Software(Enum):
    def _generate_next_value_(name, start, count, last_values):
        return name.lower()

    VINA = auto()
    PSOVINA = auto()
    QVINA = auto()
    SMINA = auto()

@dataclass(repr=True, eq=False)
class VinaMetadata(CalculationMetadata):
    """

    Attributes
    ---------
    software : Software
        the software that will be used
    center : Tuple[float, float, float]
        the center of the docking box
    size : Tuple[float, float, float]
        the x-, y-, and z-radii of the docking box
    ncpu : int
        the number of cpu cores to use during the docking calculation
    extra : Optional[List[str]]
        additional command line arguments that will be passed to the
        docking calculation
    prepared_ligand: Optional[Path]
    prepared_receptor: Optional[Path]

    Parmeters
    ---------
    software : Software
        the software that will be used
    center : Tuple[float, float, float]
        the center of the docking box
    size : Tuple[float, float, float], default=(10., 10., 10.)
        the x-, y-, and z-radii of the docking box
    ncpu : int, default=1
        the number of cpu cores to use during the docking calculation
    extra : Optional[List[str]], default=None
        additional command line arguments that will be passed to the
        docking calculation
    prepared_ligand: Optional[S] = None,
    prepared_receptor: Optional[T] = None
    """
    software: Software
    center: Tuple[float, float, float]
    size: Tuple[float, float, float] = (10., 10., 10.)
    ncpu: int = 1
    extra: Optional[str] = None
    prepared_ligand: Optional[Path] = None,
    prepared_receptor: Optional[Path] = None

    def __post_init__(self):
        self.extra = shlex.split(self.extra) if self.extra else []
