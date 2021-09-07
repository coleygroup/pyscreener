from dataclasses import dataclass
from pathlib import Path
import shlex
from typing import Optional, Union

from pyscreener.docking.metadata import CalculationMetadata
from pyscreener.docking.vina.utils import Software

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
    prepared_ligand: Optional[Union[str, Path]] = None,
    prepared_receptor: Optional[Union[str, Path]] = None
    """
    software: Software
    extra: Optional[str] = None
    prepared_ligand: Optional[Union[str, Path]] = None
    prepared_receptor: Optional[Union[str, Path]] = None

    def __post_init__(self):
        self.extra = shlex.split(self.extra) if self.extra else []
