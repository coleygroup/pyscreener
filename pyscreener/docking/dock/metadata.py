from dataclasses import dataclass
from enum import Enum, auto
from typing import Optional, Tuple

from pyscreener.docking.metadata import CalculationMetadata
from pyscreener.docking.dock.utils import SphereMode

class SphereMode(Enum):
    BOX = auto()
    LARGEST = auto()
    LIGAND = auto()
    
@dataclass(repr=True, eq=False)
class DOCKMetadata(CalculationMetadata):
    """

    Attributes
    ---------
    center : Optional[Tuple[float, float, float]]
        the center of the docking box (if known)
    size : Tuple[float, float, float]
        the x-, y-, and z-radii of the docking box

    Parmeters
    ---------
    center : Optional[Tuple[float, float, float]], default=None
    size : Tuple[float, float, float], default=(10., 10., 10.)
    """
    center: Optional[Tuple[float, float, float]] = None
    size: Tuple[float, float, float] = (10., 10., 10.)
    probe_radius: float = 1.4
    steric_clash_dist: float = 0.0,
    min_radius: float = 1.4
    max_radius: float = 4.0,
    sphere_mode: SphereMode = SphereMode.LARGEST
    docked_ligand_file: Optional[str] = None
    enclose_spheres: bool = True
    buffer: float = 10.