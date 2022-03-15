from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Optional, Tuple, Union

from pyscreener.docking.metadata import SimulationMetadata
from pyscreener.docking.dock.utils import SphereMode


@dataclass(repr=True, eq=False)
class DOCKMetadata(SimulationMetadata):
    probe_radius: float = 1.4
    steric_clash_dist: float = 0.0
    min_radius: float = 1.4
    max_radius: float = 4.0
    sphere_mode: Union[SphereMode, str] = SphereMode.BOX
    docked_ligand_file: Optional[str] = None
    enclose_spheres: bool = True
    buffer: float = 10.0
    grid_params: Optional[Mapping] = None
    dock_params: Optional[Mapping] = None
    prepared_ligand: Optional[Union[str, Path]] = None
    prepared_receptor: Optional[Tuple[str, str]] = None

    def __post_init__(self):
        if isinstance(self.sphere_mode, str):
            self.sphere_mode = SphereMode.from_str(self.sphere_more)
