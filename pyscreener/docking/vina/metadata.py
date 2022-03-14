from dataclasses import dataclass
from pathlib import Path
import shlex
from typing import Iterable, Optional, Union

from pyscreener.exceptions import UnsupportedSoftwareError
from pyscreener.docking.metadata import SimulationMetadata
from pyscreener.docking.vina.utils import Software


@dataclass(repr=True, eq=False)
class VinaMetadata(SimulationMetadata):
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
    exhaustiveness: int
        the exhaustiveness of the global search. Larger values are more exhaustive
    num_modes: int
        the number of output modes
    energy_range: float
        the maximum energy difference (in kcal/mol) between the best and worst output binding modes
    extra : List[str]
        additional arguments that will be passed to the docking calculation
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
    exhaustiveness: int, default=8
    num_modes: int, default=9
    energy_range: float, default=3.
    extra : str, default=""
        a string containing the additional command line arguments to pass to a run of a vina-type
        software for options not contained within the default metadata. E.g. for a run of Smina, extra="--force_cap ARG" or for PSOVina, extra="-w ARG"
    prepared_ligand: Optional[Union[str, Path]] = None,
    prepared_receptor: Optional[Union[str, Path]] = None
    """

    software: Union[Software, str] = Software.VINA
    exhaustiveness: int = 8
    num_modes: int = 9
    energy_range: float = 3.0
    extra: Union[str, Iterable[str]] = ""
    prepared_ligand: Optional[Union[str, Path]] = None
    prepared_receptor: Optional[Union[str, Path]] = None

    def __post_init__(self):
        if isinstance(self.software, str):
            try:
                self.software = Software.from_str(self.software)
            except KeyError:
                raise UnsupportedSoftwareError(
                    f'"{self.software}" is not a supported software for vina-type screens!'
                )

        self.extra = shlex.split(self.extra) if isinstance(self.extra, str) else self.extra
