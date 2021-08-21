from dataclasses import dataclass
from pathlib import Path
import shlex
from typing import Optional, Tuple, Union

@dataclass(repr=False, eq=False)
class VinaCalculationData:
    ligand: str
    receptor: str
    software: str
    center: Tuple[float, float, float]
    size: Tuple[float, float, float] = (10., 10., 10.)
    ncpu: int = 1
    extra: Optional[str] = None
    name: str = None
    in_path: Union[str, Path] = '.'
    out_path: Union[str, Path] = '.'
    prepared_ligand: Optional[Union[str, Path]] = None,
    prepared_receptor: Optional[Union[str, Path]] = None

    def __post_init__(self):
        # if self.software not in ('vina', 'smina', 'psovina', 'qvina'):
        #     raise ValueError(f'Invalid docking software: "{self.software}"')
        self.extra = shlex.split(self.extra) if self.extra else []
        self.in_path = Path(self.in_path)
        self.out_path = Path(self.out_path)