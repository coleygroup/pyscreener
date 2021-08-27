from dataclasses import dataclass
from enum import Enum, auto
from pathlib import Path
from pyscreener.docking.dock.metadata import DOCKMetadata
from typing import Any, Mapping, Optional, Tuple, TypeVar, Union

from pyscreener.exceptions import InvalidResultError, NotSimulatedError
from pyscreener.utils import ScoreMode
from pyscreener.docking.vina.metadata import VinaMetadata
from pyscreener.docking.metadata import CalculationMetadata

S = TypeVar('S')
T = TypeVar('T')

class CalculationType(Enum):
    def _generate_next_value_(name, start, count, last_values):
        return name.upper()

    VINA = auto()
    DOCK = auto()

@dataclass(repr=True, eq=False)
class CalculationData:
    """

    Attributes
    ---------
    smi : str
        the SMILES string of the ligand that will be docked
    receptor : Optional[str]
        the filepath of a receptor to prepare for docking
    metadata : CalculationMetadata
        the parameters with which to prepare and run the simulation
    name : str
        the name to use when creating the ligand input file and output files
    input_file : Optional[str]
    in_path: Union[str, Path]
        the path under which input will be placed
    out_path: Union[str, Path]
        the path under which output will be placed
    score_mode : str
        the mode used to calculate a score for an individual docking
        calculation given multiple output scored conformations
    k : int
        the number of top scores to use if calculating an average
    prepared_ligand: Optional[Union[str, Path]]
    prepared_receptor: Optional[Union[str, Path]]
    result : Optional[Mapping]
        the result of the docking calculation. None if the calculation has not
        been performed yet.

    Parmeters
    ---------
    smi : str
    receptor : Optional[str], default=None
    metadata : CalculationMetadata
    name : str, default='ligand'
    input_file : Optional[str], default=None
    in_path: Union[str, Path], default='.'
    out_path: Union[str, Path], default='.'
    score_mode : str, default=ScoreMode.BEST
    k : int, default=1
    prepared_ligand: Optional[Union[str, Path]], default=None
    prepared_receptor: Optional[Union[str, Path]], default=None
    result : Optional[Mapping], default=None
    """
    smi: str
    receptor: str
    metadata: CalculationMetadata
    name: str = 'ligand'
    input_file : Optional[str] = None
    in_path: Union[str, Path] = '.'
    out_path: Union[str, Path] = '.'
    score_mode : ScoreMode = ScoreMode.BEST
    k : int = 1
    prepared_ligand: Optional[S] = None,
    prepared_receptor: Optional[T] = None
    result : Optional[Mapping] = None

    def __post_init__(self):
        self.in_path = Path(self.in_path)
        self.out_path = Path(self.out_path)
    
    @property
    def calculation_type(self) -> CalculationType:
        if isinstance(self.metadata, VinaMetadata):
            return CalculationType.VINA
        elif isinstance(self.metadata, DOCKMetadata):
            return CalculationType.DOCK

        raise RuntimeError("fell through!")

    @property
    def score(self) -> Optional[float]:
        """the docking score of this calculation
        
        Raises
        ------
        NotSimulatedError
            if this calculation has not been run yet
        """
        if self.result is None:
            raise NotSimulatedError(
                'Simulation has not been run!'
            )

        try:
            return self.result['score']
        except KeyError:
            raise InvalidResultError(
                'No key: "score" in result dictionary (self.result)!'
            )