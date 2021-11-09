from dataclasses import dataclass
from os import read
from pathlib import Path
from typing import Optional, Tuple, Union
from openbabel import pybel

from pyscreener.exceptions import InvalidResultError, NotSimulatedError
from pyscreener.utils import ScoreMode
from pyscreener.docking.metadata import CalculationMetadata
from pyscreener.docking.result import Result


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
        the filepath of an arbitrary molecular input file containing a single molecule
    in_path : Union[str, Path]
        the path under which input will be placed
    out_path: Union[str, Path]
        the path under which output will be placed
    score_mode : str
        the mode used to calculate a score for an individual docking
        calculation given multiple output scored conformations
    k : int
        the number of top scores to use if calculating an average
    prepared_ligand : Optional[Union[str, Path]]
    prepared_receptor : Optional[Union[str, Path]]
    result : Optional[Mapping]
        the result of the docking calculation. None if the calculation has not
        been performed yet.
    input_file_bytes : None 
        the file bytes for distributed computing. None if there is 
        no file being sent. 

    Parmeters
    ---------
    smi : str
    receptor : Optional[str], default=None
    metadata : CalculationMetadata
    name : str, default='ligand'
    input_file : Optional[str], default=None
    in_path : Union[str, Path], default='.'
    out_path : Union[str, Path], default='.'
    score_mode : str, default=ScoreMode.BEST
    k : int, default=1
    prepared_ligand : Optional[Union[str, Path]], default=None
    prepared_receptor : Optional[Union[str, Path]], default=None
    result : Optional[Result], default=None
    input_file_bytes : Optional[str], default = None
    """

    smi: str
    receptor: str
    center: Tuple[float, float, float]
    size: Tuple[float, float, float]
    metadata: CalculationMetadata
    ncpu: int = 1
    name: str = "ligand"
    input_file: Optional[str] = None
    in_path: Union[str, Path] = "."
    out_path: Union[str, Path] = "."
    score_mode: ScoreMode = ScoreMode.BEST
    k: int = 1
    result: Optional[Result] = None
    input_file_bytes: Optional[str] = None

    def __post_init__(self):
        self.in_path = Path(self.in_path)
        self.out_path = Path(self.out_path)
        if self.input_file is not None:
            fmt = Path(self.input_file).suffix.strip(".")
            mols = list(pybel.readfile(fmt, self.input_file))
            self.input_file_bytes = mols[0]

    @property
    def score(self) -> Optional[float]:
        """the docking score of this calculation

        Raises
        ------
        NotSimulatedError
            if this calculation has not been run yet
        """
        if self.result is None:
            raise NotSimulatedError("Simulation has not been run!")

        try:
            return self.result.score
        except AttributeError:
            raise InvalidResultError(
                'No attribute: "score" in Result object (self.result)!'
            )
