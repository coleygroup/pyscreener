from dataclasses import dataclass
from pathlib import Path
import shlex
from typing import Mapping, Optional, Tuple, Union

from pyscreener.exceptions import InvalidResultException, NotSimulatedException

@dataclass(repr=True, eq=False)
class VinaCalculationData:
    """

    Attributes
    ---------
    smi : str
        the SMILES string of the ligand that will be docked
    receptor : Optional[str]
        the filepath of a receptor to prepare for docking
    software : str
        the software that will be used
    pdbids : Optiona[str]
        a PDB ID corresponding to a receptor to prepare for docking
    center : Tuple[float, float, float]
        the center of the docking box
    size : Tuple[float, float, float], default=(10., 10., 10.)
        the x-, y-, and z-radii of the docking box
    ncpu : int, default=1
        the number of cpu cores to use during the docking calculation
    extra : Optional[List[str]]
        additional command line arguments that will be passed to the
        docking calculation
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
        

    Parmeters
    ---------
    smi : str
        the SMILES string of the ligand that will be docked
    receptor : Optional[str], default=None
        the filepath of a receptor to prepare for docking
    software : str
        the software that will be used
    pdbids : Optiona[List[str]], default=None
        a PDB ID corresponding to a receptor to prepare for docking
    center : Tuple[float, float, float], default=None
        the center of the docking box
    size : Tuple[float, float, float], default=(10., 10., 10.)
        the x-, y-, and z-radii of the docking box
    ncpu : int, default=1
        the number of cpu cores to use during the docking calculation
    extra : Optional[List[str]], default=None
        additional command line arguments that will be passed to the
        docking calculation
    name : str, default='ligand'
    input_file : Optional[str]
    in_path: Union[str, Path], default='.'
        the path under which input will be placed
    out_path: Union[str, Path], default='.'
        the path under which output will be placed
    score_mode : str, default='best'
        the mode used to calculate a score for an individual docking
        calculation given multiple output scored conformations
    k : int, default=1
        the number of top scores to use if calculating an average
    prepared_ligand: Optional[Union[str, Path]], default=None
    prepared_receptor: Optional[Union[str, Path]], default=None
    result : Optional[Mapping], default=None
    """
    smi: str
    receptor: str
    software: str
    center: Tuple[float, float, float]
    size: Tuple[float, float, float] = (10., 10., 10.)
    ncpu: int = 1
    extra: Optional[str] = None
    name: str = 'ligand'
    input_file : Optional[str] = None
    in_path: Union[str, Path] = '.'
    out_path: Union[str, Path] = '.'
    score_mode : str = 'best'
    k : int = 1
    prepared_ligand: Optional[Union[str, Path]] = None,
    prepared_receptor: Optional[Union[str, Path]] = None
    result : Optional[Mapping] = None

    def __post_init__(self):
        if self.software not in ('vina', 'smina', 'psovina', 'qvina'):
            raise ValueError(f'Invalid docking software: "{self.software}"')

        self.extra = shlex.split(self.extra) if self.extra else []
        self.in_path = Path(self.in_path)
        self.out_path = Path(self.out_path)
    
    @property
    def score(self) -> Optional[float]:
        """the docking score of this calculation
        
        Raises
        ------
        NotSimulatedError
            if this calculation has not been run yet
        """
        if self.result is None:
            raise NotSimulatedException(
                'Simulation has not been run!'
            )

        try:
            return self.result['score']
        except KeyError:
            raise InvalidResultException(
                'Invalid result mapping was set. No "score" key detected!'
            )