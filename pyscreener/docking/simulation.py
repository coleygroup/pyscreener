from abc import ABC, abstractmethod, abstractproperty
from pathlib import Path
from typing import Mapping, Sequence, TypeVar, Union

S = TypeVar('S')
T = TypeVar('T')

class DockingSimulation(ABC):
    def __init__(self, path: Union[str, Path]):
        self.path = path
    
    @property
    def path(self) -> Path:
        return self.__path
    
    @path.setter
    def path(self, path):
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        
        self.__path = path

    @abstractproperty
    def receptor(self) -> S:
        pass
    
    @abstractproperty
    def smi(self) -> str:
        pass

    @abstractproperty
    def ligand(self) -> T:
        pass

    @abstractproperty
    def score(self) -> float:
        pass

    @abstractproperty
    def result(self) -> Mapping:
        pass

    @abstractmethod
    def prepare(self):
        pass

    @abstractmethod
    def run(self) -> Sequence[float]:
        pass