from abc import ABC, abstractmethod, abstractproperty
from typing import Mapping, Optional, Sequence, TypeVar

S = TypeVar("S")
T = TypeVar("T")


class DockingCalculation(ABC):
    @abstractmethod
    def __init__(self):
        pass

    @abstractmethod
    def prepare(self):
        pass

    @abstractmethod
    def run(self) -> Optional[Sequence[float]]:
        pass

    @abstractproperty
    def score(self) -> Optional[float]:
        pass

    @abstractproperty
    def result(self) -> Optional[Mapping]:
        pass
