from abc import ABC, abstractmethod
from typing import Optional, Sequence

from pyscreener.docking.data import CalculationData

class DockingRunner(ABC):
    @staticmethod
    @abstractmethod
    def prepare(data: CalculationData) -> CalculationData:
        pass

    @staticmethod
    @abstractmethod
    def run(data: CalculationData) -> Optional[Sequence[float]]:
        pass