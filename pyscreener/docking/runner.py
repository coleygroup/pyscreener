from abc import ABC, abstractmethod
from typing import Optional, Sequence

from pyscreener.docking.sim import Simulation
from pyscreener.docking.metadata import SimulationMetadata


class DockingRunner(ABC):
    @staticmethod
    @abstractmethod
    def prepare_receptor(data: Simulation) -> Simulation:
        """Prepare the receptor file(s) for the given simulation"""

    @staticmethod
    @abstractmethod
    def prepare_ligand(data: Simulation) -> Simulation:
        """Prepare the ligand file(s) for the given simulation"""

    @staticmethod
    @abstractmethod
    def run(data: Simulation) -> Optional[Sequence[float]]:
        """Run the given simulation and return the score(s) of the docked conformers"""

    @staticmethod
    @abstractmethod
    def prepare_and_run(data: Simulation) -> Simulation:
        """Prepare the receptor and ligand files then run the given simulation. Roughly equivlaent
        to `prepare_*()` followed by `run()` but doesn't return the the scores of the conformers"""

    @staticmethod
    def validate_metadata(metadata: SimulationMetadata):
        """Validate the metadata of the simulation. E.g., ensure that the specified software is
        installed for Vina-type screens."""
