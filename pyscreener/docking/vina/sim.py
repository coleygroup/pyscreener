from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union

from pyscreener.docking.simulation import DockingSimulation
from pyscreener.docking.vinarunner import VinaRunner
from pyscreener.docking.vinadata import VinaCalculationData

class VinaCalculation(DockingSimulation):
    def __init__(
        self, smi: str, receptor: str,
        software: str, center: Tuple[float, float, float],
        size: Tuple[int, int, int] = (10, 10, 10), ncpu: int = 1,
        extra: Optional[str] = None, name: Optional[str] = None,
        input_file : Optional[str] = None,
        in_path: Union[str, Path] = '.',
        out_path: Union[str, Path] = '.',
        score_mode: str = 'best', k: int = 1,
        prepared_ligand: Optional[Union[str, Path]] = None,
        prepared_receptor: Optional[Union[str, Path]] = None
    ):
        self.data = VinaCalculationData(
            smi, receptor, software, center, size, ncpu, extra, name, 
            input_file, in_path, out_path, score_mode, k,
            prepared_ligand, prepared_receptor
        )
        self.runner = VinaRunner()
    
    def prepare(self):
        return super().prepare()
    
    def run(self) -> Sequence[float]:
        return super().run()
    
    def score(self) -> Optional[float]:
        return self.data.score
    
    def result(self) -> Mapping:
        return self.data.result