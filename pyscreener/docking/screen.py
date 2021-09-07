from dataclasses import replace
from datetime import datetime
from enum import auto, Enum
from pathlib import Path
import tempfile
from typing import Iterable, List, Optional, Tuple, Union

import ray

from pyscreener.base import VirtualScreen
from pyscreener.utils import ScoreMode
from pyscreener.docking.data import CalculationData
from pyscreener.docking.metadata import CalculationMetadata
# from pyscreener.docking.utils import ScreenType
# from pyscreener.docking.dock.runner import DOCKRunner
from pyscreener.docking.vina.runner import VinaRunner

class ScreenType(Enum):
    VINA = auto()
    DOCK = auto()

class DockingVirtualScreen(VirtualScreen):
    def __init__(
        self, receptors: Iterable[str],
        center: Optional[Tuple],
        size: Optional[Tuple],
        metadata_template: CalculationMetadata,
        ncpu: int = 1,
        base_name: str = 'ligand',
        in_path: Union[str, Path] = '.',
        out_path: Union[str, Path] = '.',
        score_mode : ScoreMode = ScoreMode.BEST,
        repeat_score_mode : ScoreMode = ScoreMode.BEST,
        ensemble_score_mode : ScoreMode = ScoreMode.BEST,
        k : int = 1,
        screen_type: ScreenType = ScreenType.VINA
    ):
        # super().__init__()

        self.receptors = receptors
        self.center = center
        self.size = size
        self.metadata = metadata_template
        self.in_path = Path(in_path)
        self.out_path = Path(out_path)
        self.score_mode = score_mode
        self.repeat_score_mode = repeat_score_mode
        self.ensemble_score_mode = ensemble_score_mode
        self.k = k

        self.runner ={
            ScreenType.VINA: VinaRunner,
            # ScreenType.DOCK: DOCKRunner
        }[screen_type]
        
        self.tmp_dir = tempfile.gettempdir()
        self.prepare_and_run = ray.remote(num_cpus=ncpu)(
            self.runner.prepare_and_run
        )
        self.data_templates = [CalculationData(
            '', receptor, center, size, metadata_template, ncpu, base_name,
            None, self.tmp_in, self.tmp_out, score_mode, self.k
        ) for receptor in receptors]

        if not ray.is_initialized():
            try:
                ray.init('auto')
            except ConnectionError:
                ray.init()
    
        refs = []
        for n in ray.nodes():
            @ray.remote(resources={f'node:{n["NodeManagerAddress"]}': 0.01})
            def prepare_templates(runner, templates):
                templates = [
                    runner.prepare_receptor(template) for template in templates
                ]
                return templates

            refs.append(prepare_templates.remote(
                self.runner, self.data_templates
            ))
        self.data_templates = ray.get(refs[-1])

        self.planned_simulationss = []
        self.completed_simulationss = []

    @property
    def in_path(self):
        return self.__in_path

    @in_path.setter
    def in_path(self, path):
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)

        self.__in_path = path
    
    @property
    def out_path(self):
        return self.__out_path
    
    @out_path.setter
    def out_path(self, path):
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)

        self.__out_path = path

    @property
    def tmp_dir(self) -> Path:
        """the Screener's temp directory"""
        return self.__tmp_dir
        
    @tmp_dir.setter
    def tmp_dir(self, path: Union[str, Path]):
        timestamp = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        tmp_dir = Path(path) / 'pyscreener' / f'session_{timestamp}'
        tmp_dir.mkdir(exist_ok=True, parents=True)
        self.__tmp_dir = tmp_dir
        self.tmp_in = tmp_dir / 'inputs'
        self.tmp_out = tmp_dir / 'outputs'

    @property
    def tmp_in(self) -> Path:
        return self.__tmp_in

    @tmp_in.setter
    def tmp_in(self, path: Union[str, Path]):
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        self.__tmp_in = path
    
    @property
    def tmp_out(self) -> Path:
        return self.__tmp_out

    @tmp_out.setter
    def tmp_out(self, path: Union[str, Path]):
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)
        self.__tmp_out = path
    
    def __call__(self, *smis: Iterable[str]) -> List[float]:
        planned_simulationss = [
            [replace(data_template, smi=smi) for smi in smis]
            for data_template in self.data_templates
        ]
        refss = [
            [self.prepare_and_run.remote(s) for s in simulations]
            for simulations in planned_simulationss
        ]
        completed_simulationss = [ray.get(refs) for refs in refss]
        self.completed_simulationss.extend(completed_simulationss)

        return completed_simulationss

    def plan(
        self, smis: Iterable[str]
    ) -> List[List[CalculationData]]:
        planned_simulationss = [
            [replace(data_template, smi=smi) for smi in smis]
            for data_template in self.data_templates
        ]
        self.planned_simulationss.extend(planned_simulationss)

        return planned_simulationss

    def prepare(self):
        pass

    def run(self) -> List[List[CalculationData]]:
        completed_simulations = [
            [self.prepare_and_run(s) for s in simulations]
            for simulations in self.planned_simulationss
        ]
        return completed_simulations
    