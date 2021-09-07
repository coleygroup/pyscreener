from copy import copy
from dataclasses import replace
from datetime import datetime
from pathlib import Path
import re
import shutil
import tarfile
import tempfile
from typing import Iterable, List, Optional, Tuple, Union

import ray

# from pyscreener.base import VirtualScreen
from pyscreener.utils import ScoreMode
from pyscreener.docking.data import CalculationData
from pyscreener.docking.metadata import CalculationMetadata
from pyscreener.docking.runner import DockingRunner
from pyscreener.docking.utils import run_on_all_nodes

class DockingVirtualScreen:
    def __init__(
        self, runner: DockingRunner,
        receptors: Iterable[str],
        center: Optional[Tuple],
        size: Optional[Tuple],
        metadata_template: CalculationMetadata,
        ncpu: int = 1,
        base_name: str = 'ligand',
        path: Union[str, Path] = '.',
        score_mode : ScoreMode = ScoreMode.BEST,
        repeat_score_mode : ScoreMode = ScoreMode.BEST,
        ensemble_score_mode : ScoreMode = ScoreMode.BEST,
        k : int = 1,
    ):
        # super().__init__()

        self.runner = runner
        self.receptors = receptors
        self.center = center
        self.size = size
        self.metadata = metadata_template
        self.base_name = base_name
        self.path = path
        self.score_mode = score_mode
        self.repeat_score_mode = repeat_score_mode
        self.ensemble_score_mode = ensemble_score_mode
        self.k = k
        
        self.tmp_dir = tempfile.gettempdir()
        self.prepare_and_run = ray.remote(num_cpus=ncpu)(
            self.runner.prepare_and_run
        )
        self.data_templates = [CalculationData(
            None, receptor, center, size, copy(metadata_template),
            ncpu, base_name, None, self.tmp_in, self.tmp_out, score_mode, k
        ) for receptor in receptors]

        if not ray.is_initialized():
            try:
                ray.init('auto')
            except ConnectionError:
                ray.init()
    
        self.data_templates = self.prepare_receptors()

        self.planned_simulationss = []
        self.completed_simulationss = []
        self.num_simulations = 0

    def __len__(self):
        return self.num_simulations

    @run_on_all_nodes
    def prepare_receptors(self):
        return [
            self.runner.prepare_receptor(template)
            for template in self.data_templates
        ]

    @property
    def path(self):
        return self.__path

    @path.setter
    def path(self, path):
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)

        self.__path = path

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
        planned_simulationss = [[
            replace(
                data_template, smi=smi, name=f'{self.base_name}_{i+len(self)}'
            ) for i, smi in enumerate(smis)
        ] for data_template in self.data_templates]
        refss = [
            [self.prepare_and_run.remote(s) for s in simulations]
            for simulations in planned_simulationss
        ]
        completed_simulationss = [ray.get(refs) for refs in refss]

        self.num_simulations += sum(
            sum(1 for _ in simulations)
            for simulations in completed_simulationss
        )
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
        completed_simulationss = [
            [self.prepare_and_run(s) for s in simulations]
            for simulations in self.planned_simulationss
        ]

        self.num_simulations += sum(
            sum(1 for _ in simulations)
            for simulations in completed_simulationss
        )
        self.completed_simulationss.extend(completed_simulationss)
        
        return completed_simulationss
    
    @run_on_all_nodes
    def collect_files(self, path: Optional[Union[str, Path]] = None):
        """Collect all the files from the local disks of the respective nodes

        For I/O purposes, input and output files for each simulation are 
        created on the local disk of each node. If these files are desired at 
        the end, they must be copied over from the node's local file system to 
        the final destination.
        
        This is achieved by creating a gzipped tar file of the temp directory 
        (the one that contains all of the input and output files for 
        simulations conducted on that node) and moving these tar files under 
        the desired path. Each tar file is named according the node ID from 
        which it originates.

        This function should ideally only be called once during the lifetime
        of a Screener because it is slow and early calls will yield nothing
        over a single, final call.

        Parameters
        ----------
        out_path : Optional[Union[str, Path]], default=None
            the path under which the tar files should be collected to.
            If None, use self.path
        """
        out_path = Path(path or self.path)
        out_path.mkdir(parents=True, exist_ok=True)

        output_id = re.sub(r'[:,.]', '', ray.state.current_node_id())
        tmp_tar = (self.tmp_dir / output_id).with_suffix('.tar.gz')
        with tarfile.open(tmp_tar, 'w:gz') as tar:
            tar.add(self.tmp_in, arcname='inputs')
            tar.add(self.tmp_out, arcname='outputs')
        shutil.copy(str(tmp_tar), str(out_path))
