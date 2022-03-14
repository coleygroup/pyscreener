from copy import copy
from dataclasses import replace
from datetime import datetime
from itertools import chain
from pathlib import Path
import re
import shutil
import tarfile
import tempfile
from typing import Iterable, List, Optional, Sequence, Tuple, Union

import numpy as np
import ray
from tqdm import tqdm

from pyscreener.utils import ScoreMode, autobox, pdbfix, reduce_scores, run_on_all_nodes
from pyscreener.docking.sim import Simulation
from pyscreener.docking.metadata import SimulationMetadata
from pyscreener.docking.result import Result
from pyscreener.docking.runner import DockingRunner


class DockingVirtualScreen:
    def __init__(
        self,
        runner: DockingRunner,
        receptors: Optional[Iterable[str]],
        center: Optional[Tuple],
        size: Optional[Tuple],
        metadata_template: SimulationMetadata,
        pdbids: Optional[Sequence[str]] = None,
        docked_ligand_file: Optional[str] = None,
        buffer: float = 10.0,
        ncpu: int = 1,
        base_name: str = "ligand",
        path: Union[str, Path] = ".",
        score_mode: Union[ScoreMode, str] = ScoreMode.BEST,
        repeat_score_mode: Union[ScoreMode, str] = ScoreMode.BEST,
        ensemble_score_mode: Union[ScoreMode, str] = ScoreMode.BEST,
        repeats: int = 1,
        k: int = 1,
        verbose: int = 0,
    ):
        self.runner = runner
        self.runner.validate_metadata(metadata_template)

        self.center = center
        self.size = size
        self.metadata = metadata_template
        self.base_name = base_name
        self.path = path

        self.score_mode = (
            score_mode if isinstance(score_mode, ScoreMode) else ScoreMode.from_str(score_mode)
        )
        self.repeat_score_mode = (
            repeat_score_mode
            if isinstance(repeat_score_mode, ScoreMode)
            else ScoreMode.from_str(repeat_score_mode)
        )
        self.ensemble_score_mode = (
            ensemble_score_mode
            if isinstance(ensemble_score_mode, ScoreMode)
            else ScoreMode.from_str(ensemble_score_mode)
        )

        self.repeats = repeats
        self.k = k

        self.receptors = receptors or []
        if pdbids is not None:
            self.receptors = list(self.receptors)
            self.receptors.extend([pdbfix.get_pdb(pdbid, path=self.path) for pdbid in pdbids])

        if self.center is None:
            if docked_ligand_file is None:
                raise ValueError(
                    '"center" and "docked_ligand_file" are both None! Cannot compute docking box.'
                )

            self.center, size = autobox.docked_ligand(docked_ligand_file, buffer)
            self.size = self.size or size
            print(
                f'Autoboxed ligand from "{docked_ligand_file}" with',
                f"center={self.center} and size={self.size}",
                flush=True,
            )

        self.tmp_dir = tempfile.gettempdir()
        self.prepare_and_run = ray.remote(num_cpus=ncpu)(self.runner.prepare_and_run)

        self.data_templates = [
            Simulation(
                None,
                receptor,
                self.center,
                self.size,
                copy(metadata_template),
                ncpu,
                base_name,
                None,
                self.tmp_in,
                self.tmp_out,
                self.score_mode,
                k,
            )
            for receptor in self.receptors
        ]

        if not ray.is_initialized():
            try:
                ray.init("auto")
            except ConnectionError:
                ray.init()

        self.data_templates = self.prepare_receptors()

        self.planned_simulationsss = []
        self.resultsss = []

        self.num_ligands = 0
        self.total_simulations = 0

    def __len__(self):
        """the number of ligands that have been simulated. NOT the total number of simulations"""
        return self.num_ligands

    def __call__(
        self, *sources: Iterable[Union[str, Iterable[str]]], smiles: bool = True
    ) -> np.ndarray:
        """dock all of the ligands and return an array of their scores

        This function may be called with invidual ligand sources, lists of ligand sources, or a
        combination thereof, where a "source" is either an invidual SMILES string or the path
        to a chemical supply file. E.g.,

        >>> vs = DockingVirtualScreen(...)
        >>> vs('c1ccccc1', 'CCCC', 'CC(=O)')
        ...
        >>> vs(['c1ccccc1', 'CCCC', 'CC(=O)'])
        ...
        >>> vs('c1ccccc1', ['CCCC', 'CC(=O)'])
        ...
        >>> vs(['c1ccccc1', ...], 'CCCC', ['CC(=O)', ...])
        ...

        Parameters
        ----------
        *sources : Iterable[Union[str, Iterable[str]]]
            an Iterable of SMILES strings, individual chemical files, or iterables thereof of the
            ligands to dock
        smiles : bool, default=True
            whether the input ligand sources are all SMILES strigs. If false, treat the sources
            as input files

        Returns
        -------
        np.ndarray
            a vector of length `n` containing the output score for each ligand, where `n` is the
            total number of ligands that were supplied
        """
        sources = list(chain(*([s] if isinstance(s, str) else s for s in sources)))

        planned_simulationsss = self.plan(sources, smiles)
        resultsss = self.run(planned_simulationsss)

        self.resultsss.extend(resultsss)
        S = np.array(
            [
                [[r.score if r else None for r in results] for results in resultss]
                for resultss in resultsss
            ],
            dtype=float,
        )
        self.num_ligands += len(S)
        self.total_simulations += S.size

        return reduce_scores(S, self.repeat_score_mode, self.ensemble_score_mode, self.k)

    @property
    def path(self):
        """The default output path under which to collect files corresponding to simulations
        run by this `VirtualScreen`"""
        return self.__path

    @path.setter
    def path(self, path):
        path = Path(path)
        path.mkdir(parents=True, exist_ok=True)

        self.__path = path

    @property
    def tmp_dir(self) -> Path:
        """the temp directory of this `VirtualScreen`"""
        return self.__tmp_dir

    @tmp_dir.setter
    def tmp_dir(self, path: Union[str, Path]):
        timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        tmp_dir = Path(path) / "pyscreener" / f"session_{timestamp}"

        self.__tmp_dir = tmp_dir
        self.tmp_in = tmp_dir / "inputs"
        self.tmp_out = tmp_dir / "outputs"

        self.make_tmp_dirs()

    @run_on_all_nodes
    def make_tmp_dirs(self):
        for d in (self.tmp_dir, self.tmp_in, self.tmp_out):
            d.mkdir(parents=True, exist_ok=True)

    @run_on_all_nodes
    def prepare_receptors(self):
        """Prepare the receptor file(s) for each of the simulation templates"""
        return [self.runner.prepare_receptor(template) for template in self.data_templates]

    def all_results(self, flatten: bool = True) -> List[Result]:
        """A flattened list of results from all of the completed simulations"""
        if flatten:
            return list(chain(*(chain(*self.resultsss))))

        return self.resultsss

    def plan(self, sources: Iterable[str], smiles: bool = True) -> List[List[List[Simulation]]]:
        if smiles:
            planned_simulationsss = [
                [
                    [
                        replace(data_template, smi=smi, name=f"{self.base_name}_{i+len(self)}_{j}")
                        for j in range(self.repeats)
                    ]
                    for data_template in self.data_templates
                ]
                for i, smi in enumerate(sources)
            ]
        else:
            planned_simulationsss = [
                [
                    [
                        replace(
                            data_template,
                            input_file=filepath,
                            name=f"{self.base_name}_{i+len(self)}_{j}",
                        )
                        for j in range(self.repeats)
                    ]
                    for data_template in self.data_templates
                ]
                for i, filepath in enumerate(sources)
            ]

        return planned_simulationsss

    def run(self, simulationsss: List[List[List[Simulation]]]) -> List[List[List[Result]]]:
        refsss = [
            [[self.prepare_and_run.remote(s) for s in sims] for sims in simss]
            for simss in simulationsss
        ]
        return [
            [ray.get(refs) for refs in refss]
            for refss in tqdm(refsss, desc="Docking", unit="ligand", smoothing=0.0)
        ]

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

        output_id = re.sub(r"[:,.]", "", ray.state.current_node_id())
        tmp_tar = (self.tmp_dir / output_id).with_suffix(".tar.gz")

        with tarfile.open(tmp_tar, "w:gz") as tar:
            tar.add(self.tmp_in, arcname="inputs")
            tar.add(self.tmp_out, arcname="outputs")

        shutil.copy(str(tmp_tar), str(out_path))
