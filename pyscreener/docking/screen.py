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

from pyscreener.utils import Reduction, autobox, pdbfix, reduce_scores, run_on_all_nodes
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
        reduction: Union[Reduction, str] = Reduction.BEST,
        receptor_reduction: Union[Reduction, str] = Reduction.BEST,
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
            reduction if isinstance(reduction, Reduction) else Reduction.from_str(reduction)
        )
        self.receptor_reduction = (
            receptor_reduction
            if isinstance(receptor_reduction, Reduction)
            else Reduction.from_str(receptor_reduction)
        )

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

        ncpu = ncpu if self.runner.is_multithreaded() else 1
        self.prepare_and_run = ray.remote(num_cpus=ncpu)(self.runner.prepare_and_run)

        self.simulation_templates = [
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

        self.simulation_templates = self.prepare_receptors()

        self.planned_simulationss = []
        self.run_simulationss = []
        self.resultss = []

        self.num_ligands = 0
        self.num_simulations = 0

    def __len__(self):
        """the number of ligands that have been simulated. NOT the total number of simulations"""
        return self.num_ligands

    def __call__(
        self,
        *sources: Iterable[Union[str, Iterable[str]]],
        smiles: bool = True,
        reduction: Optional[Reduction] = None,
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

        NOTE: this function is largely for convencience and pipelines setup() -> run() -> reduce()
        together. Values of `nan` in the returned array can indicate either an invalid ligand
        (could not be parsed by rdkit) or a failed simulation. If the distinction is meaningful to
        you, then you should manually compare the List[List[Result]] from run() to the array from
        reduce()

        Parameters
        ----------
        *sources : Iterable[Union[str, Iterable[str]]]
            an Iterable of SMILES strings, individual chemical files, or iterables thereof of the
            ligands to dock
        smiles : bool, default=True
            whether the input ligand sources are all SMILES strigs. If false, treat the sources
            as input files
        reduction : Optional[Reduction], default=None
            the reduction to apply to multiple receptor scores for the same ligand. If None, use
            self.ensemble_reduction

        Returns
        -------
        np.ndarray
            a vector of length `n` where each entry is the docking score of the `i`th ligand and `n`
            total number of ligands that were supplied. If multiple receptors were supplied, this is
            either (1) an `n x r` array if `reduction` is None or (2) a length `n` vector after
            applying the given reduction.
        """
        sources = list(chain(*([s] if isinstance(s, str) else s for s in sources)))

        simss = self.setup(sources, smiles)
        resultss = self.run(simss)

        return self.reduce(resultss, reduction)

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
        return [self.runner.prepare_receptor(template) for template in self.simulation_templates]

    def results(self) -> List[Result]:
        """A flattened list of results from all of the completed simulations"""
        return list(chain(*self.resultss))

    def simulations(self) -> List[Simulation]:
        """A flattened list of simulations from all of the completed simulations"""
        return list(chain(*self.run_simulationss))

    def setup(self, sources: Iterable[str], smiles: bool = True) -> List[List[Simulation]]:
        """Set up the simulations for the input sources

        Parameters
        ----------
        sources : Iterable[str]
            an iterable of SMILES strings or filepaths
        smiles : bool, optional
            whether the inputs are SMIELS strings, by default True

        Returns
        -------
        List[List[Simulation]]
            the Simulation objects corresponding to the input ligands
        """
        if smiles:
            simss = [
                [
                    replace(template, smi=smi, name=f"{self.base_name}_{i+len(self)}")
                    for template in self.simulation_templates
                ]
                for i, smi in enumerate(sources)
            ]
        else:
            simss = [
                [
                    replace(template, input_file=filepath, name=f"{self.base_name}_{i+len(self)}")
                    for template in self.simulation_templates
                ]
                for i, filepath in enumerate(sources)
            ]

        return simss

    def run(self, simulationss: List[List[Simulation]]) -> List[List[Result]]:
        refss = [[self.prepare_and_run.remote(s) for s in sims] for sims in simulationss]

        resultss = [
            ray.get(refs) for refs in tqdm(refss, desc="Docking", unit="ligand", smoothing=0.0)
        ]

        self.run_simulationss.extend(simulationss)
        self.resultss.extend(resultss)
        self.num_ligands += len(resultss)
        self.num_simulations += len(list(chain(*resultss)))

        return resultss

    def reduce(self, resultss: List[List[Result]], reduction: Optional[Reduction]):
        reduction = reduction or self.receptor_reduction

        S = np.array(
            [[r.score if r else None for r in results] for results in resultss], dtype=float
        )

        if S.shape[1] == 1:
            return S.flatten()

        if reduction is None:
            return S

        return reduce_scores(S, reduction, self.k)

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
