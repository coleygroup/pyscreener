from functools import partial
from pathlib import Path
import timeit
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from tqdm import tqdm

from pyscreener.docking import Screener
from pyscreener.docking.vina import _docking as vina_dock
from pyscreener.docking.vina import _preparation as vina_prep

class Vina(Screener):
    def __init__(self, software: str, receptors: List[str], 
                 center: Tuple, size: Tuple = (10., 10., 10.), ncpu: int = 1,
                 extra: Optional[List[str]] = None,
                 score_mode: str = 'best', repeats: int = 1,
                 receptor_score_mode: str = 'best', 
                 ensemble_score_mode: str = 'best',
                 distributed: bool = False, num_workers: int = -1,
                 path: str = '.', verbose: int = 0, **kwargs):
        if software not in ('vina', 'qvina', 'smina', 'psovina'):
            raise ValueError(f'Unrecognized docking software: "{software}"')
        
        self.software = software
        self.receptors = receptors
        self.center = center
        self.size = size
        self.ncpu = ncpu
        self.extra = extra or []

        self.score_mode = score_mode
        self.repeats = repeats

        super().__init__(receptor_score_mode=receptor_score_mode,
                         ensemble_score_mode=ensemble_score_mode,
                         distributed=distributed, num_workers=num_workers,
                         path=path, verbose=verbose, **kwargs)

    def __call__(self, *args, **kwargs):
        return self.dock(*args, **kwargs)

    @property
    def receptors(self):
        return self.__receptors

    @receptors.setter
    def receptors(self, receptors):
        receptors = [self.prepare_receptor(receptor) for receptor in receptors]
        self.__receptors = [
            receptor for receptor in receptors if receptor is not None]

    def prepare_receptor(self, receptor):
        return vina_prep.prepare_receptor(receptor)

    @staticmethod
    def prepare_from_smi(smi: str, name: Optional[str] = None,
                         path: str = '.') -> Tuple[str, str]:
        return vina_prep.prepare_from_smi(smi, name, path)

    @staticmethod
    def prepare_from_file(filepath: str, use_3d: bool = False,
                          name: Optional[str] = None,
                          path: str = '.') -> Tuple[str, str]:
        return vina_prep.prepare_from_file(filepath, use_3d, name, path)
    
    def run_docking(self, ligands: Sequence[Tuple[str, str]]
                   ) -> List[List[List[Dict]]]:
        dock_ligand = partial(
            vina_dock.dock_ligand,
            software=self.software, receptors=self.receptors,
            center=self.center, size=self.size, ncpu=self.ncpu,
            extra=self.extra, path=self.out_path, repeats=self.repeats
        )
        CHUNKSIZE = 32
        with self.Pool(self.distributed, self.num_workers, self.ncpu) as pool:
            ligs_recs_reps = pool.map(dock_ligand, ligands, 
                                      chunksize=CHUNKSIZE)
            ligs_recs_reps = list(tqdm(ligs_recs_reps, total=len(ligands),
                                       desc='Docking', unit='ligand'))

        return ligs_recs_reps

    @staticmethod
    def parse_ligand_results(ligand_results: List[List[Dict]],
                             score_mode: str = 'best'):
        for receptor_results in ligand_results:
            for repeat_result in receptor_results:
                # vina-type log files have scoring information between this 
                # table border and the line: "Writing output ... done."
                TABLE_BORDER = '-----+------------+----------+----------'
                try:
                    with open(repeat_result['log']) as fid:
                        for line in fid:
                            if TABLE_BORDER in line:
                                break

                        score_lines = takewhile(
                            lambda line: 'Writing' not in line, fid)
                        scores = [float(line.split()[1])
                                  for line in score_lines]

                    if len(scores) == 0:
                        score = None
                    else:
                        score = Screener.calc_score(scores, score_mode)
                except OSError:
                    score = None
                repeat_result['score'] = score
