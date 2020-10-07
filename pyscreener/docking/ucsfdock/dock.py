from functools import partial
from pathlib import Path
import timeit
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from tqdm import tqdm

from pyscreener.docking import Screener
from pyscreener.docking.ucsfdock import _docking as ucsfdock_dock
from pyscreener.docking.ucsfdock import _preparation as ucsfdock_prep

class DOCK(Screener):
    def __init__(self, receptors: List[str], 
                 center: Tuple[float, float, float],
                 size: Tuple[float, float, float] = (20., 20., 20.), 
                 docked_ligand_file: Optional[str] = None,
                 use_largest: bool = False, buffer: float = 10.,
                 enclose_spheres: bool = True,
                 repeats: int = 1, score_mode: str = 'best',
                 receptor_score_mode: str = 'best', 
                 ensemble_score_mode: str = 'best',
                 distributed: bool = False, num_workers: int = -1,
                 path: str = '.', verbose: int = 0, **kwargs):
        
        super().__init__(score_mode=score_mode, 
                         receptor_score_mode=receptor_score_mode,
                         ensemble_score_mode=ensemble_score_mode,
                         distributed=distributed, num_workers=num_workers,
                         path=path, verbose=verbose, **kwargs)
    
        self.center = center
        self.size = size
        self.docked_ligand_file = docked_ligand_file
        self.use_largest = use_largest
        self.buffer = buffer
        self.enclose_spheres = enclose_spheres
        self.receptors = receptors
        self.repeats = repeats

    def __call__(self, *args, **kwargs):
        return self.dock(*args, **kwargs)

    @property
    def receptors(self):
        return self.__receptors

    @receptors.setter
    def receptors(self, receptors):
        receptors = [self.prepare_receptor(receptor) for receptor in receptors]
        receptors = [receptor for receptor in receptors if receptor is not None]
        if len(receptors) == 0:
            raise RuntimeError('Preparation failed for all receptors!')
        self.__receptors = receptors
        
    def prepare_receptor(self, receptor: str) -> Optional[Tuple[str, str]]:
        return ucsfdock_prep.prepare_receptor(
            receptor, center=self.center, size=self.size,
            docked_ligand_file=self.docked_ligand_file,
            use_largest=self.use_largest, buffer=self.buffer,
            enclose_spheres=self.enclose_spheres, path=self.in_path
        )

    @staticmethod
    def prepare_from_smi(smi: str, name: str = 'ligand',
                         path: str = '.') -> Tuple[str, str]:
        return ucsfdock_prep.prepare_from_smi(smi, name, path)

    @staticmethod
    def prepare_from_file(filepath: str, use_3d: bool = False,
                          name: Optional[str] = None, path: str = '.'):
        return ucsfdock_prep.prepare_from_file(filepath, use_3d, name, path)

    def run_docking(self, ligands: Sequence[Tuple[str, str]]
                   ) -> List[List[List[Dict]]]:
        dock_ligand = partial(
            ucsfdock_dock.dock_ligand,
            receptors=self.receptors, path=self.path, repeats=self.repeats
        )
        CHUNKSIZE = 1
        with self.Pool(self.distributed, self.num_workers) as pool:
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
                try:
                    scores = []
                    with open(repeat_result['out']) as fid:
                        for line in fid:
                            if 'Grid_Score:' in line:
                                try:
                                    scores.append(float(line.split()[2]))
                                except:
                                    continue

                    if len(scores) == 0:
                        score = None
                    else:
                        score = Screener.calc_score(scores, score_mode)
                except OSError:
                    score = None

                repeat_result['score'] = score
                p_in = repeat_result['in']
                repeat_result['in'] = Path(p_in.parent.name) / p_in.name
                p_out = repeat_result['out']
                repeat_result['out'] = Path(p_out.parent.name) / p_out.name
                p_log = repeat_result['log']
                repeat_result['log'] = Path(p_log.parent.name) / p_log.name
        return ligand_results
