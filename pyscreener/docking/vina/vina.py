from functools import partial
from itertools import takewhile
from pathlib import Path
import timeit
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from tqdm import tqdm

from pyscreener.preprocessing import autobox
from pyscreener.docking import Screener
from pyscreener.docking.vina import _docking as vina_dock
from pyscreener.docking.vina import _preparation as vina_prep

class Vina(Screener):
    def __init__(self, software: str, receptors: Optional[List[str]] = None,
                 pdbids: Optional[List[str]] = None,
                 center: Optional[Tuple[float, float, float]] = None,
                 size: Tuple = (10., 10., 10.), ncpu: int = 1,
                 extra: Optional[List[str]] = None,
                 docked_ligand_file: Optional[str] = None, buffer: float = 10.,
                 score_mode: str = 'best', repeats: int = 1,
                 receptor_score_mode: str = 'best', 
                 ensemble_score_mode: str = 'best',
                 distributed: bool = False, num_workers: int = -1,
                 path: str = '.', verbose: int = 0, **kwargs):
        if software not in ('vina', 'qvina', 'smina', 'psovina'):
            raise ValueError(f'Unrecognized docking software: "{software}"')
        if center is None and docked_ligand_file is None:
            raise ValueError(
                'Args "center" and "docked_ligand_file" were None!')
        if center is None:
            print('Autoboxing ...', end=' ', flush=True)
            center, size = autobox.docked_ligand(docked_ligand_file, buffer)
            print('Done!')
            s_center = f'({center[0]:0.1f}, {center[1]:0.1f}, {center[2]:0.1f})'
            s_size = f'({size[0]:0.1f}, {size[1]:0.1f}, {size[2]:0.1f})'
            print(f'Autoboxed ligand from "{docked_ligand_file}" with', 
                  f'center={s_center} and size={s_size}', flush=True) 

        self.software = software
        # self.receptors = receptors
        self.center = center
        self.size = size
        self.ncpu = ncpu
        self.extra = extra or []

        self.repeats = repeats

        super().__init__(receptors=receptors, pdbids=pdbids,
                         score_mode=score_mode,
                         receptor_score_mode=receptor_score_mode,
                         ensemble_score_mode=ensemble_score_mode,
                         distributed=distributed,
                         num_workers=num_workers, ncpu=ncpu,
                         path=path, verbose=verbose, **kwargs)

    def __call__(self, *args, **kwargs):
        return self.dock(*args, **kwargs)

    # @property
    # def receptors(self):
    #     return self.__receptors

    # @receptors.setter
    # def receptors(self, receptors):
    #     receptors = [self.prepare_receptor(receptor) for receptor in receptors]
    #     self.__receptors = [
    #         receptor for receptor in receptors if receptor is not None]

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
        CHUNKSIZE = 1
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
                # p_in = repeat_result['in']
                # repeat_result['in'] = Path(p_in.parent.name) / p_in.name
                # p_out = repeat_result['out']
                # repeat_result['out'] = Path(p_out.parent.name) / p_out.name
                # p_log = repeat_result['log']
                # repeat_result['log'] = Path(p_log.parent.name) / p_log.name
        return ligand_results