from pathlib import Path
from typing import Dict, Iterable, List

from ..docking import run_and_parse_docker

def build_argv():
    pass

def parse_log():
    pass

def dock_ligand(ligand, path: str = '.', repeats: int = 1,
                score_mode: str = 'best') -> List[List[Dict]]:
    if repeats <= 0:
        raise ValueError(f'Repeats must be greater than 0! ({repeats})')

    smi, ensemble_infiles = ligand

    ensemble_rowss = []
    for infile in ensemble_infiles:
        repeat_rows = []
        for repeat in range(repeats):
            name = f'{Path(infile).stem}__{repeat}'

            outfile = Path(path / Path(infile).with_suffix('.out'))
            argv = ['-i', infile, '-o', outfile]

            score = run_and_parse_docker(argv, 'dock', outfile, score_mode)

            if score:
                repeat_rows.append({
                    'smiles': smi,
                    'name': name,
                    'in': infile,
                    'out': str(outfile),
                    'score': score
                })

        if repeat_rows:
            ensemble_rowss.append(repeat_rows)

    return ensemble_rowss

def dock_inputs():
    pass