import csv
from distutils.dir_util import copy_tree
from operator import itemgetter
from pathlib import Path

# from pyscreener import docking_mpi as docking
# from pyscreener import conversion_mpi as conversion
# from pyscreener.args import gen_args

from pyscreener import (args, preprocessing, preparation, 
                        screening, postprocessing)

def main():
    print('''\
*****************************************************************
*      ____  __  ____________________  ___  ____  ___  _____    *
*     / __ \/ / / / ___/ ___/ ___/ _ \/ _ \/ __ \/ _ \/ ___/    *
*    / /_/ / /_/ (__  ) /__/ /  /  __/  __/ / / /  __/ /        *
*   / .___/\__, /____/\___/_/   \___/\___/_/ /_/\___/_/         *
*  /_/    /____/                                                *
*****************************************************************''')
    print('Welcome to Pyscreener!\n')

    params = vars(args.gen_args())

    print('Pyscreener will be run with the following arguments:')
    for param, value in sorted(params.items()):
        print(f'  {param}: {value}')
    print(flush=True)

    print('Preprocessing ...', end=' ', flush=True)
    params = preprocessing.preprocess(**params)
    print('Done!')

    name = params['name']
    base_tmp_path = Path(params['tmp']) / name
    output_dir = Path(params['root']) / name
    
    print('Preparing inputs ...', flush=True)
    inputs = preparation.prepare(path=base_tmp_path, **params)
    print('Done!')

    print(f'Screening inputs ...', flush=True)
    d_smi_score, rows = screening.screen(path=base_tmp_path,
                                         inputs=inputs, **params)
    print('Done!')

    print(f'Postprocessing ...', flush=True)
    postprocessing.postprocess(
        d_smi_score=d_smi_score, path=name, **params)
    print('Done!')

    
    copy_tree(str(base_tmp_path), str(output_dir))

    scores_filename = output_dir / f'{name}_scores.csv'
    extended_filename = output_dir / f'{name}_extended.csv'

    smis_scores = sorted(d_smi_score.items(), key=itemgetter(1))
    with open(scores_filename, 'w') as fid:
        writer = csv.writer(fid)
        writer.writerow(['smiles', 'score'])
        writer.writerows(smis_scores)
    
    rows = sorted(rows, key=lambda row: row['score'] or float('inf'))
    with open(extended_filename, 'w') as fid:
        writer = csv.writer(fid)
        writer.writerow(
            ['smiles', 'name', 'input_file', 'out_file', 'log_file', 'score'])
        writer.writerows(row.values() for row in rows)

    print(f'Scoring data has been saved to: "{scores_filename}"')
    print(f'Extended data has been saved to: "{extended_filename}"')
    print('Thanks for using Pyscreener!')

if __name__ == '__main__':
    main()
