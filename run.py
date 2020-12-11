import csv
from distutils.dir_util import copy_tree
from operator import itemgetter
from pathlib import Path
import tempfile

import pyscreener
from pyscreener import args, preprocess, postprocess

def main():
    print('''\
***************************************************************
*      ____  __  ____________________  ___  ____  ___  _____  *
*     / __ \/ / / / ___/ ___/ ___/ _ \/ _ \/ __ \/ _ \/ ___/  *
*    / /_/ / /_/ (__  ) /__/ /  /  __/  __/ / / /  __/ /      *
*   / .___/\__, /____/\___/_/   \___/\___/_/ /_/\___/_/       *
*  /_/    /____/                                              *
***************************************************************''')
    print('Welcome to Pyscreener!\n')

    params = vars(args.gen_args())

    print('Pyscreener will be run with the following arguments:')
    for param, value in sorted(params.items()):
        print(f'  {param}: {value}')
    print(flush=True)

    name = params['name']

    tmp_dir = Path(tempfile.gettempdir()) / name
    if not tmp_dir.exists():
        tmp_dir.mkdir(parents=True)
    params['path'] = tmp_dir
    
    print('Preprocessing ...', flush=True)
    params = pyscreener.preprocess(**params)
    print('Done!')

    print(f'Preparing and screening inputs ...', flush=True)
    screener = pyscreener.build_screener(**params)
    d_smi_score, rows = screener(*params['ligands'], full_results=True, 
                                 **params)
    print('Done!')

    out_dir = Path(params['root']) / name
    if not out_dir.exists():
        out_dir.mkdir(parents=True)
    params['path'] = out_dir
    print(f'Postprocessing ...', flush=True)
    pyscreener.postprocess(d_smi_score=d_smi_score, **params)
    print('Done!')

    if params['copy_all']:
        copy_tree(str(tmp_dir), str(out_dir))

    scores_filename = out_dir / f'{name}_scores.csv'
    extended_filename = out_dir / f'{name}_extended.csv'

    with open(scores_filename, 'w') as fid:
        writer = csv.writer(fid)
        writer.writerow(['smiles', 'score'])
        writer.writerows(
            sorted(d_smi_score.items(), key=lambda k_v: k_v[1] or float('inf'))
        )
    
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
