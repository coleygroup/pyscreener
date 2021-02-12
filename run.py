import csv
from distutils.dir_util import copy_tree
import os
from pathlib import Path

import ray

import pyscreener
from pyscreener import args

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

    try:
        if 'redis_password' in os.environ:
            ray.init(
                address='auto',
                #_node_ip_address=os.environ["ip_head"].split(":")[0], 
                _redis_password=os.environ['redis_password']
            )
        else:
            ray.init(address='auto')
    except ConnectionError:
        ray.init(_temp_dir=params['tmp_dir'])
    except PermissionError:
        print('Failed to create a temporary directory for ray')
        raise
    print(ray.cluster_resources())

    

    name = params['name']

    # tmp_dir = Path(params['tmp_dir']) / name
    out_dir = Path(params['root']) / name
    # if not tmp_dir.exists():
    #     tmp_dir.mkdir(parents=True)
    params['path'] = out_dir
    
    print('Preprocessing ...', flush=True)
    params = pyscreener.preprocess(**params)
    print('Done!')

    print(f'Preparing and screening inputs ...', flush=True)
    screener = pyscreener.build_screener(**params)
    d_smi_score, rows = screener(
        *params['ligands'], full_results=True, **params
    )
    print('Done!')

    out_dir = Path(params['root']) / name
    if not out_dir.exists():
        out_dir.mkdir(parents=True)
    params['path'] = out_dir
    print(f'Postprocessing ...', flush=True)
    pyscreener.postprocess(d_smi_score=d_smi_score, **params)
    print('Done!')

    if params['collect_all']:
        print('Collecting all input and output files ...', end=' ', flush=True)
        screener.collect_files(out_dir)
        print('Done!')

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
        writer.writerow(['smiles', 'name', 'node_id', 'score'])
        writer.writerows(row.values() for row in rows)

    print(f'Scoring data has been saved to: "{scores_filename}"')
    print(f'Extended data has been saved to: "{extended_filename}"')
    print('Thanks for using Pyscreener!')

if __name__ == '__main__':
    main()
