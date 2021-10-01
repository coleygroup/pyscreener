import csv
import os
from pathlib import Path

import ray

from pyscreener.args import gen_args
from pyscreener.docking import screen

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
    args = gen_args()
    params = vars(args)
    print('Pyscreener will be run with the following arguments:')
    for param, value in sorted(params.items()):
        print(f'  {param}: {value}')
    print(flush=True)

    try:
        if 'redis_password' in os.environ:
            ray.init(
                address=os.environ["ip_head"],
                _node_ip_address=os.environ["ip_head"].split(":")[0], 
                _redis_password=os.environ['redis_password']
            )
        else:
            ray.init(address='auto')
    except ConnectionError:
        ray.init()
    except PermissionError:
        print('Failed to create a temporary directory for ray')
        raise
    
    print('Ray cluster online with resources:')
    print(ray.cluster_resources())
    print(flush=True)

    out_dir = Path(args.output_dir)
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

    print(f'Postprocessing ...', flush=True)
    pyscreener.postprocess(d_smi_score=d_smi_score, **params)
    print('Done!')

    scores_filename = out_dir / 'scores.csv'
    with open(scores_filename, 'w') as fid:
        writer = csv.writer(fid)
        writer.writerow(['smiles', 'score'])
        writer.writerows(sorted(
            d_smi_score.items(),
            key=None if args.no_sort else lambda k_v: k_v[1] or float('inf')
        ))
    print(f'Scoring data has been saved to: "{scores_filename}"')

    if args.collect_all:
        print('Collecting all input and output files ...', end=' ', flush=True)
        screener.collect_files(out_dir)
        extended_filename = out_dir / 'extended.csv'
        rows = sorted(rows, key=lambda row: row['score'] or float('inf'))
        with open(extended_filename, 'w') as fid:
            writer = csv.writer(fid)
            writer.writerow(['smiles', 'name', 'node_id', 'score'])
            writer.writerows(row.values() for row in rows)
        print('Done!')
        print(f'Extended data has been saved to: "{extended_filename}"')

    print('Thanks for using Pyscreener!')

if __name__ == '__main__':
    main()
