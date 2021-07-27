import argparse
from collections import defaultdict#, namedtuple
import csv
from pathlib import Path
import tarfile

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--smis', nargs='+')
    parser.add_argument('-d', '--root-dir', required=True)
    parser.add_argument('-p', '--path')

    args = parser.parse_args()
    smis = set(args.smis)
    
    root_dir = Path(args.root_dir)
    if not root_dir.is_dir():
        raise ValueError(f'"{args.root_dir}" is not a directory!')

    d_nodeID_ligands = defaultdict(list)
    extended_csv = root_dir / 'extended.csv'
    with open(extended_csv) as fid:
        reader = csv.reader(fid)
        next(reader)

        for smi, ligand_name, node_id, _ in reader:
            if smi not in smis:
                continue
            
            d_nodeID_ligands[node_id].append(ligand_name)

    path = Path(args.path or args.root_dir)

    for node_id in d_nodeID_ligands:
        with tarfile.open(root_dir / f'{node_id}.tar.gz') as tar:
            member_names = tar.getnames()
            ligand_names = d_nodeID_ligands[node_id]
            extracted_members = [
                tar.extract(member_name, path)
                for member_name in member_names if any(
                    ligand_name in member_name for ligand_name in ligand_names
                )
            ]
            print(f'Extracted {len(extracted_members)} from node {node_id}')


if __name__ == '__main__':
    main()