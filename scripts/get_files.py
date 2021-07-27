import argparse
from collections import defaultdict#, namedtuple
import csv
from pathlib import Path
import tarfile

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--smis', nargs='+',
                        help='the individual SMILES strings you would like to extract files for. NOTE: the SMILES strings must be quoted to avoid command line parsing issues.')
    parser.add_argument('-f', '--files', nargs='+',
                        help='a file (or files) containing the desired SMILES strings, each on a separate line. NOTE: there must be no title line!')
    parser.add_argument('-d', '--root-dir', required=True,
                        help='the output directory of a pyscreener run')
    parser.add_argument('-p', '--path',
                        help='the desired output directory for the extracted input and output file. By default, use the pyscreener output directory.')

    args = parser.parse_args()
    smis = set(args.smis)
    
    for f in args.files:
        with open(f) as fid:
            smis.update(fid.read().splitlines())

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