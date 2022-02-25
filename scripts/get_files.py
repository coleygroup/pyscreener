import argparse
from collections import defaultdict
import csv
from pathlib import Path
import tarfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-s",
        "--smis",
        nargs="+",
        help="the individual SMILES strings you would like to extract files for. NOTE: the SMILES strings must be quoted to avoid command line parsing issues.",
    )
    parser.add_argument(
        "-f",
        "--files",
        nargs="+",
        help="a file (or files) containing the desired SMILES strings, each on a separate line. NOTE: there must be no title line!",
    )
    parser.add_argument(
        "-o", "--output-dir", required=True, help="the output directory of a pyscreener run"
    )
    parser.add_argument(
        "-p",
        "--path",
        help="the desired output directory for the extracted input and output file. By default, use the pyscreener output directory.",
    )

    args = parser.parse_args()
    smis = set(args.smis)

    for f in args.files:
        with open(f) as fid:
            smis.update(fid.read().splitlines())

    output_dir = Path(args.output_dir)
    if not output_dir.is_dir():
        raise ValueError(f'"{args.output_dir}" is not a directory!')

    d_nodeID_names = defaultdict(list)
    extended_csv = output_dir / "extended.csv"
    with open(extended_csv) as fid:
        reader = csv.reader(fid)
        next(reader)

        for smi, name, node_id, _ in reader:
            if smi not in smis:
                continue

            d_nodeID_names[node_id].append(name)

    path = Path(args.path or args.output_dir)

    for node_id in d_nodeID_names:
        with tarfile.open(output_dir / f"{node_id}.tar.gz") as tar:
            names = d_nodeID_names[node_id]
            targets = [member for member in tar.getnames() if any(name in member for name in names)]
            extracted_members = [tar.extract(target, path) for target in targets]
            print(f"Extracted {len(extracted_members)} from node {node_id}")


if __name__ == "__main__":
    main()
