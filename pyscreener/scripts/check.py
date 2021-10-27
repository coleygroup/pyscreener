import argparse
import json

import pyscreener as ps


def check():
    parser = argparse.ArgumentParser()
    parser.add_argument("software", help="the screening software")
    parser.add_argument("metadata_template", nargs="?", type=json.loads, default={})
    args = parser.parse_args()

    ps.check_env(args.software, args.metadata_template)


if __name__ == "__main__":
    check()
