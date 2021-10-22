import json
import sys

import pyscreener as ps


def main():
    ps.check_env(sys.argv[1], json.loads(sys.argv[2]))

if __name__ == "__main__":
    main()