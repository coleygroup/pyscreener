import sys

from pyscreener.scripts.check import check
from pyscreener.scripts.driver import driver
from pyscreener.scripts.extract import extract
from pyscreener.scripts.setup import setup


def main():
    mode = sys.argv.pop(1)
    sys.argv[0] = f"pyscreener-{mode.lower()}"
    if mode.upper() == "CHECK":
        check()
    elif mode.upper() == "RUN":
        driver()
    elif mode.upper() == "EXTRACT":
        extract()
    elif mode.upper() == "SETUP":
        setup()
    else:
        print('usage: pyscreener MODE {"check", "driver", "extract", "setup"')
        print(f"pyscreener: error: unrecognized mode: {mode}")
        exit(1)

    exit(0)


if __name__ == "__main__":
    main()
