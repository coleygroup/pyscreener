from enum import auto

from pyscreener.utils import AutoName

from .pdbqt import *


class Software(AutoName):
    def _generate_next_value_(name, start, count, last_values):
        return name.lower()

    VINA = auto()
    PSOVINA = auto()
    QVINA = auto()
    SMINA = auto()
