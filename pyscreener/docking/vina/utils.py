from enum import Enum, auto

class Software(Enum):
    def _generate_next_value_(name, start, count, last_values):
        return name.lower()

    VINA = auto()
    PSOVINA = auto()
    QVINA = auto()
    SMINA = auto()