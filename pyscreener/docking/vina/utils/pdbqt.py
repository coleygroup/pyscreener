from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Iterable, Iterator, Tuple, Union

import numpy as np


class PDBQTRecord(Enum):
    RECORD_NAME = slice(0, 6)
    X_COORD = slice(30, 38)
    Y_COORD = slice(38, 46)
    Z_COORD = slice(46, 54)
    ATOM_TYPE = slice(77, 79)


class PDBQTParser:
    ATOM_TYPE_TO_ATOM = {"A": "C", "HD": "H", "NA": "N", "OA": "O"}

    @staticmethod
    def segment(lines: Iterable[str]) -> Iterator[list[str]]:
        for line in lines:
            model_lines = []
            if "MODEL" in line:
                while "ENDMDL" not in line:
                    model_lines.append(line)
                    line = next(lines)
                model_lines.append(line)
            yield model_lines

    @staticmethod
    def parse_model(lines: Iterable[str]) -> Tuple[list[str], np.ndarray]:
        atoms = []
        coords = []
        for line in lines:
            record_name = line[PDBQTRecord.RECORD_NAME.value]
            if record_name != "HETATM" and record_name != "ATOM":
                continue

            atoms.append(line[PDBQTRecord.ATOM_TYPE.value].strip())
            coords.append(
                (
                    float(line[PDBQTRecord.X_COORD.value]),
                    float(line[PDBQTRecord.Y_COORD.value]),
                    float(line[PDBQTRecord.Z_COORD.value]),
                )
            )

        return atoms, np.array(coords)

    @staticmethod
    def parse(filepath: Union[str, Path]) -> Tuple[list[str], np.ndarray]:
        atoms = []
        coordss = []

        with open(filepath) as fid:
            for model_lines in PDBQTParser.segment(fid):
                try:
                    atoms, coords = PDBQTParser.parse_model(model_lines)
                except ValueError:
                    continue
                coordss.append(coords)

        atoms = [PDBQTParser.ATOM_TYPE_TO_ATOM.get(a) or a for a in atoms]
        C = np.stack(coordss)

        return atoms, C
