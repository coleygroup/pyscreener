from os import read
import random
import uuid

import pytest

from pyscreener.docking import CalculationData, Result
from pyscreener.exceptions import InvalidResultError, NotSimulatedError

import tempfile


@pytest.fixture(
    params=["CCCCCCC", "C1CCC1", "CC(=O)CC", "CCCCCCCC", "CCCC1CC1"]
)
def smi(request):
    return request.param

def test_notsimulated(smi):
    data = CalculationData(smi, None, None, None, None)
    with pytest.raises(NotSimulatedError):
        data.score

def test_invalid_result(smi):
    data = CalculationData(smi, None, None, None, None)
    data.result = {"score": random.random()}

    with pytest.raises(InvalidResultError):
        data.score

def test_score(smi):
    data = CalculationData(smi, None, None, None, None)
    score = random.random()
    data.result = Result(smi, 'ligand', str(uuid.uuid4()), score)

    assert data.result.score == score

def test_file_bytes(smi):
    tempf = tempfile.TemporaryFile()
    tempf.write(b'Hello world!')

    data = CalculationData(smi, None, None, None, None, input_file_bytes = tempf)

    assert data.input_file_bytes is tempf

