import random
import uuid

import pytest

from pyscreener.docking import Simulation, Result
from pyscreener.exceptions import InvalidResultError, NotSimulatedError


@pytest.fixture(params=["CCCCCCC", "C1CCC1", "CC(=O)CC", "CCCCCCCC", "CCCC1CC1"])
def smi(request):
    return request.param


def test_notsimulated(smi):
    data = Simulation(smi, None, None, None, None)
    with pytest.raises(NotSimulatedError):
        data.score


def test_invalid_result(smi):
    data = Simulation(smi, None, None, None, None)
    data.result = {"score": random.random()}

    with pytest.raises(InvalidResultError):
        data.score


def test_score(smi):
    data = Simulation(smi, None, None, None, None)
    score = random.random()
    data.result = Result(smi, "ligand", str(uuid.uuid4()), score)

    assert data.result.score == score
