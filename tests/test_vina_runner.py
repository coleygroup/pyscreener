from pathlib import Path
import random
import uuid

import pytest

from pyscreener.docking import CalculationData, vina
from pyscreener.exceptions import NotSimulatedError
from pyscreener.utils import calc_score
from pyscreener.utils.utils import ScoreMode

TEST_DIR = Path(__file__).parent

RECEPTOR_FILEPATH = TEST_DIR / '5WIU.pdb'

@pytest.fixture(
    params=["CCCC", "c1ccccc1", "CC(=O)CC", "CN=C=O", "CC(=O)C"]
)
def smi(request):
    return request.param

@pytest.fixture
def receptor():
    return TEST_DIR / '5WIU.pdb'

@pytest.fixture
def center():
    return (-18.2, 14.4, -16.1)

@pytest.fixture
def size():
    return (15.4, 13.9, 14.5)

@pytest.fixture
def in_path(tmp_path):
    p = tmp_path / 'inputs'
    p.mkdir(parents=True, exist_ok=True)

    return p

@pytest.fixture
def out_path(tmp_path):
    p = tmp_path / 'outputs'
    p.mkdir(parents=True, exist_ok=True)

    return p

@pytest.fixture
def data(smi, receptor, center, size, in_path, out_path):
    return CalculationData(
        smi, receptor, center, size, vina.VinaMetadata(),
        -1, 'ligand', None, in_path, out_path
    )

def test_prepare_ligand(data):
    data = vina.VinaRunner.prepare_ligand(data)

    assert data.metadata.prepared_ligand.exists()

def test_prepare(receptor, center, size, in_path, out_path):
    data = CalculationData(
        "c1ccccc1", receptor, center, size, vina.VinaMetadata(),
        4, 'ligand', None, in_path, out_path
    )
    data = vina.VinaRunner.prepare(data)

    assert data.metadata.prepared_receptor.exists()
    assert data.metadata.prepared_ligand.exists()

def test_run(data):
    vina.VinaRunner.prepare(data)

    with pytest.raises(NotSimulatedError):
        data.score

    scores = vina.VinaRunner.run(data)

    assert data.score == calc_score(scores, data.score_mode, data.k)