from pathlib import Path

import pytest

from pyscreener.exceptions import MissingExecutableError, NotSimulatedError
from pyscreener.utils import reduce_scores

try:
    from pyscreener.docking import Simulation, vina
except MissingExecutableError:
    pytestmark = pytest.mark.skip()

TEST_DIR = Path(__file__).parent


@pytest.fixture(params=["CCCC", "c1ccccc1", "CC(=O)CC", "CN=C=O", "CC(=O)C"])
def smi(request):
    return request.param


@pytest.fixture(params=["C1Ccc", "foo", "CCCooC"])
def bad_smi(request):
    return request.param


@pytest.fixture
def receptor():
    return TEST_DIR / "data" / "5WIU.pdb"


@pytest.fixture
def center():
    return (-18.2, 14.4, -16.1)


@pytest.fixture
def size():
    return (15.4, 13.9, 14.5)


@pytest.fixture
def in_path(tmp_path):
    p = tmp_path / "inputs"
    p.mkdir(parents=True, exist_ok=True)

    return p


@pytest.fixture
def out_path(tmp_path):
    p = tmp_path / "outputs"
    p.mkdir(parents=True, exist_ok=True)

    return p


@pytest.fixture
def sim(smi, receptor, center, size, in_path, out_path):
    return Simulation(
        smi, receptor, center, size, vina.VinaMetadata(), -1, "ligand", None, in_path, out_path
    )


@pytest.fixture
def bad_sim(bad_smi, receptor, center, size, in_path, out_path):
    return Simulation(
        bad_smi, receptor, center, size, vina.VinaMetadata(), -1, "ligand", None, in_path, out_path
    )


def test_prepare_ligand_success(sim):
    success = vina.VinaRunner.prepare_ligand(sim)

    assert success
    assert sim.metadata.prepared_ligand.exists()


def test_prepare_ligand_failure(bad_sim):
    success = vina.VinaRunner.prepare_ligand(bad_sim)

    assert not success
    assert bad_sim.metadata.prepared_ligand is None


def test_prepare(receptor, center, size, in_path, out_path):
    sim = Simulation(
        "c1ccccc1",
        receptor,
        center,
        size,
        vina.VinaMetadata(),
        4,
        "ligand",
        None,
        in_path,
        out_path,
    )
    vina.VinaRunner.prepare(sim)

    assert sim.metadata.prepared_receptor.exists()
    assert sim.metadata.prepared_ligand.exists()


def test_run(sim):
    vina.VinaRunner.prepare(sim)

    with pytest.raises(NotSimulatedError):
        sim.score

    scores = vina.VinaRunner.run(sim)

    assert sim.score == reduce_scores(scores, sim.reduction, k=sim.k)
