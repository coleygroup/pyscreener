import csv
from pathlib import Path

import pytest
from rdkit import Chem

from pyscreener.supply import LigandSupply
from pyscreener.utils import FileType


@pytest.fixture(
    params=[
        ["CCCCCCC", "C1CCC1", "CC(=O)CC", "CCCCCCCC", "CCCC1CC1"],
        ["CC(=O)NCCC1=CNc2c1cc(OC)cc2", "CN=C=O", "CN1CCC[C@H]1c2cccnc2"],
        [
            "O1C=C[C@H]([C@H]1O2)c3c2cc(OC)c4c3OC(=O)C5=C4CCC(=O)5",
            "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)1",
        ],
    ]
)
def smis(request):
    return request.param


@pytest.fixture(params=["csv", "smi", "sdf"])
def filetype(request):
    return request.param


def mols(smis):
    return [Chem.MolFromSmiles(smi) for smi in smis]


def make_csv(smis, path):
    p_csv = path / "mols.csv"
    with open(p_csv, "w") as fid:
        writer = csv.writer(fid)
        writer.writerow(["smiles"])
        writer.writerows(([smi] for smi in smis))

    return p_csv


def make_sdf(smis, path):
    p_sdf = path / "mols.sdf"
    with Chem.SDWriter(str(p_sdf)) as w:
        for m in mols(smis):
            w.write(m)

    return p_sdf


def make_smi(smis, path):
    p_smi = path / "mols.smi"
    with Chem.SmilesWriter(str(p_smi)) as w:
        for m in mols(smis):
            w.write(m)

    return p_smi


def make_file(smis, tmp_path, filetype):
    if filetype == "csv":
        f = make_csv
    elif filetype == "sdf":
        f = make_sdf
    elif filetype == "smi":
        f = make_smi

    return f(smis, tmp_path)


def test_guess_filetype():
    filepaths = ["test.csv", "test.sdf", "test.smi", "test.pdb", "test.mol2"]
    supply = LigandSupply(filepaths)

    assert supply.filetypes == [
        FileType.CSV,
        FileType.SDF,
        FileType.SMI,
        FileType.FILE,
        FileType.FILE,
    ]


def test_ligands(smis, tmp_path, filetype):
    supply = LigandSupply([make_file(smis, tmp_path, filetype)])

    assert len(supply.ligands()) == len(smis)


def test_optimize(smis, tmp_path, filetype):
    supply = LigandSupply([make_file(smis, tmp_path, filetype)], optimize=True)
    ligands = supply.ligands()

    assert len(ligands) == len(smis)

    for ligand in ligands:
        assert Path(ligand).exists()


def test_multiple_files(smis, tmp_path):
    N = 3
    filepaths = [make_csv(smis, tmp_path) for _ in range(N)]
    supply = LigandSupply(filepaths)
    ligands = supply.ligands()

    assert len(ligands) == N * len(smis)


def test_multiple_filetypes_optimize(smis, tmp_path):
    filepaths = [
        make_csv(smis, tmp_path),
        make_sdf(smis, tmp_path),
        make_smi(smis, tmp_path)
    ]
    supply = LigandSupply(filepaths, optimize=True)
    ligands = supply.ligands()

    assert len(ligands) == len(filepaths) * len(smis)
    
    for ligand in ligands:
        assert Path(ligand).exists()
