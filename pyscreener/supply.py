import csv
from pathlib import Path
from typing import Iterable, List, Optional, Union

from openbabel import pybel
from rdkit.Chem import AllChem as Chem

from pyscreener.utils import FileType


class LigandSupply:
    def __init__(
        self,
        filepaths: Iterable[str],
        filetypes: Iterable[Union[str, FileType]],
        use_3d: bool = False,
        optimize: bool = False,
        title_line: bool = True,
        smiles_col: int = 0,
        name_col: int = 1,
        id_property: Optional[str] = None,
    ):
        self.filepaths = filepaths
        self.filetypes = [
            filetype if isinstance(filetype, FileType) else FileType.from_str(filetype)
            for filetype in filetypes
        ]
        self.use_3d = use_3d
        self.optimize = optimize
        self.title_line = title_line
        self.smiles_col = smiles_col
        self.name_col = name_col
        self.id_property = id_property

    def ligands(self) -> List[str]:
        ligands = []
        for filepath, filetype in zip(self.filepaths, self.filetypes):
            if filetype == FileType.CSV:
                ligands.extend(
                    LigandSupply.get_ligands_from_csv(
                        filepath,
                        self.title_line,
                        self.smiles_col,
                        self.id_property,
                        self.optimize,
                    )
                )
            elif filetype == FileType.FILE:
                ligands.extend(
                    LigandSupply.get_ligands_from_file(filepath, self.optimize)
                )
            elif filetype == FileType.SDF:
                ligands.extend(
                    LigandSupply.get_ligands_from_sdf(
                        filepath, self.id_property, self.optimize
                    )
                )
            elif filetype == FileType.SMI:
                ligands.extend(
                    LigandSupply.get_ligands_from_smi(
                        filepath, self.id_property, self.optimize
                    )
                )

        return ligands

    @staticmethod
    def get_ligands_from_csv(
        filepath,
        title_line: bool = True,
        smiles_col: int = 0,
        name_col: int = 1,
        optimize: bool = False,
    ) -> List[str]:
        with open(filepath) as fid:
            reader = csv.reader(fid)
            if title_line:
                next(reader)
            smis = [row[smiles_col] for row in reader]

        if not optimize:
            return smis

        mols = [Chem.MolFromSmiles(smi) for smi in smis]

        p = Path(filepath)
        return LigandSupply.optimize_and_write_mols(mols, p.parent / p.stem)

    @staticmethod
    def get_ligands_from_file(
        filepath, use_3d: bool = False, optimize: bool = False
    ) -> List[str]:
        if use_3d:
            return LigandSupply.split_file(filepath)

        fmt = Path(filepath).suffix.strip(".")

        mols = list(pybel.readfile(fmt, filepath))
        smis = [mol.write() for mol in mols]

        if not optimize:
            return [mol.write() for mol in mols]

        mols = [Chem.MolFromSmiles(smi) for smi in smis]

        p = Path(filepath)
        return LigandSupply.optimize_and_write_mols(mols, p.parent / p.stem)
        # names = []
        # for mol in pybel.readfile(fmt, filepath):
        #     smis.append(mol.write())

        #     if mol.title:
        #         names.append(mol.title)
        #     else:
        #         names.append(None)

        # return smis, names

    @staticmethod
    def get_ligands_from_sdf(
        filepath,
        id_property: Optional[str] = None,
        use_3d: bool = False,
        optimize: bool = False,
    ) -> List[str]:
        if use_3d:
            return LigandSupply.split_file(filepath)

        mols = Chem.SDMolSupplier(filepath)

        if not optimize:
            return [Chem.MolToSmiles(mol) for mol in mols]

        p = Path(filepath)
        return LigandSupply.optimize_and_write_mols(mols, p.parent / p.stem)
        # if id_property is not None:
        #     for mol in mols:
        #         if mol is None:
        #             continue

        #         smis.append(Chem.MolToSmiles(mol))
        #         names.append(mol.GetProp(id_prop_name))
        # else:
        #     for mol in mols:
        #         if mol is None:
        #             continue

        #         smis.append(Chem.MolToSmiles(mol))
        #     names = [None for _ in smis]
        # return smis, names

    @staticmethod
    def get_ligands_from_smi(
        filepath, id_property: Optional[str] = None, optimize: bool = False
    ) -> List[str]:
        mols = Chem.SmilesMolSupplier(filepath)

        if not optimize:
            return [Chem.MolToSmiles(mol) for mol in mols]

        p = Path(filepath)
        return LigandSupply.optimize_and_write_mols(mols, p.parent / p.stem)

    @staticmethod
    def optimize_and_write_mols(
        mols: Iterable[Chem.Mol], base_name: Union[str, Path]
    ) -> List[str]:
        # NOTE(degraff): this can parallelized later
        mols = [Chem.AddHs(mol) for mol in mols]
        [Chem.EmbedMolecule(mol) for mol in mols]
        [Chem.MMFFOptimizeMolecule(mol) for mol in mols]
        filenames = [f"{base_name}_{i}.mol" for i in range(len(mols))]
        [Chem.MoltoMolFile(mol, filename) for mol, filename in zip(mols, filenames)]

        return filenames

    @staticmethod
    def split_file(filepath) -> List[str]:
        p = Path(filepath)
        fmt = p.suffix.strip(".")
        base_name = p.parent / p.stem

        mols = list(pybel.readfile(fmt, filepath))
        filenames = [f"{base_name}_{i}.{fmt}" for i in range(len(mols))]
        [mol.write(fmt, filename, True) for mol, filename in zip(mols, filenames)]

        return filenames
