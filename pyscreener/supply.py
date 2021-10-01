import csv
from pathlib import Path
from typing import Iterable, List, Optional, Union

from openbabel import pybel
from rdkit.Chem import AllChem as Chem

from pyscreener.utils import FileType


class LigandSupply:
    def __init__(
        self,
        filepaths: Iterable[Union[str, Path]],
        filetypes: Optional[Iterable[Union[str, FileType]]] = None,
        use_3d: bool = False,
        optimize: bool = False,
        title_line: bool = True,
        smiles_col: int = 0,
        name_col: int = 1,
        id_property: Optional[str] = None,
    ):
        self.filepaths = [Path(filepath) for filepath in filepaths]
        if filetypes is not None:
            self.filetypes = [
                filetype if isinstance(filetype, FileType) else FileType.from_str(filetype)
                for filetype in filetypes
            ]
        else:
            self.filetypes = [LigandSupply.guess_filetype(filepath) for filepath in self.filepaths]

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
                    LigandSupply.get_ligands_from_file(filepath, self.use_3d, self.optimize)
                )
            elif filetype == FileType.SDF:
                ligands.extend(
                    LigandSupply.get_ligands_from_sdf(
                        filepath, self.id_property, self.use_3d, self.optimize
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
        filepath: Path,
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

        return LigandSupply.optimize_and_write_mols(mols, filepath.parent / filepath.stem)

    @staticmethod
    def get_ligands_from_file(
        filepath: Path, use_3d: bool = False, optimize: bool = False
    ) -> List[str]:
        if use_3d:
            return LigandSupply.split_file(filepath)

        fmt = Path(filepath).suffix.strip(".")

        mols = list(pybel.readfile(fmt, filepath))
        smis = [mol.write() for mol in mols]

        if not optimize:
            return [mol.write() for mol in mols]

        mols = [Chem.MolFromSmiles(smi) for smi in smis]

        return LigandSupply.optimize_and_write_mols(mols, filepath.parent / filepath.stem)
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
        filepath: Path,
        id_property: Optional[str] = None,
        use_3d: bool = False,
        optimize: bool = False,
    ) -> List[str]:
        if use_3d:
            return LigandSupply.split_file(filepath)

        mols = Chem.SDMolSupplier(str(filepath))

        if not optimize:
            return [Chem.MolToSmiles(mol) for mol in mols]

        return LigandSupply.optimize_and_write_mols(mols, filepath.parent / filepath.stem)
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
        filepath: Path, id_property: Optional[str] = None, optimize: bool = False
    ) -> List[str]:
        mols = Chem.SmilesMolSupplier(str(filepath))

        if not optimize:
            return [Chem.MolToSmiles(mol) for mol in mols]

        return LigandSupply.optimize_and_write_mols(mols, filepath.parent / filepath.stem)

    @staticmethod
    def optimize_and_write_mols(
        mols: Iterable[Chem.Mol], base_name: Union[str, Path]
    ) -> List[str]:
        # NOTE(degraff): this can parallelized later
        mols = [Chem.AddHs(mol) for mol in mols]
        [Chem.EmbedMolecule(mol) for mol in mols]
        [Chem.MMFFOptimizeMolecule(mol) for mol in mols]
        filenames = [f"{base_name}_{i}.mol" for i in range(len(mols))]
        [Chem.MolToMolFile(mol, filename) for mol, filename in zip(mols, filenames)]

        return filenames

    @staticmethod
    def split_file(filepath: Path) -> List[str]:
        fmt = filepath.suffix.strip(".")
        base_name = filepath.parent / filepath.stem

        mols = list(pybel.readfile(fmt, filepath))
        filenames = [f"{base_name}_{i}.{fmt}" for i in range(len(mols))]
        [mol.write(fmt, filename, True) for mol, filename in zip(mols, filenames)]

        return filenames

    @staticmethod
    def guess_filetype(filepath: Path) -> List[FileType]:
        try:
            return FileType.from_str(filepath.suffix.strip('.'))
        except KeyError:
            return FileType.FILE