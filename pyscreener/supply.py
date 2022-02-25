from __future__ import annotations

import csv
from collections.abc import Iterable
from pathlib import Path
from typing import Optional, Union

from openbabel import pybel
from rdkit.Chem import AllChem as Chem

from pyscreener.utils import FileFormat


class LigandSupply(Iterable):
    """A LigandSupply is represents an abstract collection of molecular supply files, allowing for
    the iteration between all molecules contained in a variety of file formats

    Attributes
    ----------
    ligands : Iterable[str]
        the SMILES strings or files containing the molecules of the supplied inputs
    filepaths : Iterable[Path]
        the chemical supply files containing the molecules
    formats : Iterable[FileFormat]
        the corresponding format for each input file
    use_3d : bool
        whether to use the input 3D geometry of the supplied molecules. Does nothing if the file is
        a CSV or SMI file.
    optimize : bool
        whether to perform geometry optimization of all input molecules and write the geometry to
        an output MOL file for each supplied molecule. NOTE: this will ablate any input geometries.
    title_line : bool
        whether there is a title line in CSV or SMI files
    smiles_col : int
        the column containing the SMILES strings in CSV or SMI files
    name_col : int
        the column containing the molecule name in CSV or SMI files
    id_property : Optional[str]
        NOTE: Unused. the name of the molecular property containing the molecule's ID
    path : Optional[Path]
        The path under which to output all prepared files. If None, prepare all input files in the
        same directory as their parent file.

    Parameters
    ----------
    filepaths : Iterable[Union[str, Path]]
        the chemical supply files containing the molecules
    formats : Optional[Iterable[Union[str, FileFormat]]], default=None
        the corresponding filetype for each input file. Can be supplied as either a string or a
        Filetype for each file. E.g., for a csv file, filetype can be "csv" or FileFormat.CSV.
        If None, guess the type of each file
    smis : Iterable[str], default=None
        the SMILES strings corresponding to molecules to prepare
    use_3d : bool, default=False
        whether to use the input 3D geometry of the supplied molecules. Does nothing if the file is
        a CSV or SMI file.
    optimize : bool, default=False
        whether to perform geometry optimization of all input molecules and write the geometry to
        an output MOL file for each supplied molecule. NOTE: this will ablate any input geometries.
    title_line : bool, default=True
        whether there is a title line in CSV or SMI files
    smiles_col : int, default=0
        the column containing the SMILES strings in CSV or SMI files
    name_col : int, default=1
        the column containing the molecule name in CSV or SMI files
    id_property : Optional[str], default=None
        NOTE: Unused. the name molecular property containing the molecule's ID
    path : Optional[str], default=None
        The path under which to output all prepared files. If None, prepare all input files in the
        same directory as their parent file and all SMILES strings in the current directory.
    """

    def __init__(
        self,
        filepaths: Iterable[Union[str, Path]],
        formats: Optional[Iterable[Union[str, FileFormat]]] = None,
        smis: Optional[Iterable[str]] = None,
        use_3d: bool = False,
        optimize: bool = False,
        title_line: bool = True,
        smiles_col: int = 0,
        name_col: int = 1,
        id_property: Optional[str] = None,
        path: Optional[str] = None,
    ):
        self.filepaths = [Path(filepath) for filepath in filepaths]
        if formats is not None:
            self.formats = [
                filetype if isinstance(filetype, FileFormat) else FileFormat.from_str(filetype)
                for filetype in formats
            ]
        else:
            self.formats = [LigandSupply.guess_format(filepath) for filepath in self.filepaths]

        self.use_3d = use_3d
        self.optimize = optimize
        self.title_line = title_line
        self.smiles_col = smiles_col
        self.name_col = name_col
        self.id_property = id_property
        self.path = Path(path) if path is not None else None

        ligands = []
        for filepath, filetype in zip(self.filepaths, self.formats):
            if filetype == FileFormat.CSV:
                ligands.extend(
                    LigandSupply.get_ligands_from_csv(
                        filepath,
                        self.title_line,
                        self.smiles_col,
                        self.id_property,
                        self.optimize,
                        self.path,
                    )
                )
            elif filetype == FileFormat.FILE:
                ligands.extend(
                    LigandSupply.get_ligands_from_file(
                        filepath, self.use_3d, self.optimize, self.path
                    )
                )
            elif filetype == FileFormat.SDF:
                ligands.extend(
                    LigandSupply.get_ligands_from_sdf(
                        filepath, self.id_property, self.use_3d, self.optimize, self.path
                    )
                )
            elif filetype == FileFormat.SMI:
                ligands.extend(
                    LigandSupply.get_ligands_from_smi(
                        filepath, self.id_property, self.optimize, self.path
                    )
                )
        if smis is not None:
            if not self.optimize:
                ligands.extend(smis)
            else:
                mols = [Chem.MolFromSmiles(smi) for smi in smis]
                ligands.extend(
                    LigandSupply.optimize_and_write_mols(mols, Path("ligand"), self.path)
                )
        self.ligands = ligands

    def __len__(self):
        return len(self.ligands)

    def __iter__(self) -> Iterable[str]:
        return iter(self.ligands)

    def __getitem__(self, i: int) -> str:
        return self.ligands[i]

    @staticmethod
    def get_ligands_from_csv(
        filepath: Path,
        title_line: bool = True,
        smiles_col: int = 0,
        name_col: int = 1,
        optimize: bool = False,
        path: Optional[Path] = None,
    ) -> list[str]:
        with open(filepath) as fid:
            reader = csv.reader(fid)
            if title_line:
                next(reader)
            smis = [row[smiles_col] for row in reader]

        if not optimize:
            return smis

        mols = [Chem.MolFromSmiles(smi) for smi in smis]

        return LigandSupply.optimize_and_write_mols(mols, filepath, path)

    @staticmethod
    def get_ligands_from_file(
        filepath: Path, use_3d: bool = False, optimize: bool = False, path: Optional[Path] = None
    ) -> list[str]:
        if use_3d:
            return LigandSupply.split_file(filepath)

        fmt = Path(filepath).suffix.strip(".")

        mols = list(pybel.readfile(fmt, filepath))
        smis = [mol.write() for mol in mols]

        if not optimize:
            return [mol.write() for mol in mols]

        mols = [Chem.MolFromSmiles(smi) for smi in smis]

        return LigandSupply.optimize_and_write_mols(mols, filepath, path)

    @staticmethod
    def get_ligands_from_sdf(
        filepath: Path,
        id_property: Optional[str] = None,
        use_3d: bool = False,
        optimize: bool = False,
        path: Optional[Path] = None,
    ) -> list[str]:
        if use_3d:
            return LigandSupply.split_file(filepath, path)

        mols = Chem.SDMolSupplier(str(filepath))

        if not optimize:
            return [Chem.MolToSmiles(mol) for mol in mols]

        return LigandSupply.optimize_and_write_mols(mols, filepath, path)

    @staticmethod
    def get_ligands_from_smi(
        filepath: Path,
        id_property: Optional[str] = None,
        optimize: bool = False,
        path: Optional[Path] = None,
    ) -> list[str]:
        mols = Chem.SmilesMolSupplier(str(filepath))

        if not optimize:
            return [Chem.MolToSmiles(mol) for mol in mols]

        return LigandSupply.optimize_and_write_mols(mols, filepath, path)

    @staticmethod
    def optimize_and_write_mols(
        mols: Iterable[Chem.Mol], filepath: Path, path: Optional[Path] = None
    ) -> list[str]:
        base_name = (path or filepath.parent) / filepath.stem
        # NOTE(degraff): this can parallelized later
        mols = [Chem.AddHs(mol) for mol in mols]
        [Chem.EmbedMolecule(mol) for mol in mols]
        [Chem.MMFFOptimizeMolecule(mol) for mol in mols]

        filenames = [f"{base_name}_{i}.mol" for i in range(len(mols))]
        [Chem.MolToMolFile(mol, filename) for mol, filename in zip(mols, filenames)]

        return filenames

    @staticmethod
    def split_file(filepath: Path, path: Optional[Path] = None) -> list[str]:
        fmt = filepath.suffix.strip(".")
        base_name = (path or filepath.parent) / filepath.stem

        mols = list(pybel.readfile(fmt, str(filepath)))

        filenames = [f"{base_name}_{i}.{fmt}" for i in range(len(mols))]
        [mol.write(fmt, filename, True) for mol, filename in zip(mols, filenames)]

        return filenames

    @staticmethod
    def guess_format(filepath: Path) -> FileFormat:
        try:
            return FileFormat.from_str(filepath.suffix.strip("."))
        except KeyError:
            return FileFormat.FILE
