from __future__ import annotations

from pathlib import Path
import numpy as np

import pytest

from pyscreener.docking.vina.utils.pdbqt import PDBQTParser

DATA_DIR = Path(__file__).parent / "data"


@pytest.mark.parametrize(
    "line,atom,x,y,z",
    [
        (
            "HETATM    1  CAA JQ1 A   1      29.899  18.962   1.724  1.00  7.28           C",
            "C",
            29.899,
            18.962,
            1.724,
        ),
        (
            "HETATM    2  CAB JQ1 A   1      25.204  15.116   1.393  1.00  9.95           C",
            "C",
            25.204,
            15.116,
            1.393,
        ),
        (
            "HETATM    3  CAC JQ1 A   1      25.682  14.930  -1.694  1.00 10.30           C",
            "C",
            25.682,
            14.930,
            -1.694,
        ),
        (
            "HETATM    6  CAF JQ1 A   1      28.543  11.742  -3.510  1.00 15.62           C",
            "C",
            28.543,
            11.742,
            -3.510,
        ),
        (
            "HETATM    7  OAG JQ1 A   1      32.109  13.864  -4.714  1.00 20.67           O",
            "O",
            32.109,
            13.864,
            -4.714,
        ),
        (
            "ATOM     17 CLAH JQ1 A   1      24.824  18.476  -6.770  0.00  0.00    +0.000 Cl",
            "Cl",
            24.824,
            18.476,
            -6.770,
        ),
        (
            "ATOM     10  H03 JFX A 501      51.934  -6.998 -21.457  1.00  0.00     0.186 HD",
            "HD",
            51.934,
            -6.998,
            -21.457,
        ),
        (
            "HETATM   16  O16 JFX A 501      45.628  -5.094 -15.332  1.00 64.11    -0.206 OA",
            "OA",
            45.628,
            -5.094,
            -15.332,
        ),
    ],
)
def test_parse_line(line, atom, x, y, z):
    a, xyz = PDBQTParser.parse_line(line)

    assert a == atom
    for c1, c2 in zip(xyz, (x, y, z)):
        assert pytest.approx(c1) == c2


@pytest.mark.parametrize(
    "filepath,n_models_true", [(DATA_DIR / "jfx.pdbqt", 9), (DATA_DIR / "0yb.pdbqt", 2)]
)
def test_num_models(filepath, n_models_true):
    with open(filepath) as fid:
        n_models_obs = sum(1 for _ in PDBQTParser.segment(fid))

    assert n_models_obs == n_models_true


@pytest.mark.parametrize("filepath", [(DATA_DIR / "jfx.pdbqt"), (DATA_DIR / "0yb.pdbqt")])
def test_segmentation(filepath):
    with open(filepath) as fid:
        for lines in PDBQTParser.segment(fid):
            assert "MODEL" in lines[0]
            assert "ENDMDL" in lines[-1]


@pytest.mark.parametrize(
    "lines,atoms_true",
    [
        (
            """MODEL 1
REMARK VINA RESULT:      -7.5      0.000      0.000
REMARK  6 active torsions:
REMARK  status: ('A' for Active; 'I' for Inactive)
REMARK    1  A    between atoms: C08_6  and  N09_18 
REMARK    2  A    between atoms: C11_7  and  C12_8 
REMARK    3  A    between atoms: C11_7  and  S10_22 
REMARK    4  A    between atoms: C12_8  and  C13_9 
REMARK    5  A    between atoms: C13_9  and  C14_10 
REMARK    6  A    between atoms: N09_18  and  S10_22 
ROOT
HETATM    1  C01 JFX A 501      52.052  -7.816 -18.741  1.00 62.30     0.070 C 
HETATM    2  C02 JFX A 501      50.947  -7.266 -19.636  1.00 69.74     0.034 A 
HETATM    3  C05 JFX A 501      48.968  -6.525 -20.372  1.00 66.41     0.060 A 
HETATM    4  C06 JFX A 501      49.583  -7.021 -19.238  1.00 67.31     0.018 A 
HETATM    5  C07 JFX A 501      48.844  -7.172 -18.040  1.00 64.74     0.046 A 
HETATM    6  C08 JFX A 501      47.491  -6.823 -18.011  1.00 66.27     0.043 A 
HETATM    7  C18 JFX A 501      46.867  -6.322 -19.163  1.00 58.48     0.032 A 
HETATM    8  C19 JFX A 501      47.595  -6.174 -20.338  1.00 57.92     0.034 A 
HETATM    9  N03 JFX A 501      51.078  -6.931 -20.907  1.00 73.37    -0.304 N 
ATOM     10  H03 JFX A 501      51.934  -6.998 -21.457  1.00  0.00     0.186 HD
HETATM   11  N04 JFX A 501      49.895  -6.487 -21.351  1.00 70.74    -0.300 N 
ATOM     12  H04 JFX A 501      49.719  -6.164 -22.302  1.00  0.00     0.187 HD
ENDROOT
BRANCH   6  13
HETATM   13  N09 JFX A 501      46.707  -6.967 -16.798  1.00 70.80    -0.285 N 
ATOM     14  H09 JFX A 501      46.050  -7.744 -16.721  1.00  0.00     0.178 HD
BRANCH  13  15
HETATM   15  S10 JFX A 501      46.853  -5.872 -15.495  1.00 67.91     0.231 S 
HETATM   16  O16 JFX A 501      45.628  -5.094 -15.332  1.00 64.11    -0.206 OA
HETATM   17  O17 JFX A 501      47.830  -4.827 -15.809  1.00 69.85    -0.206 OA
BRANCH  15  18
HETATM   18  C11 JFX A 501      47.292  -6.762 -13.970  1.00 70.63     0.163 C 
BRANCH  18  19
HETATM   19  C12 JFX A 501      46.157  -7.680 -13.526  1.00 65.37     0.019 C 
BRANCH  19  20
HETATM   20  C13 JFX A 501      46.690  -8.854 -12.707  1.00 65.38     0.032 C 
BRANCH  20  21
HETATM   21  C14 JFX A 501      45.592  -9.858 -12.373  1.00 68.91     0.222 C 
HETATM   22  F15 JFX A 501      45.860 -10.476 -11.188  1.00 69.95    -0.255 F 
ENDBRANCH  20  21
ENDBRANCH  19  20
ENDBRANCH  18  19
ENDBRANCH  15  18
ENDBRANCH  13  15
ENDBRANCH   6  13
TORSDOF 6
ENDMDL""",
            [
                "C",
                "A",
                "A",
                "A",
                "A",
                "A",
                "A",
                "A",
                "N",
                "HD",
                "N",
                "HD",
                "N",
                "HD",
                "S",
                "OA",
                "OA",
                "C",
                "C",
                "C",
                "C",
                "F",
            ],
        ),
        (
            """MODEL 1
REMARK VINA RESULT:      -6.4      0.000      0.000
REMARK  8 active torsions:
REMARK  status: ('A' for Active; 'I' for Inactive)
REMARK    1  A    between atoms: C25_8  and  C36_17 
REMARK    2  A    between atoms: C16_28  and  C22_9 
REMARK    3  A    between atoms: C29_10  and  O28_38 
REMARK    4  A    between atoms: C26_11  and  O28_38 
REMARK    5  A    between atoms: C10_22  and  N4_32 
REMARK    6  A    between atoms: C16_28  and  C17_27 
REMARK    7  A    between atoms: C17_27  and  C18_29 
REMARK    8  A    between atoms: C16_28  and  N14_33 
ROOT
HETATM    1  C16 0YB A 301      20.690 103.282  16.649  1.00 23.40     0.256 C 
ENDROOT
BRANCH   1   2
HETATM    2  N14 0YB A 301      21.630 102.590  17.558  1.00 22.62     0.099 N 
HETATM    3  C13 0YB A 301      22.119 103.416  18.677  1.00 23.46     0.273 C 
HETATM    4  C12 0YB A 301      22.924 102.561  19.646  1.00 22.53     0.032 C 
HETATM    5  C11 0YB A 301      24.092 101.878  18.931  1.00 20.23     0.029 C 
HETATM    6  C10 0YB A 301      23.644 101.183  17.634  1.00 19.27     0.122 C 
HETATM    7  C15 0YB A 301      22.783 102.087  16.806  1.00 19.73     0.298 C 
BRANCH   6   8
HETATM    8  N4  0YB A 301      24.703 100.580  16.853  1.00 17.74    -0.303 N 
HETATM    9  C3  0YB A 301      24.708 100.670  15.478  1.00 16.77     0.108 A 
HETATM   10  C2  0YB A 301      25.688 100.139  14.725  1.00 14.10     0.022 A 
HETATM   11  C1  0YB A 301      25.717 100.223  13.236  1.00 13.72     0.056 C 
HETATM   12  C8  0YB A 301      26.766  99.428  15.396  1.00 16.75     0.256 A 
HETATM   13  O9  0YB A 301      27.719  98.910  14.817  1.00 14.17    -0.267 OA
HETATM   14  N7  0YB A 301      26.675  99.372  16.776  1.00 17.17    -0.275 N 
HETATM   15  C5  0YB A 301      25.652  99.861  17.546  1.00 18.54     0.329 A 
ATOM     16  H7  0YB A 301      27.446  98.923  17.269  1.00  0.00     0.173 HD
HETATM   17  O6  0YB A 301      25.613  99.671  18.758  1.00 18.13    -0.246 OA
ENDBRANCH   6   8
ENDBRANCH   1   2
BRANCH   1  18
HETATM   18  C22 0YB A 301      19.226 103.181  16.999  1.00 24.21    -0.022 A 
HETATM   19  C23 0YB A 301      18.714 101.932  17.308  1.00 24.57     0.010 A 
HETATM   20  C24 0YB A 301      17.400 101.784  17.678  1.00 24.53     0.022 A 
HETATM   21  C25 0YB A 301      16.561 102.873  17.699  1.00 24.20     0.074 A 
HETATM   22  C26 0YB A 301      17.047 104.127  17.358  1.00 20.03     0.093 A 
HETATM   23  C27 0YB A 301      18.371 104.273  17.012  1.00 22.14     0.047 A 
BRANCH  21  24
HETATM   24  C36 0YB A 301      15.158 102.615  18.059  1.00 25.31     0.212 C 
HETATM   25  O38 0YB A 301      14.475 102.134  16.998  1.00 25.89    -0.644 OA
HETATM   26  O37 0YB A 301      14.676 102.871  19.150  1.00 28.67    -0.644 OA
ENDBRANCH  21  24
BRANCH  22  27
HETATM   27  O28 0YB A 301      16.178 105.193  17.303  1.00 22.79    -0.310 OA
BRANCH  27  28
HETATM   28  C29 0YB A 301      15.682 105.633  16.093  1.00 22.10     0.077 A 
HETATM   29  C30 0YB A 301      14.474 105.204  15.568  1.00 23.29     0.037 A 
HETATM   30  C31 0YB A 301      14.073 105.679  14.330  1.00 20.66     0.004 A 
HETATM   31  C32 0YB A 301      14.878 106.559  13.638  1.00 19.61     0.019 A 
HETATM   32  C33 0YB A 301      16.086 106.959  14.172  1.00 23.09     0.013 A 
HETATM   33 BR3  0YB A 301      17.217 108.185  13.289  1.00 21.32    -0.051 Br
HETATM   34  C34 0YB A 301      16.507 106.496  15.403  1.00 19.58     0.056 A 
ENDBRANCH  27  28
ENDBRANCH  22  27
ENDBRANCH   1  18
BRANCH   1  35
HETATM   35  C17 0YB A 301      21.118 104.729  16.332  1.00 24.12     0.049 C 
BRANCH  35  36
HETATM   36  C18 0YB A 301      20.549 105.417  15.083  1.00 29.47    -0.049 C 
HETATM   37  C19 0YB A 301      20.310 106.872  15.428  1.00 29.39     0.014 C 
HETATM   38  C20 0YB A 301      21.502 105.399  13.908  1.00 29.63     0.014 C 
HETATM   39  C21 0YB A 301      19.214 104.862  14.652  1.00 28.88     0.014 C 
ENDBRANCH  35  36
ENDBRANCH   1  35
TORSDOF 8
ENDMDL""",
            [
                "C",
                "N",
                "C",
                "C",
                "C",
                "C",
                "C",
                "N",
                "A",
                "A",
                "C",
                "A",
                "OA",
                "N",
                "A",
                "HD",
                "OA",
                "A",
                "A",
                "A",
                "A",
                "A",
                "A",
                "C",
                "OA",
                "OA",
                "OA",
                "A",
                "A",
                "A",
                "A",
                "A",
                "Br",
                "A",
                "C",
                "C",
                "C",
                "C",
                "C",
            ],
        ),
    ],
)
def test_parse_model_atoms(lines: list[str], atoms_true):
    model_lines = lines.splitlines()
    atoms, _ = PDBQTParser.parse_model(model_lines)

    assert atoms == atoms_true


@pytest.mark.parametrize(
    "lines,coords_true",
    [
        (
            """MODEL 1
REMARK VINA RESULT:      -7.5      0.000      0.000
REMARK  6 active torsions:
REMARK  status: ('A' for Active; 'I' for Inactive)
REMARK    1  A    between atoms: C08_6  and  N09_18 
REMARK    2  A    between atoms: C11_7  and  C12_8 
REMARK    3  A    between atoms: C11_7  and  S10_22 
REMARK    4  A    between atoms: C12_8  and  C13_9 
REMARK    5  A    between atoms: C13_9  and  C14_10 
REMARK    6  A    between atoms: N09_18  and  S10_22 
ROOT
HETATM    1  C01 JFX A 501      52.052  -7.816 -18.741  1.00 62.30     0.070 C 
HETATM    2  C02 JFX A 501      50.947  -7.266 -19.636  1.00 69.74     0.034 A 
HETATM    3  C05 JFX A 501      48.968  -6.525 -20.372  1.00 66.41     0.060 A 
HETATM    4  C06 JFX A 501      49.583  -7.021 -19.238  1.00 67.31     0.018 A 
HETATM    5  C07 JFX A 501      48.844  -7.172 -18.040  1.00 64.74     0.046 A 
HETATM    6  C08 JFX A 501      47.491  -6.823 -18.011  1.00 66.27     0.043 A 
HETATM    7  C18 JFX A 501      46.867  -6.322 -19.163  1.00 58.48     0.032 A 
HETATM    8  C19 JFX A 501      47.595  -6.174 -20.338  1.00 57.92     0.034 A 
HETATM    9  N03 JFX A 501      51.078  -6.931 -20.907  1.00 73.37    -0.304 N 
ATOM     10  H03 JFX A 501      51.934  -6.998 -21.457  1.00  0.00     0.186 HD
HETATM   11  N04 JFX A 501      49.895  -6.487 -21.351  1.00 70.74    -0.300 N 
ATOM     12  H04 JFX A 501      49.719  -6.164 -22.302  1.00  0.00     0.187 HD
ENDROOT
BRANCH   6  13
HETATM   13  N09 JFX A 501      46.707  -6.967 -16.798  1.00 70.80    -0.285 N 
ATOM     14  H09 JFX A 501      46.050  -7.744 -16.721  1.00  0.00     0.178 HD
BRANCH  13  15
HETATM   15  S10 JFX A 501      46.853  -5.872 -15.495  1.00 67.91     0.231 S 
HETATM   16  O16 JFX A 501      45.628  -5.094 -15.332  1.00 64.11    -0.206 OA
HETATM   17  O17 JFX A 501      47.830  -4.827 -15.809  1.00 69.85    -0.206 OA
BRANCH  15  18
HETATM   18  C11 JFX A 501      47.292  -6.762 -13.970  1.00 70.63     0.163 C 
BRANCH  18  19
HETATM   19  C12 JFX A 501      46.157  -7.680 -13.526  1.00 65.37     0.019 C 
BRANCH  19  20
HETATM   20  C13 JFX A 501      46.690  -8.854 -12.707  1.00 65.38     0.032 C 
BRANCH  20  21
HETATM   21  C14 JFX A 501      45.592  -9.858 -12.373  1.00 68.91     0.222 C 
HETATM   22  F15 JFX A 501      45.860 -10.476 -11.188  1.00 69.95    -0.255 F 
ENDBRANCH  20  21
ENDBRANCH  19  20
ENDBRANCH  18  19
ENDBRANCH  15  18
ENDBRANCH  13  15
ENDBRANCH   6  13
TORSDOF 6
ENDMDL""",
            np.array(
                [
                    [52.052, -7.816, -18.741],
                    [50.947, -7.266, -19.636],
                    [48.968, -6.525, -20.372],
                    [49.583, -7.021, -19.238],
                    [48.844, -7.172, -18.04],
                    [47.491, -6.823, -18.011],
                    [46.867, -6.322, -19.163],
                    [47.595, -6.174, -20.338],
                    [51.078, -6.931, -20.907],
                    [51.934, -6.998, -21.457],
                    [49.895, -6.487, -21.351],
                    [49.719, -6.164, -22.302],
                    [46.707, -6.967, -16.798],
                    [46.05, -7.744, -16.721],
                    [46.853, -5.872, -15.495],
                    [45.628, -5.094, -15.332],
                    [47.83, -4.827, -15.809],
                    [47.292, -6.762, -13.97],
                    [46.157, -7.68, -13.526],
                    [46.69, -8.854, -12.707],
                    [45.592, -9.858, -12.373],
                    [45.86, -10.476, -11.188],
                ]
            ),
        ),
        (
            """MODEL 1
REMARK VINA RESULT:      -6.4      0.000      0.000
REMARK  8 active torsions:
REMARK  status: ('A' for Active; 'I' for Inactive)
REMARK    1  A    between atoms: C25_8  and  C36_17 
REMARK    2  A    between atoms: C16_28  and  C22_9 
REMARK    3  A    between atoms: C29_10  and  O28_38 
REMARK    4  A    between atoms: C26_11  and  O28_38 
REMARK    5  A    between atoms: C10_22  and  N4_32 
REMARK    6  A    between atoms: C16_28  and  C17_27 
REMARK    7  A    between atoms: C17_27  and  C18_29 
REMARK    8  A    between atoms: C16_28  and  N14_33 
ROOT
HETATM    1  C16 0YB A 301      20.690 103.282  16.649  1.00 23.40     0.256 C 
ENDROOT
BRANCH   1   2
HETATM    2  N14 0YB A 301      21.630 102.590  17.558  1.00 22.62     0.099 N 
HETATM    3  C13 0YB A 301      22.119 103.416  18.677  1.00 23.46     0.273 C 
HETATM    4  C12 0YB A 301      22.924 102.561  19.646  1.00 22.53     0.032 C 
HETATM    5  C11 0YB A 301      24.092 101.878  18.931  1.00 20.23     0.029 C 
HETATM    6  C10 0YB A 301      23.644 101.183  17.634  1.00 19.27     0.122 C 
HETATM    7  C15 0YB A 301      22.783 102.087  16.806  1.00 19.73     0.298 C 
BRANCH   6   8
HETATM    8  N4  0YB A 301      24.703 100.580  16.853  1.00 17.74    -0.303 N 
HETATM    9  C3  0YB A 301      24.708 100.670  15.478  1.00 16.77     0.108 A 
HETATM   10  C2  0YB A 301      25.688 100.139  14.725  1.00 14.10     0.022 A 
HETATM   11  C1  0YB A 301      25.717 100.223  13.236  1.00 13.72     0.056 C 
HETATM   12  C8  0YB A 301      26.766  99.428  15.396  1.00 16.75     0.256 A 
HETATM   13  O9  0YB A 301      27.719  98.910  14.817  1.00 14.17    -0.267 OA
HETATM   14  N7  0YB A 301      26.675  99.372  16.776  1.00 17.17    -0.275 N 
HETATM   15  C5  0YB A 301      25.652  99.861  17.546  1.00 18.54     0.329 A 
ATOM     16  H7  0YB A 301      27.446  98.923  17.269  1.00  0.00     0.173 HD
HETATM   17  O6  0YB A 301      25.613  99.671  18.758  1.00 18.13    -0.246 OA
ENDBRANCH   6   8
ENDBRANCH   1   2
BRANCH   1  18
HETATM   18  C22 0YB A 301      19.226 103.181  16.999  1.00 24.21    -0.022 A 
HETATM   19  C23 0YB A 301      18.714 101.932  17.308  1.00 24.57     0.010 A 
HETATM   20  C24 0YB A 301      17.400 101.784  17.678  1.00 24.53     0.022 A 
HETATM   21  C25 0YB A 301      16.561 102.873  17.699  1.00 24.20     0.074 A 
HETATM   22  C26 0YB A 301      17.047 104.127  17.358  1.00 20.03     0.093 A 
HETATM   23  C27 0YB A 301      18.371 104.273  17.012  1.00 22.14     0.047 A 
BRANCH  21  24
HETATM   24  C36 0YB A 301      15.158 102.615  18.059  1.00 25.31     0.212 C 
HETATM   25  O38 0YB A 301      14.475 102.134  16.998  1.00 25.89    -0.644 OA
HETATM   26  O37 0YB A 301      14.676 102.871  19.150  1.00 28.67    -0.644 OA
ENDBRANCH  21  24
BRANCH  22  27
HETATM   27  O28 0YB A 301      16.178 105.193  17.303  1.00 22.79    -0.310 OA
BRANCH  27  28
HETATM   28  C29 0YB A 301      15.682 105.633  16.093  1.00 22.10     0.077 A 
HETATM   29  C30 0YB A 301      14.474 105.204  15.568  1.00 23.29     0.037 A 
HETATM   30  C31 0YB A 301      14.073 105.679  14.330  1.00 20.66     0.004 A 
HETATM   31  C32 0YB A 301      14.878 106.559  13.638  1.00 19.61     0.019 A 
HETATM   32  C33 0YB A 301      16.086 106.959  14.172  1.00 23.09     0.013 A 
HETATM   33 BR3  0YB A 301      17.217 108.185  13.289  1.00 21.32    -0.051 Br
HETATM   34  C34 0YB A 301      16.507 106.496  15.403  1.00 19.58     0.056 A 
ENDBRANCH  27  28
ENDBRANCH  22  27
ENDBRANCH   1  18
BRANCH   1  35
HETATM   35  C17 0YB A 301      21.118 104.729  16.332  1.00 24.12     0.049 C 
BRANCH  35  36
HETATM   36  C18 0YB A 301      20.549 105.417  15.083  1.00 29.47    -0.049 C 
HETATM   37  C19 0YB A 301      20.310 106.872  15.428  1.00 29.39     0.014 C 
HETATM   38  C20 0YB A 301      21.502 105.399  13.908  1.00 29.63     0.014 C 
HETATM   39  C21 0YB A 301      19.214 104.862  14.652  1.00 28.88     0.014 C 
ENDBRANCH  35  36
ENDBRANCH   1  35
TORSDOF 8
ENDMDL""",
            np.array(
                [
                    [20.69, 103.282, 16.649],
                    [21.63, 102.59, 17.558],
                    [22.119, 103.416, 18.677],
                    [22.924, 102.561, 19.646],
                    [24.092, 101.878, 18.931],
                    [23.644, 101.183, 17.634],
                    [22.783, 102.087, 16.806],
                    [24.703, 100.58, 16.853],
                    [24.708, 100.67, 15.478],
                    [25.688, 100.139, 14.725],
                    [25.717, 100.223, 13.236],
                    [26.766, 99.428, 15.396],
                    [27.719, 98.91, 14.817],
                    [26.675, 99.372, 16.776],
                    [25.652, 99.861, 17.546],
                    [27.446, 98.923, 17.269],
                    [25.613, 99.671, 18.758],
                    [19.226, 103.181, 16.999],
                    [18.714, 101.932, 17.308],
                    [17.4, 101.784, 17.678],
                    [16.561, 102.873, 17.699],
                    [17.047, 104.127, 17.358],
                    [18.371, 104.273, 17.012],
                    [15.158, 102.615, 18.059],
                    [14.475, 102.134, 16.998],
                    [14.676, 102.871, 19.15],
                    [16.178, 105.193, 17.303],
                    [15.682, 105.633, 16.093],
                    [14.474, 105.204, 15.568],
                    [14.073, 105.679, 14.33],
                    [14.878, 106.559, 13.638],
                    [16.086, 106.959, 14.172],
                    [17.217, 108.185, 13.289],
                    [16.507, 106.496, 15.403],
                    [21.118, 104.729, 16.332],
                    [20.549, 105.417, 15.083],
                    [20.310, 106.872, 15.428],
                    [21.502, 105.399, 13.908],
                    [19.214, 104.862, 14.652],
                ]
            ),
        ),
    ],
)
def test_parse_model_coords(lines: list[str], coords_true):
    model_lines = lines.splitlines()
    _, coords = PDBQTParser.parse_model(model_lines)

    np.testing.assert_array_almost_equal(coords, coords_true)


@pytest.mark.parametrize(
    "filepath,n_models_true,atoms_true",
    [
        (
            DATA_DIR / "jfx.pdbqt",
            9,
            [
                "C",
                "C",
                "C",
                "C",
                "C",
                "C",
                "C",
                "C",
                "N",
                "H",
                "N",
                "H",
                "N",
                "H",
                "S",
                "O",
                "O",
                "C",
                "C",
                "C",
                "C",
                "F",
            ],
        ),
        (
            DATA_DIR / "0yb.pdbqt",
            2,
            [
                "C",
                "N",
                "C",
                "C",
                "C",
                "C",
                "C",
                "N",
                "C",
                "C",
                "C",
                "C",
                "O",
                "N",
                "C",
                "H",
                "O",
                "C",
                "C",
                "C",
                "C",
                "C",
                "C",
                "C",
                "O",
                "O",
                "O",
                "C",
                "C",
                "C",
                "C",
                "C",
                "Br",
                "C",
                "C",
                "C",
                "C",
                "C",
                "C",
            ],
        ),
    ],
)
def test_parse_atoms_and_shape(filepath, n_models_true, atoms_true):
    atoms, C = PDBQTParser.parse(filepath)

    assert atoms == atoms_true
    assert C.shape == (n_models_true, len(atoms_true), 3)
