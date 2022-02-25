import numpy as np
import pytest

from pyscreener.utils import autobox


@pytest.fixture(
    params=[
        (
            """ATOM      1  N   GLY A  34     -20.681  20.730 -35.850  1.00  0.00           N  
ATOM      2  H   GLY A  34     -20.324  21.878 -35.915  1.00  0.00           H  
ATOM      3  H2  GLY A  34     -19.812  20.282 -35.161  1.00  0.00           H  
ATOM      4  H3  GLY A  34     -20.330  20.399 -36.944  1.00  0.00           H  
ATOM      5  CA  GLY A  34     -21.957  20.459 -35.211  1.00  0.00           C  
ATOM      6  HA2 GLY A  34     -21.838  21.623 -34.956  1.00  0.00           H  
ATOM      7  HA3 GLY A  34     -23.008  20.482 -34.641  1.00  0.00           H  
ATOM      8  C   GLY A  34     -22.617  19.194 -35.727  1.00  0.00           C  
ATOM      9  O   GLY A  34     -23.111  18.385 -34.948  1.00  0.00           O  
ATOM     10  N   ALA A  35     -22.630  19.026 -37.046  1.00  0.00           N  
ATOM     11  H   ALA A  35     -22.652  19.998 -37.743  1.00  0.00           H  
ATOM     12  CA  ALA A  35     -23.227  17.847 -37.662  1.00  0.00           C  
ATOM     13  HA  ALA A  35     -24.411  17.865 -37.484  1.00  0.00           H  
ATOM     14  C   ALA A  35     -22.484  16.575 -37.266  1.00  0.00           C  
ATOM     15  O   ALA A  35     -23.094  15.521 -37.079  1.00  0.00           O  
ATOM     16  CB  ALA A  35     -23.259  17.999 -39.183  1.00  0.00           C  
ATOM     17  HB1 ALA A  35     -22.254  18.213 -39.807  1.00  0.00           H  
ATOM     18  HB2 ALA A  35     -23.979  18.863 -39.604  1.00  0.00           H  
ATOM     19  HB3 ALA A  35     -23.693  17.078 -39.820  1.00  0.00           H  
ATOM     20  N   ALA A  36     -21.167  16.679 -37.115  1.00  0.00           N  
ATOM     21  H   ALA A  36     -20.548  17.593 -37.560  1.00  0.00           H  
ATOM     22  CA  ALA A  36     -20.349  15.533 -36.742  1.00  0.00           C  
ATOM     23  HA  ALA A  36     -20.504  14.695 -37.580  1.00  0.00           H  
ATOM     24  C   ALA A  36     -20.617  15.112 -35.308  1.00  0.00           C  
ATOM     25  O   ALA A  36     -20.604  13.924 -34.987  1.00  0.00           O  
ATOM     26  CB  ALA A  36     -18.879  15.850 -36.935  1.00  0.00           C  
ATOM     27  HB1 ALA A  36     -18.173  16.457 -36.181  1.00  0.00           H  
ATOM     28  HB2 ALA A  36     -18.218  14.853 -37.047  1.00  0.00           H  
ATOM     29  HB3 ALA A  36     -18.553  16.374 -37.965  1.00  0.00           H  
ATOM     30  N   ALA A  37     -20.859  16.084 -34.437  1.00  0.00           N  
ATOM     31  H   ALA A  37     -20.261  17.093 -34.646  1.00  0.00           H  
ATOM     32  CA  ALA A  37     -21.179  15.783 -33.044  1.00  0.00           C  
ATOM     33  HA  ALA A  37     -20.343  14.987 -32.729  1.00  0.00           H  
ATOM     34  C   ALA A  37     -22.522  15.075 -32.978  1.00  0.00           C  
ATOM     35  O   ALA A  37     -22.714  14.122 -32.222  1.00  0.00           O  
ATOM     36  CB  ALA A  37     -21.205  17.057 -32.212  1.00  0.00           C  
ATOM     37  HB1 ALA A  37     -22.250  17.339 -31.725  1.00  0.00           H  
ATOM     38  HB2 ALA A  37     -20.458  16.761 -31.336  1.00  0.00           H  
ATOM     39  HB3 ALA A  37     -20.670  18.054 -32.600  1.00  0.00           H  
ATOM     40  N   LEU A  38     -23.450  15.556 -33.792  1.00  0.00           N  
ATOM     41  H   LEU A  38     -23.343  16.606 -34.321  1.00  0.00           H  
ATOM     42  CA  LEU A  38     -24.785  14.983 -33.853  1.00  0.00           C  
ATOM     43  HA  LEU A  38     -25.189  14.914 -32.735  1.00  0.00           H  
ATOM     44  C   LEU A  38     -24.756  13.558 -34.393  1.00  0.00           C  
ATOM     45  O   LEU A  38     -25.307  12.650 -33.778  1.00  0.00           O  
ATOM     46  CB  LEU A  38     -25.692  15.859 -34.709  1.00  0.00           C  
ATOM     47  HB2 LEU A  38     -25.814  17.005 -34.397  1.00  0.00           H  
ATOM     48  HB3 LEU A  38     -25.373  15.888 -35.863  1.00  0.00           H  
ATOM     49  CG  LEU A  38     -27.107  15.326 -34.925  1.00  0.00           C  
ATOM     50  HG  LEU A  38     -27.327  14.419 -35.680  1.00  0.00           H  
ATOM     51  CD1 LEU A  38     -27.740  14.903 -33.602  1.00  0.00           C  
ATOM     52 HD11 LEU A  38     -27.466  13.865 -33.095  1.00  0.00           H  
ATOM     53 HD12 LEU A  38     -28.943  14.823 -33.700  1.00  0.00           H  
ATOM     54 HD13 LEU A  38     -27.733  15.716 -32.710  1.00  0.00           H  
ATOM     55  CD2 LEU A  38     -27.940  16.395 -35.595  1.00  0.00           C  
ATOM     56 HD21 LEU A  38     -28.264  16.055 -36.704  1.00  0.00           H  
ATOM     57 HD22 LEU A  38     -29.028  16.654 -35.152  1.00  0.00           H  
ATOM     58 HD23 LEU A  38     -27.643  17.539 -35.842  1.00  0.00           H  
ATOM     59  N   VAL A  39     -24.111  13.363 -35.539  1.00  0.00           N  
ATOM     60  H   VAL A  39     -23.154  14.014 -35.774  1.00  0.00           H  
ATOM     61  CA  VAL A  39     -24.027  12.037 -36.147  1.00  0.00           C  
ATOM     62  HA  VAL A  39     -25.125  11.609 -36.360  1.00  0.00           H  
ATOM     63  C   VAL A  39     -23.284  11.055 -35.250  1.00  0.00           C  
ATOM     64  O   VAL A  39     -23.807   9.991 -34.900  1.00  0.00           O  
ATOM     65  CB  VAL A  39     -23.340  12.098 -37.523  1.00  0.00           C  
ATOM     66  HB  VAL A  39     -22.282  12.624 -37.692  1.00  0.00           H  
ATOM     67  CG1 VAL A  39     -23.012  10.698 -38.015  1.00  0.00           C  
ATOM     68 HG11 VAL A  39     -21.851  10.375 -38.046  1.00  0.00           H  
ATOM     69 HG12 VAL A  39     -23.691   9.784 -37.645  1.00  0.00           H  
ATOM     70 HG13 VAL A  39     -23.206  10.511 -39.186  1.00  0.00           H  
ATOM     71  CG2 VAL A  39     -24.232  12.811 -38.518  1.00  0.00           C  
ATOM     72 HG21 VAL A  39     -23.653  13.193 -39.500  1.00  0.00           H  
ATOM     73 HG22 VAL A  39     -24.930  13.773 -38.365  1.00  0.00           H  
ATOM     74 HG23 VAL A  39     -25.089  12.134 -39.024  1.00  0.00           H  
ATOM     75  N   GLY A  40     -22.069  11.417 -34.859  1.00  0.00           N  
ATOM     76  H   GLY A  40     -21.389  12.296 -35.254  1.00  0.00           H  
ATOM     77  CA  GLY A  40     -21.281  10.581 -33.966  1.00  0.00           C  
ATOM     78  HA2 GLY A  40     -20.924   9.657 -34.637  1.00  0.00           H  
ATOM     79  HA3 GLY A  40     -20.182  10.983 -33.700  1.00  0.00           H  
ATOM     80  C   GLY A  40     -21.989  10.324 -32.649  1.00  0.00           C  
ATOM     81  O   GLY A  40     -21.926   9.220 -32.107  1.00  0.00           O  
ATOM     82  N   GLY A  41     -22.670  11.347 -32.142  1.00  0.00           N  
ATOM     83  H   GLY A  41     -21.827  12.188 -32.137  1.00  0.00           H  
ATOM     84  CA  GLY A  41     -23.358  11.248 -30.871  1.00  0.00           C  
ATOM     85  HA2 GLY A  41     -22.609  10.897 -30.023  1.00  0.00           H  
ATOM     86  HA3 GLY A  41     -23.892  12.310 -30.836  1.00  0.00           H  
ATOM     87  C   GLY A  41     -24.517  10.275 -30.917  1.00  0.00           C  
ATOM     88  O   GLY A  41     -24.666   9.428 -30.041  1.00  0.00           O  
ATOM     89  N   VAL A  42     -25.339  10.394 -31.952  1.00  0.00           N  
ATOM     90  H   VAL A  42     -25.294  11.402 -32.555  1.00  0.00           H  
ATOM     91  CA  VAL A  42     -26.482   9.507 -32.107  1.00  0.00           C  
ATOM     92  HA  VAL A  42     -27.188   9.574 -31.148  1.00  0.00           H  
ATOM     93  C   VAL A  42     -26.003   8.071 -32.266  1.00  0.00           C  
ATOM     94  O   VAL A  42     -26.536   7.151 -31.638  1.00  0.00           O  
ATOM     95  CB  VAL A  42     -27.361   9.933 -33.311  1.00  0.00           C  
ATOM     96  HB  VAL A  42     -26.790  10.016 -34.360  1.00  0.00           H  
ATOM     97  CG1 VAL A  42     -28.311   8.834 -33.704  1.00  0.00           C  
ATOM     98 HG11 VAL A  42     -29.054   8.271 -32.943  1.00  0.00           H  
ATOM     99 HG12 VAL A  42     -27.954   7.855 -34.303  1.00  0.00           H  
ATOM    100 HG13 VAL A  42     -29.118   9.174 -34.530  1.00  0.00           H  
ATOM    101  CG2 VAL A  42     -28.123  11.205 -32.984  1.00  0.00           C  
ATOM    102 HG21 VAL A  42     -29.247  10.806 -32.786  1.00  0.00           H  
ATOM    103 HG22 VAL A  42     -28.042  11.992 -32.091  1.00  0.00           H  
ATOM    104 HG23 VAL A  42     -28.421  11.717 -34.043  1.00  0.00           H  
ATOM    105  N   LEU A  43     -24.971   7.876 -33.078  1.00  0.00           N  
ATOM    106  H   LEU A  43     -24.288   8.691 -33.591  1.00  0.00           H  
ATOM    107  CA  LEU A  43     -24.416   6.541 -33.268  1.00  0.00           C  
ATOM    108  HA  LEU A  43     -25.336   5.824 -33.532  1.00  0.00           H  
ATOM    109  C   LEU A  43     -23.861   6.000 -31.961  1.00  0.00           C  
ATOM    110  O   LEU A  43     -24.075   4.829 -31.620  1.00  0.00           O  
ATOM    111  CB  LEU A  43     -23.325   6.553 -34.334  1.00  0.00           C  
ATOM    112  HB2 LEU A  43     -22.406   6.763 -33.600  1.00  0.00           H  
ATOM    113  HB3 LEU A  43     -22.610   5.686 -34.775  1.00  0.00           H  
ATOM    114  CG  LEU A  43     -23.840   6.689 -35.762  1.00  0.00           C  
ATOM    115  HG  LEU A  43     -24.534   7.615 -36.075  1.00  0.00           H  
ATOM    116  CD1 LEU A  43     -22.671   6.753 -36.731  1.00  0.00           C  
ATOM    117 HD11 LEU A  43     -23.082   7.183 -37.776  1.00  0.00           H  
ATOM    118 HD12 LEU A  43     -22.268   5.695 -37.131  1.00  0.00           H  
ATOM    119 HD13 LEU A  43     -21.582   7.265 -36.655  1.00  0.00           H  
ATOM    120  CD2 LEU A  43     -24.773   5.542 -36.125  1.00  0.00           C  
ATOM    121 HD21 LEU A  43     -24.805   4.493 -35.535  1.00  0.00           H  
ATOM    122 HD22 LEU A  43     -24.651   5.054 -37.219  1.00  0.00           H  
ATOM    123 HD23 LEU A  43     -25.929   5.821 -36.291  1.00  0.00           H  
ATOM    124  N   LEU A  44     -23.166   6.856 -31.217  1.00  0.00           N  
ATOM    125  H   LEU A  44     -22.994   7.955 -31.593  1.00  0.00           H  
ATOM    126  CA  LEU A  44     -22.591   6.429 -29.955  1.00  0.00           C  
ATOM    127  HA  LEU A  44     -21.999   5.463 -30.346  1.00  0.00           H""",
            0,
            11,
        ),
        (
            """ATOM   5816  HB2 ALA A 462     -27.968 -14.918 -35.890  1.00  0.00           H  
ATOM   5817  HB3 ALA A 462     -27.536 -13.525 -34.715  1.00  0.00           H  
ATOM   5818  OXT ALA A 462     -29.496 -16.673 -34.979  1.00  0.00           O  
TER    5819      ALA A 462
HETATM 5820  CAB AQD A1201     -16.765  11.793 -16.344  1.00  0.00           C  
HETATM 5821  OAR AQD A1201     -16.479  13.184 -16.085  1.00  0.00           O  
HETATM 5822  CAW AQD A1201     -15.615  13.449 -14.994  1.00  0.00           C  
HETATM 5823  CAK AQD A1201     -15.231  12.434 -14.098  1.00  0.00           C  
""",
            4,
            0,
        ),
        (
            """HETATM 5820  CAB AQD A1201     -16.765  11.793 -16.344  1.00  0.00           C  
HETATM 5821  OAR AQD A1201     -16.479  13.184 -16.085  1.00  0.00           O  
HETATM 5822  CAW AQD A1201     -15.615  13.449 -14.994  1.00  0.00           C  
HETATM 5823  CAK AQD A1201     -15.231  12.434 -14.098  1.00  0.00           C  
HETATM 5824  CAV AQD A1201     -14.376  12.763 -13.031  1.00  0.00           C  
HETATM 5825  NAP AQD A1201     -13.960  11.763 -12.096  1.00  0.00           N  
HETATM 5826  CAA AQD A1201     -14.651  10.489 -12.085  1.00  0.00           C  
HETATM 5827  CAU AQD A1201     -13.891  14.071 -12.887  1.00  0.00           C  
HETATM 5828 CLA  AQD A1201     -12.777  14.490 -11.593  1.00  0.00          Cl  
HETATM 5829  CAL AQD A1201     -14.265  15.058 -13.780  1.00  0.00           C  
HETATM 5830  CAX AQD A1201     -15.118  14.748 -14.834  1.00  0.00           C  
HETATM 5831  CAS AQD A1201     -15.526  15.904 -15.744  1.00  0.00           C  
HETATM 5832  OAD AQD A1201     -14.929  16.950 -15.642  1.00  0.00           O  
HETATM 5833  NAQ AQD A1201     -16.638  15.809 -16.648  1.00  0.00           N  
HETATM 5834  CAZ AQD A1201     -16.997  16.959 -17.469  1.00  0.00           C  
HETATM 5835  CAY AQD A1201     -17.416  16.578 -18.817  1.00  0.00           C  
HETATM 5836  CAC AQD A1201     -16.699  15.290 -19.320  1.00  0.00           C  
HETATM 5837  CAM AQD A1201     -18.166  17.653 -16.797  1.00  0.00           C  
HETATM 5838  CAN AQD A1201     -19.349  17.390 -17.881  1.00  0.00           C  
HETATM 5839  NBA AQD A1201     -18.795  16.355 -18.773  1.00  0.00           N  
HETATM 5840  CAO AQD A1201     -19.331  16.401 -20.098  1.00  0.00           C  
HETATM 5841  CAT AQD A1201     -20.852  16.758 -20.138  1.00  0.00           C  
HETATM 5842  CAI AQD A1201     -21.259  18.001 -20.632  1.00  0.00           C  
HETATM 5843  CAG AQD A1201     -22.616  18.335 -20.665  1.00  0.00           C  
HETATM 5844  CAF AQD A1201     -23.565  17.421 -20.210  1.00  0.00           C  
HETATM 5845  CAH AQD A1201     -23.166  16.184 -19.721  1.00  0.00           C  
HETATM 5846  CAJ AQD A1201     -21.806  15.847 -19.682  1.00  0.00           C""",
            27,
            0,
        ),
    ]
)
def pdbfile(request, tmp_path):
    pdbblock, n_hetatm_lines, n_residue_lines = request.param
    filepath = tmp_path / "tmp.pdb"
    with open(filepath, "w") as fid:
        fid.write(pdbblock)

    return filepath, n_hetatm_lines, n_residue_lines


@pytest.fixture(
    params=[
        (
            "HETATM 5820  CAB AQD A1201     -16.765  11.793 -16.344  1.00  0.00           C  ",
            (-16.765, 11.793, -16.344),
        ),
        (
            "HETATM 5821  OAR AQD A1201     -16.479  13.184 -16.085  1.00  0.00           O  ",
            (-16.479, 13.184, -16.085),
        ),
        (
            "HETATM 5822  CAW AQD A1201     -15.615  13.449 -14.994  1.00  0.00           C  ",
            (-15.615, 13.449, -14.994),
        ),
        (
            "HETATM 5823  CAK AQD A1201     -15.231  12.434 -14.098  1.00  0.00           C  ",
            (-15.231, 12.434, -14.098),
        ),
        (
            "HETATM 5824  CAV AQD A1201     -14.376  12.763 -13.031  1.00  0.00           C  ",
            (-14.376, 12.763, -13.031),
        ),
        (
            "HETATM 5825  NAP AQD A1201     -13.960  11.763 -12.096  1.00  0.00           N  ",
            (-13.96, 11.763, -12.096),
        ),
        (
            "HETATM 5826  CAA AQD A1201     -14.651  10.489 -12.085  1.00  0.00           C  ",
            (-14.651, 10.489, -12.085),
        ),
        (
            "HETATM 5827  CAU AQD A1201     -13.891  14.071 -12.887  1.00  0.00           C  ",
            (-13.891, 14.071, -12.887),
        ),
        (
            "HETATM 5828 CLA  AQD A1201     -12.777  14.490 -11.593  1.00  0.00           Cl  ",
            (-12.777, 14.49, -11.593),
        ),
        (
            "HETATM 5829  CAL AQD A1201     -14.265  15.058 -13.780  1.00  0.00           C  ",
            (-14.265, 15.058, -13.78),
        ),
        (
            "HETATM 5830  CAX AQD A1201     -15.118  14.748 -14.834  1.00  0.00           C  ",
            (-15.118, 14.748, -14.834),
        ),
        (
            "HETATM 5831  CAS AQD A1201     -15.526  15.904 -15.744  1.00  0.00           C  ",
            (-15.526, 15.904, -15.744),
        ),
        (
            "HETATM 5832  OAD AQD A1201     -14.929  16.950 -15.642  1.00  0.00           O  ",
            (-14.929, 16.95, -15.642),
        ),
        (
            "HETATM 5833  NAQ AQD A1201     -16.638  15.809 -16.648  1.00  0.00           N  ",
            (-16.638, 15.809, -16.648),
        ),
        (
            "HETATM 5834  CAZ AQD A1201     -16.997  16.959 -17.469  1.00  0.00           C  ",
            (-16.997, 16.959, -17.469),
        ),
        (
            "HETATM 5835  CAY AQD A1201     -17.416  16.578 -18.817  1.00  0.00           C",
            (-17.416, 16.578, -18.817),
        ),
    ]
)
def line(request):
    return request.param


@pytest.fixture(params=[1.0, 2.0, 4.0, 8.0])
def length(request):
    return request.param


@pytest.fixture(params=[1e-6, 0.1, 1.0, 2.0])  # small amount for comparison stability
def buffer(request):
    return request.param


def test_extract_residue_lines(pdbfile):
    pdb_filepath, n_hetatm_lines, n_residue_lines = pdbfile
    lines = autobox.extract_residues_lines(pdb_filepath, range(0, 1000))

    assert len(lines) == n_residue_lines


def test_extract_hetatm_lines(pdbfile):
    pdb_filepath, n_hetatm_lines, n_residue_lines = pdbfile
    lines = autobox.extract_hetatm_lines(pdb_filepath)

    assert len(lines) == n_hetatm_lines


def test_parse_line(line):
    line, coords = line
    x, y, z = autobox.parse_coordinates(line)

    assert (x, y, z) == coords


def sample_spherical(npoints, ndim=3):
    vec = np.random.randn(ndim, npoints)
    vec /= np.linalg.norm(vec, axis=0)

    return vec


def coords(length):
    N = 100
    return length * sample_spherical(N)


def test_minimum_bounding_box(length, buffer):
    X = coords(length).T
    center, size = np.array(autobox.minimum_bounding_box(X, buffer))
    mins = center - (size + buffer)
    maxs = center + (size + buffer)

    assert all([[min_c <= c <= max_c for c, min_c, max_c in zip(xyz, mins, maxs)] for xyz in X])
