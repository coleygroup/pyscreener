"""
pyscreener
pythonic interface to virtual screening software
"""
import sys
from setuptools import setup, find_packages
import versioneer

short_description = __doc__.split("\n")

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:])


setup(
    name="pyscreener",
    author="david graff",
    author_email="deg711@g.harvard.edu",
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license="MIT",
    packages=find_packages(exclude=["pyscreener.dft", "pyscreener.md", "pyscreener.docking.dock.scripts"]),
    include_package_data=True,
    setup_requires=[],
    url="https://github.com/coleygroup/pyscreener",
    platforms=["Linux", "Mac OS-X", "Unix"],
    python_requires=">=3.7",
    install_requires=[
        "configargparse",
        "h5py",
        "numpy",
        "ray[default]",
        "pandas",
        "pdbfixer @ git+https://github.com/openmm/pdbfixer.git",
        "seaborn",
        "scikit_learn",
        "scipy",
        "tqdm",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Framework :: Pytest",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Typing :: Typed"
    ]
)
