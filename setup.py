"""
pyscreener
pythonic interface to virtual screening software
"""
import sys
from setuptools import setup, find_packages
import versioneer

short_description = __doc__.split("\n")

# from https://github.com/pytest-dev/pytest-runner#conditional-requirement
needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)

try:
    with open("README.md", "r") as handle:
        long_description = handle.read()
except:
    long_description = "\n".join(short_description[2:])


setup(
    name='pyscreener',
    author='david graff',
    author_email='deg711@g.harvard.edu',
    description=short_description[0],
    long_description=long_description,
    long_description_content_type="text/markdown",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='MIT',

    # Which Python importable modules should be included when your package is 
    # installed. Handled automatically by setuptools. Use 'exclude' to prevent 
    # some specific subpackage(s) from being added, if needed
    packages=find_packages(
        exclude=[
            'pyscreener.dft',
            'pyscreener.md'
            'pyscreener.docking.dock.scripts',
        ]
    ),

    # Optional include package data to ship with your package
    # Customize MANIFEST.in if the general case does not suit your needs
    # Comment out this line to prevent the files from being packaged with your software
    include_package_data=True,
    setup_requires=[],

    # url='http://www.my_package.com',  # Website
    platforms=['Linux', 'Mac OS-X', 'Unix'],
    python_requires=">=3.7",
    install_requires=[
        'configargparse',
        'h5py',
        'numpy',
        'ray[default]',
        'pandas',
        'pdbfixer @ git+https://github.com/openmm/pdbfixer.git',
        'seaborn',
        'scikit_learn',
        'scipy',
        'tqdm'
    ],

    # Manual control if final package is compressible or not, set False to prevent the .egg from being made
    # zip_safe=False,

)
