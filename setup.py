#! /usr/bin/env python
#
# Copyright (C) Hideto Mori

DESCRIPTION = "QUEEN (a Python module to universally program, QUinE, and Edit Nucleotide sequences)"
LONG_DESCRIPTION = """\
QUEEN (a Python module to universally program, quine, and edit nucleotide sequences) enables universal editing and annotation of double-stranded (ds)DNA objects for designing molecular cloning and simulation of genome editing, DNA recombination, and genetic circuits. DNA parts information can be imported from external annotated DNA files (GenBank and FASTA format). Output file (GenBank format) encodes the full information of the constructed DNA and enables to produce a quine code that self-reproduces the output file itself. In QUEEN, all of the manipulations required in DNA engineering are covered by four common operational functions, “cut,” “end-modification,” “flip” and “join,” which can represent all of the standard DNA cloning processes of today and more like super functions, “edit sequence” and “edit feature,” along with a simple set of analytical and visualization functions. A new DNA can be designed by programming a Python script or by an interactive interpreter Jupyter Notebook. The designed DNA product can be output as a GenBank format file that encodes a quine code to reproduce the file itself. The “quinable” feature of QUEEN guarantees to produce annotated DNA information whose design and construction processes are fully transparent, reproducible, and modifiable by the community.
"""

DISTNAME         = 'QUEEN'
MAINTAINER       = 'Hideto Mori'
MAINTAINER_EMAIL = 'morityunåsfc.keio.ac.jp/hidto7592ågmail.com'
URL = 'https://github.com/yachielab/QUEEN'
LICENSE = 'GNU General Public License v3.0'
DOWNLOAD_URL = 'https://github.com/yachielab/QUEEN'
VERSION = '1.0.0'
PYTHON_REQUIRES = ">=3.7"

INSTALL_REQUIRES = [
    'numpy>=1.2',
    'biopython>=1.78',
    'matplotlib>=3.4',
    'regex>=2.5',
    'graphviz>=0.17'
]


PACKAGES = [
    'QUEEN'
]

CLASSIFIERS = [
    'Intended Audience :: Science/Research',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'License :: OSI Approved :: GNU General Public License v3.0',
    'Topic :: Bioinformatics',
    'Operating System :: OS Independent',
]


if __name__ == "__main__":

    from setuptools import setup

    import sys
    if sys.version_info[:2] < (3, 7):
        raise RuntimeError("QUEEN requires python >= 3.7.")

    setup(
        name=DISTNAME,
        author=MAINTAINER,
        author_email=MAINTAINER_EMAIL,
        maintainer=MAINTAINER,
        maintainer_email=MAINTAINER_EMAIL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        license=LICENSE,
        url=URL,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        python_requires=PYTHON_REQUIRES,
        install_requires=INSTALL_REQUIRES,
        packages=PACKAGES,
        classifiers=CLASSIFIERS
    )
