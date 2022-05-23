#! /usr/bin/env python
#
# Copyright (C) Hideto Mori

DESCRIPTION      = "QUEEN (a Python module to universally program, QUinE, and Edit Nucleotide sequences)"
DISTNAME         = 'python-queen'
MAINTAINER       = 'Hideto Mori'
MAINTAINER_EMAIL = 'hidto7592@gmail.com'
URL              = 'https://github.com/yachielab/QUEEN'
LICENSE          = 'MIT'
DOWNLOAD_URL     = 'https://github.com/yachielab/QUEEN'
VERSION          = '1.0.1'
PYTHON_REQUIRES  = ">=3.7"

INSTALL_REQUIRES = [
    'numpy>=1.2',
    'biopython>=1.78',
    'matplotlib>=3.2',
    'requests~=2.23.0',
    'regex>=2.5',
    'graphviz==0.17', 
    'beautifulsoup4>=4.4'
]

PACKAGES = [
    'QUEEN'
]

CLASSIFIERS = [
    'Intended Audience :: Science/Research',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'License :: OSI Approved :: MIT License',
]

with open('README.md', 'r', encoding='utf-8') as fp:
    readme = fp.read()
LONG_DESCRIPTION = readme
LONG_DESCRIPTION_CONTENT_TYPE = 'text/markdown'

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
        long_description_content_type=LONG_DESCRIPTION_CONTENT_TYPE,
        license=LICENSE,
        url=URL,
        version=VERSION,
        download_url=DOWNLOAD_URL,
        python_requires=PYTHON_REQUIRES,
        install_requires=INSTALL_REQUIRES,
        packages=PACKAGES,
        classifiers=CLASSIFIERS
    )
