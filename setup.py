#!/usr/bin/env python
'''PNPy setup.py file'''

# import os
# import shutil
try:
    from setuptools import setup
except ImportError as ie:
    raise ie, 'please install setuptools'

setup(
    name = "PNPy",
    version = "0.0.2",
    maintainer = "Carl H Lubba",
    maintainer_email = 'c.lubba15@imperial.ac.uk',
    packages = ['PNPy'],
    url='https://github.com/caaarl/PNPy',
    download_url = 'https://github.com/caaarl/PNPy/archive/0.0.2.tar.gz',
    license='LICENSE',
    description='A peripheral nerve simulator built on NEURON',
    #long_description=,
    install_requires = ['numpy', 'scipy', 'matplotlib', 'neuron'],
    provides = ['PNPy'],
    )
