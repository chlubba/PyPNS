#!/usr/bin/env python
'''PyPNS setup.py file'''

# import os
# import shutil
try:
    from setuptools import setup
except ImportError as ie:
    raise ie, 'please install setuptools'

setup(
    name = "PyPNS",
    version = "0.0.4",
    maintainer = "Carl H Lubba",
    maintainer_email = 'c.lubba15@imperial.ac.uk',
    packages = ['PyPNS'],
    url='https://github.com/caaarl/PyPNS',
    download_url = 'https://github.com/caaarl/PyPNS/archive/0.0.4.tar.gz',
    license='LICENSE',
    description='A peripheral nerve simulator built on NEURON',
    #long_description=,
    install_requires = ['numpy', 'scipy', 'matplotlib'], # , 'neuron' exluded because we need a higher version.
    provides = ['PyPNS'],
    )
