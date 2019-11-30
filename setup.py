#!/usr/bin/env python
'''PyPNS setup.py file'''

# import os
# import shutil
try:
    from setuptools import setup
except ImportError as ie:
    print('please install setuptools')
    raise ie


setup(
    name = "PyPNS",
    version = "0.1.0",
    maintainer = "Carl H Lubba",
    maintainer_email = 'c.lubba15@imperial.ac.uk',
    packages = ['PyPNS'],
    url='https://github.com/caaarl/PyPNS',
    download_url = 'https://github.com/caaarl/PyPNS/archive/0.1.0.tar.gz',
    license='LICENSE',
    description='A peripheral nerve simulator built on NEURON',
    #long_description=,
    install_requires = ['numpy', 'scipy', 'matplotlib'], # , 'neuron' excluded because we need a higher version.
    provides = ['PyPNS'],
    classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
            "Operating System :: OS Independent",
            ],
    )
