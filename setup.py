#!/usr/bin/env python
'''PyPN setup.py file'''

import os
import shutil
try:
    from setuptools import setup
except ImportError as ie:
    raise ie, 'please install setuptools'


#try and locate the nrnivmodl script of NEURON in PATH so that the
#NEURON extension files AXNODE.mod, xtra.mod and vecevent.mod can be compiled
from distutils.spawn import find_executable, spawn
if find_executable('nrnivmodl') is not None:
    os.chdir('PyPN')
    for path in ['x86_64', 'i686', 'powerpc']:
        if os.path.isdir(path):
            shutil.rmtree(path)
    spawn([find_executable('nrnivmodl')])
    os.chdir('..')
else:
    print("nrnivmodl script not found in PATH, thus NEURON .mod files could     " +
          "not be compiled")

setup(
    name = "PyPN",
    version = "0.1.0",
    maintainer = "Carl Lubba",
    maintainer_email = 'c.lubba15@imperial.ac.uk',
    packages = ['PyPN'],
    #url='',
    # download_url = '',
    # license='LICENSE',
    description='A module for simulating bundles of axons built on NEURON and LFPy',
    #long_description=,
    install_requires = [
        'numpy', 'scipy', 'matplotlib', 'neuron', 'Cython','LFPy'
        ],
    provides = ['PyPN'],
    )
