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

# setup(
#     name = "PyPN    ",
#     version = "1.1.1",
#     maintainer = "Espen Hagen",
#     maintainer_email = 'e.hagen@fz-juelich.de',
#     packages = ['LFPy'],
#     package_data = {'LFPy' : ['stick.hoc', 'sinsyn.mod',
#                               os.path.join('i686', '*'),
#                               os.path.join('i686', '.libs', '*'),
#                               os.path.join('x86_64', '*'),
#                               os.path.join('x86_64', '.libs', '*'),
#                               os.path.join('powerpc', '*'),
#                               os.path.join('powerpc', '.libs', '*'),
#                               ]},
#     cmdclass = cmdclass,
#     ext_modules = ext_modules,
#     url='http://LFPy.github.io',
#     download_url = 'https://github.com/LFPy/LFPy/tarball/v1.1.1',
#     license='LICENSE',
#     description='A module for modeling Local Field Potentials built on NEURON',
#     long_description=long_description,
#     classifiers=[
#         'License :: OSI Approved :: GNU General Public License (GPL)',
#         'Programming Language :: Python',
#         'Programming Language :: Python :: 2.6',
#         'Programming Language :: Python :: 2.7',
#         'Programming Language :: Python :: 3.4',
#         'Programming Language :: Cython',
#         'Operating System :: OS Independent',
#         'Topic :: Scientific/Engineering',
#         'Topic :: Scientific/Engineering :: Physics',
#         'Topic :: Utilities',
#         'Intended Audience :: Developers',
#         'Intended Audience :: Science/Research',
#         'Development Status :: 5 - Production/Stable',
#         ],
#     install_requires = [
#         'numpy', 'scipy', 'matplotlib', 'neuron', 'Cython'
#         ],
#     provides = ['LFPy'],
#     )