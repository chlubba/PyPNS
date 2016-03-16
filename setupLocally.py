#!/usr/bin/env python

# this file compiles the necessary extensions of NEURON

import os
import shutil

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