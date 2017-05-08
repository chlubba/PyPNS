PNPy
====

PNPy is a Python module for the simulation of periperal nerves. Axon models are simulated in the compartmental simulator NEURON (http://www.neuron.yale.edu/neuron) over its Python interface (http://www.frontiersin.org/neuroinformatics/10.3389/neuro.11.001.2009/abstract). Extracellular potentials from membrane currents or stimulation electrodes are calculated in a resistive, electro-quasistatic approximation of the Maxwell equations from either a homogeneous analytical solution, precomputed and imported finite element model outputs or analytical functions fit to FEM results. 

PNPy was developed in the Department of Bioengineering, Centre of Neurotechnology at Imperial College London.

This scientific software is released under the GNU Public License GPLv3.


Requirements
============

To install PNPy you will need:

- Python modules numpy, scipy and matplotlib
- NEURON (from http://www.neuron.yale.edu) compiled as a Python module. See below links for help.

http://www.tc.umn.edu/~haszx010/files/vpl_dbs_docs/Installation.html
http://www.davison.webfactional.com/notes/installation-neuron-python/


Installation
============


From source:
::

    tar -xzf PNPy-x.x.tar.gz
    cd PNPy-x.x
    (sudo) python setup.py install (--user)


Usage
============

To run PNPy properly, several NEURON extensions need to be compiled. Those are located in the 'mods'-directory of PNPy. Download this directory and run 
::

	nrnivmodl
from the console. A new folder will be generated containing the compiled files. They need to be present in the working directory of your project.


Documentation
=============

To generate the html documentation, issue from the PNPy source code directory:
::
    
    sphinx-build -b html /path/to/PNPy/documentation/sources path/to/dest

The main html file is now in path/to/dest/index.html

