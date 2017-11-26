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

Using pip or easy install

	pip install PNPy

	easy_install PNPy

From source:

    tar -xzf PNPy-x.x.tar.gz
    cd PNPy-x.x
    python setup.py install


Usage
============

To run PNPy properly, several NEURON extensions need to be compiled for the myelinated axon model. Those are located in the 'mods'-directory of PNPy. Download this directory and run 

	nrnivmodl

from the console. A new folder will be generated containing the compiled files. They need to be present in the working directory of your project.

If FEM results are to be used for either recording or stimulation, the 'Fields'-directory needs to be present within the working directory as well. Different field dictionaries need to be contained within subdirectories. Subdirectory name equals field name. One example field can be downloaded from GitHub.


Documentation
============

Do compile the Sphinx documentation, locate the source folder of PNPy and run

	sphinx-build -b html path/to/documentation/source output/path

