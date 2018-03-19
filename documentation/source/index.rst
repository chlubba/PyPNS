.. PyPNS documentation master file, created by
   sphinx-quickstart on Sat May 13 17:07:34 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

The PyPNS API
============

.. module:: PyPNS

To use ``PyPNS``, you need to::
		
	import PyPNS

.. * :class:`axonClass.Myelinated` and :class:`axonClass.Unmyelinated`

This will import the classes

* :class:`bundleClass.Bundle`
* :class:`extracellularMechanismClass.homogeneous` and :class:`extracellularMechanismClass.precomputedFEM`
* :class:`recordingMechanismClass.RecordingMechanism`
* :class:`stimulusClass.StimIntra`, :class:`stimulusClass.StimField`, :class:`upstreamSpikingClass.UpstreamSpiking`

as well as the modules

* :mod:`createGeometry`
* :mod:`nameSetters`
* :mod:`signalGeneration`
* :mod:`spikeTrainGeneration`

The ``Bundle`` Class
____________________

The central object in ``PyPNS`` is :class:`bundleClass.Bundle`. 

.. autoclass:: bundleClass.Bundle()

	.. automethod:: bundleClass.Bundle.__init__

Simulating a ``Bundle``
~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: bundleClass.Bundle.simulate

Configuring the Bundle for Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: bundleClass.Bundle.add_excitation_mechanism

.. automethod:: bundleClass.Bundle.add_recording_mechanism

Recovering the Results of a Simulation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: bundleClass.Bundle.get_CAP_from_file

.. automethod:: bundleClass.Bundle.get_SFAPs_from_file

.. automethod:: bundleClass.Bundle.get_voltage_from_file_one_axon

Rerunning the CAP Calculation from Precomputed Membrane Current
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: bundleClass.Bundle.compute_CAPs_from_imem_files

.. automethod:: bundleClass.Bundle.clear_all_CAP_files

.. automethod:: bundleClass.Bundle.clear_all_recording_mechanisms

The ``ExtracellularMechanism`` Class
____________________________________

An ``ExtracellularMechanism`` is needed in ``PyPNS`` to compute the extracellular field from currents. Two variants exist:


.. autoclass:: extracellularMechanismClass.homogeneous
	:members:

.. autoclass:: extracellularMechanismClass.precomputedFEM
	:members:

The ``RecordingMechanism`` Class
________________________________

A ``RecordingMechansim`` calculates the signal that would be picked up by an electrode. One electrode is approximated as a set of points.

.. autoclass:: recordingMechanismClass.RecordingMechanism
	:members:

Stimulating the Axons
_____________________

To make the axons spike, ``ExcitationMechanisms`` can be added to the ``Bundle``. The user can choose from the following.

The ``StimulationMechanism`` Class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: stimulusClass.StimIntra
	:members:

.. autoclass:: stimulusClass.StimField
	:members:

The ``UpstreamSpiking`` Class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: upstreamSpikingClass.UpstreamSpiking
	:members:

.. ###### FROM HERE THE MODULES ARE DOCUMENTED. ######

The ``createGeometry`` Module
_____________________________

.. automodule:: createGeometry
	:members:

The ``nameSetters`` Module
__________________________

.. automodule:: nameSetters
	:members:

The ``signalGeneration`` Module
_______________________________

.. automodule:: signalGeneration
	:members:

The ``spikeTrainGeneration`` Module
___________________________________

.. automodule:: spikeTrainGeneration
	:members:

Contents:

.. toctree::
   :maxdepth: 2


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

