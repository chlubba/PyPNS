import PyPN
import matplotlib.pyplot as plt
from PyPN.takeTime import *
import numpy as np
import os

# ------------------------------------------------------------------------------
# ------------------------------- SCRIPT CONTROL -------------------------------
# ------------------------------------------------------------------------------

calculationFlag = False # run simulation or load latest bundle with this parameters (not all taken into account for identification)
electricalStimulusOn = True


# ------------------------------------------------------------------------------
# --------------------------------- DEFINITION ---------------------------------
# ------------------------------------------------------------------------------

# ----------------------------- simulation params ---------------------------

tStop=100
timeRes=0.0025

# ----------------------------- bundle params -------------------------------

# set length of bundle and number of axons
lengthOfBundle = 20000 # 20000 # 400000
numberOfAxons = 500

# set the diameter distribution or fixed value
# see http://docs.scipy.org/doc/numpy/reference/routines.random.html
# 5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.
myelinatedDiam =  {'distName' : 'normal', 'params' : (1.0, 0.7)} # (2.0, 0.7)
unmyelinatedDiam = {'distName' : 'normal', 'params' : (0.7, 0.3)}

# axon definitions
myelinatedParameters = {'fiberD': myelinatedDiam}
unmyelinatedParameters = {'fiberD': unmyelinatedDiam}

# bundle guide
segmentLengthAxon = 10
bundleGuide = PyPN.createGeometry.get_bundle_guide_corner(lengthOfBundle, segmentLengthAxon)

# set all properties of the bundle
bundleParameters = {    'radius': 300, #150, #um Radius of the bundle (typically 0.5-1.5mm)
                        'length': lengthOfBundle, # um Axon length
                        'randomDirectionComponent' : 0,
                        # 'bundleGuide': bundleGuide,

                        'numberOfAxons': numberOfAxons, # Number of axons in the bundle
                        'pMyel': .1, # Percentage of myelinated fiber type A
                        'pUnmyel': .9, #Percentage of unmyelinated fiber type C
                        'paramsMyel': myelinatedParameters, #parameters for fiber type A
                        'paramsUnmyel': unmyelinatedParameters, #parameters for fiber type C

                        'tStop' : tStop,
                        'timeRes' : timeRes,

                        # 'saveI':True,
                        'saveV': False,

                        'numberOfSavedSegments' : 50, # number of segments of which the membrane potential is saved to disk
                        'downsamplingFactor': 100
}

# ----------------------------- stimulation params ---------------------------

# parameters of signals for stimulation
rectangularSignalParams = {'amplitude': 10.,  # Pulse amplitude (mA)
                           'frequency': 20.,  # Frequency of the pulse (kHz)
                           'dutyCycle': 0.5,  # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                           'stimDur': 0.05,  # Stimulus duration (ms)
                           'waveform': 'MONOPHASIC',  # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                           'delay': 0.,  # ms
                           # 'invert': True,
                           'timeRes': timeRes,
                           }


intraParameters = {'stimulusSignal': PyPN.signalGeneration.rectangular(**rectangularSignalParams)}

# ----------------------------- recording params -------------------------------

recordingParameters = { 'radius': 200,
                        'numberOfElectrodes': 3,
                        'positionMax': .5,
                        'numberOfPoles': 1,
                        'poleDistance': 200,
                        # 'sigma': 2
}

recordingParametersBip = {'radius': 200,
                          'numberOfElectrodes': 3,
                          'positionMax': .5,
                          'numberOfPoles': 2,
                          'poleDistance': 1000,
                          }

# ------------------------------------------------------------------------------
# ---------------------------------- CALCULATION -------------------------------
# ------------------------------------------------------------------------------


if calculationFlag:

    # create the bundle with all properties of axons and recording setup
    bundle = PyPN.Bundle(**bundleParameters)

    # spiking through a single electrical stimulation
    if electricalStimulusOn:
        bundle.add_excitation_mechanism(PyPN.StimIntra(**intraParameters))

    bundle.add_recording_mechanism(PyPN.FEMRecCuff2D(**recordingParameters))
    bundle.add_recording_mechanism(PyPN.RecCuff2D(**recordingParameters))
    bundle.add_recording_mechanism(PyPN.FEMRecCuff2D(**recordingParametersBip))
    bundle.add_recording_mechanism(PyPN.RecCuff2D(**recordingParametersBip))

    # PyPN.plot.geometry_definition(bundle)
    # plt.show()

    # run the simulation
    bundle.simulate()

    # save the bundle to disk
    PyPN.save_bundle(bundle)
else:

    # try to open a bundle with the parameters set above
    # bundle = PyPN.open_recent_bundle(bundleParameters)
    # bundle = PyPN.open_bundle_from_location('/media/carl/4ECC-1C44/PyPN/dt=0.0025 tStop=100 pMyel=0.1 pUnmyel=0.9 L=20000 nAxons=500/bundle00000')
    bundle = PyPN.open_bundle_from_location('/media/carl/4ECC-1C44/PyPN/dt=0.0025 tStop=100 pMyel=0.1 pUnmyel=0.9 L=20000 nAxons=500/bundle00001')

# ------------------------------------------------------------------------------
# ---------------------------------- PLOTTING ----------------------------------
# ------------------------------------------------------------------------------

print '\nStarting to plot'
# # # # PyPN.plot.geometry_definition(bundle)
# PyPN.plot.CAP1D(bundle)
# PyPN.plot.CAP1D(bundle, recMechIndex=1)
# PyPN.plot.CAP1D(bundle, recMechIndex=2)
# PyPN.plot.CAP1D(bundle, recMechIndex=3)

plt.figure()
time, CAP = bundle.get_CAP_from_file(0)
plt.plot(time, CAP[-1,:], label='FEM')
time, CAP = bundle.get_CAP_from_file(1)
plt.plot(time, CAP[-1,:]/2*0.3, label='homogeneous')
plt.xlabel('time [ms]')
plt.ylabel('voltage [mV]')
plt.legend()
plt.tight_layout()

plt.figure()
time, CAP = bundle.get_CAP_from_file(2)
plt.plot(time, CAP[-1,:], label='FEM')
time, CAP = bundle.get_CAP_from_file(3)
plt.plot(time, CAP[-1,:]/2*0.3, label='homogeneous')
plt.xlabel('time [ms]')
plt.ylabel('voltage [mV]')
plt.legend()
plt.tight_layout()

plt.figure()
time, CAP = bundle.get_CAP_from_file(0)
plt.plot(time, CAP[-1,:], label='FEM monopolar')
time, CAP = bundle.get_CAP_from_file(2)
plt.plot(time, CAP[-1,:]/2*0.3, label='FEN bipolar')
plt.xlabel('time [ms]')
plt.ylabel('voltage [mV]')
plt.legend()
plt.tight_layout()

# t, CAP = bundle.get_CAP_from_file()
# np.save(os.path.join('/media/carl/4ECC-1C44/PyPN/FEM_CAPs/forPoster', 'homogeneous2.npy'), [t, CAP])


# PyPN.plot.CAP1D(bundle, recMechIndex=1)
# PyPN.plot.voltage(bundle)
# PyPN.plot.diameterHistogram(bundle)
plt.show()

bundle = None