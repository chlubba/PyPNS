import PyPN
import matplotlib.pyplot as plt
from PyPN.takeTime import *
import numpy as np

# ------------------------------------------------------------------------------
# ------------------------------- SCRIPT CONTROL -------------------------------
# ------------------------------------------------------------------------------

calculationFlag = True # run simulation or load latest bundle with this parameters (not all taken into account for identification)
electricalStimulusOn = True


# ------------------------------------------------------------------------------
# --------------------------------- DEFINITION ---------------------------------
# ------------------------------------------------------------------------------

# ----------------------------- simulation params ---------------------------

tStop=10
timeRes=0.0025

# ----------------------------- bundle params -------------------------------

# set length of bundle and number of axons
lengthOfBundle = 1000 # 400000
numberOfAxons = 1

# set the diameter distribution or fixed value
# see http://docs.scipy.org/doc/numpy/reference/routines.random.html
# 5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.
myelinatedDiam =  {'distName' : 'uniform', 'params' : (1.5, 4)}
unmyelinatedDiam = {'distName' : 'uniform', 'params' : (0.1, 2)}

# axon definitions
myelinatedParameters = {'fiberD': myelinatedDiam}
unmyelinatedParameters = {'fiberD': unmyelinatedDiam}

# set all properties of the bundle
bundleParameters = {    'radius': 150, #150, #um Radius of the bundle (typically 0.5-1.5mm)
                        'length': lengthOfBundle, # um Axon length
                        'randomDirectionComponent' : 0,

                        'numberOfAxons': numberOfAxons, # Number of axons in the bundle
                        'pMyel': 0., # Percentage of myelinated fiber type A
                        'pUnmyel': 1., #Percentage of unmyelinated fiber type C
                        'paramsMyel': myelinatedParameters, #parameters for fiber type A
                        'paramsUnmyel': unmyelinatedParameters, #parameters for fiber type C

                        'tStop' : tStop,
                        'timeRes' : timeRes,

                        # 'saveI':True,
                        # 'saveV':False,

                        'numberOfSavedSegments' : 50, # number of segments of which the membrane potential is saved to disk
                        'downsamplingFactor': 100
}

# ----------------------------- stimulation params ---------------------------

# parameters of signals for stimulation
rectangularSignalParams = {'amplitude': 0.05,  # Pulse amplitude (mA)
                           'frequency': 20.,  # Frequency of the pulse (kHz)
                           'dutyCycle': 0.5,  # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                           'stimDur': 0.05,  # Stimulus duration (ms)
                           'waveform': 'MONOPHASIC',  # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                           'delay': 5,  # ms
                           'invert': True,
                            'timeRes': timeRes,
                           }

asymmetricSignalParams = {  'tDelay': 0, # ms
                            'tC':     1,  # ms
                            'aC':     0.01, #mA
                            'tExp':   1,  # ms
                            'cExp':   -5, # ms
                            'tD':     2,  # ms
                            'aD':     -0.001, #mA
                            'timeRes': timeRes
}

# definition of the stimulation type of the axon
stimCuffParameters = {'stimulusSignal': PyPN.signalGeneration.rectangular(**rectangularSignalParams),
                      'radius': 50,  # um
                      }
intraParameters = {'stimulusSignal': PyPN.signalGeneration.biphasic_decaying(**asymmetricSignalParams)}

# ----------------------------- recording params -------------------------------

recordingParameters = { 'radius': 200,
                        'numberOfElectrodes': 2,
                        'positionMax': 1.,
                        'numberOfPoles': 2
}

# ------------------------------------------------------------------------------
# ---------------------------------- CALCULATION -------------------------------
# ------------------------------------------------------------------------------


if calculationFlag:

    # create the bundle with all properties of axons and recording setup
    bundle = PyPN.Bundle(**bundleParameters)

    # spiking through a single electrical stimulation
    if electricalStimulusOn:
        bundle.add_excitation_mechanism(PyPN.StimCuff(**stimCuffParameters))
        # bundle.add_excitation_mechanism(PyPN.StimIntra(**intraParameters))


    bundle.add_recording_mechanism(PyPN.RecCuff2D(**recordingParameters))

    # run the simulation
    bundle.simulate()

    # save the bundle to disk
    PyPN.save_bundle(bundle)
else:

    # try to open a bundle with the parameters set above
    bundle = PyPN.open_recent_bundle(bundleParameters)
    # bundle = PyPN.open_bundle_from_location('/media/carl/4ECC-1C44/PyPN/dt=0.0025 tStop=20 pMyel=1.0 pUnmyel=0.0 L=4000 nAxons=2/bundle00010')

# ------------------------------------------------------------------------------
# ---------------------------------- PLOTTING ----------------------------------
# ------------------------------------------------------------------------------

print '\nStarting to plot'
PyPN.plot.geometry_definition(bundle)
PyPN.plot.CAP1D(bundle)
PyPN.plot.voltage(bundle)
PyPN.plot.diameterHistogram(bundle)
plt.show()

bundle = None