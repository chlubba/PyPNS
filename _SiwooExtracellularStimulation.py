import PyPN
import matplotlib.pyplot as plt
from PyPN.takeTime import *
import numpy as np

# import cPickle as pickle
# import os

calculationFlag = True # run simulation or load latest bundle with this parameters (not all taken into account for identification)

electricalStimulusOn = True

# set simulation params
tStop=10
timeRes=0.0001#0.0025

# set length of bundle and number of axons
lengthOfBundle = 1000 # 4000 # 400000
numberOfAxons = 1

# set the diameter distribution or fixed value
# see http://docs.scipy.org/doc/numpy/reference/routines.random.html
# 5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.
myelinatedDiam =  {'distName' : 'uniform', 'params' : (1.5, 4)} # .2 #
unmyelinatedDiam = 1. # {'distName' : 'uniform', 'params' : (0.1, 2)} # .2 #

# definition of the stimulation type of the axon
cuffParameters = {      'amplitude': .08, # 0.005, # 0.016,#0.2,# .0001,#1.5, #0.2, # 0.004, # 10., #  # Pulse amplitude (mA)
                        'frequency': 5., # 20., # Frequency of the pulse (kHz)
                        'dutyCycle': .5, # 0.05, # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                        'stimDur' : 0.2, # 0.05, # Stimulus duration (ms)
                        'waveform': 'BIPHASIC', # 'MONOPHASIC', ## Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                        'radius' : 50, #um
                        'timeRes' : timeRes,
                        'delay': 5, # ms
                        'invert': True,
                        'rho':500 # S/cm
}

recordingParameters = { 'radius': 200,
                        'numberOfElectrodes': 2,
                        'positionMax': 1.,
                        'numberOfPoles': 2
}

# axon parameters
myelinatedParameters = {'fiberD': myelinatedDiam, # um Axon diameter (5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0)
}

# axon parameters
unmyelinatedParameters = {'fiberD': unmyelinatedDiam, # um Axon diameter
}

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


if calculationFlag:

    # create the bundle with all properties of axons and recording setup
    bundle = PyPN.Bundle(**bundleParameters)

    # spiking through a single electrical stimulation
    if electricalStimulusOn:
        bundle.add_excitation_mechanism(PyPN.StimCuff(**cuffParameters))


    bundle.add_recording_mechanism(PyPN.RecCuff2D(**recordingParameters))

    # run the simulation
    bundle.simulate()

    # save the bundle to disk
    PyPN.save_bundle(bundle)
else:

    # try to open a bundle with the parameters set above
    bundle = PyPN.open_recent_bundle(bundleParameters)
    # bundle = PyPN.open_bundle_from_location('/media/carl/4ECC-1C44/PyPN/dt=0.0025 tStop=20 pMyel=1.0 pUnmyel=0.0 L=4000 nAxons=2/bundle00010')


print '\nStarting to plot'
PyPN.plot.geometry_definition(bundle)
PyPN.plot.CAP1D(bundle)
PyPN.plot.voltage(bundle)
PyPN.plot.diameterHistogram(bundle)
plt.show()

bundle = None