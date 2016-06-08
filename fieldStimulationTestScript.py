import PyPN
import matplotlib.pyplot as plt
from PyPN.takeTime import *
import numpy as np

# import cPickle as pickle
# import os

calculationFlag = True # run simulation or load latest bundle with this parameters (not all taken into account for identification)

upstreamSpikingOn = False
electricalStimulusOn = True

# set simulation params
tStop=10
timeRes=0.0025#0.0025

# set length of bundle and number of axons
lengthOfBundle = 1000 # 4000 # 400000
numberOfAxons = 2

# # create a guide the axons will follow
# segmentLengthAxon = 10
# # bundleGuide = PyPN.createGeometry.get_bundle_guide_straight(lengthOfBundle, segmentLengthAxon)
# bundleGuide = PyPN.createGeometry.get_bundle_guide_straight_2radii(lengthOfBundle, segmentLengthAxon, radii=(1000, 100))

# set the diameter distribution or fixed value
# see http://docs.scipy.org/doc/numpy/reference/routines.random.html
# 5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.
myelinatedDiam =  {'distName' : 'uniform', 'params' : (1.5, 4)} # .2 #
unmyelinatedDiam = {'distName' : 'uniform', 'params' : (0.1, 2)} # .2 #




intraParameters = {     'amplitude': 3., # 0.005, # 0.016,#0.2,# .0001,#1.5, #0.2, # 0.004, # 10., #  # Pulse amplitude (mA)
                        'frequency': 20., # Frequency of the pulse (kHz)
                        'dutyCycle': .5, # 0.05, # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                        'stimDur' : 0.05, # Stimulus duration (ms)
                        'waveform': 'MONOPHASIC', # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                        'timeRes' : timeRes,
                        'delay': 5, # ms
                        # 'invert': True
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
                        # 'bundleGuide' : bundleGuide,
                        # 'randomDirectionComponent' : 0,

                        'numberOfAxons': numberOfAxons, # Number of axons in the bundle
                        'pMyel': .1, # Percentage of myelinated fiber type A
                        'pUnmyel': .9, #Percentage of unmyelinated fiber type C
                        'paramsMyel': myelinatedParameters, #parameters for fiber type A
                        'paramsUnmyel': unmyelinatedParameters, #parameters for fiber type C

                        'tStop' : tStop,
                        'timeRes' : timeRes,

                        # 'saveI':True,
                        # 'saveV':False
}

# combine parameters for the bundle creation
Parameters = dict(bundleParameters, **recordingParameters)

if calculationFlag:

    # create the bundle with all properties of axons and recording setup
    bundle = PyPN.Bundle(**bundleParameters)

    # # spiking through a single electrical stimulation
    # if electricalStimulusOn:
    #     # stimulusInstance = PyPN.Stimulus(**stimulusParameters)
    #     # plt.plot(stimulusInstance.t, stimulusInstance.stimulusSignal)
    #     # plt.title('stimulus signal without delay')
    #     # plt.show()
    #     # bundle.add_excitation_mechanism(PyPN.StimCuff(**cuffParameters))
    #     # bundle.add_excitation_mechanism(PyPN.SimpleIClamp(**stimulusParameters))
    #
    #     # bundle.add_excitation_mechanism(PyPN.StimTripolarPoint(radius=1000, poleDistance=100, stimDur=1, amplitude=10.5, frequency=1, dutyCycle=0.5, waveform='BIPHASIC', timeRes=timeRes, delay=5))
    #     # bundle.add_excitation_mechanism(PyPN.StimCuff(radius=1000, stimDur=1, amplitude=10.5, frequency=1, dutyCycle=0.5, waveform='BIPHASIC', timeRes=timeRes, delay=5))
    #
    #     # bundle.add_excitation_mechanism(PyPN.StimIntra(**intraParameters))
    #
    #     bundle.add_excitation_mechanism(
    #         PyPN.StimField('/media/carl/4ECC-1C44/PyPN/voltageFieldDummy2', 30, 0.1, bundle))



    print bundle.approx_num_segs()

    bundle.add_recording_mechanism(PyPN.RecCuff2D(**recordingParameters))

    # PyPN.plot.geometry_definition(bundle)
    # plt.show()

    # run the simulation
    bundle.simulate()


    # save the bundle to disk
    PyPN.save_bundle(bundle)
else:

    # try to open a bundle with the parameters set above
    bundle = PyPN.open_recent_bundle(Parameters)
    # bundle = PyPN.open_bundle_from_location('/media/carl/4ECC-1C44/PyPN/dt=0.0025 tStop=20 pMyel=1.0 pUnmyel=0.0 L=4000 nAxons=2/bundle00010')



PyPN.plot.CAP1D(bundle)
PyPN.plot.voltage(bundle)


plt.show()

bundle = None