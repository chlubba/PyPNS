import PyPN
import matplotlib.pyplot as plt
from PyPN.takeTime import *
import numpy as np
import os
from pprint import pprint
import sys
import cPickle as pickle

# ------------------------------------------------------------------------------
# ------------------------------- SCRIPT CONTROL -------------------------------
# ------------------------------------------------------------------------------

calculationFlag = True # run simulation or load latest bundle with this parameters (not all taken into account for identification)
electricalStimulusOn = True


# ------------------------------------------------------------------------------
# --------------------------------- DEFINITION ---------------------------------
# ------------------------------------------------------------------------------

# ----------------------------- simulation params ---------------------------

tStop=120
timeRes=0.0025

# ----------------------------- bundle params -------------------------------

# set length of bundle and number of axons
lengthOfBundle = 7000 # 20000 # 400000
numberOfAxons = 1

# bundle guide
segmentLengthAxon = 30
bundleGuide = PyPN.createGeometry.get_bundle_guide_straight(lengthOfBundle, segmentLengthAxon)

# ----------------------------- stimulation params ---------------------------

# parameters of signals for stimulation
rectangularSignalParams = {'amplitude': 50., #50,  # Pulse amplitude (mA)
                           'frequency': 20.,  # Frequency of the pulse (kHz)
                           'dutyCycle': 0.5,  # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                           'stimDur': 0.05,  # Stimulus duration (ms)
                           'waveform': 'MONOPHASIC',  # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                           'delay': 0.,  # ms
                           # 'invert': True,
                           # 'timeRes': timeRes,
                           }


intraParameters = {'stimulusSignal': PyPN.signalGeneration.rectangular(**rectangularSignalParams)}

# ----------------------------- recording params -------------------------------

recordingParametersNew = {'bundleGuide': bundleGuide,
                          'radius': 250,
                          'positionAlongBundle': 4000,
                          'numberOfPoles': 1,
                          'poleDistance': 1000,
                        }


# ------------------------------------------------------------------------------
# ---------------------------------- CALCULATION -------------------------------
# ------------------------------------------------------------------------------

diameters = np.flipud(np.arange(.2, 4., .5))
temperatures = np.arange(5, 46, 5)
Ras = np.arange(50, 300, 50)

RDCs = [0.2, 0.6, 1.0] # [0, 0.2, 0.4, 0.6, 0.8, 1.] # np.arange(0, 1., 0.15)

if calculationFlag:

    # (f, axarr) = plt.subplots(2, 3, sharey=True)

    onsetInd = np.ones((2, 3)) * 100 / 0.0025
    lengthInInd = 0

    maxAmp = 0
    minAmp = 0

    for RDCInd, RDC in enumerate(RDCs):

        plotRow = int(RDCInd/3)
        plotColumn = np.mod(RDCInd, 3)
        # axis = axarr[plotRow][plotColumn]

        # for temperatureInd, temperature in enumerate(temperatures):
        numRuns = 1
        for runInd in range(numRuns):

            vAPs = []

            LFPMech = PyPN.Extracellular.homogeneous(sigma=1)
            LFPMech2 = PyPN.Extracellular.precomputedFEM(bundleGuide)

            electrodePos = PyPN.createGeometry.circular_electrode(**recordingParametersNew)

            modularRecMech = PyPN.RecordingMechanism(electrodePos, LFPMech)

            # set the diameter distribution or fixed value
            # see http://docs.scipy.org/doc/numpy/reference/routines.random.html
            # 5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.
            myelinatedDiam = 1.  # {'distName' : 'normal', 'params' : (1.0, 0.7)} # (2.0, 0.7)
            unmyelinatedDiam = 1.  # {'distName' : 'normal', 'params' : (0.7, 0.3)}

            # axon definitions
            myelinatedParameters = {'fiberD': myelinatedDiam}
            unmyelinatedParameters = {'fiberD': unmyelinatedDiam}

            # set all properties of the bundle
            bundleParameters = {'radius': 300,  # 150, #um Radius of the bundle (typically 0.5-1.5mm)
                                'length': lengthOfBundle,  # um Axon length
                                'randomDirectionComponent': RDC,
                                # 'bundleGuide': bundleGuide,
                                'segmentLengthAxon': 30,

                                'numberOfAxons': 10,  # Number of axons in the bundle
                                'pMyel': 0.,  # Percentage of myelinated fiber type A
                                'pUnmyel': 1.,  # Percentage of unmyelinated fiber type C
                                'paramsMyel': myelinatedParameters,  # parameters for fiber type A
                                'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C

                                'tStop': tStop,
                                'timeRes': 'variable', #0.0025, #

                                # 'saveI':True,
                                # 'saveV': False,
                                # 'saveLocation': '/media/carl/4ECC-1C44/PyPN/',

                                'numberOfSavedSegments': 50,
                                # number of segments of which the membrane potential is saved to disk
                                }

            # create the bundle with all properties of axons and recording setup
            bundle = PyPN.Bundle(**bundleParameters)

            ax = PyPN.plot.geometry_definition(bundle)
            ax.set_xlim((1000, 1300))
            ax.set_ylim((-200, 200))
            ax.set_zlim((-200, 200))
            ax.view_init(90, 90)
            # plt.show()


plt.show()

bundle = None