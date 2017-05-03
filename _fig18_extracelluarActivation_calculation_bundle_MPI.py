import PyPN
import matplotlib.pyplot as plt
from PyPN.takeTime import *
import numpy as np
import os
from pprint import pprint
import sys
import cPickle as pickle
import time

from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

# ------------------------------------------------------------------------------
# ------------------------------- SCRIPT CONTROL -------------------------------
# ------------------------------------------------------------------------------

calculationFlag = True # run simulation or load latest bundle with this parameters (not all taken into account for identification)


# ------------------------------------------------------------------------------
# --------------------------------- DEFINITION ---------------------------------
# ------------------------------------------------------------------------------

# ----------------------------- simulation params ---------------------------

tStop=20
timeRes=0.0025

# ----------------------------- bundle params -------------------------------

# set length of bundle and number of axons
lengthOfBundle = 25000 # 400000
numberOfAxons = 15

# bundle guide
segmentLengthAxon = 30
bundleGuide = PyPN.createGeometry.get_bundle_guide_straight(lengthOfBundle, segmentLengthAxon)

# ----------------------------- recording params -------------------------------

recordingParametersNew = {'bundleGuide': bundleGuide,
                          'radius': 100,
                          'positionAlongBundle': 10000,
                          'numberOfPoles': 1,
                          'poleDistance': 1000,
                        }


# ------------------------------------------------------------------------------
# ---------------------------------- CALCULATION -------------------------------
# ------------------------------------------------------------------------------

legends = ['Unmyelinated', 'Myelinated']

diameters = [3.] # np.arange(1.0, 4, 0.8)
RDCs = [0, 0.05, 0.075, 0.3, 0.7]# [0.05] # [0.1, 0.2, 0.3] # [0.5, 1.]
amplitudes = np.logspace(2, 2.3, 10) # np.logspace(1.5, 5, 15) * 3 # np.logspace(1.5, 3, 12) * 3 # [150, 450, 750, 1050, 1450, 1750, 2050, 2350, 2650] # [150, 300, 450, 600, 750, 900, 1050, 1300, 1450] # np.logspace(2,4,10)
runs = 1

activationMatrix = np.zeros((len(diameters), len(RDCs), len(amplitudes)))

calculationCounter = 0 # one number unique to each parameter combination and run
for i in [1]:

    tStart = time.time()

    for run in range(runs):

        for diameterInd, diameter in enumerate(diameters):

            # axon definitions
            myelinatedParameters = {'fiberD': diameter}
            unmyelinatedParameters = {'fiberD': diameter}

            for RDCInd, RDC in enumerate(RDCs):

                for amplitudeInd, amplitude in enumerate(amplitudes):

                    if np.mod(calculationCounter, size) == rank:

                        print 'run:'
                        print run
                        print 'amplitude:'
                        print amplitude
                        print 'RDC:'
                        print RDC
                        print 'diameter:'
                        print diameter

                        # ----------------------------- stimulation params ---------------------------

                        # parameters of signals for stimulation
                        rectangularSignalParams = {'amplitude': amplitude,  # 100000.,# .0005, # #50,  # Pulse amplitude (nA)
                                                   'frequency': 1,  # Frequency of the pulse (kHz)
                                                   'dutyCycle': 0.5,
                                                   # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                                                   'stimDur': 1.,  # Stimulus duration (ms)
                                                   'waveform': 'BIPHASIC',
                                                   # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                                                   'delay': 0.,  # ms
                                                   # 'invert': True,
                                                   # 'timeRes': timeRes,
                                                   }

                        elecPosStim = PyPN.createGeometry.circular_electrode(bundleGuide, positionAlongBundle=12500, radius=235,
                                                                             numberOfPoles=2, poleDistance=1000)
                        extPotMechStim = PyPN.Extracellular.precomputedFEM(bundleGuide, 'oil190Inner50Endoneurium')

                        extraParameters = {'stimulusSignal': PyPN.signalGeneration.rectangular(**rectangularSignalParams),
                                           'electrodePositions': elecPosStim,
                                           'extPotMech': extPotMechStim}  # extPotMechStim

                        # -------------------------------- bundle definition ------------------------------------

                        # set all properties of the bundle
                        bundleParameters = {'radius': 190,  # 150, #um Radius of the bundle (typically 0.5-1.5mm)
                                            'length': lengthOfBundle,  # um Axon length
                                            'randomDirectionComponent': RDC,
                                            # 'bundleGuide': bundleGuide,

                                            'numberOfAxons': numberOfAxons,  # Number of axons in the bundle
                                            'pMyel': i,  # Percentage of myelinated fiber type A
                                            'pUnmyel': 1 - i,  # Percentage of unmyelinated fiber type C
                                            'paramsMyel': myelinatedParameters,  # parameters for fiber type A
                                            'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C
                                            # 'axonCoords': [0, 180],

                                            'tStop': tStop,
                                            'timeRes': 0.0025, #'variable', #

                                            # 'saveI':True,
                                            # 'saveV': False,
                                            # 'saveLocation': '/media/carl/SANDISK/PyPN/Results',

                                            'numberOfSavedSegments': 50,
                                            # number of segments of which the membrane potential is saved to disk
                                            }

                        # create the bundle with all properties of axons and recording setup
                        bundle = PyPN.Bundle(**bundleParameters)

                        # spiking through a single electrical stimulation
                        bundle.add_excitation_mechanism(PyPN.StimFieldQuasistatic(**extraParameters))
                        # bundle.add_excitation_mechanism(PyPN.StimCuff(stimulusSignal=PyPN.signalGeneration.rectangular(**rectangularSignalParams), radius=235))

                        # run the simulation
                        bundle.simulate()

                        for ii in range(numberOfAxons):
                            t, v = bundle.get_voltage_from_file_one_axon(ii)
                            maxPostStimPot =  np.max(v[t>1,:])

                            if maxPostStimPot > -40:
                                activationMatrix[diameterInd, RDCInd, amplitudeInd] += 1

                        # PyPN.plot.voltage(bundle)
                        # plt.show()

                        bundle.clear_all_voltage_files()
                        bundle = None

                    calculationCounter += 1



if not rank == 0:
    comm.send(activationMatrix, dest=0)
else:
    for i in range(1,size):
        activationMatrix += comm.recv(source=MPI.ANY_SOURCE)

    saveDict = {'activationMatrix': activationMatrix,
                'RDCs': RDCs,
                'amplitudes': amplitudes,
                'diameters': diameters
                }

    print 'Overall processing of took %5.2f' % (time.time() - tStart)

    pickle.dump(saveDict, open(os.path.join('activationExtracellular', 'myelinated15AxonsRDC0_005_0075_03_07amp100_200LogRadius190.dict'), "wb"))



