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

tStop=50
timeRes=0.0025

# ----------------------------- bundle params -------------------------------

# set length of bundle and number of axons
lengthOfBundle = 5000 # 20000 # 400000
numberOfAxons = 1

# bundle guide
segmentLengthAxon = 30

# ----------------------------- stimulation params ---------------------------

# parameters of signals for stimulation
rectangularSignalParams = {'amplitude': 5., #50,  # Pulse amplitude (mA)
                           'frequency': 20.,  # Frequency of the pulse (kHz)
                           'dutyCycle': 0.5,  # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                           'stimDur': 0.05,  # Stimulus duration (ms)
                           'waveform': 'MONOPHASIC',  # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                           'delay': 0.,  # ms
                           # 'invert': True,
                           # 'timeRes': timeRes,
                           }


intraParameters = {'stimulusSignal': PyPN.signalGeneration.rectangular(**rectangularSignalParams)}



# ------------------------------------------------------------------------------
# ---------------------------------- CALCULATION -------------------------------
# ------------------------------------------------------------------------------

diametersUnmyel = (0.120, 0.17,  0.21,  0.26,  0.32,  0.37,  0.41,  0.47,  0.51,  0.56,  0.62,  0.67,  0.72,  0.77,  0.84,  0.92,  0.97,  1.02,  1.07,  1.12,  1.17,  1.22,  1.27,  1.32, 1.36, 1.41, 1.48, 1.52)# [0.12, 1.52] # np.arange(0.2, 2, 0.3)
# diametersMyel = (1.01,    1.19,  1.22,  1.40,  1.41,  1.58,  1.61,  1.78,  1.81,  1.99,  2.01,  2.18,  2.22, 2.39,    2.41,  2.58,  2.61,  2.79,  2.81,  2.99,  3.01,  3.19,  3.61,  3.79,  3.81,  3.99,  4.02,  4.20) # [1.1, 4.2] # [4.2] # [2.3, 2.6, 2.9] # np.arange(0.2, 4, 0.3)]

diametersMyel = np.arange(0.4, 1.1, 0.02)

diametersBothTypes = [diametersUnmyel, diametersMyel]

tStartPlots = [0.2, 0.05]

diameters = np.arange(0.2, 4, 0.7) # [1, 2, 4] # [4, 3, 2, 1, 0.5, 0.2] # np.flipud(np.arange(.2, 4., .3))
temperatures = np.arange(5, 46, 5)
Ras = np.arange(50, 300, 50)

RDCs = [0, 0.2, 0.4, 0.6, 0.8, 1.] # np.arange(0, 1., 0.15)

simTimes = [40, 20]

saveDict = {'unmyelinatedDiameters' : diametersUnmyel,
            'unmyelinatedSFAPsHomo': [],
            'unmyelinatedSFAPsFEM': [],
            'unmyelinatedCV' : [],
            't': [],
            'myelinatedDiameters': diametersMyel,
            'myelinatedSFAPsHomo': [],
            'myelinatedSFAPsFEM': [],
            'myelinatedCV' : [],
            }

if calculationFlag:

    # variables to safe SFAPs in
    SFAPsFEM = []
    SFAPsHomo = []

    legends = ['Unmyelinated', 'Myelinated']
    bundleLengths = [5000, 30000]
    for i in [0,1]:

        vAPCollection = []

        diameters = diametersBothTypes[i]

        recMechLegends = ['homogeneous', 'FEM']
        recMechMarkers = ['o', 'v']

        SFAPsHomo = []
        SFAPsFEM = []

        for diameterInd, diameter in enumerate(diameters):

            # set the diameter distribution or fixed value
            # see http://docs.scipy.org/doc/numpy/reference/routines.random.html
            # 5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.
            myelinatedDiam = diameter  # {'distName' : 'normal', 'params' : (1.0, 0.7)} # (2.0, 0.7)
            unmyelinatedDiam = diameter  # {'distName' : 'normal', 'params' : (0.7, 0.3)}

            # axon definitions
            myelinatedParameters = {'fiberD': myelinatedDiam}
            unmyelinatedParameters = {'fiberD': unmyelinatedDiam}

            # set all properties of the bundle
            bundleParameters = {'radius': 300,  # 150, #um Radius of the bundle (typically 0.5-1.5mm)
                                'length': bundleLengths[i],  # um Axon length
                                # 'randomDirectionComponent': RDC,
                                # 'bundleGuide': bundleGuide,

                                'numberOfAxons': numberOfAxons,  # Number of axons in the bundle
                                'pMyel': i,  # Percentage of myelinated fiber type A
                                'pUnmyel': 1 - i,  # Percentage of unmyelinated fiber type C
                                'paramsMyel': myelinatedParameters,  # parameters for fiber type A
                                'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C

                                'tStop': simTimes[i],
                                'timeRes': 0.0025, #'variable', #

                                # 'saveI':True,
                                # 'saveV': False,
                                'saveLocation': '/media/carl/4ECC-1C44/PyPN/',

                                'numberOfSavedSegments': 50,
                                # number of segments of which the membrane potential is saved to disk
                                }

            # create the bundle with all properties of axons and recording setup
            bundle = PyPN.Bundle(**bundleParameters)

            LFPMech = []
            LFPMech.append(PyPN.Extracellular.homogeneous(sigma=1))
            LFPMech.append(PyPN.Extracellular.precomputedFEM(bundle.bundleCoords, fieldName='oil_sigma_0.00001_xP0'))


            # spiking through a single electrical stimulation
            if electricalStimulusOn:
                bundle.add_excitation_mechanism(PyPN.StimIntra(**intraParameters))


            if i == 1:

                recordingParametersNew = {'bundleGuide': bundle.bundleCoords,
                                          'radius': 200,
                                          'positionAlongBundle': np.floor(bundleLengths[i]*0.8 / bundle.axons[0].lengthOneCycle) *
                                                                 bundle.axons[0].lengthOneCycle,
                                          'numberOfPoles': 1,
                                          'poleDistance': 1000,
                                          }

            else:
                recordingParametersNew = {'bundleGuide': bundle.bundleCoords,
                                          'radius': 200,
                                          'positionAlongBundle': 3000,
                                          'numberOfPoles': 1,
                                          'poleDistance': 1000,
                                          }

            electrodePos = PyPN.createGeometry.circular_electrode(**recordingParametersNew)

            modularRecMechs = []
            for recMechIndex in [0,1]:
                modularRecMechs.append(PyPN.RecordingMechanism(electrodePos, LFPMech[recMechIndex]))
                bundle.add_recording_mechanism(modularRecMechs[-1])

            # run the simulation
            bundle.simulate()

            # PyPN.plot.voltage(bundle)

            # t, SFAPs = bundle.get_SFAPs_from_file()
            # plt.plot(t[t>tStartPlots[i]], SFAPs[t>tStartPlots[i]])
            # t, SFAPs = bundle.get_SFAPs_from_file(1)
            # plt.plot(t[t>tStartPlots[i]], SFAPs[t>tStartPlots[i]])
            # plt.show()

            t, v = bundle.get_voltage_from_file_one_axon(0)
            # plt.plot(t,v)
            # plt.show()


            currentNumberOfSegments = np.shape(v)[1]
            currentAxonLength = bundle.axons[0].L

            from PyPN.axonClass import *

            if not type(bundle.axons[0]) == Myelinated:
                iMaxFirstSegment = np.argmax(v[:, 0])
                iMaxLastSegment = np.argmax(v[:, -1])

                # plt.plot(v[:,0])
                # plt.show()

                tMinFirstSegment = t[iMaxFirstSegment]
                tMinLastSegment = t[iMaxLastSegment]

                distance = currentAxonLength

                vel = distance / 1000 / (tMinLastSegment - tMinFirstSegment)

            else:
                Nnodes = bundle.axons[0].axonnodes
                numberOfRecordingSites = np.shape(v)[1]
                nodePositions = range(0, (Nnodes - 1) * 11, 11)

                nodeDistance = bundle.axons[0].lengthOneCycle
                distance = Nnodes * nodeDistance

                iMaxFirstSegment = np.argmax(v[:, nodePositions[0]])
                iMaxLastSegment = np.argmax(v[:, nodePositions[-1]])

                # plt.plot(v[:,0])
                # plt.show()

                tMinFirstSegment = t[iMaxFirstSegment]
                tMinLastSegment = t[iMaxLastSegment]

                vel = distance / 1000 / (tMinLastSegment - tMinFirstSegment)

            t, SFAPs = bundle.get_SFAPs_from_file(0)
            SFAPsHomo.append(np.squeeze(SFAPs))
            _, SFAPs = bundle.get_SFAPs_from_file(1)
            SFAPsFEM.append(np.squeeze(SFAPs))

        # saveDict = {'unmyelinatedDiameters': diametersUnmyel,
        #             'unmyelinatedSFAPs': [],
        #             't': [],
        #             'myelinatedDiameters': diametersMyel,
        #             'myelinatedSFAPs': [],
        #             }
            if i == 0:
                saveDict['unmyelinatedCV'].append(vel)
            else:
                saveDict['myelinatedCV'].append(vel)

        if i == 0:
            saveDict['t'] = t
            saveDict['unmyelinatedSFAPsHomo'] = SFAPsHomo
            saveDict['unmyelinatedSFAPsFEM'] = SFAPsFEM
        else:
            saveDict['myelinatedSFAPsHomo'] = SFAPsHomo
            saveDict['myelinatedSFAPsFEM'] = SFAPsFEM


else:

    # try to open a bundle with the parameters set above
    # bundle = PyPN.open_recent_bundle(bundleParameters)
    # bundle = PyPN.open_bundle_from_location('/media/carl/4ECC-1C44/PyPN/dt=0.0025 tStop=100 pMyel=0.1 pUnmyel=0.9 L=20000 nAxons=500/bundle00000')
    # bundle = PyPN.open_bundle_from_location('/media/carl/4ECC-1C44/PyPN/dt=0.0025 tStop=100 pMyel=0.1 pUnmyel=0.9 L=20000 nAxons=500/bundle00001')
    bundle = PyPN.open_bundle_from_location(
        '/media/carl/4ECC-1C44/PyPN/dt=0.0025 tStop=60 pMyel=0.0 pUnmyel=1.0 L=5000 nAxons=13/bundle00000')

# ------------------------------------------------------------------------------
# ---------------------------------- PLOTTING ----------------------------------
# ------------------------------------------------------------------------------


pickle.dump(saveDict, open(os.path.join('/media/carl/4ECC-1C44/PyPN/SFAPs', 'SFAPsOil.dict'), "wb"))

print '\nStarting to plot'

# # pp = pprint.PrettyPrinter(indent=4)
# # pp.pprint(bundle)
# pprint (vars(bundle))
# # for axon in bundle.axons:
# #     # print axon.nbytes
# #     pprint (vars(axon))
# pprint(vars(bundle.axons[0]))
# pprint(vars(bundle.recordingMechanisms[0]))




# # first load the desired data from file
# for i in range(len(bundle.recordingMechanisms)):
#     time, CAP = bundle.get_CAP_from_file(i)
#     plt.plot(time, CAP, label='recMech'+str(i))
# plt.legend()

# plt.figure()
# plt.plot(diameters, vAPs)



# # # # PyPN.plot.geometry_definition(bundle)
# PyPN.plot.CAP1D(bundle)
# PyPN.plot.CAP1D(bundle, recMechIndex=1)
# PyPN.plot.CAP1D(bundle, recMechIndex=2)
# PyPN.plot.CAP1D(bundle, recMechIndex=3)

# plt.figure()
# time, CAP = bundle.get_CAP_from_file(0)
# plt.plot(time, CAP[-1,:], label='FEM')
# time, CAP = bundle.get_CAP_from_file(1)
# plt.plot(time, CAP[-1,:]/2*0.3, label='homogeneous')
# plt.xlabel('time [ms]')
# plt.ylabel('voltage [mV]')
# plt.legend()
# plt.tight_layout()
#
# plt.figure()
# time, CAP = bundle.get_CAP_from_file(2)
# plt.plot(time, CAP[-1,:], label='FEM')
# time, CAP = bundle.get_CAP_from_file(3)
# plt.plot(time, CAP[-1,:]/2*0.3, label='homogeneous')
# plt.xlabel('time [ms]')
# plt.ylabel('voltage [mV]')
# plt.legend()
# plt.tight_layout()
#
# plt.figure()
# time, CAP = bundle.get_CAP_from_file(0)
# plt.plot(time, CAP[-1,:], label='FEM monopolar')
# time, CAP = bundle.get_CAP_from_file(2)
# plt.plot(time, CAP[-1,:]/2*0.3, label='FEN bipolar')
# plt.xlabel('time [ms]')
# plt.ylabel('voltage [mV]')
# plt.legend()
# plt.tight_layout()

# t, CAP = bundle.get_CAP_from_file()
# np.save(os.path.join('/media/carl/4ECC-1C44/PyPN/FEM_CAPs/forPoster', 'homogeneous2.npy'), [t, CAP])


# PyPN.plot.CAP1D(bundle, recMechIndex=1)
# PyPN.plot.voltage(bundle)
# PyPN.plot.diameterHistogram(bundle)
plt.show()

bundle = None