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
bundleGuide = PyPN.createGeometry.get_bundle_guide_straight(lengthOfBundle, segmentLengthAxon)

# ----------------------------- recording params -------------------------------

recordingParametersNew = {'bundleGuide': bundleGuide,
                          'radius': 100,
                          'positionAlongBundle': 10000,
                          'numberOfPoles': 2,
                          'poleDistance': 1000,
                        }


# ------------------------------------------------------------------------------
# ---------------------------------- CALCULATION -------------------------------
# ------------------------------------------------------------------------------

diametersUnmyel = [3] # np.arange(0.2, 2, 0.3)
diametersMyel = [3] # [2.3, 2.6, 2.9] # np.arange(0.2, 4, 0.3)
diametersBothTypes = [diametersUnmyel, diametersMyel]

tStartPlots = [0.2, 0.05]

diameters = np.arange(0.2, 4, 0.7) # [1, 2, 4] # [4, 3, 2, 1, 0.5, 0.2] # np.flipud(np.arange(.2, 4., .3))
temperatures = np.arange(5, 46, 5)
Ras = np.arange(50, 300, 50)

RDCs = [0, 0.2, 0.4, 0.6, 0.8, 1.] # np.arange(0, 1., 0.15)

if calculationFlag:

    (f, axarr) = plt.subplots(1, 2, sharey=True)

    legends = ['Unmyelinated', 'Myelinated']
    bundleLengths = [5000, 50000]
    for i in [1]:

        vAPCollection = []

        diameters = diametersBothTypes[i]

        recMechLegends = ['homogeneous', 'FEM Saline', 'FEM Oil 2cm .00001', 'FEM Oil 5cm .00001'] # , 'FEM Oil only .0000001', 'FEM Oil contact']
        recMechMarkers = ['o', 'v']
        for recMechIndex in [2, 3]:

            vAPs = []
            vAPs2 = []

            LFPMech = []
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

                                    'tStop': tStop,
                                    'timeRes': 0.0025, #'variable', #

                                    # 'saveI':True,
                                    # 'saveV': False,
                                    'saveLocation': '/media/carl/4ECC-1C44/PyPN/',

                                    'numberOfSavedSegments': 50,
                                    # number of segments of which the membrane potential is saved to disk
                                    }

                # create the bundle with all properties of axons and recording setup
                bundle = PyPN.Bundle(**bundleParameters)

                LFPMech.append(PyPN.Extracellular.homogeneous(sigma=1))
                LFPMech.append(PyPN.Extracellular.precomputedFEM(bundle.bundleCoords))
                LFPMech.append(PyPN.Extracellular.precomputedFEM(bundle.bundleCoords, fieldName='oil_sigma_0.00001_xP0'))
                # LFPMech.append(PyPN.Extracellular.precomputedFEM(bundle.bundleCoords, fieldName='oil_sigma_0.0000001_xP0'))
                # LFPMech.append(PyPN.Extracellular.precomputedFEM(bundle.bundleCoords, fieldName='oil_sigma_0.0000001_contact_xP0'))
                LFPMech.append(PyPN.Extracellular.precomputedFEM(bundle.bundleCoords, fieldName='oil_sigma_0.00001_Oil5cm_xP0'))

                # ----------------------------- stimulation params ---------------------------

                amplitudes = (50., .5)

                # parameters of signals for stimulation
                rectangularSignalParams = {'amplitude': amplitudes[i], # 50.,  # 50,  # Pulse amplitude (mA)
                                           'frequency': 20.,  # Frequency of the pulse (kHz)
                                           'dutyCycle': 0.5,  # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                                           'stimDur': 0.05,  # Stimulus duration (ms)
                                           'waveform': 'MONOPHASIC',  # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                                           'delay': 0.,  # ms
                                           # 'invert': True,
                                           # 'timeRes': timeRes,
                                           }

                intraParameters = {'stimulusSignal': PyPN.signalGeneration.rectangular(**rectangularSignalParams)}

                # spiking through a single electrical stimulation
                if electricalStimulusOn:
                    bundle.add_excitation_mechanism(PyPN.StimIntra(**intraParameters))

                # bundle.add_recording_mechanism(PyPN.FEMRecCuff2D(**recordingParameters))
                # bundle.add_recording_mechanism(PyPN.RecCuff2D(**recordingParameters))
                # bundle.add_recording_mechanism(PyPN.FEMRecCuff2D(**recordingParametersBip))
                # bundle.add_recording_mechanism(PyPN.RecCuff2D(**recordingParametersBip))

                relPositions = [0] # np.arange(0, 0.5, 0.05)
                modularRecMech = [[] for jj in range(len(relPositions))]
                if i == 1:

                    for ii, relPos in enumerate(relPositions):

                        recordingParametersNew = {'bundleGuide': bundle.bundleCoords,
                                                  'radius': 200,
                                                  'positionAlongBundle': np.floor(bundleLengths[i]*0.8 / bundle.axons[0].lengthOneCycle) *
                                                                         bundle.axons[0].lengthOneCycle + bundle.axons[0].lengthOneCycle*relPos,
                                                  'numberOfPoles': 1,
                                                  'poleDistance': 2000,
                                                  }
                        electrodePos = PyPN.createGeometry.circular_electrode(**recordingParametersNew)

                        modularRecMech[ii] = PyPN.RecordingMechanism(electrodePos, LFPMech[recMechIndex])

                        bundle.add_recording_mechanism(modularRecMech[ii])

                else:
                    recordingParametersNew = {'bundleGuide': bundle.bundleCoords,
                                              'radius': 200,
                                              'positionAlongBundle': 3000,
                                              'numberOfPoles': 1,
                                              'poleDistance': 2000,
                                              }
                    electrodePos = PyPN.createGeometry.circular_electrode(**recordingParametersNew)

                    modularRecMech[0] = PyPN.RecordingMechanism(electrodePos, LFPMech[recMechIndex])

                    bundle.add_recording_mechanism(modularRecMech[0])

                # run the simulation
                bundle.simulate()

                # import matplotlib.cm as cm
                # import matplotlib.colors as colors
                #
                # jet = plt.get_cmap('jet')
                # cNorm = colors.Normalize(vmin=0, vmax=len(relPositions) - 1)
                # scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)
                #
                # for i in range(len(bundle.recordingMechanisms)):
                #
                #     colorVal = scalarMap.to_rgba(i)
                #
                #     t, SFAPs = bundle.get_SFAPs_from_file(i)
                #     plt.plot(t, SFAPs, label=str(relPositions[i]*100)+'%', color=colorVal)
                # # plt.legend(loc='best')
                # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
                # plt.xlabel('time [ms]')
                # plt.ylabel('$V_{ext}$ [mV]')
                # plt.grid()
                # plt.show()
                #
                #
                # tStartPlot = tStartPlots[i] # 0.2
                # tStopPlot = 40
                # selectionArray = np.logical_and(t>tStartPlot, t<tStopPlot)
                # tNoArt = t[selectionArray]
                # SFAPNoArtScaled = SFAPs[selectionArray]
                #
                # vAPs.append(np.max(SFAPNoArtScaled) - np.min(SFAPNoArtScaled))
                # print vAPs
                #
                # plt.plot(tNoArt, SFAPNoArtScaled)
                # # plt.title(recMechLegends[recMechIndex] + ' ' + legends[i])
                # # plt.show()

            t, SFAPs = bundle.get_SFAPs_from_file()

            axarr[i].plot(t,SFAPs, label=recMechLegends[recMechIndex])
            axarr[i].set_xlabel('time [ms]')

            if i == 0:
                axarr[i].set_xlim((1.8, 7.))
                axarr[i].set_ylabel('$V_{ext}$ [mV]')
                axarr[i].set_title('Unmyelinated')
            else:
                axarr[i].set_xlim((.4, 1.4))
                # axarr[i].legend()
                axarr[i].set_title('Myelinated')

        axarr[i].grid()
        axarr[i].set_ylim((-0.008, 0.004))
        plt.legend()
        # plt.figure()
        # plt.plot(diameters, np.divide(vAPCollection[0], vAPCollection[1]))


    # plt.xlabel('diameter [$\mu$m]')
    # plt.ylabel('$V_{ext}$ [mV]')
    # plt.legend(loc='best')
    plt.show()
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


# pickle.dump(saveDict, open(os.path.join('/media/carl/4ECC-1C44/PyPN/condVel', 'conductionVelocitiesUnmyelinatedRa.dict'), "wb"))

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