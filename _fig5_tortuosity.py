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

RDCs = [0, 0.2, 0.4, 0.6, 0.8, 1.] # np.arange(0, 1., 0.15)

if calculationFlag:

    (f, axarr) = plt.subplots(2, 3, sharey=True)

    onsetInd = np.ones((2, 3)) * 100 / 0.0025
    lengthInInd = 0

    maxAmp = 0
    minAmp = 0

    for RDCInd, RDC in enumerate(RDCs):

        plotRow = int(RDCInd/3)
        plotColumn = np.mod(RDCInd, 3)
        axis = axarr[plotRow][plotColumn]

        # for temperatureInd, temperature in enumerate(temperatures):
        numRuns = 4
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

                                'numberOfAxons': numberOfAxons,  # Number of axons in the bundle
                                'pMyel': 0.,  # Percentage of myelinated fiber type A
                                'pUnmyel': 1.,  # Percentage of unmyelinated fiber type C
                                'paramsMyel': myelinatedParameters,  # parameters for fiber type A
                                'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C

                                'tStop': tStop,
                                'timeRes': 'variable', #0.0025, #

                                # 'saveI':True,
                                # 'saveV': False,
                                'saveLocation': '/media/carl/4ECC-1C44/PyPN/',

                                'numberOfSavedSegments': 50,
                                # number of segments of which the membrane potential is saved to disk
                                }

            # create the bundle with all properties of axons and recording setup
            bundle = PyPN.Bundle(**bundleParameters)

            # spiking through a single electrical stimulation
            if electricalStimulusOn:
                bundle.add_excitation_mechanism(PyPN.StimIntra(**intraParameters))

            # bundle.add_recording_mechanism(PyPN.FEMRecCuff2D(**recordingParameters))
            # bundle.add_recording_mechanism(PyPN.RecCuff2D(**recordingParameters))
            # bundle.add_recording_mechanism(PyPN.FEMRecCuff2D(**recordingParametersBip))
            # bundle.add_recording_mechanism(PyPN.RecCuff2D(**recordingParametersBip))

            # bundle.add_recording_mechanism(modularRecMech1)
            bundle.add_recording_mechanism(modularRecMech)
            # bundle.add_recording_mechanism(modularRecMech3)

            # PyPN.plot.geometry_definition(bundle)
            # plt.show()

            # run the simulation
            bundle.simulate()

            t, SFAPs = bundle.get_SFAPs_from_file()

            # import matplotlib.cm as cm
            # import matplotlib.colors as colors
            #
            # jet = plt.get_cmap('jet')
            # cNorm = colors.Normalize(vmin=0, vmax=numRuns - 1)
            # scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)
            # colorVal = scalarMap.to_rgba(runInd)

            from scipy.interpolate import interp1d
            f = interp1d(t, np.squeeze(SFAPs))
            tReg = np.arange(0, max(t), 0.0025)
            SFAPReg = f(tReg)

            SFAPActive = np.where(np.abs(SFAPReg) > 0.01)[0]
            SFAPActive = SFAPActive[SFAPActive>1/0.0025]

            # find first
            onsetInd[plotRow, plotColumn] = np.min((np.min(SFAPActive), onsetInd[plotRow, plotColumn]))
            lengthInInd = np.max((np.max(SFAPActive) - onsetInd[plotRow, plotColumn], lengthInInd))

            maxAmp = np.max((np.max(SFAPReg[SFAPActive]), maxAmp))
            minAmp = np.min((np.min(SFAPReg[SFAPActive]), minAmp))

            # tStartPlot = 1
            # tStopPlot = 40
            # selectionArray = np.logical_and(tReg>tStartPlot, tReg<tStopPlot)
            # tNoArt = tReg[selectionArray]
            # SFAPNoArtScaled = SFAPReg[selectionArray]


            # axis.plot(tNoArt, SFAPNoArtScaled, color=colorVal)
            axis.plot(tReg, SFAPReg) # , color=colorVal)
            # axis.set_ylim((-2, 1))
            # plt.show()

            if plotColumn == 0:
                axis.set_ylabel('$V_{ext}$ [mV]')
            #
            # if plotRow == 0:
            #     axis.tick_params(
            #         axis='x',  # changes apply to the x-axis
            #         which='both',  # both major and minor ticks are affected
            #         bottom='off',  # ticks along the bottom edge are off
            #         top='off',  # ticks along the top edge are off
            #         labelbottom='off')  # labels along the bottom edge are off
            # else:
            #     axis.set_xlabel('time [ms]')
            if plotRow == 1:
                axis.set_xlabel('time [ms]')

            # PyPN.plot.geometry_definition(bundle)

            # # get SFAP
            # time, CAP = bundle.get_CAP_from_file()
            # plt.plot(time, CAP, label='RDC = ' + str(RDC))
            # plt.show()
        axis.set_title(str(RDC))
        axis.grid()

    for rowI in [0,1]:
        for columnI in [0,1,2]:
            axarr[rowI][columnI].set_xlim((onsetInd[rowI][columnI] * 0.0025 - 1, (onsetInd[rowI][columnI] + lengthInInd) * 0.0025 + 1))
            axarr[rowI][columnI].set_ylim((minAmp*1.1, maxAmp*1.1))


    # plt.legend()

    # import matplotlib.cm as cm
    # import matplotlib.colors as colors
    # jet = plt.get_cmap('jet')
    # cNorm = colors.Normalize(vmin=0, vmax=len(temperatures)-1)
    # scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)
    # colorVal = scalarMap.to_rgba(RaInd)
    # plt.plot(diameters, vAPs, label=str(Ra) + ' Ohm cm', color=colorVal)
    #
    #         # # save the bundle to disk
    #         # PyPN.save_bundle(bundle)
    #
    # plt.plot(np.arange(0,3.7,0.2), np.sqrt(np.arange(0,3.7,0.2))*2, linestyle='--', color=np.array([1, 1, 1])*0.7, label='theory')
    # plt.xlabel('diameter [um]')
    # plt.ylabel('conduction velocity [m/s]')
    # plt.title('Unmyelinated Axon')
    # plt.legend(loc='best')
    # plt.grid()

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