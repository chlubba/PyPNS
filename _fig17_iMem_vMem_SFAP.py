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

tStop=40
timeRes=0.0025

# ----------------------------- bundle params -------------------------------

# set length of bundle and number of axons
lengthOfBundle = 10000 # 20000 # 400000
numberOfAxons = 1

# bundle guide
segmentLengthAxon = 30
bundleGuide = PyPN.createGeometry.get_bundle_guide_straight(lengthOfBundle, segmentLengthAxon)

# ----------------------------- stimulation params ---------------------------

# parameters of signals for stimulation
rectangularSignalParams = {'amplitude': 10., # 100. #50,  # Pulse amplitude (mA)
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
                          'radius': 200,
                          'positionAlongBundle': 40000,
                          'numberOfPoles': 1,
                          'poleDistance': 1000,
                        }


# ------------------------------------------------------------------------------
# ---------------------------------- CALCULATION -------------------------------
# ------------------------------------------------------------------------------

diametersUnmyel = [1] # [0.2, 1, 3] # np.arange(0.2, 2, 0.3)
diametersMyel = [1] # [0.2, 1, 3] # [2.3, 2.6, 2.9] # np.arange(0.2, 4, 0.3)
diametersBothTypes = [diametersUnmyel, diametersMyel]

temperatures = np.arange(5, 46, 5)
Ras = [70] # np.arange(50, 300, 50)
Cms = [0.1, 0.3, 1.] # np.arange(0.15, 0.2, 0.01)
gkbars = np.arange(0.1, 0.2, 0.02) # [1.2] #

# saveDict = {'axonType': 'unmyelinated',
#             'diameters': diameters,
#             'temperatures': temperatures,
#             'velocityArray': []}

if calculationFlag:

    firstAP = []

    (f,axarr) = plt.subplots(1,3)

    typeLegendStrings = ['unmyelinated', 'myelinated']
    for typeInd in [0]:

        diameters = diametersBothTypes[typeInd]

        LFPMech = []
        recMechIndex = 0
        # for temperatureInd, temperature in enumerate(temperatures):
        for diameterInd, diameter in enumerate(diameters):


            # set the diameter distribution or fixed value
            # see http://docs.scipy.org/doc/numpy/reference/routines.random.html
            # 5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.
            myelinatedDiam = diameter  # {'distName' : 'normal', 'params' : (1.0, 0.7)} # (2.0, 0.7)
            unmyelinatedDiam = diameter  # {'distName' : 'normal', 'params' : (0.7, 0.3)}

            # axon definitions
            myelinatedParameters = {'fiberD': myelinatedDiam} # , 'gkbar_axnode': gkbar} # , 'temperature': temperature}
            unmyelinatedParameters = {'fiberD': unmyelinatedDiam} # , 'temperature': temperature}

            # set all properties of the bundle
            bundleParameters = {'radius': 300,  # 150, #um Radius of the bundle (typically 0.5-1.5mm)
                                'length': lengthOfBundle,  # um Axon length
                                # 'randomDirectionComponent': .9,
                                # 'bundleGuide': bundleGuide,

                                'numberOfAxons': numberOfAxons,  # Number of axons in the bundle
                                'pMyel': typeInd,  # Percentage of myelinated fiber type A
                                'pUnmyel': 1 - typeInd,  # Percentage of unmyelinated fiber type C
                                'paramsMyel': myelinatedParameters,  # parameters for fiber type A
                                'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C

                                'tStop': tStop,
                                'timeRes': 0.0025, #'variable', #

                                'saveI':True,
                                # 'saveV': False,
                                'saveLocation': '/media/carl/4ECC-1C44/PyPN/',

                                'numberOfSavedSegments': 1000,
                                # number of segments of which the membrane potential is saved to disk
                                # 'downsamplingFactor': 100
                                }

            # create the bundle with all properties of axons and recording setup
            bundle = PyPN.Bundle(**bundleParameters)

            LFPMech.append(PyPN.Extracellular.homogeneous(sigma=1))
            LFPMech.append(PyPN.Extracellular.precomputedFEM(bundle.bundleCoords))

            # spiking through a single electrical stimulation
            if electricalStimulusOn:
                bundle.add_excitation_mechanism(PyPN.StimIntra(**intraParameters))

            relPositions = np.arange(0, 0.5, 0.05)
            modularRecMech = [[] for jj in range(len(relPositions))]
            if typeInd == 1:

                for ii, relPos in enumerate(relPositions):
                    recordingParametersNew = {'bundleGuide': bundle.bundleCoords,
                                              'radius': 200,
                                              'positionAlongBundle': np.floor(8000. / bundle.axons[0].lengthOneCycle) *
                                                                     bundle.axons[0].lengthOneCycle + bundle.axons[
                                                                                                          0].lengthOneCycle * relPos,
                                              'numberOfPoles': 1,
                                              'poleDistance': 1000,
                                              }
                    electrodePos = PyPN.createGeometry.circular_electrode(**recordingParametersNew)

                    modularRecMech[ii] = PyPN.RecordingMechanism(electrodePos, LFPMech[recMechIndex])

                    bundle.add_recording_mechanism(modularRecMech[ii])

            else:
                recordingParametersNew = {'bundleGuide': bundle.bundleCoords,
                                          'radius': 200,
                                          'positionAlongBundle': 8000,
                                          'numberOfPoles': 1,
                                          'poleDistance': 1000,
                                          }
                electrodePos = PyPN.createGeometry.circular_electrode(**recordingParametersNew)

                modularRecMech[0] = PyPN.RecordingMechanism(electrodePos, LFPMech[recMechIndex])

                bundle.add_recording_mechanism(modularRecMech[0])

            # run the simulation
            bundle.simulate()

            t, imem = bundle.get_imem_from_file_axonwise(0)

            t,v = bundle.get_voltage_from_file_one_axon(0)
            # plt.plot(t,v)
            # plt.show()

            if typeInd == 1:
                nodeIndex = np.floor(bundle.axons[0].totnsegs / 11 * 0.8) * 11
                iSignal = np.sum(imem[nodeIndex:nodeIndex+3, :], axis=0)
            else:
                nodeIndex = 10
                iSignal = imem[nodeIndex, :]

            # from scipy.interpolate import interp1d
            # f = interp1d(t, iSignal)
            # tReg = np.arange(0,max(t),0.0025)
            # vReg = f(tReg)
            #
            # sp = np.fft.fft(iSignal)
            # freq = np.fft.fftfreq(tReg.shape[-1], d=0.0025/1000)
            # plt.semilogy(freq, np.abs(sp[1:]), label=typeLegendStrings[typeInd] + ', ' + str(diameter)+ ' $\mu$m') # len(freq) #  / max(np.abs(sp[1:]))

            # plt.plot(t, iOneSegment)
            if typeInd == 1:
                labelStringsSections = ['node', 'MYSA', 'FLUT']

                # first plot voltage
                axarr[0].plot(t, np.mean(v[:, nodeIndex:nodeIndex+3], axis=1), label='mean of all', linewidth=1.5)
                for ii in range(3):
                    axarr[0].plot(t, v[:, nodeIndex + ii], label=labelStringsSections[ii])
                axarr[0].set_ylabel('$V_m$ [mv]')
                axarr[0].set_title('Membrane Voltage')
                axarr[0].legend(loc='best')
                axarr[0].set_xlim((1.6, 2.8))
                axarr[0].set_ylim((-90, 40))

                # plot membrane current in second subplot
                axarr[1].plot(t, np.sum(imem[nodeIndex:nodeIndex + 3, :], axis=0), label='Sum of all',
                                    linewidth=1.5) # , color='k'
                for ii in range(3):
                    axarr[1].plot(t, imem[nodeIndex + ii, :], label=labelStringsSections[ii])

                axarr[1].set_title('Membrane Current')
                axarr[1].legend(loc='best')
                axarr[1].set_ylabel('$I_m$ [nA]')
                axarr[1].grid()
                axarr[1].set_xlim((1.6, 2.8))
                axarr[1].set_ylim((-0.7, 0.3))
                # axarr[typeInd].set_xlim((0.85, 1.7))

                # plot CAP in third one
                import matplotlib.cm as cm
                import matplotlib.colors as colors

                jet = plt.get_cmap('jet')
                cNorm = colors.Normalize(vmin=0, vmax=len(relPositions) - 1)
                scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

                for i in range(len(bundle.recordingMechanisms)):
                    colorVal = scalarMap.to_rgba(i)

                    t, SFAPs = bundle.get_SFAPs_from_file(i)
                    plt.plot(t, SFAPs*1000, label=str(relPositions[i] * 100) + '%', color=colorVal)
                # plt.legend(loc='best')
                axarr[2].legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
                axarr[2].set_ylabel('$V_{ext}$ [$\mu$V]')
                axarr[2].set_title('Extracellular Potential')
                axarr[2].set_xlim((1.9, 2.9))
                axarr[2].set_ylim((-0.035, 0.025))


            else:

                axarr[0].plot(t, v[:, nodeIndex])
                # axarr[0].plot(t, v)
                axarr[0].set_ylabel('$V_m$ [mv]')
                axarr[0].set_title('Membrane Voltage')
                # # axarr[0].legend(loc='best')
                axarr[0].set_xlim((0, 3))
                axarr[0].set_ylim((-70, -10))

                # plot membrane current in second subplot
                imem = imem.T
                axarr[1].plot(t, imem[:, nodeIndex])
                # axarr[1].plot(t, np.transpose(imem))

                axarr[1].set_title('Membrane Current')
                axarr[1].set_ylabel('$I_m$ [nA]')
                # axarr[1].grid()
                axarr[1].set_xlim((0, 3))
                axarr[1].set_ylim((-0.35, 0.15))

                t, SFAPs = bundle.get_SFAPs_from_file()
                plt.plot(t, SFAPs * 1000)
                # plt.legend(loc='best')
                # axarr[2].legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1)
                axarr[2].set_ylabel('$V_{ext}$ [$\mu$V]')
                axarr[2].set_title('Extracellular Potential')
                axarr[2].set_xlim((16, 25))
                axarr[2].set_ylim((-0.02, 0.005))

            for i in range(3):
                axarr[i].set_xlabel('time [ms]')
                axarr[i].grid()




            # plt.tight_layout(pad=1.4)




    # plt.xlim((0,20000))
    # plt.grid()
    # plt.xlabel('Frequency [Hz]')
    # plt.ylabel('Fourier coeff [nA]')
    # plt.title('FFT of the membrane current')
    # plt.xlim((0, 50000))
    # plt.legend(loc='best')

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


# pickle.dump(saveDict, open(os.path.join('/media/carl/4ECC-1C44/PyPN/condVel', 'conductionVelocitiesUnmyelinatedCm.dict'), "wb"))

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