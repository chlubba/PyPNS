import PyPN
import matplotlib.pyplot as plt
from PyPN.takeTime import *
import numpy as np
import os
from pprint import pprint
import sys
import cPickle as pickle
from scipy.interpolate import interp1d

# ------------------------------------------------------------------------------
# ------------------------------- SCRIPT CONTROL -------------------------------
# ------------------------------------------------------------------------------

calculationFlag = False # run simulation or load latest bundle with this parameters (not all taken into account for identification)
electricalStimulusOn = True


# ------------------------------------------------------------------------------
# --------------------------------- DEFINITION ---------------------------------
# ------------------------------------------------------------------------------

# ----------------------------- simulation params ---------------------------

tStop=70 # 50
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

diametersUnmyel = [1] # np.arange(0.2, 2, 0.3)
diametersMyel = [2] # [2.3, 2.6, 2.9] # np.arange(0.2, 4, 0.3)
diametersBothTypes = [diametersUnmyel, diametersMyel]

tStartPlots = [0.2, 0.05]

diameters = np.arange(0.2, 4, 0.7) # [1, 2, 4] # [4, 3, 2, 1, 0.5, 0.2] # np.flipud(np.arange(.2, 4., .3))
temperatures = np.arange(5, 46, 5)
Ras = np.arange(50, 300, 50)

RDCs = [0, 0.2, 0.4, 0.6, 0.8, 1.] # np.arange(0, 1., 0.15)

if calculationFlag:

    (f, axarr) = plt.subplots(1, 2, sharey=True)

    legends = ['Unmyelinated', 'Myelinated']
    bundleLengths = [40000, 40000] # [40000, 40000] # [60000, 40000]
    for i in [0,1]:

        vAPCollection = []

        diameters = diametersBothTypes[i]

        recMechLegends = ['homogeneous', 'FEM Oil 2cm inhomo'] # , 'FEM Oil only .0000001', 'FEM Oil contact']
        recMechMarkers = ['o', 'v']
        for recMechIndex in [1]:

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
                                    'randomDirectionComponent': 0,
                                    # 'bundleGuide': bundleGuide,

                                    'numberOfAxons': numberOfAxons,  # Number of axons in the bundle
                                    'pMyel': i,  # Percentage of myelinated fiber type A
                                    'pUnmyel': 1 - i,  # Percentage of unmyelinated fiber type C
                                    'paramsMyel': myelinatedParameters,  # parameters for fiber type A
                                    'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C

                                    'tStop': tStop,
                                    'timeRes': 0.0025, #'variable', #

                                    'saveI':True,
                                    'saveV': False,
                                    'saveLocation': '/media/carl/SANDISK/PyPN/',

                                    'numberOfSavedSegments': 50,
                                    # number of segments of which the membrane potential is saved to disk
                                    }

                # create the bundle with all properties of axons and recording setup
                bundle = PyPN.Bundle(**bundleParameters)

                # ----------------------------- stimulation params ---------------------------

                amplitudes = (50., 5.5)

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

                # run the simulation
                bundle.simulate()

                PyPN.save_bundle(bundle)
else:

    # try to open a bundle with the parameters set above
    # bundle = PyPN.open_recent_bundle(bundleParameters)
    # bundle = PyPN.open_bundle_from_location('/media/carl/4ECC-1C44/PyPN/dt=0.0025 tStop=100 pMyel=0.1 pUnmyel=0.9 L=20000 nAxons=500/bundle00000')
    # bundle = PyPN.open_bundle_from_location('/media/carl/4ECC-1C44/PyPN/dt=0.0025 tStop=100 pMyel=0.1 pUnmyel=0.9 L=20000 nAxons=500/bundle00001')
    # bundle = PyPN.open_bundle_from_location(
    #     '/media/carl/4ECC-1C44/PyPN/dt=0.0025 tStop=60 pMyel=0.0 pUnmyel=1.0 L=5000 nAxons=13/bundle00000')

    def smooth(y, box_pts):
        box = np.ones(box_pts) / box_pts
        y_smooth = np.convolve(y, box, mode='same')
        return y_smooth

    (f, axarr) = plt.subplots(3, 3, sharey='col', sharex='col')

    # bundleLocations = ['/media/carl/SANDISK/PyPN/dt=0.0025 tStop=70 pMyel=0 pUnmyel=1 L=40000 nAxons=1/bundle00013',
    #                    '/media/carl/SANDISK/PyPN/dt=0.0025 tStop=70 pMyel=1 pUnmyel=0 L=40000 nAxons=1/bundle00003']
    # bundleLocations = ['/Volumes/SANDISK/PyPN/dt=0.0025 tStop=70 pMyel=0 pUnmyel=1 L=40000 nAxons=1/bundle00016',
    #                    '/Volumes/SANDISK/PyPN/dt=0.0025 tStop=70 pMyel=1 pUnmyel=0 L=40000 nAxons=1/bundle00005']
    bundleLocations = ['/media/carl/SANDISK/PyPN/dt=0.0025 tStop=70 pMyel=0 pUnmyel=1 L=40000 nAxons=1/bundle00018',
                       '/media/carl/SANDISK/PyPN/dt=0.0025 tStop=70 pMyel=1 pUnmyel=0 L=40000 nAxons=1/bundle00006']

    legends = ['Unmyelinated', 'Myelinated']
    bundleLengths = [40000, 40000] # [40000, 40000] # [60000, 40000]
    for i in [0,1]:

        bundle = PyPN.open_bundle_from_location(bundleLocations[i])

        smoothingLength = [1, 100, 1000] # [0.01, 0.001, 0.0001]  # , 0.00001]
        xPs = [0, 0.000090, 0.000180]
        FEMFileNames = ['smoothed1_zeroed.npy', 'smoothed100_zeroed.npy', 'smoothed1000_zeroed.npy']
        profileParams = [smoothingLength, xPs, FEMFileNames]

        for profileInd in [1]: # [0, 1, 2]

            import matplotlib.cm as cm
            import matplotlib.colors as colors
            jet = plt.get_cmap('jet')
            cNorm = colors.Normalize(vmin=0, vmax=len(smoothingLength) - 1)
            scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)

            LFPMechs = []
            modularRecMechs = []
            for profileParamInd, profileParam in enumerate(profileParams[profileInd]):

                # define a function translating z-position of a 1 nA current point source to voltage
                zValues = np.linspace(-0.015, 0.015, 7500)

                if profileInd == 0:
                    profileTitle = 'idealized cuff field. smoothing boxwidth in um'
                    profileLegends = (np.array(smoothingLength) * 4).astype(str)

                    vValues = np.maximum(0, np.multiply(np.subtract(1, np.abs(zValues / 0.01)), 8.83e-5))

                    # smooth to compare to FEM output smoothed
                    vValues = smooth(vValues, profileParam)
                    # vValues = vValues / np.max(vValues) * 8.83e-5


                elif profileInd == 1:

                    profileTitle = 'idealized field close to axon segment. legend: factors.'
                    profileLegends = (np.array(profileParams[profileInd]) * 1.5).astype(str)

                    angle = 0
                    xP = profileParam

                    peakFactor = max(0, (1 - np.mod(angle, np.pi) / np.pi * 5)) * np.min([1, (xP / 0.000190) ** 5])
                    a = 1.9E-9  # 2.5E-9;
                    b = 0.00005

                    # zValues = np.arange(0, 1.05, 0.01) / 100
                    vValues = a * (1.0 / (np.abs(zValues) + b)) * peakFactor + np.maximum(0, (np.subtract(1, np.abs(zValues / 0.01)) * 8.83e-5))

                    # a = 2.5e-9
                    # b = 0.00005 # 0.00025
                    #
                    # peakFactors = profileParams[profileInd]
                    #
                    # vValues = a * (1. / (np.abs(zValues) + b)) * peakFactors[profileParamInd] + np.maximum(0, np.multiply(np.subtract(1, np.abs(zValues / 0.01)), 8.83e-5))

                    # plt.close('all')
                    # plt.plot(zValues, vValues)
                    # plt.plot(zValues, np.maximum(0, np.multiply(np.subtract(1, np.abs(zValues / 0.01)), 8.83e-5)))
                    # plt.show()


                elif profileInd == 2:

                    profileTitle = 'FEM output, smoothed over x um'
                    smoothWidth = [1, 100, 1000]
                    profileLegends = (np.array(smoothWidth) * 4).astype(str)

                    smoothedFunctionArrayPacked = np.load(
                        os.path.join(
                            '/Volumes/SANDISK/ComsolData/interpolated_function_electrode_fixed_a_z_variable',
                            profileParam))
                    smoothedFunctionArray = smoothedFunctionArrayPacked[()]
                    zValues = smoothedFunctionArray[0]
                    vValues = smoothedFunctionArray[1]

                if i == 0:
                    colorVal = scalarMap.to_rgba(profileParamInd)
                    axarr[profileInd][0].plot(zValues, vValues*1000, label=profileLegends[profileParamInd], color=colorVal)
                    axarr[profileInd][0].set_title('phi over z')
                    axarr[profileInd][0].set_xlabel('z (m)')
                    axarr[profileInd][0].set_ylabel('V (uV)')
                    axarr[profileInd][0].grid()

                zInterpolator = interp1d(zValues, vValues, bounds_error=False, fill_value=0)

                LFPMechs.append(PyPN.Extracellular.interpolator(bundle.bundleCoords, method='z')) #  interpolator=zInterpolator,

                if i == 1:

                    relPosition = 0
                    recordingParametersNew = {'bundleGuide': bundle.bundleCoords,
                                              'radius': 220,
                                              'positionAlongBundle': np.floor(bundleLengths[i]*0.6 / bundle.axons[0].lengthOneCycle) *
                                                                     bundle.axons[0].lengthOneCycle + bundle.axons[0].lengthOneCycle*relPosition,
                                              'numberOfPoles': 1,
                                              'poleDistance': 2000,
                                              }

                else:

                    recordingParametersNew = {'bundleGuide': bundle.bundleCoords,
                                              'radius': 220,
                                              'positionAlongBundle': 11000,
                                              'numberOfPoles': 1,
                                              'poleDistance': 2000,
                                              }


                electrodePos = PyPN.createGeometry.circular_electrode(**recordingParametersNew)
                modularRecMechs.append(PyPN.RecordingMechanism(electrodePos, LFPMechs[-1]))

                bundle.add_recording_mechanism(modularRecMechs[-1])

            axarr[profileInd][0].legend(loc='best')
            axarr[profileInd][0].set_title(profileTitle)


            bundle.compute_CAPs_from_imem_files()

            # ------------------------------------------------------------------------------
            # ---------------------------------- PLOTTING ----------------------------------
            # ------------------------------------------------------------------------------

            for recMechIndex in range(len(bundle.recordingMechanisms)):

                colorVal = scalarMap.to_rgba(recMechIndex)

                t, SFAPs = bundle.get_SFAPs_from_file(recMechIndex)
                axarr[profileInd][i+1].plot(t, SFAPs*1000, color=colorVal)

            axarr[profileInd][i+1].set_xlabel('time (ms)')
            axarr[profileInd][i+1].set_ylabel('$V_{ext}$ (uV)')
            axarr[profileInd][i+1].grid()

            titles = ['unmyelinated', 'myelinated']
            axarr[profileInd][i+1].set_title(titles[i])

            if i == 0:
                axarr[profileInd][i+1].set_xlim([0, 60]) # [24,32]
                axarr[profileInd][i + 1].set_ylim([-10, 3])
                # axarr[profileInd][i + 1].legend(loc='best')
            elif i == 1:
                axarr[profileInd][i+1].set_xlim([0, 7])
                # axarr[profileInd][i + 1].set_ylim([-.35, 0.2])

            # delete what has been calculated before
            # todo: make this function just parse all folders with the right name and delete them.
            bundle.clear_all_recording_mechanisms()



    plt.show()

bundle = None