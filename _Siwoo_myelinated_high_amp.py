import PyPN
import matplotlib.pyplot as plt
from PyPN.takeTime import *
import numpy as np
import os

# ------------------------------------------------------------------------------
# ------------------------------- SCRIPT CONTROL -------------------------------
# ------------------------------------------------------------------------------

calculationFlag = True # run simulation or load latest bundle with this parameters (not all taken into account for identification)
electricalStimulusOn = True


# ------------------------------------------------------------------------------
# --------------------------------- DEFINITION ---------------------------------
# ------------------------------------------------------------------------------

# ----------------------------- simulation params ---------------------------

tStop= 50 # 100
timeRes=0.0025

# ----------------------------- bundle params -------------------------------

# set length of bundle and number of axons
lengthOfBundle = 20000 # 20000 # 400000
numberOfAxons = 1

# set the diameter distribution or fixed value
# see http://docs.scipy.org/doc/numpy/reference/routines.random.html
# 5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.

# myelinatedDiam =  1.5 # {'distName' : 'normal', 'params' : (1.0, 0.7)} # (2.0, 0.7)
unmyelinatedDiam = 0.7 # {'distName' : 'normal', 'params' : (0.7, 0.3)}
# unmyelinatedDiam = {'distName' : 'manual', 'params' : {'densities' : (3, 7, 23,    31,    53,    61,    105,   115,   187,   230,   255,   333,   369,   402,   427,   458,   524,   471,   504,   499,   520,   452,   485,   415,   444,   359,   327,   260,   259,   257,   186,   176,   148,   132,   113,   79,    70,    52,    51,    44,    27,    18,    20,    8,11,  7, 9, 2, 7, 5),
#                                            'diameters': (0.7096, 0.7121,    0.7146,    0.7172,    0.7197,    0.722, 0.7248,    0.7273,    0.7298,    0.7324,    0.7349,    0.737, 0.7399,    0.74252,   0.7450,    0.7475,    0.75012,   0.7526,    0.7551,    0.7577,    0.7602,    0.7627,    0.7653,    0.7678,    0.77036,   0.7728,    0.7754,    0.7779,    0.7804,    0.7830,    0.7855,    0.78808,   0.7906,    0.7931,    0.7956,    0.7982,    0.8007,    0.8032,    0.80579,   0.8083,    0.8108,    0.8133,    0.8159,    0.8184,    0.8209,    0.823, 0.826, 0.9, 0.9, 0.9)}}

# unmyelinatedDiam = {'distName' : 'manual', 'params' : {'diameters' : (0.120, 0.17,  0.21,  0.26,  0.32,  0.37,  0.41,  0.47,  0.51,  0.56,  0.62,  0.67,  0.72,  0.77,  0.84,  0.92,  0.97,  1.02,  1.07,  1.12,  1.17,  1.22,  1.27,  1.32, 1.36, 1.41, 1.48, 1.52),
#                                            'densities': (0.0691040631732923, 0.182192465406599, 0.429980837522710, 0.632957475186409, 2.05015339910575,  3.10696898591111,  4.54590886074274,  7.22064649366380,  7.60343269800399,  8.61543655035694,  8.07683524571988,  7.15617584468796,  7.04457675416097,  6.77590492234067,  5.67583310442061,  5.20464797635635,  3.22856301277829,  2.51011904564906,  2.06140597644239,  1.50026642131635,  1.32118496258518,  0.849999834520921, 0.760773515404445, 0.312027350382088, 0.200593738933586, 0.201222559431810)}}

# myelinatedDiam = {'distName' : 'manual', 'params' : {'diameters' : (1.01,    1.19,  1.22,  1.40,  1.41,  1.58,  1.61,  1.78,  1.81,  1.99,  2.01,  2.18,  2.22, 2.39,    2.41,  2.58,  2.61,  2.79,  2.81,  2.99,  3.01,  3.19,  3.61,  3.79,  3.81,  3.99,  4.02,  4.20),
#                                                      'densities' : (0.6553,  0.658, 2.245, 2.282, 4.627, 4.665, 16.734,    16.737,    19.393,    19.396,    17.776,    17.779,    15.503,    15.506,    9.26,  9.234, 5.993, 5.961, 2.272, 2.275, 2.138, 2.106, 1.734, 1.634, 1.151, 1.189, 0.948, 0.917)}}

# unmyelinatedDiam = 0.9

myelinatedDiam = 1.5

# axon definitions
myelinatedParameters = {'fiberD': myelinatedDiam,
                        # 'temperature': 39
                        }
unmyelinatedParameters = {'fiberD': unmyelinatedDiam,
                          # 'temperature': 34
                          }

# bundle guide
segmentLengthAxon = 100
bundleGuide = PyPN.createGeometry.get_bundle_guide_straight(lengthOfBundle, segmentLengthAxon)

# ----------------------------- stimulation params ---------------------------

# parameters of signals for stimulation
rectangularSignalParams = {'amplitude': 0.1,  # Pulse amplitude (mA)
                           'frequency': 1,  #20,     # Frequency of the pulse (kHz)
                           'dutyCycle': 1,  #0.5,    # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                           'stimDur': 0.5,  # Stimulus duration (ms)
                           'waveform': 'BIPHASIC',   # 'MONOPHASIC',  # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                           'delay': 0.,  # ms
                           # 'invert': True,
                           # 'timeRes': timeRes,
                           }


intraParameters = {'stimulusSignal': PyPN.signalGeneration.rectangular(**rectangularSignalParams)}

# ----------------------------- recording params -------------------------------

recordingParametersNew = {'bundleGuide': bundleGuide,
                          'radius': 350,
                          'positionAlongBundle': 8000,
                          'numberOfPoles': 2,
                          'poleDistance': 5000,
                        }

LFPMech1 = PyPN.Extracellular.precomputedFEM(bundleGuide)
LFPMech2 = PyPN.Extracellular.homogeneous(sigma=1)

electrodePos = PyPN.createGeometry.circular_electrode(**recordingParametersNew)

modularRecMech1 = PyPN.RecordingMechanism(electrodePos, LFPMech1)
modularRecMech2 = PyPN.RecordingMechanism(electrodePos, LFPMech2)

# ------------------------------------------------------------------------------
# ---------------------------------- CALCULATION -------------------------------
# ------------------------------------------------------------------------------


if calculationFlag:

    for randomComponent in np.arange(0, 0.1, 0.1):

        # set all properties of the bundle
        bundleParameters = {'radius': 300,  # 150, #um Radius of the bundle (typically 0.5-1.5mm)
                            'length': lengthOfBundle,  # um Axon length
                            'randomDirectionComponent': 0,
                            # 'bundleGuide': bundleGuide,

                            'numberOfAxons': numberOfAxons,  # Number of axons in the bundle
                            'pMyel': 1.,  # Percentage of myelinated fiber type A
                            'pUnmyel': 0.,  # Percentage of unmyelinated fiber type C
                            'paramsMyel': myelinatedParameters,  # parameters for fiber type A
                            'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C

                            'tStop': tStop,
                            'timeRes': timeRes,

                            # 'saveI':True,
                            'saveV': False,
                            'saveLocation': '/media/carl/4ECC-1C44/PyPN/',

                            'numberOfSavedSegments': 50,
                            # number of segments of which the membrane potential is saved to disk
                            # 'downsamplingFactor': 100
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

        bundle.add_recording_mechanism(modularRecMech1)
        bundle.add_recording_mechanism(modularRecMech2)

        # PyPN.plot.geometry_definition(bundle)
        # plt.show()

        # run the simulation
        bundle.simulate()

        # # get SFAP
        # time, CAP = bundle.get_CAP_from_file()
        # plt.plot(time, CAP)

    # save the bundle to disk
    PyPN.save_bundle(bundle)
else:

    # try to open a bundle with the parameters set above
    # bundle = PyPN.open_recent_bundle(bundleParameters)
    # bundle = PyPN.open_bundle_from_location('/media/carl/4ECC-1C44/PyPN/dt=0.0025 tStop=100 pMyel=0.1 pUnmyel=0.9 L=20000 nAxons=500/bundle00000')
    bundle = PyPN.open_bundle_from_location('/home/siwoo/PyCharmProjects/PyPN/28_7/dt=0.0025 tStop=300 pMyel=0.0 pUnmyel=1.0 L=100000 nAxons=15 from Fowley/bundle00000/')

# ------------------------------------------------------------------------------
# ---------------------------------- PLOTTING ----------------------------------
# ------------------------------------------------------------------------------

print '\nStarting to plot'

# first load the desired data from file
fig = plt.figure()
for i in range(len(bundle.recordingMechanisms)):
    time, CAP = bundle.get_CAP_from_file(i)
    # ax = fig.add_subplot(2,1,i)
    plt.plot(time, CAP, label='recMech'+str(i))
    plt.xlabel('time [ms]')
    plt.ylabel('extracellular voltage [mV]')
    # plt.set_title('recMech'+str(i))
plt.legend()


# PyPN.plot.geometry_definition(bundle)
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

PyPN.plot.CAP1D_singleAxon(bundle)
# PyPN.plot.CAP1D(bundle, recMechIndex=1)
# PyPN.plot.voltage(bundle)
# PyPN.plot.diameterHistogram(bundle)
plt.show()

bundle = None