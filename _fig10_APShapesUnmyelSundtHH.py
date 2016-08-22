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

tStop=100
timeRes=0.0025

# ----------------------------- bundle params -------------------------------

# set length of bundle and number of axons
lengthOfBundle = 5000 # 20000 # 400000
numberOfAxons = 1

# bundle guide
segmentLengthAxon = 30
bundleGuide = PyPN.createGeometry.get_bundle_guide_straight(lengthOfBundle, segmentLengthAxon)

# ----------------------------- stimulation params ---------------------------

# parameters of signals for stimulation
rectangularSignalParams = {'amplitude': 100., #50,  # Pulse amplitude (mA)
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

diameters = np.flipud(np.arange(.2, 4., .5))
diameter = 0.8
temperatures = np.arange(5, 46, 5)
Ras = [70] # np.arange(50, 300, 50)
Cms = [0.1, 0.3, 1.] # np.arange(0.15, 0.2, 0.01)

# saveDict = {'axonType': 'unmyelinated',
#             'diameters': diameters,
#             'temperatures': temperatures,
#             'velocityArray': []}

if calculationFlag:

    firstAP = []


    LFPMech = PyPN.Extracellular.homogeneous(sigma=1)

    electrodePos = PyPN.createGeometry.circular_electrode(**recordingParametersNew)

    modularRecMech = PyPN.RecordingMechanism(electrodePos, LFPMech)

    # set the diameter distribution or fixed value
    # see http://docs.scipy.org/doc/numpy/reference/routines.random.html
    # 5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.
    myelinatedDiam = diameter  # {'distName' : 'normal', 'params' : (1.0, 0.7)} # (2.0, 0.7)
    unmyelinatedDiam = diameter  # {'distName' : 'normal', 'params' : (0.7, 0.3)}

    # axon definitions
    myelinatedParameters = {'fiberD': myelinatedDiam} # , 'temperature': temperature}
    unmyelinatedParameters = {'fiberD': unmyelinatedDiam} # , 'cm': Cm, 'Ra': Ra} # , 'temperature': temperature}

    # set all properties of the bundle
    bundleParameters = {'radius': 300,  # 150, #um Radius of the bundle (typically 0.5-1.5mm)
                        'length': lengthOfBundle,  # um Axon length
                        # 'randomDirectionComponent': .9,
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

    # bundle.add_recording_mechanism(modularRecMech1)
    bundle.add_recording_mechanism(modularRecMech)
    # bundle.add_recording_mechanism(modularRecMech3)

    # PyPN.plot.geometry_definition(bundle)
    # plt.show()

    # run the simulation
    bundle.simulate()


    # get the membrane voltage
    t, v = bundle.get_voltage_from_file_one_axon(0)

    # regularize time
    from scipy.interpolate import interp1d
    f = interp1d(t, v[:,30])
    tReg = np.arange(0,max(t),0.0025)
    vReg = f(tReg)

    # plt.subplot(len(Ras), len(Cms), RaInd * len(Cms) + CmInd)

    plt.plot(tReg, vReg, label='Hodgkin-Huxley')

    data = np.loadtxt("Fig5doutput_2.dat")

    numrows = len(data)
    numcols = len(data[0])

    x = data[:, 0]
    peri0 = data[:, 1]

    # adapt time step
    from scipy.interpolate import interp1d
    f = interp1d(x, peri0)
    # tRegSundt = np.arange(0, max(x), 0.0025)
    tRegSundtInt = np.arange(50, 150, 0.0025)
    tRegSundt = tRegSundtInt - 50
    vRegSundt = f(tRegSundtInt)

    correlation = np.correlate(vReg - np.mean(vReg), vRegSundt - np.mean(vRegSundt), 'full')
    # plt.plot(correlation)
    # plt.show()
    lag = len(vReg) - np.argmax(correlation)

    # plt.plot(tReg, firstAP)
    # plt.plot(tReg, vReg)
    # plt.title('lag' + str(lag*0.0025))
    # plt.show()


    # # plt.subplot(4,1,1)
    # plt.plot(x, peri0)  # distal peripheral
    # plt.ylabel('Potential (mV)')
    # plt.xlabel('Time (ms)')

    plt.plot(tRegSundt[:-lag], vRegSundt[lag:], label='Sundt et al. 2015')
    # plt.plot(tRegSundt[:-50/0.0025], vRegSundt[50/0.0025:], label='Sundt et al. 2015')
    # plt.plot(tRegSundt, vRegSundt, label='Sundt et al. 2015')
    # plt.plot(tReg, vReg)

    plt.xlim((9,15))
    plt.xlabel('time [ms]')
    plt.ylabel('$V_m$ [mV]')
    plt.legend()

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