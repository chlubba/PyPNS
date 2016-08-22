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
diameter = 5.7
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

    # for temperatureInd, temperature in enumerate(temperatures):
    for gkbarInd, gkbar in enumerate(gkbars):


        # set the diameter distribution or fixed value
        # see http://docs.scipy.org/doc/numpy/reference/routines.random.html
        # 5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.
        myelinatedDiam = diameter  # {'distName' : 'normal', 'params' : (1.0, 0.7)} # (2.0, 0.7)
        unmyelinatedDiam = diameter  # {'distName' : 'normal', 'params' : (0.7, 0.3)}

        # axon definitions
        myelinatedParameters = {'fiberD': myelinatedDiam, 'gkbar_axnode': gkbar} # , 'temperature': temperature}
        unmyelinatedParameters = {'fiberD': unmyelinatedDiam} # , 'temperature': temperature}

        # set all properties of the bundle
        bundleParameters = {'radius': 300,  # 150, #um Radius of the bundle (typically 0.5-1.5mm)
                            'length': lengthOfBundle,  # um Axon length
                            # 'randomDirectionComponent': .9,
                            # 'bundleGuide': bundleGuide,

                            'numberOfAxons': numberOfAxons,  # Number of axons in the bundle
                            'pMyel': 1.,  # Percentage of myelinated fiber type A
                            'pUnmyel': 0.,  # Percentage of unmyelinated fiber type C
                            'paramsMyel': myelinatedParameters,  # parameters for fiber type A
                            'paramsUnmyel': unmyelinatedParameters,  # parameters for fiber type C

                            'tStop': tStop,
                            'timeRes': 0.0025, #'variable', #

                            # 'saveI':True,
                            # 'saveV': False,
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

        # run the simulation
        bundle.simulate()


        # membrane voltage approach

        t, v = bundle.get_voltage_from_file_one_axon(0)

        nodeIndex = np.floor(bundle.axons[0].totnsegs/11*0.8)*11

        from scipy.interpolate import interp1d
        f = interp1d(t, v[:,nodeIndex])

        tReg = np.arange(0,max(t),0.0025)
        vReg = f(tReg)

        plt.plot(tReg, vReg, label=str(gkbar))

        # plt.subplot(len(Ras), len(Cms), RaInd * len(Cms) + CmInd)

        # saveDict['velocityArray'].append(np.array(vAPs))

        # import matplotlib.cm as cm
        # import matplotlib.colors as colors
        # jet = plt.get_cmap('jet')
        # cNorm = colors.Normalize(vmin=0, vmax=len(temperatures)-1)
        # scalarMap = cm.ScalarMappable(norm=cNorm, cmap=jet)
        # colorVal = scalarMap.to_rgba(CmInd)
        # plt.plot(diameters, vAPs, label=str(Cm) + ' $\mu$F/cm$^2$', color=colorVal)

            # # save the bundle to disk
            # PyPN.save_bundle(bundle)

    # plt.plot(np.arange(0,3.7,0.2), np.sqrt(np.arange(0,3.7,0.2))*2, linestyle='--', color=np.array([1, 1, 1])*0.7, label='theory')
    # plt.xlabel('diameter [um]')
    # plt.ylabel('conduction velocity [m/s]')
    # plt.title('Unmyelinated Axon with Ra = 50 Ohm cm')
    # plt.legend(loc='best')
    # plt.grid()

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