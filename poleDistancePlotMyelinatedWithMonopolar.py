import PyPN
import matplotlib.pyplot as plt
from PyPN.takeTime import *
import numpy as np

# import cPickle as pickle
import os

calculationFlag = True # run simulation or load latest bundle with this parameters (not all taken into account for identification)

upstreamSpikingOn = False
electricalStimulusOn = True

# set simulation params
tStop=20
timeRes=0.0025#0.0025

# set length of bundle and number of axons
lengthOfBundle = 10000
numberOfAxons = 10 # 1000

# set the diameter distribution or fixed value
# see http://docs.scipy.org/doc/numpy/reference/routines.random.html
# 5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.
myelinatedDiam =  2.75 # {'distName' : 'normal', 'params' : (2.0, 0.7)}
unmyelinatedDiam = 1. # {'distName' : 'normal', 'params' : (0.7, 0.3)}


myelinatedDiams = np.arange(0.2, 4, 0.3)

vppsMono = []

for myelinatedDiam in myelinatedDiams:

    intraParameters = {     'amplitude': 4., # 0.005, # 0.016,#0.2,# .0001,#1.5, #0.2, # 0.004, # 10., #  # Pulse amplitude (mA)
                            'frequency': 20., # Frequency of the pulse (kHz)
                            'dutyCycle': .5, # 0.05, # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                            'stimDur' : 0.05, # Stimulus duration (ms)
                            'waveform': 'MONOPHASIC', # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                            'timeRes' : timeRes,
                            'delay': 5, # ms
                            # 'invert': True
    }

    # axon parameters
    myelinatedParameters = {'fiberD': myelinatedDiam, # um Axon diameter (5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0)
    }

    # axon parameters
    unmyelinatedParameters = {'fiberD': unmyelinatedDiam, # um Axon diameter
    }

    # set all properties of the bundle
    bundleParameters = {    'radius': 150, #150, #um Radius of the bundle (typically 0.5-1.5mm)
                            'length': lengthOfBundle, # um Axon length
                            # 'bundleGuide' : bundleGuide,
                            'randomDirectionComponent' : 0,

                            'numberOfAxons': numberOfAxons, # Number of axons in the bundle
                            'pMyel': 1.,#0.01, # Percentage of myelinated fiber type A
                            'pUnmyel': 0., #.99, #Percentage of unmyelinated fiber type C
                            'paramsMyel': myelinatedParameters, #parameters for fiber type A
                            'paramsUnmyel': unmyelinatedParameters, #parameters for fiber type C

                            'tStop' : tStop,
                            'timeRes' : timeRes,

                            # 'saveI':True,
                            'saveV':False
    }


    if calculationFlag:

        # create the bundle with all properties of axons and recording setup
        bundle = PyPN.Bundle(**bundleParameters)

        # spiking through a single electrical stimulation
        if electricalStimulusOn:
            # stimulusInstance = PyPN.Stimulus(**stimulusParameters)
            # plt.plot(stimulusInstance.t, stimulusInstance.stimulusSignal)
            # plt.title('stimulus signal without delay')
            # plt.show()
            # bundle.add_excitation_mechanism(PyPN.StimCuff(**cuffParameters))
            bundle.add_excitation_mechanism(PyPN.StimIntra(**intraParameters))
            # bundle.add_excitation_mechanism(PyPN.SimpleIClamp(**stimulusParameters))

        # bundle.add_recording_mechanism(PyPN.RecCuff3D(radius=500, positionMax=0.8, sigma=1., width=2000, distanceOfRings=100, pointsPerRing=20))
        # bundle.add_recording_mechanism(PyPN.RecCuff3D(radius=500, positionMax=0.1, sigma=1., width=2000, distanceOfRings=100, pointsPerRing=20))
        # bundle.add_recording_mechanism(PyPN.RecCuff2D(radius=500, positionMax=0.2, sigma=5.))

        for poleDistance in np.arange(10, 1000, 20):
            bundle.add_recording_mechanism(PyPN.RecCuff2D(radius=200, numberOfPoles=2, poleDistance=poleDistance, positionMax=0.9, sigma=1.))

        # run the simulation
        bundle.simulate()

        # PyPN.plot.voltage(bundle)
        # plt.show()

        # save the bundle to disk
        PyPN.save_bundle(bundle)
    else:

        # try to open a bundle with the parameters set above
        bundle = PyPN.open_recent_bundle(bundleParameters)
        # bundle = PyPN.open_bundle_from_location('/media/carl/4ECC-1C44/PyPN/dt=0.0025 tStop=20 pMyel=1.0 pUnmyel=0.0 L=4000 nAxons=2/bundle00010')

    # bundle.clear_all_recording_mechanisms()
    # # bundle.add_recording_mechanism(PyPN.RecCuff2D(radius=200, positionMax=0.2, sigma=1.))
    # # investigate effect of bipolar electrode distance
    # for poleDistance in np.arange(10, 1000, 20):
    #     bundle.add_recording_mechanism(PyPN.RecCuff2D(radius=200, numberOfPoles=2, poleDistance=poleDistance, positionMax=0.2, sigma=1.))
    #
    # bundle.compute_CAPs_from_imem_files()
    #
    # # save the bundle to disk
    # PyPN.save_bundle(bundle)

    vpps = []
    distances = []
    for recMechIndex in range(len(bundle.recordingMechanisms)-1):

        recMech = bundle.add_recording_mechanism(recMechIndex)

        poleDistance = recMech.poleDistance
        distances.append(poleDistance)

        time, CAP = bundle.get_CAP_from_file(recordingMechanismIndex=recMechIndex)
        afterStim = CAP[:,time > 5.2]

        vpp = np.max(afterStim) - np.min(afterStim)
        vpps.append(vpp)


    titleString = 'Radius Bundle %i um, radius 2D electrode %i um,\n%i straight unmyelinated axons of diameter %.3f um' % (bundle.radiusBundle, bundle.recordingMechanisms[0].radius, numberOfAxons, unmyelinatedDiam)
    plt.title(titleString)
    plt.plot(distances, vpps)
    plt.xlabel('pole distance [um]')
    plt.ylabel('Vpp [mV]')

    saveString = 'Radius Bundle %i um, radius 2D electrode %i um, %i straight myelinated axons of diameter %.3f um' % (
    bundle.radiusBundle, bundle.recordingMechanisms[0].radius, numberOfAxons, myelinatedDiam)

    np.savetxt(os.path.join('/media/carl/4ECC-1C44/PyPN/poleDistanceAgainstVpp/Data', saveString+'_vpp.txt'), vpps)
    np.savetxt(os.path.join('/media/carl/4ECC-1C44/PyPN/poleDistanceAgainstVpp/Data', saveString + '_distances.txt'), distances)

    plt.savefig(os.path.join('/media/carl/4ECC-1C44/PyPN/poleDistanceAgainstVpp/Fig', saveString+'.png'))

    time, CAP = bundle.get_CAP_from_file(recordingMechanismIndex=len(bundle.recordingMechanisms)-1)
    afterStim = CAP[:, time > 5.2]

    vpp = np.max(afterStim) - np.min(afterStim)

    vppsMono.append(vpp)

plt.figure()
# titleString = 'Radius Bundle %i um, radius 2D electrode %i um,\n%i straight unmyelinated axons of diameter %.3f um' % (bundle.radiusBundle, bundle.recordingMechanisms[0].radius, numberOfAxons, unmyelinatedDiam)
plt.title('peak to peak voltage against diameter')
plt.semilogy(myelinatedDiams, vpps)
plt.xlabel('Diameter [um]')
plt.ylabel('Vpp [mV]')

filename = 'vpp_against_diameter_mono_myelinated'

np.savetxt(os.path.join('/media/carl/4ECC-1C44/PyPN/poleDistanceAgainstVpp/Data', filename+'.txt'), vppsMono)

plt.savefig(os.path.join('/media/carl/4ECC-1C44/PyPN/poleDistanceAgainstVpp/Fig', filename+'.png'))

plt.show()

bundle = None