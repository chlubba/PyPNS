import PyPN
import matplotlib.pyplot as plt

# import cPickle as pickle
# import os

calculationFlag = True # run simulation or load latest bundle with this parameters (not all taken into account for identification)

upstreamSpikingOn = False
electricalStimulusOn = True

# set simulation params
tStop=50
timeRes=0.005#0.0025

# set length of bundle and number of axons
lengthOfBundle = 10000
numberOfAxons = 1

# create a guide the axons will follow
segmentLengthAxon = 10
# bundleGuide = PyPN.createGeometry.get_bundle_guide_corner(lengthOfBundle, segmentLengthAxon)
bundleGuide = PyPN.createGeometry.get_bundle_guide_straight(lengthOfBundle, segmentLengthAxon)

# set the diameter distribution or fixed value
# see http://docs.scipy.org/doc/numpy/reference/routines.random.html
# 5.7, 7.3, 8.7, 10., 11.5, 12.8, 14., 15., 16.
myelinatedDiam =  10.2 # 0.3 # {'distName' : 'uniform', 'params' : (5.1, 16)} # .2 #
unmyelinatedDiam =  {'distName' : 'uniform', 'params' : (0.1, 20)} # .2 #

# definition of the stimulation type of the axon
stimulusParameters = {  'stimType': "INTRA", #Stimulation type either "INTRA" or "EXTRA"
                        'amplitude': 1.5, #0.2, # 0.004, # 10., #  # Pulse amplitude (nA)
                        'frequency': 20., # Frequency of the pulse (kHz)
                        'dutyCycle': 1., # 0.05, # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                        'stimDur' : 0.05, # Stimulus duration (ms)
                        'waveform': 'MONOPHASIC', # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                        'radiusBundle' : 150, #um
                        # 'tStop' : tStop,
                        'timeRes' : timeRes,
                        'delay': 5 # ms
}
# # definition of the stimulation type of the axon
# stimulusParameters = {  'delay': 5, # delay (ms)
#                         'stimDur': 0.05, # Stimulus duration (ms)
#                         'amplitude': 1.5 # 0.15 # Pulse amplitude (nA)
# }

# recording parameters of the cuff electrodes
recordingParameters = { 'numberContactPoints': 8, # Number of points on the circle constituing the cuff electrode
                        'recordingElecPos': [lengthOfBundle, lengthOfBundle + 50], #um Position of the recording electrode along axon in um, in "BIPOLAR" case the position along axons should be given as a couple [x1,x2]
                        'numberElecs': 1, # number of electrodes along the bundle
}

# axon parameters
myelinatedParametersA = {'fiberD': myelinatedDiam, # um Axon diameter (5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0)
}

# axon parameters
unmyelinatedParameters = {'fiberD': unmyelinatedDiam, # um Axon diameter
}

# set all properties of the bundle
bundleParameters = {    'radiusBundle': 150, #um Radius of the bundle (typically 0.5-1.5mm)
                        'lengthOfBundle': lengthOfBundle, # um Axon length
                        'numberOfAxons': numberOfAxons, # Number of axons in the bundle
                        'p_A': 1., # Percentage of myelinated fiber type A
                        'p_C': 0., #Percentage of unmyelinated fiber type C
                        'myelinated_A': myelinatedParametersA, #parameters for fiber type A
                        'unmyelinated': unmyelinatedParameters, #parameters for fiber type C
                        'bundleGuide' : bundleGuide,
                        'tStop' : tStop,
                        'timeRes' : timeRes,
                        'randomDirectionComponent' : 0
}

# combine parameters for the bundle creation
Parameters = dict(bundleParameters, **recordingParameters)

if calculationFlag:

    # create the bundle with all properties of axons and recording setup
    bundle = PyPN.Bundle(**Parameters)

    # spiking through a single electrical stimulation
    if electricalStimulusOn:
        stimulusInstance = PyPN.Stimulus(**stimulusParameters)
        bundle.add_excitation_mechanism(PyPN.Stimulus(**stimulusParameters))
        # bundle.add_excitation_mechanism(PyPN.SimpleIClamp(**stimulusParameters))

    # run the simulation
    bundle.simulate()

    # save the bundle to disk
    PyPN.save_bundle(bundle)
else:

    # try to open a bundle with the parameters set above
    bundle = PyPN.open_recent_bundle(Parameters)

# # plot geometry, intra and extracellular recording, axon diameters
# print '\nStarting to plot'
# PyPN.plot.geometry(bundle)
# PyPN.plot.CAP1D_singleAxon(bundle, 10)
# PyPN.plot.CAP1D(bundle)
PyPN.plot.voltage(bundle)
# PyPN.plot.voltage_one_myelinated_axon(bundle)
# PyPN.plot.diameterHistogram(bundle)

# conVelDict = bundle.conduction_velocities(saveToFile=True) # (plot=False)
# pickle.dump(conVelDict,open( os.path.join(bundle.basePath, 'conductionVelocities.dict'), "wb" ))


# import matplotlib2tikz as mtz
# mtz.save('CAP.tex')

plt.show()

bundle = None