import PyPN
import matplotlib.pyplot as plt
import cPickle as pickle
import os

calculationFlag = True # run simulation or load latest bundle with this parameters (not all taken into account for identification)

upstreamSpikingOn = False
electricalStimulusOn = True

# set simulation params
tStop=30
timeRes=0.0025

# set length of bundle and number of axons
lengthOfBundle = 2000
numberOfAxons = 50

# create a guide, the axons will follow
segmentLengthAxon = 10
bundleGuide = PyPN.createGeometry.get_bundle_guide_corner(lengthOfBundle, segmentLengthAxon)

# set the diameter distribution or fixed value
myelinatedDiam = .2 # {'distName' : 'gamma', 'params' : (6,2)}
unmyelinatedDiam = 0.1

# definition of the stimulation type of the axon
stimulusParameters = {  'stimType': "INTRA", #Stimulation type either "INTRA" or "EXTRA"
                        'amplitude': 0.05, # Pulse amplitude (nA)
                        'frequency': 0.1, # Frequency of the pulse (kHz)
                        'dutyCycle': 0.005, # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
                        'stimDur' : 10, # Stimulus duration (ms)
                        'waveform': 'MONOPHASIC', # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
                        'radiusBundle' : 150, #um
                        'tStop' : tStop,
                        'timeRes' : timeRes
}

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
                        'timeRes' : timeRes
}

# combine parameters for the bundle creation
Parameters = dict(bundleParameters, **recordingParameters)

if calculationFlag:

    # create the bundle with all properties of axons and recording setup
    bundle = PyPN.Bundle(**Parameters)

    # spiking through a single electrical stimulation
    if electricalStimulusOn:
        bundle.add_excitation_mechanism(PyPN.Stimulus(**stimulusParameters))

    # run the simulation
    bundle.simulate()

    # save the bundle definition file
    bundleSaveLocation = bundle.basePath
    pickle.dump(bundle,open( os.path.join(bundleSaveLocation, 'bundle.cl'), "wb" ))
else:
    bundleSaveLocation = PyPN.get_bundle_directory(Parameters, new = False)
    try:
        bundle = pickle.load(open(os.path.join(bundleSaveLocation, 'bundle.cl'), "rb" ))
    except:
        print 'No bundle with these parameters has been generated yet. Set calculationFlag to True.'
        quit()

# plot geometry, voltage in axon and extracellular recording
print '\nStarting to plot'
PyPN.plot.geometry(bundle)
PyPN.plot.CAP1D_singleAxon(bundle, 10)
PyPN.plot.CAP1D(bundle)
PyPN.plot.voltage(bundle)

# import matplotlib2tikz as mtz
# mtz.save('CAP.tex')

plt.show()

bundle = None