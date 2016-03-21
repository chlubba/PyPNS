import PyPN
import matplotlib.pyplot as plt
import numpy as np

# set simulation params
tStop=30
timeRes=0.0025

# set length of bundle and number of axons
lengthOfBundle = 2000
numberOfAxons = 30

# create a guide, the axons will follow
segmentLengthAxon = 10
bundleGuide = PyPN.createGeometry.get_bundle_guide_corner(lengthOfBundle, segmentLengthAxon)

myelinatedDiam = {'distName' : 'normal',
                  'params' : (6,2)}

unmyelinatedDiam = 5.

# definition of properties of spontaeous spiking on axons
upstreamSpikingDict = { 'lambd' : 500, # mean number of pulses per second
                        'correlation' : 0.1, # pairwise corrleation between neurons
                        'tStart' : 0,
                        'tStop' : tStop,
                        'nAxons' : numberOfAxons
}

# definition of the stimulation type of the axon
stimulusParameters = {  'stimType': "INTRA", #Stimulation type either "INTRA" or "EXTRA"
                        'amplitude': 1.5, # Pulse amplitude (nA)
                        'frequency': 0.1, # Frequency of the pulse (kHz)
                        'dutyCycle': 0.01, # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
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
myelinatedParametersA = {'fiberD': myelinatedDiam, # um Axon diameter (5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0) or dict
}

# axon parameters
unmyelinatedParameters = {'fiberD': unmyelinatedDiam, # um Axon diameter or dict defining the distribution
}

# set all properties of the bundle
bundleParameters = {    'radiusBundle': 150, #um Radius of the bundle (typically 0.5-1.5mm)
                        'lengthOfBundle': lengthOfBundle, # um Axon length
                        'numberOfAxons': numberOfAxons, # Number of axons in the bundle
                        'p_A': 0., # Percentage of myelinated fiber type A
                        'p_C': 1., #Percentage of unmyelinated fiber type C
                        'myelinated_A': myelinatedParametersA, #parameters for fiber type A
                        'unmyelinated': unmyelinatedParameters, #parameters for fiber type C
                        'bundleGuide' : bundleGuide,
                        'tStop' : tStop,
                        'timeRes' : timeRes
}

# combine parameters for the bundle creation
Parameters = dict(bundleParameters, **recordingParameters)

# create the bundle with all properties of axons and recording setup
bundle = PyPN.Bundle(**Parameters)

# # create the wanted mechanisms of how to cause spikes on the axons
# # continuous spiking
# bundle.add_excitation_mechanism(PyPN.UpstreamSpiking(**upstreamSpikingDict))

# spiking through a single electrical stimulation
bundle.add_excitation_mechanism(PyPN.Stimulus(**stimulusParameters))

# run the simulation
bundle.simulate()

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