import PyPN
import matplotlib.pyplot as plt
import cPickle as pickle
import os

calculationFlag = False # run simulation or load latest bundle with this parameters (not all taken into account for identification)

upstreamSpikingOn = True
electricalStimulusOn = False

# set simulation params
tStop=30
timeRes=0.0025

# set length of bundle and number of axons
lengthOfBundle = 2000
numberOfAxons = 1

# create a guide, the axons will follow
segmentLengthAxon = 10
bundleGuide = PyPN.createGeometry.get_bundle_guide_corner(lengthOfBundle, segmentLengthAxon)

# here the distributions of myelinated and unmyelinated axon diameters are defined
myelinatedDistribution = {
    'densities':[100,300,1150,2750,3650,2850,1750,900,500,250,200,150,110,100,110,100,105], #fibers densities can be given either in No/mm2 or percentage
    'diameters': [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7.,7.5, 8., 8.5, 9.],  # corresponding diameters for each densities
}
unmyelinatedDistribution = {
    'densities':[250,1250,5000,8000,9800,10200,8900,7600,5700,4000,3900,2300,2000,1300,900,750,600,600,500,250], #fibers densities can be given either in No/mm2 or percentage
    'diameters': [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.], # corresponding diameters for each densities
}

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
myelinatedParametersA = {'fiberD': myelinatedDistribution, # um Axon diameter (5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0)
}

# axon parameters
unmyelinatedParameters = {'fiberD': unmyelinatedDistribution, # um Axon diameter
}

# set all properties of the bundle
bundleParameters = {    'radiusBundle': 150, #um Radius of the bundle (typically 0.5-1.5mm)
                        'lengthOfBundle': lengthOfBundle, # um Axon length
                        'numberOfAxons': numberOfAxons, # Number of axons in the bundle
                        'p_A': 0.1, # Percentage of myelinated fiber type A
                        'p_C': 0.9, #Percentage of unmyelinated fiber type C
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

    # create the wanted mechanisms of how to cause spikes on the axons
    # continuous spiking
    if upstreamSpikingOn:
        bundle.add_excitation_mechanism(PyPN.UpstreamSpiking(**upstreamSpikingDict))

    # spiking through a single electrical stimulation
    if electricalStimulusOn:
        bundle.add_excitation_mechanism(PyPN.Stimulus(**stimulusParameters))

    # run the simulation
    bundle.simulate()

    # save the bundle definition file
    bundleSaveLocation = bundle.basePath
    pickle.dump(bundle,open( os.path.join(bundleSaveLocation, 'bundle.cl'), "wb" ))
else:
    saveParams={'elecCount': len(Parameters['recordingElecPos']), 'dt': Parameters['timeRes'], 'tStop': Parameters['tStop'], 'p_A': bundleParameters['p_A'],
                'myelinatedDiam': myelinatedParametersA['fiberD'], 'unmyelinatedDiam': unmyelinatedParameters['fiberD'],
                'L': bundleParameters['lengthOfBundle']}
    bundleSaveLocation = PyPN.get_bundle_directory(new = False, **saveParams)
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