import neuron
from neuron import h
h('load_file("noload.hoc")')
import PyPN
import math
import cPickle as pickle
import matplotlib.pyplot as plt
# from pprint import pprint
# import sys
import time
time.sleep(0) # delays for x seconds

h.celsius = 33 # set temperature in celsius
h.tstop = 30 # set simulation duration (ms)
h.dt = 0.0025 # set time step (ms)
h.finitialize(-65) # initialize voltage state

# decide what's to be done in this run
calculationFlag = True # execute the calculation or load the last run made with this set of parameters

plottingFlag = True # plot something
plotGeometry = True # plot the bundle in space
plotCAP = True # plot the compound action potential
plotCAP1D = True # plot the time course for each electrode location in one separate plot
plotCAP2D = True # if many (>100 electrodes have been used, this plots the CAP in time and space in one image
plotCAP1D_1Axon = True # plot CAP axonwise
plotVoltage = True

# Set parameters

# bundle characteristics
p_A = [.5]#[0.175,0.1,1.0, 0.0] # share of myelinated fibers
fiberD_A = 'draw' # um diameter myelinated axons 'draw' OR one of 5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0
fiberD_C = 'draw' # um diameter unmyelinated axons 'draw' OR some value below 3. (don't forget the dot to make it a float)
radiusBundle = 150.0 # um Radius of the bundle (typically 0.5-1.5mm)
numberOfAxons = 4 # 25
lengthOfBundle = 2000 # 8000

randomDirectionComponent = .1 # curviness of the axon
segmentLengthAxon = 10
bundleGuide = PyPN.createGeometry.get_bundle_guide_cordner(lengthOfBundle, segmentLengthAxon)

# recoding params
numberContactPoints=  8 # Number of points on the circle constituing the cuff electrode
recordingElecPos = [math.floor(lengthOfBundle)]#, math.floor(lengthOfBundle) + 100] #[10000], #Position of the recording electrode along axon in um, in "BIPOLAR" case the position along axons should be given as a couple [x1,x2]
numberElecs =  3 #150#150, #number of electrodes along the bundle

## how to excite axon

# stimulus characteristics
electricalStimulusOn = True
stimTypes = ["INTRA"]#, "INTRA", "EXTRA", "NONE"
waveforms = ["BIPHASIC"]#,"MONOPHASIC", "BIPHASIC"
frequencies = [0.1]#,0.1,0.1]
dutyCycles = [0.01]#[0.001]#,0.01,0.005]
amplitudes = [1.0]#,2.0,0.5]
stimDur = [10]

# upstream nerve stream activity characteristics
upstreamSpikingOn = True
spikingRate = 500. # mean number of pulses per second
spikingCorrelation = 0.1 # pairwise corrleation between neurons
startTimeSpiking = 0
stopTimeSpiking = h.tstop


myelinatedDistribution = {
    'densities':[100,300,1150,2750,3650,2850,1750,900,500,250,200,150,110,100,110,100,105], #fibers densities can be given either in No/mm2 or percentage
    'diameters': [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7.,7.5, 8., 8.5, 9.],  # corresponding diameters for each densities
}

unmyelinatedDistribution = {
    'densities':[250,1250,5000,8000,9800,10200,8900,7600,5700,4000,3900,2300,2000,1300,900,750,600,600,500,250], #fibers densities can be given either in No/mm2 or percentage
    'diameters': [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.], # corresponding diameters for each densities
}

# Do not change from here

if fiberD_A == 'draw':
    del fiberD_A
    fiberD_A = myelinatedDistribution
if fiberD_C == 'draw':
    del fiberD_C
    fiberD_C = unmyelinatedDistribution

upstreamSpikingDict = { 'lambd' : spikingRate, # mean number of pulses per second
                        'correlation' : spikingCorrelation, # pairwise corrleation between neurons
                        'tStart' : startTimeSpiking,
                        'tStop' : stopTimeSpiking,
                        'nAxons' : numberOfAxons}


for j in range(len(dutyCycles)):
    for k in range(len(p_A)):
        stimulusParameters = {
            'stimType': stimTypes[j], #Stimulation type either "INTRA" or "EXTRA"
            'amplitude': amplitudes[j], # Pulse amplitude (nA)
            'frequency': frequencies[j], # Frequency of the sin pulse (kHz)
            'dutyCycle': dutyCycles[j], # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
            'stimDur' : stimDur[j], # Stimulus duration (ms)
            'waveform': waveforms[j], # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
            'radiusBundle' : radiusBundle,
        }

        recordingParameters = {
            "numberContactPoints": numberContactPoints, #Number of points on the circle constituing the cuff electrode
            'recordingElecPos': recordingElecPos,#[10000], #Position of the recording electrode along axon in um, in "BIPOLAR" case the position along axons should be given as a couple [x1,x2]
            'numberElecs': numberElecs,#150, #number of electrodes along the bundle
        }

        myelinatedParametersA = {
            'name': "myelinated_axonA", # axon name (for neuron)
            #'Nnodes': 11, #Number of nodes
            'fiberD': fiberD_A, #fiberD_A, #Diameter of the fiber
            'layout3D': "PT3D",#"DEFINE_SHAPE",# # either "DEFINE_SHAPE" or "PT3D" using hoc corresponding function
            'rec_v': True, # set voltage recorders True or False
        }
        unmyelinatedParameters = {
            'name': "unmyelinated_axon", # axon name (for neuron)
            'fiberD': fiberD_C, #Axon diameter (micrometer)
            'cm' : 1.0, #Specific membrane capacitance (microfarad/cm2)
            'Ra': 200.0, #Specific axial resistance (Ohm cm)
            'layout3D': "PT3D",#"DEFINE_SHAPE", # either "DEFINE_SHAPE" or "PT3D" using hoc corresponding function
            'rec_v': True, # set voltage recorders True or False
        }

        bundleParameters = {         #parameters for Bundle classe
            'radiusBundle': radiusBundle, #um Radius of the bundle (typically 0.5-1.5mm)
            'lengthOfBundle': lengthOfBundle, #Axon length (micrometer)
            'numberOfAxons': numberOfAxons,#640, # Number of axons in the bundle
            'p_A': p_A[k], # Percentage of myelinated fiber type A
            'p_C': 1-p_A[k], #Percentage of unmyelinated fiber type C
            'myelinated_A': myelinatedParametersA, #parameters for fiber type A
            'unmyelinated': unmyelinatedParameters, #parameters for fiber type C
            'randomDirectionComponent' : randomDirectionComponent,
            'segmentLengthAxon' : segmentLengthAxon,
            'bundleGuide' : bundleGuide
        }


        Parameters = dict(bundleParameters, **recordingParameters)

        saveParams={'elecCount': len(recordingElecPos), 'dt': h.dt, 'tStop': h.tstop, 'p_A': bundleParameters['p_A'],
                'myelinatedDiam': myelinatedParametersA['fiberD'], 'unmyelinatedDiam': unmyelinatedParameters['fiberD'],
                'L': bundleParameters['lengthOfBundle']}
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
            bundle.simulate_bundle()

            # save the bundle definition file
            bundleSaveLocation = bundle.basePath
            pickle.dump(bundle,open( bundleSaveLocation+'bundle.cl', "wb" ))
        else:
            bundleSaveLocation = PyPN.get_bundle_directory(new = False, **saveParams)
            try:
                bundle = pickle.load(open(bundleSaveLocation+'bundle.cl', "rb" ))
            except:
                print 'No bundle with these parameters has been generated yet. Set calculationFlag to True.'
                quit()
            # bundle = None

        # pprint (vars(bundle))
        # # for axon in bundle.axons:
        # #     # print axon.nbytes
        # #     pprint (vars(axon))
        # pprint(vars(bundle.axons[0]))

        if plottingFlag:

            if plotGeometry:

                PyPN.plotBundleClass.plot_geometry(bundle)

            if plotCAP:

                if plotCAP1D_1Axon:

                    PyPN.plotBundleClass.plot_CAP1D_singleAxon(bundle, 10)

                if plotCAP1D:

                    PyPN.plotBundleClass.plot_CAP1D(bundle)

                if plotCAP2D:

                    PyPN.plotBundleClass.plot_CAP2D(bundle)

            if plotVoltage:

                PyPN.plotBundleClass.plot_voltage(bundle)

        if plottingFlag:
            # finally show result
            plt.show()

bundle = None


