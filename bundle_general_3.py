from neuron import h
h('load_file("noload.hoc")')
from bundleClass import *
import cPickle as pickle
from pprint import pprint
import sys
import time
time.sleep(0) # delays for x seconds

""" McIntyre Model parameters
celsius=37			
v_init=-80 //mV//  		
dt=0.005 //ms//         	
tstop=10"""
h.celsius = 33 # set temperature in celsius
h.tstop = 30 # set simulation duration (ms)
h.dt = 0.0025 # set time step (ms)
h.finitialize(-65) # initialize voltage state

# Set parameters
calculationFlag = True

plottingFlag = True

plotGeometry = True

plotCAP = True
plotCAP1D = True
plotCAP2D = False
plotCAP1D_1Axon = True

plotVoltage = True

# bundle characteristics
p_A = [0.]#[0.175,0.1,1.0, 0.0] # share of myelinated fibers
fiberD_A = 'draw' #5.7# 16.0 #um diameter myelinated axons 'draw' OR one of 5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0
fiberD_C = 'draw' #1.5 #'draw'
myelinatedCurviness = 0.314


radius_bundle = 150.0 #um Radius of the bundle (typically 0.5-1.5mm)
draw_distribution = True #Boolean stating the distribution of fibre should be drawn
number_of_axons = 50#25
lengthOfBundle = 5000#8000


# stimulus characteristics
stim_types = ["EXTRA"]#, "INTRA", "EXTRA", "NONE"
waveforms = ["BIPHASIC"]#,"MONOPHASIC", "BIPHASIC"
frequencies = [0.1]#,0.1,0.1]
duty_cycles = [0.01]#[0.001]#,0.01,0.005]
amplitudes = [2.0]#,2.0,0.5]
stimDur = [10]

# upstream nerve stream activity characteristics
upstreamSpikingOn = True
spikingRate = 500. # mean number of pulses per second
spikingCorrelation = 0.1 # pairwise corrleation between neurons
startTimeSpiking = 0
stopTimeSpiking = h.tstop

# recoding params
number_contact_points=  8 #Number of points on the circle constituing the cuff electrode
recording_elec_pos = [math.floor(lengthOfBundle*0.7), math.floor(lengthOfBundle*0.7) + 100] #[10000], #Position of the recording electrode along axon in um, in "BIPOLAR" case the position along axons should be given as a couple [x1,x2]
number_elecs =  1#150#150, #number of electrodes along the bundle

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

upstreamSpikingDict = { #'upstreamSpikingOn' : upstreamSpikingOn,
                        'lambd' : spikingRate, # mean number of pulses per second
                        'correlation' : spikingCorrelation, # pairwise corrleation between neurons
                        'tStart' : startTimeSpiking,
                        'tStop' : stopTimeSpiking}



for j in range(len(duty_cycles)):
    for k in range(len(p_A)):#[0]:#range(1,len(p_A)):#
        stimulusParameters = {
            'jitter_para': [0,0], #Mean and standard deviation of the delay
            'stim_type': stim_types[j], #Stimulation type either "INTRA" or "EXTRA"
            # stim_coord is NOT USED default value is used directly
            'stim_coord': [[0,50,0]], # spatial coordinates  of the stimulating electrodes, example for tripolar case=[[xe0,ye0,ze0], [xe1,ye1,ze1], [xe2,ye2,ze2]] (order is important with the middle being the cathode), INTRA case only use the position along x for IClamp
            'amplitude': amplitudes[j], # Pulse amplitude (nA)
            'freq': frequencies[j], # Frequency of the sin pulse (kHz)
            'duty_cycle': duty_cycles[j], # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
            'stim_dur' : stimDur[j], # Stimulus duration (ms)
            'waveform': waveforms[j], # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
        }

        recordingParameters = {
            "number_contact_points": number_contact_points, #Number of points on the circle constituing the cuff electrode
            'recording_elec_pos': recording_elec_pos,#[10000], #Position of the recording electrode along axon in um, in "BIPOLAR" case the position along axons should be given as a couple [x1,x2]
            'number_elecs': number_elecs,#150, #number of electrodes along the bundle
            'dur': h.tstop, # Simulation duration (ms)
            'rec_CAP': True, #If false means we avoid spending time using LFPy functions
        }
        myelinatedParametersA = {
            'name': "myelinated_axonA", # axon name (for neuron)
            'Nnodes': 11, #Number of nodes
            'fiberD': fiberD_A, #fiberD_A, #Diameter of the fiber
            'layout3D': "DEFINE_SHAPE",#"PT3D",# # either "DEFINE_SHAPE" or "PT3D" using hoc corresponding function
            'rec_v': True, # set voltage recorders True or False
            # 'nodelength' : nodelength, #set node length (um)
            # 'paralength1': paralength1, #set the length of the nearest paranode segment to the node
            # 'paralength2': paralength2_A,  #set the length of the second paranode segment followed by the internodes segments
            # 'interlength': interlength_A, #set the length of the internode part comprising the 6 segments between two paranodes2
        }


        unmyelinatedParameters = {
            'name': "unmyelinated_axon", # axon name (for neuron)
            'L': lengthOfBundle,#Axon length (micrometer)
            'diam': fiberD_C, #Axon diameter (micrometer)
            'cm' : 1.0, #Specific membrane capacitance (microfarad/cm2)
            'Ra': 200.0, #Specific axial resistance (Ohm cm)
            'layout3D': "PT3D",#"DEFINE_SHAPE", # either "DEFINE_SHAPE" or "PT3D" using hoc corresponding function
            'rec_v': True, # set voltage recorders True or False
        }
        bundleParameters = {         #parameters for Bundle classe
            'radius_bundle': radius_bundle, #um Radius of the bundle (typically 0.5-1.5mm)
            'draw_distribution': draw_distribution, #Boolean stating the distribution of fibre should be drawn
            'number_of_axons': number_of_axons,#640, # Number of axons in the bundle
            'p_A': p_A[k], # Percentage of myelinated fiber type A
            'p_C': 1-p_A[k], #Percentage of unmyelinated fiber type C
            'myelinated_A': myelinatedParametersA, #parameters for fiber type A
            'unmyelinated': unmyelinatedParameters, #parameters for fiber type C
        }


        Parameters1 = dict(bundleParameters, **stimulusParameters)
        Parameters = dict(Parameters1, **recordingParameters)

        saveParams={'elecCount': len(recording_elec_pos), 'dt': h.dt, 'tStop': h.tstop, 'p_A': bundleParameters['p_A'],
                'myelinatedDiam': myelinatedParametersA['fiberD'], 'unmyelinatedDiam': unmyelinatedParameters['diam'],
                'L': unmyelinatedParameters['L'], 'stimType': stimulusParameters['stim_type'], 'stimWaveform' : stimulusParameters['waveform'],
                'stimDutyCycle': stimulusParameters['duty_cycle'], 'stimAmplitude' : stimulusParameters['amplitude']}



        # bundleDirectory = getDirectoryName('bundle', **saveParams)
        # if not os.path.exists(bundleDirectory):
        #     os.makedirs(bundleDirectory)

        if calculationFlag:

            bundle = Bundle(**Parameters)

            if upstreamSpikingOn:
                bundle.addUpstreamSpiking(**upstreamSpikingDict)

            bundle.simulateBundle()

            # save the whole bundle
            # bundleSaveLocation = getFileName("bundle", saveParams)
            bundleSaveLocation = bundle.basePath
            pickle.dump(bundle,open( bundleSaveLocation+'bundle.cl', "wb" ))
        else:
            #directory = getDirectoryName("bundle", **saveParams)
            bundleSaveLocation = getBundleDirectory(new = False, **saveParams)
            try:
                bundle = pickle.load(open(bundleSaveLocation+'bundle.cl', "rb" ))
            except:
                print 'No bundle with these parameters has been generated yet. Set calculationFlag to True.'
                quit()
            # bundle = None

        # pprint (vars(bundle))
        # print bundle.nbytes
        # for axon in bundle.axons:
        #     # print axon.nbytes
        #     pprint (vars(axon))
        # pprint(vars(bundle.axons[0]))

        if plottingFlag:

            if plotGeometry:

                bundle.plot_geometry()

            if plotCAP:

                if plotCAP1D_1Axon:

                    bundle.plot_CAP1D_singleAxon(10)

                if plotCAP1D:

                    bundle.plot_CAP1D()

                if plotCAP2D:

                    bundle.plot_CAP2D()

            if plotVoltage:

                bundle.plot_voltage()

        if plottingFlag:
            # finally show result
            plt.show()

bundle = None


