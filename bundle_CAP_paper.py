from classes import *
from classes_plot import *
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



### !!! USING MYELINATED DISTRUTION NOT COMPATIBLE WITH VOLTAGE RECORDING ###
### ANYWAY VOLTAGE RECORDING SHOULD LOGICALLY BE PERFORMED FOR ONLY ONE AXON AT A TIME ###
myelinatedDistribution = {
    'densities':[100,300,1150,2750,3650,2850,1750,900,500,250,200,150,110,100,110,100,105], #fibers densities can be given either in No/mm2 or percentage
    'diameters': [ 1., 1.5, 2., 2.5, 3., 3.5, 4., 4.5, 5., 5.5, 6., 6.5, 7.,7.5, 8., 8.5, 9.],  # corresponding diameters for each densities
}

unmyelinatedDistribution = {
    'densities':[250,1250,5000,8000,9800,10200,8900,7600,5700,4000,3900,2300,2000,1300,900,750,600,600,500,250], #fibers densities can be given either in No/mm2 or percentage
    'diameters': [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.], # corresponding diameters for each densities
}



#fiberD from 5.7, 7.3, 8.7, 10.0, 11.5, 12.8, 14.0, 15.0, 16.0
fiberD_A = 16.0 #um

### PARAMETERS BELOW ARE ACTUALLY SET IN THE CODE BUT USEFUL FOR THE HEADER AND PLOTS ###
## Values of nodelength and paralength1 are constant
## Values of paralength2 and interlength are to be set according to the chosen fiberD value
nodelength = 1.0 #um
paralength1 = 1.3 #um
[paralength2_A, interlength_A] = fiberD_dependent_param(fiberD_A, nodelength, paralength1)








p_A = [0.175,0.1,1.0, 0.0]
for j in range(1):
    for k in range(len(p_A)):#[0]:#range(1,len(p_A)):#
        stimulusParameters = {
            'jitter_para': [0,0], #Mean and standard deviation of the delay
            'stim_type': "EXTRA", #Stimulation type either "INTRA" or "EXTRA"
            # stim_coord is NOT USED default value is used directly
            'stim_coord': [[0,50,0]], # spatial coordinates  of the stimulating electrodes, example for tripolar case=[[xe0,ye0,ze0], [xe1,ye1,ze1], [xe2,ye2,ze2]] (order is important with the middle being the cathode), INTRA case only use the position along x for IClamp
            'amplitude': 2.0, # Pulse amplitude (nA)
            'freq': 0.1, # Frequency of the sin pulse (kHz)
            'duty_cycle': 0.005, # Percentage stimulus is ON for one period (t_ON = duty_cyle*1/f)
            'stim_dur' : 10, # Stimulus duration (ms)
            'waveform': "BIPHASIC", # Type of waveform either "MONOPHASIC" or "BIPHASIC" symmetric
}

        recordingParameters = {
            "number_contact_points": 8, #Number of points on the circle constituing the cuff electrode
            'recording_elec_pos': [1000],#[10000], #Position of the recording electrode along axon in um, in "BIPOLAR" case the position along axons should be given as a couple [x1,x2]
            'number_elecs': 50,#150, #number of electrodes along the bundle
            'dur': h.tstop, # Simulation duration (ms)
            'rec_CAP': True, #If false means we avoid spending time using LFPy functions
        }
        myelinatedParametersA = {
            'name': "myelinated_axonA", # axon name (for neuron)
            'Nnodes': 11, #Number of nodes
            'fiberD': myelinatedDistribution, #fiberD_A, #Diameter of the fiber
            'layout3D': "DEFINE_SHAPE", # either "DEFINE_SHAPE" or "PT3D" using hoc corresponding function
            'rec_v': True, # set voltage recorders True or False
            # 'nodelength' : nodelength, #set node length (um)
            # 'paralength1': paralength1, #set the length of the nearest paranode segment to the node
            # 'paralength2': paralength2_A,  #set the length of the second paranode segment followed by the internodes segments
            # 'interlength': interlength_A, #set the length of the internode part comprising the 6 segments between two paranodes2
        }


        unmyelinatedParameters = {
            'name': "unmyelinated_axon", # axon name (for neuron)
            'L': 1000,#10000, #Axon length (micrometer)
            'diam': unmyelinatedDistribution, #Axon diameter (micrometer)
            'cm' : 1.0, #Specific membrane capacitance (microfarad/cm2)
            'Ra': 200.0, #Specific axial resistance (Ohm cm)
            'layout3D': "DEFINE_SHAPE", # either "DEFINE_SHAPE" or "PT3D" using hoc corresponding function
            'rec_v': True, # set voltage recorders True or False
        }
        bundleParameters = {         #parameters for Bundle classe
            'radius_bundle': 150.0, #um Radius of the bundle (typically 0.5-1.5mm)
            'draw_distribution': True, #Boolean stating the distribution of fibre should be drawn
            'number_of_axons': 50,#640, # Number of axons in the bundle
            'p_A': p_A[k], # Percentage of myelinated fiber type A
            'p_C': 1-p_A[k], #Percentage of unmyelinated fiber type C
            'myelinated_A': myelinatedParametersA, #parameters for fiber type A
            'unmyelinated': unmyelinatedParameters, #parameters for fiber type C
        }


        Parameters1 = dict(bundleParameters, **stimulusParameters)
        Parameters = dict(Parameters1, **recordingParameters)

        # saveParams={'elecCount': len(Parameters['recording_elec_pos']), 'dt': h.dt, 'p_A': Parameters['p_A'], 'p_C': Parameters['p_C'], 'L': unmyelinatedParameters['L'] }
        # directory = getDirectoryName("CAP", **saveParams)
        directory = "FOR_PAPER/CAP2D/recordings/"
        if not os.path.exists(directory):
            os.makedirs(directory)
        bundle = Bundle(**Parameters)

        # When saving voltage to file limit the number of axons to 10. If 100 unmyelinated it produces a 1Go file, and if 100 myelinated 2Go.
    
        save_CAP_tofile(bundle,Parameters,directory)
        # save_voltage_tofile(bundle,Parameters,directory)
        bundle = None
        






